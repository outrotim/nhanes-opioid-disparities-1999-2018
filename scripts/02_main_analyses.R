# ============================================================
# 02_main_analyses.R — Q1 Oaxaca-Blinder + Q2 ITS + Q4 trend
# Repository: nhanes-opioid-disparities-1999-2018
# Q1: Opioid / gabapentinoid / SNRI / TCA / NSAID-Rx prescription trends
# Q2: CDC 2016 ITS (cycle-level, 2 yr granularity)
# Q3: Opioid days x MetS RCS
# Q4: Race/SES disparity in opioid Rx (Oaxaca-Blinder + RII)
# ============================================================
source(file.path(Sys.getenv("REPO_ROOT", unset=here::here()), "environment.R"))
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(tidyr); library(survey); library(srvyr)
  library(rms); library(ggplot2); library(broom)
})

d <- arrow::read_parquet(file.path(DIR_POOLED, "analytic.parquet"))
cat("N total =", nrow(d), " weighted N =", sum(d$wt_comb, na.rm=TRUE), "\n")

d$cancer_excl <- replace_na(d$had_cancer, 0L)

des_all <- svydesign(id=~psu, strata=~strata, weights=~wt_comb,
                     nest=TRUE, data=d |> filter(!is.na(psu), !is.na(strata)))
options(survey.lonely.psu="adjust")

# ==========================================================
# Q1: Prescription trends per class over 10 cycles
# ==========================================================
cat("\n========== Q1: Trends (weighted prevalence) ==========\n")
classes <- c("opioid_any","gaba_any","snri_any","tca_any","nsaid_rx_any")
q1_tab <- purrr::map_dfr(classes, function(cls) {
  f <- as.formula(paste0("~I(", cls, "*1)"))
  byf <- as.formula("~cycle_yr")
  r <- svyby(f, byf, des_all, svymean, na.rm=TRUE)
  tibble(class=cls, cycle_yr=r$cycle_yr, prev=r[,2], se=r$se)
})
q1_tab <- q1_tab |> mutate(lo=pmax(0, prev - 1.96*se), hi=prev + 1.96*se)
readr::write_csv(q1_tab, file.path(DIR_RESULTS, "Q1_trends.csv"))

q1_trend_tests <- purrr::map_dfr(classes, function(cls) {
  ff <- as.formula(paste0(cls, " ~ cycle_yr"))
  m  <- svyglm(ff, design=des_all, family=quasibinomial())
  tidy(m) |> filter(term=="cycle_yr") |>
    mutate(class=cls, or=exp(estimate), or_lo=exp(estimate-1.96*std.error),
           or_hi=exp(estimate+1.96*std.error))
})
readr::write_csv(q1_trend_tests, file.path(DIR_RESULTS, "Q1_trend_tests.csv"))
print(q1_trend_tests)

# ==========================================================
# Q2: ITS — cycle-level segmented regression around CDC 2016
# ==========================================================
cat("\n========== Q2: ITS CDC 2016 ==========\n")
its_dat <- q1_tab |>
  filter(class %in% c("opioid_any","gaba_any","snri_any")) |>
  mutate(time  = cycle_yr - min(cycle_yr),
         post  = as.integer(cycle_yr >= 2015),
         time_post = pmax(0, cycle_yr - 2015) )

its_fits <- its_dat |>
  group_by(class) |>
  group_modify(~ tidy(lm(prev ~ time + post + time_post, data=.x,
                         weights=1/pmax(.x$se^2, 1e-6))))
readr::write_csv(its_fits, file.path(DIR_RESULTS, "Q2_ITS_coeffs.csv"))
print(its_fits |> dplyr::select(class, term, estimate, std.error, p.value), n=Inf)

# ==========================================================
# Q3: RCS — opioid days vs MetS
# ==========================================================
cat("\n========== Q3: RCS opioid_days x MetS ==========\n")
q3 <- d |> filter(!is.na(mets), !is.na(age), !is.na(sex_f), !is.na(race_f), !is.na(pir),
                  cancer_excl == 0)
cat("Q3 analytic N =", nrow(q3), "\n")

q3 <- q3 |> mutate(
  opi_cat = case_when(
    !opioid_any ~ "None",
    opioid_days < 30 ~ "Short",
    opioid_days < 90 ~ "Medium",
    TRUE ~ "Chronic90+"),
  opi_cat = factor(opi_cat, levels=c("None","Short","Medium","Chronic90+"))
)
des3 <- svydesign(id=~psu, strata=~strata, weights=~wt_comb, nest=TRUE, data=q3)

m3_cat <- svyglm(mets ~ opi_cat + age + sex_f + race_f + pir + bmi,
                 design=des3, family=quasibinomial())
q3_tab <- tidy(m3_cat) |>
  mutate(or=exp(estimate), or_lo=exp(estimate-1.96*std.error),
         or_hi=exp(estimate+1.96*std.error))
readr::write_csv(q3_tab, file.path(DIR_RESULTS, "Q3_opioid_cat_mets_OR.csv"))
print(q3_tab |> filter(grepl("opi_cat", term)) |>
        dplyr::select(term, or, or_lo, or_hi, p.value))

# ==========================================================
# Q4: Disparity in opioid Rx receipt
#   Oaxaca-Blinder + RII/SII, time-varying cycle-by-cycle
# ==========================================================
cat("\n========== Q4: Racial/SES disparity ==========\n")

race_trend <- d |> filter(!is.na(race_f)) |>
  group_by(cycle_yr, race_f) |>
  summarise(prev = weighted.mean(opioid_any, wt_comb, na.rm=TRUE),
            n = n(), .groups="drop")
readr::write_csv(race_trend, file.path(DIR_RESULTS, "Q4_opioid_by_race_cycle.csv"))

pir_trend <- d |> filter(!is.na(pir_q)) |>
  group_by(cycle_yr, pir_q) |>
  summarise(prev = weighted.mean(opioid_any, wt_comb, na.rm=TRUE),
            .groups="drop")
pir_rii <- pir_trend |> group_by(cycle_yr) |>
  summarise(rii = prev[pir_q==1] / prev[pir_q==5],
            sii = prev[pir_q==1] - prev[pir_q==5])
readr::write_csv(pir_rii, file.path(DIR_RESULTS, "Q4_pir_rii_sii.csv"))
print(pir_rii)

# Oaxaca-Blinder: Black vs White (NH), pooled across all cycles
if (requireNamespace("oaxaca", quietly=TRUE)) {
  library(oaxaca)
  ox_dat <- d |> filter(race_f %in% c("NHWhite","NHBlack"),
                        !is.na(pir), !is.na(age), !is.na(educ_f), !is.na(bmi)) |>
    mutate(group = as.integer(race_f == "NHBlack"),
           opioid_bin = as.integer(opioid_any))
  ox_fit <- tryCatch(
    oaxaca(opioid_bin ~ age + pir + I(educ_f=="Coll+") + bmi + insured | group,
           data=ox_dat, R=300),
    error=function(e) NULL)
  if (!is.null(ox_fit)) {
    saveRDS(ox_fit, file.path(DIR_RESULTS, "Q4_oaxaca.rds"))
    cat("\nOaxaca-Blinder (Black vs White), pooled:\n")
    print(ox_fit$twofold$overall)
  }

  # ---------- Time-varying Oaxaca: one fit per cycle ----------
  cat("\n---------- Time-varying Oaxaca (per cycle) ----------\n")
  cycles <- sort(unique(ox_dat$cycle_yr))
  ox_time <- purrr::map_dfr(cycles, function(y) {
    sub <- ox_dat |> filter(cycle_yr == y)
    if (nrow(sub) < 500) return(NULL)
    fit <- tryCatch(
      oaxaca(opioid_bin ~ age + pir + I(educ_f=="Coll+") + bmi + insured | group,
             data=sub, R=100),
      error=function(e) NULL)
    if (is.null(fit)) return(NULL)
    tf <- fit$twofold$overall
    i <- which.min(abs(tf[,"group.weight"] - 0.5))  # Reimers pooled
    tibble(cycle_yr=y,
           explained=tf[i,"coef(explained)"],
           explained_se=tf[i,"se(explained)"],
           unexplained=tf[i,"coef(unexplained)"],
           unexplained_se=tf[i,"se(unexplained)"])
  })
  readr::write_csv(ox_time, file.path(DIR_RESULTS, "Q4_oaxaca_by_cycle.csv"))
  print(ox_time)

  # Inverse-variance-weighted linear slope for unexplained component
  if (nrow(ox_time) >= 5) {
    w <- 1 / pmax(ox_time$unexplained_se^2, 1e-8)
    slope_fit <- lm(unexplained ~ cycle_yr, data=ox_time, weights=w)
    cat("\nIVW slope of unexplained component:\n")
    print(tidy(slope_fit))
  }
}

# Pre-/post-CDC 2016 test on unexplained component
if (exists("ox_time") && nrow(ox_time) >= 5) {
  ox_time$post2016 <- as.integer(ox_time$cycle_yr >= 2015)
  pp_fit <- lm(unexplained ~ post2016, data=ox_time,
               weights=1/pmax(ox_time$unexplained_se^2, 1e-8))
  cat("\nPre-/post-CDC-2016 change in unexplained gap:\n")
  print(tidy(pp_fit))
}

# ==========================================================
# Summary Table 1 (weighted)
# ==========================================================
cat("\n========== Table 1 (weighted) ==========\n")
t1 <- d |>
  group_by(cycle_yr) |>
  summarise(n=n(),
            age_mean = weighted.mean(age, wt_comb, na.rm=TRUE),
            female_pct = weighted.mean(sex==2, wt_comb, na.rm=TRUE)*100,
            white_pct = weighted.mean(race_f=="NHWhite", wt_comb, na.rm=TRUE)*100,
            black_pct = weighted.mean(race_f=="NHBlack", wt_comb, na.rm=TRUE)*100,
            hisp_pct  = weighted.mean(race_f=="Hispanic", wt_comb, na.rm=TRUE)*100,
            pir_mean  = weighted.mean(pir, wt_comb, na.rm=TRUE),
            opioid_pct= weighted.mean(opioid_any, wt_comb, na.rm=TRUE)*100,
            gaba_pct  = weighted.mean(gaba_any, wt_comb, na.rm=TRUE)*100,
            mets_pct  = weighted.mean(mets, wt_comb, na.rm=TRUE)*100)
readr::write_csv(t1, file.path(DIR_TABLES, "table1_weighted.csv"))
print(t1 |> dplyr::select(cycle_yr, n, opioid_pct, gaba_pct, mets_pct))

cat("\nOUTPUTS saved under results/ figures/ tables/\n")
