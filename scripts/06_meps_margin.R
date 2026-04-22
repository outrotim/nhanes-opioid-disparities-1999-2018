################################################################################
# 06_meps_margin.R
#
# MEPS 2015 + 2016 pooled extensive vs chronic vs short-course margin
# decomposition for adult opioid prescription receipt.
#
# Margins:
#   Extensive   = any opioid fill in calendar year (cumulative days > 0)
#   Chronic     = cumulative days-of-supply >= 90
#   Short-course = cumulative days-of-supply > 0 and < 90
#
# Primary estimand: share of the Black-White extensive-margin gap that is
# attributable to (a) chronic component vs (b) short-course component, with
# 1,000-replicate person-level bootstrap percentile confidence intervals.
#
# Produces numerical inputs for Exhibit 1 in the v2.0 manuscript.
#
# Requires four MEPS files in data/raw/meps/:
#   h181.ssp   2015 Full-Year Consolidated   (HC-181)
#   h178a.ssp  2015 Prescribed Medicines     (HC-178A)
#   h192.ssp   2016 Full-Year Consolidated   (HC-192)
#   h188a.ssp  2016 Prescribed Medicines     (HC-188A)
#
# Outputs:
#   results/meps_margin_decomposition.csv  margin × race summary
#   results/meps_margin_bootstrap.rds      1,000-rep bootstrap output
################################################################################

suppressMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(survey)
})
options(survey.lonely.psu = "adjust")
set.seed(20260421)

MEPS_DIR <- "data/raw/meps"

# ---- Opioid identification (case-insensitive RxName matching) ----
opioid_re <- paste0(
  "(HYDROCODONE|OXYCODONE|HYDROMORPH|MORPHINE|FENTANYL|CODEINE|METHADONE|",
  "TRAMADOL|OXYMORPH|TAPENTADOL|MEPERIDINE|BUPRENORPH|PERCOCET|VICODIN|",
  "NORCO|LORTAB|LORCET|HYDROCO/APAP|OXYC/APAP|OXYC/ACE|PENTAZOCINE|",
  "NUCYNTA|BELBUCA|BUTRANS|KADIAN|EXALGO|DILAUDID|MS CONTIN|AVINZA|",
  "EMBEDA|DURAGESIC|ZOHYDRO|OPANA|CONZIP|DARVOCET|DARVON|ULTRACET|",
  "ULTRAM|DEMEROL|BUTORPHAN|NALBUPHINE|TALWIN|STADOL|SUBOXONE|SUBUTEX)"
)
exclude_re <- "NALOXONE"

# ---- Build person-year table for one MEPS year ----
build_py <- function(hc_path, rx_path, age_var, perwt_var, year_label) {
  hc <- read_xpt(hc_path)
  rx <- read_xpt(rx_path)

  rx$is_opioid <- as.integer(
    grepl(opioid_re, rx$RXNAME, ignore.case = TRUE) &
    !grepl(exclude_re, rx$RXNAME, ignore.case = TRUE)
  )
  rx$RXDAYSUP_clean <- ifelse(rx$RXDAYSUP < 0, NA_real_, rx$RXDAYSUP)

  person_rx <- rx %>%
    filter(is_opioid == 1) %>%
    group_by(DUPERSID) %>%
    summarise(total_days = sum(RXDAYSUP_clean, na.rm = TRUE),
              .groups = "drop")

  hc_sel <- hc %>%
    rename(age_v = !!age_var, perwt_v = !!perwt_var) %>%
    filter(age_v >= 20) %>%
    select(DUPERSID, age_v, RACETHX, perwt_v, VARPSU, VARSTR) %>%
    mutate(race_f = case_when(
      RACETHX == 2 ~ "NHWhite",
      RACETHX == 3 ~ "NHBlack",
      TRUE ~ "Other"
    ))

  hc_sel %>%
    left_join(person_rx, by = "DUPERSID") %>%
    mutate(
      total_days = ifelse(is.na(total_days), 0, total_days),
      any_opioid = as.integer(total_days > 0),
      is_chronic = as.integer(total_days >= 90),
      is_short   = as.integer(any_opioid == 1 & total_days < 90),
      year       = year_label
    )
}

py15 <- build_py(file.path(MEPS_DIR, "h181.ssp"),
                 file.path(MEPS_DIR, "h178a.ssp"),
                 "AGE15X", "PERWT15F", "2015")
py16 <- build_py(file.path(MEPS_DIR, "h192.ssp"),
                 file.path(MEPS_DIR, "h188a.ssp"),
                 "AGE16X", "PERWT16F", "2016")

# Pool with combined weight scaled by number of years to preserve population
py15$wt_pool <- py15$perwt_v / 2
py16$wt_pool <- py16$perwt_v / 2
py_all <- bind_rows(py15, py16) %>%
  filter(race_f %in% c("NHBlack", "NHWhite"))

cat(sprintf("Pooled MEPS 2015+2016 NHB+NHW: %d (NHB=%d, NHW=%d)\n",
            nrow(py_all),
            sum(py_all$race_f == "NHBlack"),
            sum(py_all$race_f == "NHWhite")))

# ---- Weighted prevalences via survey design ----
des <- svydesign(ids = ~VARPSU, strata = ~VARSTR, weights = ~wt_pool,
                 nest = TRUE, data = py_all)

mean_by_race <- function(formula) {
  out <- svyby(formula, ~race_f, des, svymean, na.rm = TRUE)
  c(NHB = as.numeric(out[out$race_f == "NHBlack", 2]),
    NHB_SE = sqrt(as.numeric(vcov(out)[1, 1])),
    NHW = as.numeric(out[out$race_f == "NHWhite", 2]),
    NHW_SE = sqrt(as.numeric(vcov(out)[2, 2])))
}

ext  <- mean_by_race(~any_opioid)
chr  <- mean_by_race(~is_chronic)
sho  <- mean_by_race(~is_short)

# ---- Bootstrap shares ----
B <- 1000
boot <- replicate(B, {
  idx <- sample(nrow(py_all), replace = TRUE)
  pb  <- py_all[idx, ]
  b   <- pb %>% filter(race_f == "NHBlack")
  w   <- pb %>% filter(race_f == "NHWhite")
  g_ext <- weighted.mean(b$any_opioid, b$wt_pool) -
           weighted.mean(w$any_opioid, w$wt_pool)
  g_chr <- weighted.mean(b$is_chronic, b$wt_pool) -
           weighted.mean(w$is_chronic, w$wt_pool)
  g_sho <- weighted.mean(b$is_short, b$wt_pool) -
           weighted.mean(w$is_short, w$wt_pool)
  c(g_ext = g_ext, g_chr = g_chr, g_sho = g_sho,
    prop_chr = g_chr / g_ext, prop_sho = g_sho / g_ext)
})
boot_df <- as.data.frame(t(boot))

summary_tbl <- data.frame(
  margin                   = c("extensive", "chronic", "short_course"),
  NHB                      = c(ext["NHB"], chr["NHB"], sho["NHB"]),
  NHB_SE                   = c(ext["NHB_SE"], chr["NHB_SE"], sho["NHB_SE"]),
  NHW                      = c(ext["NHW"], chr["NHW"], sho["NHW"]),
  NHW_SE                   = c(ext["NHW_SE"], chr["NHW_SE"], sho["NHW_SE"]),
  gap_BminusW_pp           = c((ext["NHB"]-ext["NHW"])*100,
                               (chr["NHB"]-chr["NHW"])*100,
                               (sho["NHB"]-sho["NHW"])*100),
  gap_lo_pp                = c(quantile(boot_df$g_ext, 0.025, na.rm=TRUE)*100,
                               quantile(boot_df$g_chr, 0.025, na.rm=TRUE)*100,
                               quantile(boot_df$g_sho, 0.025, na.rm=TRUE)*100),
  gap_hi_pp                = c(quantile(boot_df$g_ext, 0.975, na.rm=TRUE)*100,
                               quantile(boot_df$g_chr, 0.975, na.rm=TRUE)*100,
                               quantile(boot_df$g_sho, 0.975, na.rm=TRUE)*100),
  share_of_extensive_gap   = c(1,
                               median(boot_df$prop_chr, na.rm=TRUE),
                               median(boot_df$prop_sho, na.rm=TRUE))
)

dir.create("results", showWarnings = FALSE)
write.csv(summary_tbl, "results/meps_margin_decomposition.csv",
          row.names = FALSE)
saveRDS(list(summary = summary_tbl, bootstrap = boot_df),
        "results/meps_margin_bootstrap.rds")

cat("\nResults saved to results/meps_margin_decomposition.csv\n")
cat("Bootstrap saved to results/meps_margin_bootstrap.rds\n")
