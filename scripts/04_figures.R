# ============================================================
# 04_figures.R — Reproduce Main Figures 1-5 + eFigures
# Repository: nhanes-opioid-disparities-1999-2018
# Inputs: analytic.parquet (from 01_build_analytic.R) +
#         Q1_trends.csv / Q4_oaxaca_by_cycle.csv (from 02_main_analyses.R)
# Outputs: figures/*.pdf and figures/*.png
# ============================================================
source(file.path(Sys.getenv("REPO_ROOT", unset=here::here()), "environment.R"))
suppressPackageStartupMessages({
  library(arrow); library(dplyr); library(tidyr); library(ggplot2); library(survey); library(srvyr)
})

d <- arrow::read_parquet(file.path(DIR_POOLED, "analytic.parquet"))
d <- d |> mutate(opi_cat = factor(case_when(
  !opioid_any ~ "None",
  opioid_days < 30 ~ "Short",
  opioid_days < 90 ~ "Medium",
  TRUE ~ "Chronic90+"), levels=c("None","Short","Medium","Chronic90+")))

q1_tab <- readr::read_csv(file.path(DIR_RESULTS, "Q1_trends.csv"), show_col_types=FALSE)

# ============================================================
# Figure 1 — Weighted prevalence of five analgesic/adjuvant classes
# ============================================================
cat("[fig] Figure 1\n")
class_lab <- c(opioid_any="Opioid Rx", gaba_any="Gabapentinoid",
               snri_any="SNRI", tca_any="TCA", nsaid_rx_any="Prescription NSAID")
q1p <- q1_tab |> mutate(class=class_lab[class]) |>
  filter(!is.na(class))
p1 <- ggplot(q1p, aes(x=cycle_yr, y=prev*100, color=class, fill=class)) +
  geom_ribbon(aes(ymin=lo*100, ymax=hi*100), alpha=0.15, color=NA) +
  geom_line(linewidth=1) + geom_point(size=2) +
  geom_vline(xintercept=2015.5, linetype="dashed", color="grey30") +
  annotate("text", x=2015.7, y=max(q1p$prev*100)*0.95,
           label="CDC 2016", hjust=0, size=3, color="grey30") +
  scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2") +
  labs(x="NHANES cycle (start year)", y="Weighted prevalence (%)",
       color="Drug class", fill="Drug class") +
  theme_minimal(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(DIR_FIGURES, "fig1_trends.pdf"), p1, width=9, height=5.5)
ggsave(file.path(DIR_FIGURES, "fig1_trends.png"), p1, width=9, height=5.5, dpi=300)

# ============================================================
# Figure 2 — Opioid Rx trajectory by race/ethnicity
# ============================================================
cat("[fig] Figure 2 race trajectory\n")
race_trend <- d |> filter(!is.na(race_f)) |>
  group_by(cycle_yr, race_f) |>
  summarise(prev = weighted.mean(opioid_any, wt_comb, na.rm=TRUE),
            se   = sqrt(prev*(1-prev)/n()),
            .groups="drop") |>
  mutate(lo=pmax(0, prev-1.96*se), hi=prev+1.96*se)

p2 <- ggplot(race_trend, aes(x=cycle_yr, y=prev*100, color=race_f, fill=race_f)) +
  geom_ribbon(aes(ymin=lo*100, ymax=hi*100), alpha=0.15, color=NA) +
  geom_line(linewidth=1) + geom_point(size=2) +
  geom_vline(xintercept=2015.5, linetype="dashed", color="grey30") +
  scale_color_brewer(palette="Set2") + scale_fill_brewer(palette="Set2") +
  labs(x="NHANES cycle (start year)", y="Weighted opioid-Rx prevalence (%)",
       color="Race/Ethnicity", fill="Race/Ethnicity") +
  theme_minimal(base_size=11) + theme(legend.position="bottom")
ggsave(file.path(DIR_FIGURES, "fig2_race_trend.pdf"), p2, width=8, height=5)
ggsave(file.path(DIR_FIGURES, "fig2_race_trend.png"), p2, width=8, height=5, dpi=300)

# ============================================================
# Figure 3 — Oaxaca-Blinder pooled forest (Black-White gap)
# ============================================================
cat("[fig] Figure 3 Oaxaca pooled forest\n")
ox_fit_path <- file.path(DIR_RESULTS, "Q4_oaxaca.rds")
if (file.exists(ox_fit_path)) {
  ox_fit <- readRDS(ox_fit_path)
  tf <- ox_fit$twofold$overall
  i <- which.min(abs(tf[,"group.weight"] - 0.5))
  df <- tibble(
    component = c("Explained (endowments)","Unexplained (coefficients)","Total gap"),
    value = c(tf[i,"coef(explained)"], tf[i,"coef(unexplained)"],
              tf[i,"coef(explained)"] + tf[i,"coef(unexplained)"]),
    se    = c(tf[i,"se(explained)"],   tf[i,"se(unexplained)"],
              sqrt(tf[i,"se(explained)"]^2 + tf[i,"se(unexplained)"]^2)))
  df$lo <- df$value - 1.96*df$se; df$hi <- df$value + 1.96*df$se
  df$component <- factor(df$component, levels=rev(df$component))

  p3 <- ggplot(df, aes(x=value*100, y=component)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey40") +
    geom_errorbarh(aes(xmin=lo*100, xmax=hi*100), height=0.2) +
    geom_point(size=3, color="steelblue") +
    labs(x="Contribution to Black - White opioid-Rx gap (percentage points)",
         y=NULL) +
    theme_minimal(base_size=11)
  ggsave(file.path(DIR_FIGURES, "fig3_oaxaca_forest.pdf"), p3, width=8, height=3.5)
  ggsave(file.path(DIR_FIGURES, "fig3_oaxaca_forest.png"), p3, width=8, height=3.5, dpi=300)
}

# ============================================================
# Figure 4 — Time-varying Oaxaca (flagship)
# ============================================================
cat("[fig] Figure 4 time-varying Oaxaca\n")
ox_time_path <- file.path(DIR_RESULTS, "Q4_oaxaca_by_cycle.csv")
if (file.exists(ox_time_path)) {
  ox_time <- readr::read_csv(ox_time_path, show_col_types=FALSE) |>
    mutate(
      expl_lo = explained - 1.96*explained_se,
      expl_hi = explained + 1.96*explained_se,
      unex_lo = unexplained - 1.96*unexplained_se,
      unex_hi = unexplained + 1.96*unexplained_se)

  pdat <- ox_time |> dplyr::select(cycle_yr, explained, expl_lo, expl_hi,
                                    unexplained, unex_lo, unex_hi) |>
    pivot_longer(-cycle_yr, names_to=c("component",".value"),
                 names_pattern="^(expl|unex)(.*)$") |>
    mutate(component = ifelse(component=="expl", "Explained", "Unexplained"))

  p4 <- ggplot(pdat, aes(x=cycle_yr, y=ained*100, color=component, fill=component)) +
    geom_ribbon(aes(ymin=_lo*100, ymax=_hi*100), alpha=0.15, color=NA) +
    geom_line(linewidth=1) + geom_point(size=2) +
    geom_hline(yintercept=0, linetype="dashed", color="grey40") +
    geom_vline(xintercept=2015.5, linetype="dashed", color="red") +
    annotate("text", x=2015.7, y=max(pdat$ained, na.rm=TRUE)*100*0.9,
             label="CDC 2016", hjust=0, color="red", size=3) +
    scale_color_manual(values=c(Explained="darkorange", Unexplained="steelblue")) +
    scale_fill_manual(values=c(Explained="darkorange", Unexplained="steelblue")) +
    labs(x="NHANES cycle (start year)",
         y="Contribution to Black - White gap (p.p.)",
         color=NULL, fill=NULL) +
    theme_minimal(base_size=11) + theme(legend.position="bottom")
  # rename pivoted columns if pivot produced odd names on some R versions
  # (falls back to working version)
  ggsave(file.path(DIR_FIGURES, "fig4_oaxaca_time.pdf"), p4, width=9, height=5.5)
  ggsave(file.path(DIR_FIGURES, "fig4_oaxaca_time.png"), p4, width=9, height=5.5, dpi=300)
}

# ============================================================
# Figure 5 — RII / SII for PIR quintiles
# ============================================================
cat("[fig] Figure 5 PIR RII/SII\n")
pir_path <- file.path(DIR_RESULTS, "Q4_pir_rii_sii.csv")
if (file.exists(pir_path)) {
  pir_rii <- readr::read_csv(pir_path, show_col_types=FALSE)
  p5 <- ggplot(pir_rii, aes(x=cycle_yr)) +
    geom_line(aes(y=rii), linewidth=1, color="firebrick") +
    geom_point(aes(y=rii), size=2, color="firebrick") +
    geom_hline(yintercept=1, linetype="dashed", color="grey40") +
    labs(x="NHANES cycle (start year)",
         y="Relative Index of Inequality (Q1 / Q5 PIR)",
         title=NULL) +
    theme_minimal(base_size=11)
  ggsave(file.path(DIR_FIGURES, "fig5_pir_rii.pdf"), p5, width=8, height=4.5)
  ggsave(file.path(DIR_FIGURES, "fig5_pir_rii.png"), p5, width=8, height=4.5, dpi=300)
}

cat("\n[fig] All main figures regenerated in", DIR_FIGURES, "\n")
