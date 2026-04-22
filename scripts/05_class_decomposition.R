################################################################################
# 05_class_decomposition.R
#
# NHANES 5-class Oaxaca-Blinder twofold decomposition (Reimers pooled) for
# opioid, gabapentinoid, SNRI, tricyclic antidepressant (analgesic dose),
# and prescription NSAID / COX-2 inhibitor. Fitted separately for cycles
# 2005-2013 (pre-CDC 2016) and cycles 2015-2017 (post-CDC 2016), with
# inverse-variance-weighted z-test on the pre/post difference in the
# unexplained component per class.
#
# Produces numerical inputs for Exhibit 3 in the v2.0 manuscript.
#
# Requires:
#   - data/pooled/analytic.parquet       (built by 01_build_analytic.R)
#   - R packages: arrow, oaxaca, dplyr, tidyr
#
# Outputs:
#   results/class_decomposition_long.csv     per class × period long-format
#   results/class_decomposition_prepost.csv  pre/post difference with z, p
#   results/class_decomposition.rds          all oaxaca objects
################################################################################

suppressMessages({
  library(arrow)
  library(oaxaca)
  library(dplyr)
  library(tidyr)
})

set.seed(20260421)

# ---- Data ----
d <- read_parquet("data/pooled/analytic.parquet") %>%
  filter(cycle_yr >= 2005, cycle_yr <= 2017) %>%
  filter(race_f %in% c("NHBlack", "NHWhite")) %>%
  filter(age >= 20) %>%
  mutate(
    black   = as.integer(race_f == "NHBlack"),
    college = as.integer(educ_f == "Coll+"),
    period  = ifelse(cycle_yr >= 2015, "post", "pre")
  )

cat(sprintf("Analytic N = %d (pre: %d, post: %d)\n",
            nrow(d),
            sum(d$period == "pre"),
            sum(d$period == "post")))

covars   <- c("age", "pir", "college", "bmi", "insured")
outcomes <- c(
  opioid        = "opioid_any",
  gabapentinoid = "gaba_any",
  snri          = "snri_any",
  tca           = "tca_any",
  nsaid_rx      = "nsaid_rx_any"
)

run_class_oaxaca <- function(data_in, y_var, R = 300) {
  d_cc <- data_in %>%
    select(all_of(c(y_var, "black", covars))) %>%
    filter(complete.cases(.))
  d_cc$outcome_num <- as.numeric(d_cc[[y_var]])
  formula_str <- paste("outcome_num ~", paste(covars, collapse = " + "),
                       "| black")
  tryCatch(
    list(obj = oaxaca(as.formula(formula_str),
                      data = d_cc, R = R, reg.fun = lm),
         n = nrow(d_cc)),
    error = function(e) list(obj = NULL, n = 0, error = conditionMessage(e))
  )
}

results_list <- list()
summary_rows <- list()

for (class_name in names(outcomes)) {
  y_var <- outcomes[[class_name]]
  for (p in c("pre", "post")) {
    key <- paste(class_name, p, sep = "_")
    cat(sprintf(">> %s | %s (n=%d)\n", class_name, p,
                sum(d$period == p)))
    d_sub <- d %>% filter(period == p)
    res   <- run_class_oaxaca(d_sub, y_var, R = 300)
    results_list[[key]] <- res

    if (!is.null(res$obj)) {
      tf  <- res$obj$twofold$overall
      idx <- which(abs(tf[, "group.weight"] - 0.5) < 1e-6)
      if (length(idx) == 0) idx <- 1
      row <- tf[idx, , drop = TRUE]

      summary_rows[[key]] <- data.frame(
        class         = class_name,
        period        = p,
        n             = res$n,
        y_diff_pp     = as.numeric(res$obj$y[["y.diff"]]) * 100,
        explained_pp  = as.numeric(row["coef(explained)"]) * 100,
        explained_se_pp   = as.numeric(row["se(explained)"]) * 100,
        unexplained_pp    = as.numeric(row["coef(unexplained)"]) * 100,
        unexplained_se_pp = as.numeric(row["se(unexplained)"]) * 100
      )
    }
  }
}

summary_long <- bind_rows(summary_rows)

summary_wide <- summary_long %>%
  select(class, period, unexplained_pp, unexplained_se_pp) %>%
  pivot_wider(
    names_from  = period,
    values_from = c(unexplained_pp, unexplained_se_pp),
    names_glue  = "{period}_{.value}"
  ) %>%
  mutate(
    diff_post_minus_pre_pp = post_unexplained_pp - pre_unexplained_pp,
    diff_se_pp = sqrt(pre_unexplained_se_pp^2 + post_unexplained_se_pp^2),
    z          = diff_post_minus_pre_pp / diff_se_pp,
    p_value    = 2 * pnorm(-abs(z))
  )

dir.create("results", showWarnings = FALSE)
write.csv(summary_long, "results/class_decomposition_long.csv",
          row.names = FALSE)
write.csv(summary_wide, "results/class_decomposition_prepost.csv",
          row.names = FALSE)
saveRDS(results_list, "results/class_decomposition.rds")

cat("\nResults saved to results/class_decomposition_{long,prepost}.csv\n")
cat("Oaxaca objects saved to results/class_decomposition.rds\n")
