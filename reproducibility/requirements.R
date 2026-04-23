################################################################################
# requirements.R - Version-pinned package installation for replication
# Snapshot date: 2026-04
# R version target: 4.5.1 (2025-06-13)
#
# Usage:
#   Rscript reproducibility/requirements.R
#
# Notes:
#   - CRAN snapshot pinning via `options(repos = ...)` with a fixed MRAN/RSPM
#     date is the most rigorous form of reproducibility; here we pin per-package
#     version via remotes::install_version() against current CRAN. Replace the
#     repo mirror as appropriate.
#   - All package versions correspond to the 6-step pipeline (scripts/01-06).
################################################################################

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

# Package => version pins (matched to 2026-04 snapshot)
pkgs <- list(
  arrow               = "23.0.1.1",
  oaxaca              = "0.1.5",
  dplyr               = "1.2.0",
  tidyr               = "1.3.2",
  survey              = "4.5",
  srvyr               = "1.3.1",
  ggplot2             = "4.0.2",
  PHEindicatormethods = "2.1.1",
  tmle                = "2.1.1",
  mice                = "3.19.0",
  episensr            = "2.1.0",
  mediation           = "4.5.1",
  SuperLearner        = "2.0-40",
  haven               = "2.5.5",
  readr               = "2.2.0",
  purrr               = "1.2.1",
  broom               = "1.0.12",
  patchwork           = "1.3.2",
  readxl              = "1.4.5",
  sandwich            = "3.1-1"
)

for (p in names(pkgs)) {
  if (!requireNamespace(p, quietly = TRUE) ||
      utils::packageVersion(p) != pkgs[[p]]) {
    message(sprintf("Installing %s == %s", p, pkgs[[p]]))
    remotes::install_version(
      p,
      version = pkgs[[p]],
      repos   = "https://cloud.r-project.org",
      upgrade = "never"
    )
  } else {
    message(sprintf("  %s %s already installed", p, pkgs[[p]]))
  }
}

cat("\nAll packages installed at pinned versions.\n")
cat("R version target: 4.5.1 (2025-06-13)\n")
cat("If your R is different, results may drift slightly.\n")
