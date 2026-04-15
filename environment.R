# ============================================================
# environment.R — Package setup + path configuration
# Replication repository: nhanes-opioid-disparities-1999-2018
# ============================================================
# Usage:
#   1. Set the environment variable STUDY_ROOT to the directory where
#      you have placed NHANES XPT files + derived parquet outputs, e.g.
#
#        Sys.setenv(STUDY_ROOT = "/path/to/your/local/workdir")
#
#      Or create an .Renviron file in the repo root with:
#        STUDY_ROOT=/path/to/your/local/workdir
#
#   2. source("environment.R") before running any analysis script.
#
# Tested with R 4.4.1 on macOS 14.
# ============================================================

suppressPackageStartupMessages({
  needed <- c(
    "haven", "dplyr", "tidyr", "readr", "purrr", "stringr",
    "survey", "srvyr", "rms", "mice", "oaxaca", "segmented",
    "PHEindicatormethods", "sandwich", "lmtest", "broom",
    "arrow", "fs", "glue", "cli", "here",
    "ggplot2", "EValue"
  )
  missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(needed, library, character.only = TRUE))
})

# ---------- Resolve study root (no hard-coded author paths) ----------
STUDY_ROOT <- Sys.getenv("STUDY_ROOT", unset = NA)
if (is.na(STUDY_ROOT) || !nzchar(STUDY_ROOT)) {
  # Fall back to repo root (works when sourced from the cloned repository)
  STUDY_ROOT <- tryCatch(here::here(), error = function(e) getwd())
  message("[env] STUDY_ROOT not set; defaulting to: ", STUDY_ROOT)
}
stopifnot(dir.exists(STUDY_ROOT))

DIR_RAW        <- file.path(STUDY_ROOT, "data", "raw")
DIR_HARMONIZED <- file.path(STUDY_ROOT, "data", "harmonized")
DIR_POOLED     <- file.path(STUDY_ROOT, "data", "pooled")
DIR_RESULTS    <- file.path(STUDY_ROOT, "results")
DIR_FIGURES    <- file.path(STUDY_ROOT, "figures")
DIR_TABLES     <- file.path(STUDY_ROOT, "tables")

for (d in c(DIR_RAW, DIR_HARMONIZED, DIR_POOLED,
            DIR_RESULTS, DIR_FIGURES, DIR_TABLES)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---------- NHANES cycle registry ----------
CYCLES <- tibble::tribble(
  ~code, ~years,      ~prefix, ~suffix, ~yr_mid, ~yr_start,
  "A",   "1999-2000", "",      "",      1999.5,  1999L,
  "B",   "2001-2002", "",      "_B",    2001.5,  2001L,
  "C",   "2003-2004", "",      "_C",    2003.5,  2003L,
  "D",   "2005-2006", "",      "_D",    2005.5,  2005L,
  "E",   "2007-2008", "",      "_E",    2007.5,  2007L,
  "F",   "2009-2010", "",      "_F",    2009.5,  2009L,
  "G",   "2011-2012", "",      "_G",    2011.5,  2011L,
  "H",   "2013-2014", "",      "_H",    2013.5,  2013L,
  "I",   "2015-2016", "",      "_I",    2015.5,  2015L,
  "J",   "2017-2018", "",      "_J",    2017.5,  2017L
)

MODULES <- c(
  demo = "DEMO", rx = "RXQ_RX", rx_ing = "RXQ_DRUG",
  dpq = "DPQ", bmx = "BMX", bpx = "BPX",
  tchol = "TCHOL", trigly = "TRIGLY", hdl = "HDL",
  glu = "GLU", ghb = "GHB", smq = "SMQ", alq = "ALQ",
  diq = "DIQ", bpq = "BPQ", mcq = "MCQ", hiq = "HIQ", ocq = "OCQ"
)

# ---------- Policy anchor dates ----------
CDC_OPIOID_GUIDELINE_DATE <- as.Date("2016-03-18")
TRAMADOL_SCHEDULED_DATE   <- as.Date("2014-08-18")

# ---------- Reproducibility seeds ----------
SEED_MICE      <- 20260414L
SEED_BOOTSTRAP <- 20260415L
SEED_JOINPOINT <- 20260416L

options(survey.lonely.psu = "adjust")

# ---------- Helpers ----------
compute_combined_weight <- function(df, k_cycles) {
  stopifnot("WTMEC2YR" %in% names(df), k_cycles >= 1)
  df$WTMEC_combined <- df$WTMEC2YR / k_cycles
  df
}

cdc_xpt_url <- function(module, prefix, suffix, yr_start) {
  glue::glue("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{yr_start}/DataFiles/{prefix}{module}{suffix}.XPT")
}

message("[env] environment.R sourced. STUDY_ROOT = ", STUDY_ROOT)
