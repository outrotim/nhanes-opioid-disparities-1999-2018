# ============================================================
# 01_build_analytic.R — NHANES XPT -> pooled analytic parquet
# Repository: nhanes-opioid-disparities-1999-2018
# Build pooled analytic dataset (cycles A-J, 1999-2018).
# K cycle excluded: no drug-level data in RXQ_RX_L (only summary)
# P cycle excluded: overlaps with J; avoids double-count
# ============================================================
source(file.path(Sys.getenv("REPO_ROOT", unset=here::here()), "environment.R"))
suppressPackageStartupMessages({library(arrow); library(dplyr); library(tidyr); library(stringr)})

PQ <- file.path(DIR_HARMONIZED, "raw_parquet")

CYC <- tibble::tribble(
  ~code, ~suffix, ~year,
  "A","",       1999,
  "B","_B",     2001,
  "C","_C",     2003,
  "D","_D",     2005,
  "E","_E",     2007,
  "F","_F",     2009,
  "G","_G",     2011,
  "H","_H",     2013,
  "I","_I",     2015,
  "J","_J",     2017
)

read_pq <- function(mod, suf) {
  f <- file.path(PQ, paste0(mod, suf, ".parquet"))
  if (!file.exists(f)) return(NULL)
  arrow::read_parquet(f)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

# ---------- Opioid / adjuvant drug classification ----------
OPI_REGEX <- paste(c(
  "OXYCODONE","HYDROCODONE","MORPHINE","FENTANYL","CODEINE","OXYMORPHONE",
  "HYDROMORPHONE","TAPENTADOL","TRAMADOL","METHADONE","BUPRENORPHINE",
  "MEPERIDINE","LEVORPHANOL","PROPOXYPHENE","PENTAZOCINE","NALBUPHINE",
  "BUTORPHANOL","ALFENTANIL","SUFENTANIL","REMIFENTANIL","OPIUM"
), collapse="|")
MAT_REGEX <- "NALOXONE"
GABA_REGEX <- "GABAPENTIN|PREGABALIN"
SNRI_REGEX <- "DULOXETINE|VENLAFAXINE|MILNACIPRAN|DESVENLAFAXINE|LEVOMILNACIPRAN"
TCA_REGEX  <- "AMITRIPTYLINE|NORTRIPTYLINE|DESIPRAMINE|IMIPRAMINE|DOXEPIN|CLOMIPRAMINE|PROTRIPTYLINE|TRIMIPRAMINE"
NSAID_RX_REGEX <- "CELECOXIB|DICLOFENAC|MELOXICAM|KETOROLAC|NABUMETONE|ETODOLAC|OXAPROZIN|INDOMETHACIN|PIROXICAM|SULINDAC|TOLMETIN|FENOPROFEN|FLURBIPROFEN"

classify_drug <- function(nm) {
  nm <- toupper(nm)
  is_mat  <- !is.na(nm) & grepl(MAT_REGEX, nm)
  is_opi  <- !is.na(nm) & grepl(OPI_REGEX, nm) & !is_mat
  is_gaba <- !is.na(nm) & grepl(GABA_REGEX, nm)
  is_snri <- !is.na(nm) & grepl(SNRI_REGEX, nm)
  is_tca  <- !is.na(nm) & grepl(TCA_REGEX, nm)
  is_nsaid<- !is.na(nm) & grepl(NSAID_RX_REGEX, nm)
  tibble(opioid=is_opi, mat=is_mat, gaba=is_gaba,
         snri=is_snri, tca=is_tca, nsaid=is_nsaid)
}

MME_FACTOR <- c(
  "OXYCODONE"=1.5, "HYDROCODONE"=1.0, "MORPHINE"=1.0, "FENTANYL"=75,
  "CODEINE"=0.15, "OXYMORPHONE"=3, "HYDROMORPHONE"=4,
  "TAPENTADOL"=0.4, "TRAMADOL"=0.2, "METHADONE"=4.7,
  "BUPRENORPHINE"=10, "MEPERIDINE"=0.1
)
get_mme <- function(nm) {
  nm <- toupper(nm); f <- 0
  for (k in names(MME_FACTOR)) if (grepl(k, nm) && MME_FACTOR[k] > f) f <- MME_FACTOR[k]
  f
}

# ============================================================
# Build one cycle
# ============================================================
build_one <- function(code, suf, yr) {
  message(sprintf("[build] cycle %s (%d)", code, yr))
  demo <- read_pq("DEMO", suf); if (is.null(demo)) return(NULL)
  race <- if ("RIDRETH3" %in% names(demo)) demo$RIDRETH3 else demo$RIDRETH1
  d <- tibble(
    SEQN=demo$SEQN, cycle=code, cycle_yr=yr,
    age=safe_num(demo$RIDAGEYR), sex=safe_num(demo$RIAGENDR),
    race5=safe_num(race), pir=safe_num(demo$INDFMPIR),
    educ=safe_num(demo$DMDEDUC2), pregn=safe_num(demo$RIDEXPRG),
    exam_mo=safe_num(demo$RIDEXMON), wt_mec=safe_num(demo$WTMEC2YR),
    psu=safe_num(demo$SDMVPSU), strata=safe_num(demo$SDMVSTRA)
  )
  d <- d |> filter(age >= 20, !is.na(wt_mec), wt_mec > 0)

  # Prescription module
  rx <- read_pq("RXQ_RX", suf)
  if (!is.null(rx)) {
    drug_col <- if ("RXDDRUG" %in% names(rx)) "RXDDRUG" else "RXD240B"
    days_col <- if ("RXDDAYS" %in% names(rx)) "RXDDAYS" else "RXD295"
    rx2 <- rx |> mutate(drug = toupper(.data[[drug_col]]),
                        days = safe_num(.data[[days_col]]))
    rx2 <- bind_cols(rx2, classify_drug(rx2$drug))
    rx2 <- rx2 |> mutate(mme = if_else(opioid, vapply(drug, get_mme, numeric(1)), 0))
    rx_agg <- rx2 |>
      group_by(SEQN) |>
      summarise(
        any_rx=n()>0, opioid_any=any(opioid, na.rm=TRUE),
        opioid_days=sum(if_else(opioid & !is.na(days), days, 0), na.rm=TRUE),
        opioid_mme_d=sum(if_else(opioid, mme, 0), na.rm=TRUE),
        gaba_any=any(gaba, na.rm=TRUE), snri_any=any(snri, na.rm=TRUE),
        tca_any=any(tca, na.rm=TRUE), nsaid_rx_any=any(nsaid, na.rm=TRUE),
        n_opioid_rx=sum(opioid, na.rm=TRUE), .groups="drop"
      )
    d <- d |> left_join(rx_agg, by="SEQN") |>
      mutate(across(c(opioid_any, gaba_any, snri_any, tca_any, nsaid_rx_any),
                    ~replace_na(., FALSE)),
             across(c(opioid_days, opioid_mme_d, n_opioid_rx), ~replace_na(., 0)))
  }

  # Anthropometry + BP
  bmx <- read_pq("BMX", suf)
  if (!is.null(bmx)) d <- d |> left_join(bmx |> dplyr::select(SEQN, waist=BMXWAIST, bmi=BMXBMI), by="SEQN")
  bpx <- read_pq("BPX", suf)
  if (!is.null(bpx)) {
    sy <- grep("^BPXSY[1-4]$", names(bpx), value=TRUE)
    di <- grep("^BPXDI[1-4]$", names(bpx), value=TRUE)
    if (length(sy)>0) bpx$sbp <- rowMeans(as.matrix(bpx[, sy]), na.rm=TRUE)
    if (length(di)>0) bpx$dbp <- rowMeans(as.matrix(bpx[, di]), na.rm=TRUE)
    d <- d |> left_join(bpx |> dplyr::select(SEQN, sbp, dbp), by="SEQN")
  }

  # Labs: lipids + glucose (A-C legacy names, D+ standard)
  t <- read_pq("TCHOL", suf); tg <- read_pq("TRIGLY", suf); h  <- read_pq("HDL", suf)
  if (!is.null(t)) d <- d |> left_join(tibble(SEQN=t$SEQN, tchol=safe_num(t$LBXTC)), by="SEQN")
  if (!is.null(tg)) d <- d |> left_join(tibble(SEQN=tg$SEQN, trig=safe_num(tg$LBXTR)), by="SEQN")
  if (!is.null(h)) d <- d |> left_join(tibble(SEQN=h$SEQN, hdl=safe_num(h$LBDHDD)), by="SEQN")
  g <- read_pq("GLU", suf); gh <- read_pq("GHB", suf)
  if (!is.null(g)) d <- d |> left_join(tibble(SEQN=g$SEQN, fbg=safe_num(g$LBXGLU)), by="SEQN")
  if (!is.null(gh)) d <- d |> left_join(tibble(SEQN=gh$SEQN, hba1c=safe_num(gh$LBXGH)), by="SEQN")

  # Questionnaire modules
  smq <- read_pq("SMQ", suf)
  if (!is.null(smq) && "SMQ020" %in% names(smq)) {
    d <- d |> left_join(
      tibble(SEQN=smq$SEQN, smq020=safe_num(smq$SMQ020),
             smq040=if("SMQ040" %in% names(smq)) safe_num(smq$SMQ040) else NA_real_), by="SEQN")
  }
  alq <- read_pq("ALQ", suf)
  if (!is.null(alq)) {
    alq_col <- if ("ALQ130" %in% names(alq)) "ALQ130" else if ("ALQ120Q" %in% names(alq)) "ALQ120Q" else NA
    if (!is.na(alq_col)) d <- d |> left_join(tibble(SEQN=alq$SEQN, drinks_day=safe_num(alq[[alq_col]])), by="SEQN")
  }
  mcq <- read_pq("MCQ", suf)
  if (!is.null(mcq)) {
    d <- d |> left_join(tibble(
      SEQN=mcq$SEQN,
      had_cancer = if("MCQ220" %in% names(mcq)) as.integer(safe_num(mcq$MCQ220)==1) else NA_integer_,
      had_chd    = if("MCQ160C" %in% names(mcq)) as.integer(safe_num(mcq$MCQ160C)==1) else NA_integer_,
      had_stroke = if("MCQ160F" %in% names(mcq)) as.integer(safe_num(mcq$MCQ160F)==1) else NA_integer_,
      had_arthritis = if("MCQ160A" %in% names(mcq)) as.integer(safe_num(mcq$MCQ160A)==1) else NA_integer_
    ), by="SEQN")
  }
  diq <- read_pq("DIQ", suf)
  if (!is.null(diq) && "DIQ010" %in% names(diq)) {
    d <- d |> left_join(tibble(SEQN=diq$SEQN,
      had_dm = as.integer(safe_num(diq$DIQ010)==1),
      dm_rx  = if("DID070" %in% names(diq)) as.integer(safe_num(diq$DID070)==1) else NA_integer_), by="SEQN")
  }
  bpq <- read_pq("BPQ", suf)
  if (!is.null(bpq) && "BPQ020" %in% names(bpq)) {
    d <- d |> left_join(tibble(SEQN=bpq$SEQN,
      had_htn  = as.integer(safe_num(bpq$BPQ020)==1),
      htn_rx   = if("BPQ050A" %in% names(bpq)) as.integer(safe_num(bpq$BPQ050A)==1) else NA_integer_,
      lipid_rx = if("BPQ100D" %in% names(bpq)) as.integer(safe_num(bpq$BPQ100D)==1) else NA_integer_), by="SEQN")
  }
  hiq <- read_pq("HIQ", suf)
  if (!is.null(hiq) && "HIQ011" %in% names(hiq)) {
    d <- d |> left_join(tibble(SEQN=hiq$SEQN, insured=as.integer(safe_num(hiq$HIQ011)==1)), by="SEQN")
  }
  dpq <- read_pq("DPQ", suf)
  if (!is.null(dpq)) {
    q <- intersect(paste0("DPQ0", sprintf("%02d", 10:90)), names(dpq))
    if (length(q)>=8) {
      mat <- as.matrix(dpq[, q]); mat[mat >= 7] <- NA
      phq <- rowSums(mat, na.rm=TRUE); phq[rowSums(is.na(mat))>2] <- NA
      d <- d |> left_join(tibble(SEQN=dpq$SEQN, phq9=phq), by="SEQN")
    }
  }
  d
}

# Build + pool
pool <- purrr::map_dfr(seq_len(nrow(CYC)), \(i) build_one(CYC$code[i], CYC$suffix[i], CYC$year[i]))
message(sprintf("[build] pooled rows: %d", nrow(pool)))

# Derived variables
pool <- pool |> mutate(
  n_cycles = nrow(CYC),
  wt_comb  = wt_mec / n_cycles,
  sex_f    = factor(sex, levels=c(1,2), labels=c("M","F")),
  race_f   = case_when(
    race5 %in% c(1,2) ~ "Hispanic", race5 == 3 ~ "NHWhite",
    race5 == 4 ~ "NHBlack", race5 %in% c(5,6,7) ~ "Other", TRUE ~ NA_character_),
  educ_f = case_when(
    educ %in% c(1,2) ~ "<HS", educ == 3 ~ "HS",
    educ == 4 ~ "SomeColl", educ == 5 ~ "Coll+", TRUE ~ NA_character_),
  pir_q = ntile(pir, 5),
  smoker = case_when(
    smq020 == 2 ~ "Never",
    smq020 == 1 & smq040 %in% c(1,2) ~ "Current",
    smq020 == 1 & smq040 == 3 ~ "Former", TRUE ~ NA_character_),
  bmi_c = case_when(is.na(bmi) ~ NA_character_, bmi < 25 ~ "Normal",
                    bmi < 30 ~ "Over", TRUE ~ "Obese"),
  mets_waist = as.integer((sex==1 & waist>=102) | (sex==2 & waist>=88)),
  mets_tg    = as.integer(trig >= 150 | replace_na(lipid_rx,0L)==1),
  mets_hdl   = as.integer((sex==1 & hdl<40) | (sex==2 & hdl<50)),
  mets_bp    = as.integer(sbp>=130 | dbp>=85 | replace_na(htn_rx,0L)==1),
  mets_fbg   = as.integer(fbg>=100 | replace_na(dm_rx,0L)==1),
  mets_count = rowSums(cbind(mets_waist, mets_tg, mets_hdl, mets_bp, mets_fbg), na.rm=TRUE),
  mets_na    = is.na(waist)+is.na(trig)+is.na(hdl)+is.na(sbp)+is.na(fbg),
  mets       = as.integer(mets_count >= 3),
  mets       = ifelse(mets_na >= 2, NA_integer_, mets),
  opioid_cat = case_when(
    !opioid_any ~ "None",
    n_opioid_rx == 1 & opioid_days < 30 ~ "Short",
    TRUE ~ "Chronic"),
  post_CDC2016 = as.integer(cycle_yr >= 2015)
)

# Save
out_file <- file.path(DIR_POOLED, "analytic.parquet")
arrow::write_parquet(pool, out_file, compression="zstd")
message(sprintf("[build] analytic dataset saved: %s (%d rows, %d cols)",
                out_file, nrow(pool), ncol(pool)))
