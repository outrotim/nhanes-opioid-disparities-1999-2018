# nhanes-opioid-disparities-1999-2018

> Replication code for:
> **Two decades of racial disparities in opioid prescribing in the United States: a decomposition analysis of the National Health and Nutrition Examination Survey, 1999–2018**
> *Submitted to* The Lancet Regional Health – Americas, 2026.

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![Docs CC BY 4.0](https://img.shields.io/badge/Docs-CC--BY%204.0-lightgrey.svg)](LICENSE)

---

## 1. What is in this repository

This repository contains the **minimum necessary analytic code** to reproduce the main findings, decomposition results, and key sensitivity analyses of the manuscript. It does **not** contain any individual-level data.

```
.
├── README.md                  # This file
├── LICENSE                    # MIT (code) + CC-BY 4.0 (documentation)
├── environment.R              # Package list (R 4.4+) and path config
└── scripts/
    ├── 01_build_analytic.R    # NHANES XPT  ->  pooled analytic parquet
    ├── 02_main_analyses.R     # Q1 Oaxaca-Blinder + Q2 CDC-2016 ITS + Q4 trends
    ├── 03_sensitivity.R       # MICE / E-value / NCO / MEPS / tramadol 2014
    └── 04_figures.R           # Reproduce Main Figures 1-5 + eFigures
```

---

## 2. Data availability (tiered)

| Layer | Source | Access | In this repo? |
|---|---|---|---|
| **NHANES XPT raw files** (1999–2018, 10 cycles) | CDC / NCHS | Public — https://wwwn.cdc.gov/nchs/nhanes/ | No (download yourself, it is public) |
| **MEPS 2016 HC-188A** Prescribed Medicines | AHRQ | Public — https://meps.ahrq.gov/ | No (download yourself, it is public) |
| **Multum Lexicon drug-class mapping** | NHANES RXQ_DRUG | Public (bundled with XPT) | No |
| **Derived analytic parquet** (MEC weights, PSU, strata) | Built by `01_build_analytic.R` | Reconstructible from public XPT | No *(see §3)* |
| **Intermediate result objects** (`.rds` / `.csv`) | Built by `02_main_analyses.R`, `03_sensitivity.R` | Reconstructible from derived parquet | No *(re-run scripts)* |

### Why we do not ship the derived analytic dataset

Although every input variable is from a fully public NHANES release, the derived file combines NHANES masked variance units (SDMVPSU), strata (SDMVSTRA), and MEC-weighted indicators across ten cycles for the same respondents. Redistribution of this combined object could marginally raise re-identification theoretical risk and is redundant — any reader can rebuild it by running `01_build_analytic.R` against a local NHANES XPT mirror.

---

## 3. How to reproduce

### Prerequisites
- R ≥ 4.4
- ~5 GB free disk for NHANES XPT + derived parquet
- Unix-like shell or PowerShell

### Steps

```bash
# 1. Clone the repo
git clone https://github.com/outrotim/nhanes-opioid-disparities-1999-2018.git
cd nhanes-opioid-disparities-1999-2018

# 2. Download the ten NHANES biennial cycles (1999-2018) from CDC.
#    Place XPT files under a directory of your choice, e.g. ~/nhanes_work/data/raw/
#    Module list: DEMO, RXQ_RX, RXQ_DRUG, DPQ, BMX, BPX, TCHOL, TRIGLY, HDL,
#                  GLU, GHB, SMQ, ALQ, DIQ, BPQ, MCQ, HIQ, OCQ
#    (See CYCLES and MODULES tables in environment.R.)

# 3. Point the code to your local work directory:
export STUDY_ROOT=~/nhanes_work

# 4. Launch R and source scripts in order:
Rscript -e 'source("environment.R")'
Rscript scripts/01_build_analytic.R   # builds pooled analytic.parquet
Rscript scripts/02_main_analyses.R    # runs Q1 / Q2 / Q4
Rscript scripts/03_sensitivity.R      # runs 5 key sensitivity analyses
Rscript scripts/04_figures.R          # regenerates Fig 1-5
```

Expected runtime: ~15–25 minutes on a 2020 M1 MacBook Pro (primarily MICE m = 20).

---

## 4. Key analytic choices (please read before reusing)

- **Opioid identification** — Multum generic drug-class codes 60 (opioid analgesic combinations), 191 (narcotic analgesics), 57 (cough and cold), and 58 (analgesic combinations). Tramadol is included throughout; its 2014-08 DEA rescheduling is examined as a sensitivity stratification in `03_sensitivity.R`.
- **Race-ethnicity** — non-Hispanic White vs. non-Hispanic Black per NHANES `RIDRETH1` collapsed. Other groups (Hispanic, Asian) were not the primary contrast and are reported descriptively only.
- **Decomposition method** — Oaxaca-Blinder twofold, Reimers pooled weights, time-varying (cycle-specific), pooled via inverse-variance-weighted fixed-effects meta-analysis. Implementation: `oaxaca` R package v0.1.4.
- **Survey design** — `survey::svydesign(id = ~psu, strata = ~strata, weights = ~wt_comb, nest = TRUE)`; `options(survey.lonely.psu = "adjust")`.
- **Combined MEC weights** — divided by the number of cycles pooled (NCHS guidance); see `compute_combined_weight()` in `environment.R`.

---

## 5. Caveats

- These results apply to **non-institutionalised U.S. civilian adults**. They do not recalibrate automatically to inpatient, correctional, or pediatric populations.
- The 2019-March 2020 (`P`) pre-pandemic subsample and the August 2021-August 2023 (`L`) post-pandemic release use non-standard weighting structures and were **excluded** from the primary analysis window (1999–2018). Extending the series requires re-deriving combined weights per current NCHS guidance.
- Opioid use is captured via **pill-bottle review at the MEC examination** (past-30-day reference window). This is an upper-bound proxy for chronic opioid therapy and under-captures illicit use.

---

## 6. License

- **Code** (`.R` files, `environment.R`) — [MIT License](LICENSE).
- **Documentation** (`README.md`) — CC BY 4.0.

In either case: please cite the manuscript (§7) if you reuse the code.

---

## 7. Citation

*(Please update after DOI assignment.)*

> [USER: INSERT Authors]. Two decades of racial disparities in opioid prescribing in the United States: a decomposition analysis of the National Health and Nutrition Examination Survey, 1999–2018. *The Lancet Regional Health – Americas* 2026 *(forthcoming)*. DOI: *[to be assigned on acceptance]*.

For the exact code version used in the manuscript, please cite the Zenodo archived release:

> [USER: INSERT after Zenodo release] — Zenodo DOI: 10.5281/zenodo.XXXXXXX

---

## 8. Contact

Issues and questions: please open a GitHub issue on this repository. For data-related correspondence (MEPS / NHANES interpretation), the NCHS and AHRQ public helpdesks are the authoritative sources.

---

*Repository built as the minimum-necessary reproducibility artifact for peer review. It deliberately excludes derived datasets and intermediate results to respect NHANES redistribution norms.*
