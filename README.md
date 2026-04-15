# nhanes-opioid-disparities-1999-2018

Replication code for **"Two decades of racial disparities in opioid prescribing in the United States: a decomposition analysis of the National Health and Nutrition Examination Survey, 1999–2018"** (submitted to *The Lancet Regional Health – Americas*, 2026).

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![Docs CC BY 4.0](https://img.shields.io/badge/Docs-CC--BY%204.0-lightgrey.svg)](LICENSE)

## Contents

```
README.md
LICENSE                   MIT (code) + CC BY 4.0 (documentation)
environment.R             packages and path configuration (R 4.4+)
scripts/
  01_build_analytic.R     NHANES XPT  ->  pooled analytic parquet
  02_main_analyses.R      Q1 Oaxaca-Blinder + Q2 CDC-2016 ITS + Q4 trends
  03_sensitivity.R        MICE, E-value, NCO, MEPS, tramadol 2014
  04_figures.R            Figures 1-5 and eFigures
```

No individual-level data are included.

## Data availability

| Layer | Source | Access | In this repo |
|---|---|---|---|
| NHANES XPT files (1999–2018, 10 cycles) | CDC / NCHS | https://wwwn.cdc.gov/nchs/nhanes/ | No (download directly) |
| MEPS 2016 HC-188A Prescribed Medicines | AHRQ | https://meps.ahrq.gov/ | No (download directly) |
| Multum Lexicon drug-class mapping | NHANES RXQ_DRUG | bundled with XPT | No |
| Derived analytic parquet (MEC weights, PSU, strata) | built by `01_build_analytic.R` | rebuild locally | No (see note below) |
| Intermediate `.rds` / `.csv` results | built by `02`, `03` | rebuild locally | No |

### Why the derived analytic file is not redistributed

Every input variable comes from a public NHANES release, but the derived file joins masked variance units (SDMVPSU), strata (SDMVSTRA) and MEC-weighted indicators across ten cycles of the same respondents. Redistributing that combined object adds a small theoretical re-identification concern without adding reproducibility value, because any reader can rebuild it by running `01_build_analytic.R` against a local NHANES XPT mirror.

## Reproducing the analysis

Requirements: R ≥ 4.4, about 5 GB of free disk for the XPT files and derived parquet, any POSIX shell or PowerShell.

```bash
git clone https://github.com/outrotim/nhanes-opioid-disparities-1999-2018.git
cd nhanes-opioid-disparities-1999-2018

# Download the ten NHANES biennial cycles (1999-2018) from CDC and place
# the XPT files under a local work directory, e.g. ~/nhanes_work/data/raw/.
# Modules required: DEMO, RXQ_RX, RXQ_DRUG, DPQ, BMX, BPX, TCHOL, TRIGLY,
# HDL, GLU, GHB, SMQ, ALQ, DIQ, BPQ, MCQ, HIQ, OCQ.
# CYCLES and MODULES tables are in environment.R.

export STUDY_ROOT=~/nhanes_work

Rscript -e 'source("environment.R")'
Rscript scripts/01_build_analytic.R
Rscript scripts/02_main_analyses.R
Rscript scripts/03_sensitivity.R
Rscript scripts/04_figures.R
```

End-to-end runtime is about 15–25 minutes on a 2020 M1 MacBook Pro; the MICE step dominates.

## Analytic choices worth knowing before reuse

Opioid identification uses Multum generic drug-class codes 60, 191, 57 and 58 from RXQ_DRUG. Tramadol is included throughout; its 2014 DEA rescheduling is tested as a stratification in `03_sensitivity.R`.

Race-ethnicity is non-Hispanic White versus non-Hispanic Black (NHANES `RIDRETH1` collapsed). Hispanic and Asian respondents are not the primary contrast and appear only in descriptive tables.

The decomposition is Oaxaca–Blinder twofold with Reimers pooled weights, fitted cycle by cycle and pooled by inverse-variance-weighted fixed-effects meta-analysis. Implementation: `oaxaca` v0.1.4 on top of `survey` v4.2. The survey design is `svydesign(id = ~psu, strata = ~strata, weights = ~wt_comb, nest = TRUE)` with `survey.lonely.psu = "adjust"`. MEC weights are divided by the number of pooled cycles per NCHS guidance; see `compute_combined_weight()` in `environment.R`.

## Caveats

Results apply to non-institutionalised US civilian adults; they do not transfer to inpatient, correctional or pediatric settings without recalibration.

The 2019–March 2020 pre-pandemic subsample (`P`) and the August 2021–August 2023 release (`L`) use non-standard weighting structures and are excluded from the primary 1999–2018 window. Extending the series requires reconstructing combined weights following current NCHS guidance.

Opioid use is captured by pill-bottle review at the MEC examination over a past-30-day reference window. This is an upper bound on chronic therapy and does not capture illicit use.

## License

Code (`.R` files, `environment.R`) under the [MIT License](LICENSE); documentation under CC BY 4.0. Please cite the manuscript (§ Citation) when reusing either.

## Citation

To be updated after DOI assignment:

> [USER: INSERT Authors]. Two decades of racial disparities in opioid prescribing in the United States: a decomposition analysis of the National Health and Nutrition Examination Survey, 1999–2018. *The Lancet Regional Health – Americas* 2026 (forthcoming). DOI: [to be assigned on acceptance].

For the exact code version used in the manuscript, cite the Zenodo release:

> [USER: INSERT after Zenodo release] — Zenodo DOI: 10.5281/zenodo.XXXXXXX

## Contact

Please open a GitHub issue for code-related questions. For data interpretation, the NCHS and AHRQ public helpdesks are the authoritative sources.
