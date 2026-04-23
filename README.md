# nhanes-opioid-disparities-1999-2018

Replication code for **"The Black-White Gap In U.S. Opioid Prescribing Concentrates At The Short-Course Margin: A Cross-Database NHANES-MEPS Decomposition With Implications For The 2016 CDC Guideline"** (under review at *Pharmacoepidemiology and Drug Safety*, 2026).

[![License: MIT](https://img.shields.io/badge/Code-MIT-blue.svg)](LICENSE)
[![Docs CC BY 4.0](https://img.shields.io/badge/Docs-CC--BY%204.0-lightgrey.svg)](LICENSE)

## Contents

```
README.md
LICENSE                      MIT (code) + CC BY 4.0 (documentation)
environment.R                packages and path configuration (R 4.4+)
scripts/
  01_build_analytic.R        NHANES XPT  ->  pooled analytic parquet
  02_main_analyses.R         NHANES Oaxaca-Blinder + CDC-2016 stability test
  03_sensitivity.R           MICE, E-value, NCO, MEPS ratio check, tramadol 2014
  04_figures.R               Manuscript figures and supplementary figures
  05_class_decomposition.R   NHANES 5-class Oaxaca (opioid + 4 adjuvant classes)
  06_meps_margin.R           MEPS 2015-2016 pooled extensive vs chronic vs
                             short-course margin decomposition with bootstrap
```

No individual-level data are included.

## Data availability

| Layer | Source | Access | In this repo |
|---|---|---|---|
| NHANES XPT files (1999–2018, 10 cycles) | CDC / NCHS | https://wwwn.cdc.gov/nchs/nhanes/ | No (download directly) |
| MEPS Full-Year Consolidated 2015 (HC-181), 2016 (HC-192) | AHRQ | https://meps.ahrq.gov/ | No (download directly) |
| MEPS Prescribed Medicines 2015 (HC-178A), 2016 (HC-188A) | AHRQ | https://meps.ahrq.gov/ | No (download directly) |
| Multum Lexicon drug-class mapping | NHANES RXQ_DRUG | bundled with XPT | No |
| Derived analytic parquet (MEC weights, PSU, strata) | built by `01_build_analytic.R` | rebuild locally | No (see note below) |
| Intermediate `.rds` / `.csv` results | built by `02`, `03`, `05`, `06` | rebuild locally | No |

### Why the derived analytic file is not redistributed

Every input variable comes from a public NHANES release, but the derived file joins masked variance units (SDMVPSU), strata (SDMVSTRA) and MEC-weighted indicators across ten cycles of the same respondents. Redistributing that combined object adds a small theoretical re-identification concern without adding reproducibility value, because any reader can rebuild it by running `01_build_analytic.R` against a local NHANES XPT mirror.

## Reproducing the analysis

Requirements: R ≥ 4.4, about 6 GB of free disk for the NHANES XPT files, the MEPS files, and the derived parquet, any POSIX shell or PowerShell.

```bash
git clone https://github.com/outrotim/nhanes-opioid-disparities-1999-2018.git
cd nhanes-opioid-disparities-1999-2018

# Step 1: download the ten NHANES biennial cycles (1999-2018) from CDC
# and place the XPT files under a local work directory, e.g.
#   $STUDY_ROOT/data/raw/nhanes/<cycle_letter>/
# Modules required: DEMO, RXQ_RX, RXQ_DRUG, DPQ, BMX, BPX, TCHOL, TRIGLY,
# HDL, GLU, GHB, SMQ, ALQ, DIQ, BPQ, MCQ, HIQ, OCQ.
# CYCLES and MODULES tables are in environment.R.

# Step 2: download the MEPS files (free, no registration required) and place
# under $STUDY_ROOT/data/raw/meps/:
#   h181.ssp   2015 Full-Year Consolidated
#   h178a.ssp  2015 Prescribed Medicines
#   h192.ssp   2016 Full-Year Consolidated
#   h188a.ssp  2016 Prescribed Medicines

export STUDY_ROOT=~/nhanes_work

Rscript -e 'source("environment.R")'
Rscript scripts/01_build_analytic.R       # NHANES analytic parquet
Rscript scripts/02_main_analyses.R        # NHANES Oaxaca + CDC 2016 stability
Rscript scripts/03_sensitivity.R          # MICE / E-value / NCO / MEPS ratio
Rscript scripts/04_figures.R              # Manuscript figures
Rscript scripts/05_class_decomposition.R  # NHANES 5-class Oaxaca (Figure 3)
Rscript scripts/06_meps_margin.R          # MEPS pooled margin decomposition (Figure 1)
```

End-to-end runtime is about 25–35 minutes on a 2020 M1 MacBook Pro; the MICE step in `03_sensitivity.R` and the 1,000-replicate bootstrap in `06_meps_margin.R` dominate.

## Analytic choices worth knowing before reuse

Opioid identification uses Multum generic drug-class codes 60, 191, 57 and 58 from RXQ_DRUG in NHANES; in MEPS, opioid identification is via case-insensitive regular-expression matching on the RxName field, with naloxone-containing combinations (medication-assisted treatment) excluded. Tramadol is included throughout; its 2014 DEA rescheduling is tested as a stratification in `03_sensitivity.R`.

Race-ethnicity is non-Hispanic White versus non-Hispanic Black (NHANES `RIDRETH1` collapsed; MEPS `RACETHX` codes 2 and 3). Hispanic and Asian respondents are not the primary contrast and appear only in descriptive tables.

The decomposition is Oaxaca–Blinder twofold with Reimers pooled weights, fitted cycle by cycle and pooled by inverse-variance-weighted fixed-effects meta-analysis. Implementation: `oaxaca` v0.1.4 on top of `survey` v4.2. The NHANES survey design is `svydesign(id = ~psu, strata = ~strata, weights = ~wt_comb, nest = TRUE)` with `survey.lonely.psu = "adjust"`. MEC weights are divided by the number of pooled cycles per NCHS guidance; see `compute_combined_weight()` in `environment.R`.

The MEPS margin decomposition (`06_meps_margin.R`) splits annual opioid receipt into chronic (cumulative days-of-supply ≥ 90) and short-course (> 0 and < 90 days) components, with 1,000-replicate person-level bootstrap confidence intervals around the share of the extensive-margin gap attributable to each component.

## Caveats

Results apply to non-institutionalised US civilian adults; they do not transfer to inpatient, correctional or pediatric settings without recalibration.

The 2019–March 2020 pre-pandemic subsample (`P`) and the August 2021–August 2023 release (`L`) use non-standard weighting structures and are excluded from the primary 1999–2018 window. Extending the series requires reconstructing combined weights following current NCHS guidance.

MEPS Full-Year Consolidated files for 2017 and 2018 are distributed in compressed SAS transport (CPORT) format that the open-source statistical environment used here cannot read. The MEPS analysis is therefore restricted to 2015 and 2016. Investigators with SAS access can decompress these files and extend the MEPS pooled period.

NHANES opioid use is captured by pill-bottle review at the MEC examination over a past-30-day reference window; MEPS opioid use is captured by pharmacy-confirmed any-receipt over the calendar year. The two ascertainment frames measure different aspects of prescribing exposure; the cross-database concordance reported in the manuscript is on direction and broad magnitude, not on point-prevalence values.

## License

Code (`.R` files, `environment.R`) under the [MIT License](LICENSE); documentation under CC BY 4.0. Please cite the manuscript (§ Citation) when reusing either.

## Citation

To be updated after DOI assignment:

> [USER: INSERT Authors]. The Black-White Gap In U.S. Opioid Prescribing Concentrates At The Short-Course Margin: A Cross-Database NHANES-MEPS Decomposition With Implications For The 2016 CDC Guideline. *Pharmacoepidemiology and Drug Safety* 2026 (under review). DOI: [to be assigned on acceptance].

For the exact code version used in the manuscript, cite the Zenodo release:

> [USER: INSERT after Zenodo release] — Zenodo DOI: 10.5281/zenodo.XXXXXXX

## Contact

Please open a GitHub issue for code-related questions. For data interpretation, the NCHS and AHRQ public helpdesks are the authoritative sources.
