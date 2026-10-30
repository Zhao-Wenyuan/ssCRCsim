# ssCRCsim (development skeleton)

Two independent user functions:

- `crc_similarity(expr, threshold = 90)` — CRC-specific REO score; thresholds 70–99.
- `predict_cms(expr)` — fixed balanced REO73 CMS1–CMS4 prediction.

No Spearman, Delta Spearman, or automatic duplicate gene averaging.

## IMPORTANT: before installation

The raw rule files are NOT bundled in the development source. From RStudio, open
`D:/ssCRCsim/ssCRCsim` as a project and run:

```r
source("data-raw/build_rules.R")
```

It reads the three RDS files from `D:/ssCRCsim`, compares REO70-derived REO90
against the saved CRC90 reference and **stops on inconsistency**. It creates
`R/sysdata.rda` only when checks pass.

Then:

```r
install.packages(c("devtools", "roxygen2", "testthat"))
roxygen2::roxygenise()
devtools::load_all()
devtools::test()
devtools::check()
devtools::install()

library(ssCRCsim)
# Assuming x is your numeric gene x sample expression matrix:
# crc_similarity(x)
# crc_similarity(x, threshold = 80)
# predict_cms(x)
```

Before public release, edit the author name and email in `DESCRIPTION` and
`LICENSE`; add manuscript citation and data provenance details.

## Input

Provide a **numeric gene-by-sample** matrix (`rownames` = HGNC-style gene
symbols, `colnames` = sample IDs) or a **named numeric vector** for one sample.
Duplicate gene symbols after uppercasing/trim are an **error**, not averaged.
NA/NaN/Inf and tied expression values cannot be evaluated. A minimum of 200
evaluable REO pairs applies to CRC score and separately to each CMS score.

## Notes

The CMS construction threshold of 73% is not a prediction score cutoff.
CMS scores do not represent probabilities. CRC threshold changes do not affect
CMS prediction.

This is a starter implementation and has not yet been run against your local
REO70 file or legacy scoring outputs; perform reference-result replication
before research/public release.
