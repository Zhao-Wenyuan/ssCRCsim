# ssCRCsim 


**ssCRCsim** is an R package for assessing colorectal cancer (CRC) transcriptional similarity and predicting consensus molecular subtypes (CMS) using relative expression orderings (REOs).

The package provides two independent functions:

- `crc_similarity()`: Quantifies CRC transcriptional similarity using CRC-specific gene-pair REO rules.
- `predict_cms()`: Predicts CMS1–CMS4 using a fixed, balanced set of CMS-specific REO rules.

Both functions operate on individual samples and do not require an external patient reference expression matrix.

## Installation

Install ssCRCsim directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("Zhao-Wenyuan/ssCRCsim")

library(ssCRCsim)
```

## Input requirements

The package accepts:

- A numeric expression matrix with genes as rows and samples as columns.
- A named numeric vector representing a single sample.
- Gene symbols as row names or vector names.

Input expression values may be unnormalized or normalized using methods that preserve within-sample gene expression orderings.

Duplicate gene symbols must be resolved before analysis. The package does not automatically average or merge duplicate genes.

## 1. CRC similarity scoring

The `crc_similarity()` function calculates a continuous CRC similarity score based on CRC-specific REO rules.

```r
result <- crc_similarity(expr, threshold = 90)
```

The default REO threshold is 90, and thresholds from 70 to 99 are supported.

For each sample, the score is calculated as:

**CRC similarity score = N_matched / N_available**

where:

- `N_available`: Number of gene pairs for which both genes have finite expression values and unequal expression levels.
- `N_matched`: Number of available gene pairs matching the CRC-specific REO direction.

The function returns:

- `score`: CRC similarity score.
- `n_available`: Number of available gene pairs.
- `n_matched`: Number of matched gene pairs.
- `qc_pass`: Quality-control status.

A minimum of 200 available gene pairs is required. Samples failing this criterion receive an `NA` score.

## 2. CMS prediction

The `predict_cms()` function predicts CRC consensus molecular subtypes using fixed, balanced REO rules established at a 73% consistency threshold.

```r
cms_result <- predict_cms(expr)
```

The function returns CMS1–CMS4 similarity scores and a predicted CMS subtype.

CMS prediction is independent of the threshold used for CRC similarity scoring.

## Example: Expression matrix

```r
# Read a gene-by-sample expression matrix
expr <- as.matrix(
  read.csv(
    "expression.csv",
    row.names = 1,
    check.names = FALSE
  )
)

# CRC similarity scoring
crc_result <- crc_similarity(expr)

# CMS subtype prediction
cms_result <- predict_cms(expr)
```

## Methodological notes

- CRC similarity is assessed using within-sample REO patterns rather than direct comparisons of absolute expression levels between samples.
- The default CRC scoring threshold is 90%.
- CMS classification uses fixed balanced REO73 rules.
- No Spearman correlation or delta-Spearman calculation is required for CMS prediction.
- The package includes the required REO rules internally.

## Version

Current version: **0.1.0**
