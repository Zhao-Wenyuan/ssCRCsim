#' Compute a colorectal cancer-specific REO similarity score
#'
#' Uses the precomputed CRC-specific rules. The `threshold` is applied to both
#' the CRC-direction consistency and the non-CRC reverse-direction consistency
#' in the original reference cohorts (all-sample denominators). The score
#' is the fraction of evaluable REOs matching the CRC-specific direction.
#'
#' @param expr Numeric genes-by-samples matrix/data.frame with gene-symbol row
#'   names; or a named numeric vector for one sample.
#' @param threshold Integer percentage between 70 and 99 (default 90).
#' @return A data.frame with sample-level score, matched and available counts,
#'   coverage and QC (at least 200 evaluable gene pairs). Failed QC returns NA.
#' @export
crc_similarity <- function(expr, threshold = 90L) {
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold) ||
      !is.finite(threshold) || threshold != floor(threshold) ||
      threshold < 70 || threshold > 99) {
    stop("threshold must be an integer from 70 to 99.", call. = FALSE)
  }
  expr <- .validate_expression(expr)
  kept <- 100 * .crc_rules$n_crc_support >= threshold * .crc_counts[["CRC"]] &
          100 * .crc_rules$n_non_reverse >= threshold * .crc_counts[["nonCRC"]]
  rules <- .crc_rules[kept, c("gene_a", "gene_b", "direction"), drop = FALSE]
  out <- .score_reo_core(expr, rules)
  names(out)[names(out) == "score"] <- "CRC_similarity_score"
  out$threshold <- as.integer(threshold)
  out[, c("sample_id", "CRC_similarity_score", "n_matched", "n_available",
          "n_rules_total", "n_rules_with_both_genes", "coverage_fraction",
          "qc_pass", "threshold"), drop = FALSE]
}
