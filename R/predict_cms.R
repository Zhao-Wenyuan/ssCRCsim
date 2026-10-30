#' Predict CRC consensus molecular subtype using a fixed REO73 panel
#'
#' Predicts CMS1-CMS4 from four independent, equally sized rule sets. The
#' fixed construction threshold of 73 percent is NOT a sample-level cutoff.
#' Every CMS must have at least 200 evaluable rules. A unique maximum score
#' is required for a subtype prediction. CMS scores are NOT probabilities.
#'
#' @param expr Numeric genes-by-samples matrix/data.frame with gene-symbol row
#'   names; or a named numeric vector for one sample.
#' @return A data.frame with four CMS scores, evaluable counts, predicted CMS,
#'   status and score margin. Failed QC or tied highest scores yield NA subtype.
#' @export
predict_cms <- function(expr) {
  expr <- .validate_expression(expr)
  cms <- c("CMS1", "CMS2", "CMS3", "CMS4")
  scores <- lapply(cms, function(z) .score_reo_core(expr, .cms_rules[[z]]))
  names(scores) <- cms
  out <- data.frame(sample_id = colnames(expr), stringsAsFactors = FALSE)
  for (z in cms) {
    out[[paste0(z, "_score")]] <- scores[[z]]$score
    out[[paste0(z, "_n_available")]] <- scores[[z]]$n_available
  }
  qc <- do.call(cbind, lapply(scores, function(z) z$qc_pass))
  all_qc <- rowSums(qc) == length(cms)
  score_mat <- do.call(cbind, lapply(scores, function(z) z$score))
  colnames(score_mat) <- cms
  predicted <- rep(NA_character_, nrow(out))
  status <- rep("LOW_REO_COVERAGE", nrow(out))
  margin <- rep(NA_real_, nrow(out))
  for (i in seq_len(nrow(out))) {
    if (!all_qc[i]) next
    z <- score_mat[i, ]
    if (any(!is.finite(z))) {
      status[i] <- "MISSING_SCORE"
      next
    }
    best <- which(z == max(z))
    if (length(best) != 1L) {
      status[i] <- "TIE"
      next
    }
    predicted[i] <- cms[best]
    status[i] <- "OK"
    margin[i] <- sort(z, decreasing = TRUE)[1L] -
                 sort(z, decreasing = TRUE)[2L]
  }
  out$predicted_CMS <- predicted
  out$qc_pass <- all_qc
  out$status <- status
  out$score_margin <- margin
  out
}
