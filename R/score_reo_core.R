# Internal scoring engine. direction: 1 means gene_a > gene_b, -1 means <.
.score_reo_core <- function(expr, rules, min_available = 200L,
                            block_size = 5000L) {
  n_samples <- ncol(expr)
  n_rules <- nrow(rules)
  idx_a <- match(rules$gene_a, rownames(expr))
  idx_b <- match(rules$gene_b, rownames(expr))
  both <- which(!is.na(idx_a) & !is.na(idx_b))

  n_available <- integer(n_samples)
  n_matched <- integer(n_samples)
  if (length(both) > 0L) {
    for (first in seq.int(1L, length(both), by = block_size)) {
      last <- min(length(both), first + block_size - 1L)
      r <- both[first:last]
      a <- expr[idx_a[r], , drop = FALSE]
      b <- expr[idx_b[r], , drop = FALSE]
      available <- is.finite(a) & is.finite(b) & a != b
      match_now <- matrix(FALSE, nrow = length(r), ncol = n_samples)
      gt <- which(rules$direction[r] == 1L)
      lt <- which(rules$direction[r] == -1L)
      if (length(gt)) {
        match_now[gt, ] <- a[gt, , drop = FALSE] > b[gt, , drop = FALSE]
      }
      if (length(lt)) {
        match_now[lt, ] <- a[lt, , drop = FALSE] < b[lt, , drop = FALSE]
      }
      match_now[!available] <- FALSE
      n_available <- n_available + colSums(available)
      n_matched <- n_matched + colSums(match_now)
    }
  }

  qc_pass <- n_available >= min_available
  score <- rep(NA_real_, n_samples)
  score[qc_pass] <- n_matched[qc_pass] / n_available[qc_pass]
  data.frame(
    sample_id = colnames(expr),
    score = score,
    n_matched = as.numeric(n_matched),
    n_available = as.numeric(n_available),
    n_rules_total = n_rules,
    n_rules_with_both_genes = length(both),
    coverage_fraction = if (n_rules > 0L) n_available / n_rules else NA_real_,
    qc_pass = qc_pass,
    stringsAsFactors = FALSE
  )
}
