test_that("matched, unmatched, and tied pairs are counted correctly", {
  rules <- data.frame(
    gene_a = c("A", "A", "A"),
    gene_b = c("B", "C", "D"),
    direction = c(1L, -1L, 1L)
  )
  expr <- matrix(c(10, 2, 11, 10), ncol = 1,
                 dimnames = list(c("A", "B", "C", "D"), "X"))
  z <- ssCRCsim:::.score_reo_core(expr, rules, min_available = 1L)
  expect_equal(z$n_available, 2)
  expect_equal(z$n_matched, 2)
  expect_equal(z$score, 1)
})

test_that("low coverage returns NA", {
  rules <- data.frame(gene_a = "A", gene_b = "B", direction = 1L)
  expr <- matrix(c(2, 1), ncol = 1,
                 dimnames = list(c("A", "B"), "X"))
  z <- ssCRCsim:::.score_reo_core(expr, rules, min_available = 200L)
  expect_false(z$qc_pass)
  expect_true(is.na(z$score))
})

test_that("duplicate gene symbols after standardization are rejected", {
  expect_error(ssCRCsim:::.validate_expression(c(tp53 = 2, TP53 = 3)),
               "Duplicated gene symbols")
})

test_that("CRC threshold cannot change CMS classifier", {
  expr <- c(APC = 3, KRAS = 2, TP53 = 1)
  a <- predict_cms(expr)
  crc_similarity(expr, threshold = 70)
  crc_similarity(expr, threshold = 99)
  b <- predict_cms(expr)
  expect_identical(a, b)
})
