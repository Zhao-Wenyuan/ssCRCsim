# Internal: accept named numeric vector or genes x samples numeric matrix.
.validate_expression <- function(expr) {
  if (is.numeric(expr) && is.null(dim(expr))) {
    gene_names <- names(expr)
    if (is.null(gene_names) || length(gene_names) != length(expr)) {
      stop("A single-sample vector must have gene-symbol names.", call. = FALSE)
    }
    expr <- matrix(as.numeric(expr), ncol = 1L,
                   dimnames = list(gene_names, "Sample1"))
  } else if (is.data.frame(expr)) {
    if (!all(vapply(expr, is.numeric, logical(1)))) {
      stop("All columns of a data.frame input must be numeric.", call. = FALSE)
    }
    expr <- as.matrix(expr)
  } else if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("Input must be a named numeric vector, numeric matrix, or numeric data.frame.",
         call. = FALSE)
  }

  if (nrow(expr) < 1L || ncol(expr) < 1L) {
    stop("Expression input must contain at least one gene and one sample.", call. = FALSE)
  }
  if (is.null(rownames(expr))) {
    stop("Expression matrix must have gene symbols as row names.", call. = FALSE)
  }
  gene_names <- toupper(trimws(rownames(expr)))
  if (anyNA(gene_names) || any(!nzchar(gene_names))) {
    stop("Gene symbols cannot be missing or empty.", call. = FALSE)
  }
  if (anyDuplicated(gene_names)) {
    stop("Duplicated gene symbols after case/whitespace normalization. Resolve them yourself; this package never averages duplicates.",
         call. = FALSE)
  }
  rownames(expr) <- gene_names

  sample_names <- colnames(expr)
  if (is.null(sample_names)) {
    if (ncol(expr) == 1L) {
      sample_names <- "Sample1"
    } else {
      stop("Multi-sample matrix must have column names (sample IDs).", call. = FALSE)
    }
  }
  sample_names <- trimws(sample_names)
  if (anyNA(sample_names) || any(!nzchar(sample_names)) ||
      anyDuplicated(sample_names)) {
    stop("Sample IDs must be non-empty and unique.", call. = FALSE)
  }
  colnames(expr) <- sample_names
  expr
}
