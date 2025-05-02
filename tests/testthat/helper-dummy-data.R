# tests/testthat/helper-dummy-data.R
# ------------------------------------------------------------------
# Helper available to **all** tests: generates a small synthetic pool
# ------------------------------------------------------------------

make_dummy_pool <- function(n_markers = 1000,
                            chr       = "chrI",
                            lambda_a  = 20,
                            lambda_b  = 20,
                            seed      = NULL,
                            qtl_start = 80000,
                            qtl_end   = 120000,
                            lambda_qtl_a = 35,
                            lambda_qtl_b = 10)
{
  # 1. marker positions
  if (!is.null(seed)) set.seed(seed)
  pos <- sort(sample.int(200000, n_markers, replace = TRUE))

  # 2. baseline counts
  a <- rpois(n_markers, lambda_a)
  b <- rpois(n_markers, lambda_b)

  # 3. spike-in a fake QTL (optional but useful for sanity checks)
  qtl_region <- pos > qtl_start & pos < qtl_end
  a[qtl_region] <- rpois(sum(qtl_region), lambda_qtl_a)
  b[qtl_region] <- rpois(sum(qtl_region), lambda_qtl_b)

  data.frame(chr = chr, pos = pos, a = a, b = b)
}
