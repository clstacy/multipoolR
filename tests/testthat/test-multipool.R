library(multipoolR)

test_that("max LOD stable for synthetic data", {

  pool1 <- make_dummy_pool(seed = 123)
  pool2 <- make_dummy_pool(seed = 456)  # independent replicate

  res <- multipool(pool1, pool2,
                   N      = 50,
                   res    = 500,
                   nperm  = 0,
                   mode   = "replicates",
                   plot   = FALSE)

  # Capture the value once, update the expectation below, then it stays fixed.
  max_lod <- round(max(res$results$LOD), 3)
  expect_equal(max_lod, 8.1530)  # ← replace with printed value from first run
})
