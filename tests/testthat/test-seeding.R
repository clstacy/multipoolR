test_that("permutation FDR reproducible with seed", {
  res1 <- multipool(make_dummy_pool(seed=1), make_dummy_pool(seed=2),
                    N=50, nperm=5, seed = 999, parallel = FALSE, plot = FALSE)
  res2 <- multipool(make_dummy_pool(seed=1), make_dummy_pool(seed=2),
                    N=50, nperm=5, seed = 999, parallel = FALSE, plot = FALSE)
  expect_equal(res1$fdr_lod_threshold, res2$fdr_lod_threshold)
})
