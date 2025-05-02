test_that("negative N is rejected", {
  expect_error(multipool(make_dummy_pool(seed=1),
                         make_dummy_pool(seed=2),
                         N = -5), ">= 1")
})
