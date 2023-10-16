test_that("get_variance works", {
  variance <- get_variance(distance_matrix)
  #Test class
  expect_s3_class(variance, "tbl_df")
  #Test row number
  expect_equal(nrow(variance), 1)
  #Test variable class
  expect_equal(class(variance$`1`), "numeric")
  expect_equal(class(variance$`2`), "numeric")
  expect_equal(class(variance$`3`), "numeric")
})
