test_that("dim_reduction works", {
  dim_reduced_file <- dim_reduction(distance_matrix)
  #Test class
  expect_s3_class(dim_reduced_file, "data.frame")
  #Test column number
  expect_equal(ncol(dim_reduced_file), 2)
  #Test variable classes
  expect_equal(class(dim_reduced_file$axis_1), "numeric")
  expect_equal(class(dim_reduced_file$axis_2), "numeric")
  })
