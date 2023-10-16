test_that("create_distance_matrix works", {
  matrix_file <- create_distance_matrix(score_file)

  #Test class
  expect_true(is.matrix(matrix_file))
})
