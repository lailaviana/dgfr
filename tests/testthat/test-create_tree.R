test_that("create_tree works", {
  tree <- create_tree(distance_matrix, kmeans_output)
  #Test class
  expect_s3_class(kmeans_output, "data.frame")
  expect_true(is.matrix(distance_matrix))
  expect_s3_class(tree, "ggtree")
  #Test variable class
  expect_equal(class(kmeans_output$name), "character")
  expect_equal(class(kmeans_output$axis_1), "numeric")
  expect_equal(class(kmeans_output$axis_2), "numeric")
  expect_equal(class(kmeans_output$cluster), "integer")
})
