test_that("similarity_cluster works", {
  similarity_plot <- similarity_cluster(score_file, kmeans_output, "plot")
  #Test class
  expect_s3_class(kmeans_output, "data.frame")
  expect_s3_class(score_file, "data.frame")
  #Test column number
  expect_equal(ncol(kmeans_output), 4)
  #Test variable class
  expect_equal(class(kmeans_output$name), "character")
  expect_equal(class(kmeans_output$axis_1), "numeric")
  expect_equal(class(kmeans_output$axis_2), "numeric")
  expect_equal(class(kmeans_output$cluster), "integer")
})
