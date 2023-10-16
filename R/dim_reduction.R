#' Function to perform dimension reduction
#'
#' This function takes the distance matrix as input and reduce the dimension to 2 axis
#'
#' @param distance_matrix distance matrix
#'
#' @return A dataframe containing the id, position in x axis, position in y axis
#' @export
#'
#' @examples
#'
#' dim_reduction(distance_matrix)
#'
dim_reduction <- function(distance_matrix) {
  dim_reduced_file <- stats::cmdscale(distance_matrix,
                                      k = 2,
                                      eig = FALSE,
                                      add = FALSE,
                                      x.ret = FALSE) |>
    base::as.data.frame() |>
    dplyr::rename(axis_1 = V1, axis_2 = V2)

  dim_reduced_file
}
