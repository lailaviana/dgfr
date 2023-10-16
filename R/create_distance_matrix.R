#' Function to create a distance matrix
#'
#' This function creates a distance matrix from the pairwise alignment score file
#'
#' @param score_file dataframe containing the score of pairwise alignment
#'
#' @return It returns a large matrix
#' @export
#'
#' @examples
#'
#' create_distance_matrix(score_file)

create_distance_matrix <- function(score_file) {
  distance_matrix <- dist(score_file, method = "euclidean") |> as.matrix()
  return(distance_matrix)
}
