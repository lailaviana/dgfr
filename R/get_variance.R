#' Function to get explained variance by each axis
#'
#' This function gets the variance explained by each axis in the distance
#' matrix
#'
#' @param distance_matrix distance matrix
#'
#' @return Returns a dataframe containing 1 observation and n variables, depending of the number of input sequences
#' @export
#'
#' @examples
#'
#' get_variance(distance_matrix)
#'
get_variance <- function(distance_matrix) {
  stats::prcomp(distance_matrix, scale = TRUE) |>
    factoextra::get_eig() |>
    tibble::as_tibble() |>
    dplyr::select(variance.percent) |>
    tibble::rownames_to_column(var = "dim") |>
    tidyr::pivot_wider(id_cols = tidyr::everything(),
                       names_from = dim,
                       values_from = variance.percent)
}
