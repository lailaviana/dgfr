#' Function to estimate optimal number of clusters and clustering using kmeans
#'
#' This function estimates the optimal number of clusters in the input sample based on dimension reduced file. In case there are two or more cluster numbers supported by the same index number, the smaller number of clusters will be chosen to represent the family.
#'
#'
#' @param dim_reduced_file dimension reduced file
#' @param min_clust minimum number of clusters
#' @param max_clust maximum number of clusters
#'
#' @return A dataframe containing name of the seq, axis_1, axis_2 and cluster number)
#' @export
#'
#' @examples
#'
#' kmeans_clustering(dim_reduced_file)
#'
kmeans_clustering <- function(dim_reduced_file, min_clust = 2, max_clust = 20) {
  NbClust_ <- purrr::quietly(NbClust::NbClust)
  cluster_n <- NbClust_(dim_reduced_file, diss = NULL,
                        distance = "euclidean",
                        min.nc = min_clust,
                        max.nc = max_clust,
                        method = "kmeans")$result

  optimal_cluster_n <- cluster_n$Best.nc |>
    tibble::as_tibble() |>
    dplyr::slice(1) |>
    tidyr::pivot_longer(tidyr::everything(),
                        names_to = "test",
                        values_to = "cluster_n") |>
    dplyr::count(cluster_n) |>
    dplyr::top_n(1) |>
    dplyr::pull(cluster_n)

  optimal_cluster_n <- optimal_cluster_n[1]

  kmeans <- stats::kmeans(dim_reduced_file,
                          optimal_cluster_n,
                          nstart=100,
                          iter.max=100000)

  id_cluster <- kmeans$cluster |> as.data.frame() |>
    tibble::rownames_to_column(var = "name") |>
    dplyr::rename(cluster = `kmeans$cluster`)

  dim_reduced_file |> tibble::rownames_to_column(var = "name") |>
    dplyr::left_join(id_cluster)
}
