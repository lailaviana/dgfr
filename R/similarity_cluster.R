#' Function to assess similarity within clusters
#'
#' This function assess the mean/median similarity of each cluster
#'
#' @param score_file file containing the scores of each pairwise alignment
#' @param kmeans_output k-means file
#' @param output_type type of output (table or plot)

#' @return It returns a dataframe or a ggplot figure
#' @export
#'
#' @examples
#'
#' similarity_cluster(score_file, kmeans_output, output_type = "plot")
#' similarity_cluster(score_file, kmeans_output, output_type = "table")
#'
similarity_cluster <- function(score_file, kmeans_output, output_type) {
  res <- score_file |>
    tibble::rownames_to_column(var = "gene1") |>
    tidyr::pivot_longer(cols = -1,
                        names_to = "gene2",
                        values_to = "similarity") |>
    dplyr::left_join(kmeans_output, by = dplyr::join_by(gene1 == name)) |>
    dplyr::rename(cluster_gene1 = cluster) |>
    dplyr::select(gene1, gene2, similarity, cluster_gene1) |>
    dplyr::left_join(kmeans_output, by = dplyr::join_by(gene2 == name)) |>
    dplyr::rename(cluster_gene2 = cluster) |>
    dplyr::select(gene1, gene2, similarity, cluster_gene1, cluster_gene2) |>
    dplyr::mutate(cluster = dplyr::case_when(
      cluster_gene1 == 1 & cluster_gene2 == 1 ~ "cluster_1",
      cluster_gene1 == 2 & cluster_gene2 == 2 ~ "cluster_2",
      cluster_gene1 == 3 & cluster_gene2 == 3 ~ "cluster_3",
      cluster_gene1 == 4 & cluster_gene2 == 4 ~ "cluster_4",
      cluster_gene1 == 5 & cluster_gene2 == 5 ~ "cluster_5",
      cluster_gene1 == 6 & cluster_gene2 == 6 ~ "cluster_6",
      cluster_gene1 == 7 & cluster_gene2 == 7 ~ "cluster_7",
      cluster_gene1 == 8 & cluster_gene2 == 8 ~ "cluster_8",
      cluster_gene1 == 9 & cluster_gene2 == 9 ~ "cluster_9",
      cluster_gene1 == 10 & cluster_gene2 == 10 ~ "cluster_10",
      cluster_gene1 == 11 & cluster_gene2 == 11 ~ "cluster_11",
      cluster_gene1 == 12 & cluster_gene2 == 12 ~ "cluster_12",
      cluster_gene1 == 13 & cluster_gene2 == 13 ~ "cluster_13",
      cluster_gene1 == 14 & cluster_gene2 == 14 ~ "cluster_14",
      cluster_gene1 == 15 & cluster_gene2 == 15 ~ "cluster_15",
      cluster_gene1 == 16 & cluster_gene2 == 16 ~ "cluster_16",
      cluster_gene1 == 17 & cluster_gene2 == 17 ~ "cluster_17",
      cluster_gene1 == 18 & cluster_gene2 == 18 ~ "cluster_18",
      cluster_gene1 == 19 & cluster_gene2 == 19 ~ "cluster_19",
      cluster_gene1 == 20 & cluster_gene2 == 20 ~ "cluster_20",
      TRUE ~ "no")) |>
    dplyr::filter(cluster != "no")

  if(output_type == "plot"){
    fig <- res |>
      dplyr::group_by(cluster) |>
      ggplot2::ggplot(ggplot2::aes(x = cluster, y = similarity, fill = cluster)) +
      ggplot2::geom_violin(show.legend = FALSE) +
      ggplot2::geom_boxplot(width=0.01, outlier.shape = NA, show.legend = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Cluster", y = "Similarity") +
      ggplot2::theme(text = ggplot2::element_text(size = 15))
    return(fig)
  }

  else if(output_type == "table") {
    table <- res |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(mean_similarity = mean(similarity),
                       median_similarity = stats::median(similarity),
                       number_of_genes = dplyr::n_distinct(gene1))
    return(table)
  }
  else(
    stop("The argument output_type must be 'plot' or 'table'.")
  )

}
