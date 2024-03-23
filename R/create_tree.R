#' Function to generate a tree based on the distance matrix provided
#'
#'
#' This function creates a tree based on the distance matrix provided
#' and color the sequences based on the cluster they belong to
#'
#' @param distance_matrix distance matrix
#' @param kmeans_output k-means file
#'
#' @return Returns a tree with the sequences colored based on the cluster they belong to
#' @export
#'
#' @examples
#'
#' create_tree(distance_matrix, kmeans_output)
#'
create_tree <- function(distance_matrix, kmeans_output) {
  cluster_tree <- stats::hclust(stats::as.dist(distance_matrix))
  ggtree::ggtree(cluster_tree, layout="circular") %<+%
    kmeans_output + ggtree::geom_tippoint(ggtree::aes(color=as.character(cluster))) +
    ggplot2::labs(color = "Cluster")
}
