#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


utils::globalVariables(c("pident", "qcovs", "query", "subject", "distance",
                         "genes", "gene1", "mean_distance", "variance.percent",
                         "kmeans$cluster", "V1", "V2", "dist", "data", "seq1",
                         "seq2", "value.x", "value.y", "gene2", "similarity",
                         "everything", "starts_with", "name", "cluster",
                         "cluster_gene1", "cluster_gene2", "median"
                         ))
