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

#' Pipe operator for ggtree
#'
#' @name %<+%
#' @rdname pipe_ggtree
#' @keywords internal
#' @export
#' @importFrom ggtree %<+%
NULL

utils::globalVariables(c("pident", "qcovs", "query", "subject", "distance",
                         "genes", "gene1", "mean_distance", "variance.percent",
                         "kmeans$cluster", "V1", "V2", "dist", "data", "seq1",
                         "seq2", "value.x", "value.y", "gene2", "similarity",
                         "everything", "starts_with", "name", "cluster",
                         "cluster_gene1", "cluster_gene2", "median", "myseq"
                         ))


# auxiliary functions -----------------------------------------------------

#' calculate similarity scores
#' @description
#' This function calculates the similarity scores between two sequences
#' @param seq1 A value or the magrittr placeholder.
#' @param seq2 A function call using the magrittr semantics.
#' @param alignment_method type of alignment local or global
#' @param substitutionMatrix substitution matrix to be used
#' @param myseq sequences to be used
#' @export
#' @returns a dataframe with the similarity scores

calc_similarity <- function(seq1, seq2, alignment_method, substitutionMatrix, myseq) {
  algn_12 <- Biostrings::pairwiseAlignment(
    pattern = myseq[seq1],
    subject = myseq[seq2],
    type = alignment_method,
    substitutionMatrix = substitutionMatrix,
    gapOpening = 2,
    gapExtension = 0.5,
    scoreOnly = TRUE
  )

  algn_11 <- Biostrings::pairwiseAlignment(
    pattern = myseq[seq1],
    subject = myseq[seq1],
    type = alignment_method,
    substitutionMatrix = substitutionMatrix,
    gapOpening = 2,
    gapExtension = 0.5,
    scoreOnly = TRUE
  )

  algn_22 <- Biostrings::pairwiseAlignment(
    pattern = myseq[seq2],
    subject = myseq[seq2],
    type = alignment_method,
    substitutionMatrix = substitutionMatrix,
    gapOpening = 2,
    gapExtension = 0.5,
    scoreOnly = TRUE
  )

  similarity <- algn_12 / (sqrt(algn_11) * sqrt(algn_22))

  return(similarity)
}
