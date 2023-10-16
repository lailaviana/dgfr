#' Example of fasta file that should be used as input
#'
#' This database contains sequences in fasta format
#'
#' @format Fasta file containing headers and sequences
#' \describe{
#'
#' }
"fasta"

#' Example of how should be the output of the distance matrix file
#'
#' This database contains the distance matrix generated from BLAST file
#'
#' @format Large matrix
"distance_matrix"

#' Example of how should be the output of a dimension reduced file
#'
#' This database contains the position of each sequence in 2 dimensions
#'
#' @format A dataframe containing rownames and 2 columns in the following order
#' \describe{
#'    \item{axis_1}{position in axis 1}
#'    \item{axis_2}{position in axis 2}
#' }
"dim_reduced_file"
