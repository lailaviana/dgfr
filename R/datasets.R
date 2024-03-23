#' Example of a small fasta file that can be used as input
#'
#' This database contains sequences in fasta format
#'
#' @format small fasta file containing headers and sequences
#' \describe{
#'
#' }
"fasta"

#' Example of a complete fasta file that can be used as input
#'
#' This database contains sequences in fasta format
#'
#' @format complete fasta file containing headers and sequences
#' \describe{
#'
#' }
"tcmuc_fasta"

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
