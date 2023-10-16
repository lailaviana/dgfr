#' Function to perform pairwise alignment of all sequences and retrieve the alignment score
#'
#' This function performs pairwise alignment of all sequences and compute the score of each alignment
#'
#' @param fasta a fasta file containing all sequences from a multigene family
#' @param type if the sequences are from nucleotide (nuc) or protein (prot)
#' @param alignment_method method of sequence alignment (local) or (global)

#' @return It returns a dataframe containing the scores of each alignment
#' @export
#'
#' @examples
#'
#' get_alignment_score(fasta, "prot", "global")

get_alignment_score <- function(fasta, type, alignment_method) {

  # Loading substitution matrix ---------
  message("[1/5] Loading substitution matrix...")
  BLOSUM62 <- as.matrix(utils::read.table("ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62", check.names = FALSE))

  blast_matrix <- matrix(c(
    5, -4, -4, -4,
    -4, 5, -4, -4,
    -4, -4, 5, -4,
    -4, -4, -4, 5),
    nrow = 4, byrow = TRUE,
    dimnames = list(c("A", "T", "C", "G"), c("A", "T", "C", "G")))

  # Suppress summarise inform ----------
  options(dplyr.summarise.inform = FALSE)

  # Loading fasta file ---------
  message("[2/5] Loading fasta file...")
  myseq = fasta
  # Recovering the sequences names ------------
  names <- names(fasta)

  names_df <- names |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "seq") |>
    dplyr::mutate(seq = as.integer(seq))

  # Creating all possible combinations ------------
  message("[3/5] Creating all possible sequence combinations...")
  combinations <- expand.grid(seq1 = 1:length(myseq),
                              seq2 = 1:length(myseq))

  # Similarity Score ------------

  calc_similarity <- function(seq1, seq2, type, alignment_method) {
    seq1 <- myseq[seq1]
    seq2 <- myseq[seq2]

    algn_12 <- Biostrings::pairwiseAlignment(
      pattern = seq1,
      subject = seq2,
      type = alignment_method,
      substitutionMatrix = BLOSUM62,
      gapOpening = 2,
      gapExtension = 0.5,
      scoreOnly = TRUE)

    algn_11 <- Biostrings::pairwiseAlignment(
      pattern = seq1,
      subject = seq1,
      type = alignment_method,
      substitutionMatrix = BLOSUM62,
      gapOpening = 2,
      gapExtension = 0.5,
      scoreOnly = TRUE)

    algn_22 <- Biostrings::pairwiseAlignment(
      pattern = seq2,
      subject = seq2,
      type = alignment_method,
      substitutionMatrix = BLOSUM62,
      gapOpening = 2,
      gapExtension = 0.5,
      scoreOnly = TRUE)

    similarity <- algn_12 / (sqrt(algn_11) * sqrt(algn_22))

    return(similarity)
  }

  if (!type %in% c("nuc", "prot")) {
    stop("Invalid type. Please choose 'nuc' or 'prot'.")
  }

  message(paste("[4/5] Performing pairwise alignment for", type, "sequences. This might take a while..."))

  res <- combinations |>
    dplyr::group_by(seq1, seq2) |>
    dplyr::summarise(similarity = calc_similarity(seq1, seq2, type, alignment_method)) |>
    tibble::as_tibble() |>
    dplyr::left_join(names_df, by = c("seq1" = "seq")) |>
    dplyr::left_join(names_df, by = c("seq2" = "seq")) |>
    dplyr::rename(gene1 = value.x, gene2 = value.y) |>
    dplyr::select(-c(seq1, seq2))

  score_file <- res |>
    tidyr::pivot_wider(names_from = gene2, values_from = similarity) |>
    dplyr::mutate(gene1 = names) |>
    tibble::column_to_rownames(var = "gene1") |>
    dplyr::rename_with(~ names, everything())

  message("[5/5] Pairwise alignment is done.")

  return(score_file)
}

