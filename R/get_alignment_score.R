#' Function to perform pairwise alignment of all sequences and retrieve the alignment score
#'
#' This function performs pairwise alignment of all sequences and compute the score of each alignment
#'
#' @param fasta a fasta file containing all sequences from a multigene family
#' @param type if the sequences are from nucleotide (nuc) or protein (prot)
#' @param alignment_method method of sequence alignment (local) or (global)
#' @param cores number of cores to be used in the parallel computation

#' @return It returns a dataframe containing the scores of each alignment
#' @export
#'
#' @examples
#'
#' get_alignment_score(fasta, "prot", "global", 3)

get_alignment_score <- function(fasta, type, alignment_method, cores = 2) {

  # Loading substitution matrix ---------
  message("[1/6] Loading substitution matrix...")
  values <- c(
    4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4,
    -1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4,
    -2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4,
    -2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4,
    0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,
    -1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4,
    -1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4,
    0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4,
    -2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4,
    -1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4,
    -1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4,
    -1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4,
    -1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4,
    -2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4,
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4,
    1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4,
    0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4,
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4,
    -2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4,
    0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4,
    -2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4,
    -1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4,
    0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4,
    -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1
  )

  num_linhas <- 24
  num_colunas <- 24

  BLOSUM62 <- matrix(values, nrow = num_linhas, ncol = num_colunas, byrow = TRUE)

  nomes_aminoacidos <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*")
  rownames(BLOSUM62) <- colnames(BLOSUM62) <- nomes_aminoacidos

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
  message("[2/6] Loading fasta file...")
  myseq = fasta

  # Recovering the sequences names ------------
  new_headers <- sub(" .*", "", names(myseq))  # Remove tudo após o primeiro espaço

  # New headers to sequences -----------
  names(myseq) <- new_headers
  names <- names(myseq)

  names_df <- names |>
    tibble::as_tibble() |>
    tibble::rownames_to_column(var = "seq") |>
    dplyr::mutate(seq = as.integer(seq))

  # Creating all possible combinations ------------
  message("[3/6] Creating all possible sequence combinations...")
  combinations <- expand.grid(seq1 = 1:length(myseq),
                              seq2 = 1:length(myseq))

  future::plan(future::multisession(), workers = cores)

  if (!type %in% c("nuc", "prot")) {
    stop("Invalid type. Please choose 'nuc' or 'prot'.")
  }

  message(paste("[4/6] Performing pairwise alignment for", type, "sequences. This might take a while..."))

  if(type == "prot"){
    res <- combinations %>%
      furrr::future_pmap(
        ~calc_similarity(seq1 = .x, seq2 = .y,
                         alignment_method = alignment_method,
                         substitutionMatrix = BLOSUM62,
                         myseq = myseq)
      )
  } else {
    res <- combinations %>%
      furrr::future_pmap(
        ~calc_similarity(seq1 = .x, seq2 = .y,
                         alignment_method = alignment_method,
                         substitutionMatrix = blast_matrix,
                         myseq = myseq)
      )
  }

  res_df <- tibble::tibble(
    seq1 = combinations$seq1,
    seq2 = combinations$seq2,
    similarity = unlist(res)
  )

  # Realizar as junções e manipulações finais
  message("[5/6] Finishing up...")
  res_final <- res_df %>%
    dplyr::left_join(names_df, by = c("seq1" = "seq")) %>%
    dplyr::left_join(names_df, by = c("seq2" = "seq")) %>%
    dplyr::rename(gene1 = value.x, gene2 = value.y) %>%
    dplyr::select(-c(seq1, seq2))

  score_file <- res_final %>%
    tidyr::pivot_wider(names_from = gene2, values_from = similarity) %>%
    dplyr::mutate(gene1 = names) %>%
    tibble::column_to_rownames(var = "gene1") %>%
    dplyr::rename_with(~ names, everything())

  message("[6/6] Pairwise alignment is done.")

  on.exit(closeAllConnections())

  return(score_file)
}

