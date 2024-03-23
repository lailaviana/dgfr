## code to prepare `fasta` dataset goes here
fasta <- Biostrings::readAAStringSet("data-raw/fasta")
usethis::use_data(fasta, overwrite = TRUE)

## code to prepare `tcmuc.fasta` dataset goes here
tcmuc_fasta <- Biostrings::readAAStringSet("data-raw/tcmuc_fasta")
usethis::use_data(tcmuc_fasta, overwrite = TRUE)


## code to prepare `score_file` dataset goes here
score_file <- dgfr::get_alignment_score(fasta = fasta,
                                  type = "prot",
                                  alignment_method = "global")

usethis::use_data(score_file, overwrite = TRUE)

## code to prepare `distance_matrix` dataset goes here
distance_matrix <- dgfr::create_distance_matrix(score_file)

usethis::use_data(distance_matrix, overwrite = TRUE)

## code to prepare `dim_reduced_file` dataset goes here
dim_reduced_file <- dgfr::dim_reduction(distance_matrix = distance_matrix)
usethis::use_data(dim_reduced_file, overwrite = TRUE)

## code to prepare `kmeans_output` dataset goes here
kmeans_output <- dgfr::kmeans_clustering(dim_reduced_file = dim_reduced_file)
usethis::use_data(kmeans_output, overwrite = TRUE)


