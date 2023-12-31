
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dgfr package <img src="man/figures/dgfR.png" align="right" width = "120px"/>

<!-- badges: start -->
<!-- badges: end -->

The purpose of the dgfr package is to enable the visualization of
sequence diversity within gene/protein families. In addition, it will be
possible to determine the optimal number of clusters for grouping
members of gene families based on sequence similarity and to calculate
the mean/median similarity inside each cluster.

## Installation

Installing the development version of dgfr is as follow:

``` r
# install.packages("devtools")
#devtools::install_github("lailaviana/dgfr")
library(dgfr)
```

## How to use?

1.  First, retrieve the FASTA sequences of the gene family being studied
    from the database of your choice. The fasta file must be loaded
    inside R using the Biostrings::readAAStringSet() in case of protein
    sequences or Biostrings::readDNAStringSet() in case of nucleotide
    sequences.

``` r
fasta_file <- Biostrings::readAAStringSet("data-raw/fasta")
```

2.  The object containing your fasta should be provided to the first
    function named `get_alignment_score`, then you should also provide
    the type o sequence “nuc” for nucleotide or “prot” for protein and
    the alignment method (“global” or “local”). The score is calculated
    using the following formula:

- algn_seq1seq2 / (sqrt(algn_seq1seq1) \* sqrt(algn_seq2seq2))  
  This formula was first proposed by Bleakley and Yamanashi and is a way
  to normalize the score by the length of the alignment. The output will
  be a dataframe containing the score of each pairwise alignment among
  all sequences provided.

``` r
score_file <- dgfr::get_alignment_score(fasta = fasta,
                                        type = "prot",
                                        alignment_method = "global")
#> [1/5] Loading substitution matrix...
#> [2/5] Loading fasta file...
#> [3/5] Creating all possible sequence combinations...
#> [4/5] Performing pairwise alignment for prot sequences. This might take a while...
#> [5/5] Pairwise alignment is done.
```

3.  The output of the first function of the package, should be provided
    to the second function called `create_distance_matrix`. This
    function calculates the average Euclidian distance between the score
    of each pairwise alignment.

``` r
distance_matrix <- dgfr::create_distance_matrix(score_file = score_file)
```

4.  After generating the distance matrix, we will use the `get_variance`
    function to access the variability explained by each axis. This
    requires the provision of the file containing the distance matrix.
    In general, a cumulative value of 70% or higher for the axis is
    considered adequate for explaining the data’s variance.

``` r
variance <- dgfr::get_variance(distance_matrix = distance_matrix)
variance |> dplyr::select(1:3)
#> # A tibble: 1 × 3
#>     `1`   `2`   `3`
#>   <dbl> <dbl> <dbl>
#> 1  75.3  14.3  7.05
```

5.  The subsequent required step consists of a multidimensional scaling
    (MDS) of the data matrix. The output will consist of the name of the
    sequence and its position on axes 1 and 2.

``` r
dim_reduced_file <- dgfr::dim_reduction(distance_matrix = distance_matrix)
head(dim_reduced_file)
#>                        axis_1     axis_2
#> TcYC6_0128200-RA-p1 -3.754296 -1.3045277
#> TcYC6_0128210-RA-p1 -3.110871 -0.6916588
#> TcYC6_0128220-RA-p1 -4.136747 -1.6989440
#> TcYC6_0128230-RA-p1 -3.824319 -1.5401611
#> TcYC6_0128240-RA-p1 -4.329191 -1.6484587
#> TcYC6_0128250-RA-p1 -4.187994 -1.7303893
```

6.  A possible functionality of the dgfr package is to determine the
    optimal number of clusters in the dataset. You must provide the
    dimension-reduced file and, if desired, the minimum and maximum
    number of clusters to be identified. If you do not provide these
    values, the script will automatically search between 2 and 20
    clusters.

``` r
kmeans_output <- dgfr::kmeans_clustering(dim_reduced_file = dim_reduced_file)
```

Finally, you will have a table containing the name of each sequence, its
position on axes 1 and 2, and the cluster to which the sequence belongs.

``` r
head(kmeans_output)
#>                  name    axis_1     axis_2 cluster
#> 1 TcYC6_0128200-RA-p1 -3.754296 -1.3045277       3
#> 2 TcYC6_0128210-RA-p1 -3.110871 -0.6916588       4
#> 3 TcYC6_0128220-RA-p1 -4.136747 -1.6989440       3
#> 4 TcYC6_0128230-RA-p1 -3.824319 -1.5401611       3
#> 5 TcYC6_0128240-RA-p1 -4.329191 -1.6484587       3
#> 6 TcYC6_0128250-RA-p1 -4.187994 -1.7303893       3
```

7.  The next function, allows you to visualize the sequences and its
    cluster in a tree-like structure. You must provide the
    `kmeans_clustering` function output and the `create_distance_matrix`
    output. The output will be a plot containing the sequences and its
    cluster.

``` r
tree <- dgfr::create_tree(distance_matrix = distance_matrix,
                kmeans_output = kmeans_output)

tree
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="60%" style="display: block; margin: auto;" />

8.  Finally, you can access the mean/median similarity among each
    cluster using the function `similarity_cluster`. This function takes
    as argument the `kmeans_clustering` function output, the
    `get_alignment_score` output and you must also choose between two
    possible visualizations output. The “table” output returns a table
    containing the mean and median similarity among all clusters and the
    number of genes inside each cluster.

``` r
dgfr::similarity_cluster(score_file = score_file,
                         kmeans_output = kmeans_output,
                         output_type = "table")
#> # A tibble: 5 × 4
#>   cluster   mean_similarity median_similarity number_of_genes
#>   <chr>               <dbl>             <dbl>           <int>
#> 1 cluster_1           0.301             0.173              49
#> 2 cluster_2           0.953             0.976              78
#> 3 cluster_3           0.595             0.583             160
#> 4 cluster_4           0.405             0.326             215
#> 5 cluster_5           0.587             0.654             114
```

The “plot” output returns a plot containing a violin plot highlighting
the distribution of similarity within each cluster.

``` r
dgfr::similarity_cluster(score_file = score_file,
                         kmeans_output = kmeans_output,
                         output_type = "plot")
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="60%" style="display: block; margin: auto;" />

## Data visualisation suggestions

All visualization suggestions presented here were created using the
`ggplot2` package. However, you are allowed to use your preferred
package or tool to visualize the obtained data.

### Visualizing gene family dispersion

The simplest form of visualization is the scatter plot of sequences from
a particular gene family, where each point represents a sequence, and
the distance between two points is proportional to the similarity of
these sequences.

``` r
kmeans_output |> 
  ggplot2::ggplot(
    ggplot2::aes(x = axis_1, 
                 y = axis_2)) + 
  ggplot2::geom_point(size = 3.5, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = round(variance$`1`, 2), 
                y = round(variance$`2`, 2))
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="60%" style="display: block; margin: auto;" />

### Identifying particular sequences within the gene family

The second suggested plot is similar to the first, but it allows for the
highlighting of sequences of interest. In addition to `ggplot2`,
`gghighlight` and `ggrepel` are required to perform this analysis.

``` r
seq_name <- kmeans_output |> 
  dplyr::sample_n(5) |> 
  dplyr::select(name) |> 
  dplyr::pull()

kmeans_output |> 
  ggplot2::ggplot(
    ggplot2::aes(x = axis_1, 
                 y = axis_2)) + 
  ggplot2::geom_point(size = 3.5, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = round(variance$`1`, 2), 
                y = round(variance$`2`, 2)) +
gghighlight::gghighlight(name %in% seq_name) + 
  ggrepel::geom_text_repel(label = seq_name, 
                           min.segment.length = 0, 
                           seed = 42, 
                           box.padding = 0.5)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="60%" style="display: block; margin: auto;" />

### Identifying the number of clusters within the gene family

In this third visualization suggestion, the gene family is divided into
an optimal number of groups, each of which represented by a different
color.

``` r
kmeans_output |> 
  dplyr::mutate(cluster = as.character(cluster)) |> 
  ggplot2::ggplot(
    ggplot2::aes(x = axis_1, 
                 y = axis_2,
                 color = cluster)) + 
  ggplot2::geom_point(size = 3.5, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = round(variance$`1`, 2), 
                y = round(variance$`2`, 2))
```

<img src="man/figures/README-unnamed-chunk-15-1.png" width="60%" style="display: block; margin: auto;" />

### Aggregating and visualizing RNAseq data

Another functionality is to integrate RNAseq data to visualize the
expression and dispersion of the gene family in the scatter plot. The
size of the points is proportional to the gene expression level of each
family member, with small points representing genes with low expression
and large points representing genes with high expression.

``` r
kmeans_out_l2fc <- kmeans_output |> 
  dplyr::mutate(l2fc = runif(616, min = 0.1, max = 4))

kmeans_out_l2fc |> 
  dplyr::mutate(cluster = as.character(cluster)) |> 
  ggplot2::ggplot(
    ggplot2::aes(x = axis_1, 
                 y = axis_2,
                 color = cluster)) + 
  ggplot2::geom_point(size = kmeans_out_l2fc$l2fc, alpha = 0.3) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = round(variance$`1`, 2), 
                y = round(variance$`2`, 2))
```

<img src="man/figures/README-unnamed-chunk-16-1.png" width="60%" style="display: block; margin: auto;" />

## Warranty

If you find any bugs or incorrect information, please open an issue.

## Acknowledgment

The development of this package would not have been possible without the
knowledge shared by the instructors at
[curso-R](https://github.com/curso-r), mainly Beatriz Milz, Caio Lente
and Júlio Trecenti.
