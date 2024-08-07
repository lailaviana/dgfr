% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{calc_similarity}
\alias{calc_similarity}
\title{calculate similarity scores}
\usage{
calc_similarity(seq1, seq2, alignment_method, substitutionMatrix, myseq)
}
\arguments{
\item{seq1}{A value or the magrittr placeholder.}

\item{seq2}{A function call using the magrittr semantics.}

\item{alignment_method}{type of alignment local or global}

\item{substitutionMatrix}{substitution matrix to be used}

\item{myseq}{sequences to be used}
}
\value{
a dataframe with the similarity scores
}
\description{
This function calculates the similarity scores between two sequences
}
