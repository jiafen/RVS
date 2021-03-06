% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/Raw_process.R
\name{raw_reorder}
\alias{raw_reorder}
\title{raw_reorder should be called if both case and controls are in the same dataset.}
\usage{
raw_reorder(path1, path2, file3)
}
\arguments{
\item{path1:}{input for additive coding genotype files}

\item{path2:}{the file including the list of IIDs for cases}

\item{file3:}{the file name you want your processed genotype and phenotype to be named.}
}
\value{
a list including- geno: the genotype data with the first ncase rows for cases, genotype data starting from 7th column, ncase: Num. of case
}
\description{
It reorganize the cases' and controls' genotype, and return genotype with the first ncase rows are for case followed by ncont controls
}

