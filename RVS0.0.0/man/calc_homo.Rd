% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/help_assoc_test.R
\name{calc_homo}
\alias{calc_homo}
\title{Input}
\usage{
calc_homo(subgroup, nsnp, type)
}
\arguments{
\item{subgroup:}{a matrix for conditional expected genotype or genotype calls for case or controls. dimension nsample*nsnp}

\item{nsnp:}{number of snps}

\item{type:}{specify it is conditional expected genotype or genotype call itself}
}
\value{
a vector indicate the snp columns that are monormophic.
}
\description{
Input
}

