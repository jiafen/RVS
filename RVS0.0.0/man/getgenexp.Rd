% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/help_vcf.R
\name{getgenexp}
\alias{getgenexp}
\title{Generates the expected probabilities of the genotypes E(G_ij|D_ij).}
\usage{
getgenexp(A0M, A1M, A2M, B0M, B1M, B2M, chr, Loc)
}
\arguments{
\item{A0M,}{Genotype lkelihood for major homozygous of the cases.  output of filterSNPs, extract from VCF files}

\item{A1M,}{Genotype lkelihood for minor heterozygous of the cases.}

\item{A2M,}{Genotype lkelihood for minor homozygous of the cases.}

\item{B0M,}{Genotype lkelihood for major homozygous of the controls.}

\item{B1M,}{Genotype lkelihood for minor heterozygous of the controls.}

\item{B2M,}{Genotype lkelihood for minor homozygous of the controls.}

\item{chr,}{chromosome of the snp, output from the  function filterSNPs}

\item{Loc,}{location of the snp, Output from the function filterSNPs.}
}
\value{
list including: 2 matrix (MM and P) and 1 data frame SNPs

MM -  Matrix contains expected genotypes for each person at each loci, row - individual, column - SNP

P -  Matrix for genotype frequency in all the sample at each loci, 3 columns contains P(G=0), P(G=1) and P(G=2), row - SNP

SNPs - a data frame contains chr and loc of the returned SNPs
}
\description{
Generates the expected probabilities of the genotypes E(G_ij|D_ij).
}
\examples{
geneexp <- getgenexp(SNPs$A0M,SNPs$A1M,SNPs$A2M,SNPs$B0M,SNPs$B1M,SNPs$B2M,SNPs$L11)
}

