% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locus_matrix.R
\name{getLocusMatrix}
\alias{getLocusMatrix}
\title{XIBD Binary IBD Matrix}
\usage{
getLocusMatrix(ped.genotypes, ibd.segments)
}
\arguments{
\item{ped.genotypes}{a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
See \code{Value} description in \code{\link{getGenotypes}} for more details.
The family IDs and individual IDs in \code{pedigree} must match the family IDs and individual IDs in the header of \code{genotypes}.}

\item{ibd.segments}{a named list containing the \code{ibd_segments} inferred from pairs of samples.
See \code{value} description in \code{\link{getIBDsegments}} for more details.}
}
\value{
A data frame the first four columns:
\enumerate{
\item Chromosome (type \code{"character"}, \code{"numeric"} or \code{"integer"})
\item SNP identifiers (type \code{"character"})
\item Genetic map distance (centi morgans, cM) (type \code{"numeric"})
\item Base-pair position (type \code{"integer"})
}
where each row describes a unique SNP. These columns are headed \code{chr, snp_id, pos_M} and \code{pos_bp} respectively.
Columns 5 onwards contain the binary IBD information for each sample pair, where a single column corresponds to a single pair.
These columns are labeled with merged family IDs and individual IDs separated by a slash symbol (/). For example famA/ind1/famA/ind2.
}
\description{
Returns a binary matrix of IBD (1) and non-IBD (0) for each SNP and pair combination.
}
\examples{
# generate a binary IBD matrix
my_locus_matrix <- getLocusMatrix(ped.genotypes = example_genotypes,
                                  ibd.segments = example_ibd)
}
