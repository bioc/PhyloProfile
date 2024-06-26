% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseMainInput.R
\name{fastaParser}
\alias{fastaParser}
\title{Parse multi-fasta input file}
\usage{
fastaParser(inputFile = NULL)
}
\arguments{
\item{inputFile}{input multiple fasta file. Check extdata/test.main.fasta or
https://github.com/BIONF/PhyloProfile/wiki/Input-Data#multi-fasta-format for
the supported FASTA header.}
}
\value{
A data frame of input data in long-format containing seed gene IDs (
or orthologous group IDs), their orthologous proteins together with the
corresponding taxonomy IDs and values of (up to) two additional variables.
}
\description{
Parse multi-fasta input file
}
\examples{
inputFile <- system.file(
    "extdata", "test.main.fasta", package = "PhyloProfile", mustWork = TRUE
)
fastaParser(inputFile)
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
