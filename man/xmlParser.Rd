% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseMainInput.R
\name{xmlParser}
\alias{xmlParser}
\title{Parse orthoXML input file}
\usage{
xmlParser(inputFile = NULL)
}
\arguments{
\item{inputFile}{input file in xml format}
}
\value{
A data frame of input data in long-format containing seed gene IDs (
or orthologous group IDs), their orthologous proteins together with the
corresponding taxonomy IDs and values of (up to) two additional variables.
}
\description{
Parse orthoXML input file
}
\examples{
inputFile <- system.file(
    "extdata", "test.main.xml", package = "PhyloProfile", mustWork = TRUE
)
xmlParser(inputFile)
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
