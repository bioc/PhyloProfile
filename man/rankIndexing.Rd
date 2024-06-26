% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createTaxonomyMatrix.R
\name{rankIndexing}
\alias{rankIndexing}
\title{Indexing all available ranks (including norank)}
\usage{
rankIndexing(rankListFile = NULL)
}
\arguments{
\item{rankListFile}{Input file, where each row is a rank list of a taxon
(see rankListFile in example)}
}
\value{
A dataframe containing a list of all possible ranks and their indexed
values.
}
\description{
Indexing all available ranks (including norank)
}
\examples{
rankListFile <- system.file(
    "extdata", "data/rankList.txt", package = "PhyloProfile", mustWork = TRUE
)
PhyloProfile:::rankIndexing(rankListFile)
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
