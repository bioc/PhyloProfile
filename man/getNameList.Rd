% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsePhyloProfile.R
\name{getNameList}
\alias{getNameList}
\title{Get list of pre-installed NCBI taxon names}
\usage{
getNameList(taxDB = NULL)
}
\arguments{
\item{taxDB}{Path to the taxonomy DB files}
}
\value{
List of taxon IDs, their full names, taxonomy ranks and parent IDs
obtained from "PhyloProfile/data/taxonNamesReduced.txt"
}
\description{
Get all NCBI taxon names from
"PhyloProfile/data/taxonNamesReduced.txt"
}
\examples{
getNameList()
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
