% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsePhyloProfile.R
\name{getInputTaxaName}
\alias{getInputTaxaName}
\title{Get NCBI taxon names for a selected list of taxa}
\usage{
getInputTaxaName(rankName, taxonIDs = NULL, taxDB = NULL)
}
\arguments{
\item{rankName}{taxonomy rank (e.g. "species","phylum",...)}

\item{taxonIDs}{list of taxon IDs (e.g. ncbi1234). Default = NULL}

\item{taxDB}{Path to the taxonomy DB files}
}
\value{
Data frame contains a list of full names, taxonomy ranks and parent
IDs for the input taxa.
}
\description{
Get NCBI taxon names from
"PhyloProfile/data/taxonNamesReduced.txt" for a list of input taxa
}
\examples{
taxonIDs <- c("ncbi9606", "ncbi10116")
getInputTaxaName("species", taxonIDs)
}
\seealso{
\code{\link{getInputTaxaID}} for getting input taxon IDs,
\code{\link{getNameList}} for getting the full taxon name list
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
