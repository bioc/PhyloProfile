% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getOmaBrowser.R
\name{getDataForOneOma}
\alias{getDataForOneOma}
\title{Get OMA info for a query protein and its orthologs}
\usage{
getDataForOneOma(seedID = NULL, orthoType = "OG")
}
\arguments{
\item{seedID}{OMA protein ID}

\item{orthoType}{type of OMA orthologs ("OG" or "HOG"). Default = "OG".}
}
\value{
Data frame contains info for all sequences of the input OMA group
(or HOG). That info contains the protein IDs, taxonomy IDs, sequences,
lengths, domain annotations (tab delimited) and the corresponding seed ID.
}
\description{
Get taxonomy IDs, sequences, length and annotations for an OMA
orthologous group (or OMA HOG).
}
\examples{
### Uncomment the following line to run the function
# getDataForOneOma("HUMAN29397", "OG")
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
