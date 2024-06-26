% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getOmaBrowser.R
\name{getOmaDataForOneOrtholog}
\alias{getOmaDataForOneOrtholog}
\title{Get taxonomy ID, sequence and annotation for one OMA protein}
\usage{
getOmaDataForOneOrtholog(id = NULL)
}
\arguments{
\item{id}{oma ID of one protein}
}
\value{
Data frame contains the input protein ID with its taxonomy ID,
sequence, length and domain annotations (tab delimited) for input OMA protein
}
\description{
Get taxonomy ID, sequence and annotation for one OMA protein
}
\examples{
### Uncomment the following line to run the function
# getOmaDataForOneOrtholog("HUMAN29397")
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
