% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getOmaBrowser.R
\name{getSelectedFastaOma}
\alias{getSelectedFastaOma}
\title{Get selected fasta sequences from a raw OMA dataframe}
\usage{
getSelectedFastaOma(finalOmaDf = NULL, seqID = NULL)
}
\arguments{
\item{finalOmaDf}{raw OMA data for a list of proteins (see ?getDataForOneOma)}

\item{seqID}{OMA ID of selected protein}
}
\value{
Required protein sequence in fasta format.
}
\description{
Get selected fasta sequences from a raw OMA dataframe
}
\examples{
### Uncomment the following line to run the function
# omaData <- getDataForOneOma("HUMAN29397", "OG")
# getSelectedFastaOma(omaData, "HUMAN29397")
}
\seealso{
\code{\link{getDataForOneOma}}
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
