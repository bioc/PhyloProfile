% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parseDomainInput.R
\name{getDomainFolder}
\alias{getDomainFolder}
\title{Get domain file from a folder for a seed protein}
\usage{
getDomainFolder(seed, domainPath)
}
\arguments{
\item{seed}{seed protein ID}

\item{domainPath}{path to domain folder}
}
\value{
Domain file and its complete directory path for the selected protein.
}
\description{
Get domain file from a folder for a seed protein
}
\examples{
domainPath <- paste0(
    path.package("PhyloProfile", quiet = FALSE), "/extdata/domainFiles"
)
PhyloProfile:::getDomainFolder("101621at6656", domainPath)
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
