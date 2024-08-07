% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umapClustering.R
\name{createUmapPlotData}
\alias{createUmapPlotData}
\title{Create UMAP cluster plot}
\usage{
createUmapPlotData(umapData = NULL, data4umap = NULL, labelNr = 5, 
    excludeTaxa = "None")
}
\arguments{
\item{umapData}{data contains UMAP cluster (output from umapClustering())}

\item{data4umap}{data for UMAP clustering (output from prepareUmapData())}

\item{labelNr}{maximal number of labels. Default: 5}

\item{excludeTaxa}{hide taxa from plot. Default: "None"}
}
\value{
A plot as ggplot object
}
\description{
Create UMAP cluster plot
}
\examples{
rawInput <- system.file(
   "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
)
longDf <- createLongMatrix(rawInput)
data4umap <- prepareUmapData(longDf, "phylum")
umapData <- umapClustering(data4umap)
createUmapPlotData(umapData, data4umap)
}
\seealso{
\code{\link{prepareUmapData}}, \code{\link{umapClustering}}
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
