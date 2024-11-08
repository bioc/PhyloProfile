% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umapClustering.R
\name{plotUmap3D}
\alias{plotUmap3D}
\title{Create UMAP cluster 3D plot}
\usage{
plotUmap3D(plotDf = NULL, legendPos = "bottom", 
    colorPalette = "Set2", transparent = 0,highlightTaxa = NULL, 
    dotZoom = 0)
}
\arguments{
\item{plotDf}{data for UMAP plot}

\item{legendPos}{position of legend. Default: "right"}

\item{colorPalette}{color palette. Default: "Set2"}

\item{transparent}{transparent level (from 0 to 1). Default: 0}

\item{highlightTaxa}{list of taxa to be highlighted}

\item{dotZoom}{dot size zooming factor. Default: 0}
}
\value{
A plot as ggplot object
}
\description{
Create UMAP cluster 3D plot
}
\examples{
rawInput <- system.file(
   "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
)
longDf <- createLongMatrix(rawInput)
umapData <- prepareUmapData(longDf, "phylum")
data.umap <- umapClustering3D(umapData)
plotDf <- createUmapPlotData(data.umap, umapData)
plotUmap3D(plotDf)
}
\seealso{
\code{\link{prepareUmapData}}, \code{\link{umapClustering}},
\code{\link{createUmapPlotData}}
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
