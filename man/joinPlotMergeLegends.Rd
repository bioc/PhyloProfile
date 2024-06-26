% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createDomainPlot.R
\name{joinPlotMergeLegends}
\alias{joinPlotMergeLegends}
\title{Join multiple plots and merge legends}
\usage{
joinPlotMergeLegends(
  df1 = NULL,
  df2 = NULL,
  plot1 = NULL,
  plot2 = NULL,
  position = c("bottom", "right"),
  font = "Arial"
)
}
\arguments{
\item{df1}{Data frame for plot 1}

\item{df2}{Data frame for plot 2}

\item{plot1}{ggplot object of plot 1}

\item{plot2}{ggplot object of plot 2}

\item{position}{position of legend (bottom or right)}

\item{font}{font of text}
}
\value{
joined plots with merged legend as a grid object
}
\description{
Join multiple plots and merge legends
}
\examples{
seed <- "101621at6656"
ortho <- "101621at6656|AGRPL@224129@0|224129_0:001955|1"
ortho <- gsub("\\\\|", ":", ortho)
grepID <- paste(seed, "#", ortho, sep = "")
domainFile <- system.file(
    "extdata", "domainFiles/101621at6656.domains",
    package = "PhyloProfile", mustWork = TRUE
)
domainDf <- parseDomainInput(seed, domainFile, "file")
domainDf$feature_id_mod <- domainDf$feature_id
subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
subdomainDf$feature <- as.character(subdomainDf$feature)
orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
minStart <- min(subdomainDf$start)
maxEnd <- max(c(subdomainDf$end, subdomainDf$length))
# resolve overlapping domains
seedDf <- PhyloProfile:::resolveOverlapFeatures(seedDf)
orthoDf <- PhyloProfile:::resolveOverlapFeatures(orthoDf)
# add feature colors
featureColorDf <- PhyloProfile:::addFeatureColors(seedDf, orthoDf)
seedDf <- featureColorDf[[1]]
orthoDf <- featureColorDf[[2]]
# generate plots
plotSeed <- PhyloProfile:::singleDomainPlotting(
    seedDf, seed, minStart = minStart, maxEnd = maxEnd, font = "sans"
)
plotOrtho <- PhyloProfile:::singleDomainPlotting(
    orthoDf, ortho, minStart = minStart, maxEnd = maxEnd, font = "sans"
)
# merge plots
PhyloProfile:::joinPlotMergeLegends(
    seedDf, orthoDf, plotSeed, plotOrtho, "bottom", font = "sans")
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
