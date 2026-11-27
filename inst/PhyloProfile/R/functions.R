#' Function to keep user defined geneID order
#' @param data data frame contains gene ID column
#' @param orderType either "none", "alphabetically", "by profile similarity" or
#' "by a sorted list" (from input$orderGenes)
#' @param geneOrder named list of gene IDs. Either "more" (input list has more
#' genes than the main input), "missing" (some gene IDs are missing from the
#' input list) or "sortedGenes" (genes should be ordered by this list)
#' @return data either sorted or non-sorted
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

sortGeneIDs <- function(data, orderType, geneOrder) {
    # Check if geneID is already a factor; if not, convert it
    if (!is.factor(data$geneID)) {
        data$geneID <- as.factor(data$geneID)
    }

    if (orderType == "none") {
        # Set factor levels to the current order in the data
        levels(data$geneID) <- unique(data$geneID)
    } else if (orderType == "user defined") {
        # Apply user-defined order only if `geneOrder` is valid
        if (!is.null(geneOrder) && "sortedGenes" %in% names(geneOrder)) {
            sortedGenes <- geneOrder$sortedGenes
            data$geneID <- factor(data$geneID, levels = sortedGenes)
        }
    }
    return(data)
}

#' Check installed packages
#' and install missing packages automatically
#' @param packages list of packages need to be checked
#' @return none
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

installPackages <- function(packages){
    missingPackages <-
        packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(missingPackages))
        install.packages(
            missingPackages,
            dependencies = TRUE,
            repos = "http://cran.us.r-project.org"
        )
}

installPackagesBioconductor <- function(packages){
    missingPackages <-
        packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(missingPackages)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install(missingPackages, ask = FALSE)
    }
}

#' Check internet connection
#' @return status of internet connection
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

hasInternet <- function(){
    !is.null(curl::nslookup("r-project.org", error = FALSE))
}

#' Create link to public database
#' @return link to public database
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

createDBlink <- function(id, source, type = "", version = ""){
    linkText <- ""
    url <- ""
    if (source == "NCBI") {
        url <- paste0("https://www.ncbi.nlm.nih.gov/protein/", id)
    } else if (source == "UniProt") {
        url <- paste0("https://www.uniprot.org/uniprot/", id)
    } else if (source == "OrthoDB") {
        if (version == "") {
            url <- paste0("https://www.orthodb.org/?query=", id)
        } else {
            version <- gsub("\\.", "-", version)
            url <- paste0("https://v", version, ".orthodb.org/?query=", id)
        }
        if (type == "gene") {
            idMod <- gsub(":", "%3A", id)
            if (version == "") {
                url <- paste0("https://www.orthodb.org/?gene=", idMod)
            } else {
                url <- paste0("https://v", version, ".orthodb.org/?gene=", idMod)
            }
        }
    } else if (source == "OMA") {
        url <- paste0("https://omabrowser.org/oma/omagroup/", id, "/members/")
        if (type == "gene") {
            url <- paste0("https://omabrowser.org/oma/info/", id)
        }
    }

    if (length(url > 0)) {
        if (RCurl::url.exists(url)) {
            linkText <- paste0(
                "<p><a href='", url, "' target='_blank'>",
                source, " entry for <strong>", id, "</strong></a></p>"
            )
        }
    }
    return(linkText)
}

#' Get the most specific taxonomy rank from input taxa
#' @return the most specific taxonomy rank
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getLowestRank <- function(longInput, taxDB) {
    inputTaxa <- gsub("ncbi","",levels(as.factor(longInput$ncbiID)))
    allDB <- getNameList(taxDB)
    allRanks <- getTaxonomyRanks()
    inputRanks <- levels(as.factor(allDB[allDB$ncbiID %in% inputTaxa,]$rank))
    if ("subspecies" %in% inputRanks) inputRanks <- c(inputRanks, "strain")
    inputRanks <- inputRanks[inputRanks %in% allRanks]
    return(allRanks[min(match(inputRanks, allRanks))][1])
}

#' Get colors for gene categoies
#' @param geneCategoryFile gene category file
#' @param type Either from manually uploaded file or predefined in config file
#' @return named vector of gene categories and their colors
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getCatColors <- function(geneCategoryFile, type = "file"){
    catColors <- NULL
    if (type == "file") {
        inputCatDt <- read.table(
            file = geneCategoryFile$datapath,
            sep = "\t",
            header = FALSE,
            check.names = FALSE,
            comment.char = "",
            fill = TRUE
        )

    } else if (type == "config"){
        inputCatDt <- read.table(
            file = geneCategoryFile,
            sep = "\t",
            header = FALSE,
            check.names = FALSE,
            comment.char = "",
            fill = TRUE
        )
    }
    if (ncol(inputCatDt) == 3) {
        catColorDf <- unique(inputCatDt[,c(2,3)])
        catColors <- as.character(catColorDf$V3)
        names(catColors) <- catColorDf$V2
    }
    return(catColors)
}

#' Get if ortho ID is in BIONF format (e.g. Q6PCB6|SACCE@4932@qfo|PROTID|1)
#' @param orthoID one ortholog ID
#' @param seedID seed ID of that ortholog
#' @param ncbiID ncbi ID of that ortholog
#' @return TRUE (if it is in BIONF format) or FALSE
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

checkBionfFormat <- function(orthoID, seedID, ncbiID) {
    ortho <- strsplit(as.character(orthoID),'|',fixed = TRUE)[[1]]
    if (length(ortho) >= 3 && ortho[1] == seedID && grepl(ncbiID, ortho[2]))
        return(TRUE)
    return(FALSE)
}

#' Adapt plot size based on number of genes and taxa
#' @param nrTaxa number of taxa
#' @param nrGene number of genes
#' @param xAxis type of x-axis, either "taxa" or "genes"
#' @param dotZoom zoom factor for dots
#' @return a vector contains number of itemes in y-axis, plot height and width
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

adaptPlotSize <- function(nrTaxa = 0, nrGene = 0, xAxis = "taxa", dotZoom = 0) {
    if (nrTaxa < 10000 && nrGene < 10000) {
        # adapte to axis type
        if (xAxis == "taxa") {
            h <- nrGene
            w <- nrTaxa
        } else {
            w <- nrGene
            h <- nrTaxa
        }
        # adapt to dot zoom factor
        if (dotZoom < -0.5){
            hv <- (200 + 12 * h) * (1 + dotZoom) + 500
            wv <- (200 + 12 * w) * (1 + dotZoom) + 500
        }  else if ((dotZoom < 0)) {
            hv <- (200 + 12 * h) * (1 + dotZoom) + 200
            wv <- (200 + 12 * w) * (1 + dotZoom) + 200
        } else {
            hv <- (200 + 12 * h) * (1 + dotZoom)
            wv <- (200 + 12 * w) * (1 + dotZoom)
        }
        # minimum size
        if (hv < 300) hv <- 300
        if (wv < 300) wv <- 300
        # update plot size based on number of genes/taxa
        hv <- hv + 300
        wv <- wv + 300
        return(c(h, hv, wv))
    } else return(c())
}

#' Estimate profile plot size (in px) from the data AND the space needed by
#' axis labels, supertaxon (division) labels and the legend.
#' @param nrTaxa number of (super)taxa
#' @param nrGene number of genes
#' @param xAxis "taxa" or "genes"
#' @param dotZoom dot zoom factor (-1 to 3)
#' @param taxaLabelLen max number of characters of taxon labels
#' @param geneLabelLen max number of characters of gene labels
#' @param xSize,ySize font size of x/y axis labels (pt)
#' @param xAngle angle of x-axis labels (degree)
#' @param groupLabelLen max number of characters of supertaxon labels
#'     (0 = no supertaxon labels)
#' @param groupLabelSize,groupLabelDist size/distance of supertaxon labels
#'     (as used in addRankDivisionPlot)
#' @param legendPos legend position ("right", "left", "top", "bottom", "none")
#' @param legendSize font size of the legend (pt)
#' @param fastMode TRUE for heatmapPlottingFast (no axis tick labels)
#' @return c(number of items on y-axis, height, width); unlike adaptPlotSize()
#'     this never returns an empty vector, so callers only need to cap it
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

estimatePlotSize <- function(
    nrTaxa = 0, nrGene = 0, xAxis = "taxa", dotZoom = 0,
    taxaLabelLen = 0, geneLabelLen = 0, xSize = 14, ySize = 14, xAngle = 60,
    groupLabelLen = 0, groupLabelSize = 7, groupLabelDist = 1.5,
    legendPos = "right", legendSize = 12, fastMode = FALSE
) {
    orDefault <- function(x, d) {
        if (is.null(x) || length(x) == 0 || is.na(x[1])) d else x[1]
    }
    dotZoom <- orDefault(dotZoom, 0)
    taxaLabelLen <- orDefault(taxaLabelLen, 0)
    geneLabelLen <- orDefault(geneLabelLen, 0)
    xSize <- orDefault(xSize, 14); ySize <- orDefault(ySize, 14)
    xAngle <- orDefault(xAngle, 60)
    groupLabelLen <- orDefault(groupLabelLen, 0)
    groupLabelSize <- orDefault(groupLabelSize, 7)
    groupLabelDist <- orDefault(groupLabelDist, 1.5)
    legendPos <- orDefault(legendPos, "right")
    legendSize <- orDefault(legendSize, 12)
    charW <- 0.55 # average character width relative to font size
    # items and label lengths on each axis
    if (xAxis == "taxa") {
        nY <- nrGene; nX <- nrTaxa
        xLabLen <- taxaLabelLen; yLabLen <- geneLabelLen
    } else {
        nY <- nrTaxa; nX <- nrGene
        xLabLen <- geneLabelLen; yLabLen <- taxaLabelLen
    }
    # panel: 12 px per gene/taxon at dotZoom = 0, plus ~1 cell of padding
    cell <- 12 * max(1 + dotZoom, 0.25)
    panelH <- max(100, (nY + 1.2) * cell)
    panelW <- max(100, (nX + 1.2) * cell)
    # axis titles + outer plot margin
    top <- 10; right <- 10
    bottom <- 10 + 1.6 * xSize
    left <- 10 + 1.6 * ySize
    # axis tick labels (hidden in fast mode)
    if (!fastMode) {
        a <- xAngle * pi / 180
        xLabW <- xLabLen * xSize * charW
        bottom <- bottom + xLabW * abs(sin(a)) + xSize * abs(cos(a)) + 5
        left <- left + yLabLen * ySize * charW + 5
        # first rotated x label sticks out to the left of the panel
        left <- max(left, xLabW * abs(cos(a)) + 10)
    }
    # supertaxon labels: same margin as set in addRankDivisionPlot()
    if (groupLabelLen > 0) {
        groupMargin <- (groupLabelLen + 1) * groupLabelSize * groupLabelDist
        if (xAxis == "taxa") top <- max(top, groupMargin)
        else {
            right <- max(right, groupMargin)
            legendPos <- "top" # forced by addRankDivisionPlot()
        }
    }
    # legend
    if (legendPos %in% c("right", "left")) {
        right <- right + legendSize * 9 + 30
    } else if (legendPos %in% c("top", "bottom")) {
        top <- top + legendSize * 4 + 20
    }
    hv <- ceiling(panelH + top + bottom)
    wv <- ceiling(panelW + left + right)
    return(c(nY, hv, wv))
}

#' Get max length of supertaxon (division) labels
#' @param taxonIDs NCBI IDs of the (super)taxa at the working rank
#' @param workingRank working taxonomy rank (e.g. species)
#' @param superRank rank used for the division lines (e.g. phylum)
#' @param taxDB path to taxonomy database
#' @return max number of characters of the labels (0 if no superRank), incl.
#'     room for the "*" / ",..." marks added by addRankDivisionPlot()
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getSuperRankLabelLen <- function(
    taxonIDs = NULL, workingRank = NULL, superRank = NULL, taxDB = NULL
) {
    if (is.null(superRank) || superRank == "" || is.null(workingRank) ||
        length(taxonIDs) == 0) return(0)
    tryCatch({
        taxMatrix <- getTaxonomyMatrix(taxDB)
        superIDs <- unique(
            taxMatrix[taxMatrix[[workingRank]] %in% taxonIDs, superRank]
        )
        nameList <- getNameList(taxDB)
        labels <- nameList$fullName[nameList$ncbiID %in% superIDs]
        if (length(labels) == 0) return(15)
        max(nchar(as.character(labels))) + 4
    }, error = function(e) 15)
}

#' Scale a list of values into 0-1
#' @param x a numeric vector
#' @return a vector with values scaled from 0 to 1

scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

#' Helper function to read single-column files
#' @param fileInput fileInput object
#' @return vector of unique values in input file

readSingleColFile <- function(fileInput) {
    if (!is.null(fileInput)) {
        labels <- readLines(fileInput$datapath)
        return(unique(labels))
    }
    return(NULL)
}

idWithSuffix <- function(base, suffix) {
    if (suffix == "") base else paste0(base, suffix)
}


# FUNCTIONS FOR RENDER UI ELEMENTS ============================================
createSliderCutoff <- function(id, title, start, stop, varID){
    if (is.null(varID)) return()
    if (varID == "") {
        sliderInput(id, title,
                    min = 1,
                    max = 1,
                    step = 0.025,
                    value = 1,
                    width = 200)
    } else {
        sliderInput(id, title,
                    min = 0,
                    max = 1,
                    step = 0.025,
                    value = c(start, stop),
                    width = 200)
    }
}

updateSliderCutoff <- function(session, id, title, newVar, varID){
    if (is.null(varID) || varID == "") return()
    updateSliderInput(session, id, title,
                      value = newVar,
                      min = 0,
                      max = 1,
                      step = 0.025)
}

createPlotSize <- function(id, title, value, width = 100) {
    numericInput(id,
                 title,
                 min = 100,
                 max = 3200,
                 step = 50,
                 value = value,
                 width = width)
}

createTextSize <- function(id, title, value, width = 100) {
    numericInput(id,
                 title,
                 min = 3,
                 max = 99,
                 step = 1,
                 value = value,
                 width = width)
}

#' Replace ~ symbol by the name of home folder
#' @export
#' @param fullPath complete path to input file or folder
#' @return complete path of input file or folder without ~ symbol
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @import stringr
#' @examples
#' replaceHomeCharacter("~/path/to/something")
replaceHomeCharacter <- function (fullPath = NULL) {
    if (!(Sys.info()['sysname'] == "Windows")){
        homeName <- system("echo $HOME", intern = TRUE)
        return(stringr::str_replace(fullPath, "~", homeName))
    } else return(fullPath)
}

substrRight <- function(x, n) {
    substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n) {
    substr(x, 1, n)
}
