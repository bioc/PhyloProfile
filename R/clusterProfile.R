#' Get data for calculating distance matrix from phylogenetic profiles
#' @export
#' @usage getDataClustering(data, profileType = "binary", var1AggBy = "max",
#'     var2AggBy = "max")
#' @param data a data frame contains processed and filtered profiles (see
#' ?fullProcessedProfile and ?filterProfileData, ?fromInputToProfile)
#' @param profileType type of data used for calculating the distance matrix.
#' Either "binary" (consider only the presence/absence status of orthlogs),
#' "orthoID" (consider ortholog IDs as values for clustering),
#' "var1"/"var2" for taking values of the additional variables into account.
#' Default = "binary".
#' @param var1AggBy aggregate method for VAR1 (min, max, mean or median).
#' Default = "max".
#' @param var2AggBy aggregate method for VAR2 (min, max, mean or median).
#' Default = "max".
#' @return A wide dataframe contains values for calculating distance matrix.
#' @author Carla Mölbert (carla.moelbert@gmx.de), Vinh Tran
#' (tran@bio.uni-frankfurt.de)
#' @importFrom data.table dcast setDT
#' @seealso \code{\link{fromInputToProfile}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' var1AggregateBy <- "max"
#' var2AggregateBy <- "mean"
#' getDataClustering(data, profileType, var1AggregateBy, var2AggregateBy)

getDataClustering <- function(
    data = NULL, profileType = "binary", var1AggBy = "max", var2AggBy = "max"
) {
    if (is.null(data)) stop("Input data cannot be NULL!")
    supertaxon <- presSpec <- orthoID <- NULL
    # remove lines where there is no found ortholog
    subDataHeat <- subset(data, data$presSpec > 0)
    # predict if ortho ID in BIONF format
    idFormat <- "other"
    firstOrtho <- strsplit(
        as.character(subDataHeat[1,]$orthoID), '|', fixed = TRUE)[[1]]
    if (
        length(firstOrtho) >= 3 && firstOrtho[1] == subDataHeat[1,]$geneID &&
        grepl(subDataHeat[1,]$supertaxonID, firstOrtho[2])
    ) idFormat <- "bionf"
    # transform data into wide matrix
    if (profileType == "binary") {
        subDataHeat <- subDataHeat[, c("geneID", "supertaxon", "presSpec")]
        subDataHeat$presSpec[subDataHeat$presSpec > 0] <- 1
        subDataHeat <- subDataHeat[!duplicated(subDataHeat), ]
        wideData <- data.table::dcast(
            data.table::setDT(subDataHeat),
            geneID ~ supertaxon, value.var = "presSpec")
    } else if(profileType == "orthoID") {
        subDataHeat <- subDataHeat[, c("geneID", "supertaxon", "orthoID")]
        if (idFormat == "bionf") {
            subDataHeat <- within(
                subDataHeat,
                orthoMod <- data.frame(
                    do.call(
                        'rbind', strsplit(as.character(orthoID),'|',fixed=TRUE)
                    )
                )
            )
            subDataHeat$orthoID <- paste(
                subDataHeat$orthoMod$X2, subDataHeat$orthoMod$X3, sep = "#"
            )
            subDataHeat <- subDataHeat[,!(names(subDataHeat) %in% ("orthoMod"))]
        } else {
            subDataHeat$orthoID <- paste(
                subDataHeat$supertaxon, subDataHeat$orthoID, sep = "#"
            )
        }
        subDataHeat <- subDataHeat[!duplicated(subDataHeat), ]
        l <- unique(subDataHeat$orthoID)
        subDataHeat$orthoID <- as.numeric(factor(subDataHeat$orthoID, levels=l))
        wideData <- data.table::dcast(
            data.table::setDT(subDataHeat),
            geneID ~ supertaxon, value.var = "orthoID", fun.aggregate = sum)
    } else {
        var <- profileType
        subDataHeat <- subDataHeat[, c("geneID", "supertaxon", var)]
        subDataHeat <- subDataHeat[!duplicated(subDataHeat), ]
        # aggreagte the values by the selected method
        if (var == "var1") aggregateBy <- var1AggBy
        else aggregateBy <- var2AggBy
        subDataHeat <- stats::aggregate(
            subDataHeat[, var],
            list(subDataHeat$geneID, subDataHeat$supertaxon),
            FUN = aggregateBy
        )
        colnames(subDataHeat) <- c("geneID", "supertaxon", var)
        wideData <- data.table::dcast(
            data.table::setDT(subDataHeat),
            geneID ~ supertaxon, value.var = var)
    }
    # set name for wide matrix as gene IDs
    wideData <- as.data.frame(wideData)
    dat <- wideData[, 2:ncol(wideData)]
    rownames(dat) <- wideData[, 1]
    dat[is.na(dat)] <- 0
    return(dat)
}

#' Calculate the distance matrix
#' @export
#' @param profiles dataframe contains profile data for distance calculating
#' (see ?getDataClustering)
#' @param method distance calculation method ("euclidean", "maximum",
#' "manhattan", "canberra", "binary", "distanceCorrelation",
#' "mutualInformation" or "pearson" for binary data; "distanceCorrelation" or
#' "mutualInformation" for non-binary data). Default = "mutualInformation".
#' @return A calculated distance matrix for input phylogenetic profiles.
#' @importFrom bioDist mutualInfo
#' @importFrom bioDist cor.dist
#' @importFrom stats dist as.dist
#' @importFrom energy dcor
#' @author Carla Mölbert (carla.moelbert@gmx.de), Vinh Tran
#' (tran@bio.uni-frankfurt.de)
#' @seealso \code{\link{getDataClustering}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' method <- "mutualInformation"
#' getDistanceMatrix(profiles, method)

getDistanceMatrix <- function(profiles = NULL, method = "mutualInformation") {
    if (is.null(profiles)) stop("Profile data cannot be NULL!")
    profiles <-  profiles[, colSums(profiles != 0) > 0]
    profiles <-  profiles[rowSums(profiles != 0) > 0, ]
    distMethods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
    if (method %in% distMethods) {
        distanceMatrix <- stats::dist(profiles, method = method)
    } else if (method == "distanceCorrelation") {
        n <- seq_len(nrow(profiles))
        matrix <- matrix(0L, nrow = nrow(profiles), ncol = nrow(profiles))
        for (i in seq_len(nrow(profiles))) { # rows
            p_i <- unlist(profiles[i,])
            for (j in seq_len(nrow(profiles))) { # columns
                if (i == j) break
                matrix[i, j] <- energy::dcor(p_i, unlist(profiles[j,]))
            }
        }
        # Swich the value so that the profiles with a high correlation
        # are clustered together
        matrix <- 1 - matrix
        matrix <- as.data.frame(matrix)
        profileNames <- rownames(profiles)
        colnames(matrix) <- profileNames[seq_len(length(profileNames)) - 1]
        rownames(matrix) <- profileNames
        distanceMatrix <- stats::as.dist(matrix)
    } else if (method == "mutualInformation") {
        distanceMatrix <- bioDist::mutualInfo(as.matrix(profiles))
        distanceMatrix <- max(distanceMatrix, na.rm = TRUE) - distanceMatrix
    } else if (method == "pearson") {
        distanceMatrix <-  bioDist::cor.dist(as.matrix(profiles))
    }
    return(distanceMatrix)
}

#' Create a hclust object from the distance matrix
#' @export
#' @param distanceMatrix calculated distance matrix (see ?getDistanceMatrix)
#' @param clusterMethod clustering method ("single", "complete",
#' "average" for UPGMA, "mcquitty" for WPGMA, "median" for WPGMC,
#' or "centroid" for UPGMC). Default = "complete".
#' @return An object class hclust generated based on input distance matrix and
#' a selected clustering method.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom stats hclust
#' @seealso \code{\link{getDataClustering}},
#' \code{\link{getDistanceMatrix}}, \code{\link{hclust}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' distMethod <- "mutualInformation"
#' distanceMatrix <- getDistanceMatrix(profiles, distMethod)
#' clusterMethod <- "complete"
#' clusterDataDend(distanceMatrix, clusterMethod)

clusterDataDend <- function(distanceMatrix = NULL, clusterMethod = "complete") {
    if (is.null(distanceMatrix)) stop("Distance matrix cannot be NULL!")
    dd.col <- stats::hclust(distanceMatrix, method = clusterMethod)
    return(dd.col)
}

#' Plot dendrogram tree
#' @export
#' @param dd dendrogram object (see ?clusterDataDend)
#' @return A dendrogram plot for the genes in the input phylogenetic profiles.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom ape plot.phylo
#' @seealso \code{\link{clusterDataDend}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' data <- finalProcessedProfile
#' profileType <- "binary"
#' profiles <- getDataClustering(
#'     data, profileType, var1AggregateBy, var2AggregateBy)
#' distMethod <- "mutualInformation"
#' distanceMatrix <- getDistanceMatrix(profiles, distMethod)
#' clusterMethod <- "complete"
#' dd <- clusterDataDend(distanceMatrix, clusterMethod)
#' getDendrogram(dd)

getDendrogram <- function(dd = NULL) {
    if (is.null(dd)) stop("Input dendrogram cannot be NULL!")
    p <- ape::plot.phylo(as.phylo(dd))
    return(p)
}
