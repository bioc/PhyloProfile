# Get OMA information functions
# using OmaDB package (https://github.com/DessimozLab/OmaDB)

#' Check the validity of input OMA IDs
#' @description Check if input IDs are valid OMA IDs for OMA Browser
#' @export
#' @param ids list of ids needs to be checked
#' @return List of invalid IDs (not readable for OMA)
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' ### Uncomment the following line to run the function
#' # checkOmaID("HUMAN29398")

checkOmaID <- function(ids) {
    if (!requireNamespace("OmaDB", quietly = TRUE)) {
        stop(
            "Package \"OmaDB\" must be installed to use this function.",
            call. = FALSE
        )
    }
    ids <- as.character(ids)
    getOma <- OmaDB::getProtein(ids)
    return(setdiff(ids, names(getOma)))
}

#' Get OMA members
#' @description Get OMA ortholog group, OMA HOG or OMA pair's members for a
#' seed protein from OMA Browser.
#' @export
#' @param id ID of the seed protein (OMA or UniProt ID)
#' @param orthoType type of OMA orthologs: either "HOG", "OG"
#' (orthologous group) or "PAIR" (orthologous pair - CURRENTLY NOT WORKING).
#' Default = "OG".
#' @return List of OMA orthologs for an input seed protein.
#' @author Carla Mölbert carla.moelbert@gmx.de
#' @examples
#' ### Uncomment the following line to run the function
#' # getOmaMembers("HUMAN29397", "OG")

getOmaMembers <- function(id = NULL, orthoType = "OG") {
    if (!requireNamespace("OmaDB", quietly = TRUE)) {
        stop(
            "Package \"OmaDB\" must be installed to use this function.",
            call. = FALSE
        )
    }
    if (is.null(id)) stop("No OMA ID given!")
    # get the members of the Hierarchical Orthologous Group
    if (orthoType == "HOG") {
        members <- suppressWarnings(
            OmaDB::getHOG(id = id, level = "root", members = TRUE)$members$omaid
        )
    }
    # get the members of the Ortholoug group
    else if (orthoType == "OG") {
        members <- suppressWarnings(
            OmaDB::getOMAGroup(id = id)$members$omaid
        )
    }
    return(members)
}

#' Get domain annotation from OMA Browser
#' @export
#' @description Get domain annotation from OMA Browser based on a URL or a
#' raw data frame contains annotation info from OMA
#' @param domainURL URL address for domain annotation of ONE OMA id or a
#' raw data frame contains annotation info from OMA
#' @return Data frame contains feature names with their start and end positions
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' ### Uncomment the following line to run the function
#' # getOmaDomainFromURL("https://omabrowser.org/api/protein/7916808/domains/")

getOmaDomainFromURL <- function(domainURL = NULL) {
    if (!requireNamespace("OmaDB", quietly = TRUE)) {
        stop(
            "Package \"OmaDB\" must be installed to use this function.",
            call. = FALSE
        )
    }
    if (is.null(domainURL)) stop("OMA domain URL is NULL!")
    if (grepl("https://", domainURL[1])) {
        domains <- OmaDB::resolveURL(domainURL)$regions
    } else {
        domains <- domainURL$regions
    }

    pos <- strsplit(domains$location, ":")
    allPos <- unlist(pos)
    start <- allPos[c(TRUE, FALSE)]
    end <- allPos[c(FALSE, TRUE)]

    feature <- ifelse(
        nchar(domains$name) > 0,
        paste0(domains$source, "_", domains$domainid, " (",domains$name,")"),
        paste0(domains$source, "_", domains$domainid)
    )
    feature <- gsub("#", "-", feature)
    featureRep <- vapply(pos, FUN.VALUE = numeric(1), function (x) length(x)/2)

    return(data.frame(
        feature = rep(feature, featureRep), start, end, stringsAsFactors = FALSE
    ))
}

#' Get taxonomy ID, sequence and annotation for one OMA protein
#' @export
#' @param id oma ID of one protein
#' @return Data frame contains the input protein ID with its taxonomy ID,
#' sequence, length and domain annotations (tab delimited) for input OMA protein
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' ### Uncomment the following line to run the function
#' # getOmaDataForOneOrtholog("HUMAN29397")

getOmaDataForOneOrtholog <- function(id = NULL) {
    if (!requireNamespace("OmaDB", quietly = TRUE)) {
        stop(
            "Package \"OmaDB\" must be installed to use this function.",
            call. = FALSE
        )
    }
    if (is.null(id)) stop("No OMA ID given!")
    # get ncbi taxonomy id
    specName <- substr(id, 1, 5)
    taxonID <- paste0(
        "ncbi", OmaDB::getTaxonomy(members = specName, newick = FALSE)$id
    )
    # get raw data
    raw <- OmaDB::getProtein(id)

    # get sequence
    seq <- as.character(raw$sequence)

    # get sequence length
    length <- raw$sequence_length

    # get annotation
    rawDomains <- suppressWarnings(raw$domains)
    domainDf <- suppressWarnings(getOmaDomainFromURL(rawDomains))
    domainDfJoin <- c(domainDf, sep = "#")
    domains <- paste(unlist(do.call(paste, domainDfJoin)), collapse = "\t")

    # return data frame contains all info
    omaDf <- data.frame(
        orthoID = id,
        taxonID = taxonID,
        seq = seq,
        length = length,
        domains = domains,
        stringsAsFactors = FALSE
    )
    return(omaDf)
}

#' Get OMA info for a query protein and its orthologs
#' @description Get taxonomy IDs, sequences, length and annotations for an OMA
#' orthologous group (or OMA HOG).
#' @export
#' @param seedID OMA protein ID
#' @param orthoType type of OMA orthologs ("OG" or "HOG"). Default = "OG".
#' @return Data frame contains info for all sequences of the input OMA group
#' (or HOG). That info contains the protein IDs, taxonomy IDs, sequences,
#' lengths, domain annotations (tab delimited) and the corresponding seed ID.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' ### Uncomment the following line to run the function
#' # getDataForOneOma("HUMAN29397", "OG")

getDataForOneOma <- function(seedID = NULL, orthoType = "OG"){
    if (!requireNamespace("OmaDB", quietly = TRUE)) {
        stop(
            "Package \"OmaDB\" must be installed to use this function.",
            call. = FALSE
        )
    }
    if (is.null(seedID)) stop("No OMA seed ID given!")
    # get members
    idList <- getOmaMembers(seedID, orthoType)
    specName <- substr(idList, 1, 5)

    # get taxonomy IDs
    taxonID <- paste0(
        "ncbi",
        lapply(
            specName,
            function (x) OmaDB::getTaxonomy(members = x, newick = FALSE)$id
        )
    )

    # get sequences, protein lengths and their domain annotations
    raw <- OmaDB::getProtein(idList)
    seq <- lapply(
        seq(length(raw)), function (x) as.character(raw[[x]]$sequence)
    )
    length <- lapply(
        seq(length(raw)), function (x) raw[[x]]$sequence_length
    )
    domains <- lapply(
        seq(length(raw)),
        function (x) {
            rawDomains <- suppressWarnings(raw[[x]]$domains)
            domainDf <- suppressWarnings(getOmaDomainFromURL(rawDomains))
            domainDfJoin <- c(domainDf, sep = "#")
            paste(unlist(do.call(paste, domainDfJoin)), collapse = "\t")
        }
    )
    # return data frame
    return(data.frame(
        orthoID = idList,
        taxonID = taxonID,
        seq = unlist(seq),
        length = unlist(length),
        domains = unlist(domains),
        seed = rep(seedID, length(idList)),
        stringsAsFactors = FALSE
    ))
}

#' Create a phylogenetic profile from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of proteins (see ?getDataForOneOma)
#' @return Dataframe of the phylogenetic profiles in long format, which
#' contains the seed protein IDs, their orthologous proteins and the
#' corresponding taxononmy IDs of the orthologs.
#' @seealso \code{\link{getDataForOneOma}}
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' ### Uncomment the following lines to run the function
#' # omaData <- getDataForOneOma("HUMAN29397", "OG")
#' # createProfileFromOma(omaData)

createProfileFromOma <- function(finalOmaDf = NULL) {
    if (is.null(finalOmaDf)) stop("Raw OMA data cannot be NULL!")
    profileDf <- finalOmaDf[, c("seed", "taxonID", "orthoID")]
    colnames(profileDf) <- c("geneID", "ncbiID", "orthoID")
    return(profileDf[!duplicated(profileDf), ])
}

#' Create domain annotation dataframe from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of proteins (see ?getDataForOneOma)
#' @return Dataframe of the domain annotation used for PhyloProfile, which
#' contains seed IDs, ortholog IDs, ortholog lengths, annotated features, start
#' and end positions of those features.
#' @seealso \code{\link{getDataForOneOma}}
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table rbindlist
#' @examples
#' ### Uncomment the following line to run the function
#' # omaData <- getDataForOneOma("HUMAN29397", "OG")
#' # getAllDomainsOma(omaData)

getAllDomainsOma <- function(finalOmaDf = NULL) {
    if (is.null(finalOmaDf)) stop("Raw OMA data cannot be NULL!")

    seedID <- finalOmaDf$seed
    orthoID <- finalOmaDf$orthoID
    length <- finalOmaDf$length
    feature <- finalOmaDf$domains

    featureCount <- vapply(
        feature,
        function (x) length(unlist(strsplit(x, "\t"))),
        FUN.VALUE = numeric(1)
    )
    domainDf <- data.frame(
        do.call(rbind, strsplit(unlist(strsplit(feature, "\t")), "#")),
        stringsAsFactors = FALSE
    )
    colnames(domainDf) <- c("feature", "start", "end")
    domainDf$seedID <- rep(paste0(seedID, "#", orthoID), featureCount)
    domainDf$orthoID <- rep(orthoID, featureCount)
    domainDf$length <- rep(length, featureCount)

    seedDf <- domainDf[domainDf$orthoID %in% levels(as.factor(seedID)),]
    seedDfFull <- data.frame(
        feature = rep(seedDf$feature, nlevels(as.factor(domainDf$seedID))),
        start = rep(seedDf$start, nlevels(as.factor(domainDf$seedID))),
        end = rep(seedDf$end, nlevels(as.factor(domainDf$seedID))),
        seedID = rep(levels(as.factor(domainDf$seedID)), each = nrow(seedDf)),
        orthoID = rep(seedDf$orthoID, nlevels(as.factor(domainDf$seedID))),
        length = rep(seedDf$length, nlevels(as.factor(domainDf$seedID))),
        stringsAsFactors = FALSE
    )

    outDf <- data.table::rbindlist(list(domainDf, seedDfFull))
    outDf$length <- as.numeric(outDf$length)
    outDf$start <- as.numeric(outDf$start)
    outDf$end <- as.numeric(outDf$end)

    return(
        outDf[, c("seedID", "orthoID", "length", "feature", "start", "end")]
    )
}

#' Get all fasta sequences from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of proteins (see ?getDataForOneOma)
#' @return A list contains all protein sequences in fasta format.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{getDataForOneOma}}
#' @examples
#' ### Uncomment the following line to run the function
#' # omaData <- getDataForOneOma("HUMAN29397", "OG")
#' # getAllFastaOma(omaData)

getAllFastaOma <- function(finalOmaDf = NULL) {
    if (is.null(finalOmaDf)) stop("Raw OMA data cannot be NULL!")
    fastaDf <- finalOmaDf[, c("orthoID", "seq")]
    fastaOut <- paste(paste0(">", fastaDf$orthoID), fastaDf$seq, sep = "\n")
    return(unique(fastaOut))
}

#' Get selected fasta sequences from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of proteins (see ?getDataForOneOma)
#' @param seqID OMA ID of selected protein
#' @return Required protein sequence in fasta format.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{getDataForOneOma}}
#' @examples
#' ### Uncomment the following line to run the function
#' # omaData <- getDataForOneOma("HUMAN29397", "OG")
#' # getSelectedFastaOma(omaData, "HUMAN29397")

getSelectedFastaOma <- function(finalOmaDf = NULL, seqID = NULL) {
    if (is.null(finalOmaDf)) stop("Raw OMA data cannot be NULL!")
    selectedDf <- subset(
        finalOmaDf[, c("orthoID", "seq")],
        finalOmaDf$orthoID == seqID
    )
    header <- paste0(">", selectedDf$orthoID[1])
    return(paste(header, selectedDf$seq[1], sep = "\n"))
}
