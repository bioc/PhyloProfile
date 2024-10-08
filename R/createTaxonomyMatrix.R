#' Get all NCBI taxonomy rank names
#' @export
#' @return A list of all available NCBI taxonomy rank names.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' mainTaxonomyRank()

mainTaxonomyRank <- function() {
    return(
        c(
            "isolate","strain","pathogroup","serotype","serogroup",
            "forma","formaspecialis","varietas",
            "genotype","morph","subvariety","biotype",
            "subspecies","species","speciessubgroup",
            "speciesgroup","series","section","subgenus","genus",
            "subtribe","tribe", "subfamily","family","superfamily",
            "parvorder","infraorder","suborder","order","superorder",
            "subcohort","cohort","infraclass","subclass","class","superclass",
            "subphylum","phylum","superphylum",
            "subkingdom","kingdom","superkingdom"
        )
    )
}

#' Pre-processing NCBI taxonomy data
#' @description Download NCBI taxonomy database and parse information that are
#' needed for PhyloProfile, including taxon IDs, their scientific names,
#' systematic ranks, and parent (next higher) rank IDs.
#' @return A dataframe contains NCBI taxon IDs, taxon names, taxon ranks and the
#' next higher taxon IDs (parent's IDs) of all taxa in the NCBI taxonomy
#' database.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @export
#' @examples
#' \dontrun{
#' ?processNcbiTaxonomy
#' preProcessedTaxonomy <- PhyloProfile:::processNcbiTaxonomy()
#' # save to text (tab-delimited) file
#' write.table(
#'     preProcessedTaxonomy,
#'     file = "preProcessedTaxonomy.txt",
#'     col.names = TRUE,
#'     row.names = FALSE,
#'     quote = FALSE,
#'     sep = "\t"
#' )
#' # save to rdata file
#' save(
#'     preProcessedTaxonomy, file = "preProcessedTaxonomy.RData", compress='xz'
#' )
#' }

processNcbiTaxonomy <- function() {
    temp <- tempfile()
    utils::download.file(
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip", temp
    )
    names <- utils::read.table(
        unz(temp, "names.dmp"), header = FALSE, fill = TRUE, sep = "\t",
        quote = "", comment.char = "", stringsAsFactors = FALSE
    )
    nodes <- utils::read.table(
        unz(temp, "nodes.dmp"), header = FALSE, fill = TRUE, sep = "\t",
        quote = "", comment.char = "", stringsAsFactors = FALSE
    )
    unlink(temp)
    message("Download NCBI taxonomy done!")
    # Create data frame containing taxon ID, the scientific name, its taxonomy
    # rank and the taxon ID of the higer rank (parent's ID)
    outTax <- merge(
        names[names$V7 == "scientific name", c("V1", "V3")],
        nodes[,c("V1", "V5", "V3")],
        by = "V1"
    )
    colnames(outTax) <- c("ncbiID", "fullName", "rank","parentID")
    # Remove "'" from taxon names and ranks, remove space from taxon ranks
    outTax$fullName <- gsub("'", "",outTax$fullName)
    outTax$rank <- gsub("'", "", outTax$rank)
    outTax$rank <- gsub(" ", "", outTax$rank)
    # Add Holozoa and Holomycota ranks
    opiMod <- data.frame(
        ncbiID = c(999999001, 999999002, 999999003, 999999004, 999999005),
        fullName = c(
            "Choanozoa", "Filozoa", "Pluriformea", "Holozoa", "Holomycota"),
        rank = rep("norank", 5),
        parentID = c(999999002, 999999003, 999999004, 33154, 33154),
        stringsAsFactors = FALSE
    )
    outTax[outTax$ncbiID %in% c(28009, 33208),]$parentID <- 999999001
    outTax[outTax$ncbiID == 2687318,]$parentID <- 999999002
    outTax[outTax$ncbiID == 42461,]$parentID <- 999999003
    outTax[outTax$ncbiID == 127916,]$parentID <- 999999004
    outTax[outTax$ncbiID %in% c(4751, 2686024),]$parentID <- 999999005
    outTax <- do.call(rbind, list(outTax, opiMod))
    # Return
    message("Parsing NCBI taxonomy done!")
    return(outTax)
}

#' Get taxonomy info for a list of input taxa
#' @param inputTaxa NCBI taxonomy IDs of input taxa.
#' @param currentNCBIinfo table/dataframe of the pre-processed NCBI taxonomy
#' data (/PhyloProfile/data/preProcessedTaxonomy.txt)
#' @return A list of NCBI taxonomy info for input taxa, including the taxonomy
#' IDs, full scientific names, taxonomy ranks and the parent IDs.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table fread
#' @export
#' @examples
#' inputTaxa <- c("272557", "176299")
#' ncbiFilein <- system.file(
#'     "extdata", "data/preProcessedTaxonomy.txt",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' currentNCBIinfo <- as.data.frame(data.table::fread(ncbiFilein))
#' getTaxonomyInfo(inputTaxa, currentNCBIinfo)
getTaxonomyInfo <- function(inputTaxa = NULL, currentNCBIinfo = NULL) {
    tmp <- list()
    outList <- list()
    k <- 1
    for (refID in inputTaxa) {
        # get info for this taxon
        refEntry <- currentNCBIinfo[currentNCBIinfo$ncbiID == refID, ]
        lastID <- refEntry$parentID
        inputTaxaInfo <- refEntry
        if(length(lastID) == 0) {
            msg <- paste("PROBLEM with", refID)
            stop(msg)
        }
        while (lastID != 1) {
            if (lastID %in% names(tmp)) {
                inputTaxaInfo <- rbindlist(
                    list(inputTaxaInfo, tmp[[toString(lastID)]]),
                    use.names = TRUE, fill = TRUE, idcol = NULL)
                lastID <- 1
            } else {
                nextEntry <- currentNCBIinfo[currentNCBIinfo$ncbiID == lastID, ]
                inputTaxaInfo <- rbindlist(
                    list(inputTaxaInfo, nextEntry), use.names = TRUE,
                    fill = TRUE, idcol = NULL)
                lastID <- nextEntry$parentID
            }
        }
        for (i in seq(inputTaxaInfo$ncbiID)) {
            if (!(inputTaxaInfo$ncbiID[i] %in% names(tmp))) {
                tmp[[toString(inputTaxaInfo$ncbiID[i])]] <-
                    inputTaxaInfo[i:length(inputTaxaInfo$ncbiID)]
            }
        }
        outList[[k]] <- inputTaxaInfo
        k <- k + 1
    }
    return(outList)
}

#' Get taxonomy info for a list of taxa
#' @description Get NCBI taxonomy IDs, ranks and names for an input taxon list.
#' @param inputTaxa NCBI ID list of input taxa.
#' @param currentNCBIinfo table/dataframe of the pre-processed NCBI taxonomy
#' data (/PhyloProfile/data/preProcessedTaxonomy.txt)
#' @return A list of 3 dataframes: idList, rankList and reducedInfoList. The
#' "rankList" contains taxon names and all taxonomy ranks of the input taxa
#' including also the noranks from the input rank to the taxonomy root. The
#' "idList" contains input taxon IDs, taxon names, all the ranks from current
#' rank to the taxonomy root together with their IDs (with the format
#' "id#rank"). The reducedInfoList is a subset of preProcessedTaxonomy.txt file,
#' containing the NCBI IDs, taxon fullnames, their current rank and their
#' direct parent ID.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table fread rbindlist
#' @export
#' @examples
#' inputTaxa <- c("272557", "176299")
#' ncbiFilein <- system.file(
#'     "extdata", "data/preProcessedTaxonomy.txt",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' currentNCBIinfo <- as.data.frame(data.table::fread(ncbiFilein))
#' getIDsRank(inputTaxa, currentNCBIinfo)

getIDsRank <- function(inputTaxa = NULL, currentNCBIinfo = NULL){
    if (is.null(currentNCBIinfo)) stop("Pre-processed NCBI tax data is NULL!")
    inputTaxaInfo <- getTaxonomyInfo(inputTaxa, currentNCBIinfo)
    allMainRank <- mainTaxonomyRank()
    ## get reduced taxonomy info (subset of preProcessedTaxonomy.txt)
    reducedDf <- unique(rbindlist(inputTaxaInfo))
    ## get list of all ranks and rank#IDs
    rankMod <- ncbiID <- NULL
    inputRankIDDf <- lapply(
        seq_len(length(inputTaxaInfo)),
        function (x) {
            inputTaxaInfo[[x]]$rank[
                !(inputTaxaInfo[[x]]$rank %in% allMainRank)] <- "norank"
            inputTaxaInfo[[x]]$rankMod <- inputTaxaInfo[[x]]$rank
            if (inputTaxaInfo[[x]]$rank[1] == "norank")
                inputTaxaInfo[[x]]$rankMod[1] <-
                    paste0("strain_", inputTaxaInfo[[x]]$ncbiID[1])
            inputTaxaInfo[[x]]$rankMod <- with(
                inputTaxaInfo[[x]],
                ifelse(rankMod == "norank", paste0("norank_", ncbiID), rankMod))
            inputTaxaInfo[[x]]$id <- paste0(
                inputTaxaInfo[[x]]$ncbiID, "#", inputTaxaInfo[[x]]$rankMod)
            return(inputTaxaInfo[[x]][ , c("rankMod", "id")])
        }
    )
    inputRankList <- lapply(
        seq_len(length(inputRankIDDf)),
        function (x) {
            ll <- c(
                paste0("ncbi", inputTaxa[x]),
                reducedDf$fullName[reducedDf$ncbiID == inputTaxa[x]],
                inputRankIDDf[[x]]$rank, "norank_1")
            ll <- gsub("strain_[[:digit:]]+", "strain", ll)
            return(data.frame(
                matrix(ll, nrow = 1, byrow = TRUE), stringsAsFactors = FALSE))
        }
    )
    inputRankDf <- data.frame(data.table::rbindlist(inputRankList, fill = TRUE))
    inputIDList <- lapply(
        seq_len(length(inputRankIDDf)),
        function (x) {
            ll <- c(paste0(
                reducedDf$fullName[reducedDf$ncbiID==inputTaxa[x]],"#name"),
                inputRankIDDf[[x]]$id, "1#norank_1")
            return(data.frame(
                matrix(ll, nrow = 1, byrow = TRUE), stringsAsFactors = FALSE))
        }
    )
    inputIDDf <- data.frame(data.table::rbindlist(inputIDList, fill = TRUE))
    newCol <- seq(ncol(inputIDDf) + 1, ncol(inputRankDf))
    inputIDDf[paste0("X", newCol)] <- NA
    reducedDf$rank[!(reducedDf$rank %in% allMainRank)] <- "norank"
    return(list(inputIDDf, inputRankDf, as.data.frame(reducedDf)))
}

#' Indexing all available ranks (including norank)
#' @param rankListFile Input file, where each row is a rank list of a taxon
#' (see rankListFile in example)
#' @return A dataframe containing a list of all possible ranks and their indexed
#' values.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' rankListFile <- system.file(
#'     "extdata", "data/rankList.txt", package = "PhyloProfile", mustWork = TRUE
#' )
#' PhyloProfile:::rankIndexing(rankListFile)

rankIndexing <- function (rankListFile = NULL) {
    if (is.null(rankListFile)) stop("Rank list file is NULL!")
    . <- NULL
    rankList <- utils::read.table(
        rankListFile, sep = '\t', header = FALSE,fill = TRUE,
        stringsAsFactors = TRUE, na.strings = c("", "NA")
    )
    ### remove duplicated and "normal" taxa
    rankList <- rankList[,-(seq_len(2))]
    rankList <- rankList[!duplicated(rankList), ]
    mainRank <- mainTaxonomyRank()
    filterRanks <- c(mainRank, "norank_1", NA)
    rankList <- rankList %>% filter_all(any_vars(!(. %in% filterRanks)))
    ### get list of all inpu ranks
    allInputRank <- as.character(unique(unlist(rankList)))
    allInputRank <- allInputRank[!is.na(allInputRank)]
    ### initial index for main ranks
    mainRank <- mainTaxonomyRank()
    rank2index <- new.env(hash = TRUE)
    getHash <- Vectorize(get, vectorize.args = "x")
    assignHash <- Vectorize(assign, vectorize.args = c("x", "value"))
    for (i in seq_len(length(mainRank))) rank2index[[mainRank[i]]] <- i
    ### the magic happens here
    for (k in seq_len(nrow(rankList))) {
        subList <- rankList[k,][!is.na(rankList[k,])]
        ## indexing
        tmpEnv <- new.env(hash = TRUE)
        flag <- 0
        for (i in seq_len(length(subList))) {
            iRank <- subList[i]
            if (is.null(rank2index[[iRank]])) {
                for (j in seq_len(length(subList))) {
                    if (j < i) {
                        if (!is.null(tmpEnv[[subList[i - j]]])) {
                            tmpEnv[[iRank]] <- tmpEnv[[subList[i - j]]] + 1
                            break
                        }
                    } else j <- j - 1
                }
            } else {
                # for old rank
                if (i > 1) {
                    if (flag == 0) {
                        if (!(iRank %in% ls(rank2index)))
                            stop(iRank," not found!")
                        currentIndex <- rank2index[[iRank]]
                    } else {
                        if (!(iRank %in% ls(tmpEnv))) {
                            tmpEnv[[iRank]] <- tmpEnv[[subList[i-1]]]
                        }
                        currentIndex <- tmpEnv[[iRank]]
                    }

                    if (currentIndex <= tmpEnv[[subList[i-1]]]) {
                        if (flag == 0) {
                            tmpEnv[[iRank]] <- tmpEnv[[subList[i-1]]] + 1
                            # list of curr ranks whose index should be increased
                            candidateList <- unlist(
                                mget(
                                    ls(rank2index)[
                                        !(ls(rank2index) %in% ls(tmpEnv))
                                    ],
                                    rank2index
                                )
                            )
                            candidateList <- candidateList[
                                order(unlist(candidateList))
                            ]
                            for (c in seq_len(length(candidateList))) {
                                r <- names(candidateList)[c]
                                fromIndex <- rank2index[[iRank]]
                                if(subList[i-1] %in% ls(rank2index)) {
                                    fromIndex <- rank2index[[subList[i-1]]]
                                }

                                if (rank2index[[r]] > fromIndex) {
                                    tmpEnv[[r]] <-
                                        rank2index[[r]] +
                                        (tmpEnv[[iRank]] - rank2index[[iRank]])
                                    flag <- 1
                                }
                            }
                        } else {
                            step <- tmpEnv[[subList[i-1]]] -
                                rank2index[[iRank]] + 1
                            tmpEnv[[iRank]] <- tmpEnv[[subList[i-1]]] + 1
                            for (t in ls(tmpEnv)) {
                                if (tmpEnv[[t]] >= tmpEnv[[subList[i-1]]]){
                                    if (!(t == iRank) & !(t == subList[i-1])){
                                        tmpEnv[[t]] <- tmpEnv[[t]] + step
                                    }
                                }
                            }
                        }
                        assignHash(
                            ls(tmpEnv), getHash(ls(tmpEnv), tmpEnv), rank2index
                        )
                    } else {
                        if (is.null(tmpEnv[[iRank]])) {
                            tmpEnv[[iRank]] <- rank2index[[iRank]]
                        }
                    }
                } else {
                    tmpEnv[[iRank]] <- rank2index[[iRank]]
                }
            }
        }
        assignHash(ls(tmpEnv), getHash(ls(tmpEnv), tmpEnv), rank2index)
    }
    # convert env into dataframe and return
    index2RankList <- lapply(
        seq_len(length(allInputRank)), function (x) {
            data.frame(
                index = rank2index[[allInputRank[x]]],
                rank = allInputRank[x], stringsAsFactors = FALSE
            )
        }
    )
    index2RankDf <- do.call(rbind, index2RankList)
    index2RankDf <- index2RankDf[with(index2RankDf, order(index2RankDf$index)),]
    return(index2RankDf)
}

#' Align NCBI taxonomy IDs of list of taxa into a sorted rank list.
#' @export
#' @param idListFile a text file whose each row is a rank+ID list of a taxon
#' (see idListFile in example)
#' @param rankListFile a text file whose each row is a rank list of a taxon
#' (see rankListFile in example)
#' @return An aligned taxonomy dataframe which contains all the available
#' taxonomy ranks from the id and rank list file. This dataframe can be used for
#' creating a well resolved taxonomy tree (see ?createUnrootedTree) and sorting
#' taxa based on a selected reference taxon (see ?sortInputTaxa).
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table melt transpose setDT
#' @importFrom pbapply pblapply
#' @importFrom zoo na.locf
#' @seealso \code{\link{rankIndexing}}, \code{\link{createUnrootedTree}},
#' \code{\link{sortInputTaxa}}
#' @examples
#' idListFile <- system.file(
#'     "extdata", "data/idList.txt", package = "PhyloProfile", mustWork = TRUE
#' )
#' rankListFile <- system.file(
#'     "extdata", "data/rankList.txt", package = "PhyloProfile", mustWork = TRUE
#' )
#' taxonomyTableCreator(idListFile, rankListFile)

taxonomyTableCreator <- function(idListFile = NULL, rankListFile = NULL) {
    if (is.null(idListFile) | is.null(rankListFile))
        stop("Taxonomy ID list or rank list file is NULL!")
    index <- NULL
    index2RankDf <- rankIndexing(rankListFile)
    ncol <- max(utils::count.fields(rankListFile, sep = '\t'))
    idList <- utils::read.table(
        idListFile, sep = '\t', header = FALSE, check.names = FALSE,
        comment.char = "", fill = TRUE, stringsAsFactors = TRUE,
        na.strings = c("","NA"), col.names = paste0('X', seq_len(ncol)))
    colnames(idList)[1] <- "tip"
    orderedRank <- factor(index2RankDf$rank, levels = index2RankDf$rank)
    ### create a dataframe containing ordered ranks
    fullRankIDdf <- data.frame(
        rank = matrix(
            unlist(orderedRank), nrow = length(orderedRank), byrow = TRUE
        ), stringsAsFactors = FALSE)
    fullRankIDdf$index <- as.numeric(rownames(fullRankIDdf))
    fullRankIDdf <- data.table(fullRankIDdf)
    setkey(fullRankIDdf, rank)
    mTaxonDf <- pbapply::pblapply(
        seq_len(nrow(idList)),
        function (x) {
            ### get list of all IDs for this taxon & convert into long format
            taxonDf <- data.frame(idList[x,])
            taxonName <- unlist(
                strsplit(as.character(idList[x,]$tip), "#", fixed = TRUE))
            mTaxonDf <- suppressWarnings(
                data.table::melt(setDT(taxonDf), id = "tip")
            )
            ### get rank names and corresponding IDs, then remove NA cases
            splitCol <- data.frame(
                do.call(
                    'rbind',
                    strsplit(as.character(mTaxonDf$value), '#', fixed = TRUE)
                )
            )
            mTaxonDf <- cbind(mTaxonDf, splitCol)
            mTaxonDf <- mTaxonDf[stats::complete.cases(mTaxonDf),]
            ### subselect mTaxonDf to keep only 2 column rank id and rank name
            mTaxonDf <- mTaxonDf[, c("X1","X2")]
            if (mTaxonDf$X2[1] != index2RankDf$rank[1])
                mTaxonDf <- rbind(
                    data.frame(
                        "X1" = mTaxonDf$X1[1], "X2" = index2RankDf$rank[1]
                    ), mTaxonDf)
            ### rename columns & return
            colnames(mTaxonDf) <- c(taxonName[1], "rank")
            mTaxonDf <- data.table(mTaxonDf)
            setkey(mTaxonDf, rank)
            return(mTaxonDf)
        }
    )
    ### merge into data frame contains all available ranks from input
    mTaxonDfFull <- c(list(fullRankIDdf), mTaxonDf)
    fullRankIDdf <- Reduce(
        function (x, y) merge(x , y, all.x = TRUE, allow.cartesian = TRUE),
        mTaxonDfFull
    )
    ### reorder ranks & replace NA id by id of previous rank
    fullRankIDdf <- fullRankIDdf[order(fullRankIDdf$index),]
    fullRankIDdf <- zoo::na.locf(fullRankIDdf)
    ### remove index column & transpose into wide format
    fullRankIDdf <- subset(fullRankIDdf, select = -c(index))
    tFullRankIDdf <- data.table::transpose(fullRankIDdf)
    ### set first row to column names
    colnames(tFullRankIDdf) <- as.character(unlist(tFullRankIDdf[1,]))
    tFullRankIDdf <- tFullRankIDdf[-1,]
    ### add "abbrName  ncbiID  fullName" columns
    fullRankIDdfOut <- data.frame(
        abbrName = paste0("ncbi", unlist(tFullRankIDdf[,1])),
        ncbiID = unlist(tFullRankIDdf[,1]),
        fullName = colnames(fullRankIDdf)[-c(1)], stringsAsFactors = FALSE)
    fullRankIDdfOut <- cbind(fullRankIDdfOut, tFullRankIDdf)
    ### rename last column to "root"
    names(fullRankIDdfOut)[ncol(fullRankIDdfOut)] <- "root"
    numericCol <- c(seq(4, ncol(fullRankIDdfOut)))
    fullRankIDdfOut[,numericCol] <-
        as.numeric(as.character(unlist(fullRankIDdfOut[,numericCol])))
    return(fullRankIDdfOut)
}

#' Get taxonomy hierarchy for a list of taxon IDs
#' @description Get NCBI taxonomy hierarchy and URLs for an input taxon list.
#' @param inputTaxa NCBI ID list of input taxa.
#' @param currentNCBIinfo table/dataframe of the pre-processed NCBI taxonomy
#' data (/PhyloProfile/data/preProcessedTaxonomy.txt)
#' @return A list of dataframs containing taxonomy hierarchy and its URL to
#' NCBI database for input taxon IDs
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table fread
#' @examples
#' inputTaxa <- c("272557", "176299")
#' ncbiFilein <- system.file(
#'     "extdata", "data/preProcessedTaxonomy.txt",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' currentNCBIinfo <- as.data.frame(data.table::fread(ncbiFilein))
#' PhyloProfile:::getTaxHierarchy(inputTaxa, currentNCBIinfo)

getTaxHierarchy <- function(inputTaxa = NULL, currentNCBIinfo = NULL){
    if (is.null(currentNCBIinfo)) stop("Pre-processed NCBI tax data is NULL!")
    inputTaxaInfo <- PhyloProfile::getTaxonomyInfo(inputTaxa, currentNCBIinfo)
    ## get reduced taxonomy info (subset of preProcessedTaxonomy.txt)
    reducedDf <- unique(rbindlist(inputTaxaInfo))
    ## get list of all rank#IDs
    rankMod <- ncbiID <- NULL
    inputRankIDDf <- lapply(
        seq_len(length(inputTaxaInfo)),
        function (x) {
            inputTaxaInfo[[x]]$url <- paste0(
                "https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=",
                inputTaxaInfo[[x]]$ncbiID
            )
            inputTaxaInfo[[x]]$link <- paste0(
                "<p><a href='", inputTaxaInfo[[x]]$url, "' target='_blank'>",
                inputTaxaInfo[[x]]$fullName, "</a>",
                " (", inputTaxaInfo[[x]]$rank, ")", "</p>"
            )
            return(inputTaxaInfo[[x]][ , c("fullName", "ncbiID", "link")])
        }
    )
    return(inputRankIDDf)
}

#' Get taxon names for a list of taxon IDs
#' @export
#' @param idList list of taxonomy IDs
#' @param currentNCBIinfo table/dataframe of the pre-processed NCBI taxonomy
#' data (/PhyloProfile/data/preProcessedTaxonomy.txt)
#' @return A dataframe contains input taxon Ids and their full names.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table fread
#' @examples
#' ncbiFilein <- system.file(
#'     "extdata", "data/preProcessedTaxonomy.txt",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' currentNCBIinfo <- as.data.frame(data.table::fread(ncbiFilein))
#' idList <- c("9606", "5207", "40674", "4751")
#' id2name(idList, currentNCBIinfo)

id2name <- function(idList = NULL, currentNCBIinfo = NULL) {
    if (is.null(idList)) stop("No list of taxon IDs given!")
    nameList <- data.frame(
        ncbiID = idList,
        fullName = paste0("NCBI taxID ", idList)
    )
    if (is.null(currentNCBIinfo)) {
        ncbiFilein <- system.file(
            "PhyloProfile", "data/preProcessedTaxonomy.txt",
            package = "PhyloProfile", mustWork = FALSE
        )
        if (file.exists(ncbiFilein)) {
            currentNCBIinfo <- as.data.frame(data.table::fread(ncbiFilein))
            # get taxon names
            nameList <- currentNCBIinfo[
                currentNCBIinfo$ncbiID %in% idList, c("ncbiID","fullName")
            ]
        }
    } else {
        nameListTmp <- currentNCBIinfo[
            currentNCBIinfo$ncbiID %in% idList, c("ncbiID","fullName")
        ]
        nameList <- merge(nameList, nameListTmp, by = "ncbiID", all.x = TRUE)
        nameList$fullName.y[is.na(nameList$fullName.y)] <-
            nameList$fullName.x[is.na(nameList$fullName.y)]
        nameList <- nameList[, c("ncbiID","fullName.y")]
        colnames(nameList) <- c("ncbiID", "fullName")
    }
    return(nameList)
}
