# Parse main input functions
# These functions check the validity of the main input file
# and convert different input format into long format.

#' Check the validity of the input phylogenetic profile file
#' @description Check if input file has one of the following format: orthoXML,
#' multiple FASTA, tab-delimited matrix (wide or long), or list of OMA IDs.
#' @export
#' @param filein input file
#' @return The format of the input file format, or type of error
#' @importFrom data.table fread
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{checkOmaID}}
#' @examples
#' filein <- system.file(
#'     "extdata", "test.main.wide", package = "PhyloProfile", mustWork = TRUE
#' )
#' checkInputValidity(filein)

checkInputValidity <- function(filein) {
    if (missing(filein) | is.null(filein)) return("No input file given!")
    inputDt <- data.frame(fread(
        file = filein, sep = "\t", header = FALSE, check.names = FALSE,
        fill = TRUE, stringsAsFactors = FALSE))

    if (is.na(inputDt[1, ncol(inputDt)])) {
        return("moreCol")
    } else {
        names(inputDt) <- as.character(unlist(inputDt[1, ]))
        # XML format (starts with <?xml)
        if (grepl("<?xml", colnames(inputDt)[1])) {
            return("xml")
        }
        # FASTA format (starts with ">" )
        else if (grepl(">", colnames(inputDt)[1]) == TRUE) {
            return("fasta")
        }
        # LONG or WIDE format (starts with "geneID")
        else {
            if (grepl("geneID", colnames(inputDt)[1])) {
                # LONG format
                if (is.na(pmatch("ncbi", colnames(inputDt)[3])) ||
                    is.na(pmatch("ncbi", colnames(inputDt)[4])) ||
                    is.na(pmatch("ncbi", colnames(inputDt)[5]))) {
                    return("long")
                }
                # WIDE format
                else {
                    tmp <- inputDt[inputDt == ""][1]
                    if (!is.na(tmp) & tmp == "") {
                        return("emptyCell")
                    } else return("wide")
                }
            }
            else {
                if (grepl("\\s", inputDt[1, 1])) return("invalidFormat")
                # OMA ids
                ids <- unique(inputDt[,1])
                if (length(ids) < 2) ids <- rep(ids, 2)
                invalidOma <- checkOmaID(ids)
                if (length(invalidOma) == 0) {
                    return("oma")
                } else return(invalidOma)
            }
        }
    }
}

#' Parse orthoXML input file
#' @export
#' @param inputFile input file in xml format
#' @return A data frame of input data in long-format containing seed gene IDs (
#' or orthologous group IDs), their orthologous proteins together with the
#' corresponding taxonomy IDs and values of (up to) two additional variables.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_find_all
#' @importFrom xml2 xml_attr
#' @importFrom data.table dcast setDT
#' @examples
#' inputFile <- system.file(
#'     "extdata", "test.main.xml", package = "PhyloProfile", mustWork = TRUE
#' )
#' xmlParser(inputFile)

xmlParser <- function(inputFile = NULL){
    if (is.null(inputFile)) stop("Input file is NULL!")
    scoreType <- scoreValue <- NULL
    data <- read_xml(inputFile)
    # get all genes for each taxon
    species <- xml_find_all(data, ".//species")
    geneSpec <- xml_find_all(species, ".//gene")
    refGeneID <- xml_attr(geneSpec, "id")
    if (length(geneSpec) == 0)
        stop("Unsupported XML format! Please check the file header!")
    if (grepl("protId", geneSpec[1])) {
        orthoID <- as.data.frame(unlist(strsplit(
            xml_attr(geneSpec, "protId"), split = '|', fixed = TRUE
        )), stringsAsFactors = FALSE)[,1]
    } else if (grepl("geneId", geneSpec[1])) orthoID <- xml_attr(geneSpec, "id")
    speciesID <- xml_attr(species, "NCBITaxId")
    speciesIDrep <- unlist(
        lapply(species, function (x) length(xml_find_all(x, ".//gene"))))
    speciesDf <- data.frame(
        ncbiID = rep(paste0("ncbi", speciesID), speciesIDrep),
        refGeneID = refGeneID, orthoID = orthoID, stringsAsFactors=FALSE)
    # get orthologs and their scores
    orthoGroup <- xml_find_all(data, ".//orthologGroup")
    groupID <- xml_attr(orthoGroup, "id")
    genes <- xml_find_all(orthoGroup, ".//geneRef")
    refGeneID <- xml_attr(genes, "id")
    score <- xml_find_all(genes, ".//score")
    if (length(score) > 0) {
        scorePair <- lapply(
            score, function (x) c(xml_attr(x, "id"), xml_attr(x, "value")))
        scorePair <- as.data.frame(do.call(rbind, scorePair))
        groupIDrep <- unlist(lapply(
            orthoGroup,
            function (x) ncol(scorePair)*length(xml_find_all(x, ".//geneRef")))
        )
        orthoDf <- data.frame(
            geneID = rep(groupID, groupIDrep),
            refGeneID = rep(refGeneID, each = ncol(scorePair)),
            scoreType = scorePair$V1, scoreValue = scorePair$V2,
            stringsAsFactors = FALSE)
        orthoDf <- data.table::dcast(
            setDT(orthoDf), geneID+refGeneID ~ scoreType, value.var="scoreValue"
        )
    } else {
        groupIDrep <- unlist(lapply(
            orthoGroup,
            function (x) length(xml_find_all(x, ".//geneRef")))
        )
        orthoDf <- data.frame(
            geneID = rep(groupID, groupIDrep),
            refGeneID = refGeneID,
            stringsAsFactors = FALSE)
    }
    # merge into final dataframe
    finalDf <- merge(speciesDf, orthoDf, all.y = TRUE, by = "refGeneID")
    # remove refGeneID column and reorder columns
    finalDf <- finalDf[, -1]
    refcols <- c("geneID", "ncbiID", "orthoID")
    finalDf <- finalDf[, c(refcols, setdiff(names(finalDf), refcols))]
    # remove columns that contains only NA
    finalDf <- finalDf[, colSums(is.na(finalDf)) < nrow(finalDf)]
    return(finalDf)
}

#' Parse multi-fasta input file
#' @export
#' @param inputFile input multiple fasta file. Check extdata/test.main.fasta or
#' https://github.com/BIONF/PhyloProfile/wiki/Input-Data#multi-fasta-format for
#' the supported FASTA header.
#' @return A data frame of input data in long-format containing seed gene IDs (
#' or orthologous group IDs), their orthologous proteins together with the
#' corresponding taxonomy IDs and values of (up to) two additional variables.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom Biostrings readAAStringSet
#' @examples
#' inputFile <- system.file(
#'     "extdata", "test.main.fasta", package = "PhyloProfile", mustWork = TRUE
#' )
#' fastaParser(inputFile)

fastaParser <- function(inputFile = NULL){
    if (is.null(inputFile)) stop("Input file is NULL!")
    # split sequence IDs into columns
    fastaFile <- Biostrings::readAAStringSet(inputFile)
    seqID <- names(fastaFile)
    tmpDf <- data.frame(
        do.call(rbind, strsplit(seqID, "\\|")), stringsAsFactors = FALSE
    )
    orthoID <- tmpDf[,3]
    if (!is.null(nrow(tmpDf[,-c(1,2, ncol(tmpDf) - 1, ncol(tmpDf))])))
        orthoID <- do.call(paste, tmpDf[,-c(1,2, ncol(tmpDf) - 1, ncol(tmpDf))])
    faDf <- data.frame(
        geneID = tmpDf[,1],
        ncbiID = tmpDf[,2],
        orthoID = orthoID,
        var1 = tmpDf[, ncol(tmpDf) - 1],
        var2 = tmpDf[, ncol(tmpDf)],
        stringsAsFactors = FALSE
    )
    faDf$orthoID <- gsub(" ", "|", faDf$orthoID)

    # remove columns that contains only NA
    faDf <- faDf[, colSums(is.na(faDf)) < nrow(faDf)]
    return(faDf)
}

#' Transform input file in wide matrix into long matrix format
#' @export
#' @param inputFile input file in wide matrix format
#' @return A data frame of input data in long-format containing seed gene IDs (
#' or orthologous group IDs), their orthologous proteins together with the
#' corresponding taxonomy IDs and values of (up to) two additional variables.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table fread
#' @examples
#' inputFile <- system.file(
#'     "extdata", "test.main.wide", package = "PhyloProfile", mustWork = TRUE
#' )
#' wideToLong(inputFile)

wideToLong <- function(inputFile = NULL){
    if (is.null(inputFile)) stop("Input file is NULL!")
    wideDataframe <- data.frame(fread(
        file = inputFile, sep = "\t", header = TRUE, check.names = FALSE,
        stringsAsFactors = FALSE
    ))
    ncbiIDs <- colnames(wideDataframe[, c(-1)])
    orthoInfo <- data.frame(
        do.call(
            rbind, strsplit(as.character(unlist(wideDataframe[, c(-1)])), "#")
        ),
        stringsAsFactors = FALSE
    )

    longDataframe <- data.frame(
        geneID = rep(wideDataframe$geneID, time = ncol(wideDataframe) - 1),
        ncbiID = rep(ncbiIDs, time = 1, each = nrow(wideDataframe)),
        orthoID = orthoInfo$X1,
        var1 = suppressWarnings(as.numeric(orthoInfo$X2)),
        var2 = suppressWarnings(as.numeric(orthoInfo$X3)),
        stringsAsFactors = FALSE
    )
    return(longDataframe)
}

#' Create a long matrix format for all kinds of input phylogenetic profiles
#' @export
#' @param inputFile input profile file in orthoXML, multiple FASTA,
#' tab-delimited matrix format (wide or long).
#' @return A data frame of input data in long-format containing seed gene IDs (
#' or orthologous group IDs), their orthologous proteins together with the
#' corresponding taxonomy IDs and values of (up to) two additional variables.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table fread
#' @seealso \code{\link{xmlParser}}, \code{\link{fastaParser}},
#' \code{\link{wideToLong}}
#' @examples
#' inputFile <- system.file(
#'     "extdata", "test.main.wide", package = "PhyloProfile", mustWork = TRUE
#' )
#' createLongMatrix(inputFile)

createLongMatrix <- function(inputFile = NULL){
    if (is.null(inputFile)) stop("Input file is NULL!")

    inputType <- checkInputValidity(inputFile)
    # XML
    if (inputType == "xml") longDataframe <- xmlParser(inputFile)
    # FASTA
    else if (inputType == "fasta") longDataframe <- fastaParser(inputFile)
    # LONG
    else if (inputType == "long") {
        longDataframe <- data.frame(fread(
            file = inputFile, sep = "\t", header = TRUE, check.names = FALSE, 
            stringsAsFactors = FALSE
        ))
    }
    # WIDE
    else if (inputType == "wide") longDataframe <- wideToLong(inputFile)
    else return(NULL)

    # convert geneID, ncbiID and orthoID into factor and var1, var2 into numeric
    for (i in seq_len(3)) {
        longDataframe[, i] <- as.factor(longDataframe[, i])
    }
    if (ncol(longDataframe) > 3 & ncol(longDataframe) < 6) {
        for (j in seq(4, ncol(longDataframe))){
            longDataframe[,j] <- suppressWarnings(
                as.numeric(as.character(longDataframe[,j]))
            )
        }
    }
    if (ncol(longDataframe) == 6)
        longDataframe[, 6] <- as.factor(longDataframe[, 6])

    # remove duplicated lines
    longDataframe <- longDataframe[!duplicated(longDataframe),]
    return(longDataframe)
}
