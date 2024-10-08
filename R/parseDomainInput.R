#' Parse domain input file
#' @description Get all domain annotations for one seed protein IDs.
#' @export
#' @param seed seed protein ID
#' @param inputFile name of input file (file name or path to folder contains
#' individual domain files)
#' @param type type of data (file" or "folder"). Default = "file".
#' @return A dataframe for protein domains including seed ID, its orthologs IDs,
#' sequence lengths, feature names, start and end positions, feature weights
#' (optional) and the status to determine if that feature is important for
#' comparison the architecture between 2 proteins* (e.g. seed protein vs
#' ortholog) (optional).
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom stringr str_split_fixed
#' @importFrom data.table fread
#' @seealso \code{\link{getDomainFolder}}
#' @examples
#' seed <- "101621at6656"
#' inputFile <- system.file(
#'     "extdata", "domainFiles/101621at6656.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' type <- "file"
#' parseDomainInput(seed, inputFile, type)

parseDomainInput <- function(seed = NULL, inputFile = NULL, type = "file") {
    file <- NULL
    # check parameters
    if (type == "folder" & is.null(seed)) stop("Seed ID cannot be NULL!")
    if (is.null(inputFile)) stop("Input domain file cannot be NULL!")
    # get domain from single file
    if (type == "file") {
        file <- inputFile
    }
    # or from a domain folder
    else {
        file <- getDomainFolder(seed, inputFile)
        if (file == "noSelectHit") {return("noSelectHit")}
        else if (file == "noFileInFolder") {return("noFileInFolder")}
    }
    if (file != FALSE) {
        exeptions <- c("noFileInput", "noSelectHit", "noFileInFolder")
        if (!(file %in% exeptions)) {
            domains <- data.frame(fread(file, sep = "\t", header = FALSE))
        }
    }
    if (ncol(domains) == 5) {
        colnames(domains) <-
            c("seedID", "orthoID", "feature", "start", "end")
    } else if (ncol(domains) == 6) {
        colnames(domains) <-
            c("seedID", "orthoID", "length", "feature", "start", "end")
    } else if (ncol(domains) == 7) {
        colnames(domains) <-
            c("seedID", "orthoID", "length", "feature", "start", "end","weight")
    } else if (ncol(domains) == 8) {
        colnames(domains) <- c(
            "seedID", "orthoID", "length", "feature", "start", "end","weight",
            "path")
    } else if (ncol(domains) == 14) {
        domains <- domains[-1,]
        colnames(domains) <- c(
            "seedID", "orthoID", "length", "feature", "start", "end","weight",
            "path","acc","evalue","bitscore","pStart","pEnd","pLen")
        domains$length <- as.integer(domains$length)
        domains$start <- as.integer(domains$start)
        domains$end <- as.integer(domains$end)
        domains$pStart <- as.integer(domains$pStart)
        domains$pEnd <- as.integer(domains$pEnd)
        domains$pLen <- as.integer(domains$pLen)
        domains$evalue <- as.numeric(domains$evalue)
        domains$bitscore <- as.numeric(domains$bitscore)
    } else {
        return("ERR")
    }
    if (nrow(domains) == 0) return("ERR-0")
    domains$seedID <- as.character(domains$seedID)
    domains$orthoID <- as.character(domains$orthoID)
    domains$seedID <- gsub("\\|",":",domains$seedID)
    domains$orthoID <- gsub("\\|",":",domains$orthoID)
    domains[c("feature_type","feature_id")] <-
        stringr::str_split_fixed(domains$feature, '_', 2)
    domains$feature_id[domains$feature_type == "smart"] <-
        paste0(domains$feature_id[domains$feature_type == "smart"], "_smart")
    return(domains)
}

#' Get domain file from a folder for a seed protein
#' @param seed seed protein ID
#' @param domainPath path to domain folder
#' @return Domain file and its complete directory path for the selected protein.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' domainPath <- paste0(
#'     path.package("PhyloProfile", quiet = FALSE), "/extdata/domainFiles"
#' )
#' PhyloProfile:::getDomainFolder("101621at6656", domainPath)

getDomainFolder <- function(seed, domainPath){
    if (is.null(seed)) {
        fileDomain <- "noSelectHit"
    } else {
        # check file extension
        allExtension <- c("txt", "csv", "list", "domains", "architecture")
        fileDomain <- paste0(domainPath, "/", seed, ".", allExtension)

        checkExistance <- lapply(fileDomain, function (x) file.exists(x))
        fileDomain <- fileDomain[match(TRUE, checkExistance)]

        if (is.na(fileDomain)) fileDomain <- "noFileInFolder"
    }
    return(fileDomain)
}
