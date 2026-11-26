#' Calculate the phylogenetic gene age from the phylogenetic profiles
#' @export
#' @usage estimateGeneAge(processedProfileData, taxaCount, rankName, refTaxon,
#'     var1CO, var2CO, percentCO, taxDB = NULL)
#' @param processedProfileData dataframe contains the full processed
#' phylogenetic profiles (see ?fullProcessedProfile or ?parseInfoProfile)
#' @param taxaCount dataframe counting present taxa in each supertaxon
#' @param rankName working taxonomy rank (e.g. "strain", "species", "genus",
#' "family", "order")
#' @param refTaxon reference taxon name (e.g. "Homo sapiens", "Homo" or
#' "Hominidae")
#' @param var1CO cutoff for var1. Default: c(0, 1)
#' @param var2CO cutoff for var2. Default: c(0, 1)
#' @param percentCO cutoff for percentage of species present in each
#' supertaxon. Default: c(0, 1)
#' @param taxDB Path to the taxonomy DB files
#' @return A dataframe contains estimated gene ages for the seed proteins.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom data.table setDT setnames
#' @importFrom dplyr count
#' @seealso \code{\link{parseInfoProfile}} for creating a full processed
#' profile dataframe; \code{\link{getNameList}} and
#' \code{\link{getTaxonomyMatrix}} for getting taxonomy info,
#' \code{\link{fullProcessedProfile}} for a demo input dataframe
#' @examples
#' library(dplyr)
#' data("fullProcessedProfile", package="PhyloProfile")
#' rankName <- "class"
#' refTaxon <- "Mammalia"
#' processedProfileData <- fullProcessedProfile
#' taxonIDs <- levels(as.factor(processedProfileData$ncbiID))
#' sortedInputTaxa <- sortInputTaxa(
#'     taxonIDs, rankName, refTaxon, NULL, NULL
#' )
#' taxaCount <- sortedInputTaxa %>% dplyr::count(supertaxon)
#' var1Cutoff <- c(0, 1)
#' var2Cutoff <- c(0, 1)
#' percentCutoff <- c(0, 1)
#' estimateGeneAge(
#'     processedProfileData,
#'     taxaCount,
#'     rankName,
#'     refTaxon,
#'     var1Cutoff, var2Cutoff, percentCutoff
#' )

estimateGeneAge <- function(
    processedProfileData, taxaCount, rankName, refTaxon,
    var1CO = c(0, 1), var2CO = c(0, 1), percentCO = c(0, 1), taxDB = NULL
){
    supportedRanks <- c(
        "strain", "species", "genus", "family", "order", "class", "phylum",
        "kingdom", "superkingdom"
    )
    if (!(rankName %in% supportedRanks))
        stop(
            "Invalid taxonomy rank! Supported ranks: ",
            paste(supportedRanks, collapse = ", ")
        )
    taxList <- getNameList(taxDB)
    getName <- function(id) taxList$fullName[taxList$ncbiID %in% id][1]
    # full non-duplicated taxonomy data
    Dt <- getTaxonomyMatrix(taxDB)
    # ranks used to compare each taxon with the reference taxon, ordered from
    # the lowest to the highest. Each position of a category string ("cath")
    # corresponds to one element of c("abbrName", rankList), so the length of
    # cath and every category are derived from rankList (never hardcoded)
    rankList <- c(
        "strain", "species", "genus", "family", "order", "class", "phylum",
        "kingdom", "norank_33154", "superkingdom", "root"
    )
    rankList <- rankList[rankList %in% colnames(Dt)]
    if (!(rankName %in% rankList))
        stop("Rank ", rankName, " not found in the taxonomy matrix!")
    catCols <- c("abbrName", rankList)
    nCat <- length(catCols)
    # category of taxa sharing the reference's taxon from rank r upwards,
    # e.g. for 12 positions getCat("genus") = "000111111111"
    getCat <- function(r) {
        nZero <- match(r, catCols) - 1
        paste0(strrep("0", nZero), strrep("1", nCat - nZero))
    }
    # two-digit prefix of the age labels (see geneAgePlotDf)
    ageCode <- c(
        species = "01", genus = "02", family = "03", order = "04",
        class = "05", phylum = "06", kingdom = "07", norank_33154 = "08",
        superkingdom = "10", root = "11"
    )
    
    # get selected (super)taxon ID. Strains often have no "strain" rank in
    # the NCBI taxonomy (e.g. "norank"), so fall back to the name only
    superID <- taxList$ncbiID[
        taxList$fullName == refTaxon & taxList$rank == rankName]
    if (length(superID) == 0 || !any(Dt[, rankName] %in% superID))
        superID <- taxList$ncbiID[taxList$fullName == refTaxon]
    superID <- superID[superID %in% Dt[, rankName]]
    if (length(superID) == 0)
        stop("Reference taxon ", refTaxon, " not found at rank ", rankName)
    # get (super)taxa IDs for one of representative taxa
    firstLine <- Dt[Dt[, rankName] %in% superID, ][1, ]
    supFirstLine <- firstLine[, catCols]
    # subset of taxonomy data, containing only ranks from rankList
    subDt <- Dt[, catCols]
    # compare each taxon IDs with selected taxon & create a "category" DF
    subDtTmp <- subDt
    if ("norank_33154" %in% rankList && "superkingdom" %in% rankList) {
        # merge Unikonta (33154, 554915) or Bikonta into one group
        isUnikont <- subDtTmp$norank_33154 %in% c(554915, 33154)
        isEuk <- subDtTmp$superkingdom == 2759
        if (supFirstLine$norank_33154 %in% c(554915, 33154)) {
            subDtTmp$norank_33154[isUnikont & isEuk] <-
                supFirstLine$norank_33154
        } else {
            subDtTmp$norank_33154[!isUnikont & isEuk] <-
                supFirstLine$norank_33154
        }
    }
    # compare as character: abbrName is a factor, all other columns are IDs
    refIDs <- vapply(supFirstLine, function(x) as.character(x), character(1))
    catMatrix <- vapply(
        subDtTmp, function(col) as.integer(as.character(col) %in% refIDs),
        integer(nrow(subDtTmp))
    )
    if (is.null(dim(catMatrix))) catMatrix <- matrix(catMatrix, nrow = 1)
    cath <- apply(catMatrix, 1, paste0, collapse = "")
    # keep only the uninterrupted run of 1s at the highest ranks, so that every
    # cath is one of getCat(<rank>) (e.g. "000101111" -> "000001111")
    ones <- nchar(sub("^.*?(1*)$", "\\1", cath))
    cath <- paste0(strrep("0", nCat - ones), strrep("1", ones))
    catDf <- data.frame(
        ncbiID = as.character(subDtTmp$abbrName), cath = cath,
        stringsAsFactors = FALSE
    )
    # get main input data
    mdData <- droplevels(processedProfileData)
    # filter by var1, var2 ..
    mdData <- subset(
        mdData, mdData$var1 >= var1CO[1] & mdData$var1 <= var1CO[2]
        & mdData$var2 >= var2CO[1] & mdData$var2 <= var2CO[2])
    # calculate % present taxa
    finalPresSpecDt <- calcPresSpec(mdData, taxaCount)
    mdData <- mdData[
        , c("geneID", "supertaxon","ncbiID","orthoID","var1","var2")
    ]
    mdDataFull <- merge(
        mdData, finalPresSpecDt, by = c("geneID", "supertaxon"), all.x = TRUE
    )
    # add "category" into mdData
    mdDtExt <- merge(mdDataFull, catDf, by = "ncbiID", all.x = TRUE)
    mdDtExt$var1[mdDtExt$var1 == "NA" | is.na(mdDtExt$var1)] <- 0
    mdDtExt$var2[mdDtExt$var2 == "NA" | is.na(mdDtExt$var2)] <- 0
    # remove cat for "NA" orthologs and also for orthologs that dont fit cutoffs
    mdDtExt$cath[mdDtExt$orthoID == "NA" | is.na(mdDtExt$orthoID)] <- NA
    mdDtExt <- mdDtExt[stats::complete.cases(mdDtExt), ]
    # filter by %specpres
    mdDtExt <- mdDtExt[
        mdDtExt$presSpec >= percentCO[1] & mdDtExt$presSpec <= percentCO[2],
    ]
    
    # get the furthest common taxon with selected taxon for each gene
    minCat <- tapply(mdDtExt$cath, mdDtExt$geneID, min)
    geneAgeDf <- data.frame(
        geneID = names(minCat), cath = as.vector(minCat),
        stringsAsFactors = FALSE
    )
    ### genes without any valid ortholog and ranks that are lower than the
    ### ref rank are moved to the working (ref) rank
    refCat <- getCat(rankName)
    geneAgeDf$cath[is.na(geneAgeDf$cath) | geneAgeDf$cath > refCat] <- refCat
    
    ### convert cat into geneAge
    getOrthos <- function(gene, cat) {
        mdDtExt$ncbiID[mdDtExt$geneID == gene & mdDtExt$cath == cat]
    }
    # root: shared by 2 or 3 superkingdoms
    rootAge <- function(orthoList) {
        orthoDomainID <- unique(
            subDt$superkingdom[subDt$abbrName %in% orthoList]
        )
        orthoDomainName <- sort(unique(c(
            getName(supFirstLine$superkingdom),
            taxList$fullName[taxList$ncbiID %in% orthoDomainID]
        )))
        code <- if (length(orthoDomainName) >= 3) "12" else "11"
        paste0(code, "_", paste(orthoDomainName, collapse = "-"))
    }
    # norank_33154: Opisthokonta/Amoebozoa..., Unikonta/Bikonta or Eukaryote
    konAge <- function(orthoList) {
        if (supFirstLine$superkingdom != 2759)
            return(paste0("10_", getName(supFirstLine$superkingdom)))
        orthoKonID <- unique(c(
            as.character(subDt$norank_33154[subDt$abbrName %in% orthoList]),
            as.character(supFirstLine$norank_33154)
        ))
        if (length(orthoKonID) == 1)
            return(paste0("08_", getName(orthoKonID)))
        ukon <- unique(orthoKonID %in% c(33154, 554915))
        if (length(ukon) == 1)
            return(ifelse(ukon, "09_Unikonta", "09_Bikonta"))
        return("10_Eukaryote")
    }
    # kingdom: for Sar, check for a common subgroup of Sar
    kingdomAge <- function(orthoList) {
        kingdomName <- getName(supFirstLine$kingdom)
        sarCols <- intersect(
            colnames(Dt), c("norank_33630", "norank_543769", "norank_33634"))
        if (!identical(kingdomName, "Sar") || length(sarCols) == 0)
            return(paste0("07_", kingdomName))
        sarTaxDt <- Dt[Dt$abbrName %in% orthoList, sarCols, drop = FALSE]
        common <- unlist(lapply(
            sarTaxDt, function(x) if (length(unique(x)) == 1) unique(x)))
        if (is.null(common)) return(paste0("07_", kingdomName))
        return(paste0("07_", getName(common)))
    }
    
    geneAgeDf$age <- NA_character_
    refCode <- if (rankName == "strain") "00" else ageCode[[rankName]]
    geneAgeDf$age[geneAgeDf$cath == refCat] <- paste0(refCode, "_", refTaxon)
    higherRanks <- rankList[-seq_len(match(rankName, rankList))]
    for (r in higherRanks) {
        rCat <- getCat(r)
        idx <- which(geneAgeDf$cath == rCat)
        if (length(idx) == 0) next
        if (r %in% c("root", "norank_33154", "kingdom")) {
            ageFn <- switch(
                r, root = rootAge, norank_33154 = konAge, kingdom = kingdomAge)
            geneAgeDf$age[idx] <- vapply(
                geneAgeDf$geneID[idx],
                function(g) ageFn(getOrthos(g, rCat)), character(1)
            )
        } else {
            geneAgeDf$age[idx] <- paste0(
                ageCode[[r]], "_", getName(supFirstLine[[r]]))
        }
    }
    # return geneAge data frame
    geneAgeDf <- geneAgeDf[, c("geneID", "cath", "age")]
    return(geneAgeDf)
}

#' Create data for plotting gene ages
#' @param geneAgeDf data frame containing estimated gene ages for seed proteins
#' @return A dataframe for plotting gene age plot containing the absolute number
#' and percentage of genes for each calculated evolutionary ages and the
#' corresponding position for writting those number on the plot.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom dplyr count
#' @export
#' @seealso \code{\link{estimateGeneAge}}
#' @examples
#' geneAgeDf <- data.frame(
#' geneID = c("100136at6656", "100265at6656", "101621at6656", "103479at6656"),
#' cath = c("0000001", "0000011", "0000001", "0000011"),
#' age = c("12_LUCA", "10_Eukaryota", "12_LUCA", "10_Eukaryota")
#' )
#' geneAgePlotDf(geneAgeDf)

geneAgePlotDf <- function(geneAgeDf){
    if (is.null(geneAgeDf)) stop("Gene age data is NULL!")
    age <- NULL
    plotDf <- geneAgeDf %>% dplyr::count(age)
    plotDf$level <- substring(plotDf$age, 1,2)
    levelDf <- data.frame(
        level = c("00", "01", "02", "03", "04", "05", "06", "07", "10", "12"),
        rank = c(
            " (Strain)", " (Species)", " (Genus)", " (Family)", " (Order)",
            " (Class)", " (Phylum)", " (Kingdom)", " (Superkingdom)", " (LUCA)"
        ),
        stringsAsFactors = FALSE
    )
    plotDf <- merge(plotDf, levelDf, by = "level", all.x = TRUE)
    plotDf[is.na(plotDf)] <- ""
    plotDf$name <- paste0(substring(plotDf$age, 4), plotDf$rank)
    plotDf$name <- factor(plotDf$name, levels = plotDf$name)
    plotDf$percentage <- round(plotDf$n / sum(plotDf$n) * 100, 1)
    outDf <- plotDf[, c("name", "n", "percentage")]
    colnames(outDf) <- c("name", "count", "percentage")
    return(outDf)
}

#' Create gene age plot
#' @param geneAgePlotDf data frame required for plotting gene age (see
#' ?geneAgePlotDf)
#' @param textFactor increase factor of text size
#' @param font font of text. Default = Arial"
#' @return A gene age distribution plot as a ggplot2 object
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @export
#' @seealso \code{\link{estimateGeneAge}} and \code{\link{geneAgePlotDf}}
#' @examples
#' geneAgePlotDf <- data.frame(
#'     name = c("Streptophyta (Phylum)", "Bikonta", "Eukaryota (Superkingdom)"),
#'     count = c(7, 1, 30),
#'     percentage = c(18, 3, 79)
#' )
#' createGeneAgePlot(geneAgePlotDf, 1, "sans")

createGeneAgePlot <- function (geneAgePlotDf, textFactor = 1, font = "Arial"){
    name <- count <- percentage <- x <- y <- NULL
    arrowDf <- data.frame(time = c(0, 0), y = c(0, nrow(geneAgePlotDf) + 1))
    n <- nlevels(as.factor(geneAgePlotDf$name))
    mycolors <- grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(11, "Spectral"))(n)
    p <- ggplot(data = geneAgePlotDf, aes(x = name, y = count, fill = name)) +
        geom_bar(stat = "identity", width = 0.5) +
        geom_text(
            aes(label = paste0(percentage, "%")),
            position = position_dodge(0.9),
            size = 3.5 + 1 * textFactor, hjust = 0
        ) +
        geom_line(
            data = data.frame(x = c(0, 0), y = c(0, nrow(geneAgePlotDf) + 1)),
            aes(y = x, x = y),
            arrow = arrow(
                length = unit(0.30,"cm"), ends = "last", type = "closed"
            )
        ) +
        geom_text(
            data = NULL,
            aes(y = 0, x = (nrow(geneAgePlotDf) + 1) * 0.92, label = "time"),
            vjust = -0.5,
            angle = 90, alpha = 0.2
        ) +
        coord_flip() +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = 12 * textFactor),
            axis.text = element_text(size = 12 * textFactor),
            text = element_text(family = font)
        ) +
        scale_fill_manual(values = mycolors) +
        scale_x_discrete(limits = rev(levels(geneAgePlotDf$name))) +
        labs(y = "Number of genes")
    return(p)
}
