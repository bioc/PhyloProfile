#' Calculate the phylogenetic gene age from the phylogenetic profiles
#' @export
#' @usage estimateGeneAge(processedProfileData, taxaCount, rankName, refTaxon,
#'     var1CO, var2CO, percentCO, taxDB = NULL)
#' @param processedProfileData dataframe contains the full processed
#' phylogenetic profiles (see ?fullProcessedProfile or ?parseInfoProfile)
#' @param taxaCount dataframe counting present taxa in each supertaxon
#' @param rankName working taxonomy rank (e.g. "species", "genus", "family")
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
    if (rankName == "strain") stop("Please select a higher rank than STRAIN!")
    rankList <- c(
        "genus", "family", "class", "phylum", "kingdom", "norank_33154",
        "superkingdom", "root"
    )
    # get selected (super)taxon ID
    taxList <- getNameList(taxDB)
    superID <- taxList[
        taxList$fullName == refTaxon & taxList$rank == rankName,]$ncbiID
    # full non-duplicated taxonomy data
    Dt <- getTaxonomyMatrix(taxDB)
    # subset of taxonomy data, containing only ranks from rankList
    subDt <- Dt[, c("abbrName", rankList)]
    # get (super)taxa IDs for one of representative species
    firstLine <- Dt[Dt[, rankName] == superID, ][1, ]
    supFirstLine <- firstLine[, c("abbrName", rankList)]
    # compare each taxon IDs with selected taxon & create a "category" DF
    subDtTmp <- subDt
    if (supFirstLine$norank_33154 %in% c(554915, 33154)) {
        subDtTmp$norank_33154[
            subDtTmp$norank_33154 %in% c(554915, 33154) &
                subDtTmp$superkingdom == 2759
            ] <- supFirstLine$norank_33154
    } else {
        subDtTmp$norank_33154[
            !(subDtTmp$norank_33154 %in% c(554915, 33154)) &
                subDtTmp$superkingdom == 2759
            ] <- supFirstLine$norank_33154
    }
    catList <- lapply(
        seq(nrow(subDtTmp)), function (x) {
            cath <- subDtTmp[x, ] %in% supFirstLine
            cath <- paste0(cath, collapse = "")})
    catDf <- data.frame(ncbiID = as.character(subDtTmp$abbrName),
                        cath = do.call(rbind, catList), stringsAsFactors=FALSE)
    catDf$cath <- gsub("TRUE", "1", catDf$cath)
    catDf$cath <- gsub("FALSE", "0", catDf$cath)
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
    mdDataFull <- Reduce(
        function(x, y) merge(x, y, by = c("geneID", "supertaxon"), all.x=TRUE),
        list(mdData, finalPresSpecDt))
    # add "category" into mdData
    mdDtExt <- merge(mdDataFull, catDf, by = "ncbiID", all.x = TRUE)
    mdDtExt$var1[mdDtExt$var1 == "NA" | is.na(mdDtExt$var1)] <- 0
    mdDtExt$var2[mdDtExt$var2 == "NA" | is.na(mdDtExt$var2)] <- 0
    # remove cat for "NA" orthologs and also for orthologs that dont fit cutoffs
    if (nrow(mdDtExt[mdDtExt$orthoID == "NA" | is.na(mdDtExt$orthoID), ]) > 0)
        mdDtExt[mdDtExt$orthoID == "NA" | is.na(mdDtExt$orthoID),]$cath <- NA
    mdDtExt <- mdDtExt[stats::complete.cases(mdDtExt), ]
    # filter by %specpres
    mdDtExt <- mdDtExt[
        mdDtExt$presSpec >=percentCO[1] & mdDtExt$presSpec<= percentCO[2],
        ]

    # get the furthest common taxon with selected taxon for each gene
    geneAgeDf <- as.data.frame(tapply(mdDtExt$cath, mdDtExt$geneID, min))
    data.table::setDT(geneAgeDf, keep.rownames = TRUE)[]
    data.table::setnames(geneAgeDf, seq_len(2), c("geneID", "cath"))  #col names
    row.names(geneAgeDf) <- NULL   # remove row names

    ### move NA cat to working taxonomy rank
    if (rankName == "genus")
        geneAgeDf$cath[is.na(geneAgeDf$cath)] <- "011111111"
    else if (rankName == "family")
        geneAgeDf$cath[is.na(geneAgeDf$cath)] <- "001111111"
    else if (rankName == "class")
        geneAgeDf$cath[is.na(geneAgeDf$cath)] <- "000111111"
    else if (rankName == "phylum")
        geneAgeDf$cath[is.na(geneAgeDf$cath)] <- "000011111"
    else if (rankName == "kingdom")
        geneAgeDf$cath[is.na(geneAgeDf$cath)] <- "000001111"
    else if (rankName == "superkingdom")
        geneAgeDf$cath[is.na(geneAgeDf$cath)] <- "000000011"
    else
        geneAgeDf$cath[is.na(geneAgeDf$cath)] <- "111111111"
    ### merge ranks that are lower than ref rank to ref rank
    index <- 0
    if (rankName != "species")
        index <- match(rankName, rankList)
    maxCat <- paste0(
        paste(rep(0,index), collapse = ""), paste(rep(1,9-index), collapse = "")
    )
    geneAgeDf$cath[geneAgeDf$cath > maxCat] <- maxCat
    ### convert cat into geneAge
    domainDfList <- lapply(
        geneAgeDf[geneAgeDf$cath == "000000001",]$geneID,
        function (x) {
            orthoList <- mdDtExt[
                mdDtExt$geneID == x & mdDtExt$cath == "000000001",]$ncbiID
            orthoDomainID <- levels(as.factor(
                subDt[subDt$abbrName %in% orthoList,]$superkingdom))
            orthoDomainName <- levels(as.factor(c(
                taxList$fullName[taxList$ncbiID == supFirstLine$superkingdom],
                taxList$fullName[taxList$ncbiID %in% orthoDomainID]
            )))
            if (length(orthoDomainName) == 3) {
                age <- paste0("10_", paste(orthoDomainName, collapse = "-"))
            } else {
                age <- paste0("09_", paste(orthoDomainName, collapse = "-"))
            }
            return(data.frame(geneID = x, age, stringsAsFactors = FALSE))
        }
    )
    konList <- lapply(
        geneAgeDf[geneAgeDf$cath == "000000111",]$geneID,
        function (x) {
            orthoList <- mdDtExt[
                mdDtExt$geneID == x & mdDtExt$cath == "000000111",]$ncbiID
            orthoKonID <- levels(as.factor(
                subDt[subDt$abbrName %in% orthoList,]$norank_33154))
            superKingdom <- taxList$ncbiID[
                taxList$ncbiID == supFirstLine$superkingdom
                & taxList$rank == "superkingdom"]
            if (superKingdom == 2759) {
                orthoKonID <- unique(c(orthoKonID, supFirstLine$norank_33154))
                if (length(orthoKonID) == 1) {
                    age <- paste0(
                        "06_", taxList$fullName[taxList$ncbiID == orthoKonID])
                } else {
                    ukon <- levels(as.factor(orthoKonID%in%c(33154, 554915)))
                    if (length(ukon) == 1) {
                        if (ukon == TRUE) {
                            age <- "07_Unikonta"
                        } else {
                            age <- "07_Bikonta"
                        }
                    } else {
                        age <- "08_Eukaryote"
                    }
                }
            } else {
                age <- paste0(
                    "08_", taxList$fullName[
                        taxList$ncbiID == supFirstLine$superkingdom])
            }
            return(data.frame(geneID = x, age, stringsAsFactors = FALSE))
        }
    )
    sarList <- lapply(
        geneAgeDf[geneAgeDf$cath == "000001111",]$geneID,
        function (x) {
            if (
                taxList$fullName[taxList$ncbiID == supFirstLine$kingdom] =="Sar"
            ) {
                orthoList <- mdDtExt[
                    mdDtExt$geneID == x & mdDtExt$cath == "000001111",]$ncbiID
                sarTaxDt <- Dt[
                    Dt$abbrName %in% orthoList,
                    colnames(Dt)[
                        colnames(Dt) %in%
                            c("norank_33630", "norank_543769", "norank_33634")]]
                check <- lapply(
                    sarTaxDt,
                    function(x) {
                        if (length(unique(x)) == 1) return(unique(x))
                        else return(NULL)
                    }
                )
                if (is.null(unlist(check))) {
                    age <- paste0(
                        "05_",
                        taxList$fullName[taxList$ncbiID == supFirstLine$kingdom]
                    )
                } else {
                    age <- paste0(
                        "05_",
                        taxList$fullName[taxList$ncbiID == unlist(check)])
                }
            } else {
                age <- paste0(
                    "05_",
                    taxList$fullName[taxList$ncbiID == supFirstLine$kingdom])
            }
            return(data.frame(geneID = x, age, stringsAsFactors = FALSE))
        }
    )
    geneAgeDfPre <- do.call(rbind, c(domainDfList, konList, sarList))
    if (!(is.null(geneAgeDfPre)))
        geneAgeDf <- merge(geneAgeDf, geneAgeDfPre, by = "geneID", all.x = TRUE)
    geneAgeDf$age[geneAgeDf$cath == "000000011"] <- paste0(
        "08_", taxList$fullName[taxList$ncbiID == supFirstLine$superkingdom])
    geneAgeDf$age[geneAgeDf$cath == "000011111"] <- paste0(
        "04_", taxList$fullName[taxList$ncbiID == supFirstLine$phylum])
    geneAgeDf$age[geneAgeDf$cath == "000111111"] <- paste0(
        "03_", taxList$fullName[taxList$ncbiID == supFirstLine$class])
    geneAgeDf$age[geneAgeDf$cath == "001111111"] <- paste0(
        "02_", taxList$fullName[taxList$ncbiID == supFirstLine$family])
    geneAgeDf$age[geneAgeDf$cath == "011111111"] <- paste0(
        "01_", taxList$fullName[taxList$ncbiID == supFirstLine$genus])
    geneAgeDf$age[geneAgeDf$cath == "111111111"] <- paste0(
        "00_", taxList$fullName[
            taxList$fullName == refTaxon & taxList$rank == rankName])
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
#' age = c("07_LUCA", "06_Eukaryota", "07_LUCA", "06_Eukaryota")
#' )
#' geneAgePlotDf(geneAgeDf)

geneAgePlotDf <- function(geneAgeDf){
    if (is.null(geneAgeDf)) stop("Gene age data is NULL!")
    age <- NULL
    plotDf <- geneAgeDf %>% dplyr::count(age)
    plotDf$level <- substring(plotDf$age, 1,2)
    levelDf <- data.frame(
        level = c("00", "01", "02", "03", "04", "05", "08", "10"),
        rank = c(
            " (Species)", " (Genus)", " (Family)", " (Class)", " (Phylum)",
            " (Kingdom)", " (Superkingdom)"," (LUCA)"
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
