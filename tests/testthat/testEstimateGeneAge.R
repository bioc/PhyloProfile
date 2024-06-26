context("test estimation of gene ages")

test_that("test estimation of gene ages", {
    library(dplyr)
    data("fullProcessedProfile", package="PhyloProfile")
    rankName <- "class"
    refTaxon <- "Mammalia"
    var1Cutoff <- c(0,1)
    var2Cutoff <- c(0,1)
    percentCutoff <- c(0,1)
    taxonIDs <- levels(as.factor(fullProcessedProfile$ncbiID))
    sortedInputTaxa <- sortInputTaxa(
        taxonIDs, rankName, refTaxon, NULL
    )
    taxaCount <- sortedInputTaxa %>% dplyr::count(supertaxon)
    geneAge <- estimateGeneAge(
        fullProcessedProfile, taxaCount,
        rankName, refTaxon,
        var1Cutoff, var2Cutoff, percentCutoff
    )
    expect_true(
        geneAge$age[geneAge$geneID == "103479at6656"] == "08_Eukaryota"
    )
})

test_that("test plotting gene age plot", {
    geneAgeDf <- data.frame(
        geneID = c("OG_1017", "OG_1019"),
        cat = c("0000001", "0000001"),
        age = c("07_LUCA", "07_LUCA"),
        stringsAsFactors = FALSE
    )
    plotDf <- geneAgePlotDf(geneAgeDf)
    geneAgeText <- 1
    p <- createGeneAgePlot(plotDf, geneAgeText)
    expect_true(nrow(p$data) == 1)
})
