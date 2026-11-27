#' Process ortholog IDs
#' @description Process ortholog IDs to identify duplicated IDs
#' @param dataHeat a data frame contains processed profiles (see
#' ?fullProcessedProfile, ?filterProfileData)
#' @return the same dataframe as input, but the ortholog IDs are changed into
#' <taxID:orthoID>. New column "orthoFreq" specifies if the ortholog IDs are
#' single or duplicated
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @export
#' @examples
#' ?processOrthoID
#' data("finalProcessedProfile", package="PhyloProfile")
#' processOrthoID(finalProcessedProfile)

processOrthoID <- function(dataHeat = NULL) {
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    orthoID <- geneID <- orthoFreqCount <- orthoFreqNew <- supertaxonID <- NULL
    orthoIDNew <- orthoFreq <- NULL
    # Check if idFormat is "bionf"
    idFormat <- "other"
    firstOrtho <- strsplit(dataHeat$orthoID[1],'\\|')[[1]]
    firstGeneID <- dataHeat$geneID[1]
    firstSupertaxonID <- dataHeat$supertaxonID[1]
    if (
        length(firstOrtho) >= 3 && firstOrtho[1] == firstGeneID &&
        grepl(firstSupertaxonID, firstOrtho[2])
    ) { idFormat <- "bionf" }

    # parse orthoID
    if (idFormat == "bionf") {
        dataHeat <- within(
            dataHeat,
            orthoMod <- data.frame(
                do.call('rbind', strsplit(orthoID,'|',fixed=TRUE))
            )
        )
        dataHeat$orthoIDNew <- paste(
            dataHeat$orthoMod$X2, dataHeat$orthoMod$X3, sep = "#"
        )
    } else {
        dataHeat$orthoIDNew <- paste(
            dataHeat$supertaxonID, dataHeat$orthoID, sep = "#"
        )
    }
    dataHeat <- dataHeat[ , !(names(dataHeat) %in% ("orthoMod"))]
    # Count occurrences of orthoIDNew and label as "Single" or "Multiple"
    countOrthoDf <- dataHeat %>%
        count(orthoIDNew, name = "orthoFreq") %>%
        mutate(orthoFreq = ifelse(orthoFreq > 1, "Multiple", "Single"))
    # Merge the frequency labels back to the main data
    dataHeat <- dataHeat %>%
        left_join(countOrthoDf, by = "orthoIDNew")
    # Check for "Multiple" freq in co-orthologs at geneID and supertaxonID level
    freqDt <- dataHeat %>%
        distinct(geneID, supertaxonID, orthoFreq) %>%
        group_by(geneID, supertaxonID) %>%
        mutate(
            orthoFreq = ifelse("Multiple" %in% orthoFreq, "Multiple", orthoFreq)
        ) %>% ungroup()
    # Merge the updated frequency labels back to the main data
    dataHeat <- dataHeat %>% left_join(
        freqDt, by = c("geneID", "supertaxonID"), suffix = c("", "New"),
        relationship = "many-to-many"
    )
    # Step 3: Update orthoFreq with the new frequency labels and clean up
    dataHeat <- dataHeat %>% mutate(orthoFreq = orthoFreqNew) %>%
        select(-orthoFreqNew)
    return(dataHeat)
}

#' Create data for main profile plot
#' @export
#' @param dataHeat a data frame contains processed profiles (see
#' ?fullProcessedProfile, ?filterProfileData)
#' @return A dataframe for plotting the phylogenetic profile, containing seed
#' protein IDs (geneID), ortholog IDs (orthoID) together with their ncbi
#' taxonomy IDs (ncbiID and abbrName), full names (fullName), indexed supertaxa
#' (supertaxon), values for additional variables (var1, var2) and the aggregated
#' values of those additional variables for each supertaxon (mVar1, mVar2),
#' number of original and filtered co-orthologs in each supertaxon (paralog and
#' paralogNew), number of species in each supertaxon (numberSpec) and the % of
#' species that have orthologs in each supertaxon (presSpec).
#' @importFrom stats na.omit
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{filterProfileData}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' dataMainPlot(finalProcessedProfile)

dataMainPlot <- function(dataHeat = NULL){
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    paralogNew <- orthoID <- NULL
    # reduce number of inparalogs based on filtered dataHeat
    dataHeatTb <- data.table(stats::na.omit(dataHeat))
    dataHeatTb[, paralogNew := .N, by = c("geneID", "supertaxon")]
    dataHeatTb <- data.frame(dataHeatTb[, c(
        "geneID", "supertaxon", "paralogNew"
    )])

    dataHeat <- dataHeat %>% left_join(dataHeatTb, by=c("geneID","supertaxon"))
    dataHeat$paralog <- dataHeat$paralogNew
    dataHeat <- dataHeat[!duplicated(dataHeat), ]

    # remove unneeded dots
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    dataHeat$paralog[dataHeat$presSpec < 1] <- NA
    dataHeat$paralog[dataHeat$paralog == 1] <- NA
    # rescale numbers of paralogs
    if (length(unique(stats::na.omit(dataHeat$paralog))) > 0) {
        maxParalog <- max(stats::na.omit(dataHeat$paralog))
        dataHeat$paralogSize <- (dataHeat$paralog / maxParalog) * 3
    }
    # remove prefix number of taxa names but keep the order
    dataHeat$supertaxon <- factor(
        substr(
            as.character(dataHeat$supertaxon), 8 ,
            nchar(as.character(dataHeat$supertaxon))),
        levels = substr(
            levels(as.factor(dataHeat$supertaxon)), 8,
            nchar(levels(as.factor(dataHeat$supertaxon)))))
    # count ortholog IDs for each taxon
    dataHeat <- processOrthoID(dataHeat)
    return(dataHeat)
}

#' Create data for customized profile plot
#' @description Create data for customized profile plot based on a selected
#' list of genes and/or taxa, containing seed protein IDs (geneID), ortholog IDs
#' (orthoID) together with their ncbi taxonomy IDs (ncbiID and abbrName), full
#' names (fullName), indexed supertaxa (supertaxon), values for additional
#' variables (var1, var2) and the aggregated values of those additional
#' variables for each supertaxon (mVar1, mVar2), number of original and filtered
#' co-orthologs in each supertaxon (paralog and paralogNew), number of species
#' in each supertaxon (numberSpec) and the % of species that have orthologs in
#' each supertaxon (presSpec).
#' @export
#' @usage dataCustomizedPlot(dataHeat = NULL, selectedTaxa = "all",
#'     selectedSeq = "all")
#' @param dataHeat a data frame contains processed profiles (see
#' ?fullProcessedProfile, ?filterProfileData)
#' @param selectedTaxa selected subset of taxa. Default = "all".
#' @param selectedSeq selected subset of genes. Default = "all".
#' @return A dataframe contains data for plotting the customized profile.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom stats na.omit
#' @seealso \code{\link{filterProfileData}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' selectedTaxa <- c("Mammalia", "Saccharomycetes", "Insecta")
#' selectedSeq <- "all"
#' dataCustomizedPlot(finalProcessedProfile, selectedTaxa, selectedSeq)

dataCustomizedPlot <- function(
    dataHeat = NULL, selectedTaxa = "all", selectedSeq = "all"
){
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    geneID <- supertaxonMod <- paralogNew <- NULL
    dataHeat$supertaxonMod <- {
        substr(
            dataHeat$supertaxon, 8, nchar(as.character(dataHeat$supertaxon))
        )
    }
    if (selectedTaxa[1] == "all" & selectedSeq[1] != "all") {
        dataHeat <- subset(dataHeat, geneID %in% selectedSeq)
    } else if (selectedSeq[1] == "all" & selectedTaxa[1] != "all") {
        dataHeat <- subset(dataHeat, supertaxonMod %in% selectedTaxa)
    } else {
        dataHeat <- subset(
            dataHeat, geneID %in% selectedSeq & supertaxonMod %in% selectedTaxa
        )
    }
    # reduce number of inparalogs based on filtered dataHeat
    dataHeatTb <- data.table(stats::na.omit(dataHeat))
    dataHeatTb[, paralogNew := .N, by = c("geneID", "supertaxon")]
    dataHeatTb <- data.frame(
        dataHeatTb[, c("geneID", "supertaxon", "paralogNew")]
    )
    dataHeat <- merge(
        dataHeat, dataHeatTb, by = c("geneID", "supertaxon"), all.x = TRUE
    )
    dataHeat$paralog <- dataHeat$paralogNew
    dataHeat <- dataHeat[!duplicated(dataHeat), ]
    # remove unneeded dots
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    dataHeat$paralog[dataHeat$presSpec < 1] <- NA
    dataHeat$paralog[dataHeat$paralog == 1] <- NA
    # rescale numbers of paralogs
    if (length(unique(stats::na.omit(dataHeat$paralog))) > 0) {
        maxParalog <- max(stats::na.omit(dataHeat$paralog))
        dataHeat$paralogSize <- (dataHeat$paralog / maxParalog) * 3
    }
    # remove prefix number of taxa names but keep the order
    dataHeat$supertaxon <- factor(
        substr(
            as.character(dataHeat$supertaxon), 8,
            nchar(as.character(dataHeat$supertaxon))),
        levels = substr(
            levels(dataHeat$supertaxon), 8,
            nchar(levels(dataHeat$supertaxon))))
    # count ortholog IDs for each taxon
    if (nrow(dataHeat) == 0) return(dataHeat)
    dataHeat <- processOrthoID(dataHeat)
    return(dataHeat)
}

#' Create profile heatmap plot
#' @export
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param parm plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2) display gene IDs (default) or gene names;
#' (3+4) names of 2 variables var1ID and var2ID - default = "var1" & "var2";
#' (5+6) mid value and color for mid value of var1 -
#' default is 0.5 and #FFFFFF; (7) color for lowest var1 - default = "#FF8C00";
#' (8) color for highest var1 - default = "#4682B4"; (9+10) mid value and color
#' for mid value of var2 - default is 1 and #FFFFFF;(11) color for lowest var2 -
#' default = "#FFFFFF", (12) color for highest var2 - default = "#F0E68C", (13)
#' color of co-orthologs - default = "#07D000"; (14+15+16) text sizes for x, y
#' axis and legend - default = 9 for each; (17) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (18) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (19) angle of x-axis from 0 to
#' 90 - default = 60; (20) show/hide separate line for reference taxon 1/0 -
#' default = 0; (21) enable/disable coloring gene categories TRUE/FALSE -
#' default = FALSE; (22) enable/disable coloring duplicated ortholog IDs
#' TRUE/FALSE - default=FALSE). NOTE: Leave blank or NULL to use default values.
#' @return A profile heatmap plot as a ggplot object.
#' @import ggplot2
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE,
#'     "catColors" = NULL,
#'     "colorByOrthoID" = FALSE
#' )
#'
#' heatmapPlotting(plotDf, plotParameter)

heatmapPlotting <- function(data = NULL, parm = NULL){
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (is.null(parm))
        parm <- list(
            "xAxis" = "taxa", "geneIdType" = "geneName",
            "var1ID" = "var1", "var2ID"  = "var2",
            "midVar1" = 0.5, "midColorVar1" =  "#FFFFFF",
            "lowColorVar1" =  "#FF8C00", "highColorVar1" = "#4682B4",
            "midVar2" = 1, "midColorVar2" =  "#FFFFFF",
            "lowColorVar2" = "#CB4C4E", "highColorVar2" = "#3E436F",
            "paraColor" = "#07D000", "xSize" = 8, "ySize" = 8, "legendSize" = 8,
            "mainLegend" = "top", "dotZoom" = 0, "xAngle" = 60, "guideline" = 0,
            "colorByGroup" = FALSE,"catColors" = NULL,"colorByOrthoID" = FALSE)
    geneID<-geneName<-supertaxon<-category<-var1<-var2<- presSpec<-paralog<-NULL
    orthoFreq <- xmin <- xmax <- ymin <- ymax <- NULL
    ### create heatmap plot
    if (!(is.null(parm$geneIdType))) {
        if (parm$geneIdType == "geneName") data$geneID <- data$geneName
    }
    # create a canvas
    if (parm$xAxis == "genes") {
        p <- ggplot(data,aes(x=geneID, y=supertaxon)) +
            labs(x = "Gene ID", y = "Taxon")
    } else {
        p <- ggplot(data, aes(y = geneID, x = supertaxon)) +
            labs(y = "Gene ID", x = "Taxon")
    }
    # create geom_tile & scale_fill_gradient for var2 OR gene category
    if (parm$colorByGroup == TRUE) {
        p <- p + geom_raster(aes(fill = factor(category)), alpha = 0.3)
        if (!is.null(parm$catColors))
            p <- p + scale_fill_manual(values = parm$catColors)
    } else {
        if (length(unique(stats::na.omit(data$var2))) != 1)
            p <- p + scale_fill_gradient2(
                low = parm$lowColorVar2, high = parm$highColorVar2,
                mid = parm$midColorVar2, midpoint = parm$midVar2,
                na.value = "gray95", limits = c(0, 1)) +
                geom_raster(aes(fill = var2))
    }

    # create geom_point for found ortho; coloring by var1 or
    # orthoIDs (only when working on lowest taxonomy rank)
    # and inparalogs (co-orthologs)
    if (length(unique(stats::na.omit(data$presSpec))) < 2) {
        # working on the lowest taxonomy rank
        if (parm$colorByOrthoID == TRUE) {
            # color by ortho IDs
            p <- p + geom_point(
                aes(colour = factor(orthoFreq)), na.rm = TRUE,
                size = data$presSpec*5*(1+parm$dotZoom), show.legend = TRUE
            )
        } else {
            if (length(unique(stats::na.omit(data$var1))) == 1) {
                p <- p + geom_point(
                    aes(colour = var1), na.rm = TRUE,
                    size = data$presSpec*5*(1+parm$dotZoom), show.legend = FALSE
                )
            } else {
                p <- p +
                    geom_point(
                        aes(colour = var1), na.rm = TRUE,
                        size = data$presSpec * 5 * (1 + parm$dotZoom)) +
                    scale_color_gradient2(
                        low = parm$lowColorVar1, high = parm$highColorVar1,
                        mid = parm$midColorVar1, midpoint = parm$midVar1,
                        limits = c(0,1)
                    )
            }
            # plot inparalogs (if available)
            if (length(unique(stats::na.omit(data$paralog))) > 0) {
                # Define min and max for your desired size range
                min_size <- 1
                max_size <- 3*(1+parm$dotZoom)
                # Manually rescale paralog size
                min_val <- min(data$paralog, na.rm = TRUE)
                max_val <- max(data$paralog, na.rm = TRUE)
                # Avoid division by zero (when min == max)
                if (min_val == max_val) {
                    data$paralogSize <- max_size
                } else {
                    data$paralogSize <-
                        ((data$paralog - min_val) / (max_val - min_val)) *
                        (max_size - min_size) + min_size
                }
                # create breaks for legend
                paralog_vals <- sort(unique(na.omit(data$paralog)))
                n_vals <- length(paralog_vals)
                custom_breaks <- if (n_vals <= 5) {
                    paralog_vals
                } else {
                    paralog_vals[round(seq(2, n_vals, length.out = 5))]
                }
                p <- p +
                    geom_point(
                        data = data, aes(size = paralog), color=parm$paraColor,
                        na.rm = TRUE, show.legend = TRUE) +
                    guides(size = guide_legend(title = "# of co-orthologs")) +
                scale_size_continuous(
                    breaks = custom_breaks, range = c(min_size, max_size))
            }
        }
    } else {
        # scale dot size based on % preSpec in each super taxon
        if (length(unique(stats::na.omit(data$var1))) == 1) {
            p <- p + geom_point(aes(size=presSpec),color="#336a98",na.rm = TRUE)
        } else {
            p <- p + geom_point(aes(colour=var1, size = presSpec),na.rm = TRUE)+
                scale_color_gradient2(
                    low = parm$lowColorVar1, high = parm$highColorVar1,
                    mid = parm$midColorVar1, midpoint = parm$midVar1,
                    limits = c(0,1)
                )
        }
        # remain the scale of point while filtering % presSpec
        presentVl <- data$presSpec[!is.na(data$presSpec)]
        p <- p + scale_size_continuous(range = c(
            (floor(min(presentVl) * 10) / 10 * 5) * (1 + parm$dotZoom),
            (floor(max(presentVl) * 10) / 10 * 5) * (1 + parm$dotZoom)))
    }

    # create legend
    if (parm$colorByGroup == FALSE) {
        if (parm$colorByOrthoID == FALSE) {
            p <- p + guides(fill = guide_colourbar(title = parm$var2ID),
                            color = guide_colourbar(title = parm$var1ID))
        } else {
            p <- p + guides(fill = guide_colourbar(title = parm$var2ID),
                            color = guide_legend(title = 'OrthoID copy number'))
        }
    } else {
        p <- p + guides(fill = guide_legend("Category"),
                        color = guide_colourbar(title = parm$var1ID))
    }

    # text size, legend position
    p <- p + theme_minimal(base_size = 9)
    vjustValue <- 1
    if (parm$xAngle == 90) vjustValue <- 0.5
    p <- p + theme(
        axis.text.x = element_text(
            angle = parm$xAngle, hjust = 1, size = parm$xSize,
            vjust = vjustValue
        ),
        axis.text.y = element_text(size = parm$ySize),
        axis.title.x = element_text(size = parm$xSize),
        axis.title.y = element_text(size = parm$ySize),
        legend.title = element_text(size = parm$legendSize),
        legend.text = element_text(size = parm$legendSize),
        legend.position = parm$mainLegend)
    return(p)
}


#' Create profile heatmap plot using scattermore
#' @export
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param parm plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2) display gene IDs (default) or gene names;
#' (3+4) names of 2 variables var1ID and var2ID - default = "var1" & "var2";
#' (5+6) mid value and color for mid value of var1 -
#' default is 0.5 and #FFFFFF; (7) color for lowest var1 - default = "#FF8C00";
#' (8) color for highest var1 - default = "#4682B4"; (9+10) mid value and color
#' for mid value of var2 - default is 1 and #FFFFFF;(11) color for lowest var2 -
#' default = "#FFFFFF", (12) color for highest var2 - default = "#F0E68C", (13)
#' color of co-orthologs - default = "#07D000"; (14+15+16) text sizes for x, y
#' axis and legend - default = 9 for each; (17) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (18) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (19) color dots based on either
#' "var1" or "var2". NOTE: Leave blank or NULL to
#' use default values.
#' @return A profile heatmap plot as a ggplot object.
#' @import ggplot2
#' @import scattermore
#' @importFrom grDevices colorRampPalette
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "colorVar" = "var1"
#' )
#'
#' heatmapPlottingFast(plotDf, plotParameter)

heatmapPlottingFast <- function(data = NULL, parm = NULL) {
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (is.null(parm))
        parm <- list(
            "xAxis" = "taxa", "geneIdType" = "geneName",
            "var1ID" = "var1", "var2ID"  = "var2",
            "midVar1" = 0.5, "midColorVar1" =  "#FFFFFF",
            "lowColorVar1" =  "#FF8C00", "highColorVar1" = "#4682B4",
            "midVar2" = 1, "midColorVar2" =  "#FFFFFF",
            "lowColorVar2" = "#CB4C4E", "highColorVar2" = "#3E436F",
            "paraColor" = "#07D000", "xSize" = 8, "ySize" = 8, "legendSize" = 8,
            "mainLegend" = "top", "dotZoom" = 0, "colorVar" = "var1")
    geneID<-geneName<-supertaxon<-category<-var1<-var2<- presSpec<-paralog<-NULL
    orthoFreq <- xmin <- xmax <- ymin <- ymax <- x <- y <- NULL
    if (is.null(parm$colorVar)) parm$colorVar <- "var1"
    ### create heatmap plot
    if (!(is.null(parm$geneIdType))) {
        if (parm$geneIdType == "geneName") data$geneID <- data$geneName
    }
    # create a canvas
    if (parm$xAxis == "genes") {
        p <- ggplot(data,aes(x = geneID, y = supertaxon)) +
            labs(x = "Gene ID", y = "Taxon")
    } else {
        p <- ggplot(data, aes(y = geneID, x = supertaxon)) +
            labs(y = "Gene ID", x = "Taxon")
    }
    # format theme
    p <- p + theme_minimal(base_size = 9)
    p <- p + theme(
        axis.text = element_blank(),
        axis.title.x = element_text(size = parm$xSize),
        axis.title.y = element_text(size = parm$ySize),
        legend.title = element_text(size = parm$legendSize),
        legend.text = element_text(size = parm$legendSize),
        legend.position = parm$mainLegend)
    # Extract x, y, and color values (swap when genes are on the x-axis)
    taxIdx  <- as.numeric(data$supertaxon)
    geneIdx <- as.numeric(as.factor(data$geneID))
    if (parm$xAxis == "genes") {
        x_values <- geneIdx
        y_values <- taxIdx
    } else {
        x_values <- taxIdx
        y_values <- geneIdx
    }
    # Create color vector: green if paralog > 1, else gradient by var1 or var2
    if (parm$colorVar == "var1") {
        cl <- colorRampPalette(c(parm$lowColorVar1, parm$highColorVar1))(100)[
            cut(data$var1, breaks = 100)
        ]
    } else
        cl <- colorRampPalette(c(parm$lowColorVar2, parm$highColorVar2))(100)[
            cut(data$var2, breaks = 100)
        ]
    cl[data$paralog > 1] <- parm$paraColor
    # Add scatter plot. The raster resolution follows the plot size (if
    # given), so that dots are not squeezed when the plot is small
    rasterPx <- c(1000, 1000)
    if (!is.null(parm$width) && !is.null(parm$height) &&
        isTRUE(parm$width > 0) && isTRUE(parm$height > 0))
        rasterPx <- round(pmin(c(parm$width, parm$height), 2000))
    p <- p +
        geom_scattermost(
            cbind(x_values, y_values), color = cl,
            pointsize = 2 * (1 + parm$dotZoom), pixels = rasterPx
        )
    p <- p +
        scale_x_continuous(
            limits = c(0.5, max(x_values, na.rm = TRUE) + 0.5),
            expand = c(0, 0)
        ) +
        scale_y_continuous(
            limits = c(0.5, max(y_values, na.rm = TRUE) + 0.5),
            expand = c(0, 0)
        )
    if (parm$colorVar == "var1") {
        p <- p +
            geom_point(data=data.frame(x = double(0)), aes(x, x, color = x)) +
            scale_color_gradientn(
                colors = colorRampPalette(
                    c(parm$lowColorVar1, parm$highColorVar1))(100),
                limits = c(0,1), name = parm$var1ID
            )
    } else {
        p <- p +
            geom_point(data=data.frame(x = double(0)), aes(x, x, color = x)) +
            scale_color_gradientn(
                colors = colorRampPalette(
                    c(parm$lowColorVar2, parm$highColorVar2))(100),
                limits = c(0,1), name = parm$var2ID
            )
    }
    # Add only one real co-ortholog point with fill = green to generate legend
    if(any(data$paralog > 1, na.rm = TRUE)) {
        co_point <- data[data$paralog > 1, ][1, ]
        coIdx <- which(data$paralog > 1)[1]
        co_point$x <- x_values[coIdx]
        co_point$y <- y_values[coIdx]
        p <- p + geom_point(
            data = co_point, aes(x = x, y = y, fill = "Co-orthologs"),
            shape = 21, size = 3
        ) + scale_fill_manual(
            values = c("Co-orthologs" = parm$paraColor), name = NULL
        )
    }
    return(p)
}

#' Highlight gene and/or taxon of interest on the phylogenetic profile plot
#' @export
#' @usage highlightProfilePlot(profilePlot = NULL, plotDf = NULL,
#'     taxonHighlight = "none", workingRank = "none", geneHighlight = NULL,
#'     taxDB = NULL, xAxis = "taxa")
#' @param profilePlot initial (highlighted) profile plot
#' @param plotDf dataframe for plotting the heatmap phylogentic profile
#' @param taxonHighlight taxon of interst. Default = "none".
#' @param workingRank working taxonomy rank (needed only for highlight taxon).
#' @param geneHighlight gene of interest. Default = NULL.
#' @param taxDB Path to the taxonomy DB files
#' @param xAxis type of x-axis (either "genes" or "taxa")
#' @return A profile heatmap plot with highlighted gene and/or taxon of interest
#' as ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}},
#' \code{\link{heatmapPlotting}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE,
#'     "colorByOrthoID" = FALSE
#' )
#' profilePlot <- heatmapPlotting(plotDf, plotParameter)
#' taxonHighlight <- "none"
#' workingRank <- "class"
#' geneHighlight <- "100265at6656"
#' highlightProfilePlot(
#'     profilePlot, plotDf, taxonHighlight, workingRank, geneHighlight,
#'     NULL, plotParameter$xAxis
#' )

highlightProfilePlot <- function(
        profilePlot = NULL, plotDf = NULL, taxonHighlight = "none",
        workingRank = "none", geneHighlight = NULL, taxDB = NULL, xAxis = "taxa"
){
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    if (is.null(profilePlot)) stop("Profile plot cannot be NULL!")
    xmin <- xmax <- ymin <- ymax <- NULL
    if (is.null(taxonHighlight) && is.null(geneHighlight)) return(profilePlot)
    # highlight taxon
    if (length(taxonHighlight) > 0 && !("none" %in% taxonHighlight)) {
        # get selected highlight taxon ID
        taxName <- getNameList(taxDB)
        taxonHighlightID <- taxName$ncbiID[
            taxName$fullName %in% taxonHighlight & taxName$rank == workingRank]
        if (length(taxonHighlightID) == 0L) {
            taxonHighlightID <- taxName$ncbiID[
                taxName$fullName %in% taxonHighlight
            ]
        } else if(length(taxonHighlightID) < length(taxonHighlight)) {
            missingID <- taxName$ncbiID[
                taxName$fullName %in% taxonHighlight
            ]
            missingID <- missingID[!(missingID %in% taxonHighlightID)]
            taxonHighlightID <- c(taxonHighlightID, missingID)
        }
        # get taxonID together with it sorted index
        taxonHighlightID <- intersect(taxonHighlightID, plotDf$supertaxonID)
        if (length(taxonHighlightID) > 0) {
            selTaxon <- unique(as.character(
                plotDf$supertaxon[plotDf$supertaxonID %in% taxonHighlightID]
            ))
            selIndex <- match(selTaxon, levels(as.factor(plotDf$supertaxon)))
            if (xAxis == "taxa") {
                rect <- data.frame(
                    xmin=selIndex-0.5, xmax=selIndex+0.5, ymin=-Inf, ymax=Inf)
            } else
                rect <- data.frame(
                    ymin=selIndex-0.5, ymax=selIndex+0.5, xmin=-Inf, xmax=Inf)
            profilePlot <- profilePlot + geom_rect(
                data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
        }
    }
    # highlight gene
    if (length(geneHighlight) > 0 && !("none" %in% geneHighlight)) {
        geneHighlight <- intersect(geneHighlight, unique(plotDf$geneID))
        if (length(geneHighlight) > 0) {
            selIndex <- match(geneHighlight, levels(as.factor(plotDf$geneID)))
            if (xAxis == "taxa") {
                rect <- data.frame(
                    ymin=selIndex-0.5, ymax=selIndex+0.5, xmin=-Inf, xmax=Inf)
            } else
                rect <- data.frame(
                    xmin=selIndex-0.5, xmax=selIndex+0.5, ymin=-Inf, ymax=Inf)
            profilePlot <- profilePlot + geom_rect(
                data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
        }
    }
    return(profilePlot)
}

#' Add taxonomy rank division lines to the heatmap plot
#' @export
#' @usage addRankDivisionPlot(profilePlot = NULL, plotDf = NULL,
#'     taxDB = NULL, workingRank = NULL, superRank = NULL, xAxis = "taxa",
#'     font = "Arial", groupLabelSize = 14, groupLabelDist = 2,
#'     groupLabelAngle = 90, refLine = TRUE, lineColor = "dodgerblue4",
#'     plotSize = NULL, labelMode = "lca")
#' @param profilePlot initial (highlighted) profile plot
#' @param plotDf dataframe for plotting the heatmap phylogentic profile
#' @param taxDB path to taxonomy database (taxonomyMatrix.txt file required!)
#' @param workingRank working taxonomy rank (e.g. species)
#' @param superRank taxonomy rank for division lines (e.g. superkingdom)
#' @param xAxis type of x-axis (either "genes" or "taxa")
#' @param font font of text. Default = Arial"
#' @param groupLabelSize size of rank labels
#' @param groupLabelDist size of the plot area for rank labels
#' @param groupLabelAngle angle of rank labels
#' @param refLine add vertical line to separate reference taxon
#' @param lineColor color of the division lines. Default = "dodgerblue4"
#' @param plotSize plot size in px along the division axis (width if xAxis =
#' "taxa", height if "genes"). Used to merge overlapping labels into their
#' common ancestor (marked with "*"). NULL = never merge
#' @param labelMode how to handle overlapping labels: "lca" (merge them
#' into their lowest common ancestor, marked with "*") or "first" (show
#' only the first label of the overlapping ones, followed by ",...").
#' Default = "lca"
#' @return A profile heatmap plot with highlighted gene and/or taxon of interest
#' as ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{heatmapPlotting}}, \code{\link{highlightProfilePlot}},
#' \code{\link{getTaxonomyMatrix}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE,
#'     "colorByOrthoID" = FALSE
#' )
#' profilePlot <- heatmapPlotting(plotDf, plotParameter)
#' workingRank <- "class"
#' superRank <- "superkingdom"
#' addRankDivisionPlot(
#'     profilePlot, plotDf, NULL, workingRank, superRank, "taxa", font = "sans"
#' )

addRankDivisionPlot <- function(
        profilePlot = NULL, plotDf = NULL, taxDB = NULL,
        workingRank = NULL, superRank = NULL, xAxis = "taxa", font = "Arial",
        groupLabelSize = 14, groupLabelDist = 2, groupLabelAngle = 90,
        refLine = TRUE, lineColor = "dodgerblue4", plotSize = NULL,
        labelMode = "lca"
) {
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    if (is.null(profilePlot)) stop("Profile plot cannot be NULL!")
    if (is.null(lineColor) || lineColor == "") lineColor <- "dodgerblue4"
    
    profilePlot <- profilePlot + theme(text = element_text(family = font))
    if (is.null(workingRank) || is.null(superRank)) {
        if (refLine == TRUE) {
            # guideline for separating ref species
            if (xAxis == "genes") {
                profilePlot <- profilePlot + labs(y = "Taxon") +
                    geom_hline(yintercept = 0.5, colour = lineColor) +
                    geom_hline(yintercept = 1.5, colour = lineColor)
            } else
                profilePlot <- profilePlot + labs(x = "Taxon") +
                    geom_vline(xintercept = 0.5, colour = lineColor) +
                    geom_vline(xintercept = 1.5, colour = lineColor)
        }
        return(profilePlot)
    } else {
        # sort supertaxonID based on the sorted supertaxon
        plotDf <- plotDf[with(plotDf,order(supertaxon)),]
        plotDf$supertaxonID <- factor(
            plotDf$supertaxonID, levels = unique(plotDf$supertaxonID)
        )
        # get input (super)taxa
        inputTax <- levels(as.factor(plotDf$supertaxonID))
        # subset taxonomy matrix to contain only superrank for division line
        taxMatrix <- getTaxonomyMatrix(taxDB)
        subTaxMatrix <- subset(
            taxMatrix[taxMatrix[[workingRank]] %in% inputTax,],
            select = c(as.character(workingRank), as.character(superRank))
        )
        subTaxMatrix <- subTaxMatrix[!duplicated(subTaxMatrix),]
        # remove IDs that do not have real NCBI superRank
        nameList <- getNameList(taxDB)
        subNameList <- nameList[
            nameList$ncbiID %in% subTaxMatrix[[superRank]],
            c("ncbiID", "rank", "fullName")
        ]
        colnames(subNameList) <- c(superRank, "rank", "name")
        mergedDf <- merge(subTaxMatrix, subNameList, by = superRank, all.x=TRUE)
        mergedDf <- mergedDf[mergedDf$rank == superRank,]
        subTaxMatrix <- mergedDf
        # group input taxa based on superrank and get their index
        groupedList <- lapply(
            levels(as.factor(subTaxMatrix[[superRank]])),
            function(x) {
                tmp <- subTaxMatrix[subTaxMatrix[[superRank]] == x,]
                tmp$index <- which(
                    levels(as.factor(plotDf$supertaxonID))
                    %in% tmp[[workingRank]]
                )
                return(tmp)
            }
        )
        # collect groups (one per superRank taxon), ordered by position
        groupDf <- do.call(rbind, lapply(groupedList, function(i) {
            data.frame(
                id = i[[superRank]][1], name = as.character(i$name[1]),
                min = min(i$index), max = max(i$index),
                stringsAsFactors = FALSE
            )
        }))
        groupDf <- groupDf[order(groupDf$min), ]
        nG <- nrow(groupDf)
        # minimum distance (in cells) between two labels
        needCells <- 0
        if (!is.null(plotSize) && isTRUE(plotSize > 0)) {
            nCells <- nlevels(as.factor(plotDf$supertaxonID))
            pxPerCell <- plotSize * 0.85 / nCells        # ~85% of plot = panel
            labelPx <- groupLabelSize * 72 / 25.4 * 1.3  # text height + gap
            needCells <- labelPx / pxPerCell
        }
        if (identical(labelMode, "first")) {
            # keep the first label of a run of overlapping labels, add ",..."
            cl <- integer(nG)
            cl[1] <- 1
            start <- groupDf$min[1]
            for (k in seq_len(nG)[-1]) {
                if (groupDf$min[k] - start >= needCells) {
                    cl[k] <- cl[k - 1] + 1
                    start <- groupDf$min[k]
                } else cl[k] <- cl[k - 1]
            }
            labelDf <- do.call(rbind, lapply(
                split(seq_len(nG), cl),
                function(idx) {
                    lab <- groupDf$name[idx[1]]
                    if (length(idx) > 1) lab <- paste0(lab, ",...")
                    data.frame(
                        pos = min(groupDf$min[idx]), label = lab,
                        from = min(groupDf$min[idx]),
                        to = max(groupDf$max[idx]), merged = FALSE,
                        stringsAsFactors = FALSE
                    )
                }
            ))
        } else {
            # ancestors of each group above superRank (columns go low -> high)
            cols <- colnames(taxMatrix)
            rIdx <- match(superRank, cols)
            higher <- if (rIdx < length(cols)) cols[(rIdx + 1):
                        length(cols)] else character(0)
            subTax <- as.data.frame(
                taxMatrix[taxMatrix[[superRank]] %in% groupDf$id, ]
            )
            rowIdx <- match(groupDf$id, subTax[[superRank]])
            ancMat <- matrix(
                paste0("self_", groupDf$id), nrow = nG, 
                ncol = max(1, length(higher))
            )
            for (j in seq_along(higher)) {
                v <- as.character(subTax[[higher[j]]][rowIdx])
                ancMat[, j] <- ifelse(
                    is.na(v) | v == "", paste0("NA_", groupDf$id), v
                )
            }
            # lift overlapping groups one level at a time; neighbours sharing
            # the same ancestor at their level are merged into one label.
            # All other groups belonging to that ancestor are lifted as well,
            # even if their own labels do not overlap
            lv <- integer(nG)
            getClusters <- function(lv) {
                key <- ifelse(
                    lv == 0, paste0("0:", groupDf$id),
                    paste0(lv, ":", ancMat[cbind(seq_len(nG), pmax(lv, 1))])
                )
                cumsum(c(TRUE, key[-1] != key[-nG]))
            }
            cl <- getClusters(lv)
            for (L in seq_along(higher)) {
                starts <- tapply(groupDf$min, cl, min)
                bad <- which(diff(starts) < needCells)
                if (length(bad) == 0) break
                hit <- cl %in% c(bad, bad + 1)
                lv[ancMat[, L] %in% unique(ancMat[hit, L])] <- L
                cl <- getClusters(lv)
            }
            # one label per cluster: own name, or ancestor name + "*"
            labelDf <- do.call(
                rbind, lapply(split(seq_len(nG), cl), function(idx) {
                    if (length(idx) == 1) {
                        lab <- groupDf$name[idx]
                    } else {
                        anc <- ancMat[idx[1], lv[idx[1]]]
                        lab <- nameList$fullName[
                            as.character(nameList$ncbiID) == anc][1]
                        if (is.na(lab)) lab <- groupDf$name[idx[1]]
                        lab <- paste0(lab, "*")
                    }
                    data.frame(
                        pos = min(groupDf$min[idx]), label = as.character(lab),
                        from = min(groupDf$min[idx]),
                        to = max(groupDf$max[idx]), merged = length(idx) > 1,
                        stringsAsFactors = FALSE
                    )
            }))
        }
        # labels still overlapping at the top level are hidden
        keep <- logical(nrow(labelDf))
        last <- -Inf
        for (k in seq_len(nrow(labelDf))) {
            if (labelDf$pos[k] - last >= needCells) {
                keep[k] <- TRUE
                last <- labelDf$pos[k]
            }
        }
        labelDf <- labelDf[keep, , drop = FALSE]
        # draw division lines (all groups) and labels
        linePos <- unique(c(groupDf$min - 0.5, groupDf$max + 0.5))
        # bars along the plot edge showing the extent of merged (LCA) groups;
        # slightly shortened so that neighbouring bars stay separated
        barDf <- labelDf[labelDf$merged, , drop = FALSE]
        barFrom <- barDf$from - 0.5 + 0.2
        barTo <- barDf$to + 0.5 - 0.2
        # leading space keeps labels clear of the bars
        labelDf$label <- paste0(" ", labelDf$label)
        if (xAxis == "taxa") {
            if (nrow(barDf) > 0)
                profilePlot <- profilePlot + annotate(
                    geom = "segment", x = barFrom, xend = barTo,
                    y = Inf, yend = Inf, colour = lineColor,
                    linewidth = 2, lineend = "butt"
                )
            profilePlot <- profilePlot +
                geom_vline(xintercept = linePos, colour = lineColor) +
                annotate(
                    geom = "text", angle = groupLabelAngle,
                    hjust = 0, vjust = 1.2, size = groupLabelSize,
                    x = labelDf$pos - 0.5, y = Inf, label = labelDf$label
                )
        } else {
            if (nrow(barDf) > 0)
                profilePlot <- profilePlot + annotate(
                    geom = "segment", y = barFrom, yend = barTo,
                    x = Inf, xend = Inf, colour = lineColor,
                    linewidth = 2, lineend = "butt"
                )
            profilePlot <- profilePlot +
                geom_hline(yintercept = linePos, colour = lineColor) +
                annotate(
                    geom = "text", angle = groupLabelAngle - 90,
                    hjust = -0.05, vjust = -0.3, size = groupLabelSize,
                    y = labelDf$pos - 0.5, x = Inf, label = labelDf$label
                )
        }
        longestLabelLen <- max(nchar(labelDf$label))
        marginSize <- longestLabelLen * groupLabelSize * groupLabelDist
        if (xAxis == "taxa") {
            return(
                profilePlot +
                    coord_cartesian(clip = "off") +
                    theme(
                        plot.margin = margin(
                            t = marginSize, r = 10, b = 10, l = 10
                        )
                    )
            )
        } else
            return(
                profilePlot +
                    coord_cartesian(clip = "off") +
                    theme(
                        plot.margin = margin(
                            t = 10, r = marginSize, b = 10, l = 10
                        ),
                        legend.position = "top"
                    )
            )
    }
}

#' Add a species tree aligned to the taxon columns (rows) of a profile plot
#' @description Draw a (cladogram) species tree outside the plot panel so that
#' each tip sits exactly at the position of its (super)taxon. The tree is
#' added as an annotation (no scale training), therefore the coordinates of
#' the data points stay unchanged and click/brush events keep working.
#' @export
#' @usage addSpeciesTreePlot(profilePlot = NULL, plotDf = NULL, tree = NULL,
#'     taxaMap = NULL, xAxis = "taxa", height = 25,
#'     lineColor = "gray30", lineWidth = 0.5, returnLayout = FALSE)
#' @param profilePlot initial (highlighted) profile plot
#' @param plotDf dataframe used for the profile plot (needs a "supertaxon"
#' column)
#' @param tree a "phylo" object. Tip labels are either the supertaxon names
#' of the plot or the abbrNames (e.g. ncbi9606) given in taxaMap
#' @param taxaMap optional dataframe with columns "abbrName" and "supertaxon"
#' to map tree tips to supertaxa (e.g. the output of sortInputTaxa). The first
#' tip of each supertaxon is used as its representative
#' @param xAxis type of x-axis (either "genes" or "taxa"). The tree is drawn
#' below the taxon labels of the x-axis for "taxa" and left of the taxon
#' labels of the y-axis for "genes", so the labels are the leaf names
#' @param height height of the tree in mm
#' @param lineColor color of the tree branches
#' @param lineWidth line width of the tree branches
#' @param returnLayout if TRUE, return a list of the plot and the tree layout
#' (branch segments in data units / mm, tips below each node), e.g. for
#' showing node information on hover. Default = FALSE
#' @return A profile plot with the species tree as a ggplot object. If the
#' tree does not fit the taxon order of the plot, the input plot is returned
#' unchanged (with a message)
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{heatmapPlotting}}, \code{\link{addRankDivisionPlot}},
#' \code{\link{sortInputTaxa}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' profilePlot <- heatmapPlotting(plotDf)
#' taxa <- levels(droplevels(as.factor(plotDf$supertaxon)))
#' tree <- ape::rtree(length(taxa), tip.label = taxa)
#' addSpeciesTreePlot(profilePlot, plotDf, tree)

addSpeciesTreePlot <- function(
        profilePlot = NULL, plotDf = NULL, tree = NULL, taxaMap = NULL,
        xAxis = "taxa", height = 25,
        lineColor = "gray30", lineWidth = 0.5, returnLayout = FALSE
) {
    if (is.null(profilePlot)) stop("Profile plot cannot be NULL!")
    out <- function(p, layout = NULL) {
        if (isTRUE(returnLayout)) return(list(plot = p, layout = layout))
        return(p)
    }
    if (is.null(plotDf) || is.null(tree)) {
        message("Species tree not drawn: no plot data or no tree")
        return(out(profilePlot))
    }
    if (is.null(height) || !isTRUE(height > 0)) height <- 25
    # taxon columns (or rows) in the order of the plot
    cols <- levels(droplevels(as.factor(plotDf$supertaxon)))
    nCol <- length(cols)
    if (nCol < 2) {
        message("Species tree not drawn: less than 2 taxa in the plot")
        return(out(profilePlot))
    }
    # map tree tips to supertaxa (one representative tip per supertaxon)
    if (!is.null(taxaMap)) {
        taxaMap <- data.frame(
            abbrName = as.character(taxaMap$abbrName),
            supertaxon = as.character(taxaMap$supertaxon),
            stringsAsFactors = FALSE
        )
        # dataMainPlot removes the 7-character index prefix of the supertaxon
        # names (e.g. "100006_"), so remove it here as well if needed
        if (!any(taxaMap$supertaxon %in% cols)) {
            stripped <- substr(
                taxaMap$supertaxon, 8, nchar(taxaMap$supertaxon)
            )
            if (any(stripped %in% cols)) taxaMap$supertaxon <- stripped
        }
        m <- taxaMap[
            as.character(taxaMap$abbrName) %in% tree$tip.label &
                as.character(taxaMap$supertaxon) %in% cols, ]
        m <- m[!duplicated(as.character(m$supertaxon)), ]
        if (nrow(m) < 2) {
            message(
                "Species tree not drawn: taxa of the plot not found in the ",
                "tree (e.g. plot: ", paste(head(cols, 3), collapse = ", "),
                "; tree tips: ", paste(head(tree$tip.label, 3), collapse=", "),
                "; taxaMap: ", paste(head(taxaMap$abbrName, 3), collapse=", "),
                " / ", paste(head(taxaMap$supertaxon, 3), collapse=", "), ")"
            )
            return(out(profilePlot))
        }
        tree <- ape::keep.tip(tree, as.character(m$abbrName))
        repAbbr <- stats::setNames(
            as.character(m$abbrName), as.character(m$supertaxon))
        tree$tip.label <- as.character(m$supertaxon)[
            match(tree$tip.label, as.character(m$abbrName))]
    }
    missingTaxa <- setdiff(cols, tree$tip.label)
    if (length(missingTaxa) > 0) {
        message(
            "Species tree not drawn: ", length(missingTaxa),
            " taxa of the plot not found in the tree"
        )
        return(out(profilePlot))
    }
    if (length(setdiff(tree$tip.label, cols)) > 0)
        tree <- ape::keep.tip(tree, cols)
    tree <- ape::collapse.singles(tree)
    # root using the first taxon (reference), the same way as sortInputTaxa
    tree <- ape::root(tree, outgroup = cols[1], resolve.root = TRUE)
    tree <- ape::reorder.phylo(tree, "postorder")
    if (!exists("repAbbr", inherits = FALSE))
        repAbbr <- stats::setNames(rep(NA_character_, nCol), cols)
    # rectangular cladogram layout: tip i at position i (= its column)
    e <- tree$edge
    nTip <- length(tree$tip.label)
    nNode <- nTip + tree$Nnode
    rootNode <- nTip + 1
    tipPos <- match(tree$tip.label, cols)
    pos <- c(tipPos, rep(NA_real_, tree$Nnode))
    depth <- rep(0, nNode)
    lo <- c(tipPos, rep(Inf, tree$Nnode))
    hi <- c(tipPos, rep(-Inf, tree$Nnode))
    nTips <- c(rep(1, nTip), rep(0, tree$Nnode))
    minChild <- rep(Inf, nNode)
    maxChild <- rep(-Inf, nNode)
    descTips <- c(as.list(tree$tip.label), rep(list(character(0)), tree$Nnode))
    for (i in seq_len(nrow(e))) {
        p <- e[i, 1]
        ch <- e[i, 2]
        # postorder: all edges below ch are already processed
        if (ch > nTip) pos[ch] <- (minChild[ch] + maxChild[ch]) / 2
        lo[p] <- min(lo[p], lo[ch])
        hi[p] <- max(hi[p], hi[ch])
        nTips[p] <- nTips[p] + nTips[ch]
        minChild[p] <- min(minChild[p], pos[ch])
        maxChild[p] <- max(maxChild[p], pos[ch])
        depth[p] <- max(depth[p], depth[ch] + 1)
        descTips[[p]] <- c(descTips[[p]], descTips[[ch]])
    }
    pos[rootNode] <- (minChild[rootNode] + maxChild[rootNode]) / 2
    # every clade must cover a contiguous block of columns, otherwise the
    # tree cannot be drawn without crossing branches
    inner <- (nTip + 1):nNode
    if (any(hi[inner] - lo[inner] + 1 != nTips[inner])) {
        message(
            "Species tree not drawn: taxon order of the plot does not ",
            "match the tree"
        )
        return(out(profilePlot))
    }
    # styles of the taxon labels (axis text): the tree is placed below the
    # x-axis labels (xAxis = "taxa") or left of the y-axis labels ("genes")
    th <- tryCatch(theme_get() + profilePlot$theme, error = function(e) NULL)
    getEl <- function(name) tryCatch(
        calc_element(name, th), error = function(e) NULL
    )
    isGenes <- xAxis == "genes"
    txt <- getEl(if (isGenes) "axis.text.y.left" else "axis.text.x.bottom")
    tickLen <- getEl(
        if (isGenes) "axis.ticks.length.y.left" else"axis.ticks.length.x.bottom"
    )
    toMm <- function(u) tryCatch(
        grid::convertUnit(u, "mm", valueOnly = TRUE), error = function(e) 1
    )
    labelStyle <- list(
        size = if (!is.null(txt$size)) txt$size else 8,
        family = if (!is.null(txt$family)) txt$family else "",
        angle = if (!is.null(txt$angle)) txt$angle else 0,
        lineheight = if (!is.null(txt$lineheight)) txt$lineheight else 0.9,
        margin = if (is.null(txt$margin)) 1 else toMm(
            if (isGenes) txt$margin[2] else txt$margin[1]),
        tick = if (is.null(tickLen)) 1 else toMm(tickLen)
    )
    gap <- 1
    # filled at drawing time: offset of the tree (mm) and panel position (px)
    drawInfo <- new.env()
    treeGrob <- grid::gTree(
        cl = "ppSpeciesTree", labels = cols, labelStyle = labelStyle,
        isGenes = isGenes, nCol = nCol, gap = gap, height = height,
        maxDepth = max(depth), lineColor = lineColor, lineWidth = lineWidth,
        seg = data.frame(
            a0 = c(pos[e[, 2]], minChild[inner]),
            a1 = c(pos[e[, 2]], maxChild[inner]),
            d0 = c(depth[e[, 2]], depth[inner]),
            d1 = c(depth[e[, 1]], depth[inner])
        ),
        drawInfo = drawInfo
    )
    if (isGenes) {
        profilePlot <- profilePlot +
            annotation_custom(
                treeGrob, xmin = -Inf, xmax = -Inf,
                ymin = 0.5, ymax = nCol + 0.5
            ) +
            theme(axis.title.y = element_text(
                margin = margin(r = gap + height + 2, unit = "mm")
            ))
    } else {
        profilePlot <- profilePlot +
            annotation_custom(
                treeGrob, xmin = 0.5, xmax = nCol + 0.5,
                ymin = -Inf, ymax = -Inf
            ) +
            theme(axis.title.x = element_text(
                margin = margin(t = gap + height + 2, unit = "mm")
            ))
    }
    if (!identical(profilePlot$coordinates$clip, "off"))
        profilePlot <- profilePlot + coord_cartesian(clip = "off")
    # layout of the drawn branches (e.g. for hover info): position along the
    # taxon axis in data units, depth scaled to 0..height (mm, to be added to
    # drawInfo$offset), and the node each segment belongs to (vertical branch
    # -> child node, horizontal connector -> parent node)
    toAway <- function(d) gap + d / max(depth) * height
    layout <- list(
        xAxis = xAxis, nCol = nCol,
        segments = data.frame(
            node = c(e[, 2], inner),
            along0 = treeGrob$seg$a0, along1 = treeGrob$seg$a1,
            away0 = toAway(treeGrob$seg$d0), away1 = toAway(treeGrob$seg$d1)
        ),
        isTip = seq_len(nNode) <= nTip,
        descTips = descTips,
        repAbbr = repAbbr,
        drawInfo = drawInfo
    )
    return(out(profilePlot, layout))
}

#' Draw the species tree grob (created by addSpeciesTreePlot)
#' @description The size of the taxon labels is only known at drawing time,
#' therefore the tree is placed here: below (or left of) the axis labels.
#' The offset and the panel position on the device are stored in
#' x$drawInfo (used e.g. for the hover info in the Shiny app).
#' @param x a ppSpeciesTree grob
#' @return the grob with its branch segments
#' @importFrom grid makeContent
#' @exportS3Method grid::makeContent
makeContent.ppSpeciesTree <- function(x) {
    st <- x$labelStyle
    # extent of the rotated label block: longest label (w) and text height (h)
    gpText <- grid::gpar(
        fontsize = st$size, fontfamily = st$family, lineheight = st$lineheight
    )
    w <- max(vapply(x$labels, function(l) grid::convertWidth(
        grid::grobWidth(grid::textGrob(l, gp = gpText)), "mm", valueOnly = TRUE
    ), numeric(1)))
    h <- grid::convertHeight(
        grid::grobHeight(grid::textGrob("Ag", gp = gpText)), "mm",
        valueOnly = TRUE
    )
    rad <- st$angle * pi / 180
    labelMm <- if (x$isGenes) {
        w * abs(cos(rad)) + h * abs(sin(rad))
    } else w * abs(sin(rad)) + h * abs(cos(rad))
    offset <- st$tick + st$margin + labelMm
    away <- function(d) grid::unit(
        -(offset + x$gap + d / x$maxDepth * x$height), "mm"
    )
    along <- function(a) grid::unit((a - 0.5) / x$nCol, "npc")
    gp <- grid::gpar(
        col = x$lineColor, lwd = x$lineWidth * .pt, lineend = "square"
    )
    sg <- x$seg
    seg <- if (x$isGenes) {
        grid::segmentsGrob(
            x0 = away(sg$d0), x1 = away(sg$d1),
            y0 = along(sg$a0), y1 = along(sg$a1), gp = gp
        )
    } else {
        grid::segmentsGrob(
            x0 = along(sg$a0), x1 = along(sg$a1),
            y0 = away(sg$d0), y1 = away(sg$d1), gp = gp
        )
    }
    # remember where the tree was drawn (device inches, origin bottom-left)
    tryCatch({
        ll <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"))
        ur <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"))
        x$drawInfo$offset <- offset
        x$drawInfo$left <- as.numeric(ll$x)
        x$drawInfo$right <- as.numeric(ur$x)
        x$drawInfo$bottom <- as.numeric(ll$y)
        x$drawInfo$top <- as.numeric(ur$y)
        x$drawInfo$devHeight <- grDevices::dev.size("in")[2]
    }, error = function(e) NULL)
    grid::setChildren(x, grid::gList(seg))
}
