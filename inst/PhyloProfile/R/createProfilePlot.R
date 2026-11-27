#' Profile plot
#' @param data data for heatmap plot (from reactive fn "dataHeat")
#' @param clusteredDataHeat clustered data (from reactive fn "clusteredDataHeat"
#' @param applyCluster choose clustered data or not (from input$applyCluster)
#' @param parameters plot parameters (colors, size, variable names, ...)
#' @param inSeq subset sequences for customized profile (input$inSeq)
#' @param inTaxa subset taxa for customized profile (input$inTaxa)
#' @param rankSelect selected taxonomy rank (input$rankSelect)
#' @param inSelect selected taxon name (input$inSelect)
#' @param taxonHighlight highlighted taxon (input$taxonHighlight)
#' @param geneHighlight highlighted gene (input$geneHighlight)
#' @param typeProfile either "mainProfile" or "customizedProfile"
#' @param taxDB path to taxonomy database (taxonomyMatrix.txt file required!)
#' @param superRank taxonomy rank for division lines (e.g. superkingdom)
#' @param allTaxa list of all input (super)taxa and their ncbi IDs
#' @param mode plot mode, either "fast" or "normal"
#' @param speciesTree NULL or list(tree = phylo, taxaMap = df(abbrName,
#' supertaxon), height = tree height in mm) for drawing a species tree aligned
#' to the taxa (main plot in normal mode, subplot in fast mode)
#' @taxDB Path to the taxonomy DB files
#' @return info for selected point on the profile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")
rtCheck <- FALSE

createProfilePlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("plot.ui")),
        br(),
        downloadButton(ns("profileDownload"),"Download profile",
                       class = "butDL"),
        br(),br(),hr(),
        radioButtons(
            ns("showSettings"), "Show subplot settings",
            choices = c("show", "hide"), selected = "hide", inline = TRUE
        ),
        conditionalPanel(
            condition = "input.showSettings == 'show'",
            column(
                2,
                radioButtons(
                    inputId = ns("xAxisSub"),
                    label = "X-axis",
                    choices = list("taxa", "genes"),
                    selected = "taxa",
                    inline = TRUE
                )
            ),
            column(
                2,
                createPlotSize(ns("widthSub"), "Width (px)", 1200)
            ),
            column(
                2,
                createPlotSize(ns("heightSub"), "Height (px)", 1200)
            ),
            column(
                2,
                createTextSize(ns("fontSizeSub"), "Axis size", 14,100)
            ),
            column(
                2,
                createTextSize(ns("legendSizeSub"), "Legend size", 8,100)
            ),
            column(
                2,
                selectInput(
                    ns("legendPosSub"), label = "Legend position:",
                    choices = list(
                        "Right" = "right", "Left" = "left", "Top" = "top",
                        "Bottom" = "bottom", "Hide" = "none"
                    ),
                    selected = "right"
                )
            ),
            ns = NS(id)
        ),
        uiOutput(ns("brushedInfo.ui")),
        br(),
        downloadButton(
            ns("profileSubDownload"),"Download subplot", class = "butDL"
        ),
        strong(""),
        downloadButton(
            ns("dataSubDownload"),"Download selected data", class = "butDL"
        ),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        )
    )
}

createProfilePlot <- function(
        input, output, session, data, clusteredDataHeat, applyCluster, parameters,
        inSeq, inTaxa, rankSelect, inSelect, taxonHighlight, geneHighlight,
        typeProfile, taxDB, superRank, allTaxa, mode = "normal",
        labelMode = reactive("lca"), speciesTree = reactive(NULL)
) {
    # data for heatmap ---------------------------------------------------------
    dataHeat <- reactive({
        if (is.null(data())) stop("Profile data is NULL!")
        if (rtCheck) {
            checkpointp101 <- Sys.time()
            print(paste("checkpoint-P101 - before data heatmap",checkpointp101))
        }
        if (typeProfile() == "customizedProfile") {
            if (is.null(inTaxa()) | is.null(inSeq())) return()
            if (applyCluster() == TRUE) {
                dataHeat <- dataCustomizedPlot(
                    clusteredDataHeat(), inTaxa(), inSeq()
                )
            } else { dataHeat <- dataCustomizedPlot(data(), inTaxa(), inSeq()) }
        } else {
            if (applyCluster() == TRUE) {
                dataHeat <- dataMainPlot(clusteredDataHeat())
            } else { dataHeat <- dataMainPlot(data()) }
        }
        # add all input taxa (if available)
        if (!(is.null(allTaxa()))) {
            mergedDf <- dataHeat %>% full_join(
                allTaxa(), by = c("supertaxonID", "supertaxon")
            ) %>% mutate(
                geneID = ifelse(is.na(geneID), dataHeat$geneID[1], geneID),
                supertaxon = factor(
                    supertaxon, levels = levels(allTaxa()$supertaxon)
                )
            )
            return(mergedDf)
        }
        if (rtCheck) {
            checkpointp102 <- Sys.time()
            print(paste(
                "checkpoint-P102 - data heatmap done", checkpointp102,
                " --- ",  checkpointp102 - checkpointp101
            ))
        }
        return(dataHeat)
    })
    
    # basic profile plot -------------------------------------------------------
    basicProfile <- reactive({
        if (is.null(dataHeat())) stop("Profile data is NULL!")
        if (typeProfile() == "customizedProfile") {
            if (length(inSeq()) == 0 || length(inTaxa()) == 0) return()
            if ("all" %in% inSeq() & "all" %in% inTaxa()) return()
        }
        if (rtCheck) {
            checkpointp201 <- Sys.time()
            print(paste("checkpoint-P201 - before basic plot", checkpointp201))
        }
        if (mode() == "fast" & typeProfile() != "customizedProfile") {
            hpf <- heatmapPlottingFast(dataHeat(), parameters())
            if (rtCheck) {
                checkpointp202 <- Sys.time()
                print(paste(
                    "checkpoint-P202 - basic plot done", checkpointp202,
                    " --- ",  checkpointp202 - checkpointp201
                ))
            }
            return(hpf)
        } else {
            hp <- heatmapPlotting(dataHeat(), parameters())
            if (rtCheck) {
                checkpointp202 <- Sys.time()
                print(paste(
                    "checkpoint-P202 - basic plot done", checkpointp202,
                    " --- ",  checkpointp202 - checkpointp201
                ))
            }
            return(hp)
        }
    })
    
    # get superRank ------------------------------------------------------------
    getSuperRank <- reactive({
        if (!is.null(superRank()) && superRank() == "") return(NULL)
        else return(superRank())
    })
    
    # add species tree aligned to the taxa ----------------------------------
    # (main plot only in normal mode; subplot of fast mode with sub = TRUE)
    # store: reactiveVal that keeps the tree layout for the hover info
    treeLayoutMain <- reactiveVal(NULL)
    treeLayoutSub <- reactiveVal(NULL)
    addTreeToPlot <- function(
        p, df, xAxis = parameters()$xAxis, sub = FALSE, store = NULL
    ) {
        if (!is.null(store)) store(NULL)
        if (!sub && mode() == "fast" && typeProfile() == "mainProfile")
            return(p)
        st <- speciesTree()
        if (is.null(st) || is.null(st$tree)) return(p)
        res <- addSpeciesTreePlot(
            p, df, st$tree, st$taxaMap, xAxis, st$height,
            returnLayout = TRUE
        )
        if (!is.null(store)) store(res$layout)
        return(res$plot)
    }

    # hover info for the internal nodes of the species tree -------------------
    hoverTaxMatrix <- reactive({ getTaxonomyMatrix(taxDB()) })
    hoverNameList <- reactive({ getNameList(taxDB()) })

    # find the tree node (branch) next to the mouse pointer. mouse: position
    # (css px) on the plot image, sent by the JS in treeHoverContainer().
    # Shiny's own hover is not used, because it is not sent outside the panel
    getTreeNodeAtHover <- function(mouse, layout, tolerance = 1.5) {
        if (is.null(mouse) || is.null(layout)) return(NULL)
        di <- layout$drawInfo
        if (is.null(di$offset) || is.null(di$left)) return(NULL)
        # renderPlot draws with res = 72, i.e. 1 inch on the device = 72 css px
        pxPerMm <- 72 / 25.4
        left <- di$left * 72
        right <- di$right * 72
        bottom <- (di$devHeight - di$bottom) * 72
        top <- (di$devHeight - di$top) * 72
        if (layout$xAxis == "genes") {
            along <- 0.5 + (bottom - mouse$y) / (bottom - top) * layout$nCol
            away <- (left - mouse$x) / pxPerMm
            mmPerUnit <- (bottom - top) / layout$nCol / pxPerMm
        } else {
            along <- 0.5 + (mouse$x - left) / (right - left) * layout$nCol
            away <- (mouse$y - bottom) / pxPerMm
            mmPerUnit <- (right - left) / layout$nCol / pxPerMm
        }
        away <- away - di$offset
        if (!is.finite(away) || away < 0) return(NULL)
        # distance (mm) from the pointer to each branch segment
        sg <- layout$segments
        px <- along * mmPerUnit
        x0 <- sg$along0 * mmPerUnit; x1 <- sg$along1 * mmPerUnit
        y0 <- sg$away0; y1 <- sg$away1
        dx <- x1 - x0; dy <- y1 - y0
        len2 <- dx^2 + dy^2
        t <- ifelse(len2 > 0, ((px - x0) * dx + (away - y0) * dy) / len2, 0)
        t <- pmin(pmax(t, 0), 1)
        dist <- sqrt((x0 + t * dx - px)^2 + (y0 + t * dy - away)^2)
        if (min(dist) > tolerance) return(NULL)
        return(sg$node[which.min(dist)])
    }

    # container for a plot with species tree: sends the mouse position on the
    # plot image to input[[inputId]] (throttled) and holds the tooltip
    treeHoverContainer <- function(plotUI, inputId, tooltipId) {
        div(
            class = "pp-tree-hover", `data-input` = inputId,
            style = "position: relative;",
            singleton(tags$head(tags$script(HTML(
                "if (!window.ppTreeHoverInit) {
                    window.ppTreeHoverInit = true;
                    var ppLast = 0;
                    $(document).on('mousemove', '.pp-tree-hover img',
                        function(e) {
                            var now = Date.now();
                            if (now - ppLast < 80) return;
                            ppLast = now;
                            var r = this.getBoundingClientRect();
                            var id = $(this).closest('.pp-tree-hover')
                                .data('input');
                            Shiny.setInputValue(id, {
                                x: e.clientX - r.left, y: e.clientY - r.top
                            }, {priority: 'event'});
                        });
                    $(document).on('mouseleave', '.pp-tree-hover img',
                        function() {
                            var id = $(this).closest('.pp-tree-hover')
                                .data('input');
                            Shiny.setInputValue(id, null, {priority: 'event'});
                        });
                }"
            )))),
            plotUI,
            uiOutput(tooltipId)
        )
    }

    # name of a tree node: lowest taxon shared by all its tips
    getTreeNodeLabel <- function(node, layout) {
        tips <- layout$descTips[[node]]
        if (layout$isTip[node]) return(tips)
        label <- "Unknown"
        abbr <- layout$repAbbr[tips]
        abbr <- abbr[!is.na(abbr)]
        if (length(abbr) > 0) {
            tm <- hoverTaxMatrix()
            rows <- tm[as.character(tm$abbrName) %in% abbr, , drop = FALSE]
            rankCols <- setdiff(
                colnames(tm), c("abbrName", "ncbiID", "fullName")
            )
            for (rc in rankCols) {
                v <- unique(as.character(rows[[rc]]))
                if (length(v) == 1 && !is.na(v) && v != "") {
                    nl <- hoverNameList()
                    hit <- nl[as.character(nl$ncbiID) == v, ]
                    if (nrow(hit) > 0) label <- paste0(
                        hit$fullName[1], " (", hit$rank[1], ")"
                    )
                    break
                }
            }
        }
        return(paste0(label, " - ", length(tips), " taxa"))
    }

    makeTreeTooltip <- function(hover, layout) {
        node <- tryCatch(
            getTreeNodeAtHover(hover, layout), error = function(e) NULL
        )
        if (is.null(node)) return(NULL)
        label <- tryCatch(
            getTreeNodeLabel(node, layout), error = function(e) NULL
        )
        if (is.null(label)) return(NULL)
        div(
            style = paste0(
                "position: absolute; z-index: 100; pointer-events: none; ",
                "left: ", hover$x + 12, "px; ",
                "top: ", hover$y + 12, "px; ",
                "background: rgba(255,255,255,0.95); border: 1px solid #999; ",
                "border-radius: 3px; padding: 2px 6px; font-size: 12px; ",
                "white-space: nowrap;"
            ),
            label
        )
    }

    output$treeHoverInfo <- renderUI({
        makeTreeTooltip(input$treeMouse, treeLayoutMain())
    })
    output$subTreeHoverInfo <- renderUI({
        makeTreeTooltip(input$subTreeMouse, treeLayoutSub())
    })

    # render heatmap profile ---------------------------------------------------
    output$plot <- renderPlot({
        if (is.null(basicProfile())) return()
        if (rtCheck) {
            checkpointp301 <- Sys.time()
            print(paste("checkpoint-P301 - before final plot", checkpointp301))
        }
        withProgress(message = 'PLOTTING...', value = 0.5, {
            p <- highlightProfilePlot(
                basicProfile(), dataHeat(), taxonHighlight(), rankSelect(),
                geneHighlight(), taxDB(), parameters()$xAxis
            )
            refLine <- TRUE
            if (mode() == "fast") refLine <- FALSE
            p <- addRankDivisionPlot(
                p, dataHeat(), taxDB(), rankSelect(), getSuperRank(),
                parameters()$xAxis, parameters()$font,
                parameters()$groupLabelSize, parameters()$groupLabelDist,
                parameters()$groupLabelAngle, refLine,
                parameters()$divisionLineColor,
                if (parameters()$xAxis == "taxa") parameters()$width
                else parameters()$height,
                labelMode()
            )
            p <- addTreeToPlot(p, dataHeat(), store = treeLayoutMain)
            if (rtCheck) {
                checkpointp302 <- Sys.time()
                print(paste(
                    "checkpoint-P302 - final plot done", checkpointp302,
                    " --- ",  checkpointp302 - checkpointp301
                ))
            }
            return(p)
        })
    })
    
    output$plot.ui <- renderUI({
        ns <- session$ns
        if (typeProfile() == "customizedProfile") {
            if (is.null(inSeq()[1]) | is.null(inTaxa()[1])) return()
            else if (inSeq()[1] == "all" & inTaxa()[1] == "all")  return()
        }
        
        if (mode() == "fast" & typeProfile() == "mainProfile") {
            shinycssloaders::withSpinner(
                plotOutput(
                    ns("plot"),
                    width = parameters()$width,
                    height = parameters()$height,
                    brush = ns("plotBrush")
                )
            )
        } else {
            # container with hover tooltip for the species tree
            treeHoverContainer(
                shinycssloaders::withSpinner(
                    plotOutput(
                        ns("plot"),
                        width = parameters()$width,
                        height = parameters()$height,
                        click = ns("plotClick")
                    )
                ),
                ns("treeMouse"), ns("treeHoverInfo")
            )
        }
    })
    
    output$profileDownload <- downloadHandler(
        filename = function() {
            c("profile.svg")
        },
        content = function(file) {
            ggsave(
                file,
                plot = addTreeToPlot(addRankDivisionPlot(
                    basicProfile(), dataHeat(), taxDB(), rankSelect(),
                    getSuperRank(), parameters()$xAxis, parameters()$font,
                    parameters()$groupLabelSize, parameters()$groupLabelDist,
                    parameters()$groupLabelAngle, FALSE,
                    parameters()$divisionLineColor,
                    if (parameters()$xAxis == "taxa") parameters()$width
                    else parameters()$height,
                    labelMode()
                ), dataHeat()),
                width = parameters()$width * 0.035,
                height = parameters()$height * 0.035,
                units = "cm", dpi = 300, device = "svg", limitsize = FALSE
            )
        }
    )
    
    # get info of brushed area -------------------------------------------------
    brushedData <- reactive({
        if (is.null(input$plotBrush)) return(data.frame())
        if (parameters()$xAxis == "genes") {
            dt <- brushedPoints(
                dataHeat(), input$plotBrush, xvar = "geneID", yvar = "supertaxon"
            )
        } else {
            dt <- brushedPoints(
                dataHeat(), input$plotBrush, xvar = "supertaxon", yvar = "geneID"
            )
        }
        dt$geneID <- factor(dt$geneID)
        dt$geneID <- droplevels(dt$geneID)
        dt$supertaxon <- factor(dt$supertaxon)
        dt$supertaxon <- droplevels(dt$supertaxon)
        return(dt)
    })
    
    # create subplot -----------------------------------------------------------
    # ** deactivate download buttons -------------------------------------------
    observe({
        shinyjs::toggleState(
            id = "profileSubDownload",
            condition = nrow(brushedData()) > 0
        )
    })
    observe({
        shinyjs::toggleState(
            id = "dataSubDownload",
            condition = nrow(brushedData()) > 0
        )
    })
    
    # ** show subplot settings whenever subplot data is available -------------
    # switch "Show subplot settings" to "show" when a selection with data is
    # made (and back to "hide" when the selection is cleared). Only reacts to
    # changes in availability, so the user can still hide the settings manually
    # for the current selection.
    subDataAvailable <- reactiveVal(FALSE)
    observe({
        available <- tryCatch(
            nrow(brushedData()) > 0, error = function(e) FALSE
        )
        if (!identical(available, isolate(subDataAvailable()))) {
            subDataAvailable(available)
            updateRadioButtons(
                session, "showSettings",
                selected = if (available) "show" else "hide"
            )
        }
    })

    # ** sync subplot x-axis with main profile x-axis --------------------------
    # mainXAxis only invalidates when the main x-axis value really changes, so
    # changing other main parameters (font, colors, ...) does not reset a
    # x-axis that the user has switched manually for the subplot
    mainXAxis <- reactiveVal(NULL)
    observe({
        newXAxis <- parameters()$xAxis
        if (!is.null(newXAxis) && !identical(newXAxis, isolate(mainXAxis())))
            mainXAxis(newXAxis)
    })
    observeEvent(mainXAxis(), {
        if (!identical(mainXAxis(), input$xAxisSub)) {
            updateRadioButtons(session, "xAxisSub", selected = mainXAxis())
        }
    })
    
    # ** max length of supertaxon (division) labels of the subplot -----------
    subSuperRankLabelLen <- reactive({
        req(nrow(brushedData()) > 0)
        if (is.null(getSuperRank())) return(0)
        getSuperRankLabelLen(
            unique(as.character(brushedData()$supertaxonID)), rankSelect(),
            getSuperRank(), taxDB()
        )
    })

    # ** auto update subplot size ----------------------------------------------
    # size is based on the CURRENT subplot x-axis (input$xAxisSub), so when the
    # user switches xAxisSub, widthSub and heightSub are switched as well.
    # The size includes the space for axis labels, supertaxon labels and the
    # legend (see estimatePlotSize()).
    # This observer only writes widthSub/heightSub/legendPosSub and never reads
    # them reactively, so there is no update loop with the x-axis sync above.
    observe({
        req(nrow(brushedData()) > 0)
        req(input$xAxisSub)
        df <- brushedData()
        nrTaxa <- length(unique(as.character(df$supertaxon)))
        nrGene <- length(unique(as.character(df$geneID)))
        geneLabels <- df$geneID
        if (identical(parameters()$geneIdType, "geneName") &&
            "geneName" %in% colnames(df)) geneLabels <- df$geneName
        groupLabelLen <- subSuperRankLabelLen()
        nrY <- if (input$xAxisSub == "taxa") nrGene else nrTaxa
        # legend: on top for small profiles without division lines
        legendPos <- isolate(input$legendPosSub)
        if (is.null(legendPos) || legendPos == "") legendPos <- "right"
        if (nrY <= 20) {
            legendPos <- if (groupLabelLen == 0) "top" else "right"
            updateSelectInput(
                session, "legendPosSub", label = "Legend position:",
                choices = list("Right" = "right",
                               "Left" = "left",
                               "Top" = "top",
                               "Bottom" = "bottom",
                               "Hide" = "none"),
                selected = legendPos
            )
        }
        fontSize <- if (is.null(input$fontSizeSub)) parameters()$xSize
            else input$fontSizeSub
        legendSize <- if (is.null(input$legendSizeSub)) parameters()$legendSize
            else input$legendSizeSub
        adaptedSize <- estimatePlotSize(
            nrTaxa, nrGene, input$xAxisSub, dotZoom = parameters()$dotZoom,
            taxaLabelLen = max(nchar(as.character(df$supertaxon))),
            geneLabelLen = max(nchar(as.character(geneLabels)), na.rm = TRUE),
            xSize = fontSize, ySize = fontSize,
            xAngle = parameters()$xAngle,
            groupLabelLen = groupLabelLen,
            groupLabelSize = parameters()$groupLabelSize,
            groupLabelDist = parameters()$groupLabelDist,
            legendPos = legendPos, legendSize = legendSize
        )
        updateNumericInput(
            session, "widthSub", value = min(25000, adaptedSize[3]))
        updateNumericInput(
            session, "heightSub", value = min(25000, adaptedSize[2]))
    })
    
    # ** get settings for subplot; replace xAxis, font size and legend pos -----
    subPlotParameters <- reactive({
        plotParameters <- parameters()
        plotParameters$xAxis <- input$xAxisSub
        plotParameters$mainLegend <- input$legendPosSub
        plotParameters$xSize <- input$fontSizeSub
        plotParameters$ySize <- input$fontSizeSub
        plotParameters$legendSize <- input$legendSizeSub
        # vertical taxon labels when the species tree is shown below them
        treeShown <- tryCatch(
            !is.null(speciesTree()), error = function(e) FALSE
        )
        if (identical(input$xAxisSub, "taxa") && isTRUE(treeShown))
            plotParameters$xAngle <- 90
        return(plotParameters)
    })
    
    # ** render subplot --------------------------------------------------------
    output$brushedInfo.ui <- renderUI({
        ns <- session$ns
        req(nrow(brushedData()) > 0)
        tagList(
            h4("SELECTED DATA"),
            br(),
            treeHoverContainer(
                plotOutput(
                    ns("subPlot"),
                    width = input$widthSub,
                    height = input$heightSub,
                    click = ns("plotClick")
                ),
                ns("subTreeMouse"), ns("subTreeHoverInfo")
            )
        )
    })
    
    output$subPlot <- renderPlot({
        req(nrow(brushedData()) > 0)
        p <- heatmapPlotting(brushedData(), subPlotParameters())
        p <- highlightProfilePlot(
            p, brushedData(), taxonHighlight(), rankSelect(),
            geneHighlight(), taxDB(), subPlotParameters()$xAxis
        )
        p <- addRankDivisionPlot(
            p, brushedData(), taxDB(), rankSelect(),
            getSuperRank(), subPlotParameters()$xAxis, parameters()$font,
            parameters()$groupLabelSize, parameters()$groupLabelDist,
            parameters()$groupLabelAngle, FALSE,
            parameters()$divisionLineColor,
            if (subPlotParameters()$xAxis == "taxa") input$widthSub
            else input$heightSub,
            labelMode()
        )
        p <- addTreeToPlot(
            p, brushedData(), subPlotParameters()$xAxis, sub = TRUE,
            store = treeLayoutSub
        )
        p
    })
    
    # ** download subplot ------------------------------------------------------
    output$profileSubDownload <- downloadHandler(
        filename = function() {
            c("subProfile.svg")
        },
        content = function(file) {
            p <- heatmapPlotting(brushedData(), subPlotParameters())
            p <- addRankDivisionPlot(
                p, brushedData(), taxDB(), rankSelect(),
                getSuperRank(), subPlotParameters()$xAxis, parameters()$font,
                parameters()$groupLabelSize, parameters()$groupLabelDist,
                parameters()$groupLabelAngle, FALSE,
                parameters()$divisionLineColor,
                if (subPlotParameters()$xAxis == "taxa") input$widthSub
                else input$heightSub,
                labelMode()
            )
            p <- addTreeToPlot(
                p, brushedData(), subPlotParameters()$xAxis, sub = TRUE
            )
            ggsave(
                file, plot = p,
                width = input$widthSub * 0.035, height = input$heightSub *0.035,
                units = "cm", dpi = 300, device = "svg", limitsize = FALSE
            )
        }
    )
    
    # ** download data of subplot ----------------------------------------------
    output$dataSubDownload <- downloadHandler(
        filename = function() {
            c("subData.txt")
        },
        content = function(file){
            dataOut <- brushedData() %>% select(
                geneID, supertaxonID, supertaxon, orthoID, var1, var2, geneName
            )
            if (!(identical(dataOut$geneID, dataOut$geneName)))
                dataOut <- dataOut %>% select(-c(geneName))
            contains_id <- mapply(function(id, name) {
                grepl(paste0("@", id, "@"), name)
            }, dataOut$supertaxonID, dataOut$orthoID)
            if (all(contains_id == TRUE))
                dataOut <- dataOut %>% mutate(
                    supertaxonID = paste0("ncbi", supertaxonID)
                ) %>% select(
                    geneID, ncbiID = supertaxonID, orthoID, var1, var2
                )
            write.table(
                dataOut, file, sep = "\t", row.names = FALSE, quote = FALSE
            )
        }
    )
    
    # get info of clicked point on heatmap plot --------------------------------
    selectedpointInfo <- reactive({
        req(dataHeat())
        if (typeProfile() == "customizedProfile") {
            req(inSeq())
            req(inTaxa())
        }
        # get selected supertaxon name
        taxaList <- getNameList(taxDB())
        rankName <- rankSelect()
        inSelect <- taxaList$ncbiID[taxaList$fullName == inSelect()]
        
        dataHeat <- dataHeat()
        if (mode() == "fast" & typeProfile() == "mainProfile")
            dataHeat <- brushedData()
        
        if (is.null(dataHeat) | nrow(dataHeat) == 0) {
            message("WARNING: Data for heatmap is NULL!")
            return()
        }
        
        if (typeProfile() == "customizedProfile") {
            # get sub-dataframe of selected taxa and sequences
            if (is.null(inSeq()[1]) | is.null(inTaxa()[1])) {
                message("WARNING: Subset taxa or genes is NULL!")
                return()
            }
            if (inTaxa()[1] == "all" & inSeq()[1] != "all") {
                # select data from dataHeat for selected sequences only
                dataHeat <- subset(dataHeat, geneID %in% inSeq())
            } else if (inSeq()[1] == "all" & inTaxa()[1] != "all") {
                # select data from dataHeat for selected taxa only
                dataHeat <- subset(dataHeat, supertaxon %in% inTaxa())
            } else {
                # select data from dataHeat for selected sequences and taxa
                dataHeat <- subset(dataHeat, geneID %in% inSeq()
                                   & supertaxon %in% inTaxa())
            }
            
            # drop all other supertaxon that are not in sub-dataframe
            dataHeat$supertaxon <- factor(dataHeat$supertaxon)
            dataHeat$geneID <- factor(dataHeat$geneID)
        }
        
        # get values
        if (is.null(input$plotClick$x)) return()
        else {
            # get cooridiate point
            # in fast mode, clicks come from the subplot, which can have its
            # own x-axis (input$xAxisSub)
            clickXAxis <- parameters()$xAxis
            if (mode() == "fast" & typeProfile() == "mainProfile")
                clickXAxis <- input$xAxisSub
            if (clickXAxis == "genes") {
                corX <- round(input$plotClick$y);
                corY <- round(input$plotClick$x)
            } else {
                corX <- round(input$plotClick$x);
                corY <- round(input$plotClick$y)
            }
            # get geneID
            genes <- levels(dataHeat$geneID)
            geneID <- toString(genes[corY])
            # get supertaxon (spec)
            supertaxa <- levels(dataHeat$supertaxon)
            spec <- toString(supertaxa[corX])
            # get var1 and var2 score
            var1 <- NA
            if (!is.na(dataHeat$var1[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var1 <- max(
                    na.omit(dataHeat$var1[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec])
                )
            }
            var2 <- NA
            if (!is.na(dataHeat$var2[dataHeat$geneID == geneID
                                     & dataHeat$supertaxon == spec][1])) {
                var2 <- {
                    max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]))
                }
            }
            # get percentage of present species and total of taxa
            Percent <- NA
            if (!is.na(dataHeat$presSpec[dataHeat$geneID == geneID
                                         & dataHeat$supertaxon == spec][1])) {
                Percent <- {
                    max(
                        na.omit(
                            dataHeat$presSpec[dataHeat$geneID == geneID
                                              & dataHeat$supertaxon == spec]
                        )
                    )
                }
            }
            presentTaxa <- NA
            if (!is.na(dataHeat$presentTaxa[dataHeat$geneID == geneID
                                            & dataHeat$supertaxon == spec][1])){
                presentTaxa <- max(
                    na.omit(dataHeat$presentTaxa[dataHeat$geneID == geneID
                                                 & dataHeat$supertaxon == spec])
                )
            }
            totalTaxa <- NA
            if (
                !is.na(dataHeat$totalTaxa[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec][1])
            ) {
                totalTaxa <- max(
                    na.omit(
                        dataHeat$totalTaxa[dataHeat$geneID == geneID
                                           & dataHeat$supertaxon == spec]
                    )
                )
            }
            
            
            # get ortholog ID
            orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]
            totalOrtho <- length(orthoID)
            
            # get working taxonomy level
            strain <- "Y"
            if (nlevels(as.factor(dataHeat$totalTaxa)) > 1) strain <- "N"
            
            if (is.na(Percent)) return()
            else {
                info <- c(
                    geneID,
                    list(orthoID),
                    totalOrtho,
                    spec,
                    round(var1, 2),
                    round(Percent, 2),
                    round(var2, 2),
                    presentTaxa,
                    totalTaxa,
                    strain
                )
                return(info)
            }
        }
    })
    
    
    return(selectedpointInfo)
}
