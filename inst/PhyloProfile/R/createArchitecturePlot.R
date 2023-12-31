#' Protein domain architecture plot
#' @param pointInfo() info of clicked point
#' (from reactive fn "pointInfoDetail")
#' @param domainInfo() domain information
#' (from reactive fn "getDomainInformation")
#' @param labelArchiSize lable size (from input$labelArchiSize)
#' @param titleArchiSize title size (from input$titleArchiSize)
#' @param archiHeight plot height (from input$archiHeight)
#' @param archiWidth plot width (from input$archiWidth)
#' @param seqIdFormat sequence ID format (either bionf or unknown)
#' @param currentNCBIinfo dataframe of the pre-processed NCBI taxonomy data
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

createArchitecturePlotUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            4,
            style = "padding:0px;",
            selectInput(
                ns("showFeature"),
                "Show",
                choices = c(
                    "All features" = "all",
                    "Shared features" = "common",
                    "Unique features" = "unique"
                ),
                selected = "all"
            )
        ),
        column(
            4,
            selectizeInput(
                ns("excludeFeature"),
                "Exclude feature type(s)",
                choices = c(
                    "flps","seg","coils","signalp","tmhmm","smart","pfam"
                ),
                multiple = TRUE, options=list(placeholder = 'None')
            )
        ),
        column(4, uiOutput(ns("featureList.ui"))),
        column(12, uiOutput(ns("archiPlot.ui"))),
        downloadButton(ns("archiDownload"), "Download plot", class = "butDL"),
        tags$head(
            tags$style(HTML(
                ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
        ),
        br(),
        br(),
        h4(strong("LINKS TO ONLINE DATABASE")),
        textOutput(ns("selectedDomain")),
        tableOutput(ns("domainTable")),
        HTML(
            paste0(
                "<p><em><strong>Disclaimer:</strong> ",
                "External links are automatically generated and may point to ",
                "a wrong target (see <a ",
                "href=\"https://github.com/BIONF/PhyloProfile/wiki/FAQ",
                "#wrong-info-from-public-databases\" ",
                "target=\"_blank\">FAQ</a>)</em></p>"
            )
        )
    )
}

createArchitecturePlot <- function(
    input, output, session, pointInfo, domainInfo, labelArchiSize, 
    titleArchiSize, archiHeight, archiWidth, seqIdFormat, currentNCBIinfo
){
    # * filter domain features -------------------------------------------------
    filterDomainDf <- reactive({
        df <- domainInfo()
        if (is.null(nrow(df))) stop("Domain info is NULL!")
        # filter domain df by feature type
        if (!("feature_type" %in% colnames(df))) {
            df[c("feature_type","feature_id")] <- 
                stringr::str_split_fixed(df$feature, '_', 2)
            df$feature_id[df$feature_type == "smart"] <-
                paste0(df$feature_id[df$feature_type == "smart"], "_smart")
        }
        df <- df[!(df$feature_type %in% input$excludeFeature),]
        return(df)
    })
    
    getSeqIdFormat <- reactive({
        if (seqIdFormat() == 1) return("bionf")
        return("unknown")
    })
    
    # * render plot ------------------------------------------------------------
    output$archiPlot <- renderPlot({
        if (is.null(nrow(filterDomainDf()))) stop("Domain info is NULL!")
        # remove user specified features (from input$featureList)
        df <- filterDomainDf()
        df <- df[!(df$feature_id %in% input$featureList),]
        # generate plot
        g <- createArchiPlot(
            pointInfo(), df, labelArchiSize(), titleArchiSize(),
            input$showFeature, getSeqIdFormat(), currentNCBIinfo()
        )
        if (any(g == "No domain info available!")) {
            msgPlot()
        } else {
            grid.draw(g)
        }
    })

    output$archiPlot.ui <- renderUI({
        ns <- session$ns
        if (is.null(nrow(filterDomainDf()))) {
            msg <- paste0(
                "<p><em>Wrong domain file has been uploaded!
        Please check the correct format in
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/",
                "Input-Data#ortholog-annotations-eg-domains\"
        target=\"_blank\" rel=\"noopener\">our PhyloProfile wiki</a>.</em></p>"
            )
            HTML(msg)
        } else {
            shinycssloaders::withSpinner(
                plotOutput(
                    ns("archiPlot"),
                    height = archiHeight(),
                    width = archiWidth(),
                    click = ns("archiClick")
                )
            )
        }
    })

    output$archiDownload <- downloadHandler(
        filename = function() {
            c("domains.pdf")
        },
        content = function(file) {
            # remove user specified features (from input$featureList)
            df <- filterDomainDf()
            df <- df[!(df$feature_id %in% input$featureList),]
            # generate plot
            g <- createArchiPlot(
                pointInfo(), filterDomainDf(),labelArchiSize(),titleArchiSize(),
                input$showFeature, getSeqIdFormat(), currentNCBIinfo()
            )
            grid.draw(g)
            # save plot to file
            ggsave(
                file, plot = g,
                width = archiWidth() * 0.056458333,
                height = archiHeight() * 0.056458333,
                units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
            )
        }
    )

    # output$selectedDomain <- renderText({
    #     if (is.null(input$archiClick$y)) return("No domain selected!")
    #     y <- input$archiClick$y
    #     # paste(y, round(y), convertY(unit(y, "npc"), "px"))
    #     
    # })
    
    output$featureList.ui <- renderUI({
        ns <- session$ns
        allFeats <- getAllFeatures(pointInfo(), filterDomainDf())
        selectizeInput(
            ns("featureList"), "Exclude individual feature(s)", multiple = TRUE,
            choices = allFeats,
            options=list(placeholder = 'None')
        )
    })
    
    output$domainTable <- renderTable({
        if (is.null(nrow(filterDomainDf()))) return("No domain info available!")
        features <- getDomainLink(pointInfo(), filterDomainDf(), getSeqIdFormat())
        features
    }, sanitize.text.function = function(x) x)
}

#' plot error message
#' @return error message in a ggplot object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
msgPlot <- function() {
    msg <- paste(
        "No information about domain architecture!",
        "Please check:","if you uploaded the correct domain file/folder; or ",
        "if the selected genes (seed & ortholog) do exist in the uploaded file",
        "(please search for the corresponding seedID and hitID)",
        sep = "\n"
    )
    x <- c(1,2,3,4,5)
    y <- c(1,2,3,4,5)
    g <- ggplot(data.frame(x, y), aes(x,y)) +
        geom_point(color = "white") +
        annotate(
            "text", label = msg, x = 3.5, y = 0.5, size = 5, colour = "red"
        ) +
        theme(axis.line = element_blank(), axis.text = element_blank(),
              axis.ticks = element_blank(), axis.title = element_blank(),
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid = element_blank(),
              plot.background = element_blank()) +
        ylim(0,1)
    return(g)
}

#' Get list of all features
#' @return A list of all features
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getAllFeatures <- function(info, domainDf) {
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    subdomainDf <- domainDf[grep(grepID, domainDf$seedID), ]
    orthoID <- feature <- NULL
    if (nrow(subdomainDf) < 1) return(paste0("No domain info available!"))
    else {
        # ortho & seed domains df
        orthoDf <- subdomainDf[subdomainDf$orthoID == ortho,]
        seedDf <- subdomainDf[subdomainDf$orthoID != ortho,]
        feature <- c(
            levels(as.factor(orthoDf$feature_id)), 
            levels(as.factor(seedDf$feature_id))
        )
    }
    return(levels(as.factor(feature)))
}

#' get pfam and smart domain info (domain name, acc, profile HMM,...)
#' @return dataframe for each type of database
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getDomainInfo <- function(info, domainDf, type) {
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    # get sub dataframe based on selected groupID and orthoID
    group <- gsub("\\|", ":", group)
    ortho <- gsub("\\|", ":", ortho)
    grepID <- paste(group, "#", ortho, sep = "")
    domainDf <- domainDf[grep(grepID, domainDf$seedID),]
    domainDf <- domainDf[domainDf$feature_type == type,]
    domainInfoDf <- data.frame()
    if (nrow(domainDf) > 0) {
        # get all features of for this pair proteins
        feature <- getAllFeatures(info, domainDf)
        # filter domain info
        if ("acc" %in% colnames(domainDf)) {
            domainInfoDf <- domainDf[
                , c("orthoID", "feature_id", "acc", "evalue", "bitscore")
            ]
            domainInfoDf$evalue <- format(domainInfoDf$evalue, scientific=TRUE)
        } else {
            domainInfoDf <- domainDf[, c("orthoID", "feature_id")]
            domainInfoDf <- domainInfoDf[!(duplicated(domainInfoDf)),]
            domainInfoDf$acc <- rep(NA, nrow(domainInfoDf))
            domainInfoDf$evalue <- rep(NA, nrow(domainInfoDf))
            domainInfoDf$bitscore <- rep(NA, nrow(domainInfoDf))
        }
    }
    return(domainInfoDf)
}

#' get pfam and smart domain links
#' @return dataframe with domain IDs and their database links
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getDomainLink <- function(info, domainDf, seqIdFormat) {
    featurePfam <- getDomainInfo(info, domainDf, "pfam")
    pfamDf <- createLinkTable(featurePfam, "pfam")
    featureSmart <- getDomainInfo(info, domainDf, "smart")
    smartDf <- createLinkTable(featureSmart, "smart")
    
    featDf <- rbind(pfamDf, smartDf)
    featDf <- subset(featDf, select=c(orthoID,feature_id,link,evalue,bitscore))
    featDf <-featDf[order(featDf$orthoID),]
    if (seqIdFormat == "bionf") {
        featDf[c("groupID", "spec", "geneID", "misc")] <- 
            stringr::str_split_fixed(featDf$orthoID, ":", 4)
        featDf <- subset(featDf,select=c(geneID,feature_id,link,evalue,bitscore))
    }
    colnames(featDf) <- c("Gene ID", "Domain ID", "URL", "E-value", "Bit-score")
    return(featDf)
}

#' plot error message
#' @param featureDf dataframe contains 2 columns feature_id and acc
#' @param featureType pfam or smart
#' @return dataframe contains 2 columns feature_id and link
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
createLinkTable <- function(featureDf, featureType) {
    featureDf <- featureDf[!(duplicated(featureDf)),]
    if (nrow(featureDf) > 0) {
        if (featureType == "pfam") {
            featureDf$link[is.na(featureDf$acc)] <- paste0(
                "<a href='http://pfam-legacy.xfam.org/family/", featureDf$feature_id,
                "' target='_blank'>", "PFAM", "</a>"
            )
            featureDf$link[!(is.na(featureDf$acc))] <- paste0(
                "<a href='https://www.ebi.ac.uk/interpro/entry/pfam/", featureDf$acc, 
                "' target='_blank'>", "INTERPRO", "</a>"
            )
        } else {
            featureDf$link <- paste0(
                "<a href='http://smart.embl-heidelberg.de/smart/",
                "do_annotation.pl?BLAST=DUMMY&DOMAIN=",
                # featureDf$feature_id, "' target='_blank'>",
                gsub("_smart", "",featureDf$feature_id), "' target='_blank'>",
                "SMART", "</a>"
            )
        }
    }
    return(featureDf)
}
 

