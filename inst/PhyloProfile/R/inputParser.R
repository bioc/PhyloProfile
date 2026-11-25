#' Input parser module: convert other ortholog formats into PhyloProfile format
#'
#' @description UI of the input parser module. The format-specific inputs
#' (OrthoFinder, SonicParanoid, OMA) are rendered by the server part depending
#' on the selected input type; taxonomy ID source (PhyloProfile taxonomy DB
#' or own mapping file), options, preview and download are shared by all
#' formats.
#' @param id ID of the module
#' @return UI elements of the input parser
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
inputParserUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("title.ui")),
        parserSection(
            "1. Input data",
            uiOutput(ns("formatInput.ui"))
        ),
        parserSection(
            "2. NCBI taxonomy IDs of the taxa",
            radioButtons(
                ns("taxIdSource"), NULL,
                choices = list(
                    "From PhyloProfile taxonomy database" = "db",
                    "Upload own mapping file" = "file"
                ),
                selected = "db", inline = TRUE
            ),
            conditionalPanel(
                condition = "input.taxIdSource == 'db'", ns = ns,
                actionButton(
                    ns("getTaxIds"), "Get taxonomy IDs",
                    icon = icon("magnifying-glass"), class = "btn btn-primary"
                ),
                br(), br(),
                uiOutput(ns("taxIdResult.ui"))
            ),
            conditionalPanel(
                condition = "input.taxIdSource == 'file'", ns = ns,
                fileInput(
                    ns("mappingFile"), NULL,
                    accept = c(".tsv", ".txt", "text/tab-separated-values")
                ),
                helpText(
                    "Example mapping file (tab-delimited, with a header",
                    "line). Taxon names must match the taxa in the input;",
                    "IDs can be given as 7029 or ncbi7029:"
                ),
                uiOutput(ns("mappingExample.ui"))
            )
        ),
        parserSection(
            "3. Output options",
            fluidRow(
                column(
                    4,
                    radioButtons(
                        ns("orthoIdFormat"), "Ortholog ID style:",
                        choices = list(
                            "Sequence ID only" = "plain",
                            "geneID|ABBR@taxID@1|seqID|1" = "full"
                        ),
                        selected = "plain"
                    )
                ),
                column(4, uiOutput(ns("unmapped.ui"))),
                column(4, uiOutput(ns("omaIdType.ui")))
            )
        ),
        parserSection(
            "4. Convert",
            actionButton(
                ns("doParse"), "Do convert",
                icon = icon("wrench"), class = "btn btn-warning"
            ),
            br(), br(),
            verbatimTextOutput(ns("parseStatus")),
            tableOutput(ns("parsedPreview")),
            uiOutput(ns("download.ui"))
        )
    )
}

#' Boxed section of the input parser UI with a bold title
#' @param title section title
#' @param ... UI elements of the section
#' @return UI elements
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
parserSection <- function(title, ...) {
    div(
        class = "well",
        style = "padding: 10px 15px 5px 15px; margin-bottom: 15px;",
        h4(strong(title), style = "margin-top: 0; margin-bottom: 12px;"),
        ...
    )
}

#' Input parser module: server part
#'
#' @param input,output,session module parameters
#' @param inputType reactive with the selected input format ("orthofinder",
#' "sonicparanoid", "omahog" or "omagroup")
#' @param ncbiTaxDb reactive with the local NCBI taxonomy DB (data.frame
#' with ncbiID, fullName, rank, parentID), used to find taxonomy IDs for the
#' input taxa when the user does not upload a mapping file
#' @param roots roots for the folder selection (shinyFiles); default home
#' folder and system volumes
#' @return reactive with the converted data (NULL if nothing converted yet)
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
inputParser <- function(
        input, output, session, inputType, ncbiTaxDb = NULL,
        roots = c(Home = fs::path_home(), shinyFiles::getVolumes()())
) {
    ns <- session$ns
    isOma <- reactive(inputType() %in% c("omahog", "omagroup"))
    getTaxDb <- function() {
        if (is.function(ncbiTaxDb)) ncbiTaxDb() else ncbiTaxDb
    }
    
    output$unmapped.ui <- renderUI({
        choices <- list(
            "Stop with an error" = "error",
            "Skip them with a warning" = "warn"
        )
        if (identical(input$taxIdSource, "db"))
            choices <- c(choices, list("Use suggested pseudo IDs" = "pseudo"))
        radioButtons(
            ns("unmapped"), "Taxa without taxonomy ID:",
            choices = choices, selected = "error"
        )
    })
    
    # format specific inputs ---------------------------------------------------
    output$title.ui <- renderUI({
        title <- switch(
            inputType(),
            "orthofinder" = "OrthoFinder",
            "sonicparanoid" = "SonicParanoid",
            "omahog" = "OMA HOGs",
            "omagroup" = "OMA groups"
        )
        h3(
            paste("Convert", title, "output to PhyloProfile format"),
            style = "margin-top: 0; margin-bottom: 15px;"
        )
    })
    
    output$formatInput.ui <- renderUI({
        switch(
            inputType(),
            "orthofinder" = orthoTableInputUI(
                ns, "OrthoFinder orthogroup table (e.g. Orthogroups.tsv)"
            ),
            "sonicparanoid" = orthoTableInputUI(
                ns, paste(
                    "SonicParanoid flat ortholog group table",
                    "(e.g. flat.ortholog_groups.tsv)"
                )
            ),
            "omahog" = omaInputUI(ns, "HOG"),
            "omagroup" = omaInputUI(ns, "OMA group")
        )
    })
    
    output$omaIdType.ui <- renderUI({
        req(isOma())
        radioButtons(
            ns("omaIdType"), "Use as ortholog ID:",
            choices = list(
                "Protein ID" = "proteinID",
                "Cross reference" = "crossRef"
            ),
            selected = "proteinID"
        )
    })
    
    output$mappingExample.ui <- renderUI({
        example <- switch(
            inputType(),
            "orthofinder" = c(
                "species\tncbiID", "Acyrthosiphon_pisum\t7029",
                "Myzus_persicae\t13164", "Rhodnius_prolixus\t13249"
            ),
            "sonicparanoid" = c(
                "proteome\tncbiID", "GCF_000005845.2.faa\t511145",
                "GCF_000006765.1.faa\t208964", "GCF_000009045.1.faa\t224308"
            ),
            c(
                "taxon\tncbiID", "Homo sapiens\t9606",
                "MOUSE\t10090", "Bos taurus\t9913"
            )
        )
        note <- switch(
            inputType(),
            "orthofinder" = NULL,
            "sonicparanoid" = paste(
                "The FASTA extension (.faa, .fa, .fasta) is optional,",
                "e.g. GCF_000005845.2 also works."
            ),
            paste(
                "Taxa can be given by their name (column Taxon of the OMA",
                "file) or by their 5-letter OMA species code",
                "(e.g. MOUSE for MOUSE29524)."
            )
        )
        tagList(
            tags$pre(
                style = "font-size: 85%; max-width: 400px;",
                paste(example, collapse = "\n")
            ),
            if (!is.null(note)) helpText(note)
        )
    })
    
    # OMA folder selection -----------------------------------------------------
    shinyFiles::shinyDirChoose(input, "omaDir", roots = roots, session = session)
    omaDirPath <- reactive({
        req(is.list(input$omaDir))
        dirPath <- shinyFiles::parseDirPath(roots, input$omaDir)
        if (length(dirPath) == 0) return(NULL)
        as.character(dirPath)
    })
    omaDirFiles <- reactive({
        req(omaDirPath())
        list.files(
            omaDirPath(), pattern = "\\.(csv|tsv|txt)$",
            full.names = TRUE, ignore.case = TRUE
        )
    })
    output$omaDirInfo <- renderText({
        req(omaDirPath())
        nFiles <- length(omaDirFiles())
        paste0(
            omaDirPath(), "\n", nFiles, " file(s) found",
            if (nFiles > 0) paste0(
                ": ",
                paste(utils::head(basename(omaDirFiles()), 5), collapse = ", "),
                if (nFiles > 5) ", ..."
            )
        )
    })
    
    # default group name = name of the uploaded OMA file
    observeEvent(input$omaFile, {
        updateTextInput(
            session, "omaGroupName",
            value = tools::file_path_sans_ext(input$omaFile$name)
        )
    })
    
    # OMA input files (without group names)
    omaFilePaths <- function() {
        if (input$omaInputMode == "file") {
            if (is.null(input$omaFile)) stop("Please upload an OMA file!")
            return(input$omaFile$datapath)
        }
        dirFiles <- tryCatch(omaDirFiles(), error = function(e) NULL)
        if (is.null(dirFiles)) stop("Please select a folder with OMA files!")
        if (length(dirFiles) == 0)
            stop("No .csv, .tsv or .txt files found in the folder!")
        dirFiles
    }
    
    # taxa of the current input: data.frame with taxon and code (OMA only)
    getInputTaxa <- function() {
        if (isOma()) {
            omaDf <- readOmaInput(omaFilePaths())
            taxa <- omaDf[!duplicated(omaDf$taxon), c("taxon", "code")]
            return(taxa)
        }
        if (is.null(input$inputFile))
            stop("Please upload an ortholog group table!")
        header <- readLines(input$inputFile$datapath, n = 1, warn = FALSE)
        taxa <- trimws(strsplit(sub("\r$", "", header), "\t")[[1]][-1])
        if (inputType() == "sonicparanoid")
            taxa <- sub("\\.(faa|fa|fasta|fas|pep)$", "", taxa)
        data.frame(taxon = taxa, code = NA_character_, stringsAsFactors = FALSE)
    }
    
    # search taxonomy IDs in the PhyloProfile taxonomy DB (only on click) ------
    taxIdResult <- reactiveVal(NULL)
    taxIdMsg <- reactiveVal("")
    
    # reset the search result when the input changes
    observeEvent(
        list(inputType(), input$inputFile, input$omaInputMode,
             input$omaFile, input$omaDir),
        {
            taxIdResult(NULL)
            taxIdMsg("")
        },
        ignoreInit = TRUE
    )
    
    observeEvent(input$getTaxIds, {
        taxIdResult(NULL)
        taxIdMsg("")
        res <- tryCatch({
            taxDb <- getTaxDb()
            if (is.null(taxDb) || nrow(taxDb) == 0)
                stop("PhyloProfile taxonomy database is not available!")
            taxa <- getInputTaxa()
            withProgress(message = "Searching taxonomy IDs...", value = 0.5, {
                hits <- getTaxonIdsFromDb(taxa$taxon, taxDb)
            })
            hits$code <- taxa$code
            hits
        }, error = function(e) {
            taxIdMsg(paste("ERROR:", conditionMessage(e)))
            NULL
        })
        taxIdResult(res)
    })
    
    output$taxIdResult.ui <- renderUI({
        if (taxIdMsg() != "")
            return(div(class = "alert alert-danger", taxIdMsg()))
        taxMap <- taxIdResult()
        req(taxMap)
        nType <- table(factor(
            taxMap$matchType,
            levels = c("exact", "strain", "species name", "parent name",
                       "not found")
        ))
        nType <- nType[nType > 0]
        nMissing <- sum(is.na(taxMap$ncbiID))
        tagList(
            div(
                class = if (nMissing > 0) "alert alert-warning" else
                    "alert alert-success",
                style = "padding: 6px 8px; margin: 5px 0;",
                paste0(
                    nrow(taxMap), " taxa in the input. Matched by: ",
                    paste0(names(nType), " (", nType, ")", collapse = ", "),
                    "."
                ),
                if (nMissing > 0) tagList(
                    br(),
                    paste(
                        nMissing, "taxa have no taxonomy ID. You can use the",
                        "suggested pseudo IDs, skip these taxa, or add their",
                        "IDs to the downloaded taxon mapping and upload it",
                        "as your own mapping file."
                    )
                )
            ),
            DT::dataTableOutput(ns("taxonMappingTable")),
            downloadButton(ns("downloadTaxonMapping"), "Download taxon mapping"),
            downloadButton(ns("downloadTaxonTable"), "Download table"),
            helpText(
                "Taxon mapping: mapping file for PhyloProfile (taxon name,",
                "NCBI ID). Taxa without ID get their suggested pseudo ID if",
                "'Use suggested pseudo IDs' is selected, otherwise an empty ID."
            )
        )
    })
    
    # taxon table as shown in the app (taxa without ID first)
    taxonTable <- reactive({
        taxMap <- taxIdResult()
        req(taxMap)
        out <- data.frame(
            Taxon = taxMap$taxon,
            "NCBI ID" = taxMap$ncbiID,
            "Matched name" = taxMap$matchedName,
            "Matched by" = taxMap$matchType,
            "Suggested pseudo ID" = taxMap$pseudoID,
            check.names = FALSE, stringsAsFactors = FALSE
        )
        if (isOma())
            out <- cbind(out[1], "OMA code" = taxMap$code, out[-1])
        out[order(is.na(taxMap$pseudoID)), , drop = FALSE]
    })
    
    output$taxonMappingTable <- DT::renderDataTable({
        DT::datatable(
            taxonTable(), rownames = FALSE,
            options = list(
                pageLength = 10, lengthChange = FALSE, ordering = TRUE
            )
        )
    })
    
    # mapping file as PhyloProfile requires: taxon name and NCBI ID
    output$downloadTaxonMapping <- downloadHandler(
        filename = function() "taxon_mapping.tsv",
        content = function(file) {
            taxMap <- taxIdResult()
            ids <- taxMap$ncbiID
            if (identical(input$unmapped, "pseudo"))
                ids[is.na(ids)] <- taxMap$pseudoID[is.na(ids)]
            out <- data.frame(
                taxon = taxMap$taxon,
                ncbiID = ifelse(is.na(ids), "", ids),
                stringsAsFactors = FALSE
            )
            utils::write.table(
                out, file, sep = "\t", quote = FALSE, row.names = FALSE
            )
        }
    )
    
    # the taxon table as shown in the app
    output$downloadTaxonTable <- downloadHandler(
        filename = function() "taxon_mapping_table.tsv",
        content = function(file) {
            out <- taxonTable()
            out[is.na(out)] <- ""
            utils::write.table(
                out, file, sep = "\t", quote = FALSE, row.names = FALSE
            )
        }
    )
    
    # conversion ---------------------------------------------------------------
    parsedData <- reactiveVal(NULL)
    parseMsg <- reactiveVal("")
    
    # reset result when inputs change, so an old result is never downloaded
    observeEvent(
        list(inputType(), input$inputFile, input$mappingFile,
             input$orthoIdFormat, input$unmapped, input$omaInputMode,
             input$omaFile, input$omaGroupName, input$omaDir,
             input$omaIdType, input$taxIdSource, taxIdResult()),
        {
            parsedData(NULL)
            parseMsg("")
        },
        ignoreInit = TRUE
    )
    
    # taxonomy ID mapping for the parsers (data.frame or mapping file path)
    getMapping <- function() {
        if (input$taxIdSource == "file") {
            if (is.null(input$mappingFile))
                stop("Please upload a taxon to taxonomy ID mapping file!")
            return(list(mapping = input$mappingFile$datapath,
                        unmapped = input$unmapped))
        }
        taxMap <- taxIdResult()
        if (is.null(taxMap))
            stop("Please click on 'Get taxonomy IDs' first!")
        unmapped <- input$unmapped
        if (unmapped == "pseudo") {
            toPseudo <- is.na(taxMap$ncbiID)
            if (any(toPseudo)) {
                taxMap$ncbiID[toPseudo] <- taxMap$pseudoID[toPseudo]
                warning(
                    sum(toPseudo), " taxa use pseudo taxonomy IDs: ",
                    paste0(
                        taxMap$taxon[toPseudo], " (",
                        taxMap$pseudoID[toPseudo], ")", collapse = ", "
                    ),
                    call. = FALSE
                )
            }
            unmapped <- "error"
        }
        list(
            mapping = data.frame(
                taxon = taxMap$taxon, ncbiID = taxMap$ncbiID,
                stringsAsFactors = FALSE
            ),
            unmapped = unmapped
        )
    }
    
    runParser <- function() {
        if (!isOma()) {
            if (is.null(input$inputFile))
                stop("Please upload an ortholog group table!")
            mapping <- getMapping()
            parser <- switch(
                inputType(),
                "orthofinder" = orthofinder2pp,
                "sonicparanoid" = sonicparanoid2pp
            )
            return(parser(
                orthoFile = input$inputFile$datapath,
                mappingFile = mapping$mapping,
                orthoIdFormat = input$orthoIdFormat,
                unmapped = mapping$unmapped
            ))
        }
        omaFiles <- omaFilePaths()
        groupNames <- NULL
        if (input$omaInputMode == "file") {
            groupNames <- trimws(input$omaGroupName)
            if (groupNames == "")
                groupNames <- tools::file_path_sans_ext(input$omaFile$name)
        }
        mapping <- getMapping()
        oma2pp(
            omaFiles = omaFiles,
            mappingFile = mapping$mapping,
            groupNames = groupNames,
            orthoIdType = input$omaIdType,
            orthoIdFormat = input$orthoIdFormat,
            unmapped = mapping$unmapped
        )
    }
    
    observeEvent(input$doParse, {
        parsedData(NULL)
        warnMsg <- character(0)
        result <- withCallingHandlers(
            tryCatch(
                runParser(),
                error = function(e) {
                    parseMsg(paste("ERROR:", conditionMessage(e)))
                    NULL
                }
            ),
            warning = function(w) {
                warnMsg <<- c(warnMsg, paste("WARNING:", conditionMessage(w)))
                invokeRestart("muffleWarning")
            }
        )
        if (is.null(result)) return()
        parsedData(result)
        parseMsg(paste(
            c(
                paste0(
                    "Done! ", nrow(result), " orthologs in ",
                    length(unique(result$geneID)), " groups from ",
                    length(unique(result$ncbiID)), " taxa."
                ),
                warnMsg,
                if (nrow(result) > 0)
                    paste("First", min(10, nrow(result)), "lines:")
            ),
            collapse = "\n"
        ))
    })
    
    output$parseStatus <- renderText(parseMsg())
    
    output$parsedPreview <- renderTable({
        req(parsedData())
        utils::head(parsedData(), 10)
    })
    
    # download -----------------------------------------------------------------
    output$download.ui <- renderUI({
        req(parsedData())
        req(nrow(parsedData()) > 0)
        downloadButton(ns("download"), "Save PhyloProfile file")
    })
    
    output$download <- downloadHandler(
        filename = function() {
            inName <- if (!isOma()) {
                tools::file_path_sans_ext(input$inputFile$name)
            } else if (input$omaInputMode == "file") {
                unique(parsedData()$geneID)[1]
            } else {
                basename(omaDirPath())
            }
            paste0(inName, ".phyloprofile.long")
        },
        content = function(file) {
            utils::write.table(
                parsedData(), file, sep = "\t", quote = FALSE,
                row.names = FALSE
            )
        }
    )
    
    return(parsedData)
}

#' UI for table based inputs (OrthoFinder, SonicParanoid)
#' @param ns namespace function of the module
#' @param description short description of the expected input file
#' @return UI elements
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
orthoTableInputUI <- function(ns, description) {
    tagList(
        helpText(description),
        fileInput(
            ns("inputFile"), NULL,
            accept = c(".tsv", ".txt", "text/tab-separated-values")
        )
    )
}

#' UI for OMA inputs (HOGs, OMA groups)
#' @param ns namespace function of the module
#' @param groupType "HOG" or "OMA group", used in labels
#' @return UI elements
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
omaInputUI <- function(ns, groupType) {
    tagList(
        helpText(paste0(
            groupType, " file(s) downloaded from the OMA browser (CSV or TSV)"
        )),
        radioButtons(
            ns("omaInputMode"), NULL,
            choices = list(
                "Single file" = "file",
                "Folder with several files" = "folder"
            ),
            selected = "file", inline = TRUE
        ),
        conditionalPanel(
            condition = "input.omaInputMode == 'file'", ns = ns,
            fluidRow(
                column(
                    6,
                    fileInput(
                        ns("omaFile"), paste(groupType, "file:"),
                        accept = c(".csv", ".tsv", ".txt")
                    )
                ),
                column(
                    6,
                    textInput(
                        ns("omaGroupName"),
                        paste(groupType, "name (default: file name):"),
                        value = ""
                    )
                )
            )
        ),
        conditionalPanel(
            condition = "input.omaInputMode == 'folder'", ns = ns,
            shinyFiles::shinyDirButton(
                ns("omaDir"), "Select input folder",
                title = paste("Please select a folder with", groupType, "files"),
                buttonType = "default", class = NULL
            ),
            br(), br(),
            div(
                class = "alert alert-info",
                style = "padding: 6px 8px; margin: 5px 0;",
                paste0(
                    "All .csv, .tsv and .txt files in the folder will be ",
                    "used. Each file is one ", groupType, ", named after its ",
                    "file name (without extension)."
                )
            ),
            verbatimTextOutput(ns("omaDirInfo"))
        )
    )
}
