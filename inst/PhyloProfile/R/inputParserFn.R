#' Convert an OrthoFinder orthogroup table into PhyloProfile long format
#'
#' @description Parse an OrthoFinder orthogroup table (e.g. Orthogroups.tsv,
#' first column = orthogroup ID, one column per species, each cell holding a
#' comma-separated list of sequence IDs) and convert it into the PhyloProfile
#' long format with the columns geneID, ncbiID and orthoID.
#' @param orthoFile path to the OrthoFinder table, or a data.frame with the
#' same layout.
#' @param mappingFile path to a tab-delimited species-to-taxonomy mapping file,
#' or a data.frame. It needs two columns: species name (exactly as in the
#' table header) and NCBI taxonomy ID (either "9606" or "ncbi9606").
#' Additional columns are ignored. A header line is expected. Species with
#' an empty or NA ID are treated as unmapped (see the unmapped argument).
#' @param orthoIdFormat "plain" (default) uses the sequence ID as orthoID;
#' "full" builds a PhyloProfile-style ID "geneID|ABBR@taxID@1|seqID|1",
#' where ABBR is the species abbreviation from getSpeciesAbbr()
#' (e.g. Acanthosoma_haemorrhoidale -> ACAHA, Tuberolachnus -> TUBER).
#' @param unmapped what to do with species that are missing in the mapping
#' file: "error" (default) or "warn" (drop these species with a warning).
#' @param outFile optional path; if given, the result is written as a
#' tab-delimited file.
#' @return A data.frame with the columns geneID, ncbiID and orthoID.
#' @examples
#' \dontrun{
#' pp <- orthofinder2pp("Orthogroups.tsv", "species2taxid.tsv")
#' }
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
orthofinder2pp <- function(
        orthoFile, mappingFile, orthoIdFormat = c("plain", "full"),
        unmapped = c("error", "warn"), outFile = NULL
) {
    orthoIdFormat <- match.arg(orthoIdFormat)
    unmapped <- match.arg(unmapped)
    ortho <- readGroupTable(orthoFile)
    mapping <- readTaxMapping(mappingFile)
    pp <- groupTable2pp(ortho, mapping, orthoIdFormat, unmapped)
    writePPOutput(pp, outFile)
    return(pp)
}

#' Convert a SonicParanoid ortholog group table into PhyloProfile long format
#'
#' @description Parse the SonicParanoid flat ortholog group table
#' (flat.ortholog_groups.tsv: first column = group_id, one column per input
#' proteome, each cell holding a comma-separated list of sequence IDs or "*"
#' if the proteome has no ortholog in the group) and convert it into the
#' PhyloProfile long format with the columns geneID, ncbiID and orthoID.
#' @param orthoFile path to the SonicParanoid flat ortholog group table, or a
#' data.frame with the same layout.
#' @param mappingFile path to a tab-delimited mapping file, or a data.frame,
#' with the proteome names (column names of the SonicParanoid table) and their
#' NCBI taxonomy IDs. Names can be given with or without the FASTA file
#' extension (.faa, .fa, .fasta, .fas, .pep), e.g. "GCF_000069205.1.faa"
#' or "GCF_000069205.1".
#' @param orthoIdFormat "plain" (default) or "full", see orthofinder2pp().
#' @param unmapped "error" (default) or "warn", see orthofinder2pp().
#' @param groupPrefix text added in front of the numeric SonicParanoid group
#' IDs to create the geneIDs (default "SP", giving SP1, SP2, ...). Use "" to
#' keep the plain numbers.
#' @param outFile optional path; if given, the result is written as a
#' tab-delimited file.
#' @return A data.frame with the columns geneID, ncbiID and orthoID.
#' @examples
#' \dontrun{
#' pp <- sonicparanoid2pp("flat.ortholog_groups.tsv", "proteome2taxid.tsv")
#' }
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
sonicparanoid2pp <- function(
        orthoFile, mappingFile, orthoIdFormat = c("plain", "full"),
        unmapped = c("error", "warn"), groupPrefix = "SP", outFile = NULL
) {
    orthoIdFormat <- match.arg(orthoIdFormat)
    unmapped <- match.arg(unmapped)
    ortho <- readGroupTable(orthoFile, strict = TRUE)
    ## "*" marks a proteome without orthologs in the group
    ortho[-1] <- lapply(ortho[-1], function(x) {
        x[x == "*"] <- ""
        x
    })
    ortho[[1]] <- paste0(groupPrefix, ortho[[1]])
    ## proteome names without FASTA extension, in table and mapping
    fastaExt <- "\\.(faa|fa|fasta|fas|pep)$"
    colnames(ortho)[-1] <- sub(fastaExt, "", colnames(ortho)[-1])
    mapping <- readTaxMapping(mappingFile)
    mapping$species <- sub(fastaExt, "", mapping$species)
    pp <- groupTable2pp(ortho, mapping, orthoIdFormat, unmapped)
    writePPOutput(pp, outFile)
    return(pp)
}

#' Convert OMA HOG / OMA group files into PhyloProfile long format
#'
#' @description Parse one or more OMA HOG or OMA group files downloaded from
#' the OMA browser (CSV or TSV, with the columns "Protein ID",
#' "Cross reference" and "Taxon"; other columns are ignored) and convert them
#' into the PhyloProfile long format with the columns geneID, ncbiID and
#' orthoID. Each file is one ortholog group (geneID).
#' @param omaFiles path to a single file, a vector of file paths, or a folder.
#' For a folder, all .csv, .tsv and .txt files in it are used.
#' @param mappingFile optional (NULL) path to a tab-delimited mapping file, or
#' a data.frame, with taxa and their NCBI taxonomy IDs ("9606" or
#' "ncbi9606"). Taxa can be given by their name as in the "Taxon" column or
#' by their 5-letter OMA species code (e.g. HUMAN). The mapping file has
#' priority over the search in taxDb.
#' @param taxDb optional (NULL) data.frame of the local NCBI taxonomy DB with
#' (at least) the columns ncbiID and fullName. Taxa not found in the mapping
#' file are searched in it with matchTaxaNcbi().
#' @param groupNames optional names for the groups, one per file. By default
#' (NULL) the file names without extension are used.
#' @param orthoIdType "proteinID" (default) uses the OMA protein ID as
#' ortholog ID, "crossRef" uses the cross reference. Proteins without a
#' cross reference keep their OMA protein ID (with a warning).
#' @param orthoIdFormat "plain" (default) or "full". "full" builds
#' "geneID|CODE@taxID@1|seqID|1" with CODE = the OMA species code.
#' @param unmapped what to do with taxa without taxonomy ID: "error"
#' (default), "warn" (drop them with a warning) or "pseudo" (use pseudo IDs
#' above the highest ID in taxDb, see getOmaTaxonMapping()).
#' @return A data.frame with the columns geneID, ncbiID and orthoID.
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @examples
#' \dontrun{
#' pp <- oma2pp("hog_F0848137.csv", taxDb = currentNCBIinfo())
#' pp <- oma2pp("oma_hogs/", "taxon2taxid.tsv", orthoIdType = "crossRef")
#' }
oma2pp <- function(
        omaFiles, mappingFile = NULL, taxDb = NULL, groupNames = NULL,
        orthoIdType = c("proteinID", "crossRef"),
        orthoIdFormat = c("plain", "full"),
        unmapped = c("error", "warn", "pseudo")
) {
    orthoIdType <- match.arg(orthoIdType)
    orthoIdFormat <- match.arg(orthoIdFormat)
    unmapped <- match.arg(unmapped)
    if (is.null(mappingFile) && is.null(taxDb))
        stop("Please provide a mapping file and/or the NCBI taxonomy DB")
    
    omaDf <- readOmaInput(omaFiles, groupNames)
    if (nrow(omaDf) == 0) {
        warning("No orthologs found in the input files", call. = FALSE)
        return(data.frame(
            geneID = character(0), ncbiID = character(0),
            orthoID = character(0), stringsAsFactors = FALSE
        ))
    }
    
    ## taxonomy IDs for all taxa
    taxMap <- getOmaTaxonMapping(omaDf, mappingFile, taxDb)
    missingTax <- taxMap$taxon[is.na(taxMap$ncbiID)]
    if (length(missingTax) > 0) {
        msg <- paste0(
            length(missingTax), " taxa without NCBI taxonomy ID: ",
            paste(missingTax, collapse = ", ")
        )
        if (unmapped == "error") stop(msg)
        if (unmapped == "warn") {
            warning(msg, call. = FALSE)
        } else {
            if (all(is.na(taxMap$pseudoID)))
                stop("Pseudo IDs need the NCBI taxonomy DB (taxDb)")
            toPseudo <- is.na(taxMap$ncbiID)
            taxMap$ncbiID[toPseudo] <- taxMap$pseudoID[toPseudo]
            warning(
                length(missingTax), " taxa use pseudo taxonomy IDs: ",
                paste0(
                    taxMap$taxon[toPseudo], " (", taxMap$pseudoID[toPseudo],
                    ")", collapse = ", "
                ),
                call. = FALSE
            )
        }
    }
    omaDf$taxId <- taxMap$ncbiID[match(omaDf$taxon, taxMap$taxon)]
    omaDf <- omaDf[!is.na(omaDf$taxId), , drop = FALSE]
    
    ## ortholog IDs
    seqIds <- omaDf$proteinID
    if (orthoIdType == "crossRef") {
        noRef <- omaDf$crossRef == ""
        if (any(noRef))
            warning(
                sum(noRef), " proteins without cross reference keep their ",
                "OMA protein ID: ",
                paste(utils::head(omaDf$proteinID[noRef], 10), collapse = ", "),
                if (sum(noRef) > 10) ", ...", call. = FALSE
            )
        seqIds[!noRef] <- omaDf$crossRef[!noRef]
    }
    if (orthoIdFormat == "full") {
        seqIds <- paste0(
            omaDf$geneID, "|", omaDf$code, "@", omaDf$taxId, "@1|",
            seqIds, "|1"
        )
    }
    pp <- data.frame(
        geneID = omaDf$geneID,
        ncbiID = paste0("ncbi", omaDf$taxId),
        orthoID = seqIds,
        stringsAsFactors = FALSE
    )
    rownames(pp) <- NULL
    return(pp)
}

#' Read OMA HOG / OMA group input files
#'
#' @param omaFiles path to a single file, a vector of file paths, or a folder
#' (all .csv, .tsv and .txt files in it are used)
#' @param groupNames optional group names, one per file (default: file names
#' without extension)
#' @return data.frame with geneID, proteinID, crossRef, taxon and code (the
#' 5-letter OMA species code)
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
readOmaInput <- function(omaFiles, groupNames = NULL) {
    if (length(omaFiles) == 1 && dir.exists(omaFiles)) {
        omaFiles <- list.files(
            omaFiles, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE,
            ignore.case = TRUE
        )
        if (length(omaFiles) == 0)
            stop("No .csv, .tsv or .txt files found in the folder")
    }
    notFound <- omaFiles[!file.exists(omaFiles)]
    if (length(notFound) > 0)
        stop("File not found: ", paste(notFound, collapse = ", "))
    if (is.null(groupNames))
        groupNames <- tools::file_path_sans_ext(basename(omaFiles))
    groupNames <- trimws(as.character(groupNames))
    if (length(groupNames) != length(omaFiles))
        stop("Number of group names does not match the number of files")
    if (any(groupNames == "") || anyDuplicated(groupNames))
        stop("Group names must be non-empty and unique")
    omaDf <- do.call(rbind, lapply(seq_along(omaFiles), function(i) {
        df <- readOmaFile(omaFiles[i])
        if (nrow(df) == 0) return(NULL)
        df$geneID <- groupNames[i]
        df
    }))
    if (is.null(omaDf)) {
        omaDf <- data.frame(
            proteinID = character(0), crossRef = character(0),
            taxon = character(0), geneID = character(0),
            stringsAsFactors = FALSE
        )
    }
    omaDf$code <- substr(omaDf$proteinID, 1, 5)
    return(omaDf)
}

#' Get NCBI taxonomy IDs for the taxa of OMA data
#'
#' @description Taxa are first looked up in the user's mapping file (by
#' taxon name or OMA species code), then in the local NCBI taxonomy DB with
#' matchTaxaNcbi(). Taxa that are still missing get a suggested pseudo ID,
#' counting up from the highest ID in taxDb.
#' @param omaDf data.frame from readOmaInput()
#' @param mappingFile optional mapping file path or data.frame (see oma2pp())
#' @param taxDb optional local NCBI taxonomy DB (columns ncbiID, fullName)
#' @return data.frame with one row per taxon: taxon, code, nProteins,
#' ncbiID (NA if not found), matchedName, matchType ("mapping file",
#' "exact", "strain", "species name", "parent name" or "not found") and
#' pseudoID (only for taxa not found, NA otherwise)
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
getOmaTaxonMapping <- function(omaDf, mappingFile = NULL, taxDb = NULL) {
    taxa <- unique(omaDf[, c("taxon", "code")])
    taxa <- taxa[!duplicated(taxa$taxon), , drop = FALSE]
    taxa$nProteins <- as.integer(table(omaDf$taxon)[taxa$taxon])
    taxa$ncbiID <- NA_character_
    taxa$matchedName <- NA_character_
    taxa$matchType <- "not found"
    
    ## 1. user mapping file: by taxon name, then by OMA code
    if (!is.null(mappingFile)) {
        mapping <- readTaxMapping(mappingFile)
        mapping$species <- gsub("\\s+", " ", mapping$species)
        taxIds <- stats::setNames(mapping$taxId, mapping$species)
        ids <- unname(taxIds[taxa$taxon])
        ids[is.na(ids)] <- unname(taxIds[taxa$code[is.na(ids)]])
        found <- !is.na(ids)
        taxa$ncbiID[found] <- ids[found]
        taxa$matchType[found] <- "mapping file"
    }
    
    ## 2. local NCBI taxonomy DB
    if (!is.null(taxDb) && any(is.na(taxa$ncbiID))) {
        toSearch <- is.na(taxa$ncbiID)
        hits <- matchTaxaNcbi(taxa$taxon[toSearch], taxDb)
        taxa$ncbiID[toSearch] <- hits$ncbiID
        taxa$matchedName[toSearch] <- hits$matchedName
        taxa$matchType[toSearch] <- hits$matchType
    }
    
    ## 3. pseudo IDs for the remaining taxa
    taxa$pseudoID <- NA_character_
    missingTax <- is.na(taxa$ncbiID)
    if (!is.null(taxDb) && any(missingTax)) {
        maxId <- max(suppressWarnings(as.numeric(taxDb$ncbiID)), na.rm = TRUE)
        taxa$pseudoID[missingTax] <- format(
            maxId + seq_len(sum(missingTax)), scientific = FALSE, trim = TRUE
        )
    }
    rownames(taxa) <- NULL
    return(taxa)
}

#' Get NCBI taxonomy IDs for taxon names from the local taxonomy DB
#'
#' @description Search the taxon names in the local NCBI taxonomy DB with
#' matchTaxaNcbi() ("_" in names are treated as spaces, e.g. for OrthoFinder
#' species names) and suggest pseudo IDs for taxa that are not found,
#' counting up from the highest ID in the DB.
#' @param taxonNames character vector of taxon names
#' @param taxDb data.frame with the columns ncbiID and fullName (and
#' optionally rank)
#' @return data.frame with taxon, ncbiID (NA if not found), matchedName,
#' matchType and pseudoID (only for taxa not found, NA otherwise)
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
getTaxonIdsFromDb <- function(taxonNames, taxDb) {
    hits <- matchTaxaNcbi(gsub("_", " ", taxonNames), taxDb)
    hits$taxon <- taxonNames
    hits$pseudoID <- NA_character_
    missingTax <- is.na(hits$ncbiID)
    if (any(missingTax)) {
        maxId <- max(suppressWarnings(as.numeric(taxDb$ncbiID)), na.rm = TRUE)
        hits$pseudoID[missingTax] <- format(
            maxId + seq_len(sum(missingTax)), scientific = FALSE, trim = TRUE
        )
    }
    rownames(hits) <- NULL
    return(hits)
}

#' Find NCBI taxonomy IDs for taxon names in the local NCBI taxonomy DB
#'
#' @description Names are compared case-insensitively and with collapsed
#' whitespace, in this order:
#' 1. "exact": the whole name, e.g. "Homo sapiens".
#' 2. "strain": for OMA names with strain info in brackets, e.g.
#' "Myxococcus stipitatus (strain DSM 14675 / JCM 12634 / Mx s8)", each
#' designation is tried as "<species> <designation>" (also with "str." or
#' "strain" in between), ignoring spaces and punctuation, so that
#' "Myxococcus stipitatus DSM 14675" is found.
#' 3. "species name": the name without the part in brackets
#' ("Myxococcus stipitatus").
#' 4. "parent name": the name with the last word(s) removed, down to two
#' words (e.g. "Saimiri boliviensis boliviensis" -> "Saimiri boliviensis").
#' If a name has several hits, the one with the lowest rank is preferred
#' (strain/subspecies before species before others).
#' @param taxonNames character vector of taxon names
#' @param taxDb data.frame with the columns ncbiID and fullName (and
#' optionally rank)
#' @return data.frame with taxon, ncbiID (NA if not found), matchedName and
#' matchType
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
matchTaxaNcbi <- function(taxonNames, taxDb) {
    normName <- function(x) tolower(gsub("\\s+", " ", trimws(x)))
    alnum <- function(x) gsub("[^a-z0-9]", "", tolower(x))
    ## only DB entries of the same genera (first words) are needed
    genera <- unique(sub(" .*", "", normName(taxonNames)))
    dbNames <- tolower(as.character(taxDb$fullName))
    keep <- sub(" .*", "", dbNames) %in% genera
    taxDb <- taxDb[keep, , drop = FALSE]
    dbNames <- dbNames[keep]
    dbIds <- as.character(taxDb$ncbiID)
    rankOrder <- c("strain", "isolate", "serotype", "serogroup",
                   "subspecies", "varietas", "forma", "species")
    dbRank <- if (!is.null(taxDb$rank)) {
        r <- match(as.character(taxDb$rank), rankOrder)
        ifelse(is.na(r), length(rankOrder) + 1, r)
    } else rep(1, nrow(taxDb))
    
    ## best hit (lowest rank) among the DB rows with one of these names
    lookup <- function(names) {
        idx <- which(dbNames %in% names)
        if (length(idx) == 0) return(NULL)
        idx[which.min(dbRank[idx])]
    }
    
    res <- lapply(taxonNames, function(taxon) {
        name <- normName(taxon)
        out <- function(idx, type) data.frame(
            taxon = taxon, ncbiID = dbIds[idx],
            matchedName = as.character(taxDb$fullName[idx]),
            matchType = type, stringsAsFactors = FALSE
        )
        ## 1. exact
        idx <- lookup(name)
        if (!is.null(idx)) return(out(idx, "exact"))
        ## 2. strain designations in brackets
        base <- trimws(sub("\\s*\\(.*\\)\\s*$", "", name))
        if (grepl("\\(.+\\)", name) && base != "") {
            inside <- sub("^.*\\((.*)\\)\\s*$", "\\1", name)
            inside <- sub(
                "^(strain|isolate|substrain|serovar|serotype)\\s+", "", inside
            )
            designations <- trimws(strsplit(inside, "/", fixed = TRUE)[[1]])
            designations <- designations[designations != ""]
            cand <- which(startsWith(dbNames, paste0(base, " ")))
            if (length(cand) > 0 && length(designations) > 0) {
                candKeys <- alnum(dbNames[cand])
                queryKeys <- alnum(c(
                    paste(base, designations),
                    paste(base, "str", designations),
                    paste(base, "strain", designations)
                ))
                hit <- cand[candKeys %in% queryKeys]
                if (length(hit) > 0)
                    return(out(hit[which.min(dbRank[hit])], "strain"))
            }
        }
        ## 3. species name without bracket part
        if (base != name && base != "") {
            idx <- lookup(base)
            if (!is.null(idx)) return(out(idx, "species name"))
        }
        ## 4. remove last words, down to two words
        words <- strsplit(base, " ", fixed = TRUE)[[1]]
        while (length(words) > 2) {
            words <- words[-length(words)]
            idx <- lookup(paste(words, collapse = " "))
            if (!is.null(idx)) return(out(idx, "parent name"))
        }
        data.frame(
            taxon = taxon, ncbiID = NA_character_, matchedName = NA_character_,
            matchType = "not found", stringsAsFactors = FALSE
        )
    })
    do.call(rbind, res)
}

#' Create species abbreviations
#'
#' @description Abbreviate species names as the first 3 letters of the genus
#' plus the first 2 letters of the species epithet, in upper case
#' (Acanthosoma_haemorrhoidale or Acanthosoma-haemorrhoidale -> ACAHA). A
#' one-word name gives its first 5 letters (Tuberolachnus -> TUBER). Words can
#' be separated by "_", "-" or spaces; only the first two words are used.
#' @param speciesNames character vector of species names.
#' @return character vector of abbreviations.
#' @examples
#' getSpeciesAbbr(c("Acanthosoma_haemorrhoidale", "Tuberolachnus"))
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
getSpeciesAbbr <- function(speciesNames) {
    words <- strsplit(trimws(as.character(speciesNames)), "[_ -]+")
    abbr <- vapply(words, function(w) {
        w <- gsub("[^A-Za-z0-9]", "", w)
        w <- w[w != ""]
        if (length(w) == 0) return(NA_character_)
        if (length(w) == 1) return(substr(w[1], 1, 5))
        paste0(substr(w[1], 1, 3), substr(w[2], 1, 2))
    }, character(1))
    toupper(abbr)
}

# ---- internal helpers --------------------------------------------------------

#' Read a group table (group ID column + one column per taxon)
#' @param orthoFile file path or data.frame
#' @param strict TRUE: every line must have as many columns as the header
#' @return data.frame of trimmed character columns
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @keywords internal
readGroupTable <- function(orthoFile, strict = FALSE) {
    if (is.data.frame(orthoFile)) {
        ortho <- orthoFile
    } else {
        if (!file.exists(orthoFile)) stop("File not found: ", orthoFile)
        lines <- readLines(orthoFile, warn = FALSE)
        lines <- sub("\r$", "", lines)
        lines <- lines[trimws(lines) != ""]
        if (length(lines) < 2) stop("Input table has no data lines")
        if (!any(grepl("\t", lines, fixed = TRUE)))
            stop("Input table is not tab-delimited")
        fields <- strsplit(lines, "\t", fixed = TRUE)
        ## strsplit drops trailing empty fields; pad them back
        nCol <- length(fields[[1]])
        nField <- lengths(fields)
        badRows <- if (strict) which(nField != nCol) else which(nField > nCol)
        if (length(badRows) > 0)
            stop(
                "These lines have ", if (strict) "a different number of" else "more",
                " columns than the header line: ",
                paste(utils::head(badRows, 10), collapse = ", "),
                if (length(badRows) > 10) ", ..."
            )
        fields <- lapply(fields, function(x) c(x, rep("", nCol - length(x))))
        ortho <- as.data.frame(
            do.call(rbind, fields[-1]), stringsAsFactors = FALSE
        )
        colnames(ortho) <- fields[[1]]
    }
    if (ncol(ortho) < 2) stop("Input table needs at least 2 columns")
    ortho[] <- lapply(ortho, function(x) trimws(as.character(x)))
    colnames(ortho) <- trimws(colnames(ortho))
    ortho[is.na(ortho)] <- ""
    if (anyDuplicated(colnames(ortho)[-1]))
        stop("Duplicated taxon names in the header line")
    return(ortho)
}

#' Read a taxon-to-NCBI-taxonomy-ID mapping
#' @param mappingFile file path or data.frame (taxon name, taxonomy ID)
#' @return data.frame with columns species and taxId (numeric strings)
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @keywords internal
readTaxMapping <- function(mappingFile) {
    if (is.data.frame(mappingFile)) {
        mapping <- mappingFile
    } else {
        if (!file.exists(mappingFile)) stop("File not found: ", mappingFile)
        mapping <- utils::read.table(
            mappingFile, header = TRUE, sep = "\t", quote = "",
            comment.char = "", colClasses = "character", check.names = FALSE
        )
    }
    if (ncol(mapping) < 2) stop("Mapping file needs at least 2 columns")
    mapping <- data.frame(
        species = trimws(as.character(mapping[[1]])),
        taxId = sub("^ncbi", "", trimws(as.character(mapping[[2]]))),
        stringsAsFactors = FALSE
    )
    ## empty / NA IDs count as unmapped taxa
    mapping <- mapping[
        !is.na(mapping$taxId) & !(mapping$taxId %in% c("", "NA")), ,
        drop = FALSE
    ]
    if (any(!grepl("^[0-9]+$", mapping$taxId)))
        stop("Non-numeric taxonomy IDs in mapping file: ",
             paste(mapping$taxId[!grepl("^[0-9]+$", mapping$taxId)],
                   collapse = ", "))
    if (anyDuplicated(mapping$species))
        stop("Duplicated species in mapping file: ",
             paste(unique(mapping$species[duplicated(mapping$species)]),
                   collapse = ", "))
    return(mapping)
}

#' Convert a group table into PhyloProfile long format
#' @param ortho data.frame from readGroupTable() (missing cells = "")
#' @param mapping data.frame from readTaxMapping()
#' @param orthoIdFormat "plain" or "full"
#' @param unmapped "error" or "warn"
#' @return data.frame with geneID, ncbiID, orthoID
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @keywords internal
groupTable2pp <- function(ortho, mapping, orthoIdFormat, unmapped) {
    speciesCols <- colnames(ortho)[-1]
    taxIds <- stats::setNames(mapping$taxId, mapping$species)
    
    ## check that every taxon in the table has a taxonomy ID
    missingSpec <- setdiff(speciesCols, names(taxIds))
    if (length(missingSpec) > 0) {
        msg <- paste0(
            length(missingSpec), " species not found in mapping file: ",
            paste(missingSpec, collapse = ", ")
        )
        if (unmapped == "error") stop(msg)
        warning(msg, call. = FALSE)
        speciesCols <- setdiff(speciesCols, missingSpec)
    }
    
    ## wide -> long: one row per (group, taxon, sequence)
    geneIDs <- ortho[[1]]
    res <- lapply(speciesCols, function(sp) {
        seqs <- strsplit(ortho[[sp]], ",", fixed = TRUE)
        n <- lengths(seqs)
        if (sum(n) == 0) return(NULL)
        out <- data.frame(
            geneID = rep(geneIDs, n),
            ncbiID = paste0("ncbi", taxIds[[sp]]),
            orthoID = trimws(unlist(seqs, use.names = FALSE)),
            species = sp,
            stringsAsFactors = FALSE
        )
        out[out$orthoID != "", , drop = FALSE]
    })
    pp <- do.call(rbind, res)
    if (is.null(pp) || nrow(pp) == 0) {
        warning("No orthologs found in the input table", call. = FALSE)
        return(data.frame(
            geneID = character(0), ncbiID = character(0),
            orthoID = character(0), stringsAsFactors = FALSE
        ))
    }
    ## keep the original group order
    pp <- pp[order(match(pp$geneID, geneIDs)), , drop = FALSE]
    if (orthoIdFormat == "full") {
        spec <- unique(pp$species)
        abbr <- stats::setNames(getSpeciesAbbr(spec), spec)
        dupAbbr <- unique(abbr[duplicated(abbr)])
        if (length(dupAbbr) > 0)
            warning("Species sharing the same abbreviation: ",
                    paste(names(abbr)[abbr %in% dupAbbr], collapse = ", "),
                    call. = FALSE)
        pp$orthoID <- paste0(
            pp$geneID, "|", abbr[pp$species], "@",
            sub("^ncbi", "", pp$ncbiID), "@1|", pp$orthoID, "|1"
        )
    }
    pp <- pp[, c("geneID", "ncbiID", "orthoID")]
    rownames(pp) <- NULL
    return(pp)
}

#' Read one OMA HOG / OMA group file (CSV or TSV)
#' @param file path to the file
#' @return data.frame with proteinID, crossRef and taxon
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @keywords internal
readOmaFile <- function(file) {
    firstLine <- readLines(file, n = 1, warn = FALSE, encoding = "UTF-8")
    sep <- if (grepl("\t", firstLine, fixed = TRUE)) "\t" else ","
    df <- utils::read.table(
        file, header = TRUE, sep = sep, quote = "\"", comment.char = "",
        colClasses = "character", check.names = FALSE, fill = TRUE,
        fileEncoding = "UTF-8-BOM", na.strings = character(0)
    )
    colnames(df) <- trimws(colnames(df))
    required <- c("Protein ID", "Cross reference", "Taxon")
    missingCols <- setdiff(required, colnames(df))
    if (length(missingCols) > 0)
        stop(basename(file), " is not an OMA HOG / group file (missing ",
             "column: ", paste(missingCols, collapse = ", "), ")")
    df <- data.frame(
        proteinID = trimws(df[["Protein ID"]]),
        crossRef = trimws(df[["Cross reference"]]),
        taxon = gsub("\\s+", " ", trimws(df[["Taxon"]])),
        stringsAsFactors = FALSE
    )
    df[df$proteinID != "", , drop = FALSE]
}

#' Write PhyloProfile long format output
#' @param pp data.frame to write
#' @param outFile output path or NULL (nothing is written)
#' @author Vinh Tran (tran@bio.uni-frankfurt.de)
#' @keywords internal
writePPOutput <- function(pp, outFile) {
    if (!is.null(outFile)) {
        utils::write.table(
            pp, outFile, sep = "\t", quote = FALSE, row.names = FALSE
        )
    }
}
