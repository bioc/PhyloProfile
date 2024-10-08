# CREATE MANUAL COLOR SCHEME FOR A LIST

#' Create qualitative colours
#' @param n number of colors
#' @param light light colors TRUE or FALSE
#' @return list of n different colors
#' @source Modified based on https://gist.github.com/peterk87/6011397
#' @importFrom RColorBrewer brewer.pal
#' @examples
#' PhyloProfile:::qualitativeColours(5)

qualitativeColours <- function(n, light = FALSE) {
    # For more than 21 colours needed
    rich12equal <- c(
        "#000040", "#000093", "#0020E9", "#0076FF", "#00B8C2", "#04E466",
        "#49FB25", "#E7FD09", "#FEEA02", "#FFC200", "#FF8500", "#FF3300")
    # Qualitative colour schemes by Paul Tol
    ifelse (n >= 19 & n <= 21, return(grDevices::colorRampPalette(c(
        "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
        "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA",
        "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
        "#771122", "#AA4455", "#DD7788"))(n)),
        ifelse (n >= 16 & n <= 18, return(grDevices::colorRampPalette(c(
            "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD",
            "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77",
            "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455",
            "#DD7788"))(n)),
            ifelse (n == 15, return(grDevices::colorRampPalette(c(
                "#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA",
                "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711",
                "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
                "#771122", "#AA4455", "#DD7788"))(n)),
                ifelse(n > 12 & n <= 14, return(grDevices::colorRampPalette(c(
                    "#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7",
                    "#7BAFDE", "#4EB265", "#90C987", "#CAE0AB", "#F7EE55",
                    "#F6C141", "#F1932D", "#E8601C", "#DC050C"))(n)),
                    ifelse(n > 9 & n <= 12,
                        ifelse(
                            light,
                            return(RColorBrewer::brewer.pal(
                                n = n, name = 'Set3')),
                            return(RColorBrewer::brewer.pal(
                                n = n, name = 'Paired'))),
                        ifelse(n <= 9,
                            ifelse(
                                light,
                                return(RColorBrewer::brewer.pal(
                                    n = n, name = 'Pastel1')),
                                return(RColorBrewer::brewer.pal(
                                    n = n, name = 'Set1'))),
                            return(grDevices::colorRampPalette(rich12equal)(n))
                        )
                    )
                )
            )
        )
    )
}

#' Get color for a list of items
#' @return list of colors for each element (same elements will have the same
#' color)
#' @param x input list
#' @export
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{qualitativeColours}}
#' @examples
#' items <- c("a", "b", "c")
#' getQualColForVector(items)

getQualColForVector <- function(x = NULL) {
    if (is.null(x)) stop("Input list is NULL!")
    types <- unique(x)
    types <- types[!is.na(types)]
    typeColors <- qualitativeColours(length(types))
    colorsTypes <- as.vector(x)
    countTypes <- countColors <- 1
    while (countTypes <= length(types)) {
        if (countColors > length(typeColors)) countColors <- 1
        colorsTypes[colorsTypes == types[countTypes]] <- typeColors[countColors]
        countColors <- countColors + 1
        countTypes <- countTypes + 1
    }
    return(unlist(colorsTypes))
}


#' Check if a color pallete has enough colors for a list of items
#' @export
#' @param items vector contains list of items
#' @param pallete name of color palette
#' @return TRUE if color pallete has enough colors, otherwise FALSE
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom RColorBrewer brewer.pal.info
#' @examples
#' myItems <- rep("a",3)
#' checkColorPalette(myItems, "Set1")

checkColorPalette <- function(items, pallete = "Paired") {
    if (is.null(items)) stop("No item given!")
    colorDf <- data.frame(RColorBrewer::brewer.pal.info)
    if (!(pallete %in% row.names(colorDf)))
        stop("Given color pallete not found")
    colorDf$name <- row.names(colorDf)
    if (length(items) < colorDf$maxcolors[colorDf$name == pallete])
        return(TRUE)
    return(FALSE)
}
