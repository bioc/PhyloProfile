% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsePhylotree.R
\name{checkNewick}
\alias{checkNewick}
\title{Check the validity of input newick tree}
\usage{
checkNewick(tree, inputTaxonID = NULL)
}
\arguments{
\item{tree}{input newick tree}

\item{inputTaxonID}{list of all input taxon IDs for the phylogenetic profiles}
}
\value{
Possible formatting error of input tree. 0 = suitable tree for using
with PhyloProfile, 1 = missing parenthesis; 2 = missing comma;
3 = tree has singleton; or a list of taxa that do not exist in the input
phylogenetic profile.
}
\description{
Check the validity of input newick tree
}
\examples{
data("ppTree", package="PhyloProfile")
checkNewick(ppTree, c("ncbi3702", "ncbi3711", "ncbi7029"))
}
\seealso{
\code{\link{getInputTaxaID}} for getting input taxon IDs,
\code{\link{ppTree}} for an example of input tree
}
\author{
Vinh Tran tran@bio.uni-frankfurt.de
}
