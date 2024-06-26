#' NCBI Taxonomy reduced data set
#'
#' @docType data
#' @name taxonNamesReduced
#' @description A list of NCBI taxonomy info (including taxon IDs, taxon names,
#' their systematic taxonomy rank and IDs of their next rank - parent IDs) for
#' 95 taxa in two experimental sets included in PhyloProfilData package.
#' @format Dataframe
#' @return A data frame with 4 columns:
#' \itemize{
#'     \item ncbiID e.g. "10090"
#'     \item fullName e.g. "Mus musculus"
#'     \item rank e.g. "species"
#'     \item parentID e.g. "862507"
#' }
#' @usage data(taxonNamesReduced)
NULL

#' NCBI ID list for experimental data sets
#'
#' @docType data
#' @name idList
#' @description Data frame, in which each row contains the complete taxonomy
#' ranks from the lowest systematic level (strain/species) upto the taxonomy
#' root and the corresponding IDs for one taxon in the experimental data sets.
#' @format Dataframe
#' @return A data frame with up to 41 columns and 95 rows corresponding to 95
#'     taxa in the 2 experimental data sets
#' @usage data(idList)
NULL

#' NCBI rank list for experimental data sets
#'
#' @docType data
#' @name rankList
#' @description Data frame, in which each row contains the complete taxonomy
#' ranks from the lowest systematic level (strain/species) upto the taxonomy
#' root for one taxon in the experimental data sets.
#' @format Dataframe
#' @return A data frame with up to 41 columns and 95 rows corresponding to 95
#'     taxa in the 2 experimental data sets
#' @usage data(rankList)
NULL

#' Taxonomy matrix for experimental data sets
#'
#' @docType data
#' @name taxonomyMatrix
#' @description Data frame containing the fully aligned taxonomy IDs of 95 taxa
#' in the experimental data sets. By talking into account both the defined
#' ranks (e.g. strain, This data is used for clustering and then creating a
#' taxon tree. It is used also for cross-linking between different taxonomy
#' ranks within a taxon.
#' @format Dataframe
#' @return A data frame with up to 149 columns and 95 rows corresponding to 95
#'     taxa in the 2 experimental data sets
#' @usage data(taxonomyMatrix)
NULL

#' An example of a taxonomy tree in newick format
#'
#' @docType data
#' @name ppTree
#' @format Dataframe
#' @return A data frame with only one entry
#' \itemize{
#'     \item V1 tree in newick format
#' }
#' @usage data(ppTree)
NULL

#' An example of a taxonomy matrix
#'
#' @docType data
#' @name ppTaxonomyMatrix
#' @format Dataframe
#' @return A data frame with 10 rows and 162 variables:
#' \itemize{
#'     \item abbrName  e.g. "ncbi10090"
#'     \item ncbiID e.g. "10090"
#'     \item fullName e.g. "Mus musculus"
#'     \item strain e.g. "10090"
#'     ...
#' }
#' @usage data(ppTaxonomyMatrix)
NULL

#' An example of a raw long input file
#'
#' @docType data
#' @name mainLongRaw
#' @format Dataframe
#' @return A data frame with 168 rows and 5 variables:
#' \itemize{
#'     \item geneID  Seed or ortholog group ID, e.g. "100136at6656"
#'     \item ncbiID  Taxon ID, e.g. "ncbi36329"
#'     \item orthoID Ortholog ID, e.g. "100136at6656|PLAF7@36329@1|Q8ILT8|1"
#'     \item FAS_F First additional variable
#'     \item FAS_B Second additional variable
#' }
#' @usage data(mainLongRaw)
NULL

#' An example of a raw long input file together with the taxonomy info
#'
#' @docType data
#' @name profileWithTaxonomy
#' @format Dataframe
#' @return A data frame with 20 rows and 12 variables:
#' \itemize{
#'     \item geneID Seed or ortholog group ID, e.g. "OG_1017"
#'     \item ncbiID Taxon ID, e.g. "ncbi176299"
#'     \item orthoID Ortholog ID, e.g. "A.fabrum@176299@1582"
#'     \item var1 First additional variable
#'     \item var2 Second additional variable
#'     \item paralog Number of co-orthologs in the current taxon
#'     \item abbrName e.g. "ncbi176299"
#'     \item taxonID Taxon ID, e.g. "176299"
#'     \item fullName Full taxon name, e.g. "Agrobacterium fabrum str. C58"
#'     \item supertaxonID Supertaxon ID (only different than ncbiID in case
#'     working with higher taxonomy rank than input's)
#'     \item supertaxon Name of the corresponding supertaxon
#'     \item rank Rank of the supertaxon
#' }
#' @usage data(profileWithTaxonomy)
NULL

#' An example of a fully processed phylogenetic profile
#'
#' @docType data
#' @name fullProcessedProfile
#' @format Dataframe
#' @return A data frame with 168 rows and 14 variables:
#' \itemize{
#'     \item supertaxon Supertaxon name together with its ordered index, e.g.
#'     "1001_Mammalia"
#'     \item ncbiID Taxon ID, e.g. "ncbi10116"
#'     \item geneID Seed or ortholog group ID, e.g. "100136at6656"
#'     \item orthoID Ortholog ID, e.g. "100136at6656|HUMAN@9606@1|Q9UNQ2|1"
#'     \item var1 First additional variable
#'     \item var2 Second additional variable
#'     \item paralog Number of co-orthologs in the current taxon
#'     \item abbrName NCBI ID of the ortholog, e.g. "ncbi9606"
#'     \item taxonID Taxon ID of the ortholog, in this case: "0"
#'     \item fullName Full taxon name of the ortholog, e.g. "Homo sapiens"
#'     \item supertaxonID Supertaxon ID (only different than ncbiID in case
#'     working with higher taxonomy rank than input's). e.g. "40674"
#'     \item rank Rank of the supertaxon, e.g. "class"
#'     \item category "cat
#'     \item numberSpec Total number of species in each supertaxon
#' }
#' @usage data(fullProcessedProfile)
NULL

#' An example of a filtered phylogenetic profile
#'
#' @docType data
#' @name filteredProfile
#' @format Dataframe
#' @return A data frame with 168 rows and 20 variables:
#' \itemize{
#'     \item geneID Seed or ortholog group ID, e.g. "100136at6656"
#'     \item supertaxon Supertaxon name together with its ordered index, e.g.
#'     "1001_Mammalia"
#'     \item ncbiID Taxon ID, e.g. "ncbi10116"
#'     \item orthoID Ortholog ID, e.g. "100136at6656|HUMAN@9606@1|Q9UNQ2|1"
#'     \item var1 First additional variable
#'     \item var2 Second additional variable
#'     \item paralog Number of co-orthologs in the current taxon
#'     \item abbrName NCBI ID of the ortholog, e.g. "ncbi9606"
#'     \item taxonID Taxon ID of the ortholog, in this case: "0"
#'     \item fullName Full taxon name of the ortholog, e.g. "Homo sapiens"
#'     \item supertaxonID Supertaxon ID (only different than ncbiID in case
#'     working with higher taxonomy rank than input's). e.g. "40674"
#'     \item rank Rank of the supertaxon, e.g. "class"
#'     \item category "cat
#'     \item numberSpec Total number of species in each supertaxon
#'     \item taxonMod Name of supersupertaxon w/o its index, e.g. "Mammalia"
#'     \item presSpec Percentage of taxa having orthologs in each supertaxon
#'     \item presentTaxa Number of taxa that have ortho in each supertaxon
#'     \item totalTaxa Total number of taxa in each supertaxon
#'     \item mVar1 Value of the 1. variable after grouping into supertaxon
#'     \item mVar2 Value of the 2. variable after grouping into supertaxon
#' }
#' @usage data(filteredProfile)
NULL

#' An example of a final processed & filtered phylogenetic profile
#'
#' @docType data
#' @name finalProcessedProfile
#' @format Dataframe
#' @return A data frame with 88 rows and 11 variables:
#' \itemize{
#'     \item geneID Seed or ortholog group ID, e.g. "100136at6656"
#'     \item supertaxon Supertaxon name together with its ordered index, e.g.
#'     "1001_Mammalia"
#'     \item supertaxonID Supertaxon ID (only different than ncbiID in case
#'     working with higher taxonomy rank than input's). e.g. "40674"
#'     \item var1 First additional variable
#'     \item presSpec The percentage of species presenting in each supertaxon
#'     \item category "cat"
#'     \item orthoID Ortholog ID, e.g. "100136at6656|RAT@10116@1|G3V7R8|1"
#'     \item var2 Second additional variable
#'     \item paralog Number of co-orthologs in the current taxon
#'     \item presentTaxa Number of taxa that have ortho in each supertaxon
#'     \item totalTaxa Total number of taxa in each supertaxon
#' }
#' @usage data(finalProcessedProfile)
NULL
