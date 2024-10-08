\name{PhyloProfilenews}
\title{PhyloProfile news}
\encoding{UTF-8}

\section{Version 1.20.0}{
  \itemize{
    \item added UMAP clustering
    \item improved domain plot #110
    \item fast mode for large profile plot
    \item option to change font
    \item option to show gene names
    \item option to sort genes by defined list
    \item option to upload preprocessed data from an input folder
  }
}

\section{Version 1.18.0}{
  \itemize{
    \item option to specify version for OrthoDB
  }
}

\section{Version 1.16.5}{
  \itemize{
    \item fixed detailed plot when name of a taxon is substring of another
    \item show supertaxon names when rotate the profile plot
  }
}

\section{Version 1.16.3}{
  \itemize{
    \item fixed clustering method not applied
  }
}

\section{Version 1.16.2}{
  \itemize{
    \item check if taxonomy DB exists in cwd instead of r/library
  }
}

\section{Version 1.16.1}{
  \itemize{
    \item fixed bug wrong order in detailed plot
    \item fixed bug mVar for co-orthologs
  }
}

\section{Version 1.16.0}{
  \itemize{
    \item extended options in config file
    \item option to specify IP and port
    \item able to generate domain plot from main/customized profile tab
    \item highlight multiple genes/taxa; division lines for taxon group (#110)
    \item simplify gene IDs in detailed plot; show taxonomy hierarchy (#110)
    \item option to filter features in architecture plot (#110)
    \item option to show all input taxa (#110)
    \item added more info to domain plot (e-value, bit-score, pHMM)
    \item option to specify host and port for runPhyloProfile() fn
    \item fixed bug when lowest rank is subspecies
  }
}

\section{Version 1.14.5}{
  \itemize{
    \item improved highlighted duplicated ortho IDs
  }
}

\section{Version 1.14.4}{
  \itemize{
    \item fixed bug ordering genes
    \item auto identify lowest rank
  }
}

\section{Version 1.14.2}{
  \itemize{
    \item option to use user defined taxonomy DB (#124)
  }
}

\section{Version 1.14.1}{
  \itemize{
    \item option for highlight duplicated ortholog IDs
    \item option for cluster profiles based on ortholog IDs
    \item speedup by pre-calculating taxonomy tree (#123)
    \item disable rank selection after plotting
    \item turn on profile clustering by default
    \item changed "ftp://" to "https://" for taxdmp.zip
    \item fixed #126
  }
}

\section{Version 1.14.0}{
  \itemize{
    \item option for uploading sorted taxon list
  }
}

\section{Version 1.12.6}{
  \itemize{
    \item fixed bug ordering gene IDs when using gene categories
  }
}

\section{Version 1.12.5}{
  \itemize{
    \item fixed bug multiple entries for one (super)taxon ID
  }
}

\section{Version 1.12.4}{
  \itemize{
    \item fixed bug download data
  }
}

\section{Version 1.12.3}{
  \itemize{
    \item fixed bug customized profile not showed
  }
}

\section{Version 1.12.1}{
  \itemize{
    \item added number of co-orthologs and number of taxa in each supertaxon
    \item replaced pfam link by interpro url
  }
}

\section{Version 1.10.3}{
  \itemize{
    \item fixed bug parsing taxon names for input with more than 9999 taxa
  }
}
 
\section{Version 1.10.1}{
  \itemize{
    \item fixed bug button for tree uploading disappear #119
    \item improved parsing seq ID to create links for external DB #117
    \item fixed bug grepping domains for group comparison fn
    \item fixed bug parsing data for a list of input gene IDs #122
  }
}

\section{Version 1.8.6}{
  \itemize{
    \item option to identify sequence source
  }
}

\section{Version 1.8.5}{
  \itemize{
    \item added functions for exporting plot settings
  }
}

\section{Version 1.8.4}{
  \itemize{
    \item fixed error reading taxonNamesReduced.txt that contains "#"
    \item added functions for import and export taxonomy DB
  }
}

\section{Version 1.8.3}{
  \itemize{
    \item fixed bug of group comparison function
  }
}

\section{Version 1.8.2}{
  \itemize{
    \item fixed loading cluster from config file
  }
}

\section{Version 1.8.1}{
  \itemize{
    \item do not show pfam links for smart domains and vice versa
  }
}

\section{Version 1.7.8}{
  \itemize{
    \item only mainInput is required in config file
  }
}

\section{Version 1.7.7}{
  \itemize{
    \item use config file for input files, refspec selection and othes settings
  }
}

\section{Version 1.6.6}{
  \itemize{
    \item turn off auto sizing for large number of taxa or genes (>= 10.000)
  }
}

\section{Version 1.6.5}{
  \itemize{
    \item fixed filter for "species relation"
    \item added midpoint colors
  }
}

\section{Version 1.6.2}{
  \itemize{
    \item increase default font size for profile and domain plot
    \item make links to DBs: ncbi taxonomy, ncbi protein, uniprot, pfam, smart
  }
}

\section{Version 1.6.1}{
  \itemize{
    \item modified taxonomy ranks (ropensci#875)
    \item improved rank indexing function (ropensci#874)
    \item improved x-axis label (#116)
  }
}

\section{Version 1.4.12}{
  \itemize{
    \item domain plot works with group ID containing pipe
  }
}

\section{Version 1.4.11}{
  \itemize{
    \item added check for invalid ids in newTaxa file
  }
}

\section{Version 1.4.10}{
  \itemize{
    \item fixed bug comparing unk taxon IDs with max ncbi taxID
  }
}

\section{Version 1.4.9}{
  \itemize{
    \item fixed bug comparing taxon IDs with max ncbi taxID
  }
}

\section{Version 1.4.8}{
  \itemize{
    \item fixed bug filterProfileData (some filtered data could not be plotted)
  }
}

\section{Version 1.4.5}{
  \itemize{
    \item Add auto sizing feature
  }
}

\section{Version 1.4.4}{
  \itemize{
    \item Fixed orthoxml parsing function
  }
}

\section{Version 1.4.3}{
  \itemize{
    \item Plot background for NA orthologs
    \item Fixed highlight by categories
  }
}

\section{Version 1.4.2}{
  \itemize{
    \item Fixed error when parsing large amount of new taxa
    \item Improved the speed of parsing new taxa
    \item Apply filter to refspec orthologs when not working with strain/species
  }
}
 
\section{Version 1.4.1}{
  \itemize{
    \item Fixed bug data lost after filtering var2
    \item Add filter data to detailed plot
  }
}

\section{Version 1.3.11}{
  \itemize{
    \item Add filter number of co-orthologs to parseInfoProfile()
  }
}

\section{Version 1.3.10}{
  \itemize{
    \item Calculate percentage of present taxa after filtering of var1 and var2
    \item Fixed filter when working with high taxonomy ranks
  }
}

\section{Version 1.2.8}{
  \itemize{
    \item Added new NCBI taxonomy ranks (e.g. biotype, isolate, pathogroup, ...)
    \item Added function to reset taxonomy data
  }
}

\section{Version 1.2.7}{
  \itemize{
    \item Solved problem with new NCBI taxonomy rank "clade" by replace them with "norank"
  }
}

\section{Version 1.2.6}{
  \itemize{
    \item Fixed bug customized profile of subset of taxa not clickable
  }
}

\section{Version 1.2.4}{
  \itemize{
    \item Fixed bug checking invalid taxon IDs
  }
}

\section{Version 1.2.1}{
  \itemize{
    \item Fixed bug in rankIndexing and processNcbiTaxonomy
    \item Improved check for invalid input taxon IDs
  }
}

\section{Version 1.1.12}{
  \itemize{
    \item Added pseudo IDs for Holozoa and Holomycota clades and their non-NCBI sub-clades
  }
}

\section{Version 1.1.11}{
  \itemize{
    \item Increased the resolution of the phylostratigraphy function
  }
}

\section{Version 1.1.10}{
  \itemize{
    \item Improved the rank indexing function
  }
}

\section{Version 1.1.6}{
  \itemize{
    \item Removed grDevices, grid, stats, utils from Imports list in DESCRIPTION
  }
}

\section{Version 1.1.5}{
  \itemize{
    \item Added progress indicators for data loading and processing steps
    \item Removed shinycssloaders
  }
}

\section{Version 1.1.3}{
  \itemize{
    \item Fixed bug listing taxa by selecting a supertaxon from a higher level
    \item Added message for downloading preprocessed data
  }
}

\section{Version 1.1.2}{
  \itemize{
    \item Fixed bug sorting taxa when taxonomy rank is set as strain
  }
}

\section{Version 1.1.1}{
  \itemize{
    \item Commented examples of OMA parser functions
    \item Turned off tests for OMA parser functions
  }
}

\section{Version 0.99.31}{
  \itemize{
    \item Submitted to Bioconductor
  }
}