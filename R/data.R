#' Pathway path data
#'
#' Pathway path used as background data, generated using STRINGdb PPI.
#'
#' @docType data
#'
#' @usage data(pathway.path)
#'
#' @format An object of list containing pathway paths that is used as background data.
#'
#' @keywords datasets
#'
#' @examples
#' head(summary(pathway.path)
"pathway.path"


#' ROR1 data
#'
#' Data from the ROR1+ cells differentiated from human pluripotent stem cells,
#' contains transcriptomic profiles data of lens epithelial cell (LEC).
#' The cell has 2 biological replicates.
#'
#' @docType data
#'
#' @usage data(ROR1.data)
#'
#' @format An object of matrix containing Count Per Million (CPM) and log2 normalized gene expression value of LEC.
#'
#' @keywords datasets
#'
#' @references Murphy P et al, Development 2018 Jan 9; 145(1)
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/29217756}{PubMed})
#'
#' @examples
#' head(ROR1.data)
"ROR1.data"


#' Housekeeping gene data
#'
#' Housekeeping gene used as background data, generated using gene expression profiles of different cells/tissues from ENCODE project.
#'
#' @docType data
#'
#' @usage data(housekeeping.gene)
#'
#' @format An object of vector containing housekeeping genes that is used as background data.
#'
#' @keywords datasets
#'
#' @examples
#' head(housekeeping.gene)
"housekeeping.gene"


#' Receptor protein data
#'
#' Receptor protein, used as source components of pathway paths.
#'
#' @docType data
#'
#' @usage data(RP.protein)
#'
#' @format An object of vector containing receptor proteins.
#'
#' @keywords datasets
#'
#' @examples
#' head(RP.protein)
"RP.protein"


#' Kinase protein data
#'
#' Kinase protein, used as middle components of pathway paths.
#'
#' @docType data
#'
#' @usage data(KN.protein)
#'
#' @format An object of vector containing kinase proteins.
#'
#' @keywords datasets
#'
#' @examples
#' head(KN.protein)
"KN.protein"


#' Transcription factor protein data
#'
#' Transcription factor protein, used as target components of pathway paths.
#'
#' @docType data
#'
#' @usage data(TF.protein)
#'
#' @format An object of vector containing Transcription factor proteins.
#'
#' @keywords datasets
#'
#' @examples
#' head(TF.protein)
"TF.protein"
