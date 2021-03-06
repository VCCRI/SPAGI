% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spagi_master.R
\name{combine_mm_hs_ppi}
\alias{combine_mm_hs_ppi}
\title{combine_mm_hs_ppi}
\usage{
combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
}
\arguments{
\item{mm.ppi}{The PPI data for "mmusculus" species generated by the function get_ppi_for_molecules.}

\item{hs.ppi}{The PPI data for "hsapiens" species generated by the function get_ppi_for_molecules.}

\item{RP.protein}{A vector containg the same receptor (RP) proteins that are used in the function get_ppi_for_molecules.}

\item{KN.protein}{A vector containg the same kinase (KN) proteins that are used in the function get_ppi_for_molecules.}

\item{TF.protein}{A vector containg the same transcription factor (TF) proteins that are used in the function get_ppi_for_molecules.}
}
\value{
This function returns a list consisting of the combined filtered PPI data, the RP proteins and the TF proteins of the combined filtered PPI data to generate the pathway path.
}
\description{
This function combines the PPI data for both "mmusculus" and "hsapiens" species created by the get_ppi_for_molecules function.
}
\details{
This function combines the PPI data for both "mmusculus" and "hsapiens" species created by the get_ppi_for_molecules function.
}
\examples{
## Need two folder at working directory for downloading stringdb files for each species - stringdb_mouse, stringdb_human.
## It takes some time to download the data, and then can reuse the downloaded data.
## Here we will use RP.protein, KN.protein, TF.protein protein parameters. These data are automatically loaded with the package. You can modify these parameters.
## We will generate PPI data for two species - "mmusculus" and "hsapiens" by calling the function get_ppi_for_molecules two times.
## Then we will combine these two PPI data sets by using the combine_mm_hs_ppi function that will be used later on to generate the pathway path data.

## Get PPI data for the protein molecules of species "mmusculus".
mm.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="mmusculus")
## Get PPI data for the protein molecules of species "hsapiens".
hs.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="hsapiens")
## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
comb.ppi.result<-combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
head(summary(comb.ppi.result))

}
