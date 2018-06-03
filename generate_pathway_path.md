Pathway path generation
================
Md Humayun Kabir

<!-- README.md is generated from README.Rmd -->
Instruction to generate pathway path background data for SPAGI
==============================================================

Introduction
------------

This file contains the necessary instruction and R code of how to utilize the *SPAGI* package functions to generate the background pathway path data.

To generate the background pathway path data it is needed to create two folder (stringdb\_mouse and stringdb\_human) at working directory for downloading stringdb files for each species "mmusculus" and "hsapiens". It takes some time to download the data, and then can reuse the downloaded data. Here we will use RP.protein, KN.protein, TF.protein protein as molecules for the pathway path data. These data are collected from public data sets. These data sets are automatically loaded with the package. We will generate PPI data for two species - "mmusculus" and "hsapiens" by calling the function get\_ppi\_for\_molecules two times for the same set of protein molecules. Then we will combine these two PPI data sets by using the combine\_mm\_hs\_ppi function. Finally we will generate the pathway path data using the combined PPI data. In the final step we will use a list of housekeeping genes to filter out the paths where all molecules are housekeeping genes. The list of housekeeping genes are generated based on gene expression profiles of a large number of cell types or tissues obtained from the ENCODE project.

To generate the background pathway path data the *SPAGI* package relies on two other packages *STRINGdb* for the protein-protein interaction of the molecules and *igraph* for implementation of graph algorithm.

Installation
------------

*SPAGI* depends on one package *data.table*. To generate the background pathway path data it relies on packages *STRINGdb* and *igraph*. Make sure you have installed all the packages. The commands to do so are as follows:

Open an R session and do the following:

``` r
source("http://bioconductor.org/biocLite.R")
biocLite("STRINGdb")

install.packages('data.table')
install.packages('igraph')
```

Make sure you have *devtools* installed, then install the *SPAGI* package from github repository:

``` r
install.packages('devtools')
devtools::install_github('VCCRI/SPAGI')
```

Finally load the packages with:

``` r
library(STRINGdb)
library(igraph)
library(spagi)
```

Package installation and loading is done!

Now, you can run the following code to generate the pathway path data:

``` r
## Get PPI data for the protein molecules of species "mmusculus".
mm.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="mmusculus")
## Get PPI data for the protein molecules of species "hsapiens".
hs.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="hsapiens")
## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
comb.ppi.result<-combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
##Generate the pathway path data using the comb.ppi.result and housekeeping.gene data sets
pathway.path<-generate_pathway_path(ppi.result=comb.ppi.result, housekeeping.gene)
head(summary(pathway.path))
```

Pathway path data generation is complete.
-----------------------------------------
