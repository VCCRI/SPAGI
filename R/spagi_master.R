#########################################################################################
#' @title SPAGI: Signalling Pathway Analysis for putative Gene regulatory network Identification
#'
#' @description SPAGI is an R package for active pathway identification for RNA-seq gene expression profiles. This package contains the neccessary R code to run SPAGI as described in "SPAGI: Identification of active signalling pathways using RNA-seq gene expression profiles". SPAGI is implemented within a helpful framework to identify active pathway paths from RNA-seq gene expression profiles.
#'
#'
#' @author Md Humayun Kabir <humayun.mkabir@gmail.com>
#'
#' @docType package
#' @name spagi-package
#' @aliases spagi SPAGI
#'
#' @examples
#' ## Do a sample analysis using human ocular lens epithelial cell (LEC) RNA-seq gene expression data.
#'
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#'
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format.
#' ROR1.processed.data<-preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway<-identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' ## Get active pathway ranking metric
#' ROR1.active.pathway.ranking.metric<-get_active_pathway_ranking_metric(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)
#' ## Display top n pathways
#' display_top_ranked_pathways(pathway.ranking.metric = ROR1.active.pathway.ranking.metric)
#'
NULL
########################################################################################






##################################################################################
# Copyright Victor Chang Cardiac Research Institute 2018
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################






######################################################################################
# Users should also check the names of functions loaded by this script
#
#
###############################################################
# Please note that it is necessary to pre-process the RNA-seq query data to calculate RPKM/FPKM/CPM of raw read count data and make log2 normalization before utilizing these data with the SPAGI package.
# Also note that the background pathway path and housekeeping gene data are in official gene symbol format. So please make your gene ids as official gene symbols before using with SPAGI package.
# The SPAGI package does not perform any normalizations. It assumes that all query data are in RPKM/FPKM/CPM and log2 normalized format.
# To utilize the SPAGI package, you have to provide an expression cut-off threshold and a high expression threshold (i.e., an expression value that is high enough) for your query data.
###############################################################
######################################################################################








#################################################################################################
# calculate mean of the cells based on replicates

#' @title compute_mean
#'
#' @description
#' This function calculates average expression values of multiple replicate samples
#'
#' @rdname compute_mean
#' @name compute_mean
#'
#' @details
#' This function calculates average expression values of multiple replicate samples using rowMeans. If the data has no replicates (a single vector or a matrix with only one column) the original values are returned.
#'
#' @return This function returns a vector with the mean expression values.
#'
#' @param data.matrix A matrix with gene expression values, rows denote genes and columns denote replicate samples. Alternatively a vector of gene expression values.
#'
#' @export
#'
#' @examples
#' dd<-matrix(sample(1:10, 30, replace=TRUE), 10, 3)
#' rownames(dd)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' compute_mean(dd)
#'

compute_mean <- function(data.matrix){
  if((is.vector(data.matrix)==TRUE) || (ncol(data.matrix)<2))
    data.matrix<-cbind(val1=data.matrix,vall2=data.matrix)
  return(rowMeans(data.matrix))
}
################################################################################################








####################################################################################################
# format matrix data of gene expression profile
# generate a matrix from user input (make averages of the replicates)

#' @title format_matrix_data
#'
#' @description
#' This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
#'
#' @rdname format_matrix_data
#' @name format_matrix_data
#'
#' @details
#' This function re-formats gene expression matrix data by averaging replicate samples. The cell type identifiers should be the column names of the expression matrix, otherwise the experiment.descriptor parameter is provided to specify cell-types of the samples (i.e. from a separate metadata file).
#'
#' @return This function returns a matrix with mean value of the replicates for each cell type.
#'
#' @param matrix.data A matrix with gene expression profiles where rows denote the genes and the columns denote the cells and/or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param experiment.descriptor A vector corresponding to the columns of matrix.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to the column names of the matrix.data.
#'
#' @export
#'
#' @examples
#' mm<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(mm)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(mm)<-c("cell2", "cell3", "cell1", "cell2", "cell1", "cell3", "cell3", "cell2", "cell3", "cell1")
#' format_matrix_data(mm)
#'

format_matrix_data<-function(matrix.data, experiment.descriptor = colnames(matrix.data)){
  colnames(matrix.data)<-experiment.descriptor

  #sorting the matrix data sets according to the column names
  ordered.column.names<-colnames(matrix.data)[order(colnames(matrix.data))]
  matrix.data.ordered<-matrix.data[,order(colnames(matrix.data))]
  colnames(matrix.data.ordered)<-ordered.column.names

  repCellTissue=table(colnames(matrix.data.ordered))
  nCellTissue<-length(repCellTissue)
  nameCellTissue<-names(repCellTissue)


  #processing the matrix data to make a new matrix with number of rows same as before, but
  #number of columns equal to the number of cells/tissues by calculating the average value of the replicates
  if(ncol(matrix.data.ordered)==nCellTissue){
    f.mat<-matrix.data.ordered  #here each cell/tissue has only 1 replicate, with existing cell names
  }
  else{                         #here each cell/tissue has one or more replicates
    f.mat<-matrix(0, nrow=nrow(matrix.data.ordered), ncol=nCellTissue)
    r<-repCellTissue
    s<-1
    e<-r[1]
    for(i in 1:nCellTissue){
      mean.data<-compute_mean(matrix.data.ordered[,s:e])
      f.mat[,i]<-mean.data
      if(i != nCellTissue){
        s<-s+r[i]
        e<-e+r[i+1]
      }
    }
    rownames(f.mat)<-rownames(matrix.data.ordered)
    colnames(f.mat)<-nameCellTissue
  }

  return(f.mat)
}
##############################################################################################








###########################################################################################################
# For average value calculation of the replicates and then get expressed genes based on expression cut-off;
# return a list of expressed genes for each cell type or tissue

#' @title preprocess_querydata
#'
#' @description
#' This function preprocesses the query data to calculate average value of the replicates and then get expressed genes based on expression cut-off.
#'
#' @rdname preprocess_querydata
#' @name preprocess_querydata
#'
#' @details
#' This function preprocesses the query data to calculate average value of the replicates and then get expressed genes based on expression cut-off.
#'
#' @return This function returns a list with specifically expressed genes for each cell type / tissue
#'
#' @param cell.tissue.data A matrix containing cell type / tissue RNA-seq gene expression data. It is assumed that all query data are in RPKM/FPKM/CPM and log normalized form. Also assume that gene ids are official gene symbols. For the matrix, rows denote the genes and the columns denote the cell types or tissues. Duplicate column names are expected in this case denoting replicate samples. All the replicate samples for a specific cell or tissue should have identical column names, otherwise the experiment.descriptor parameter should be used to identify replicate samples of a specific cell type or tissue.
#' @param exp.cutoff.th An expression cut-off threshold value for the query data.
#' @param species The species abbreviation of the query data (cell.tissue.data). Default is "hsapiens".
#' @param data.format Format of cell.tissue.data. Default is "matrix".
#' @param experiment.descriptor A vector corresponding to the matrix column names of cell.tissue.data, containing the cell type or tissues of each sample. The names should be identical for a specific cell or tissue.  Defaults to NULL.
#'
#' @export
#'
#' @examples
#' query.data<-matrix(sample(1:10, 100, replace=TRUE),10,10)
#' rownames(query.data)<-c("CRYAA", "CRYAB", "CRYBB3", "PAX6", "SOX2", "PROX1", "SIX3", "CRIM1", "CRYBB2", "BMP7")
#' colnames(query.data)<-c("cell1", "cell1", "cell1", "cell2", "cell2", "cell2", "cell3", "cell3", "cell3", "cell3")
#' preprocess_querydata(cell.tissue.data=query.data, exp.cutoff.th=5)
#'

preprocess_querydata<-function(cell.tissue.data, exp.cutoff.th, species = "hsapiens", data.format = "matrix", experiment.descriptor = NULL){
  if(data.format=="matrix"){
    ##to make all genes uppercase for mmusculus as PPI data are all uppercase
    if(species=="mmusculus"){
      rownames(cell.tissue.data)<-toupper(rownames(cell.tissue.data))
    }
    else if(species=="hsapiens"){
      cell.tissue.data<-cell.tissue.data
    }
    else{
      print("ERROR: could not support other species at this moment!!")
      return(NULL)
    }
    ##


    ##process the query data
    #calculate average value for each gene based on replicates for each cell type
    if(is.null(experiment.descriptor)) experiment.descriptor<-colnames(cell.tissue.data)
    cell.data.formatted<-format_matrix_data(cell.tissue.data, experiment.descriptor)

    #get expressed genes for each cell based on cut-off threshold
    expressed.cell.data<-list()
    cell.names<-colnames(cell.data.formatted)
    for(i in 1:ncol(cell.data.formatted)){
      each.cell.data<-cell.data.formatted[,i]
      names(each.cell.data)<-rownames(cell.data.formatted)
      #now get expressed genes based on cut-off value
      each.cell.exp.data<-each.cell.data[each.cell.data>=exp.cutoff.th]
      expressed.cell.data[[cell.names[i]]]<-each.cell.exp.data
    }

    #return the list with expressed genes data for each cell
    return(expressed.cell.data)
    ##
  }
  else{
    print("ERROR: could not support other data format at this moment!!")
    return(NULL)
  }
}
##########################################################################################








###########################################################################################
# identify active pathway path for RNA-seq gene expression profile
# generate a list of pathways where each sublist denote a path of the pathway

#' @title identify_active_pathway_path
#'
#' @description
#' This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
#'
#' @rdname identify_active_pathway_path
#' @name identify_active_pathway_path
#'
#' @details
#' This function identifies active pathway paths for RNA-seq gene expression profile. It utilizes background pathway path data to identify the active pathway paths.
#'
#' @return This function returns a list of pathways where each sublist denote a path of the pathway.
#'
#' @param pathway.path A list with pathway path data where each sublist denotes a path of the pathway. This is used as background data.
#' @param processed.query.data A list with expressed query data where each sublist corresponds for each cell/tissue type.
#'
#' @importFrom data.table chmatch
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#'
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format.
#' ROR1.processed.data<-preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway<-identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' head(ROR1.active.pathway$ROR1_LEC$FGFR1)
#'

identify_active_pathway_path<-function(pathway.path, processed.query.data){
  ##process separately for each cell or tissue type
  active_pathway_path<-lapply(processed.query.data, function(each.query.cell.exp.data){
    ##get the pathway paths in which all elements are expressed in query input data
    #each.cell.exp.gene.name<-names(each.query.cell.exp.data)
    each.cell.exp.gene.name<-sort(names(each.query.cell.exp.data))
    pathway.path.exist<-lapply(pathway.path, function(y){
      tmp.path.exist<-lapply(y, function(x){
        #for taking only the path where all molecules are expressed in gene expression data
        #if(all(unlist(x) %chin% each.cell.exp.gene.name == "TRUE")){
        if(!(anyNA(chmatch(x, each.cell.exp.gene.name)))){
          return(x)
        }
      })
      return(tmp.path.exist)
    })
    ##


    ##take only the existing pathway paths without null paths
    pathway.path.exist.clean<-lapply(pathway.path.exist, function(x){
      return(x[!(sapply(x, is.null))])
    })
    ##


    ##take only the pathways that have at least one complete path
    pathway.path.exist.clean.2<-list()
    for(i in 1:length(pathway.path.exist.clean)){
      if(length(pathway.path.exist.clean[[i]])!=0){
        pathway.path.exist.clean.2[[names(pathway.path.exist.clean)[i]]]<-pathway.path.exist.clean[[i]]
      }
    }
    ##


    ##return the active pathway path for each cell or tissue
    return(pathway.path.exist.clean.2)
    ##
  })


  ##Finally return the whole result for active pathway path of every cell/tissue
  return(active_pathway_path)
  ##
}
##########################################################################################









################################################################################################
# rank the active pathways based on their active (i.e., highly expressed) gene count proportion
# generate a list of pathways' ranking metric where each sublist corresponds for each cell/tissue type

#' @title get_active_pathway_ranking_metric
#'
#' @description
#' This function generates active pathway ranking metric for each cell/tissue type. It uses active pathway path and processed query data with a high expression threshold to generate the ranking metric.
#'
#' @rdname get_active_pathway_ranking_metric
#' @name get_active_pathway_ranking_metric
#'
#' @details
#' This function generates active pathway ranking metric for each cell/tissue type. It uses active pathway path and processed query data with a high expression threshold to generate the ranking metric.
#'
#' @return This function returns a list of pathways' ranking metric for each cell/tissue type.
#'
#' @param active.pathway.path A list of active pathway path data for each cell/tissue type.
#' @param processed.query.data A list with expressed query data where each sublist corresponds for each cell/tissue type.
#' @param high.exp.th A high expression threshold value for the processed query data. It is used to get active (i.e., highly expressed) molecule proportion for each pathway path.
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#'
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format.
#' ROR1.processed.data<-preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway<-identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' ## Get active pathway ranking metric
#' ROR1.active.pathway.ranking.metric<-get_active_pathway_ranking_metric(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)
#' head(sort(unlist(ROR1.active.pathway.ranking.metric$ROR1_LEC), decreasing=T))
#'

get_active_pathway_ranking_metric<-function(active.pathway.path, processed.query.data, high.exp.th){
  ##process separately each cell/tissue to get active pathway ranking metric
  active.pathway.ranking<-list()
  for(i in 1:length(active.pathway.path)){
    #get each cell/tissue active pathway paths
    each.cell.active.pathway.path<-active.pathway.path[[i]]

    #take the cell/tissue name from the pathway to get that cell/tissue processed.query.data and high.exp.th
    tmp.cell.name<-names(active.pathway.path[i])

    #take the respective cell/tissue processed data
    tmp.cell.processed.data<-processed.query.data[[tmp.cell.name]]

    ##get the pathway paths average active gene count proportion
    pathway.path.active.gene.count.proportion<-lapply(each.cell.active.pathway.path, function(y){
      individual.pathway.active.gene.count<-lapply(y, function(x){
        #get active (i.e., highly expressed) gene count proportion for each path and return
        tmp.path.gene.exp<-tmp.cell.processed.data[unlist(x)]
        tmp.path.active.genes<-names(tmp.path.gene.exp[tmp.path.gene.exp>=high.exp.th])
        return(length(tmp.path.active.genes)/length(unlist(x)))
      })
      #calcaulate average active gene count proportion for each pathway and return
      return(sum(unlist(individual.pathway.active.gene.count)) / length(y))
    })
    ##

    #assign average active gene count proportion for each cell/tissue
    active.pathway.ranking[[tmp.cell.name]]<-pathway.path.active.gene.count.proportion
  }
  #return average active gene count proportion for all cell/tissue
  return(active.pathway.ranking)
  ##
}
################################################################################################








################################################################################################
# display the sorted top n ranked pathways for each cell type or tissue

#' @title display_top_ranked_pathways
#'
#' @description
#' This function displays the sorted top n ranked pathways for each cell type or tissue in a barplot.
#'
#' @rdname display_top_ranked_pathways
#' @name display_top_ranked_pathways
#'
#' @details
#' This function displays the sorted top n ranked pathways for each cell type or tissue in a barplot.
#'
#' @return NULL
#'
#' @param pathway.ranking.metric The ranking metric result returned by 'get_active_pathway_ranking_metric' function.
#' @param top.n.pathway An integer number of how many top ranked pathways will be displayed. Default is 50.
#'
#' @export
#'
#' @examples
#' ## Here we will use "pathway.path" as background data from the SPAGI repository.
#' ## Also we will use "ROR1.data" as query RNA-seq gene expression data. This data is for ocular lens epithelial cell differentiated from human pluripotent stem cells.
#' ## These data sets are loaded automatically with the package.
#'
#' ## Pre-process the query data (ROR1.data), the data has already been made in CPM and log2 normalized format.
#' ROR1.processed.data<-preprocess_querydata(cell.tissue.data = ROR1.data, exp.cutoff.th = 1.8)
#' ## Identify active pathway paths of the processed query data
#' ROR1.active.pathway<-identify_active_pathway_path(pathway.path = pathway.path, processed.query.data = ROR1.processed.data)
#' ## Get active pathway ranking metric
#' ROR1.active.pathway.ranking.metric<-get_active_pathway_ranking_metric(active.pathway.path = ROR1.active.pathway, processed.query.data = ROR1.processed.data, high.exp.th = 7)
#' ## Display top n pathways
#' display_top_ranked_pathways(pathway.ranking.metric = ROR1.active.pathway.ranking.metric)
#'

display_top_ranked_pathways<-function(pathway.ranking.metric, top.n.pathway=50){
  #in each loop, the result of one query cell/tissue type is processed and showed in the bar plot
  for(i in 1:length(pathway.ranking.metric)){
    cell.tissue.names<-names(pathway.ranking.metric[i])
    #sort the pathway based on the ranking value
    individual.sorted.pathway<-sort(unlist(pathway.ranking.metric[[i]]), decreasing = T)

    #for setting the title
    if(nchar(cell.tissue.names)>40){
      title<-paste("Pathway ranking of:\n", cell.tissue.names, sep = " ")
    }
    else{
      title<-paste("Pathway ranking of", cell.tissue.names, sep = " ")
    }

    #taking reverse of the negative log10 p-values and displayed in the bar plot
    barplot(individual.sorted.pathway[1:top.n.pathway], main = title, xlab = "Ranking metric", las=2, horiz=TRUE, cex.names = 0.5)

    #for console message for each cell/tissue type ploting
    consol.msg<-paste(cell.tissue.names, "-- result plotting done!!", sep = " ")
    print(consol.msg)
  }
}
#####################################################################################################








##################################################################################################
#############Need two folder for downloading stringdb files for each species - stringdb_mouse, stringdb_human
#############It takes some time to download the data, and then can reuse the downloaded data

#' @title get_ppi_for_molecules
#'
#' @description
#' This function gets the PPI data from STRINGdb for the protein molecules provided.
#'
#' @rdname get_ppi_for_molecules
#' @name get_ppi_for_molecules
#'
#' @details
#' This function gets the PPI data from STRINGdb for the protein molecules provided.
#'
#' @return This function returns a data frame of the PPI for the molecules.
#'
#' @param RP.protein A vector containg the receptor (RP) proteins.
#' @param KN.protein A vector containg the kinase (KN) proteins.
#' @param TF.protein A vector containg the transcription factor (TF) proteins.
#' @param species The species name, either "hsapiens" or "mmusculus".
#' @param score An interger value for the STRINGdb PPI score threshold cutoff. Default is 700.
#'
#' @importFrom STRINGdb STRINGdb
#'
#' @export
#'
#' @examples
#' ## Need two folder at working directory for downloading stringdb files for each species - stringdb_mouse, stringdb_human.
#' ## It takes some time to download the data, and then can reuse the downloaded data.
#' ## Here we will use RP.protein, KN.protein, TF.protein protein parameters. These data are automatically loaded with the package. You can modify these parameters.
#' ## And we will use the species as "mmusculus".
#'
#' ## Get PPI data for the protein molecules of species "mmusculus".
#' mm.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="mmusculus")
#' head(mm.ppi)
#'

get_ppi_for_molecules<-function(RP.protein, KN.protein, TF.protein, species, score=700){
  ##get ppi interactions for molecules
  if(species=="mmusculus"){
    #initiate the connection, id  10090 for mouse
    string_db_mouse <- STRINGdb$new(version="10", species=10090, score_threshold=0, input_directory="stringdb_mouse" )
    #now combine all the protein
    all.protein<-unique(c(RP.protein, KN.protein, TF.protein))
    #make a data frame from all the protein
    all.protein.df<-data.frame("gene"=all.protein)
    #mapping gene names to string ids
    all.protein.mapped <- string_db_mouse$map(all.protein.df, "gene", takeFirst = T, removeUnmappedRows = TRUE)
    #get interactions information
    all.protein.mapped.interactions<-string_db_mouse$get_interactions(all.protein.mapped$STRING_id)
    #get only interactions and score
    all.protein.mapped.interactions.score<-all.protein.mapped.interactions[,c(1,2,16)]
  }
  else if(species=="hsapiens"){
    #initiate the connection, id  9606 for human
    string_db_human <- STRINGdb$new(version="10", species=9606, score_threshold=0, input_directory="stringdb_human" )
    #now combine all the protein and make uppercase
    all.protein<-toupper(unique(c(RP.protein, KN.protein, TF.protein)))
    #make a data frame from all the protein
    all.protein.df<-data.frame("gene"=all.protein)
    #mapping gene names to string ids
    all.protein.mapped <- string_db_human$map(all.protein.df, "gene", takeFirst = T, removeUnmappedRows = TRUE)
    #get interactions information
    all.protein.mapped.interactions<-string_db_human$get_interactions(all.protein.mapped$STRING_id)
    #get only interactions and score
    all.protein.mapped.interactions.score<-all.protein.mapped.interactions[,c(1,2,16)]
  }
  else{
    print("ERROR: Do not support other species at this moment.")
    return(NULL)
  }
  ##


  ##from STRING_id to gene name conversion
  all.factor.M<-all.protein.mapped.interactions.score
  all.factor.N<-all.protein.mapped
  all.factor.M[,1]<-all.factor.N[match(all.factor.M$from, all.factor.N$STRING_id),1]
  all.factor.M[,2]<-all.factor.N[match(all.factor.M$to, all.factor.N$STRING_id),1]
  all.factor.PPI<-all.factor.M
  ##


  ##get only the significant interactions, here by default combined score >= 700
  all.factor.PPI.significant<-all.factor.PPI[all.factor.PPI$combined_score>=score,]
  ##



  #########To get all interactions without considering the directions
  #########Here, we will take the highest score value for duplicates
  ##1st get the original interactions
  all.ppi.sig.1<-all.factor.PPI.significant
  rownames(all.ppi.sig.1)<-NULL
  ##


  #####
  ##combine the neighboring factors to treat as a single vector - original order
  comb.ppi.1<-list()
  for(i in 1:nrow(all.ppi.sig.1)){
    comb.ppi.1[[i]]<-paste(all.ppi.sig.1[i,1], all.ppi.sig.1[i,2], sep="*")
  }
  ##

  ##make the first df (original order) with the combined_score
  comb.ppi.1.df<-data.frame("interaction"=unlist(comb.ppi.1), "score"=all.ppi.sig.1$combined_score)
  ##
  #####


  #####
  ##combine the neighboring factors to treat as a single vector - reverse order
  comb.ppi.2<-list()
  for(j in 1:nrow(all.ppi.sig.1)){
    comb.ppi.2[[j]]<-paste(all.ppi.sig.1[j,2], all.ppi.sig.1[j,1], sep="*")
  }
  ##

  ##make the second df (reverse order) with the combined_score
  comb.ppi.2.df<-data.frame("interaction"=unlist(comb.ppi.2), "score"=all.ppi.sig.1$combined_score)
  ##
  #####


  ##Now add both the interactions' data frame - original order and reverse order
  comb.ppi.df<-rbind(comb.ppi.1.df, comb.ppi.2.df)
  ##

  ##order according to the score value - highest to lowest
  comb.ppi.df.ordered<-comb.ppi.df[order(comb.ppi.df$score, decreasing = T),]
  ##

  ##take PPIs with the highest score valued unique one from the duplicates
  comb.ppi.df.ordered.unique <- comb.ppi.df.ordered[!duplicated(comb.ppi.df.ordered$interaction),]
  rownames(comb.ppi.df.ordered.unique)<-NULL
  ##

  ##Finally return the combined PPI data frame with score value
  return(comb.ppi.df.ordered.unique)
  ##
  ##########
}
################################################################################################








#########################################################################################
###############Function for combining both mm.ppi and hs.ppi#############################
##########Also to get filtered PPI (RP-RP-KN-...-KN-TF) and their list of RPs and TFs####

#' @title combine_mm_hs_ppi
#'
#' @description
#' This function combines the PPI data for both "mmusculus" and "hsapiens" species created by the get_ppi_for_molecules function.
#'
#' @rdname combine_mm_hs_ppi
#' @name combine_mm_hs_ppi
#'
#' @details
#' This function combines the PPI data for both "mmusculus" and "hsapiens" species created by the get_ppi_for_molecules function.
#'
#' @return This function returns a list consisting of the combined filtered PPI data, the RP proteins and the TF proteins of the combined filtered PPI data to generate the pathway path.
#'
#' @param mm.ppi The PPI data for "mmusculus" species generated by the function get_ppi_for_molecules.
#' @param hs.ppi The PPI data for "hsapiens" species generated by the function get_ppi_for_molecules.
#' @param RP.protein A vector containg the same receptor (RP) proteins that are used in the function get_ppi_for_molecules.
#' @param KN.protein A vector containg the same kinase (KN) proteins that are used in the function get_ppi_for_molecules.
#' @param TF.protein A vector containg the same transcription factor (TF) proteins that are used in the function get_ppi_for_molecules.
#'
#' @export
#'
#' @examples
#' ## Need two folder at working directory for downloading stringdb files for each species - stringdb_mouse, stringdb_human.
#' ## It takes some time to download the data, and then can reuse the downloaded data.
#' ## Here we will use RP.protein, KN.protein, TF.protein protein parameters. These data are automatically loaded with the package. You can modify these parameters.
#' ## We will generate PPI data for two species - "mmusculus" and "hsapiens" by calling the function get_ppi_for_molecules two times.
#' ## Then we will combine these two PPI data sets by using the combine_mm_hs_ppi function that will be used later on to generate the pathway path data.
#'
#' ## Get PPI data for the protein molecules of species "mmusculus".
#' mm.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="mmusculus")
#' ## Get PPI data for the protein molecules of species "hsapiens".
#' hs.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="hsapiens")
#' ## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
#' comb.ppi.result<-combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
#' head(summary(comb.ppi.result))
#'

combine_mm_hs_ppi<-function(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein){
  #####combine, order and take the unique PPIs with highest score
  #Combine the both PPI data
  comb.ppi<-rbind(mm.ppi, hs.ppi)
  #order according to the score value - highest to lowest
  comb.ppi.ordered<-comb.ppi[order(comb.ppi$score, decreasing = T),]
  #take PPIs with the highest score valued unique one from the duplicates
  comb.ppi.ordered.unique <- comb.ppi.ordered[!duplicated(comb.ppi.ordered$interaction), ]
  #####



  #####
  #now separating the links using a list that contains all the links as vectors
  comb.ppi.interaction.split<-lapply(as.vector(comb.ppi.ordered.unique$interaction), function(x) {return(unlist(strsplit(x, split = "[*]")))})
  #making data frame from the unique split lists
  comb.ppi.interaction.split.df<-as.data.frame(do.call(rbind, lapply(comb.ppi.interaction.split, rbind)))
  #set the column names of the data frame
  colnames(comb.ppi.interaction.split.df)<-c("from", "to")
  #now add the score value as a 3rd column
  all.factor.PPI.significant<-data.frame(comb.ppi.interaction.split.df, "score" = as.vector(comb.ppi.ordered.unique$score))
  #####



  #####First make all protein symbols as uppercase
  RP.protein<-toupper(RP.protein)
  KN.protein<-toupper(KN.protein)
  TF.protein<-toupper(TF.protein)
  #####



  #####To get only the significant links exist from RP - KN - TF
  #####FOr RP - RP, we have allowed maximum of 2 layers according to our design,
  #####If you need different design you should change in this section according to your design.
  ##get interactions from RP to KN
  RP.to.KN.significant.ppi<-all.factor.PPI.significant[((all.factor.PPI.significant$from %in% RP.protein) &
                                                          (all.factor.PPI.significant$to %in% KN.protein)),]

  ##get interactions from KN to KN - for all KNs
  KN.to.KN.significant.ppi<-all.factor.PPI.significant[((all.factor.PPI.significant$from %in% KN.protein) &
                                                          (all.factor.PPI.significant$to %in% KN.protein)),]

  ##get interactions from KN to TF - for all KNs
  KN.to.TF.significant.ppi<-all.factor.PPI.significant[((all.factor.PPI.significant$from %in% KN.protein) &
                                                          (all.factor.PPI.significant$to %in% TF.protein)),]

  ##get the RPs that have no direct interaction with the KNs
  RP.not.connected.with.KN<-setdiff(RP.protein, unique(RP.to.KN.significant.ppi$from))

  ##get the ppi from 'RP.not.connected.with.KN' to 'unique(RP.to.KN.significant.ppi$from)'
  #this will give us interaction for 2 RP layers
  #get interactions from from RP not connected with KN to RP connected with KN
  #these combined RPs will act as source to finding the paths
  RP.to.RP.significant.ppi<-all.factor.PPI.significant[((all.factor.PPI.significant$from %in% RP.not.connected.with.KN) &
                                                          (all.factor.PPI.significant$to %in% unique(RP.to.KN.significant.ppi$from))),]

  ##And finally combine all the interactions from RP-RP-KN-...-KN-TF
  all.significant.filtered.ppi<-rbind(RP.to.RP.significant.ppi, RP.to.KN.significant.ppi,
                                      KN.to.KN.significant.ppi, KN.to.TF.significant.ppi)
  rownames(all.significant.filtered.ppi)<-NULL
  #####



  #####
  ##Now get the RP and TF of the interactions
  RPs<-unique(c(unique(as.vector(RP.to.RP.significant.ppi$from)),
                unique(as.vector(RP.to.KN.significant.ppi$from))))

  TFs<-unique(as.vector(KN.to.TF.significant.ppi$to))
  #####



  #####
  ##Finally make a list of all.significant.filtered.ppi, RPs and TFs and then return
  comb.ppi.result<-list()
  comb.ppi.result[["PPI"]]<-all.significant.filtered.ppi
  comb.ppi.result[["RPs"]]<-RPs
  comb.ppi.result[["TFs"]]<-TFs
  return(comb.ppi.result)
  ##
  #####
}
###########################################################################################








##############################################################################################
#####Here may be you will get some warnings, but these are only for not reachable paths for source to destination

#' @title generate_pathway_path
#'
#' @description
#' This function generates the background pathway path data using the list object returned by the function combine_mm_hs_ppi.
#'
#' @rdname generate_pathway_path
#' @name generate_pathway_path
#'
#' @details
#' This function generates the background pathway path data using the list object returned by the function combine_mm_hs_ppi.
#'
#' @return This function returns a list consisting of the pathway path data that will be used as background data for the SPAGI package.
#'
#' @param ppi.result The list object returned by the function combine_mm_hs_ppi.
#' @param housekeeping.gene A vector consisting of housekeeping genes. This data is loaded automatically with the package. This data was generated using the gene expression profiles of different cell types and/or tissues from the ENCODE human and mouse project.
#' @param max.path.length An integer number indicating the maximum path length. Default is 7.
#'
#' @importFrom igraph graph.data.frame shortest_paths
#'
#' @export
#'
#' @examples
#' ## Need two folder at working directory for downloading stringdb files for each species - stringdb_mouse, stringdb_human.
#' ## It takes some time to download the data, and then can reuse the downloaded data.
#' ## Here we will use RP.protein, KN.protein, TF.protein protein parameters. These data are automatically loaded with the package. You can modify these parameters.
#' ## We will generate PPI data for two species - "mmusculus" and "hsapiens" by calling the function get_ppi_for_molecules two times.
#' ## Then we will combine these two PPI data sets by using the combine_mm_hs_ppi function that will be used to generate the pathway path data.
#' ## Finally we will generate the pathway path data using the combined PPI data
#'
#' ## Get PPI data for the protein molecules of species "mmusculus".
#' mm.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="mmusculus")
#' ## Get PPI data for the protein molecules of species "hsapiens".
#' hs.ppi<-get_ppi_for_molecules(RP.protein, KN.protein, TF.protein, species="hsapiens")
#' ## Now combine and get the filtered PPI and the RP and TF proteins of the combined filtered PPI
#' comb.ppi.result<-combine_mm_hs_ppi(mm.ppi, hs.ppi, RP.protein, KN.protein, TF.protein)
#' ##Generate the pathway path data using the comb.ppi.result and housekeeping.gene data sets
#' pathway.path<-generate_pathway_path(ppi.result=comb.ppi.result, housekeeping.gene)
#' head(summary(pathway.path))
#'

generate_pathway_path<-function(ppi.result, housekeeping.gene, max.path.length=7){
  #####preprocess the ppi.result data
  ##Assign the result data to the objects
  all.significant.filtered.ppi<-ppi.result$PPI
  RPs<-ppi.result$RPs
  TFs<-ppi.result$TFs
  ##

  ##NOTE
  ##First, it is good to know that when looking up paths, igraph understands weights as costs,
  ##i.e. on edges with higher weight it costs more to travel,
  ##so it will consider shorter the paths with lower sum weight.
  ##It is easy to turn this into the opposite, here we will do the score (range 0 to 999) as weight by 1000-score.

  ##First calculate the weight from the score as (1000-score)
  edge.weight<- 1000 - all.significant.filtered.ppi$score
  ##

  ##Now, add the weight column to the data frame
  all.edges<-data.frame(all.significant.filtered.ppi[,1:2], "weight"=edge.weight)
  #make the rownames null
  rownames(all.edges)<-NULL
  ##
  #####



  #####Now create a graph and generate the pathway paths
  ##make the graph data frame from the "all.edges"
  g1 <- graph.data.frame(d = all.edges, directed = TRUE)
  ##

  ##Find all shortest paths for all source nodes ("RPs") to destination ("TFs")
  ##Here by default uses the Dijkstra's algorithm for weighted directed graph
  ll.all.path<-list()
  for(i in 1:length(RPs)){
    ll.all.path[[RPs[i]]]<-shortest_paths(g1, from=RPs[i], to = TFs, mode = "out")
    print(i)
  }
  ##

  ##Finding all the complete (RP-KN-...-KN-TF) paths
  ##Here we considered for at least 3 layers and maximum 7 layers by default
  ll.all.path.complete<-list()
  for(i in 1:length(ll.all.path)){
    #this tmp variable is used for all filtered paths for a pathway
    tmp.individual.pathway.paths.clean.ind<-NULL
    #this loop to check each individual path for a pathway
    for(j in 1:length(ll.all.path[[i]]$vpath)){
      if((length(unlist(ll.all.path[[i]]$vpath[j])) >= 3) & (length(unlist(ll.all.path[[i]]$vpath[j])) <= max.path.length))
        tmp.individual.pathway.paths.clean.ind<-c(tmp.individual.pathway.paths.clean.ind, j)
    }
    #this combines all the filtered paths for a pathway to the respective pathway
    ll.all.path.complete[[names(ll.all.path)[i]]]<-ll.all.path[[i]]$vpath[tmp.individual.pathway.paths.clean.ind]
  }
  ##

  ##take only the pathways that have at least one complete path
  ll.all.path.complete.exist<-list()
  for(i in 1:length(ll.all.path.complete)){
    if(length(ll.all.path.complete[[i]])!=0){
      ll.all.path.complete.exist[[names(ll.all.path.complete)[i]]]<-ll.all.path.complete[[i]]
    }
  }
  ##

  ##get only the pathway path elements, i.e., the names
  pathway.path.all<-list()
  for(i in 1:length(ll.all.path.complete.exist)){
    tmp.pathway.path<-lapply(ll.all.path.complete.exist[[i]], function(x){return(names(x))})
    pathway.path.all[[names(ll.all.path.complete.exist)[i]]]<-tmp.pathway.path
  }
  ##
  #####



  #####This section is for removing the paths where all elements are housekeeping genes
  ##get the pathway paths in which all elements are not hk genes
  pathway.path.specific<-list()
  for(i in 1:length(pathway.path.all)){
    tmp.path.spec<-lapply(pathway.path.all[[i]], function(x){
      if(!(all(x %in% housekeeping.gene == "TRUE")))
        return(x)
    })
    pathway.path.specific[[names(pathway.path.all)[i]]]<-tmp.path.spec
  }
  ##

  ##take only the existing pathway paths without null paths
  pathway.path.specific.clean<-lapply(pathway.path.specific, function(x){
    return(x[!(sapply(x, is.null))])
  })
  ##

  ##take only the pathways that have at least one path
  pathway.path.specific.clean.2<-list()
  for(i in 1:length(pathway.path.specific.clean)){
    if(length(pathway.path.specific.clean[[i]])!=0){
      pathway.path.specific.clean.2[[names(pathway.path.specific.clean)[i]]]<-pathway.path.specific.clean[[i]]
    }
  }
  ##

  ##Finally return the pathway.path.specific.clean.2 data
  ##This data will be used as background pathway path data
  return(pathway.path.specific.clean.2)
  ##
  #####
}
################################################################################################





