#' Convert Data to Appropriate pmartRseq Class
#'
#' Converts a list object or several data.frames of rRNA (16S/ITS/18S),
#' metatranscript, or metagenomic data to an object of the class 'seqData'.
#' Objects of the class 'seqData' are lists with two obligatory components
#' \code{e_data} and \code{f_data}. An optional list component \code{e_meta}
#' is used if analysis or visualization at other levels (e.g. taxonomy) is
#' also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of features observed and \eqn{n} is the
#'   number of samples (an additional feature identifier/name column should also
#'   be present anywhere in the data.frame). Each row corresponds to data for
#'   each feature.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample.
#' @param e_meta an optional data.frame with \eqn{p} rows. Each row corresponds
#'   to a feature with one column giving identifiers (must be named the same as
#'   the column in \code{e_data}) and other columns giving meta information
#'   (e.g. mappings of OTU identification to taxonomy).
#' @param e_tree an optional NEXUS or Newick formatted phylogenetic tree file,
#'   imported using ape::read.tree(tree_path). The OTU labels in the tree file
#'   should match the OTU identifiers in the preceeding data fields.
#' @param e_seq an optional fasta formatted representation of biological
#'   sequences imported using Biostrings::readDNAStringSet(fasta_path, ...).
#'   Each OTU in the fasta maps to at least one sequence in the preceeding data
#'   fields.
#' @param data_type character string specifying if this is 'rRNA' (for
#'   16S/ITS/18S), 'metagenomic', or 'metatranscriptomic' data.
#' @param edata_cname character string specifying the name of the column
#'   containing the identifiers in \code{e_data} and \code{e_meta} (if
#'   applicable).
#' @param taxa_cname optional character string specifying the name of the column
#'   containing the taxonomy in \code{e_meta} (if applicable).
#'   Defaults to NULL. If \code{e_meta} is NULL, then specify \code{taxa_cname}
#'   as NULL.
#' @param ec_cname optional character string specifying the name of the column
#'   containing the EC numbers in \code{e_meta} (if applicable).
#'   Defaults to NULL. If \code{e_meta} is NULL, then specify \code{ec_cname}
#'   as NULL.
#' @param gene_cname optional character string specifying the name of the column
#'   containing the gene names in \code{e_meta} (if applicable).
#'   Defaults to NULL. If \code{e_meta} is NULL, then specify \code{gene_cname}
#'   as NULL.
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param ... further arguments
#'
#' @details Objects of class 'seqData' contain some attributes that are
#'   referenced by downstream functions. These attributes can be changed from
#'   their default value by manual specification. A list of these attributes as
#'   well as their default values are as follows: \tabular{ll}{ data_scale \tab
#'   Scale of the data provided in \code{e_data}. Acceptable values are 'log2',
#'   'log10', 'log', 'count', and 'abundance', which indicate data is log base
#'   2, base 10, natural log transformed, raw count data, and raw abundance,
#'   respectively. Default values is 'count'. \cr \tab \cr data_norm \tab A
#'   logical argument, specifying whether the data has been normalized or not.
#'   Default value is 'FALSE'. \cr \tab \cr norm_method \tab Null if data_norm
#'   is FALSE. If data_norm is TRUE, character string defining which
#'   normalization method was used. Default value is 'NULL'. \cr \tab \cr
#'   location_param \tab NULL if there are no location parameters from
#'   normalization, otherwise a vector detailing the normalization location
#'   parameters for each sample. \cr \tab \cr scale_param \tab NULL if there are
#'   no scale parameters from normalization, otherwise a vector detailing the
#'   normalization scale parameters for each sample. \cr \tab \cr seq_type
#'   \tab Character string describing the type of sequencer (e.g. 'HiSeq').
#'   Default value is 'NULL'. \cr \tab \cr db \tab Character
#'   string describing which database was used to process the data (e.g.
#'   "TIGR"). Default value is 'NULL'. \cr \tab \cr db_version \tab Character
#'   string describing which version of the database was used. Default value is
#'   'NULL'. If db is NULL, then db_version will default to a NULL value. }
#'
#'   Computed values included in the \code{data_info} attribute are as follows:
#'   \tabular{ll}{ num_edata \tab The number of unique \code{edata_cname}
#'   entries.\cr \tab \cr num_na \tab The number of NA observations in the
#'   dataset.\cr \tab \cr frac_na \tab The prportion of \code{e_data} values
#'   that are NA. \cr \tab \cr num_zero \tab The number of observations that
#'   equal 0 in the dataset. \cr \tab \cr frac_zero \tab The proportion of
#'   \code{e_data} values that are 0. \cr \tab \cr num_taxa \tab The number of
#'   unique \code{taxa_cname} entries. \cr \tab \cr num_ec \tab The number of
#'   unique \code{ec_cname} entries. \cr \tab \cr num_gene \tab The number of
#'   unique \code{gene_cname} entries. \cr \tab \cr num_samps \tab The number
#'   of samples that make up the columns of \code{e_data}.\cr \tab \cr meta_info
#'   \tab A logical argument, specifying whether \code{e_meta} is provided.\cr
#'   \tab \cr }
#'
#' @author Allison Thompson, Lisa Bramer
#'
#'
#' @export

as.seqData <- function(e_data, f_data, e_meta=NULL, edata_cname, fdata_cname, data_type, taxa_cname=NULL, ...){
  .as.seqData(e_data, f_data, e_meta, e_tree, e_seq, edata_cname, fdata_cname, data_type, taxa_cname, ...)
}

## seqData ##
.as.seqData <- function(e_data, f_data, e_meta=NULL, e_tree, e_seq,
                         edata_cname, fdata_cname, data_type,
                         taxa_cname, ec_cname=NULL, gene_cname=NULL,
                         data_scale="count", data_norm=FALSE, norm_method=NULL,
                         location_param=NULL, scale_param=NULL,
                         seq_type=NULL, db=NULL, db_version=NULL, ...){
  # initial checks #
  #browser()

  # names checks #
  colnames(e_data) <- gsub("[^A-Za-z0-9_]","\\.",colnames(e_data))
  colnames(f_data) <- gsub("[^A-Za-z0-9_]","\\.",colnames(f_data))
  colnames(e_meta) <- gsub("[^A-Za-z0-9_]","\\.",colnames(e_meta))

  colnames(f_data) <- sapply(colnames(f_data), function(x) ifelse(!is.na(as.numeric(substr(x,1,1))), paste("X",x,sep=""),x))
  colnames(e_data) <- sapply(colnames(e_data), function(x) ifelse(!is.na(as.numeric(substr(x,1,1))), paste("X",x,sep=""),x))
  colnames(e_meta) <- sapply(colnames(e_meta), function(x) ifelse(!is.na(as.numeric(substr(x,1,1))), paste("X",x,sep=""),x))

  if(!is.na(as.numeric(substr(fdata_cname, 1, 1)))){
    fdata_cname <- paste("X",fdata_cname,sep="")
  }

  fdata_cname <- gsub("[^A-Za-z0-9_]","\\.",fdata_cname)
  edata_cname <- gsub("[^A-Za-z0-9_]","\\.",edata_cname)
  taxa_cname <- if(!is.null(taxa_cname)) gsub("[^A-Za-z0-9_]","\\.",taxa_cname)
  ec_cname <- if(!is.null(ec_cname)) gsub("[^A-Za-z0-9_]","\\.",ec_cname)
  gene_cname <- if(!is.null(gene_cname)) gsub("[^A-Za-z0-9_]","\\.",gene_cname)

  f_data[,fdata_cname] <- gsub("[^A-Za-z0-9_]","\\.",f_data[,fdata_cname])
  f_data[,fdata_cname] <- sapply(f_data[,fdata_cname], function(x) ifelse(!is.na(as.numeric(substr(x,1,1))), paste("X",x,sep=""),x))

  # check that the OTU column exists in e_data and e_meta (if applicable) #
  if (!(edata_cname %in% names(e_data))) stop(paste("OTU column ", edata_cname," not found in e_data. See details of as.seqData for specifying column names.", sep = ""))

  # check that data_type is correct type #
  if(!(tolower(data_type) %in% c("rrna","metag","metagenome","metagenomic","metat","metatranscript","metatranscriptome","metatranscriptomic"))){
    stop(paste("Data type ",data_type," is not supported - must be one of 'rRNA', 'metagenomic', or 'metatranscriptomic'.",sep=""))
  }

  # standardize data_type for metagenomic data #
  if(tolower(data_type) %in% c("metag","metagenome","metagenomic")){
    data_type = "metagenomic"
  }

  # standardize data_type for metatranscriptomic data #
  if(tolower(data_type) %in% c("metat","metatranscript","metatranscriptome","metatranscriptomic")){
    data_type = "metatranscriptomic"
  }

  # check that edata_cname is in e_meta, if e_meta is provided #
  if (!is.null(e_meta)) {
    if (!(edata_cname %in% names(e_meta))) stop(paste("OTU column ", edata_cname," not found in e_meta. Column names for OTUs must match for e_data and e_meta. See details of as.seqData for specifying column names.", sep = ""))
  }

  # if e_meta is NULL set emeta_cnames to NULL #
  if (is.null(e_meta)){
    taxa_cname = NULL
    ec_cname = NULL
    gene_cname = NULL
  }

  # if e_meta is not NULL check that the emeta_cnames (taxa_cname, ec_cname, gene_cname) columns are found, if given #
  if (!is.null(e_meta)){
    # check taxa_cname
    if (!is.null(taxa_cname)){
      if (!(taxa_cname %in% names(e_meta))) stop(paste("Taxonomy column ", taxa_cname, " not found in e_meta. See details of as.seqData for specifying column names.", sep = "") )
    }
    # check ec_cname
    if (!is.null(ec_cname)){
      if (!(ec_cname %in% names(e_meta))) stop(paste("Enzyme commission number column ", ec_cname, " not found in e_meta. See details of as.seqData for specifying column names.", sep = "") )
    }
    # check gene_cname
    if (!is.null(gene_cname)){
      if (!(gene_cname %in% names(e_meta))) stop(paste("Gene column ", gene_cname, " not found in e_meta. See details of as.seqData for specifying column names.", sep = "") )
    }
  }

  # check that the Sample column name is in f_data column names #
  if (!(fdata_cname %in% names(f_data))) stop(paste("Sample column ", fdata_cname, " not found in f_data. See details of as.seqData for specifying column names.", sep = ""))

  # check that all samples in e_data are present in f_data #
  edat_sampid = which(names(e_data) == edata_cname)
  samps.miss = sum(!(names(e_data[,-edat_sampid]) %in% f_data[,fdata_cname]))
  if( samps.miss > 0) stop(paste( samps.miss, " samples from e_data not found in f_data", sep = ""))

  # check for any extra samples in f_data than in e_data - necessary to remove before group_designation function #
  if(any(!(f_data[,fdata_cname] %in% names(e_data)))){
    f_data <- f_data[-which(!(f_data[,fdata_cname] %in% names(e_data))),]
  }

  # check that f_data has at least 2 columns #
  if(ncol(f_data) < 2) stop("f_data must contain at least 2 columns")

  # if e_meta is provided, check that all OTUs in e_data occur in e_meta #
  if(!is.null(e_meta)){
    if(sum(!(e_data[,edata_cname] %in% e_meta[,edata_cname])) > 0 ) stop("Not all OTUs in e_data are present in e_meta")
  }

  # if e_meta is provided, remove any extra features that were provided #
  if(!is.null(e_meta)){
    if(any(!(e_meta[,edata_cname] %in% e_data[,edata_cname]))){
      e_meta <- e_meta[-which(!(e_meta[,edata_cname] %in% e_data[,edata_cname])),]
    }
  }

  # check that data_scale is one of the acceptable options #
  if(!(data_scale %in% c('log2', 'log10', 'log', 'count', 'abundance'))) stop(paste(data_scale, " is not a valid option for 'data_scale'. See details of as.seqData for specifics.", sep=""))

  # store results #
  seqObj = list(e_data = e_data, f_data = f_data, e_meta = e_meta)
  # #----check for unhallowed characters in column names----#
  # names(seqObj$f_data) <- make.names(names(seqObj$f_data))
  # names(seqObj$e_data) <- make.names(names(seqObj$e_data))

  # set column name attributes #
  attr(seqObj, "cnames") = list(edata_cname = edata_cname, fdata_cname = fdata_cname, taxa_cname = taxa_cname, ec_cname = ec_cname, gene_cname = gene_cname)

  # count missing values in e_data #
  #num_miss_obs = sum(is.na(e_data[,-which(names(e_data)==edata_cname)])) + length(which(e_data[,-which(names(e_data)==edata_cname)] == 0))
  #frac_missing = mean(is.na(e_data[,-which(names(e_data)==edata_cname)]))
  num_na = sum(is.na(e_data[,-which(names(e_data)==edata_cname)]))
  num_zero = length(which(e_data[,-which(names(e_data)==edata_cname)] == 0))
  frac_na = mean(is.na(e_data[,-which(names(e_data)==edata_cname)]))
  frac_zero = length(which(e_data[,-which(names(e_data)==edata_cname)]==0)) / length(which(e_data[,-which(names(e_data)==edata_cname)]>=0))

  # number of unique OTUs #
  num_edata = length(unique(e_data[, edata_cname]))

  # number of samples #
  num_samps = ncol(e_data) - 1

  if(!is.null(e_meta)){
    # number of unique taxonomy that map to an OTU in e_data #
    if(!is.null(taxa_cname)){
      num_taxa = length(unique(e_meta[which(as.character(e_meta[, edata_cname]) %in% as.character(e_data[, edata_cname])), taxa_cname]))
    }else{num_taxa = NULL}
  }else{
    num_taxa = NULL
  }

  if(!is.null(e_meta)){
    # number of unique ec that map to an OTU in e_data #
    if(!is.null(ec_cname)){
      num_ec = length(unique(e_meta[which(as.character(e_meta[, ec_cname]) %in% as.character(e_data[, edata_cname])), ec_cname]))
    }else{num_ec = NULL}
  }else{
    num_ec = NULL
  }

  if(!is.null(e_meta)){
    # number of unique genes that map to an OTU in e_data #
    if(!is.null(gene_cname)){
      num_gene = length(unique(e_meta[which(as.character(e_meta[, edata_cname]) %in% as.character(e_data[, edata_cname])), gene_cname]))
    }else{num_gene = NULL}
  }else{
    num_gene = NULL
  }

  # set data information attributes #
  attr(seqObj, "data_info") = list(data_type = data_type, data_scale = data_scale, data_norm = data_norm,
                                norm_method = norm_method, location_param = location_param, scale_param = scale_param,
                                num_samps = num_samps, num_edata = num_edata, num_taxa = num_taxa, num_ec = num_ec, num_gene = num_gene,
                                num_na = num_na, frac_na = frac_na, num_zero = num_zero, frac_zero = frac_zero,
                                #num_miss_obs = num_miss_obs, frac_missing = frac_missing,
                                seq_type = seq_type)

  # set meta data attributes #
  if(!is.null(e_meta)){
    attr(seqObj, "meta_info") = TRUE
  }else{ attr(seqObj, "meta_info") = FALSE}

  # set database attributes #
  if(is.null(db)){
    db_version = NULL
  }else{
    db_version = db_version
  }
  attr(seqObj, "database") = list(db = db, db_version = db_version)

  # set group dataframe attribute to NULL, will be filled in after running group_designation function #
  attr(seqObj, "group_DF") = NULL

  # set class of list #
  class(seqObj) = "seqData"

  return(seqObj)

}
