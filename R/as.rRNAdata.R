#' Convert Data to Appropriate mintR Class
#'
#' Converts a list object or several data.frames of rRNA-level data (16s) to an
#' object of the class 'rRNAdata'. Objects of the class 'rRNAdata' are lists
#' with two obligatory components \code{e_data} and \code{f_data}. An optional
#' list component \code{e_meta} is used if analysis or visualization at other
#' levels (e.g. taxonomy) is also desired.
#'
#' @param e_data a \eqn{p \times n + 1} data.frame of expression data, where
#'   \eqn{p} is the number of accession numbers observed and \eqn{n} is the
#'   number of samples (an additional feature identifier/name column should also
#'   be present anywhere in the data.frame). Each row corresponds to data for
#'   each accession number.
#' @param f_data a data.frame with \eqn{n} rows. Each row corresponds to a
#'   sample with one column giving the unique sample identifiers found in e_data
#'   column names and other columns providing qualitative and/or quantitative
#'   traits of each sample.
#' @param e_meta an optional data.frame with \eqn{p} rows. Each row corresponds
#'   to an OTU with one column giving OTU identifiers (must be named the same as
#'   the column in \code{e_data}) and other columns giving meta information
#'   (e.g. mappings of OTU identification to taxonomy).
#' @param tree_path an optional quoted file path to either a NEXUS or Newick
#'   formatted phylogenetic tree file. The OTU labels in the tree file should
#'   match the OTU identifiers in the preceeding data fields.
#' @param fasta_path an optional quoted file path to a fasta formatted
#'   representation of biological sequences. Each OTU in the fasta maps to at
#'   least one sequence in the preceeding data fields.
#' @param edata_cname character string specifying the name of the column
#'   containing the identifiers in \code{e_data} and \code{e_meta} (if
#'   applicable).
#' @param emeta_cname character string specifying the name of the column
#'   containing the mapped identifiers in \code{e_meta} (if applicable).
#'   Defaults to NULL. If \code{e_meta} is NULL, then specify \code{emeta_cname}
#'   as NULL.
#' @param fdata_cname character string specifying the name of the column
#'   containing the sample identifiers in \code{f_data}.
#' @param ... further arguments some of which can be passed to \code{read.tree()}
#'   and \code{readDNAStringSet}
#'
#' @details Objects of class 'rRNAdata' contain some attributes that are
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
#'   normalization scale parameters for each sample. \cr \tab \cr data_types
#'   \tab Character string describing the type of data (e.g. 'HiSeq' or
#'   'Positive ion'). Default value is 'NULL'. \cr \tab \cr db \tab Character
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
#'   \code{e_data} values that are 0. \cr \tab \cr num_emeta \tab The number of
#'   unique \code{emeta_cname} entries. \cr \tab \cr num_samps \tab The number
#'   of samples that make up the columns of \code{e_data}.\cr \tab \cr meta_info
#'   \tab A logical argument, specifying whether \code{e_meta} is provided.\cr
#'   \tab \cr }
#'
#' @author Allison Thompson, Lisa Bramer
#' @seealso \code{\link{as.gDNAdata}}
#' @seealso \code{\link{as.cDNAdata}}
#'
#'
#' @export

as.rRNAdata <- function(e_data, f_data, e_meta=NULL, tree_path = NULL, fasta_path = NULL, edata_cname, fdata_cname, emeta_cname=NULL, ...){
  .as.rRNAdata(e_data, f_data, e_meta, tree_path, fasta_path, edata_cname, fdata_cname, emeta_cname, ...)
}

## rRNA data ##
.as.rRNAdata <- function(e_data, f_data, e_meta=NULL, tree_path, fasta_path,
                         edata_cname, fdata_cname, emeta_cname,
                         data_scale="count", data_norm=FALSE, norm_method=NULL,
                         location_param=NULL, scale_param=NULL,
                         data_types=NULL, db=NULL, db_version=NULL, ...){
  # initial checks #

  # check that e_data and f_data are data.frames #
  if(class(e_data) != "data.frame") stop("e_data must be of the class 'data.frame'")
  if(class(f_data) != "data.frame") stop("f_data must be of the class 'data.frame'")

  # check to see if e_meta is NULL, if not check that it is a data.frame #
  if(!is.null(e_meta)){ if(class(e_meta) != "data.frame") stop("e_meta must be of the class 'data.frame'")}

  # check to see if tree_path is NULL, if not check that it is a character #
  browser()
  if(!is.null(tree_path)){
    if(class(tree_path) != "character") stop("tree_path must be of the class 'character'")
    e_tree <- ape::read.tree(tree_path)
    }

  # check to see if fasta_path is NULL, if not check that it is a character #
  if(!is.null(fasta_path)){
    if(class(fasta_path) != "character") stop("fasta_path must be of the class 'character'")
    e_fasta <- Biostrings::readDNAStringSet(fasta_path, ...)
    }

  # check that the OTU column exists in e_data and e_meta (if applicable) #
  if(!(edata_cname %in% names(e_data))) stop(paste("OTU column ", edata_cname," not found in e_data. See details of as.rRNAdata for specifying column names.", sep = ""))

  if(!is.null(e_meta)){
    if(!(edata_cname %in% names(e_meta))) stop(paste("OTU column ", edata_cname," not found in e_meta. Column names for OTUs must match for e_data and e_meta. See details of as.rRNAdata for specifying column names.", sep = ""))
  }

  # if e_meta is NULL set emeta_cname to NULL #
  if(is.null(e_meta)) emeta_cname = NULL

  # if e_meta is not NULL check that the taxonomy column is found #
  if(!is.null(e_meta)){
    if(!is.null(emeta_cname)){
      if(!(emeta_cname %in% names(e_meta))) stop(paste("Taxonomy column ", emeta_cname, " not found in e_meta. See details of as.rRNAdata for specifying column names.", sep = "") )
    }
  }

  # check that the Sample column name is in f_data column names #
  if(!(fdata_cname %in% names(f_data))) stop(paste("Sample column ", fdata_cname, " not found in f_data. See details of as.rRNAdata for specifying column names.", sep = ""))

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
  if(!(data_scale %in% c('log2', 'log10', 'log', 'count', 'abundance'))) stop(paste(data_scale, " is not a valid option for 'data_scale'. See details of as.rRNAdata for specifics.", sep=""))

  # if e_meta is NULL, set emeta_cname to NULL #
  if(is.null(e_meta)){
    emeta_cname = NULL
  }

  # store results #
  rRNAobj = list(e_data = e_data, f_data = f_data, e_meta = e_meta, e_tree = e_tree, e_fasta = e_fasta)

  # set column name attributes #
  attr(rRNAobj, "cnames") = list(edata_cname = edata_cname, emeta_cname = emeta_cname, fdata_cname = fdata_cname)

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
    if(!is.null(emeta_cname)){
      num_emeta = length(unique(e_meta[which(as.character(e_meta[, edata_cname]) %in% as.character(e_data[, edata_cname])), emeta_cname]))
    }else{num_emeta = NULL}
  }else{
    num_emeta = NULL
  }

  # set data information attributes #
  attr(rRNAobj, "data_info") = list(data_scale = data_scale, data_norm = data_norm, norm_method = norm_method,
                                location_param = location_param, scale_param = scale_param,
                                num_samps = num_samps, num_edata = num_edata, num_emeta = num_emeta,
                                num_na = num_na, frac_na = frac_na, num_zero = num_zero, frac_zero = frac_zero,
                                #num_miss_obs = num_miss_obs, frac_missing = frac_missing,
                                data_types = data_types)

  # set meta data attributes #
  if(!is.null(e_meta)){
    attr(rRNAobj, "meta_info") = TRUE
  }else{ attr(rRNAobj, "meta_info") = FALSE}

  # set database attributes #
  if(is.null(db)){
    db_version = NULL
  }else{
    db_version = db_version
  }
  attr(rRNAobj, "database") = list(db = db, db_version = db_version)

  # set group dataframe attribute to NULL, will be filled in after running group_designation function #
  attr(rRNAobj, "group_DF") = NULL

  # set class of list #
  class(rRNAobj) = "rRNAdata"

  return(rRNAobj)

}
