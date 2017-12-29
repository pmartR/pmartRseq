#' Metadata filter object
#'
#' This function returns a metaFilter object, performing all of the
#' calculations needed to filter the data based off a specified
#' e_meta level
#'
#' @param omicsData An object of the classes "seqData"
#'
#' @param criteria Specify which omicsData$e_meta column name to filter on
#'
#' @return An object of class metaFilter (also a data.frame) that contains the
#'   sample identifier, the values in the criteria column, and the corresponding
#'   sum abundances.
#'
#' @author Allison Thompson and Sarah Reehl
#'
#' @examples
#'
#' @export
metadata_based_filter <- function(omicsData, criteria) {

  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if (!class(omicsData) %in% c("seqData")) stop("omicsData must be of class 'seqData'")

  if (attr(omicsData, "data_info")$data_scale!='count') {
    warning("This function is meant for count data like 'rRNA', 'gDNA' or 'cDNA' data.")
  }

  if(!(tolower(criteria) %in% tolower(colnames(omicsData$e_meta)))){
    stop("criteria must be a column name in omicsData$e_meta")
  }

  ## end initial checks ##

  # extract relevant data
  emeta <- omicsData$e_meta
  edata_cname <- attr(omicsData, "cnames")$edata_cname

  # get feature data
  infrequent_OTUs <- data.frame(emeta[,which(tolower(colnames(emeta)) %in% c(tolower(edata_cname), tolower(criteria)))])

  # get count of each feature
  sums <- data.frame(OTU=omicsData$e_data[,which(colnames(omicsData$e_data) == edata_cname)], Sum=apply(omicsData$e_data[,-which(colnames(omicsData$e_data) == edata_cname)], 1, function(x) sum(x, na.rm=TRUE)))

  # format output
  infrequent_OTUs <- merge(infrequent_OTUs, sums, by=edata_cname)

  class(infrequent_OTUs) <- c("metaFilter",class(infrequent_OTUs))

  attr(infrequent_OTUs, "group_DF") <- attr(omicsData, "group_DF")
  attr(infrequent_OTUs, "criteria") <- criteria

  threshold <- quantile(table(infrequent_OTUs[,2]), 0.95)
  attr(infrequent_OTUs, "threshold") <- threshold

  return(infrequent_OTUs)

}
