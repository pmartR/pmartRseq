#' Calculates Abundance
#'
#' This function calculates the total number of features seen in each sample.
#'
#' @param omicsData an object of the class 'seqData' usually created by \code{\link{as.seqData}}.
#'
#' @details Abundance is the total number of features (individuals) seen in each sample.
#'
#' @return An object of class abunRes (also a data.frame) containing the abundance value for every sample in the data object.
#'
#' @examples
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_abundance <- abundance_calc(omicsData = rRNA_data)
#' rRNA_abundance
#' summary(rRNA_abundance)
#' plot(rRNA_abundance)
#'
#' @author Allison Thompson
#'
#' @export

abundance_calc <- function(omicsData){

  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if(!class(omicsData) %in% c("seqData")) stop("omicsData must be of class 'seqData'")

  if(attr(omicsData, "data_info")$data_scale!='count'){
    warning("This function is meant for count data like 'rRNA', 'gDNA' or 'cDNA' data.")
  }

  if(!attr(omicsData, "data_info")$data_norm){
    warning("We suggest normalizing before running this function.")
  }

  ## end initial checks ##

  edata_cname <- attr(omicsData, "cnames")$edata_cname

  # calculate abundance
  abundance <- apply(omicsData$e_data[,-which(colnames(omicsData$e_data)==edata_cname)], 2, function(y) sum(y, na.rm=TRUE))
  abundance <- as.data.frame(abundance)

  # make an abundance object
  attr(abundance, "group_DF") <- attr(omicsData, "group_DF")
  attr(abundance, "cnames") <- attr(omicsData, "cnames")
  class(abundance) <- c("abunRes", class(abundance))

  return(abundance)

}
