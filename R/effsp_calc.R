#' Calculates Effective Species
#'
#' This function calculates the number of effective species for every sample in a count dataset
#'
#' @param omicsData an object of the class 'seqData' reated by \code{\link{as.seqData}}.
#'
#' @details Calculates effective species of count data
#'
#' @return An object of class effspRes (also a data.frame) containing the number of effective species for every sample in the data object.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_effsp <- effsp_calc(omicsData = rRNA_data)
#' rRNA_effsp
#' summary(rRNA_effsp)
#' plot(rRNA_effsp)
#' }
#'
#' @author Allison Thompson
#'
#' @export

effsp_calc <- function(omicsData){

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

  # change 0 to NA, makes for easier calculation
  omicsData$e_data[omicsData$e_data==0] <- NA

  edata_cname <- attr(omicsData, "cnames")$edata_cname

  # calculate effective species
  effectiveSpecies <- apply(omicsData$e_data[,-which(colnames(omicsData$e_data)==edata_cname)], 2, function(x) 1/sum((x/sum(x, na.rm=TRUE))^2, na.rm=TRUE))
  effectiveSpecies <- as.data.frame(effectiveSpecies)

  # make an effective species object
  attr(effectiveSpecies, "group_DF") <- attr(omicsData, "group_DF")
  attr(effectiveSpecies, "cnames") <- attr(omicsData, "cnames")
  class(effectiveSpecies) <- c("effspRes", class(effectiveSpecies))

  return(effectiveSpecies)

}
