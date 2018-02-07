#' Calculates Alpha Diversity Indices
#'
#' This function calculates Shannon's, Simpson's, and/or Inverse Simpson's alpha diversity indices.
#'
#' @param omicsData an object of the class 'seqData' created by \code{\link{as.seqData}}.
#' @param index a character vector stating which of the calculations to perform - "shannon" for Shannon's diversity index, "simpson" for Simpson's diversity index, and/or "invsimpson" for the inverse Simpson's diversity index. Default is to perform all 3 calculations.
#'
#' @details Alpha diversity is calculated for each sample in the data, using Shannon's diversity index (Shannon's H), Simpson's diversity index (Simpson's D), and/or inverse Simpson's diversity index.
#'
#' Shannon or Shannon–Weaver (or Shannon–Wiener) index is defined as H = -sum(p_i*log(p_i)), where p_i is the proportional abundance of species i.
#'
#' Both variants of Simpson's index are based on D = sum(p_i^2). Choice simpson returns 1-D and invsimpson returns 1/D.
#'
#' @return An object of class alphaRes (also a data.frame) containing the diversity value(s) for every sample in the data object.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_diversity <- alphaDiv_calc(omicsData = rRNA_data)
#' rRNA_diversity
#' summary(rRNA_diversity)
#' plot(rRNA_diversity)
#'}
#'
#' @author Allison Thompson
#'
#' @export
alphaDiv_calc <- function(omicsData, index=c("shannon","simpson","invsimpson")){

  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if(!class(omicsData) %in% c("seqData")) stop("omicsData must be of class 'seqData'")

  if(attr(omicsData, "data_info")$data_scale!='count'){
    warning("This function is meant for count data like 'rRNA', 'gDNA' or 'cDNA' data.")
  }

  # Make sure that index matches
  if(!(all(tolower(index) %in% c("shannon","simpson","invsimpson")))){
    stop("Error in index. Must contain only 'shannon', 'simpson', and/or 'invsimpson'.")
  }

  ## end initial checks ##

  # change 0 to NA, makes for easier calculation
  omicsData$e_data[omicsData$e_data == 0] <- NA #mintR:::helper_edata_replace(omicsData=omicsData, x=0 , y=NA)
  rownames(omicsData$e_data) <- omicsData$e_data[,which(colnames(omicsData$e_data) == attr(omicsData,"cnames")$edata_cname)]
  omicsData$e_data <- omicsData$e_data[,-which(colnames(omicsData$e_data) == attr(omicsData,"cnames")$edata_cname)]

  res <- list()
  results <- list()

  # Calculate Shannon's H diversity index
  if("shannon" %in% tolower(index)){
    shannon <- apply(omicsData$e_data, 2, function(x){
      pd <- sum(x, na.rm=TRUE)
      -1 * sum((x/pd) * log(x/pd), na.rm=TRUE)
    })
    res[["shannon"]] <- shannon
  }

  # Calculate Simpson's D diversity index
  if("simpson" %in% tolower(index)){
    simpson <- apply(omicsData$e_data, 2, function(x){
      pd <- sum(x, na.rm=TRUE)
      1 - sum((x/pd)^2, na.rm=TRUE)
    })
    res[["simpson"]] <- simpson
  }

  # Calculate inverse Simpson's diversity index
  if("invsimpson" %in% tolower(index)){
    invsimpson <- apply(omicsData$e_data, 2, function(x){
      pd <- sum(x, na.rm=TRUE)
      1 / sum((x/pd)^2, na.rm=TRUE)
    })
    res[["invsimpson"]] <- invsimpson
  }

  # Combine results
  res <- do.call(rbind, res)
  res <- as.data.frame(res)

  results <- res

  attr(results,"index") <- index
  attr(results, "group_DF") <- attr(omicsData, "group_DF")
  attr(results, "cnames") <- attr(omicsData, "cnames")
  class(results) <- c("alphaRes", class(res))

  return(results)

}
