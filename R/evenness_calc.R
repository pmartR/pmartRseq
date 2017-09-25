#' Calculates Evenness Indices
#'
#' This function calculates Shannon's, Simpson's, and/or Inverse Simpson's evenness indices.
#'
#' @param omicsData an object of the class 'seqData' reated by \code{\link{as.seqData}}.
#' @param index a character vector stating which of the calculations to perform - "shannon" for Shannon's evenness index and/or "simpson" for Simpson's evenness index. Default is to perform both calculations.
#'
#' @details Evenness is calculated for each sample in the data. Evenness is calculated using Shannon's evenness index or Simpson's evenness index.
#'
#' @return An object of class evenRes (also a data.frame) containing the evenness value(s) for every sample in the data object.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_even <- evenness_calc(omicsData = rRNA_data)
#' rRNA_even
#' summary(rRNA_even)
#' plot(rRNA_even)
#' }
#'
#' @author Allison Thompson
#'
#' @export
evenness_calc <- function(omicsData, index=c("shannon","simpson")){

  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if(!class(omicsData) %in% c("seqData")) stop("omicsData must be of class 'seqData'")

  if(attr(omicsData, "data_info")$data_scale!='count'){
    warning("This function is meant for count data like 'rRNA', 'gDNA' or 'cDNA' data.")
  }

  ## end initial checks ##

  # change 0 to NA, makes for easier calculation
  #omicsData <- helper_edata_replace(omicsData=omicsData, x=0 , y=NA)
  omicsData$e_data[omicsData$e_data==0] <- NA
  rownames(omicsData$e_data) <- omicsData$e_data[,which(colnames(omicsData$e_data) == attr(omicsData,"cnames")$edata_cname)]
  omicsData$e_data <- omicsData$e_data[,-which(colnames(omicsData$e_data) == attr(omicsData,"cnames")$edata_cname)]

  res <- list()
  results <- list()

  # Calculate Shannon's evenness index
  if("shannon" %in% tolower(index)){
    shannon <- apply(omicsData$e_data, 2, function(x){
      pd <- sum(x, na.rm=TRUE)
      S <- length(which(!is.na(x) & x!=0))
      (-1 * sum((x/pd) * log(x/pd), na.rm=TRUE)) / log(S)
    })
    res[["shannon"]] <- shannon
    #attr(results, "indexValueRange")$shannon <- data.frame(Min=min(res[["shannon"]], na.rm=TRUE), Max=max(res[["shannon"]], na.rm=TRUE))
  }

#   # Calculate Simpson's evenness index
#   if("simpson" %in% tolower(index)){
#     simpson <- apply(omicsData$e_data, 2, function(x){
#       pd <- sum(x, na.rm=TRUE)
#       1 - sum((x/pd)^2, na.rm=TRUE)
#     })
#     res[["simpson"]] <- simpson
#     attr(results, "indexValueRange")$simpson <- data.frame(Min=min(res[["simpson"]], na.rm=TRUE), Max=max(res[["simpson"]], na.rm=TRUE))
#   }

  # Calculate inverse Simpson's evenness index
  if("simpson" %in% tolower(index)){
    simpson <- apply(omicsData$e_data, 2, function(x){
      pd <- sum(x, na.rm=TRUE)
      S <- length(which(!is.na(x) & x!=0))
      (1 / sum((x/pd)^2, na.rm=TRUE)) * (1/S)
    })
    res[["simpson"]] <- simpson
    #attr(results, "indexValueRange")$invsimpson <- data.frame(Min=min(res[["invsimpson"]], na.rm=TRUE), Max=max(res[["invsimpson"]], na.rm=TRUE))
  }

  # Combine results
  res <- do.call(rbind, res)
  res <- as.data.frame(res)

  results <- res

  attr(results,"index") <- index
  attr(results, "group_DF") <- attr(omicsData, "group_DF")
  attr(results, "cnames") <- attr(omicsData, "cnames")
  class(results) <- c("evenRes", class(res))

  return(results)

}
