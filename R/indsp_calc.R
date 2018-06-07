#' Computes Indicator Species Analysis for Count Data
#'
#' This function is a wrapper for an indicator species analysis using the package \code{indicspecies}, studying the association between species patterns and combinations of group sites.
#'
#' @param omicsData an object of the class 'seqData' created by \code{\link{as.seqData}}.
#' @param within If desired, can run indicator species analysis between one effect and within another effect.
#' @param pval_thresh P-value threshold for creating flags for significance. The default is 0.05.
#' @param max_grp Integer indicating the maximum number of group combinations to be considered. If NULL, will use n-1 groups. Max is n-1 groups.
#'
#' @details This function uses the \code{multipatt} function in the package \code{indicspecies} to create combinations of input clusters and compare each combination with the species in \code{e_data}. Please refer to the \code{multipatt} function in the package \code{indicspecies} for more information.
#'
#' @return An object of class indspRes (also a data.frame) containing the indicator species analysis of the data.
#'
#' @references De Cáceres, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_data <- group_designation(omicsData = rRNA_data, main_effects = c("site", "treatment"), time_course=NULL)
#'
#' rRNA_indsp <- indsp_calc(omicsData = rRNA_data, within=NULL, pval_thresh=0.05)
#' head(rRNA_indsp)
#' summary(rRNA_indsp)
#' plot(rRNA_indsp)
#' }
#'
#' @author Allison Thompson
#'
#' @export
indsp_calc <- function(omicsData, within=NULL, pval_thresh=0.05, max_grp=NULL){
  library(indicspecies)

    ## initial checks ##

  # check that omicsData is of appropriate class #
  if(!class(omicsData) %in% c("seqData")) stop("omicsData must be of class 'seqData'")

  if(attr(omicsData, "data_info")$data_scale!='count'){
    warning("This function is meant for count data like 'rRNA', 'gDNA' or 'cDNA' data.")
  }

  if(!attr(omicsData, "data_info")$data_norm){
    warning("We suggest normalizing before running this function.")
  }

  if(is.null(attr(omicsData, "group_DF"))){
    stop("Need to run group_designation function first.")
  }
  
  if(!is.null(max_grp)){
    if(!is.numeric(max_grp) | max_grp < 1){
      stop("max_grp needs to be either NULL or a numeric integer.")
    }
    if(max_grp > length(unique(attr(omicsData, "group_DF")$Group)) - 1){
      stop("max_grp must be no greater than n - 1, where n is the number of groups.")
    }
  }

  ## end of initial checks ##

  # format e_data
  #indspec <- mintR:::helper_edata_replace(omicsData, x=NA , y=0)$e_data
  indspec <- omicsData$e_data
  indspec[is.na(indspec)] <- 0

  rownames(indspec) <- omicsData$e_data[,attr(omicsData, "cnames")$edata_cname]
  indspec <- indspec[,-which(colnames(indspec) == attr(omicsData,"cnames")$edata_cname)]

  indspec <- t(indspec)
  indspec <- as.data.frame(indspec)

  # format groups
  group <- attr(omicsData, "group_DF")$Group[match(rownames(indspec), attr(omicsData,"group_DF")[,attr(omicsData,"cnames")$fdata_cname])]

  if(is.null(within)){
    # indicator species analysis from indicspecies
    IS <- indicspecies::multipatt(indspec, group, func="IndVal.g", control=how(nperm=999), max.order = max_grp)

    # format results
    results <- IS$sign
    results$Flag <- ifelse(results$p.value <= pval_thresh, 1, 0)

  }else{
    attr(omicsData,"group_DF")[,within] <- as.factor(attr(omicsData,"group_DF")[,within])
    fdata_cname <- attr(omicsData, "cnames")$fdata_cname
    IS <- lapply(levels(attr(omicsData, "group_DF")[,within]), function(x){
      # subset grouping information to only those samples within a specified group
      meta.data <- subset(attr(omicsData, "group_DF"), get(within) == x )

      # subset data to only those samples within a specified group
      indspec2 <- indspec[which(rownames(indspec) %in% meta.data[,fdata_cname]), ]

      # re-order meta.data
      meta.data <- meta.data[match(rownames(indspec2),meta.data[,fdata_cname]),]

      # indicator species analysis
      IS <- indicspecies::multipatt(indspec2, as.character(meta.data$Group), func="IndVal.g", control=how(nperm=999), max.order = max_grp)
      results <- IS$sign
      results$Flag <- ifelse(results$p.value <= pval_thresh, 1, 0)

      # format results
      idx1 <- grep("index", colnames(results))
      idx2 <- grep("stat", colnames(results))
      idx3 <- grep("p.value", colnames(results))
      idx4 <- grep("Flag", colnames(results))
      idx <- c(idx1, idx2, idx3, idx4)
      colnames(results)[idx] <- unlist(lapply(idx, function(y) paste(colnames(results)[y], x, sep="_")))

      # order results so cbind doesn't mess anything up
      results <- results[order(rownames(results)),]
      return(results)
    })
    results <- do.call(cbind, IS)
  }

  # results <- cbind(rownames(results), results)
  # colnames(results)[1] <- attr(omicsData, "cnames")$edata_cname
  # rownames(results) <- NULL
  #
  res <- results

  # make an indicator species object
  attr(res, "group_DF") <- attr(omicsData, "group_DF")
  attr(res, "cnames") <- attr(omicsData, "cnames")
  attr(res, "Threshold") <- pval_thresh
  attr(res, "within") <- within
  class(res) <- c("indspRes", class(results))

  return(res)

}
