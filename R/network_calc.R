#' Calculates Correlations
#'
#' This function calculates correlations for use in network analysis.
#'
#' @param omicsData an object of the class 'seqData' reated by \code{\link{as.seqData}}.
#' @param type Character, which correlation metric to use. Must be one of 'pearson' or 'spearman'. Default is 'spearman'.
#' @param group Logical, should correlation analysis be performed on sub-groups of the data. Default is FALSE, will use all of the samples in the data.
#' @param group_var Character, if group==TRUE, which variable to use as the grouping variable.
#' @param fdr_method a character vector stating which cutoff method to use, one of 'fndr', 'pct0', or 'locfdr'. Default is 'fndr'.
#' @param missing_val Use 0 or NA for any 'missing' values? Default is 0.
#'
#' @details Correlation is calculated between every pair of features.
#'
#' @return An object of class corrRes (also a data.frame) containing the correlation values.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' mynetwork <- network_calc(omicsData = rRNA_data)
#' head(mynetwork)
#' }
#'
#' @author Allison Thompson
#'
#' @references
#'
#' @export
network_calc <- function(omicsData, type="spearman", group=FALSE, group_var=NULL, fdr_method="fndr", missing_val=0){

 library(Hmisc)
 library(fdrtool)

  ### Initial Checks ###

  if(class(omicsData)[1] != "seqData"){
    stop("omicsData must be an object of class 'seqData'")
  }

  if(!(type %in% c("spearman","pearson"))){
    stop("type must be one of 'spearman' or 'pearson'.")
  }

  if(!(fdr_method %in% c("fndr","pot0","locfdr"))){
    stop("fdr_method must be one of 'fndr', 'pot0', or 'locfdr'.")
  }

  if(missing_val != 0 & !is.na(missing_val)){
    stop("missing_val must be one of 0 or NA.")
  }

  if(!is.logical(group)){
    stop("group must be TRUE or FALSE")
  }

  if(!(group_var %in% c(colnames(attr(omicsData, "group_DF")), colnames(omicsData$f_data)))){
    stop("group_var must be a column name found in either omicsData$f_data or attr(omicsData, 'group_DF')")
  }

  ### End Initial Checks ###

  # Extract data
  cordata <- omicsData$e_data

  if(missing_val == 0){
    cordata[is.na(cordata)] <- 0
  }else if(is.na(missing_val)){
    cordata[cordata==0] <- NA
  }else{
    stop("missing_val must be 0 or NA")
  }

  # Extract cnames
  edata_cname <- attr(omicsData, "cnames")$edata_cname
  fdata_cname <- attr(omicsData, "cnames")$fdata_cname

  # If want to run network on subgroups
  if(group){

    # Get samples in each group
    if(is.null(group_var) & is.null(attr(omicsData, "group_DF"))){
      stop("In order to perform network calculation across groups, must provide a grouping variable.")
    }else if(is.null(group_var) & !is.null(attr(omicsData, "group_DF"))){
      warning("No grouping variable specified, will use 'Group' found in group_DF.")
      groupings <- unique(attr(omicsData, "group_DF")$Group)
      samps <- lapply(groupings, function(x) colnames(cordata)[which(colnames(cordata) %in% subset(attr(omicsData, "group_DF"), Group == x)[,edata_cname])])
      names(samps) <- groupings
    }else{
      if(!is.null(attr(omicsData, "group_DF")) & group_var %in% colnames(attr(omicsData, "group_DF"))){
        groupings <- unique(attr(omicsData, "group_DF")[,group_var])
        samps <- lapply(groupings, function(x) colnames(cordata)[which(colnames(cordata) %in% attr(omicsData, "group_DF")[which(attr(omicsData, "group_DF")[,group_var]==x),fdata_cname])])
        names(samps) <- groupings
      }else{
        groupings <- unique(omicsData$f_data[,group_var])
        samps <- lapply(groupings, function(x) colnames(cordata)[which(colnames(cordata) %in% omicsData$f_data[which(omicsData$f_data[,group_var] == x), fdata_cname])])
        names(samps)  <- groupings
      }
    }

    if(any(sapply(samps, length) < 4)){
      warning("Grouping leads to a network with less than 4 samples - cannot have less than 4 samples - will remove this group")
      samps <- samps[-which(sapply(samps, length) < 4)]
    }

    if(length(samps) < 1){
      stop("No groups remaining, each group must have at least 4 samples and there must be at least one group.")
    }

    # Run correlation calculation within each group
    res <- lapply(c(1:length(samps)), function(x){
      # Get group
      grp <- names(samps)[x]

      # Subset to samples within group
      cord <- cordata[,which(colnames(cordata) %in% c(edata_cname, samps[[x]]))]

      # Make this a matrix
      rownames(cord) <- cord[,which(colnames(cord) == edata_cname)]
      cord <- cord[,-which(colnames(cord) == edata_cname)]
      if(any(rowSums(cord, na.rm=TRUE) == 0)){
        cord <- cord[-which(rowSums(cord, na.rm=TRUE) == 0),]
      }
      cord <- as.matrix(t(cord))

      # Run correlations
      res <- Hmisc::rcorr(cord, type=type)

      # Exctract correlation coefficient, p-value, and which features are compared
      pairs <- flattenCorrMatrix(res$r, res$P)

      # Calculate FDR
      pairs$q.value <- fdrtool::fdrtool(pairs$p, statistic="pvalue", plot=FALSE, verbose=FALSE, cutoff.method=fdr_method)$qval

      # State which group this is for
      pairs$Group <- grp

      return(pairs)
    })

    # Combine results for all groups
    results <- do.call(rbind, res)
    colnames(results) <- c("Row","Column","cor.coeff","p.value","q.value","Group")

    attr(results, "group_var") <- group_var

  }else{
    # Make this a matrix
    cord <- cordata
    rownames(cord) <- cord[,which(colnames(cord) == edata_cname)]
    cord <- cord[,-which(colnames(cord) == edata_cname)]
    if(any(rowSums(cord, na.rm=TRUE) == 0)){
      cord <- cord[-which(rowSums(cord, na.rm=TRUE) == 0),]
    }
    cord <- as.matrix(t(cord))

    # Run correlation
    res <- Hmisc::rcorr(cord, type=type)

    # Exctract correlation coefficient, p-value, and which features are compared
    pairs <- flattenCorrMatrix(res$r, res$P)

    # Calculate FDR
    pairs$q.value <- fdrtool::fdrtool(pairs$p, statistic="pvalue", plot=FALSE, verbose=FALSE, cutoff.method=fdr_method)$qval

    results <- pairs
    colnames(results) <- c("Row","Column","cor.coeff","p.value","q.value")
  }

  # Format results
  attr(results,"type") <- type
  attr(results, "group_DF") <- attr(omicsData, "group_DF")
  attr(results, "cnames") <- attr(omicsData, "cnames")
  attr(results, "e_meta") <- omicsData$e_meta
  class(results) <- c("corrRes", class(results))

  return(results)

}


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

