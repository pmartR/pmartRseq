#' Calculates Jaccard Index
#'
#' This function calculates Jaccard's dis/similarity coefficient between group replicates.
#'
#' @param omicsData an object of the class 'seqData' created by \code{\link{as.seqData}}
#' @param sim TRUE/FALSE indicating if you want a value for Jaccard similarity (TRUE) or dissimilarity (FALSE). Default is TRUE.
#'
#' @details Jaccard dis/similarity index is calculated across replicates within groups. The index is then averaged for each replicate.
#'
#' @return An object of class jaccardRes (also a data.frame) containing the Jaccard dis/similarity value for every sample in the data object.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_jaccard <- jaccard_calc(omicsData = rRNA_data, sim = TRUE)
#' rRNA_jaccard
#' summary(rRNA_jaccard)
#' plot(rRNA_jaccard)
#' }
#'
#' @author Allison Thompson
#'
#' @export
jaccard_calc <- function(omicsData, sim=TRUE){

  if(is.null(attr(omicsData, "group_DF"))){
    stop("Must run group designation function first.")
  }

  # Format expression data
  edata <- omicsData$e_data
  rownames(edata) <- edata[,attr(omicsData, "cnames")$edata_cname]
  edata <- edata[,-which(colnames(omicsData$e_data) == attr(omicsData, "cnames")$edata_cname)]
  edata[edata==0] <- NA

  # Pull out group information
  groups <- unique(attr(omicsData, "group_DF")$Group)

  jac <- lapply(groups, function(g){
    # Subset data to only those samples in the group of interest
    samples <- subset(attr(omicsData, "group_DF"), Group==g)[,attr(omicsData,"cnames")$fdata_cname]
    temp <- edata[,colnames(omicsData$e_data)[which(colnames(omicsData$e_data) %in% samples)]]

    # A matrix showing all pairs of samples within the group
    pairs <- combn(colnames(temp), 2)

    # Perform for each pair of samples
    res <- apply(pairs, 2, function(p){
      temp1 <- temp[,p[1]]
      names(temp1) <- rownames(temp)
      temp1 <- temp1[-which(is.na(temp1))]

      temp2 <- temp[,p[2]]
      names(temp2) <- rownames(temp)
      temp2 <- temp2[-which(is.na(temp2))]

      if(sim){
        # Jaccard similarity
        jac_idx <- length(intersect(names(temp1),names(temp2))) / length(union(names(temp1),names(temp2)))
      }else{
        # Jaccard dissimilarity
        jac_idx <- 1 - length(intersect(names(temp1),names(temp2))) / length(union(names(temp1),names(temp2)))
      }

      # Return Jaccard value with group and sample information
      res <- data.frame(Group=g, Sample1=p[1], Sample2=p[2], Jaccard=jac_idx)
      return(res)
    })
    do.call(rbind, res)
  })

  all_res <- do.call(rbind, jac)

  # Median Jaccard score for each sample
  results <- lapply(colnames(edata), function(x){
    rep <- subset(all_res, Sample1==x | Sample2==x)
    med <- median(rep$Jaccard, na.rm=TRUE)
    iqr <- IQR(rep$Jaccard, na.rm=TRUE)
    avg <- mean(rep$Jaccard, na.rm=TRUE)
    sd <- sd(rep$Jaccard, na.rm=TRUE)
    data.frame(Group=unique(rep$Group), Sample=x, Median=med, InterQuartileRange=iqr, Average=avg, StdDev=sd)
  })

  jac_res <- do.call(rbind, results)
  colnames(jac_res)[which(colnames(jac_res) == "Sample")] <- attr(omicsData,"cnames")$fdata_cname

  jac_results <- jac_res

  attr(jac_results, "group_DF") <- attr(omicsData, "group_DF")
  attr(jac_results, "cnames") <- attr(omicsData, "cnames")
  attr(jac_results, "similarity") <- sim
  attr(jac_results, "allRes") <- all_res
  class(jac_results) <- c("jaccardRes", class(jac_res))


  return(jac_results)

}
