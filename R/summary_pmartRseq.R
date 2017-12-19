#' Produce a basic summary of an pmartRseq S3 Object
#'
#' This function will provide basic summary statistics about an object of the class \code{seqData}.
#'
#' @param omicsData an object of the class 'seqData' created by \code{\link{as.seqData}}
#' @param omicsDataFilter an optional filter object for the respective \code{omicsData} class.
#'
#' @return a list of elements
#' \itemize{
#' \item{num.samps} number of samples in \code{e_data}
#' \item{num.edata} number of unique identifying features in \code{e_data}
#' \item{num.emeta} number of unique identifying features in \code{e_meta}
#' \item{num_na} number of observations equal to NA in \code{e_data}
#' \item{frac_na} fraction of observations equal to NA in \code{e_data}
#' \item{num_zero} number of observations equal to 0 in \code{e_data}
#' \item{frac_zero} fraction of observations equal to 0 in \code{e_data}
#' }
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' summary(rRNA_data)
#' }
#'
#' @author Lisa Bramer, Kelly Stratton, Allison Thompson
#'
#' @export


#'@export
#'@rdname summary-pmartRseq
#'@name summary-pmartRseq
summary.seqData <- function(omicsData, omicsDataFilter = NULL){
  # if no filter is provided, just pull info from the pepdata object #
  #if(is.null(omicsDataFilter)){
    res = list(num_samps = attr(omicsData, "data_info")$num_samps,
               num_edata = attr(omicsData, "data_info")$num_edata,
               num_taxa = attr(omicsData, "data_info")$num_taxa,
               num_na = attr(omicsData, "data_info")$num_na,
               frac_na = attr(omicsData, "data_info")$frac_na,
               num_zero = attr(omicsData, "data_info")$num_zero,
               frac_zero = attr(omicsData, "data_info")$frac_zero)

    # construct output #
    newres <- lapply(res, function(x) ifelse(is.null(x), "NA", as.character(x)))
    catmat <- data.frame(unlist(newres, use.names=FALSE))
    colnames(catmat) <- NULL
    rownames(catmat) <- c("Samples ", "Rows (e_data) ", "Rows (taxa) ", "Missing Observations (NA) ", "Fraction Missing (NA) ", "Missing Observations (0) ", "Fraction Missing (0) ")

    cat("\nSummary of 'seqData' Object\n----------------------------")
    cat(capture.output(catmat), sep="\n")
    cat("\n")
    return(invisible(res))
  # }else{
  #   # implement filter to get counts after application #
  #   temp_data = pmartRseq_filter(omicsData, omicsDataFilter)
  #
  #   res = list(num_samps = attr(temp_data, "data_info")$num_samps,
  #              num_edata = attr(temp_data, "data_info")$num_edata,
  #              num_emeta = attr(temp_data, "data_info")$num_emeta,
  #              num_miss_obs = attr(temp_data, "data_info")$num_miss_obs,
  #              frac_missing = attr(temp_data, "data_info")$num_miss_obs/(nrow(temp_data$e_data)*ncol(temp_data$e_data)))
  #
  # }
  #
  # return(res)
}



#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.countSTAT_results <- function(pmartRseq_results){
  fidx <- grep("Flag",colnames(pmartRseq_results$allResults))
  if(length(fidx) > 1){
    num_sig <- apply(pmartRseq_results$allResults[,fidx], 2, function(x) length(which(x != 0)))
  }else{
    num_sig <- length(which(pmartRseq_results$allResults[,fidx] != 0))
  }
  num_sig <- as.data.frame(num_sig)
  res = list()

  res[[1]] <- rbind(PValue_Threshold=attr(pmartRseq_results,"Threshold"),
                    Tests = paste(attr(pmartRseq_results, "Tests")$Test,collapse="; "),
                    Adjustment = attr(pmartRseq_results, "Adjustment"))
  res[[1]] <- as.data.frame(res[[1]])
  colnames(res[[1]]) <- NULL

  res[[2]] <- num_sig
  res[[2]] <- rbind("NumSig",res[[2]])
  rownames(res[[2]])[1] <- "Comparison"
  rownames(res[[2]])[-1] <- unlist(lapply(colnames(pmartRseq_results$allResults)[fidx], function(x) gsub("Flag_","",x)))
  colnames(res[[2]]) <- NULL
  #res <- rbind(res[[1]], res[[2]])

  cat("\nSummary of 'countSTAT' Object\n-----------------------------")
  cat(capture.output(res[[1]]), sep="\n")
  cat(capture.output(res[[2]]), sep="\n")
  cat("\n")
  return(invisible(res))
}


#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.alphaRes <- function(pmartRseq_results){
   res <- t(do.call(rbind,apply(pmartRseq_results, 1, function(x){
     data.frame(Min=min(x), Median=median(x), Mean=mean(x), Max=max(x))
   })))

   res <- as.data.frame(rbind(t(data.frame(Test=colnames(res))),res))
   colnames(res) <- NULL
   cat("\nSummary of 'alphaDiversity' Object\n-----------------------------------")
   cat(capture.output(res), sep="\n")
   cat("\n")

   return(invisible(res))
}


#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.evenRes <- function(pmartRseq_results){
  res <- t(do.call(rbind,apply(pmartRseq_results, 1, function(x){
    data.frame(Min=min(x), Median=median(x), Mean=mean(x), Max=max(x))
  })))

  res <- as.data.frame(rbind(t(data.frame(Test=colnames(res))),res))
  colnames(res) <- NULL
  cat("\nSummary of 'evenness' Object\n-----------------------------")
  cat(capture.output(res), sep="\n")
  cat("\n")

  return(invisible(res))
}

#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.richRes <- function(pmartRseq_results){
  res <- t(do.call(rbind,apply(pmartRseq_results, 1, function(x){
    data.frame(Min=min(x), Median=median(x), Mean=mean(x), Max=max(x))
  })))

  res <- as.data.frame(rbind(t(data.frame(Test=colnames(res))),res))
  colnames(res) <- NULL
  cat("\nSummary of 'richness' Object\n-----------------------------------")
  cat(capture.output(res), sep="\n")
  cat("\n")

  return(invisible(res))
}

#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.abunRes <- function(pmartRseq_results){
  res <- data.frame(Min=min(pmartRseq_results$abundance), Median=median(pmartRseq_results$abundance),
                    Mean=mean(pmartRseq_results$abundance), Max=max(pmartRseq_results$abundance))

  res <- t(res)
  colnames(res) <- "Abundance"
  cat("\nSummary of 'abundance' Object\n-----------------------------\n")
  cat(capture.output(res), sep="\n")
  cat("\n")

  return(invisible(res))
}

#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.jaccardRes <- function(pmartRseq_results){
  res <- apply(pmartRseq_results[,-c(1:2)], 2, function(x){
    data.frame(Min=min(x, na.rm=TRUE), Median=median(x, na.rm=TRUE),
               Mean=mean(x, na.rm=TRUE), Max=max(x, na.rm=TRUE))
  })
  res <- do.call(rbind, res)
  #res <- rbind(colnames(res),res)
  #rownames(res)[1] <- "Statistic"

  #colnames(res) <- NULL

  cat("\nSummary of 'jaccard' Object\n----------------------------\n")
  cat(capture.output(as.data.frame(t(res))), sep="\n")
  cat("\n")
}


#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.effspRes <- function(pmartRseq_results){
  res <- data.frame(Min=min(pmartRseq_results$effectiveSpecies), Median=median(pmartRseq_results$effectiveSpecies),
                    Mean=mean(pmartRseq_results$effectiveSpecies), Max=max(pmartRseq_results$effectiveSpecies))

  res <- t(res)
  colnames(res) <- "EffectiveSpecies"
  cat("\nSummary of 'effectiveSpecies' Object\n-------------------------------------\n")
  cat(capture.output(res), sep="\n")
  cat("\n")

  return(invisible(res))
}


#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
#' @param min_num An integer value specifying the minimum number of fn counts a feature must have to be kept for statistical analysis. If fn="percent", give the decimal number, not the percentage.
summary.countFilter <- function(pmartRseq_results, min_num=NULL){

  if(!is.null(min_num)){
    # check that min_num is not a vector #
    if(length(min_num) > 1){ stop("min_num must be of length 1")}
    # check that min_num is numeric >= 0 #
    if(!is.numeric(min_num) | min_num < 0){ stop("min_num must be an integer >= 0")}
  }

  if(!is.null(attr(pmartRseq_results, "group_var"))){
    # return numeric version of plot, the threshold used, the number that would be tested and the number that would not be tested #
    filt_fn <- attr(pmartRseq_results, "function")

    if(!is.null(min_num)){
      res <- lapply(unique(pmartRseq_results$Group), function(x){
        temp <- pmartRseq_results[which(pmartRseq_results$Group == x),]

        # get number molecules tested
        num_not_tested <- length(which(temp[,paste(filt_fn,"OTUs",sep="")] <= min_num))
        # get number molecules not tested
        num_tested <- length(which(temp[,paste(filt_fn,"OTUs",sep="")] > min_num))

        data.frame(Group=x, Num.Not.Filtered=num_tested, Num.Filtered=num_not_tested)
      })

      res <- do.call(rbind, res)

    }else{
      num_tested = "NULL"
      num_not_tested = "NULL"
      min_num = "NULL"
      res <- list(Num.Not.Filtered=num_tested, Num.Filtered=num_not_tested)
    }

    # create output #
    cat("\nSummary of Count Filter\n-----------------------\n")
    catmat <- data.frame(c("Filter Function:"=filt_fn,"Minimum Number:"=min_num, res))
    #colnames(catmat) <- NULL
    cat(capture.output(print(catmat, row.names=TRUE)), sep="\n")

    return(invisible(res))
  }else{
    # return numeric version of plot, the threshold used, the number that would be tested and the number that would not be tested #
    filt_fn <- attr(pmartRseq_results, "function")

    if(!is.null(min_num)){
      # get number molecules tested
      num_not_tested <- length(which(pmartRseq_results[,paste(filt_fn,"OTUs",sep="")] <= min_num))

      # get number molecules not tested
      num_tested <- length(which(pmartRseq_results[,paste(filt_fn,"OTUs",sep="")] > min_num))

    }else{
      num_tested = "NULL"
      num_not_tested = "NULL"
      min_num = "NULL"
    }

    res <- list(min_num=min_num, num_not_filtered=num_tested, num_filtered=num_not_tested)

    # create output #
    cat("\nSummary of Count Filter\n-----------------------")
    catmat <- data.frame(c("Filter Function:"=filt_fn,"Minimum Number:"=min_num, "Num Filtered:"=res$num_filtered, "Num Not Filtered:"=res$num_not_filtered))
    colnames(catmat) <- NULL
    cat(capture.output(print(catmat, row.names=TRUE)), sep="\n")

    return(invisible(res))
  }
}


#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.indspRes <- function(pmartRseq_results){
  fidx <- grep("Flag",colnames(pmartRseq_results))
  if(length(fidx) > 1){
    num_sig <- apply(pmartRseq_results[,fidx], 2, function(x) length(which(x != 0)))
    num_sig <- data.frame(Comparison=unlist(lapply(names(num_sig),function(x) strsplit(x,"Flag_")[[1]][2])),as.data.frame(num_sig))
  }else{
    num_sig <- length(which(pmartRseq_results[,fidx] != 0))
  }
  res = list(pval_threshold = attr(pmartRseq_results,"Threshold"),
             num_sig = num_sig
  )

  res <- data.frame(pval_threshold=attr(pmartRseq_results,"Threshold"), num_sig)
  #res <- rbind(colnames(res), res)
  rownames(res) <- NULL

  #colnames(res) <- NULL
  cat("\nSummary of 'indicatorSpecies' Object\n-------------------------------------\n")
  cat(capture.output(res), sep="\n")
  cat("\n")

  return(invisible(res))
}


#' @export
#' @rdname summary-pmartRseq
#' @name summary-pmartRseq
summary.paRes <- function(pmartRseq_results){
  library(dplyr)

  num_sig <- pmartRseq_results$results %>% dplyr::group_by(term) %>% dplyr::summarise(NumSig=length(which(p.value <= attr(pmartRseq_results, "pval_thresh")$pval_thresh)))
  #
  # fidx <- grep("Flag",colnames(pmartRseq_results$allResults))
  # if(length(fidx) > 1){
  #   num_sig <- apply(pmartRseq_results$allResults[,fidx], 2, function(x) length(which(x != 0)))
  # }else{
  #   num_sig <- length(which(pmartRseq_results$allResults[,fidx] != 0))
  # }
  # num_sig <- as.data.frame(num_sig)
  # res = list()
  #
  # res[[1]] <- rbind(PValue_Threshold=attr(pmartRseq_results,"Threshold"),
  #                   Tests = paste(attr(pmartRseq_results, "Tests")$Test,collapse="; "),
  #                   Adjustment = attr(pmartRseq_results, "Adjustment"))
  # res[[1]] <- as.data.frame(res[[1]])
  # colnames(res[[1]]) <- NULL
  #
  # res[[2]] <- num_sig
  # res[[2]] <- rbind("NumSig",res[[2]])
  # rownames(res[[2]])[1] <- "Comparison"
  # rownames(res[[2]])[-1] <- unlist(lapply(colnames(pmartRseq_results$allResults)[fidx], function(x) gsub("Flag_","",x)))
  # colnames(res[[2]]) <- NULL
  # #res <- rbind(res[[1]], res[[2]])

  cat("\nSummary of 'paRes' Object\n-----------------------------\n")
  cat(capture.output(as.data.frame(num_sig)), sep="\n")
  cat("\n")
  return(invisible(num_sig))
}
