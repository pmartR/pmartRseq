#' Splits an e_meta column into separate columns
#'
#' This function can take one large string (like a taxonomy) and split it into a desired number of different columns. For instance, split a large string containing a taxonomy into columns for Kingdom, Phylum, etc.
#'
#' @param omicsData an object of one of the class "seqData"
#' @param cname the column in e_meta to split. If NULL, will use attr(omicsData,"cnames")$emeta_cname
#' @param split1 variable used to split string into columns. Default is ";".
#' @param numcol number of columns to split given column into. If NULL, will calculate from data. Default is NULL.
#' @param split2 variable used to split string again, if needed. For example, can split "k__Bacteria" into "Bacteria" and will remove the "k__" part. Default is "__".
#' @param num which part to keep after splitting for the second time. Default is 2.
#' @param newnames new column names after splitting is done. If NULL, will use c("Kingdom","Phylum","Class","Order","Family","Genus","Species"). Default is NULL.
#'
#' @return seqData object with a split e_meta
#'
#' @author Allison Thompson
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' load(rRNA_data)
#'
#' omicsData <- rRNA_data
#' head(omicsData$e_meta)
#'
#' omicsData_split <- split_emeta(omicsData, emeta_cname="Taxonomy", split1=";", numcol=7, split2="__", num=2, newnames=NULL)
#' head(omicsData_split$e_meta)
#' }
#'
#' @export
split_emeta <- function(omicsData, cname=NULL, split1=";", numcol=NULL, split2="__", num=2, newnames=NULL){
  library(stringr)

  ## inital checks ##

  if(is.null(omicsData$e_meta)){
    stop("e_meta must contain data")
  }

  ## end initial checks ##

  # pull out relevant information
  emeta <- omicsData$e_meta

  if(is.null(cname)){
    cname <- attr(omicsData,"cnames")$taxa_cname
  }

  # calculate number of taxonomy columns
  if(!is.null(split1)){
    if(is.null(numcol)){
      numcol <- max(unlist(lapply(lapply(emeta[,which(colnames(emeta)==cname)], function(x){
        strsplit(as.character(x), split1)
        }), function(y){
          length(y[[1]])
          }
      )))
    }
  }else{
    numcol <- ncol(emeta) - 1
  }

  # split taxonomy column into calculated number of columns
  if(!is.null(split1)){
    split_data <- stringr::str_split_fixed(emeta[,cname], split1, numcol)
  }else{
    split_data <- emeta[,-which(colnames(emeta)==cname)]
  }

  # if necessary, can split again on another character, this time will remove other part (eg k__)
  if(!is.null(split2)){
    split_data <- apply(split_data, 1:2, function(x) strsplit(as.character(x), split2)[[1]][num])
  }

  # put the split data back together
  split_data <- as.data.frame(split_data)
  split_data <- apply(split_data, 1:2, as.character)
  split_data[is.na(split_data)] <- "Unknown"
  split_data <- data.frame(OTU=emeta[,which(colnames(emeta) == cname)], split_data)
  colnames(split_data)[1] <- cname

  # add split data to e_meta
  res <- omicsData
  res$e_meta <- merge(split_data, emeta, by=cname)

  # update column names
  if(is.null(newnames)){
    names.vec <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    colnames(res$e_meta) <- c(cname,names.vec[1:numcol],colnames(emeta)[-which(colnames(emeta)==cname)])
  }else{
    colnames(res$e_meta) <- c(cname,newnames,colnames(emeta)[-which(colnames(emeta)==cname)])
  }

  return(res)
}
