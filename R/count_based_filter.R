#' Count filter object
#'
#' This function returns a countFilter object, performing all of the
#' calculations needed to filter the data based off a specified function and
#' limit
#'
#' @param omicsData An object of one of the classes "seqData"
#'
#' @param fn Specify "mean" to use the mean count of each OTU, "percent" to use
#'   mean counts lower than a certain percent, "max" to use the max count across
#'   all samples, "sum" to use the total count of each OTU, "nonmiss" to use
#'   presence/absence counts, or "ka" to use k over a filtering (need at least
#'   k counts of OTUs seen in at least a samples).
#'
#' @return An object of class countFilter (also a data.frame) that contains the
#'   molecule identifier and the mean/percent/max/sum/nonmissing count across
#'   all samples.
#'
#' @author Allison Thompson, Bryan Stanfill
#'
#' @references Arumugam, Manimozhiyan, et al. "Enterotypes of the human gut
#' microbiome." nature 473.7346 (2011): 174-180.
#'
#' https://bioinformatics.oxfordjournals.org/content/early/2013/07/15/bioinformatics.btt350.full
#'
#' @examples
#'\dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' omicsData <- rRNA_data
#'
#' #Find mean count of OTUs
#' mean_lim <- count_based_filter(omicsData, fn="mean")
#' head(mean_lim)
#' summary(mean_lim)
#' plot(mean_lim)
#'
#' #Find percentage of each OTU
#' perc_lim <- count_based_filter(omicsData, fn="percent")
#' head(perc_lim)
#' summary(perc_lim)
#' plot(perc_lim)
#'
#' #Find maximum count of OTUs
#' max_lim <- count_based_filter(omicsData, fn="max")
#' head(max_lim)
#' summary(max_lim)
#' plot(max_lim)
#'
#' #Find total count of OTUs
#' sum_lim <- count_based_filter(omicsData, fn="sum")
#' head(sum_lim)
#' summary(sum_lim)
#' plot(sum_lim)
#'
#' #Find number of nonmissing OTUs
#' nonmiss_lim <- count_based_filter(omicsData, fn="nonmiss")
#' head(nonmiss_lim)
#' summary(nonmiss_lim)
#' plot(nonmiss_lim)
#'
#' #Find order of samples for k/a filtering
#' ka_lim <- count_based_filter(omicsData, fn="ka")
#' head(ka_lim)
#' summary(ka_lim)
#' plot(ka_lim)
#' }
#'
#' @export
count_based_filter <- function(omicsData, fn="sum"){

  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if(!class(omicsData) %in% c("seqData")) stop("omicsData must be of class 'seqData'")

  if(attr(omicsData, "data_info")$data_scale!='count'){
    warning("This function is meant for count data like 'rRNA', 'gDNA' or 'cDNA' data.")
  }

  if(!(tolower(fn) %in% c("mean","percent","max","sum","nonmiss", "ka"))){
    stop("fn must only be 'mean', 'percent', 'max', 'sum', 'nonmiss', or 'ka'.")
  }

  ## end initial checks ##

  edata <- omicsData$e_data
  edata_cname <- attr(omicsData,"cnames")$edata_cname

  if(fn == "mean"){
    # Mean count of each OTU
    mean_OTUs <- apply(edata[,-which(names(omicsData$e_data) == edata_cname)],1,function(x){return(mean(x,na.rm=TRUE))})
    infrequent_OTUs <- data.frame(omicsData$e_data[,edata_cname], mean_OTUs)
    colnames(infrequent_OTUs) <- c(edata_cname, "meanOTUs")

  }else if(fn == "percent"){
    # Percent of population each OTU makes up
    perc_OTUs <- rowSums(edata[,-which(names(omicsData$e_data) == edata_cname)],na.rm=TRUE)/sum(edata[,-which(names(omicsData$e_data) == edata_cname)],na.rm=TRUE)
    infrequent_OTUs <- data.frame(omicsData$e_data[,edata_cname], perc_OTUs)
    colnames(infrequent_OTUs) <- c(edata_cname, "percentOTUs")

  }else if(fn == "max"){
    # Max count of each OTU
    max_OTUs <- apply(edata[,-which(colnames(edata)==edata_cname)],1,function(x){return(max(x,na.rm=T))})
    infrequent_OTUs <- data.frame(omicsData$e_data[,edata_cname], max_OTUs)
    colnames(infrequent_OTUs) <- c(edata_cname, "maxOTUs")

  }else if(fn == "sum"){
    # Total number of times each OTU is seen
    sum_OTUs <- rowSums(edata[, -which(colnames(edata) == edata_cname)], na.rm=TRUE)
    infrequent_OTUs <- data.frame(omicsData$e_data[, edata_cname], sum_OTUs)
    colnames(infrequent_OTUs) <- c(edata_cname, "sumOTUs")

  }else if(fn == "nonmiss"){
    # Presence/absence of each OTU
    nonmiss_OTUs <- !is.na(omicsData$e_data[, -which(colnames(edata) == edata_cname)])
    infrequent_OTUs <- data.frame(omicsData$e_data[, edata_cname], rowSums(nonmiss_OTUs))
    colnames(infrequent_OTUs) <- c(edata_cname, "nonmissOTUs")

  }else if(fn == "ka"){
    # k/a filtering - an OTU must be seen at least k times in at least a samples
    ka_edata <- edata
    rownames(ka_edata) <- ka_edata[,which(colnames(ka_edata) == edata_cname)]
    ka_edata <- ka_edata[,-which(colnames(ka_edata) == edata_cname)]
    ka_edata <- as.matrix(ka_edata)
    ka_edata[which(is.na(ka_edata))] <- 0

    ka_OTUs <- lapply(c(1:nrow(ka_edata)), function(x) as.vector(ka_edata[x,])[order(as.vector(ka_edata[x,]),decreasing=TRUE)])
    ka_OTUs <- lapply(ka_OTUs, unname)
    ka_OTUs <- do.call(rbind, ka_OTUs)
    colnames(ka_OTUs) <- sapply(c(1:ncol(ka_edata)), function(x) paste("NumSamples_",x,sep=""))

    infrequent_OTUs <- data.frame(omicsData$e_data[, edata_cname], ka_OTUs)
    colnames(infrequent_OTUs)[1] <- edata_cname

  }

  class(infrequent_OTUs) <- c("countFilter",class(infrequent_OTUs))

  attr(infrequent_OTUs, "sample_names") <- names(omicsData$e_data)[-which(names(omicsData$e_data) == edata_cname)]
  attr(infrequent_OTUs, "group_DF") <- attr(omicsData, "group_DF")
  attr(infrequent_OTUs, "function") <- fn

  threshold <- quantile(reshape2::melt(infrequent_OTUs)$value, 0.95)
  attr(infrequent_OTUs, "threshold") <- threshold

  return(infrequent_OTUs)

}
