#' Plot Indicator species Analysis Results
#'
#' This function plots the results from an indicator species analysis.
#'
#' @param indsp an indicator species analysis results object, created by \code{\link{indsp_calc}}
#' @param omicsData an object of the class 'seqData' created by \code{\link{as.seqData}}.
#' @param x_axis which grouping value to put on the x-axis. Default is "Group".
#' @param group Which taxonomic level to group the results by. Default is "Phylum".
#'
#' @details This function creates a heatmap-type plot showing the results of an indicator species analysis.
#'
#' @return A plot of indicator species analysis results.
#'
#' @references De Cáceres, M. and Legendre, P. 2009. Associations between species and groups of sites: indices and statistical inference. Ecology 90(12): 3566-3574.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_data <- group_designation(omicsData = rRNA_data, main_effects = c("site", "treatment"), time_course=NULL)
#' rRNA_norm <- normalize_data(omicsData = rRNA_data, norm_fn = "css", normalize = TRUE)
#'
#' rRNA_indsp <- indsp_calc(omicsData = rRNA_norm, within=NULL, pval_thresh=0.05)
#' head(rRNA_indsp)
#' summary(rRNA_indsp)
#' plot_indsp(indsp = rRNA_indsp, omicsData = rRNA_data)
#' }
#'
#' @author Allison Thompson
#'
#' @export

plot_indsp <- function(indsp, omicsData, x_axis = "Group", group = "Phylum"){

  library(reshape)
  library(dplyr)
  library(ggplot2)

  # format normalized data
  attr(omicsData, "group_DF")$Group <- gsub("-","\\.",attr(omicsData, "group_DF")$Group)
  groupDF <- unique(attr(omicsData, "group_DF")[,-which(colnames(attr(omicsData, "group_DF")) == attr(omicsData,"cnames")$fdata_cname)])
  normdata <- omicsData$e_data
  normdata <- melt(normdata)
  normdata <- merge(normdata, attr(omicsData,"group_DF"), by.x="variable", by.y=attr(omicsData,"cnames")$fdata_cname)

  # get mean normalized abundance for every group for every feature
  vars <- c(attr(omicsData,"cnames")$edata_cname, x_axis)
  vars <- lapply(vars, as.symbol)
  normgp <- normdata %>%
    dplyr::group_by_(.dots=vars) %>%
    dplyr::summarise(MeanNorm=mean(value, na.rm=TRUE))

  # format indicator species results
  colnames(indsp) <- gsub("-","\\.",colnames(indsp))
  myinds <- data.frame(OTU=rownames(indsp), indsp)
  colnames(myinds)[1] <- attr(omicsData, "cnames")$edata_cname

  # separate out the necessary information
  if(!is.null(attr(indsp,"within"))){
    alldata <- lapply(unique(attr(omicsData,"group_DF")[,attr(indsp,"within")]), function(x){
                    temp <- myinds[,c(1,grep(x, colnames(myinds)))]
                    ids <- grep("s\\.",colnames(temp))
                    temp <- merge(melt(temp[,c(1,ids)]), temp[,-ids], by=attr(omicsData,"cnames")$edata_cname)
                    colnames(temp)[grep("index",colnames(temp))] <- "index"
                    colnames(temp)[grep("stat",colnames(temp))] <- "stat"
                    colnames(temp)[grep("p.value",colnames(temp))] <- "p.value"
                    colnames(temp)[grep("Flag",colnames(temp))] <- "Flag"
                    return(temp)
                  })
    alldata <- do.call(rbind, alldata)
  }else{
    ids <- grep("s\\.",colnames(myinds))
    alldata <- merge(melt(myinds[,c(1,inds)]), myinds[,-inds], by=attr(omicsData,"cnames")$edata_cname)
  }

  # comine indicator species results with e_meta information
  alldata <- merge(alldata, omicsData$e_meta, by=attr(omicsData, "cnames")$edata_cname)
  alldata$variable <- gsub("s\\.","",alldata$variable)
  alldata <- merge(alldata, groupDF, by.x = "variable", by.y = "Group")
  colnames(alldata)[which(colnames(alldata)=="variable")] = "Group"

  # separate out only the significant features
  alldata$indsp <- sapply(c(1:nrow(alldata)), function(x){
    ifelse(alldata$value[x]  == 1 & alldata$p.value[x] <= attr(indsp, "Threshold") & !is.na(alldata$p.value[x]), 1, 0)
  })
  alldata$indsp[is.na(alldata$indsp)] <- 0

  # combine indicator species results with normalized abundances
  alldata <- merge(alldata, normgp, by=c(attr(omicsData,"cnames")$edata_cname, "Group"))

  # separate out only the significant features
  sig <- unique(alldata[which(alldata$indsp == 1),attr(omicsData,"cnames")$edata_cname])
  alldata <- alldata[which(alldata[,attr(omicsData,"cnames")$edata_cname] %in% sig),]

  # order the important variables
  alldata[,group] <- as.factor(alldata[,group])
  alldata <- alldata[order(alldata[,group]),]
  alldata[,attr(omicsData,"cnames")$edata_cname] <- factor(alldata[,attr(omicsData,"cnames")$edata_cname],
                                                           levels=unique(alldata[order(alldata[,group]),attr(omicsData,"cnames")$edata_cname]))

  # make better colors
  TaxonBarPallette <- c("#FF00BB","#CC00FF","#F2BFFF","#7A0099","#0022FF","#8091FF","#001499","#00F2FF","#CCFCFF","#009199","#00D90E","#BFFFC4","#007308","#FFFF00","#DDFF00","#B3B300","#FF9100","#FFC880","#995700","#FF0000","#FFABAB","#990000","#BFBFBF","#636363","#000000")

  map <- aes_string(x=x_axis, y=attr(omicsData,"cnames")$edata_cname, fill=group)
  ggplot(alldata, map)+
    geom_tile(alpha=0.25)+
    geom_tile(aes(alpha=ifelse(indsp==1,1,0)))+
    guides(alpha=FALSE)+
    scale_fill_manual(values=rep(TaxonBarPallette,3))+
    #geom_text(aes(label=round(MeanNorm, digits=1)), cex=3)+
    theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title=element_text(size=18, face="bold"),
          title=element_text(size=18, face="bold"),
          legend.title=element_text(size=18, face="bold"),
          axis.text.x=element_text(size=9,angle=45,vjust=1,hjust=1),
          axis.text.y=element_text(size=7),
          legend.text=element_text(size=12))+
    labs(title="Indicator Species", y=attr(omicsData,"cnames")$edata_cname, x=x_axis)+
    scale_y_discrete(limits=levels(alldata[,attr(omicsData,"cnames")$edata_cname]))+
    #coord_flip()+
    theme(aspect.ratio=7/3)


}
