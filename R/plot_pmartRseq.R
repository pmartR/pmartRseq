#' Produce a plot of an omicsData Object
#'
#' This function will proivde a plot for a \code{omicsData} or \code{omicsData_results} object.
#'
#' @param results_object is an object of class countSTAT_results (\code{\link{countSTAT}}), alphaDiversity (\code{\link{alphaDiversity}}), evenness (\code{\link{evenness}}), jaccard (\code{\link{mint_jaccard}}), countFilter (\code{\link{count_based_filter}}), seqData (\code{\link{as.seqData}}), richness (\code{\link{richness}}), abundance (\code{\link{abundance}}), effectiveSpecies (\code{\link{effectiveSpecies}}), or indicatorSpecies (\code{\link{indicatorSpecies}}).


#' @export
#' @rdname plot_pmartRseq
#' @name plot_pmartRseq
#' @param type Character vector specifying which type of plot to create. "pvals" for a line plot showing the number of significantly expresed biomolecules at different p-value thresholds, "flag" for a bar plot showing the number of significantly expressed biomolecules at the previously specified threshold, "logfc" for a heatmap showing the log2 fold changes of the significantly expressed biomolecules, and/or "volcano" for a volcano plot showing the log2 fold changes versus -log10(pvalue). The default is "pvals".
#' @param test Optional, a character vector specifying which differential expression test to plot. The default is NULL, which will create plots for all tests.
#' @param plot_title Optional, a character vector to use as the plot title
#' @param leglab Optional, a character vector to use as the legend label
#' @param x_lab Optional, a character vector to use as the x-axis label
#' @param y_lab Optional, a character vector to use as the y-axis label
#' @param FCThresh Optional, numeric threshold for the Log2FoldChanges for the volcano plot. Default is log2(1.5)
plot.countSTAT_results <- function(results_object, type="pvals", test=NULL, x_lab=NULL, y_lab=NULL, plot_title=NULL, leglab=NULL, FCThresh=log2(1.5), comparison=NULL, ...){
  .plot.countSTAT_results(results_object, type, test, x_lab, y_lab, plot_title, leglab, FCThresh, comparison, ...)
}

.plot.countSTAT_results <- function(results_object, type="pvals", test=NULL, x_lab=NULL, y_lab=NULL, plot_title=NULL, leglab=NULL, FCThresh=log2(1.5), comparison=NULL){
  library(ggplot2)
  library(reshape2)
  library(gplots)
  library(RColorBrewer)

  ## initial checks ##
  if(!is.null(plot_title)){
    if(!is.character(plot_title)){ stop("plot_title must be a character vector")}
  }
  if(!is.null(x_lab)){
    if(!is.character(x_lab)){ stop("x_lab must be a character vector")}
  }
  if(!is.null(y_lab)){
    if(!is.character(y_lab)){ stop("y_lab must be a character vector")}
  }
  if(!is.null(test)){
    if(!is.character(test)){ stop("test must be a character vector")}
  }
  if(!is.null(type)){
    if(!is.character(type)){ stop("type must be a character vector")}
    if(any(!(type %in% c("pvals","flag","logfc","volcano")))){ stop("type must be at least one of 'pvals', 'flag', 'logfc', or 'volcano'")}
  }
  if(!("countSTAT_results" %in% class(results_object))){
    stop("results_object must be of class countSTAT_results")
  }
  ## end of initial checks ##

  if(is.null(test)){
    test <- attr(results_object,"Tests")$Test
  }

  if(is.null(comparison)){
    comparison <- attr(results_object,"comparisons")$comparison
  }

  PThresh <- attr(results_object, "Threshold")

  ## p-value plot ##
  if("pvals" %in% tolower(type)){
    thresholds <- seq(0,0.4,0.05)

    invisible(lapply(test, function(x){
      idx1 <- grep("padj", names(results_object$allResults))
      idx2 <- grep(x, names(results_object$allResults))
      idx <- intersect(idx1, idx2)

      data <- results_object$allResults[,idx]

      if(length(idx) > 1){
        Counts = lapply(thresholds, function(x){
          counts =apply(data, 2, function(y){
            count = length(which(y < x))
            return(count)
          })
          return(counts)
        })

        Counts <- do.call(rbind, Counts)
        Counts <- data.frame(Threshold=thresholds, Counts)
        Counts$Threshold <- as.factor(Counts$Threshold)
        Counts <- melt(Counts)
        Counts$Threshold <- as.numeric(as.character(Counts$Threshold))
        Counts$Comparison <- unlist(lapply(as.character(Counts$variable), function(x) strsplit(x,"padj_")[[1]][2]))

      }else{
        Counts = lapply(thresholds, function(x) length(which(data < x)))

        Counts <- do.call(rbind, Counts)
        Counts <- data.frame(Threshold=thresholds, Counts)
        Counts$Comparison <- strsplit(colnames(results_object$allResults)[idx],"padj_")[[1]][2]
        names(Counts) <- gsub("Counts","value",names(Counts))
      }


      p1 <- ggplot(Counts, aes(x=Threshold, y=value, colour=Comparison)) +
        geom_line(cex=1.5) +
        theme_bw() +
        theme(axis.line.x = element_line(colour = "black"),
              axis.line.y = element_line(colour="black"),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_blank())
      #labs(x="P-Value Threshold",y="Number of differentially expressed genes") +
      #ggtitle(paste(x," Test",sep=""))

      if(!is.null(x_lab)){
        p1 <- p1 + labs(x=x_lab)
      }else{
        p1 <- p1 + labs(x="P-Value Threshold")
      }

      if(!is.null(y_lab)){
        p1 <- p1 + labs(y=y_lab)
      }else{
        p1 <- p1 + labs(y="Number of differentially expressed genes")
      }

      if(!is.null(plot_title)){
        p1 <- p1 + ggtitle(plot_title)
      }else{
        p1 <- p1 + ggtitle(paste(x," Test",sep=""))
      }

      if(!is.null(leglab)){
        p1 <- p1 + guides(colour=guide_legend(title=leglab))
      }

      print(p1)
    }))
  }

  ## flags barplot ##
  if("flag" %in% tolower(type)){
    pal <- brewer.pal(9, "Set1")

    invisible(lapply(test, function(r){
      idx1 <- grep("Flag", colnames(results_object$allResults))
      idx2 <- grep(r, colnames(results_object$allResults))
      idx <- intersect(idx1, idx2)

      data <- results_object$allResults[,idx]
      data <- melt(data)
      if(length(idx) == 1){data$variable = colnames(results_object$allResults)[idx]}
      data$Pairs <- unlist(lapply(as.character(data$variable), function(x) strsplit(x,"Flag_")[[1]][2]))
      data <- data[-which(data$value==0),]
      data$Direction <- unlist(lapply(data$value, function(x) ifelse(x==1,"Up","Down")))
      data$Direction <- factor(data$Direction, levels=c("Up","Down"))

      plotData <- table(data$Pairs, data$Direction)
      plotData <- melt(plotData)
      names(plotData) <- c("Comparison","Direction","Count")

      p2 <- ggplot(plotData, aes(x=Comparison,y=Count,fill=Direction)) +
        geom_bar(stat="identity", position="dodge") +
        scale_fill_manual(values=pal[c(3,1)], name="Regulatory\nDirection") +
        theme_bw() +
        theme(axis.line.x = element_line(colour = "black"),
              axis.line.y = element_line(colour="black"),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_blank()) +
        geom_text(aes(label=Count), position=position_dodge(width=0.9), vjust=1)
      #labs(x="Comparisons",y="Number of differentially expressed genes") +
      #ggtitle(paste(r, " Test",sep=""))

      if(!is.null(x_lab)){
        p2 <- p2 + labs(x=x_lab)
      }else{
        p2 <- p2 + labs(x="Comparisons")
      }

      if(!is.null(y_lab)){
        p2 <- p2 + labs(y=y_lab)
      }else{
        p2 <- p2 + labs(y="Number of differentially expressed genes")
      }

      if(!is.null(plot_title)){
        p2 <- p2 + ggtitle(plot_title)
      }else{
        p2 <- p2 + ggtitle(paste(r," Test",sep=""))
      }

      if(!is.null(leglab)){
        p2 <- p2 + guides(fill=guide_legend(title=leglab))
      }

      print(p2)
    }))
  }

  ## logfc heatmap ##
  if("logfc" %in% tolower(type)){
    invisible(lapply(test, function(r){
      idx1 <- grep("Flag", colnames(results_object$allResults))
      idx2 <- grep(r, colnames(results_object$allResults))
      idx <- intersect(idx1, idx2)

      sig.genes <- lapply(idx, function(x) which(results_object$allResults[,x] != 0))
      sig.genes <- unlist(sig.genes)

      idx3 <- grep("logFC", colnames(results_object$allResults))
      idx4 <- grep(r, colnames(results_object$allResults))
      idx <- intersect(idx3, idx4)

      data <- results_object$allResults[unique(sig.genes),idx]

      data <- as.matrix(data)

      numCols <- length(seq(from=0, to=max(abs(data)), by=0.1))
      cols1 <- colorRampPalette(colors=rev(brewer.pal(11, "RdYlGn")[6:11]))(numCols)
      cols2 <- colorRampPalette(colors=rev(brewer.pal(11, "RdYlGn")[1:6]))(numCols)

      cols <- c(cols1, cols2)
      breaks <- numCols + numCols + 1

      if(length(sig.genes) > 1 & length(idx) > 1){
        rowLabels <- unlist(lapply(colnames(data), function(x) strsplit(as.character(x),"logFC_")[[1]][2]))

        res <- quantile(data, c(0.02, 0.98))
        a <- as.numeric(res[1])
        b <- as.numeric(res[2])

        data[data < a] <- a # truncate fold change values at 2nd and 98th percentiles
        data[data > b] <- b # for higher resolution in the bulk of the plot

        data <- t(data)

        par(cex.main=.75)
        heatmap.2(data, dendrogram="none", col=cols, breaks=breaks, trace="none", Rowv=FALSE, Colv=TRUE,
                  key=TRUE, keysize=2, labCol="", labRow=rowLabels, cexRow=.8, xlab=x_lab, ylab=y_lab,
                  main=ifelse(is.null(plot_title),paste(r," Test\nLog2FoldChanges",sep=""),plot_title))
      }else if(length(sig.genes) > 1 & length(idx) == 1){
        data <- as.matrix(data[order(data[,1]),])
        rowLabels = unlist(lapply(names(results_object$allResults)[idx], function(x) strsplit(as.character(x),"logFC_")[[1]][2]))
        image(data, col = cols, breaks = seq(min(data),max(data),length.out=breaks), ylab=rowLabels, xaxt="n", yaxt="n", main=ifelse(is.null(plot_title),paste(r," Test\nLog2FoldChanges",sep=""),plot_title))
      }else if(length(sig.genes) == 1 & length(idx) > 1){
        data <- as.matrix(t(data))
        rowLabels = unlist(lapply(names(results_object$allResults)[idx], function(x) strsplit(as.character(x),"logFC_")[[1]][2]))
        image(data, col = cols, breaks = seq(min(data),max(data),length.out=breaks), ylab=rowLabels, xaxt="n", yaxt="n", main=ifelse(is.null(plot_title),paste(r," Test\nLog2FoldChanges",sep=""),plot_title))
      }

    }))
  }

  ## volcano plot ##
  if("volcano" %in% tolower(type)){

    invisible(lapply(test, function(r){

      data <- lapply(comparison, function(c){

        idx1 <- grep(c, colnames(results_object$allResults))
        idx2 <- grep(r, colnames(results_object$allResults))
        idx <- intersect(idx1, idx2)

        data <- results_object$allResults[,idx]
        data$Threshold <- as.factor(abs(data[,grep("logFC",colnames(data))]) > FCThresh & data[,grep("padj",colnames(data))] < PThresh)
        data$logFC <- data[,grep("logFC",colnames(data))]
        #         res <- quantile(data$logFC, c(0.02, 0.98))
        #         a <- as.numeric(res[1])
        #         b <- as.numeric(res[2])
        #         data$logFC[data$logFC < a] <- a # truncate fold change values at 2nd and 98th percentiles
        #         data$logFC[data$logFC > b] <- b # for higher resolution in the bulk of the plot
        data$padj <- data[,grep("padj",colnames(data))]
        #         res2 <- quantile(data$padj, 0.02, na.rm=TRUE)
        #         c <- as.numeric(res2)
        #         data$padj[data$padj < c] <- c
        data$Comp <- c
        data <- data[,c("Comp","logFC","padj","Threshold")]
        return(data)
      })

      plotData <- do.call(rbind,data)
      plotData$Comp <- as.factor(plotData$Comp)

      p3 <- ggplot(plotData, aes(x=logFC, y=-log10(padj), colour=Threshold)) +
        geom_point(alpha=0.4, size=1.75) +
        facet_wrap(~Comp) +
        theme_bw() +
        theme(axis.line.x = element_line(colour = "black"),
              axis.line.y = element_line(colour="black"),
              plot.background=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              panel.border=element_blank())

      if(!is.null(x_lab)){
        p3 <- p3 + labs(x=x_lab)
      }else{
        p3 <- p3 + labs(x="Log2FC")
      }

      if(!is.null(y_lab)){
        p3 <- p3 + labs(y=y_lab)
      }else{
        p3 <- p3 + labs(y="-log10(padj)")
      }

      if(!is.null(plot_title)){
        p3 <- p3 + ggtitle(plot_title)
      }else{
        p3 <- p3 + ggtitle(paste(r," Test", sep=""))
      }

      if(!is.null(leglab)){
        p3 <- p3 + guides(colour=guide_legend(title=leglab))
      }

      print(p3)
    }))
  }

}


#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param x_axis Required, a character vector specifying which variable to put on the x-axis, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param color Optional, a character vector specifying which variable to map to colors, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param shape Optional, a character vector specifying which variable to map to shape, must be one of the column names in attr(results_object, "group_DF"). Default is NULL.
#'@param scales Optional, a character vector to describe if any of the axes should be free when faceting. Default is "free_y" for a free y-axis.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.alphaRes <- function(results_object, x_axis="Group", color="Group", shape=NULL, plot_title=NULL, scales="free_y",
                          x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.alphaRes(results_object, x_axis, color, shape, plot_title, scales, x_lab, y_lab, leglab, ...)
}

.plot.alphaRes <- function(results_object, x_axis="Group", color=NULL, shape=NULL, plot_title=NULL, scales="free_y",
                           x_lab=NULL, y_lab=NULL, leglab=NULL) {

  library(ggplot2)

  if(x_axis %in% c("Group","Groups","group","groups","G","g")){
    x_axis <- "Group"
  }else{
    if(!(x_axis %in% names(attr(results_object, "group_DF")))){
      stop("x_axis must be one of the columns in group_DF")
    }
  }

  if(!is.null(color)){
    if(!(color %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("color must be one of the columns in group_DF")
    }
  }

  if(!is.null(shape)){
    if(!(shape %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("shape must be one of the columns in group_DF")
    }
  }

  if(!is.null(scales)){
    if(!(scales %in% c("free_y","free_x","free","fixed"))){
      stop("scales must be one of 'free_y', 'free_x', 'free', or 'fixed'.")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }


  results_object2 <- cbind(rownames(results_object), results_object)
  plot.data <- melt(results_object2)
  names(plot.data)[1] <- "Test"
  names(plot.data)[2] <- attr(results_object, "cnames")$fdata_cname

  plot.data <- merge(plot.data, attr(results_object, "group_DF"), by=attr(results_object, "cnames")$fdata_cname)

  map <- aes_string(x=x_axis, y="value", colour=color, shape=shape)
  p <- ggplot(plot.data, map) +
    geom_jitter(size=3, width=0.1, height=0) +
    facet_wrap(~Test, scales=scales) +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(x_lab)){
    p <- p + labs(x=x_lab)
  }else{
    p <- p + labs(x="Group")
  }

  if(!is.null(y_lab)){
    p <- p + labs(y=y_lab)
  }else{
    p <- p + labs(y="Alpha Diversity Measure")
  }

  if(!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }else{
    p <- p + ggtitle("Alpha Diversity")
  }

  if(!is.null(leglab)){
    p <- p + guides(colour=guide_legend(title=leglab))
  }

  return(p)
}




#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param x_axis Required, a character vector specifying which variable to put on the x-axis, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param color Optional, a character vector specifying which variable to map to colors, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param shape Optional, a character vector specifying which variable to map to shape, must be one of the column names in attr(results_object, "group_DF"). Default is NULL.
#'@param scales Optional, a character vector to describe if any of the axes should be free when faceting. Default is "free_y" for a free y-axis.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.evenRes <- function(results_object, x_axis="Group", color="Group", shape=NULL, plot_title=NULL, scales="free_y",
                         x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.evenRes(results_object, x_axis, color, shape, plot_title, scales, x_lab, y_lab, leglab, ...)
}

.plot.evenRes <- function(results_object, x_axis="Group", color=NULL, shape=NULL, plot_title=NULL, scales="free_y",
                          x_lab=NULL, y_lab=NULL, leglab=NULL) {

  library(ggplot2)

  if(x_axis %in% c("Group","Groups","group","groups","G","g")){
    x_axis <- "Group"
  }else{
    if(!(x_axis %in% names(attr(results_object, "group_DF")))){
      stop("x_axis must be one of the columns in group_DF")
    }
  }

  if(!is.null(color)){
    if(!(color %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("color must be one of the columns in group_DF")
    }
  }

  if(!is.null(shape)){
    if(!(shape %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("shape must be one of the columns in group_DF")
    }
  }

  if(!is.null(scales)){
    if(!(scales %in% c("free_y","free_x","free","fixed"))){
      stop("scales must be one of 'free_y', 'free_x', 'free', or 'fixed'.")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  results_object2 <- cbind(rownames(results_object), results_object)
  plot.data <- melt(results_object2)
  names(plot.data)[1] <- "Test"
  names(plot.data)[2] <- attr(results_object, "cnames")$fdata_cname

  plot.data <- merge(plot.data, attr(results_object, "group_DF"), by=attr(results_object, "cnames")$fdata_cname)

  map <- aes_string(x=x_axis, y="value", colour=color, shape=shape)
  p <- ggplot(plot.data, map) +
    geom_jitter(size=3, width=0.1, height=0) +
    facet_wrap(~Test, scales=scales) +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(x_lab)){
    p <- p + labs(x=x_lab)
  }else{
    p <- p + labs(x="Group")
  }

  if(!is.null(y_lab)){
    p <- p + labs(y=y_lab)
  }else{
    p <- p + labs(y="Evenness Measure")
  }

  if(!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }else{
    p <- p + ggtitle("Evenness")
  }

  if(!is.null(leglab)){
    p <- p + guides(colour=guide_legend(title=leglab))
  }

  return(p)
}




#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param variable Required, character vector specifying Which Jaccard variable to plot - options are "Median", "InterQuartileRange", "Average", and "StdDev". Default is "Median".
#'@param x_axis Required, a character vector specifying which variable to put on the x-axis, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param color Optional, a character vector specifying which variable to map to colors, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param shape Optional, a character vector specifying which variable to map to shape, must be one of the column names in attr(results_object, "group_DF"). Default is NULL.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.jaccardRes <- function(results_object, variable="Median", x_axis="Group", color="Group", shape=NULL, plot_title=NULL,
                            x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.jaccardRes(results_object, variable, x_axis, color, shape, plot_title, x_lab, y_lab, leglab, ...)
}

.plot.jaccardRes <- function(results_object, variable="Median", x_axis="Group", color="Group", shape=NULL, plot_title=NULL,
                             x_lab=NULL, y_lab=NULL, leglab=NULL) {

  library(ggplot2)

  ## initial checks ##
  if(x_axis %in% c("Group","Groups","group","groups","G","g")){
    x_axis <- "Group"
  }else{
    if(!(x_axis %in% names(attr(results_object, "group_DF")))){
      stop("x_axis must be one of the columns in group_DF")
    }
  }

  if(!is.null(color)){
    if(!(color %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("color must be one of the columns in group_DF")
    }
  }

  if(!is.null(shape)){
    if(!(shape %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("shape must be one of the columns in group_DF")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  if(length(variable) > 1){
    stop("Can only use one parameter for variable")
  }

  if(!(variable %in% c("Median","InterQuartileRange","Average","StdDev"))){
    stop("variable must be one of 'Median', 'InterQuartileRange', 'Average', or 'StdDev'")
  }
  ## end initial checks ##

  plot.data <- melt(results_object)
  names(plot.data)[1] <- "Grp"
  names(plot.data)[2] <- attr(results_object, "cnames")$fdata_cname
  plot.data <- plot.data[which(plot.data$variable == variable),]

  plot.data <- merge(plot.data, attr(results_object, "group_DF"), by=attr(results_object, "cnames")$fdata_cname)

  map <- aes_string(x=x_axis, y="value", colour=color, shape=shape)
  p <- ggplot(plot.data, map) +
    geom_jitter(size=3, width=0.1, height=0) +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(x_lab)){
    p <- p + labs(x=x_lab)
  }else{
    p <- p + labs(x="Group")
  }

  if(!is.null(y_lab)){
    p <- p + labs(y=y_lab)
  }else{
    if(attr(results_object, "similarity")){
      p <- p + labs(y=paste(variable," Jaccard Similarity", sep=""))
    }else{
      p <- p + labs(y=paste(variable," Jaccard Dissimilarity", sep=""))
    }
  }

  if(!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }else{
    p <- p + ggtitle("Jaccard Index")
  }

  if(!is.null(leglab)){
    p <- p + guides(colour=guide_legend(title=leglab))
  }

  return(p)
}




#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param breaks Required, a number specifying the number of breaks to have in the cumulative graph. Default is 100.
#'@param max_count Optional, a number specifying the maximum count number to show on the graph. Default is NULL.
#'@param min_num Optional, a number specifying the desired cut point in order to visualize how many OTUs would be lost if that cut point was used. sum_based_filter uses strictly less than the desired min_num, so this will show the valus for strictly less than the desired min_num. Default is NULL.
#'@param min_samp Optional, for k/a filtering, a number specifying that OTUs must be seen in at least this many samples. Default is 2 for ka filters and NULL for everything else.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.countFilter <- function(results_object, breaks=100, max_count=NULL, min_num=NULL, min_samp=NULL, plot_title=NULL, x_lab=NULL, y_lab=NULL, ...) {
  .plot.countFilter(results_object, breaks, max_count, min_num, min_samp, plot_title, x_lab, y_lab, ...)
}

.plot.countFilter <- function(results_object, breaks=100, max_count=NULL, min_num=NULL, min_samp=NULL, plot_title=NULL, x_lab=NULL, y_lab=NULL) {

  library(ggplot2)

  ## initial checks ##
  if(!is.null(min_num)){
    # check that min_num is numeric >= 0 #
    if(!(class(min_num) %in% c("numeric","integer")) || min_num < 0){
      stop("min_num must be an integer >= 0")
    }
  }

  # Check if min_samp is provided if fn="ka"
  if(is.null(min_samp)){
    if(attr(results_object, "function") == "ka"){
      warning("No minimum sample size specified, defaulting to 2")
      min_samp = 2
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  ## end of initial checks ##

  fn <- attr(results_object, "function")

  # format k/a data
  if(fn == "ka"){
    results_object <- results_object[,c(1,grep(paste("NumSamples_",min_samp,sep=""), colnames(results_object)))]
    colnames(results_object)[2] <- "kaOTUs"
  }

  # limit data, no need to look at all of it #
  if(!is.null(max_count)){
    if(fn == "percent"){
      max_count <- max_count / 100
    }
    results_object <- results_object[which(results_object[,paste(fn,"OTUs",sep="")] < max_count),]
  }else{
    if(fn != "ka"){
      results_object <- results_object[which(results_object[,paste(fn,"OTUs",sep="")] < attr(results_object, "threshold")),]
    }else{
      results_object <- results_object[which(results_object[,paste(fn,"OTUs",sep="")] < quantile(results_object[,2],0.95)),]
    }
  }

  # get abundance counts #
  brkpt <- max(results_object[,paste(fn,"OTUs",sep="")])/breaks

  break_data <- seq(from=0, to=ifelse(is.null(max_count),max(results_object[,paste(fn,"OTUs",sep="")]),max_count), by=brkpt)

  all_counts <- lapply(break_data, function(x){
    data.frame(point=x, sumcount=length(which(results_object[,paste(fn,"OTUs",sep="")] <= x)))
  })
  all_counts <- do.call(rbind, all_counts)

  fn_lab <- paste(toupper(substring(fn, 1, 1)),substring(fn, 2), sep="", collapse=" ")

  # make labels #
  if(fn == "ka"){
    xlabel <- ifelse(is.null(x_lab), paste("Max Count of Biomolecule in each of ", min_samp, " Samples", sep=""), x_lab)
    ylabel <- ifelse(is.null(y_lab), "Cumulative Frequency", y_lab)
    plot_title <- ifelse(is.null(plot_title), paste("Cumulative Frequency of OTUs seen in ", min_samp, " Samples", sep=""), plot_title)
  }else{
    xlabel <- ifelse(is.null(x_lab), paste(fn_lab," Count of Biomolecule Across All Samples", sep=""), x_lab)
    ylabel <- ifelse(is.null(y_lab), "Cumulative Frequency", y_lab)
    plot_title <- ifelse(is.null(plot_title), paste("Cumulative Frequency of ", fn_lab, " of Biomolecules in Samples",sep=""), plot_title)
  }

  p <- ggplot(all_counts) +
    geom_rect(aes(xmin=point-brkpt/2, xmax=point+brkpt/2,
                  ymin=0, ymax=sumcount), fill="royalblue1", col="black") +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(plot_title) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(min_num)) {
    # mark on graph where min_num is #
    if(fn == "percent"){
      min_num = min_num / 100
    }
    num_tested <- length(which(results_object[,paste(fn,"OTUs",sep="")] <= min_num))
    p <- p + ggplot2::annotate(geom = "rect", xmin = min_num - brkpt/2, xmax = min_num + brkpt/2,
                               ymin = 0, ymax = num_tested, fill="black", col="black", size=1.5)
  }

  return(p)

}

#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param breaks Required, a number specifying the number of breaks to have in
#'  the cumulative graph. Default is 100.
#'@param max_count Optional, a number specifying the maximum count number to
#'  show on the graph. Default is NULL.
#'@param min_num Optional, a number specifying the desired cut point in order to
#'  visualize how many Sampless would be lost if that cut point was used.
#'  sample_based_filter uses strictly less than the desired min_num, so this
#'  will show the valus for strictly less than the desired min_num. Default is
#'  NULL.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.sampleFilter <- function(results_object, breaks=100, max_count=NULL, min_num=NULL, plot_title=NULL, x_lab=NULL, y_lab=NULL, ...) {
  .plot.sampleFilter(results_object, breaks, max_count, min_num, plot_title, x_lab, y_lab, ...)
}

.plot.sampleFilter <- function(results_object, breaks=100, max_count=NULL, min_num=NULL, plot_title=NULL, x_lab=NULL, y_lab=NULL) {

  library(ggplot2)

  ## initial checks ##
  if(!is.null(min_num)){
    # check that min_num is numeric >= 0 #
    if(!(class(min_num) %in% c("numeric","integer")) || min_num < 0){
      stop("min_num must be an integer >= 0")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  ## end of initial checks ##

  fn <- attr(results_object, "function")

  # limit data, no need to look at all of it #
  if(!is.null(max_count)){
    results_object <- results_object[which(results_object[,paste(fn,"Samps",sep="")] < max_count),]
  }else{
    results_object <- results_object[which(results_object[,paste(fn,"Samps",sep="")] < attr(results_object, "threshold")),]
  }

  # get abundance counts #
  brkpt <- max(results_object[,paste(fn,"Samps",sep="")])/breaks

  break_data <- seq(from=0, to=ifelse(is.null(max_count),max(results_object[,paste(fn,"Samps",sep="")]),max_count), by=brkpt)

  all_counts <- lapply(break_data, function(x){
    data.frame(point=x, sumcount=length(which(results_object[,paste(fn,"Samps",sep="")] <= x)))
  })
  all_counts <- do.call(rbind, all_counts)

  fn_lab <- paste(toupper(substring(fn, 1, 1)),substring(fn, 2), sep="", collapse=" ")
  # make labels #
  xlabel <- ifelse(is.null(x_lab), paste(fn_lab," Count of Samples Across All Biomolecules", sep=""), x_lab)
  ylabel <- ifelse(is.null(y_lab), "Cumulative Frequency", y_lab)
  plot_title <- ifelse(is.null(plot_title), paste("Cumulative Frequency of ", fn_lab, " of Sample Biomolecules",sep=""), plot_title)

  p <- ggplot(all_counts) +
    geom_rect(aes(xmin=point-brkpt/2, xmax=point+brkpt/2,
                  ymin=0, ymax=sumcount), fill="royalblue1", col="black") +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(plot_title) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(min_num)) {
    # mark on graph where min_num is #
    num_tested <- length(which(results_object[,paste(fn,"Samps",sep="")] <= min_num))
    p <- p + ggplot2::annotate(geom = "rect", xmin = min_num - brkpt/2, xmax = min_num + brkpt/2,
                               ymin = 0, ymax = num_tested, fill="black", col="black", size=1.5)
  }

  return(p)

}



#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param x_axis Required, a character vector specifying which variable to group data by and put on the x-axis, must be one of the column names in results_object$e_data, results_object$e_meta, or attr(results_object, "group_DF"). Default is "Group".
#'@param class Required, a character vector specifying which variable to group data by (e.g. "Phylum") and to color data by, must be one of the column names in results_object$e_data, results_object$e_meta, or attr(results_object, "group_DF"). Default is "Phylum".
#'@param grp_fn Required, a character vector specifying which function to use to summarise values across replicates. Can use one of "median", "mean", or "sum". Default is "median".
#'@param subset_by Optional, which variable to subset the data by. NULL will not subset the data. Default is NULL.
#'@param subset_val Optional, which values to include in the subset.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.seqData <- function(results_object, x_axis="Group", class="Phylum", grp_fn="median", subset_by=NULL, subset_val=NULL,
                          plot_title=NULL, x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.seqData(results_object, x_axis, class, grp_fn, subset_by, subset_val, plot_title, x_lab, y_lab, leglab, ...)
}

.plot.seqData <- function(results_object, x_axis="Group", class="Phylum", grp_fn="median", subset_by=NULL, subset_val=NULL, plot_title=NULL,
                           x_lab=NULL, y_lab=NULL, leglab=NULL) {

  # normalized or not-normalized counts???
  # remove 0/NA from data, so med(0,0,4,0,0) is 4, not 0???

  library(reshape2)
  library(ggplot2)
  library(dplyr)

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  if(!is.null(leglab)){
    if(!is.character(leglab)){
      stop("leglab must be a character vector")
    }
  }

  if(!is.null(subset_by)){
    if(!is.character(subset_by) | !(subset_by %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"groupDF"))))){
      stop("subset_by must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
    }
  }

  if(!is.character(x_axis) | !(x_axis %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"group_DF"))))){
    stop("x_axis must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
  }

  if(!is.character(class) | !(class %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"groupDF"))))){
    stop("class must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
  }

  if(!is.null(subset_val)){
    if(!is.character(subset_val)){
      stop("subset_val must be a character vector")
    }
  }

  ## Check group designation ##
  if(is.null(attr(results_object,"group_DF"))){
    stop("Run group designation function first")
  }

  ## Check subsetting ##
  if(!is.null(subset_by) & is.null(subset_val)){
    stop("Must specify values to subset to when specifying variable to subset with")
  }else if(is.null(subset_by) & !is.null(subset_val)){
    stop("Must specify which variable to use for subsetting when specifying values for subset")
  }

  ## Merge e_data with e_meta to get feature meta data ##
  if(!is.null(results_object$e_meta)){
    data <- merge(results_object$e_data, results_object$e_meta, by=attr(results_object,"cnames")$edata_cname)
  }else{
    #data <- results_object$e_data
    stop("Data must have an e_meta mapping")
  }

  ## Subset ##
  if(!is.null(subset_by) & !is.null(subset_val)){
    data <- subset(data, get(subset_by) %in% subset_val)
  }

  ## Format data for plotting ##
  data_melt <- melt(data)
  colnames(data_melt)[which(colnames(data_melt)=="variable")] <- attr(results_object, "cnames")$fdata_cname

  ## Merge data with group_DF to get groupings ##
  data_melt <-  merge(data_melt, attr(results_object,"group_DF"), by=attr(results_object, "cnames")$fdata_cname)

  vars <- c(x_axis, class, attr(results_object, "cnames")$fdata_cname)
  vars <- lapply(vars, as.symbol)

  data_grp1 <- data_melt %>% dplyr::group_by_(.dots=vars) %>% dplyr::summarise(sum=sum(value, na.rm=TRUE))
  data_grp1 <- data_grp1[-which(data_grp1$sum==0),]

  vars1 <- c(x_axis, class)
  vars1 <- lapply(vars1, as.symbol)

  if(grp_fn %in% c("median","med")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=median(sum, na.rm=TRUE))
  }else if(grp_fn %in% c("mean","average","avg")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=mean(sum, na.rm=TRUE))
  }else if(grp_fn %in% c("sum")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=sum(sum, na.rm=TRUE))
  }

  data_grp <- data_grp[order(-data_grp$value),]

  TaxonBarPallette <- c("#FF00BB","#CC00FF","#F2BFFF","#7A0099","#0022FF","#8091FF","#001499","#00F2FF","#CCFCFF","#009199","#00D90E","#BFFFC4","#007308","#FFFF00","#DDFF00","#B3B300","#FF9100","#FFC880","#995700","#FF0000","#FFABAB","#990000","#BFBFBF","#636363","#000000")

  ## Bar plot ##
  map <- aes_string(x=x_axis, y="value", fill=class)
  p <- ggplot(data_grp, map) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_manual(values=rep(TaxonBarPallette,5))+
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(x_lab)){
    p <- p + labs(x=x_lab)
  }else{
    p <- p + labs(x=x_axis)
  }

  if(!is.null(y_lab)){
    p <- p + labs(y=y_lab)
  }else{
    p <- p + labs(y="Abundance")
  }

  if(!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }else{
    p <- p + ggtitle("Taxonomy")
  }

  if(!is.null(leglab)){
    p <- p + guides(fill=guide_legend(title=leglab))
  }

  return(p)

}



#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param x_axis Required, a character vector specifying which variable to group data by and put on the x-axis, must be one of the column names in results_object$e_data, results_object$e_meta, or attr(results_object, "group_DF"). Default is "Group".
#'@param class Required, a character vector specifying which variable to group and color data by, must be one of the column names in results_object$e_data, results_object$e_meta, or attr(results_object, "group_DF"). Default is "ECNum".
#'@param grp_fn Required, a character vector specifying which function to use to summarise values across replicates. Can use one of "median", "mean", or "sum". Default is "median".
#'@param subset_by Optional, which variable to subset the data by. NULL will not subset the data. Default is NULL.
#'@param subset_val Optional, which values to include in the subset.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.cDNAdata <- function(results_object, x_axis="Group", class="ECNum", grp_fn="median", subset_by=NULL, subset_val=NULL,
                          plot_title=NULL, x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.cDNAdata(results_object, x_axis, class, grp_fn, subset_by, subset_val, plot_title, x_lab, y_lab, leglab, ...)
}

.plot.cDNAdata <- function(results_object, x_axis="Group", class="ECNum", grp_fn="median", subset_by=NULL, subset_val=NULL, plot_title=NULL,
                           x_lab=NULL, y_lab=NULL, leglab=NULL) {

  # normalized or not-normalized counts???
  # remove 0/NA from data, so med(0,0,4,0,0) is 4, not 0???

  library(reshape2)
  library(ggplot2)
  library(dplyr)

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  if(!is.null(leglab)){
    if(!is.character(leglab)){
      stop("leglab must be a character vector")
    }
  }

  if(!is.null(subset_by)){
    if(!is.character(subset_by) | !(subset_by %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"groupDF"))))){
      stop("subset_by must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
    }
  }

  if(!is.character(x_axis) | !(x_axis %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"group_DF"))))){
    stop("x_axis must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
  }

  if(!is.character(class) | !(class %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"groupDF"))))){
    stop("class must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
  }

  if(!is.null(subset_val)){
    if(!is.character(subset_val)){
      stop("subset_val must be a character vector")
    }
  }

  ## Check group designation ##
  if(is.null(attr(results_object,"group_DF"))){
    stop("Run group designation function first")
  }

  ## Check subsetting ##
  if(!is.null(subset_by) & is.null(subset_val)){
    stop("Must specify values to subset to when specifying variable to subset with")
  }else if(is.null(subset_by) & !is.null(subset_val)){
    stop("Must specify which variable to use for subsetting when specifying values for subset")
  }

  ## Merge e_data with e_meta to get feature meta data ##
  if(!is.null(results_object$e_meta)){
    data <- merge(results_object$e_data, results_object$e_meta, by=attr(results_object,"cnames")$edata_cname)
  }else{
    #data <- results_object$e_data
    stop("Data must have an e_meta mapping")
  }

  ## Subset ##
  if(!is.null(subset_by) & !is.null(subset_val)){
    data <- subset(data, get(subset_by) %in% subset_val)
  }

  ## Format data for plotting ##
  data_melt <- melt(data)
  colnames(data_melt)[which(colnames(data_melt)=="variable")] <- attr(results_object, "cnames")$fdata_cname

  ## Merge data with group_DF to get groupings ##
  data_melt <-  merge(data_melt, attr(results_object,"group_DF"), by=attr(results_object, "cnames")$fdata_cname)

  vars <- c(x_axis, class, attr(results_object, "cnames")$fdata_cname)
  vars <- lapply(vars, as.symbol)

  data_grp1 <- data_melt %>% dplyr::group_by_(.dots=vars) %>% dplyr::summarise(sum=sum(value, na.rm=TRUE))
  data_grp1 <- data_grp1[-which(data_grp1$sum==0),]

  vars1 <- c(x_axis, class)
  vars1 <- lapply(vars1, as.symbol)

  if(grp_fn %in% c("median","med")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=median(sum, na.rm=TRUE))
  }else if(grp_fn %in% c("mean","average","avg")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=mean(sum, na.rm=TRUE))
  }else if(grp_fn %in% c("sum")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=sum(sum, na.rm=TRUE))
  }

  data_grp <- data_grp[order(-data_grp$value),]

  ## Bar plot ##
  map <- aes_string(x=x_axis, y="value", fill=class)
  p <- ggplot(data_grp, map) +
    geom_bar(stat="identity", position="stack") +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(x_lab)){
    p <- p + labs(x=x_lab)
  }else{
    p <- p + labs(x=x_axis)
  }

  if(!is.null(y_lab)){
    p <- p + labs(y=y_lab)
  }else{
    p <- p + labs(y="Abundance")
  }

  if(!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }else{
    p <- p + ggtitle("cDNAdata")
  }

  if(!is.null(leglab)){
    p <- p + guides(fill=guide_legend(title=leglab))
  }

  return(p)

}


#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param x_axis Required, a character vector specifying which variable to group data by and put on the x-axis, must be one of the column names in results_object$e_data, results_object$e_meta, or attr(results_object, "group_DF"). Default is "Group".
#'@param class Required, a character vector specifying which variable to group and color data by, must be one of the column names in results_object$e_data, results_object$e_meta, or attr(results_object, "group_DF"). Default is "ECNum".
#'@param grp_fn Required, a character vector specifying which function to use to summarise values across replicates. Can use one of "median", "mean", or "sum". Default is "median".
#'@param subset_by Optional, which variable to subset the data by. NULL will not subset the data. Default is NULL.
#'@param subset_val Optional, which values to include in the subset.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.gDNAdata <- function(results_object, x_axis="Group", class="ECNum", grp_fn="median", subset_by=NULL, subset_val=NULL,
                          plot_title=NULL, x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.gDNAdata(results_object, x_axis, class, grp_fn, subset_by, subset_val, plot_title, x_lab, y_lab, leglab, ...)
}

.plot.gDNAdata <- function(results_object, x_axis="Group", class="ECNum", grp_fn="median", subset_by=NULL, subset_val=NULL, plot_title=NULL,
                           x_lab=NULL, y_lab=NULL, leglab=NULL) {

  # normalized or not-normalized counts???
  # remove 0/NA from data, so med(0,0,4,0,0) is 4, not 0???

  library(reshape2)
  library(ggplot2)
  library(dplyr)

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  if(!is.null(leglab)){
    if(!is.character(leglab)){
      stop("leglab must be a character vector")
    }
  }

  if(!is.null(subset_by)){
    if(!is.character(subset_by) | !(subset_by %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"groupDF"))))){
      stop("subset_by must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
    }
  }

  if(!is.character(x_axis) | !(x_axis %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"group_DF"))))){
    stop("x_axis must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
  }

  if(!is.character(class) | !(class %in% c(colnames(results_object$e_data), colnames(results_object$e_meta), colnames(attr(results_object,"groupDF"))))){
    stop("class must be a character vector and a column name found in results_object$e_data, results_object$e_meta, or group_DF")
  }

  if(!is.null(subset_val)){
    if(!is.character(subset_val)){
      stop("subset_val must be a character vector")
    }
  }

  ## Check group designation ##
  if(is.null(attr(results_object,"group_DF"))){
    stop("Run group designation function first")
  }

  ## Check subsetting ##
  if(!is.null(subset_by) & is.null(subset_val)){
    stop("Must specify values to subset to when specifying variable to subset with")
  }else if(is.null(subset_by) & !is.null(subset_val)){
    stop("Must specify which variable to use for subsetting when specifying values for subset")
  }

  ## Merge e_data with e_meta to get feature meta data ##
  if(!is.null(results_object$e_meta)){
    data <- merge(results_object$e_data, results_object$e_meta, by=attr(results_object,"cnames")$edata_cname)
  }else{
    #data <- results_object$e_data
    stop("Data must have an e_meta mapping")
  }

  ## Subset ##
  if(!is.null(subset_by) & !is.null(subset_val)){
    data <- subset(data, get(subset_by) %in% subset_val)
  }

  ## Format data for plotting ##
  data_melt <- melt(data)
  colnames(data_melt)[which(colnames(data_melt)=="variable")] <- attr(results_object, "cnames")$fdata_cname

  ## Merge data with group_DF to get groupings ##
  data_melt <-  merge(data_melt, attr(results_object,"group_DF"), by=attr(results_object, "cnames")$fdata_cname)

  vars <- c(x_axis, class, attr(results_object, "cnames")$fdata_cname)
  vars <- lapply(vars, as.symbol)

  data_grp1 <- data_melt %>% dplyr::group_by_(.dots=vars) %>% dplyr::summarise(sum=sum(value, na.rm=TRUE))
  data_grp1 <- data_grp1[-which(data_grp1$sum==0),]

  vars1 <- c(x_axis, class)
  vars1 <- lapply(vars1, as.symbol)

  if(grp_fn %in% c("median","med")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=median(sum, na.rm=TRUE))
  }else if(grp_fn %in% c("mean","average","avg")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=mean(sum, na.rm=TRUE))
  }else if(grp_fn %in% c("sum")){
    data_grp <- data_grp1 %>% dplyr::group_by_(.dots=vars1) %>% dplyr::summarise(value=sum(sum, na.rm=TRUE))
  }

  data_grp <- data_grp[order(-data_grp$value),]

  ## Bar plot ##
  map <- aes_string(x=x_axis, y="value", fill=class)
  p <- ggplot(data_grp, map) +
    geom_bar(stat="identity", position="stack") +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(x_lab)){
    p <- p + labs(x=x_lab)
  }else{
    p <- p + labs(x=x_axis)
  }

  if(!is.null(y_lab)){
    p <- p + labs(y=y_lab)
  }else{
    p <- p + labs(y="Abundance")
  }

  if(!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }else{
    p <- p + ggtitle("gDNAdata")
  }

  if(!is.null(leglab)){
    p <- p + guides(fill=guide_legend(title=leglab))
  }

  return(p)

}


#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param abun Optional, an object of class 'abundance', used for plotting richness vs abundance
#'@param x_axis Required, a character vector specifying which variable to put on the x-axis, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param color Optional, a character vector specifying which variable to map to colors, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param scales Optional, a character vector specifying if any/both of the axes be free. Default is "free_y", for a free y-axis.
#'@param shape Optional, a character vector specifying which variable to map to shape, must be one of the column names in attr(results_object, "group_DF"). Default is NULL.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.richRes <- function(results_object, abun=NULL, x_axis="Group", color="Group", scales="free_y", shape=NULL, plot_title=NULL,
                         x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.richRes(results_object, abun, x_axis, color, scales, shape, plot_title, x_lab, y_lab, leglab, ...)
}

.plot.richRes <- function(results_object, abun=NULL, x_axis="Group", color="Group", scales="free_y", shape=NULL, plot_title=NULL,
                          x_lab=NULL, y_lab=NULL, leglab=NULL) {

  library(ggplot2)

  if(x_axis %in% c("Group","Groups","group","groups","G","g")){
    x_axis <- "Group"
  }else{
    if(!(x_axis %in% names(attr(results_object, "group_DF")))){
      stop("x_axis must be one of the columns in group_DF")
    }
  }

  if(!is.null(color)){
    if(!(color %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("color must be one of the columns in group_DF")
    }
  }

  if(!is.null(shape)){
    if(!(shape %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("shape must be one of the columns in group_DF")
    }
  }

  if(!is.null(scales)){
    if(!(scales %in% c("free_y","free_x","free","fixed"))){
      stop("scales must be one of 'free_y', 'free_x', 'free', or 'fixed'.")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }

  if(is.null(abun)){
    results_object2 <- cbind(rownames(results_object), results_object)
    plot.data <- melt(results_object2)
    names(plot.data)[1] <- "Test"
    names(plot.data)[2] <- attr(results_object, "cnames")$fdata_cname

    plot.data <- merge(plot.data, attr(results_object, "group_DF"), by=attr(results_object, "cnames")$fdata_cname)

    map <- aes_string(x=x_axis, y="value", colour=color, shape=shape)
    p <- ggplot(plot.data, map) +
      geom_jitter(size=3, width=0.1, height=0) +
      facet_wrap(~Test, scales=scales) +
      theme_bw() +
      theme(axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour="black"),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank())

    if(!is.null(x_lab)){
      p <- p + labs(x=x_lab)
    }else{
      p <- p + labs(x="Group")
    }

    if(!is.null(y_lab)){
      p <- p + labs(y=y_lab)
    }else{
      p <- p + labs(y="Richness")
    }

    if(!is.null(plot_title)){
      p <- p + ggtitle(plot_title)
    }else{
      p <- p + ggtitle("Richness")
    }

    if(!is.null(leglab)){
      p <- p + guides(colour=guide_legend(title=leglab))
    }
  }else{
    if(!("abunRes" %in% class(abun))){stop("abun must be of class 'abunRes'")}
    rich <- data.frame(Samples=colnames(results_object), Richness=as.data.frame(t(results_object))$observed)
    abun <- data.frame(Samples=rownames(abun), Abundance=abun$abundance)

    data <- merge(rich, abun, by="Samples")
    data <- merge(data, attr(results_object, "group_DF"), by.x="Samples", by.y=attr(results_object, "cnames")$fdata_cname)

    map <- aes_string(x="Richness", y="Abundance", colour=color, shape=shape)
    p <- ggplot(data, map) +
      geom_jitter(size=3, width=0.1, height=0) +
      theme_bw() +
      theme(axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour="black"),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank())

    if(!is.null(x_lab)){
      p <- p + labs(x=x_lab)
    }else{
      p <- p + labs(x="Richness")
    }

    if(!is.null(y_lab)){
      p <- p + labs(y=y_lab)
    }else{
      p <- p + labs(y="Abundance")
    }

    if(!is.null(plot_title)){
      p <- p + ggtitle(plot_title)
    }else{
      p <- p + ggtitle("Richness vs Abundance")
    }

    if(!is.null(leglab)){
      p <- p + guides(colour=guide_legend(title=leglab))
    }
  }

  return(p)
}



#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param rich Optional, an object of class 'richness', used for plotting richness vs abundance
#'@param x_axis Required, a character vector specifying which variable to put on the x-axis, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param color Optional, a character vector specifying which variable to map to colors, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param shape Optional, a character vector specifying which variable to map to shape, must be one of the column names in attr(results_object, "group_DF"). Default is NULL.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.abunRes <- function(results_object, rich=NULL, x_axis="Group", color="Group", shape=NULL, plot_title=NULL,
                         x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.abunRes(results_object, rich, x_axis, color, shape, plot_title, x_lab, y_lab, leglab, ...)
}

.plot.abunRes <- function(results_object, rich=NULL, x_axis="Group", color=NULL, shape=NULL, plot_title=NULL,
                          x_lab=NULL, y_lab=NULL, leglab=NULL) {

  library(ggplot2)

  if(x_axis %in% c("Group","Groups","group","groups","G","g")){
    x_axis <- "Group"
  }else{
    if(!(x_axis %in% names(attr(results_object, "group_DF")))){
      stop("x_axis must be one of the columns in group_DF")
    }
  }

  if(!is.null(color)){
    if(!(color %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("color must be one of the columns in group_DF")
    }
  }

  if(!is.null(shape)){
    if(!(shape %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("shape must be one of the columns in group_DF")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }


  if(is.null(rich)){
    results_object2 <- cbind(rownames(results_object), results_object)
    plot.data <- melt(results_object2)
    names(plot.data)[2] <- "Abundance"
    names(plot.data)[1] <- attr(results_object, "cnames")$fdata_cname

    plot.data <- merge(plot.data, attr(results_object, "group_DF"), by=attr(results_object, "cnames")$fdata_cname)

    map <- aes_string(x=x_axis, y="value", colour=color, shape=shape)
    p <- ggplot(plot.data, map) +
      geom_jitter(size=3, width=0.1, height=0) +
      theme_bw() +
      theme(axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour="black"),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank())

    if(!is.null(x_lab)){
      p <- p + labs(x=x_lab)
    }else{
      p <- p + labs(x="Group")
    }

    if(!is.null(y_lab)){
      p <- p + labs(y=y_lab)
    }else{
      p <- p + labs(y="Abundance")
    }

    if(!is.null(plot_title)){
      p <- p + ggtitle(plot_title)
    }else{
      p <- p + ggtitle("Abundance")
    }

    if(!is.null(leglab)){
      p <- p + guides(colour=guide_legend(title=leglab))
    }
  }else{
    if(!("richRes" %in% class(rich))){stop("rich must be of class 'richRes'")}
    abun <- data.frame(Samples=rownames(results_object), Abundance=results_object$abundance)
    rich <- data.frame(Samples=colnames(rich), Richness=as.data.frame(t(rich))$observed)

    data <- merge(rich, abun, by="Samples")
    data <- merge(data, attr(results_object, "group_DF"), by.x="Samples", by.y=attr(results_object, "cnames")$fdata_cname)

    map <- aes_string(x="Richness", y="Abundance", colour=color, shape=shape)
    p <- ggplot(data, map) +
      geom_jitter(size=3, width=0.1, height=0) +
      theme_bw() +
      theme(axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour="black"),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank())

    if(!is.null(x_lab)){
      p <- p + labs(x=x_lab)
    }else{
      p <- p + labs(x="Richness")
    }

    if(!is.null(y_lab)){
      p <- p + labs(y=y_lab)
    }else{
      p <- p + labs(y="Abundance")
    }

    if(!is.null(plot_title)){
      p <- p + ggtitle(plot_title)
    }else{
      p <- p + ggtitle("Richness vs Abundance")
    }

    if(!is.null(leglab)){
      p <- p + guides(colour=guide_legend(title=leglab))
    }
  }

  return(p)
}



#'@export
#'@rdname plot_pmartRseq
#'@name plot_pmartRseq
#'@param x_axis Required, a character vector specifying which variable to put on the x-axis, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param color Optional, a character vector specifying which variable to map to colors, must be one of the column names in attr(results_object, "group_DF"). Default is "Group".
#'@param shape Optional, a character vector specifying which variable to map to shape, must be one of the column names in attr(results_object, "group_DF"). Default is NULL.
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.effspRes <- function(results_object, x_axis="Group", color="Group", shape=NULL, plot_title=NULL,
                          x_lab=NULL, y_lab=NULL, leglab=NULL, ...) {
  .plot.effspRes(results_object, x_axis, color, shape, plot_title, x_lab, y_lab, leglab, ...)
}

.plot.effspRes <- function(results_object, x_axis="Group", color=NULL, shape=NULL, plot_title=NULL,
                           x_lab=NULL, y_lab=NULL, leglab=NULL) {

  library(ggplot2)

  if(x_axis %in% c("Group","Groups","group","groups","G","g")){
    x_axis <- "Group"
  }else{
    if(!(x_axis %in% names(attr(results_object, "group_DF")))){
      stop("x_axis must be one of the columns in group_DF")
    }
  }

  if(!is.null(color)){
    if(!(color %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("color must be one of the columns in group_DF")
    }
  }

  if(!is.null(shape)){
    if(!(shape %in% c("Group","Groups","group","groups","G","g",names(attr(results_object, "group_DF"))))){
      stop("shape must be one of the columns in group_DF")
    }
  }

  if(!is.null(plot_title)){
    if(!is.character(plot_title)){
      stop("plot_title must be a character vector")
    }
  }

  if(!is.null(x_lab)){
    if(!is.character(x_lab)){
      stop("x_lab must be a character vector")
    }
  }

  if(!is.null(y_lab)){
    if(!is.character(y_lab)){
      stop("y_lab must be a character vector")
    }
  }


  results_object2 <- cbind(rownames(results_object), results_object)
  plot.data <- melt(results_object2)
  names(plot.data)[2] <- "EffSpecies"
  names(plot.data)[1] <- attr(results_object, "cnames")$fdata_cname

  plot.data <- merge(plot.data, attr(results_object, "group_DF"), by=attr(results_object, "cnames")$fdata_cname)

  map <- aes_string(x=x_axis, y="value", colour=color, shape=shape)
  p <- ggplot(plot.data, map) +
    geom_jitter(size=3, width=0.1, height=0) +
    theme_bw() +
    theme(axis.line.x = element_line(colour = "black"),
          axis.line.y = element_line(colour="black"),
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank())

  if(!is.null(x_lab)){
    p <- p + labs(x=x_lab)
  }else{
    p <- p + labs(x="Group")
  }

  if(!is.null(y_lab)){
    p <- p + labs(y=y_lab)
  }else{
    p <- p + labs(y="Effective Species")
  }

  if(!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }else{
    p <- p + ggtitle("Effective Species")
  }

  if(!is.null(leglab)){
    p <- p + guides(colour=guide_legend(title=leglab))
  }

  return(p)
}



#' @export
#' @rdname plot_pmartRseq
#' @name plot_pmartRseq
#' @param type Required, a character vector specifying which type of plot to create. "pvals" for a line plot showing the number of significantly expresed biomolecules at different p-value thresholds; or "flag" for a bar plot showing the number of indicated species in each group, at the previously specified threshold. The default is "pvals".
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
plot.indspRes <- function(results_object, type="pvals", ...){
  .plot.indspRes(results_object, type, ...)
}

.plot.indspRes <- function(results_object, type="pvals", x_lab=NULL, y_lab=NULL, plot_title=NULL, leglab=NULL){
  library(ggplot2)
  library(reshape2)
  library(gplots)
  library(RColorBrewer)

  ## initial checks ##
  if(!is.null(plot_title)){
    if(!is.character(plot_title)){ stop("plot_title must be a character vector")}
  }
  if(!is.null(x_lab)){
    if(!is.character(x_lab)){ stop("x_lab must be a character vector")}
  }
  if(!is.null(y_lab)){
    if(!is.character(y_lab)){ stop("y_lab must be a character vector")}
  }
  if(!is.null(type)){
    if(!is.character(type)){ stop("type must be a character vector")}
    if(any(!(type %in% c("pvals","flag")))){ stop("type must be at least one of 'pvals' or 'flag'")}
  }
  if(!("indspRes" %in% class(results_object))){
    stop("results_object must be of class indspRes")
  }
  ## end of initial checks ##

  ## p-value plot ##
  if("pvals" %in% tolower(type)){
    thresholds <- seq(0,0.4,0.05)

    data <- results_object[,grep("p.value", colnames(results_object))]

    Counts <-  lapply(thresholds, function(x) length(which(data < x)))
    Counts <- do.call(rbind, Counts)

    Counts <- data.frame(Threshold=thresholds, Counts)
    Counts$Threshold <- as.factor(Counts$Threshold)
    Counts <- melt(Counts)
    Counts$Threshold <- as.numeric(as.character(Counts$Threshold))
    names(Counts) <- gsub("Counts","value",names(Counts))

    p1 <- ggplot(Counts, aes(x=Threshold, y=value)) +
      geom_line(cex=1.5) +
      theme_bw() +
      theme(axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour="black"),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank())
    #labs(x="P-Value Threshold",y="Number of differentially expressed genes") +
    #ggtitle(paste(x," Test",sep=""))

    if(!is.null(x_lab)){
      p1 <- p1 + labs(x=x_lab)
    }else{
      p1 <- p1 + labs(x="P-Value Threshold")
    }

    if(!is.null(y_lab)){
      p1 <- p1 + labs(y=y_lab)
    }else{
      p1 <- p1 + labs(y="Number of indicator species")
    }

    if(!is.null(plot_title)){
      p1 <- p1 + ggtitle(plot_title)
    }else{
      p1 <- p1 + ggtitle("Indicator Species")
    }

    if(!is.null(leglab)){
      p1 <- p1 + guides(colour=guide_legend(title=leglab))
    }

    print(p1)
  }

  ## flags barplot ##
  if("flag" %in% tolower(type)){
    pal <- brewer.pal(9, "Set1")

    idx <- grep("Flag", colnames(results_object))

    numsig <- lapply(unique(attr(results_object, "group_DF")$Group), function(y){
      idx1 <- grep(y, colnames(results_object))
      idx2 <- idx[which(idx > idx1)[1]]

      sig <- length(which(results_object[,idx1] == 1 & results_object[,idx2] == 1))
      data.frame(Sample=y,NumIndSp=sig)
    })
    plotData <- do.call(rbind,numsig)

    # data <- results_object[,idx]
    # #data <- cbind(rownames(data),data)
    # #colnames(data)[1] <- attr(results_object, "cnames")$edata_cname
    # data <- melt(data)
    # if(length(idx) == 1){
    #   data$variable = colnames(results_object)[idx]
    #   data <- data[-which(data$value==0),]
    #   data$value <- gsub("1","Flags",data$value)
    # # data$Direction <- unlist(lapply(data$value, function(x) ifelse(x==1,"Up","Down")))
    # # data$Direction <- factor(data$Direction, levels=c("Up","Down"))
    # #
    #   plotData <- table(data$value)
    #   plotData <- melt(plotData)
    #   names(plotData) <- c("Flags","Count")
    #   plotData$Comparison <- "AllComparisons"
    # }else{
    #   data$Comparison <- unlist(lapply(as.character(data$variable), function(x) strsplit(x,"Flag_")[[1]][2]))
    #   data <- data[-which(data$value == 0),]
    #   data$value <- gsub("1","Flags",data$value)
    #   plotData <- table(data$Comparison, data$value)
    #   plotData <- melt(plotData)
    #   names(plotData) <- c("Comparison", "Flags", "Count")
    # }

    p2 <- ggplot(plotData, aes(x=Sample,y=NumIndSp,fill="NumFlags")) +
      geom_bar(stat="identity", position="dodge") +
      scale_fill_manual(values=pal[c(3,1)], name="NumIndSp") +
      theme_bw() +
      theme(axis.line.x = element_line(colour = "black"),
            axis.line.y = element_line(colour="black"),
            plot.background=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_blank()) +
      geom_text(aes(label=NumIndSp), position=position_dodge(width=0.9), vjust=1)
    #labs(x="Comparisons",y="Number of differentially expressed genes") +
    #ggtitle(paste(r, " Test",sep=""))

    if(!is.null(x_lab)){
      p2 <- p2 + labs(x=x_lab)
    }else{
      p2 <- p2 + labs(x="Comparisons")
    }

    if(!is.null(y_lab)){
      p2 <- p2 + labs(y=y_lab)
    }else{
      p2 <- p2 + labs(y="Number of indicator species")
    }

    if(!is.null(plot_title)){
      p2 <- p2 + ggtitle(plot_title)
    }else{
      p2 <- p2 + ggtitle("Indicator Species")
    }

    if(!is.null(leglab)){
      p2 <- p2 + guides(fill=guide_legend(title=leglab))
    }

    print(p2)
  }

}
