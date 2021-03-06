#' Plot all differential abundance results
#'
#' This function takes the differential abundance results of an omicsData object and creates a plot showing all the log2 fold changes of all features and the corresponding p-values. This looks at the entire dataset to help determine shifts in the whole community.
#'
#' @param countSTAT_results Required, an object of the class \code{countSTAT} created by \code{\link{countSTAT}}.
#' @param omicsData Required, omicsData ,an object of the class 'seqData' created by \code{\link{as.seqData}}.
#' @param x_axis Required, a character vector stating which variable to group data by and put on the x-axis, must be one of "Comparison" or one of the column names in omicsData$e_meta. Default is "Phylum".
#'@param plot_title Optional, a character vector to use as the plot title
#'@param leglab Optional, a character vector to use as the legend label
#'@param x_lab Optional, a character vector to use as the x-axis label
#'@param y_lab Optional, a character vector to use as the y-axis label
#'
#' @return A plot showing log2 fold change and p-values for all comparisons and all features
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' library(mintR)
#' data("rRNAdata")
#'
#'
#' rRNAdata <- group_designation(omicsData=rRNA_data, main_effects="Treatment")
#' norm_factors <- normalize_data(omicsData=rRNAdata, subset_fn="none", norm_fn="percentile", normalize=FALSE)
#' norm_data <- normalize_data(omicsData=rRNAdata, subset_fn="none", norm_fn="percentile", normalize=TRUE)
#' diffexp <- countSTAT(omicsData=rRNAdata, norm_factors=norm_factors, comparisons="all", test="dw", pval_adjust="none", pval_thresh=0.05)
#'
#' plot_all_diffabun(countSTAT_results=diffexp, omicsData=norm_data, x_axis="Phylum")
#'}
#'
#' @author Allison Thompson
#'
#' @export
#'
plot_all_diffabun <- function(countSTAT_results, omicsData, x_axis="Phylum", scales="fixed",
                                x_lab=NULL, y_lab=NULL, leglab=NULL, plot_title=NULL){

  library(RColorBrewer)
  library(ggplot2)
  library(reshape2)

  # initial checks #

  if(is.null(x_axis)){
    stop("x_axis cannot be NULL")
  }

  if(!(x_axis %in% c("Comparison", colnames(omicsData$e_meta)))){
    stop("x_axis must be either 'Comparison' or one of the column names in omicsData$e_meta")
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

  # end initial checks #

  # Extract differential expression results
  data <- countSTAT_results$allResults

  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  # Create a new plot for every test that was run
 lapply(attr(countSTAT_results, "Tests")$Test, function(t){
    # Extract only data for specific test
    data <- data[,grep(t, colnames(data))]

    # Extract and format adjusted p-values
    data.padj <- data[,grep("padj", colnames(data))]
    data.padj <- data.frame(Temp=rownames(data), data.padj)
    colnames(data.padj)[1] <- attr(omicsData, "cnames")$edata_cname
    colnames(data.padj)[-1] <- colnames(data)[grep("padj",colnames(data))]
    data.padj <- reshape2::melt(data.padj)
    colnames(data.padj)[3] <- "padj"
    data.padj$Comparison <- unlist(lapply(data.padj$variable, function(x) strsplit(as.character(x),"padj_")[[1]][2]))

    # Extract and format log2 fold change values
    data.logfc <- data[,grep("logFC", colnames(data))]
    data.logfc <- data.frame(Temp=rownames(data), data.logfc)
    colnames(data.logfc)[1] <- attr(omicsData, "cnames")$edata_cname
    colnames(data.logfc)[-1] <- colnames(data)[grep("logFC",colnames(data))]
    data.logfc <- reshape2::melt(data.logfc)
    colnames(data.logfc)[3] <- "logFC"
    data.logfc$Comparison <- unlist(lapply(data.logfc$variable, function(x) strsplit(as.character(x),"logFC_")[[1]][2]))

    # Merge adjusted p-values and log2fc data
    plot.data <- merge(data.padj, data.logfc, by=c(attr(omicsData, "cnames")$edata_cname, "Comparison"))
    # Merge data with e_meta
    plot.data <- merge(plot.data, omicsData$e_meta, by=attr(omicsData, "cnames")$edata_cname)

    lapply(attr(countSTAT_results, "comparisons")$comparison, function(c){
      c <- gsub("[^A-Za-z0-9_]","\\.",c)
      map <- ggplot2::aes_string(x=x_axis, y="logFC", colour="padj")
      p <- ggplot2::ggplot(subset(plot.data, Comparison==c), map) +
        ggplot2::geom_point(size=ifelse(subset(plot.data, Comparison==c)$padj <= 0.1, 3, 1))+
        ggplot2::scale_colour_gradientn(colours=rev(myPalette(100)), limits=c(0,0.1), name="P-value\n")+
        ggplot2::geom_hline(yintercept=log2(2),lty="dashed")+
        ggplot2::geom_hline(yintercept=-log2(2),lty="dashed")+
        ggplot2::geom_hline(yintercept=log2(10),lty="dotted")+
        ggplot2::geom_hline(yintercept=-log2(10),lty="dotted")+
        ggplot2::theme_bw()+
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle=90),
              axis.line.x = ggplot2::element_line(colour = "black"),
              axis.line.y = ggplot2::element_line(colour="black"),
              plot.background = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.border = ggplot2::element_blank())

      if(!is.null(x_lab)){
        p <- p + ggplot2::labs(x=x_lab)
      }else{
        p <- p + ggplot2::labs(x=x_axis)
      }

      if(!is.null(y_lab)){
        p <- p + ggplot2::labs(y=y_lab)
      }else{
        p <- p + ggplot2::labs(y="Log2 Fold Change")
      }

      if(!is.null(plot_title)){
        p <- p + ggplot2::ggtitle(plot_title)
      }else{
        p <- p + ggplot2::ggtitle(paste(t," ",c," Results"))
      }

      if(!is.null(leglab)){
        p <- p + ggplot2::guides(colour=guide_legend(title=leglab))
      }

      return(p)
   })

  })

}
