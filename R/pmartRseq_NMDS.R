#' NMDS plot for beta diversity
#'
#' This function creates an NMDS plot for a beta diversity object from vegan.
#'
#' @param res an object created by vegan::metaMDS
#' @param omicsData a seqData object
#' @param grp name of column to use for grouping variables
#' @param k number of dimensions
#' @param x_axis which dimension to put on the x-axis
#' @param y_axis which dimension to put on the y-axis
#' @param ellipses logical indicating whether or not to show ellipses on the plot
#'
#' @details After beta diversity is calculated, this function will create a plot of the results, along with ellipses, to show the different groupings.
#'
#' @return An NMDS plot of a beta diversity index.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' library(vegan)
#' data(rRNA_data)
#' rRNA_data <- group_designation(omicsData = rRNA_data, main_effects = c("treatment"))
#' rRNA_norm <- normalize_data(omicsData = rRNA_data, norm_fn = "css", normalize = TRUE)
#' rRNA_veg <- pmartRseq_to_vegan(omicsData = rRNA_norm)
#' rRNA_metamds <- vegan::metaMDS(rRNA_veg, distance = "bray", k = 4, autotransform = FALSE, na.rm = TRUE)
#' pmartRseq_NMDS(rRNA_metamds, as.factor(attr(rRNA_norm,"group_DF")[match(rownames(rRNA_veg),attr(rRNA_norm,"group_DF")[,attr(rRNA_norm,"cnames")$fdata_cname]),"Group]))
#'}
#'
#' @author Allison Thompson
#'
#' @export
pmartRseq_NMDS <- function(res, omicsData, grp, k, x_axis="NMDS1", y_axis="NMDS2", ellipses=TRUE){

  library(ggplot2)

  NMDS <- data.frame(SampleID=rownames(scores(res)),scores(res))
  colnames(NMDS)[1] <- attr(omicsData, "cnames")$fdata_cname
  NMDS <- merge(NMDS, attr(omicsData, "group_DF"), by=attr(omicsData, "cnames")$fdata_cname)

  testgrp <- table(NMDS[,grp])
  if(any(testgrp < 3)){
    names <- names(which(testgrp < 3))
    NMDS <- NMDS[-which(NMDS[,grp] %in% names),]
    NMDS[,grp] <- droplevels(NMDS[,grp])
  }

  # Format treatment, "group"
  if(any(levels(NMDS[,grp]) == "")){
    NMDS[,grp] <- as.character(NMDS[,grp])
    NMDS[,grp] <- as.factor(NMDS[,grp])
  }

  if(ellipses){
    # Function for drawing ellipses
    veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100){
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }

    df_ell <- data.frame()
    for(g in levels(NMDS[,grp])){
      df_ell <- rbind(df_ell,
                      cbind(as.data.frame(with(NMDS[NMDS[,grp]==g,],
                        veganCovEllipse(cov.wt(NMDS[which(NMDS[,grp]==g),c(x_axis,y_axis)],
                                               wt=rep(1/length(NMDS[which(NMDS[,grp]==g),x_axis])
                                                      ,length(NMDS[which(NMDS[,grp]==g),x_axis])))$cov,
                                        center=c(mean(NMDS[which(NMDS[,grp]==g),x_axis]),
                                                 mean(NMDS[which(NMDS[,grp]==g),y_axis]))))),
                            group=g))
    }

    # Plot
    map1 <- ggplot2::aes_string(x=x_axis, y=y_axis, color=grp)
    map2 <- ggplot2::aes_string(x=x_axis, y=y_axis, color="group")
    X1 <- ggplot2::ggplot(data = NMDS) +
      ggplot2::geom_point(map1, size=2.5, alpha=1) +
      ggplot2::geom_path(data = df_ell, map2, size=1.5, linetype=5, alpha=0.7)+
      ggplot2::theme_bw()+
      tggplot2::heme(aspect.ratio = 1,
            axis.text.x = ggplot2::element_text(size=20),
            axis.text.y = ggplot2::element_text(size=20),
            axis.title.x = ggplot2::element_text(size=20),
            axis.title.y = ggplot2::element_text(size=20),
            legend.title = ggplot2::element_text(size=15),
            legend.text = ggplot2::element_text(size=15),
            panel.grid = ggplot2::element_blank())
  }else{
    # Plot
    map1 <- ggplot2::aes_string(x=x_axis, y=y_axis, color=grp)
    X1 <- ggplot2::ggplot(data = NMDS) +
      ggplot2::geom_point(map1, size=2.5, alpha=1) +
      ggplot2::theme_bw()+
      ggplot2::theme(aspect.ratio = 1,
            axis.text.x = ggplot2::element_text(size=20),
            axis.text.y = ggplot2::element_text(size=20),
            axis.title.x = ggplot2::element_text(size=20),
            axis.title.y = ggplot2::element_text(size=20),
            legend.title = ggplot2::element_text(size=15),
            legend.text = ggplot2::element_text(size=15),
            panel.grid = ggplot2::element_blank())
  }
  X1
}

# mead_PCA <- function(res, grp){
#
#   plotdata <- data.frame(SampleID = rownames(scores(res)), PC1 = data.frame(scores(res))$NMDS1, PC2 = data.frame(scores(res))$NMDS2, Group = grp)
#
#   ggplot(plotdata, aes(x=PC1, y=PC2)) +
#     geom_point(aes(colour=Group), size=2.5, alpha=0.75)+
#     theme_bw()+
#     theme(aspect.ratio=1,
#           axis.text.x=element_text(size=20),
#           axis.text.y=element_text(size=20),
#           axis.title.x=element_text(size=20),
#           axis.title.y=element_text(size=20),
#           legend.title=element_text(size=15),
#           legend.text=element_text(size=15),
#           panel.grid=element_blank())
#
# }

