#' NMDS plot for beta diversity
#'
#' This function creates an NMDS plot for a beta diversity object from vegan.
#'
#' @param res an object created by vegan::metaMDS
#' @param grp vector of grouping variables
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
pmartRseq_NMDS <- function(res,grp,ellipses){

  library(ggplot2)

  # Extract component scores
  NMDS1 <- data.frame(scores(res))$NMDS1
  NMDS2 <- data.frame(scores(res))$NMDS2

  testgrp <- table(grp)
  if(any(testgrp < 3)){
    names <- names(which(testgrp < 3))
    ids <- which(grp %in% names)
    grp <- grp[-ids]
    grp <- droplevels(grp)
    NMDS1 <- NMDS1[-ids]
    NMDS2 <- NMDS2[-ids]
  }

  # Format treatment, "group"
  Treatment <- grp
  if(any(levels(Treatment) == "")){
    Treatment <- as.character(Treatment)
    Treatment <- as.factor(Treatment)
  }

  # Aggregate data using mean
  NMDS <- data.frame(NMDS1, NMDS2, Treatment)
  NMDS.mean = aggregate(NMDS[,1:2], list(group=Treatment), mean)

  if(ellipses){
    # Function for drawing ellipses
    veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100){
      theta <- (0:npoints) * 2 * pi/npoints
      Circle <- cbind(cos(theta), sin(theta))
      t(center + scale * t(Circle %*% chol(cov)))
    }

    # Create a new dataframe with data and ellipses
    df_ell <- data.frame()
    for(g in levels(NMDS$Treatment)){
      df_ell <- rbind(df_ell,
                      cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                      veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))),
                      group=g))
    }

    # Plot
    X1 <- ggplot(data = NMDS, aes(NMDS1, NMDS2)) +
      geom_point(aes(color = Treatment), size=2.5, alpha=1) +
      geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2, colour=group), size=1.5, linetype=5, alpha=0.7)+
      theme_bw()+
      theme(aspect.ratio=1,
            axis.text.x=element_text(size=20),
            axis.text.y=element_text(size=20),
            axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            legend.title=element_text(size=15),
            legend.text=element_text(size=15),
            panel.grid=element_blank())
  }else{
    df_ell <- data.frame()
    for(g in levels(NMDS$Treatment)){
      df_ell <- rbind(df_ell,
                      cbind(as.data.frame(NMDS[NMDS$Treatment==g,],group=g)))
    }

    # Plot
    X1 <- ggplot(data = NMDS, aes(NMDS1, NMDS2)) +
      geom_point(aes(color = Treatment), size=2.5, alpha=0.75) +
      theme_bw()+
      theme(aspect.ratio=1,
            axis.text.x=element_text(size=20),
            axis.text.y=element_text(size=20),
            axis.title.x=element_text(size=20),
            axis.title.y=element_text(size=20),
            legend.title=element_text(size=15),
            legend.text=element_text(size=15),
            panel.grid=element_blank())
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

