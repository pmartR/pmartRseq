#' Generate Network Plot
#'
#' This function generates a network plot for the network data.
#'
#' @param omicsData an object of the class 'seqData' usually created by \code{\link{as.seqData}}.
#' @param netData an object of class 'netData', created by \code{link{network_calc}}
#' @param coeff Optional, cutoff value to use for the correlation coefficient
#' @param pval Optional, cutoff value to use for p-values
#' @param qval Optional, cutoff value to use for q-values
#' @param colour Optional, if desired, can colour vertices by a taxonomic level
#' @param legend.show Logical, should a legend be shown
#' @param legend.pos Optional, if legend==TRUE, where to position the legend. Default is 'bottomleft'.
#'
#' @details A network graph is created for the network(s) that were generated.
#'
#' @return A network graph.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' mynetwork <- network_calc(omicsData = rRNA_data)
#' network_plot(mynetwork)
#' }
#'
#' @author Allison Thompson
#'
#' @export

network_plot <- function(omicsData, netData, coeff=NULL, pval=NULL, qval=NULL, colour="Phylum", legend.show=TRUE, legend.pos="bottomleft"){
  library(igraph)

  #if you want to associate the taxonomy with your taxa import a taxonomy key now
  taxa <- omicsData$e_meta

  net <- netData

  # Correlation Coefficient cutoff
  if(!is.null(coeff)){
    net <- subset(net, abs(cor.coeff) >= coeff)
  }
  # p value cutoff
  if(!is.null(pval)){
    net <- subset(net, p.value <= pval)
  }
  # q value cutoff
  if(!is.null(qval)){
    net <- subset(net, q.value <= qval)
  }

  if(!is.null(attr(netData, "group_var"))){

    gN <- lapply(unique(net$Group), function(x){
      tmp <- subset(net, Group == x)

      gN <- simplify(graph.edgelist(as.matrix(tmp[,c("Row","Column")]), directed = FALSE))

      if(!is.null(taxa)){
        #make taxonomy as vertex attribute
        vgn <- as.data.frame(as.matrix(V(gN)))
        vgn$Feature <- rownames(vgn)
        colnames(vgn)[which(colnames(vgn) == "Feature")] <- attr(omicsData, "cnames")$edata_cname
        vgn <- merge(vgn, taxa, by=attr(omicsData, "cnames")$edata_cname)
        vgn <- droplevels(vgn)
      }else{
        vgn <- NULL
      }

      #In our graph we will base the color and thickness of lines on edge attributes, here we add these based on the stregth of the correlation and pos/neg
      E(gN)$strength <- abs(tmp$cor.coeff)
      E(gN)$dirct <- tmp$cor.coeff

      #you can alter the layout algorithm to create a more visually appealling plot- Fruchterman-Reingold & Kamada-Kawai layouts are common ones to try
      #Fruchterman-Reingold layout and edit the size
      l <- layout_with_fr(gN)
      l <- norm_coords(l, ymin= -1, ymax= 1, xmin=-1, xmax=1)

      if(!is.null(colour) & length(vgn[,colour]) > 0){
        #Get colors
        cols <- iwanthue(length(levels(vgn[,colour])), random=TRUE)

        my_colour <- cols[as.numeric(as.factor(vgn[,colour]))]

        plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.color=my_colour, vertex.label=NA, vertex.size=5, edge.width=(E(gN)$strength)*4, edge.color= ifelse(E(gN)$dirct > 0, "blue", "red3"))
        if(legend.show){
          legend(legend.pos, legend=levels(as.factor(vgn[,colour])), col=cols, bty="n", pch=20, pt.cex=3, cex=0.5, horiz=FALSE, ncol=2, text.width=.1, x.intersp=.25)
        }
        title(paste(x," Network",sep=""))
      }else{
        plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.label=NA, vertex.size=5, edge.width=(E(gN)$strength)*4, edge.color= ifelse(E(gN)$dirct > 0, "blue", "red3"))
        title(paste(x," Network",sep=""))
      }

      return(gN)
    })


    #Preliminary code to look at consistencies between networks- note that this won't run currently as we've only made gN1 in the code above. The graph.intersection function will keep only the edges that occur in both (all) graphs
    g_intersection <- Reduce(function(x, y) graph.intersection(x, y, byname=TRUE, keep.all.vertices=FALSE), gN)

    if(!is.null(taxa)){
      #make taxonomy as vertex attribute
      ivgn <- as.data.frame(as.matrix(V(g_intersection)))
      ivgn$Feature <- rownames(ivgn)
      colnames(ivgn)[which(colnames(ivgn) == "Feature")] <- attr(omicsData, "cnames")$edata_cname
      ivgn <- merge(ivgn, taxa, by=attr(omicsData, "cnames")$edata_cname)
      ivgn <- droplevels(ivgn)
    }

    if(!is.null(colour) & length(ivgn[,colour]) > 0){
      #Get colors
      icols <- iwanthue(length(levels(ivgn[,colour])), random=TRUE)

      imy_colour <- icols[as.numeric(as.factor(ivgn[,colour]))]


      plot(g_intersection, vertex.color=imy_colour, vertex.label=NA, vertex.size=5)
      if(legend.show){
        legend(legend.pos, legend=levels(as.factor(ivgn[,colour])), col=icols, bty="n", pch=20, pt.cex=3, cex=0.5, horiz=FALSE, ncol=2, text.width=.1, x.intersp=.5)
      }
      title("Graph Intersection")
    }else{
      plot(g_intersection, vertex.label=NA, vertex.size=5)
      title("Graph Intersection")
    }

  }else{
    tmp <- net

    gN <- simplify(graph.edgelist(as.matrix(tmp[,c("Row","Column")]), directed = FALSE))

    if(!is.null(taxa)){
      #make taxonomy as vertex attribute
      vgn <- as.data.frame(as.matrix(V(gN)))
      vgn$Feature <- rownames(vgn)
      colnames(vgn)[which(colnames(vgn) == "Feature")] <- attr(omicsData, "cnames")$edata_cname
      vgn <- merge(vgn, taxa, by=attr(omicsData, "cnames")$edata_cname)
      vgn <- droplevels(vgn)
    }

    #In our graph we will base the color and thickness of lines on edge attributes, here we add these based on the stregth of the correlation and pos/neg
    E(gN)$strength <- abs(tmp$cor.coeff)*0.6
    E(gN)$dirct <- tmp$cor.coeff

    #you can alter the layout algorithm to create a more visually appealling plot- Fruchterman-Reingold & Kamada-Kawai layouts are common ones to try
    #Fruchterman-Reingold layout and edit the size
    l <- layout_with_fr(gN)
    l <- norm_coords(l, ymin= -1, ymax= 1, xmin=-1, xmax=1)

    if(!is.null(colour) & length(vgn[,colour]) > 0){
      #Get colors
      cols <- iwanthue(length(levels(vgn[,colour])), random=TRUE)

      my_colour <- cols[as.numeric(as.factor(vgn[,colour]))]

      plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.color=my_colour, vertex.label=NA, vertex.size=5, edge.width=(E(gN)$strength)*6, edge.color= ifelse(E(gN)$dirct > 0, "blue", "red3"))
      if(legend.show){
        legend(legend.pos, legend=levels(as.factor(vgn[,colour])), col=cols, bty="n", pch=20, pt.cex=3, cex=0.5, horiz=FALSE, ncol=2, text.width=.1, x.intersp=.25)
      }
      title("Network")
    }else{
      plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.label=NA, vertex.size=5, edge.width=(E(gN)$strength)*6, edge.color= ifelse(E(gN)$dirct > 0, "blue", "red3"))
      title("Network")
    }
  }

}


iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100,
                     plot=FALSE, random=FALSE) {
  # Presently doesn't allow hmax > hmin (H is circular)
  # n: number of colours
  # hmin: lower bound of hue (0-360)
  # hmax: upper bound of hue (0-360)
  # cmin: lower bound of chroma (0-180)
  # cmax: upper bound of chroma (0-180)
  # lmin: lower bound of luminance (0-100)
  # lmax: upper bound of luminance (0-100)
  # plot: plot a colour swatch?
  # random: should clustering be random? (if FALSE, seed will be set to 1,
  #         and the RNG state will be restored on exit.)
  require(colorspace)
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0,
            hmax <= 360, cmax <= 180, lmax <= 100,
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            n > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  set.seed(64)
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1),
                                   seq(-100, 100, 5),
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin &
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    #hcl <- hcl[-which(is.na(coords(hcl)[, 2]))]
    lab <- as(hcl, 'LAB')
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), n, iter.max=50)
  if (isTRUE(plot)) {
    swatch(hex(LAB(clus$centers)))
  }
  hex(LAB(clus$centers))
}

