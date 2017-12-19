#' Generate Network Plot
#'
#' This function generates a network plot for the network data.
#'
#' @param netGraph an object of class 'networkGraph', created by \code{\link{pmartRseq_igraph}}
#' @param omicsData Optional, an object of the class 'seqData' usually created by \code{\link{as.seqData}}, if want to colour by taxonomy and/or scale vertices by abundance
#' @param modData Optional, an object of class 'modData', created by \code{\link{detect_modules}}, if want to colour by modules.
#' @param colour Optional, if desired, can colour vertices by a taxonomic level or 'Module' for module. Use 'NULL' if no colour is desired.
#' @param vsize Logical, should vertices be scaled by median abundance of taxa
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
#' mygraph <- pmartRseq_igraph(netData = mynetwork, coeff=0.6, pval=NULL, qval=0.05)
#' network_plot(omicsData = rRNA_data, netGraph = mygraph, colour = "Phylum", vsize = TRUE, legend.show = TRUE, legend.pos = "bottomleft")
#' }
#'
#' @author Allison Thompson
#'
#' @export

network_plot <- function(netGraph, omicsData=NULL, modData=NULL, colour="Phylum", vsize=FALSE, legend.show=TRUE, legend.pos="bottomleft"){
  library(igraph)

  ### Initial Checks ###

  if(!is.null(omicsData) & class(omicsData) != "seqData"){
    stop("omicsData must be an object of class 'seqData'")
  }

  if(!is.null(modData) & class(modData) != "modData"){
    stop("modData must be an object of class 'modData'")
  }

  if(class(netGraph) != "networkGraph"){
    stop("netGraph must be an object of class 'networkGraph'.")
  }

  if(!(colour %in% c(colnames(omicsData$e_meta), "Module"))){
    stop("colour must be either a column name in omicsData$e_meta or 'Module' if coloring by module.")
  }

  if(!(legend.pos %in% c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right","center"))){
    stop("legend.pos must be one of 'bottomright', 'bottom', 'bottomleft', 'left', 'topleft', 'top', 'topright', 'right', 'center'")
  }

  if(!is.logical(vsize)){
    stop("vsize must be one of TRUE or FALSE")
  }

  if(vsize & is.null(omicsData)){
    stop("omicsData must be provided in order to scale vertices by abundance")
  }

  if(!is.logical(legend.show)){
    stop("legend.show must be one of TRUE or FALSE")
  }

  ### End Initial Checks ###

  #if you want to associate the taxonomy with your taxa import a taxonomy key now
  taxa <- omicsData$e_meta

  net <- netGraph

  if(!is.null(attr(netGraph, "group_var"))){

    gN <- lapply(names(netGraph), function(x){
      gN <- netGraph[[x]]

      #gN <- simplify(graph.edgelist(as.matrix(tmp[,c("Row","Column")]), directed = FALSE))


      if(!is.null(taxa) & (colour %in% colnames(omicsData$e_meta))){
        #make taxonomy as vertex attribute
        vgn <- as.data.frame(as.matrix(V(gN)))
        vgn$Feature <- rownames(vgn)
        colnames(vgn)[which(colnames(vgn) == "Feature")] <- attr(netGraph, "cnames")$edata_cname

        vgn <- merge(vgn, taxa, by=attr(netGraph, "cnames")$edata_cname)
        vgn <- droplevels(vgn)
        vgn <- vgn[match(names(V(gN)), vgn[,attr(netGraph, "cnames")$edata_cname]),]
      }else if(!is.null(modData) & (colour == "Module")){
        vgn <- as.data.frame(as.matrix(V(gN)))
        vgn$Feature <- rownames(vgn)
        colnames(vgn)[which(colnames(vgn) == "Feature")] <- attr(netGraph, "cnames")$edata_cname

        vgn <- merge(vgn, mods[[x]], by=attr(netGraph, "cnames")$edata_cname)
        vgn$Module <- as.factor(vgn$Module)
        vgn <- droplevels(vgn)
        vgn <- vgn[match(names(V(gN)), vgn[,attr(netGraph, "cnames")$edata_cname]),]
      }else{
        vgn <- NULL
      }

      #you can alter the layout algorithm to create a more visually appealling plot- Fruchterman-Reingold & Kamada-Kawai layouts are common ones to try
      #Fruchterman-Reingold layout and edit the size
      l <- layout_with_fr(gN)
      l <- norm_coords(l, ymin= -1, ymax= 1, xmin=-1, xmax=1)

      if(vsize){
        if(attr(netGraph, "group_var") %in% colnames(attr(omicsData, "group_DF"))){
          samps <- attr(omicsData, "group_DF")[which(attr(omicsData, "group_DF")[,attr(netGraph, "group_var")] == x), attr(netGraph, "cnames")$fdata_cname]
        }else if(attr(netGraph, "group_var") %in% colnames(omicsData$f_data)){
          samps <- omicsData$f_data[which(omicsData$f_data[,attr(netGraph, "group_var")] == x), attr(netGraph, "cnames")$fdata_cname]
        }else{
          stop("Something went wrong, please double check group var in network data, group_DF in omics data, and f_data in omics data.")
        }

        size <- omicsData$e_data[,which(colnames(omicsData$e_data) %in% samps)]
        size <- apply(size, 1, function(x) median(x, na.rm=TRUE))
        size <- data.frame(Features=omicsData$e_data[,which(colnames(omicsData$e_data) == attr(netGraph, "cnames")$edata_cname)], Median=size)
        colnames(size)[which(colnames(size) == "Features")] <- attr(netGraph, "cnames")$edata_cname

        v.size <- as.data.frame(as.matrix(V(gN)))
        v.size$Feature <- rownames(v.size)
        colnames(v.size)[which(colnames(v.size) == "Feature")] <- attr(netGraph, "cnames")$edata_cname
        v.size <- merge(v.size, size, by=attr(netGraph, "cnames")$edata_cname)
        v.size <- v.size[match(names(V(gN)), v.size[,attr(netGraph, "cnames")$edata_cname]), "Median"]
        sizex <- max(v.size, na.rm=TRUE)/20
        v.size <- v.size/sizex
      }else{
        v.size = 5
      }

      if(!is.null(colour) & length(vgn[,colour]) > 0){
        #Get colors
        cols <- iwanthue(length(levels(vgn[,colour])), random=TRUE)

        my_colour <- cols[as.numeric(as.factor(vgn[,colour]))]

        plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.color=my_colour, vertex.label=NA, vertex.size=v.size, edge.width=(E(gN)$strength)*4, edge.color= ifelse(E(gN)$dirct > 0, "black", "red3"))
        #plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.color=my_colour, vertex.size=v.size, edge.width=(E(gN)$strength)*4, edge.color= ifelse(E(gN)$dirct > 0, "black", "red3"))
        if(legend.show){
          legend(legend.pos, legend=levels(as.factor(vgn[,colour])), col=cols, bty="n", pch=20, pt.cex=3, cex=0.5, horiz=FALSE, ncol=2, text.width=.1, x.intersp=.25)
        }
        title(paste(x," Network",sep=""))
      }else{
        plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.label=NA, vertex.size=v.size, edge.width=(E(gN)$strength)*4, edge.color= ifelse(E(gN)$dirct > 0, "black", "red3"))
        title(paste(x," Network",sep=""))
      }

      return(gN)
    })


    #Preliminary code to look at consistencies between networks- note that this won't run currently as we've only made gN1 in the code above. The graph.intersection function will keep only the edges that occur in both (all) graphs
    g_intersection <- Reduce(function(x, y) graph.intersection(x, y, byname=TRUE, keep.all.vertices=FALSE), gN)

    if(!is.null(taxa) & (colour %in% colnames(omicsData$e_meta))){
      #make taxonomy as vertex attribute
      ivgn <- as.data.frame(as.matrix(V(g_intersection)))
      ivgn$Feature <- rownames(ivgn)
      colnames(ivgn)[which(colnames(ivgn) == "Feature")] <- attr(netGraph, "cnames")$edata_cname

      ivgn <- merge(ivgn, taxa, by=attr(netGraph, "cnames")$edata_cname)
      ivgn <- droplevels(ivgn)
      ivgn <- ivgn[match(names(V(g_intersection)), ivgn[,attr(netGraph, "cnames")$edata_cname]),]
    }else if(!is.null(modData) & (colour == "Module")){
      ivgn <- as.data.frame(as.matrix(V(g_intersection)))
      ivgn$Feature <- rownames(ivgn)
      colnames(ivgn)[which(colnames(ivgn) == "Feature")] <- attr(netGraph, "cnames")$edata_cname

      ivgn <- merge(ivgn, do.call(rbind, mods), by=attr(netGraph, "cnames")$edata_cname)
      ivgn$Module <- as.factor(ivgn$Module)
      ivgn <- droplevels(ivgn)
      ivgn <- ivgn[match(names(V(g_intersection)), ivgn[,attr(netGraph, "cnames")$edata_cname]),]
    }else{
      ivgn <- NULL
    }

    if(vsize){
      # if(attr(netGraph, "group_var") %in% colnames(attr(omicsData, "group_DF"))){
      #   samps <- attr(omicsData, "group_DF")[which(attr(omicsData, "group_DF")[,attr(netGraph, "group_var")] == x), attr(netGraph, "cnames")$fdata_cname]
      # }else if(attr(netGraph, "group_var") %in% colnames(omicsData$f_data)){
      #   samps <- omicsData$f_data[which(omicsData$f_data[,attr(netGraph, "group_var")] == x), attr(netGraph, "cnames")$fdata_cname]
      # }else{
      #   stop("Something went wrong, please double check group var in network data, group_DF in omics data, and f_data in omics data.")
      # }

      size <- omicsData$e_data[,-which(colnames(omicsData$e_data) == attr(netGraph, "cnames")$edata_cname)]
      size <- apply(size, 1, function(x) median(x, na.rm=TRUE))
      size <- data.frame(Features=omicsData$e_data[,which(colnames(omicsData$e_data) == attr(netGraph, "cnames")$edata_cname)], Median=size)
      colnames(size)[which(colnames(size) == "Features")] <- attr(netGraph, "cnames")$edata_cname

      v.size <- as.data.frame(as.matrix(V(g_intersection)))
      v.size$Feature <- rownames(v.size)
      colnames(v.size)[which(colnames(v.size) == "Feature")] <- attr(netGraph, "cnames")$edata_cname
      v.size <- merge(v.size, size, by=attr(netGraph, "cnames")$edata_cname)
      v.size <- v.size[match(names(V(g_intersection)), v.size[,attr(netGraph, "cnames")$edata_cname]), "Median"]
      sizex <- max(v.size, na.rm=TRUE)/20
      v.size <- v.size/sizex
    }else{
      v.size = 5
    }

    if(!is.null(colour) & length(ivgn[,colour]) > 0){
      #Get colors
      icols <- iwanthue(length(levels(ivgn[,colour])), random=TRUE)

      imy_colour <- icols[as.numeric(as.factor(ivgn[,colour]))]


      plot(g_intersection, vertex.color=imy_colour, vertex.label=NA, vertex.size=v.size)
      if(legend.show){
        legend(legend.pos, legend=levels(as.factor(ivgn[,colour])), col=icols, bty="n", pch=20, pt.cex=3, cex=0.5, horiz=FALSE, ncol=2, text.width=.1, x.intersp=.5)
      }
      title("Graph Intersection")
    }else{
      plot(g_intersection, vertex.label=NA, vertex.size=v.size)
      title("Graph Intersection")
    }

  }else{
    gN <- net

    #gN <- simplify(graph.edgelist(as.matrix(tmp[,c("Row","Column")]), directed = FALSE))

    if(!is.null(taxa) & (colour %in% colnames(omicsData$e_meta))){
      #make taxonomy as vertex attribute
      vgn <- as.data.frame(as.matrix(V(gN)))
      vgn$Feature <- rownames(vgn)
      colnames(vgn)[which(colnames(vgn) == "Feature")] <- attr(netGraph, "cnames")$edata_cname

      vgn <- merge(vgn, taxa, by=attr(netGraph, "cnames")$edata_cname)
      vgn <- droplevels(vgn)
      vgn <- vgn[match(names(V(gN)), vgn[,attr(netGraph, "cnames")$edata_cname]),]
    }else if(!is.null(modData) & (colour == "Module")){
      vgn <- as.data.frame(as.matrix(V(gN)))
      vgn$Feature <- rownames(vgn)
      colnames(vgn)[which(colnames(vgn) == "Feature")] <- attr(netGraph, "cnames")$edata_cname

      vgn <- merge(vgn, mods, by=attr(netGraph, "cnames")$edata_cname)
      vgn$Module <- as.factor(vgn$Module)
      vgn <- droplevels(vgn)
      vgn <- vgn[match(names(V(gN)), vgn[,attr(netGraph, "cnames")$edata_cname]),]
    }else{
      vgn <- NULL
    }

    if(vsize){

      size <- omicsData$e_data[,-which(colnames(omicsData$e_data) == attr(netGraph, "cnames")$edata_cname)]
      size <- apply(size, 1, function(x) median(x, na.rm=TRUE))
      size <- data.frame(Features=omicsData$e_data[,which(colnames(omicsData$e_data) == attr(netGraph, "cnames")$edata_cname)], Median=size)
      colnames(size)[which(colnames(size) == "Features")] <- attr(netGraph, "cnames")$edata_cname

      v.size <- as.data.frame(as.matrix(V(gN)))
      v.size$Feature <- rownames(v.size)
      colnames(v.size)[which(colnames(v.size) == "Feature")] <- attr(netGraph, "cnames")$edata_cname
      v.size <- merge(v.size, size, by=attr(netGraph, "cnames")$edata_cname)
      v.size <- v.size[match(names(V(gN)), v.size[,attr(netGraph, "cnames")$edata_cname]), "Median"]
      sizex <- max(v.size, na.rm=TRUE)/20
      v.size <- v.size/sizex
    }else{
      v.size = 5
    }

    #you can alter the layout algorithm to create a more visually appealling plot- Fruchterman-Reingold & Kamada-Kawai layouts are common ones to try
    #Fruchterman-Reingold layout and edit the size
    l <- layout_with_fr(gN)
    l <- norm_coords(l, ymin= -1, ymax= 1, xmin=-1, xmax=1)

    if(!is.null(colour) & length(vgn[,colour]) > 0){
      #Get colors
      cols <- iwanthue(length(levels(vgn[,colour])), random=TRUE)

      my_colour <- cols[as.numeric(as.factor(vgn[,colour]))]

      plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.color=my_colour, vertex.label=NA, vertex.size=v.size, edge.width=(E(gN)$strength)*6, edge.color= ifelse(E(gN)$dirct > 0, "black", "red3"))
      if(legend.show){
        legend(legend.pos, legend=levels(as.factor(vgn[,colour])), col=cols, bty="n", pch=20, pt.cex=3, cex=0.5, horiz=FALSE, ncol=2, text.width=.1, x.intersp=.25)
      }
      title("Network")
    }else{
      plot(gN, rescale = FALSE, ylim=c(-1,1),xlim=c(-1,1), asp = 0, rescale=T, layout=l, vertex.label=NA, vertex.size=v.size, edge.width=(E(gN)$strength)*6, edge.color= ifelse(E(gN)$dirct > 0, "black", "red3"))
      title("Network")
    }
  }

  return(NULL)

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

