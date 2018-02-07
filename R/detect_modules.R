#' Detect modules for Network
#'
#' This function detects modules for a network.
#'
#' @param netGraph an object of class 'networkGraph', created by \code{link{pmartRseq_igraph}}
#' @param cluster which clustering method to use. Must be one of 'edge_betweenness', 'fast_greedy', 'infomap', 'label_prop', 'leading_eigen', 'louvain', 'optimal', 'spinglass', or 'walktrap'. Default is 'louvain'.
#' @param cutoff Any modules which come back with a total of less than this number of members will all be grouped together.
#'
#' @details A network graph is created for the network(s) that were generated.
#'
#' @return An object of class 'modData', with group designations for every vertex in the graph.
#'
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' mynetwork <- network_calc(omicsData = rRNA_data)
#' mygraph <- pmartRseq_igraph(netData = mynetwork, coeff=0.6, pval=NULL, qval=0.05)
#' mymods <- detect_modules(netGraph = mygraph)
#' }
#'
#' @author Allison Thompson
#'
#' @export

detect_modules <- function(netGraph, cluster="louvain", cutoff=5){

  library(igraph)

  ### Initial Checks ###

  if(cutoff < 0 | length(cutoff) > 1 | !is.numeric(cutoff)){
    stop("cutoff must be a non-negative integer of length 1.")
  }

  if(class(netGraph)[1] != "networkGraph"){
    stop("netGraph must be an object of class 'networkGraph'.")
  }

  ### End Initial Checks ###

  # if per group, detect modules per group
  if(!is.null(attr(netGraph, "group_var"))){

    # use appropriate clustering method
    membs <- lapply(names(netGraph), function(x){
      if(cluster == "edge_betweenness"){
        mods <- igraph::cluster_edge_betweenness(netGraph[[x]])
      }else if(cluster == "fast_greedy"){
        mods <- igraph::cluster_fast_greedy(netGraph[[x]])
      }else if(cluster == "infomap"){
        mods <- igraph::cluster_infomap(netGraph[[x]])
      }else if(cluster == "label_prop"){
        mods <- igraph::cluster_label_prop(netGraph[[x]])
      }else if(cluster == "leading_eigen"){
        mods <- igraph::cluster_leading_eigen(netGraph[[x]])
      }else if(cluster == "louvain"){
        mods <- igraph::cluster_louvain(netGraph[[x]])
      }else if(cluster == "optimal"){
        mods <- igraph::cluster_optimal(netGraph[[x]])
      }else if(cluster == "spinglass"){
        mods <- igraph::cluster_spinglass(netGraph[[x]])
      }else if(cluster == "walktrap"){
        mods <- igraph::cluster_walktrap(netGraph[[x]])
      }else{
        stop("cluster must be one of 'edge_betweenness', 'fast_greedy', 'infomap', 'label_prop', 'leading_eigen', 'louvain', 'optimal', 'spinglass', or 'walktrap'.")
      }

      # calculate membership of each module
      membs <- igraph::membership(mods)
      membs <- data.frame(Features=names(membs), Module=as.matrix(membs))

      # change modules which have less than cutoff members
      if(any(table(membs$Module) <= cutoff)){
        ids <- table(membs$Module)
        ids <- names(ids)[which(ids <= cutoff)]
        other <- subset(membs, Module %in% ids)

        membs <- membs[-which(membs$Features %in% other$Features),]

        other$Module <- "Non-Modular"

        membs <- rbind(membs, other)
      }

      colnames(membs) <- c(attr(netGraph, "cnames")$edata_cname, "Module")

      attr(membs, "modularity") <- modularity(mods)
      attr(membs, "sizes") <- sizes(mods)

      return(membs)

    })

    attr(membs, "group_var") <- attr(netGraph, "group_var")
    names(membs) <- names(netGraph)

  }else{
    # use appropriate clustering method
    if(cluster == "edge_betweenness"){
      mods <- igraph::cluster_edge_betweenness(netGraph)
    }else if(cluster == "fast_greedy"){
      mods <- igraph::cluster_fast_greedy(netGraph)
    }else if(cluster == "infomap"){
      mods <- igraph::cluster_infomap(netGraph)
    }else if(cluster == "label_prop"){
      mods <- igraph::cluster_label_prop(netGraph)
    }else if(cluster == "leading_eigen"){
      mods <- igraph::cluster_leading_eigen(netGraph)
    }else if(cluster == "louvain"){
      mods <- igraph::cluster_louvain(netGraph)
    }else if(cluster == "optimal"){
      mods <- igraph::cluster_optimal(netGraph)
    }else if(cluster == "spinglass"){
      mods <- igraph::cluster_spinglass(netGraph)
    }else if(cluster == "walktrap"){
      mods <- igraph::cluster_walktrap(netGraph)
    }else{
      stop("cluster must be one of 'edge_betweenness', 'fast_greedy', 'infomap', 'label_prop', 'leading_eigen', 'louvain', 'optimal', 'spinglass', or 'walktrap'.")
    }

    # calculate membership of each module
    membs <- membership(mods)
    membs <- data.frame(Features=names(membs), Module=as.matrix(membs))

    # change modules which have less than cutoff members
    if(any(table(membs$Module) <= cutoff)){
      ids <- table(membs$Module)
      ids <- names(ids)[which(ids <= cutoff)]
      other <- subset(membs, Module %in% ids)

      membs <- membs[-which(membs$Features %in% other$Features),]

      other$Module <- "Non-Modular"

      membs <- rbind(membs, other)
    }

    colnames(membs) <- c(attr(netGraph, "cnames")$edata_cname, "Module")

    attr(membs, "modularity") <- modularity(mods)
    attr(membs, "sizes") <- sizes(mods)

  }

  attr(membs, "cnames") <- attr(netGraph, "cnames")
  attr(membs, "e_meta") <- attr(netGraph, "e_meta")

  class(membs) <- c("modData",class(membs))

  return(membs)
}
