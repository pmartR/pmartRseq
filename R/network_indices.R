#' Calculate Indices (metrics) for Network
#'
#' This function calculates a variety of indices for a network
#'
#' @param netGraph an object of class 'networkGraph', created by \code{link{pmartRseq_igraph}}
#'
#' @details A network graph is created for the network(s) that were generated.
#'
#' @return A list containing 1) Metrics - metric indices for the graph, and 2) Random - metric indices for a random graph to compare against.
#'
#' @references  Ju et al., 2014. Taxonomic relatedness shapes bacterial assembly in activated sludge of globally distributed wastewater treatment plants. Environ Microbiol. 2014 Aug;16(8):2421-32)
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' mynetwork <- network_calc(omicsData = rRNA_data)
#' mygraph <- pmartRseq_igraph(netData = mynetwork, coeff=0.6, pval=NULL, qval=0.05)
#' myindices <- network_indices(netGraph = mygraph)
#' }
#'
#' @author Allison Thompson
#'
#' @export

network_indices <- function(netGraph){

  library(igraph)

  ### Initial Checks ###

  if(class(netGraph)[1] != "networkGraph"){
    stop("netGraph must be an object of class 'networkGraph'.")
  }

  ### End Initial Checks ###

  if(!is.null(attr(netGraph, "group_var"))){
    indices <- lapply(names(netGraph), function(x){

      # vertices
      btwn <- betweenness(netGraph[[x]])
      dgr <- degree(netGraph[[x]])
      vertices <- data.frame(Betweenness=btwn, Degree=dgr)
      vertices <- data.frame(Features=rownames(vertices), vertices)
      colnames(vertices)[which(colnames(vertices) == "Features")] = attr(netGraph, "cnames")$edata_cname

      # degree
      dgr_dist <- degree_distribution(netGraph[[x]])

      # edge
      edge_btwn <- edge_betweenness(netGraph[[x]])
      edge_btwn <- data.frame(EdgeVertices=attr(E(netGraph[[x]]),"vnames"), EdgeBetweenness=edge_btwn)

      # Transitivity
      trans <- transitivity(netGraph[[x]])

      # Network wide metrics
      meandist <- mean_distance(netGraph[[x]])

      # Other metrics
      other <- list(DegreeDistribution=dgr_dist, Transitivity=trans, MeanDistance=meandist)

      # Put it all together
      res <- list(Vertices=vertices, Edges=edge_btwn, Other=other)

      ## Create a random netowrk to compare against

      #Makes a single random network of same size as real data set- has to be set manually-
      v <- length(V(netGraph[[x]]))   #number of vertices
      e <- length(E(netGraph[[x]]))  #number of edges
      rand <- erdos.renyi.game(v, e, type="gnm")

      #Vetrex metrics - betweenness and degree
      rand_btwn <- betweenness(rand)
      rand_dgr <- degree(rand)
      rand_vertices <- data.frame(Betweenness=rand_btwn, Degree=rand_dgr)

      # Can also look at the distribution of degree in network
      rand_dgr_dist <- degree_distribution(rand)

      # edge
      rand_edge_btwn <- edge_betweenness(rand)
      rand_edge_btwn <- data.frame(EdgeBetweenness=rand_edge_btwn)

      # Transitivity
      rand_trans <- transitivity(rand)

      # Network wide metrics
      rand_meandist <- mean_distance(rand)

      # Other metrics
      rand_other <- list(DegreeDistribution=rand_dgr_dist, Transitivity=rand_trans, MeanDistance=rand_meandist)

      # Put it all together
      rand_res <- list(Vertices=rand_vertices, Edges=rand_edge_btwn, Other=rand_other)

      # Put graph and random metrics together
      results <- list(Metrics=res, Random=rand_res)
      return(results)

    })

    names(indices) <- names(netGraph)
    attr(indices, "group_var") <- attr(netGraph, "group_var")

  }else{
    # vertices
    btwn <- betweenness(netGraph)
    dgr <- degree(netGraph)
    vertices <- data.frame(Betweenness=btwn, Degree=dgr)
    vertices <- data.frame(Features=rownames(vertices), vertices)
    colnames(vertices)[which(colnames(vertices) == "Features")] = attr(netGraph, "cnames")$edata_cname

    # degree
    dgr_dist <- degree_distribution(netGraph)

    # edge
    edge_btwn <- edge_betweenness(netGraph)
    edge_btwn <- data.frame(EdgeVertices=attr(E(netGraph),"vnames"), EdgeBetweenness=edge_btwn)

    # Transitivity
    trans <- transitivity(netGraph)

    # Network wide metrics
    meandist <- mean_distance(netGraph)

    # Other metrics
    other <- list(DegreeDistribution=dgr_dist, Transitivity=trans, MeanDistance=meandist)

    # Put it all together
    res <- list(Vertices=vertices, Edges=edge_btwn, Other=other)

    ## Create a random netowrk to compare against

    #Makes a single random network of same size as real data set- has to be set manually-
    v <- length(V(netGraph))   #number of vertices
    e <- length(E(netGraph))  #number of edges
    rand <- erdos.renyi.game(v, e, type="gnm")

    #Vetrex metrics - betweenness and degree
    rand_btwn <- betweenness(rand)
    rand_dgr <- degree(rand)
    rand_vertices <- data.frame(Betweenness=rand_btwn, Degree=rand_dgr)

    # Can also look at the distribution of degree in network
    rand_dgr_dist <- degree_distribution(rand)

    # edge
    rand_edge_btwn <- edge_betweenness(rand)
    rand_edge_btwn <- data.frame(EdgeBetweenness=rand_edge_btwn)

    # Transitivity
    rand_trans <- transitivity(rand)

    # Network wide metrics
    rand_meandist <- mean_distance(rand)

    # Other metrics
    rand_other <- list(DegreeDistribution=rand_dgr_dist, Transitivity=rand_trans, MeanDistance=rand_meandist)

    # Put it all together
    rand_res <- list(Vertices=rand_vertices, Edges=rand_edge_btwn, Other=rand_other)

    # Put graph and random metrics together
    results <- list(Metrics=res, Random=rand_res)

    indices <- results
  }

  attr(indices, "cnames") <- attr(netGraph, "cnames")
  attr(indices, "thresholds") <- attr(netGraph, "thresholds")
  attr(indices, "e_meta") <- attr(netGraph, "e_meta")

  class(indices) <- c("netIndices",class(indices))

  return(indices)

}
