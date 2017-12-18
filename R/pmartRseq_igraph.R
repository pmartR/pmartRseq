#' Calculates a network graph
#'
#' This function calculates a network graph, using igraph, in order to create a network plot and/or calcuate network indices.
#'
#' @param netData an object of class 'netData', created by \code{link{network_calc}}
#' @param coeff Optional, cutoff value to use for the correlation coefficient
#' @param pval Optional, cutoff value to use for p-values
#' @param qval Optional, cutoff value to use for q-values
#'
#' @details Create an igraph object, which can then be used to visualize a network plot and/or calculate network indices.
#'
#' @return An object of class networkGraph (also a data.frame or list) containing the network.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' mynetwork <- network_calc(omicsData = rRNA_data)
#' myNetGraph <- pmartRseq_igraph(mynetwork)
#' }
#'
#' @author Allison Thompson
#'
#' @export

pmartRseq_igraph <- function(netData, coeff=NULL, pval=NULL, qval=NULL){
  library(igraph)

  ### Initial Checks ###

  if(class(netData) != "corrRes"){
    stop("netGraph must be an object of class 'corrRes'.")
  }

  if(!is.null(coeff) & !is.numeric(coeff)){
    stop("coeff must be a numeric value")
  }

  if(!is.null(pval) & (pval < 0 | pval > 1 | !is.numeric(pval))){
    stop("pval must be a numeric value between 0 and 1.")
  }

  if(!is.null(qval) & (qval < 0 | qval > 1 | !is.numeric(qval))){
    stop("qval must be a numeric value between 0 and 1.")
  }

  ### End Initial Checks ###

 # Function to make an igraph object from a network object

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

  # If network was created in groups, create object in groups.
  if(!is.null(attr(netData, "group_var"))){

    gN <- lapply(unique(net$Group), function(x){
      tmp <- subset(net, Group == x)

      gN <- simplify(graph.edgelist(as.matrix(tmp[,c("Row","Column")]), directed = FALSE))

      #In our graph, carry through correlation coefficients
      E(gN)$strength <- abs(tmp$cor.coeff)
      E(gN)$dirct <- tmp$cor.coeff

      return(gN)
    })

    names(gN) <- unique(net$Group)

    attr(gN, "group_var") <- attr(netData, "group_var")

  }else{
    # Otherwise, create object over all data.
    tmp <- net

    gN <- simplify(graph.edgelist(as.matrix(tmp[,c("Row","Column")]), directed = FALSE))

    #In our graph, carry through correlation coefficients
    E(gN)$strength <- abs(tmp$cor.coeff)
    E(gN)$dirct <- tmp$cor.coeff


  }

  # Attach thresholds as attribute
  attr(gN, "thresholds") <- list(CorCoeff=coeff, PValue=pval, QValue=qval)
  attr(gN, "cnames") <- attr(netData, "cnames")
  attr(gN, "e_meta") <- attr(netData, "e_meta")

  # Make class networkGraph
  class(gN) <- c("networkGraph",class(gN))

  return(gN)

}
