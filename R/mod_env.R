#' Function to calculate module eigen vectors and correlate against environmental variables
#'
#' This function calculates module eigen vectors and then allows for correlation against environmental variables.
#'
#' @param omicsData An object of the class 'seqData' usually created by \code{\link{as.seqData}}. omicsData$f_data must contain environmental variables to compare against.
#' @param modData An object of class 'modData', created by \code{\link{detect_modules}}, if want to colour by modules.
#' @param envVars A character vector with the names of the environmental variables to compare against. Must all be column names in omicsData$f_data. If NULL, will not run correlation with environmental variables, but will still calculate pca for modules.
#' @param pca.method A string specifying the pca method to use. Default is 'svd'. Other options include 'nipals', 'rnipals', 'bpca', 'ppca', 'svdImpute', 'robusePca', 'nlpca', 'llsImpute', and 'llsImputeAll'.
#' @param cor.method A string specifying the correlation method to use. Default is 'spearman', which is nonparametric and doesn't make any assumptions about the data. Other options include 'kendall' and 'pearson'. 'Pearson' assumes normality - check assumptions and make any necessary transformations before using.
#' @param use A string specifying either 'pairwise' or 'complete'. 'pairwise' does pairwise deletion of cases and 'complete' selects just the complete cases. Default is 'pairwise'.
#' @param padjust A string specifying which adjustment for multiple test should be used. Default is 'BH'. Options are 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', and 'none'.
#'
#' @details A network graph is created for the network(s) that were generated.
#'
#' @return A list (or list of lists, if running in groups) containing pca - the pca values for every module and corr - the correlation values for every pca/module against environmental variables.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' mynetwork <- network_calc(omicsData = rRNA_data)
#' mygraph <- pmartRseq_igraph(netData = mynetwork, coeff=0.6, pval=NULL, qval=0.05)
#' mymods <- detect_modules(netGraph = mygraph)
#' myeigen <- mod_eigen(omicsData = rRNA_data, modData=mymods, envVars=c("MBC","MBN","SOC"), method="spearman", use="pairwise", padjust="BH")
#' myeigen
#' }
#'
#' @references corr.test from psych: Revelle, W. (2017) psych: Procedures for Personality and Psychological Research, Northwestern University, Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.7.8. ;   pca from pcaMethods: Stacklies, W., Redestig, H., Scholz, M., Walther, D. and Selbig, J.  pcaMethods -- a Bioconductor package providing PCA methods for incomplete data. Bioinformatics, 2007, 23, 1164-1167
#'
#' @author Allison Thompson
#'
#' @export

mod_env <- function(omicsData, modData, envVars, pca.method="svd", cor.method="spearman", use="pairwise", padjust="BH"){

  library(pcaMethods)
  library(psych)

  ### Initial Checks ###

  if(!is.null(omicsData) & class(omicsData)[1] != "seqData"){
    stop("omicsData must be an object of class 'seqData'")
  }

  if(!is.null(modData) & class(modData)[1] != "modData"){
    stop("modData must be an object of class 'modData'")
  }

  if(!all(envVars %in% colnames(omicsData$f_data)) & !is.null(envVars)){
    stop("all envVars must be found in colnames(omicsData$f_data)")
  }

  if(!(pca.method %in% c("svd", "nipals", "rnipals", "bpca", "ppca", "svdImpute", "robusePca", "nlpca", "llsImpute", "llsImputeAll")) | length(pca.method) != 1){
    stop("pca.method must be one (and only one) of 'svd', 'nipals', 'rnipals', 'bpca', 'ppca', 'svdImpute', 'robusePca', 'nlpca', 'llsImpute', or 'llsImputeAll'")
  }

  if(!(cor.method %in% c("spearman","pearson","kendall")) | length(cor.method) != 1){
    stop("cor.method must be one (and only one) of 'spearman', 'pearson', or 'kendall'")
  }

  if(!(use %in% c("pairwise","complete")) | length(use) != 1){
    stop("use must be one (and only one) of 'pairwise' or 'complete'")
  }

  if(!(padjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr","none")) | length(padjust) != 1){
    stop("padjust must be one (and only one) of holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', or 'none'")
  }

  ### End Initial Checks ###

  if(!is.null(attr(modData, "group_var"))){

    eigen_res <- lapply(names(modData), function(x){

      # Get module information for specified group
      mods <- modData[[x]]

      # Obtain samples in specified group
      if(attr(modData, "group_var") %in% colnames(attr(omicsData, "group_DF"))){
        samps <- attr(omicsData, "group_DF")[which(attr(omicsData, "group_DF")[,attr(modData, "group_var")] == x), attr(modData, "cnames")$fdata_cname]
      }else if(attr(netGraph, "group_var") %in% colnames(omicsData$f_data)){
        samps <- omicsData$f_data[which(omicsData$f_data[,attr(modData, "group_var")] == x), attr(modData, "cnames")$fdata_cname]
      }else{
        stop("Something went wrong, please double check group var in module data, group_DF in omics data, and f_data in omics data.")
      }

      # Subset omicsData down to samples in group
      abunData <- omicsData$e_data[,which(colnames(omicsData$e_data) %in% c(attr(omicsData, "cnames")$edata_cname, as.character(samps)))]

      pcas <- lapply(unique(mods$Module), function(y){
        # Subset module data down to specific module
        module <- subset(mods, Module == y)

        # Merge module information with abundance information and format
        module <- merge(module, abunData, by=attr(modData, "cnames")$edata_cname, all.x=TRUE)
        rownames(module) <- module[,which(colnames(module) == attr(modData, "cnames")$edata_cname)]
        module <- module[,-which(colnames(module) %in% c(attr(modData, "cnames")$edata_cname, "Module"))]
        module <- t(module)

        # pca <- dudi.pca(module, center=FALSE, scale=FALSE, scannf=FALSE, nf=2)

        # Use pca to calculate principal component axes and format
        pca <- pcaMethods::pca(object=t(module), nPcs=2, scale="none", center=FALSE, method=pca.method)
        pca <- data.frame(Samples=rownames(pca@loadings), pca@loadings)
        colnames(pca) <- c(attr(modData, "cnames")$fdata_cname, paste(colnames(pca)[2],"_Module",y,"_Group",x,sep=""), paste(colnames(pca)[3],"_Module",y,"_Group",x,sep=""))
        return(pca)
      })

      # Combine data from all modules in group
      pcas <- Reduce(function(x, y) merge(x, y, by=attr(modData, "cnames")$fdata_cname, all=TRUE), pcas)

      # If envirnmental variables, correlate module PCAs to environmental variables
      if(!is.null(envVars)){
        # Format data
        pcas <- merge(pcas, omicsData$f_data[,which(colnames(omicsData$f_data) %in% c(attr(omicsData, "cnames")$fdata_cname, envVars))], by=attr(omicsData, "cnames")$fdata_cname)
        pcas$group_var <- x

        # Correlation test of environmental variables with module PCAs
        env.cor <- psych::corr.test(x=pcas[,grep("PC[12]_Module",colnames(pcas))], y=pcas[,which(colnames(pcas) %in% envVars)], method=cor.method, use=use, adjust=padjust)

        # Format correlation coefficient and p-value results
        env.r <- reshape2::melt(env.cor$r)
        colnames(env.r) <- c("Module","EnvVar","CorrCoeff")
        env.p <- reshape2::melt(env.cor$p)
        colnames(env.p) <- c("Module","EnvVar","p.value")
        env.cor <- merge(env.r, env.p, by=c("Module","EnvVar"))

        res <- list(pca=pcas, corr=env.cor)
      }else{
        res <- list(pca=pcas)
      }

      return(res)

    })

    names(eigen_res) <- names(modData)

    attr(eigen_res, "group_var") <- attr(modData, "group_var")

  }else{

    mods <- modData

    abunData <- omicsData$e_data

    pcas <- lapply(unique(mods$Module), function(y){
      # Subset module data down to specific module
      module <- subset(mods, Module == y)

      # Merge module information with abundance information and format
      module <- merge(module, abunData, by=attr(modData, "cnames")$edata_cname, all.x=TRUE)
      rownames(module) <- module[,which(colnames(module) == attr(modData, "cnames")$edata_cname)]
      module <- module[,-which(colnames(module) %in% c(attr(modData, "cnames")$edata_cname, "Module"))]
      module <- t(module)

      # pca <- dudi.pca(module, center=FALSE, scale=FALSE, scannf=FALSE, nf=2)

      # Use pca to calculate principal component axes and format
      pca <- pcaMethods::pca(object=t(module), nPcs=2, scale="none", center=FALSE, method=pca.method)
      pca <- data.frame(Samples=rownames(pca@loadings), pca@loadings)
      colnames(pca) <- c(attr(modData, "cnames")$fdata_cname, paste(colnames(pca)[2],"_Module",y,sep=""), paste(colnames(pca)[3],"_Module",y,sep=""))
      return(pca)
    })

    # Combine data from all modules
    pcas <- Reduce(function(x, y) merge(x, y, by=attr(modData, "cnames")$fdata_cname, all=TRUE), pcas)

    # If envirnmental variables, correlate module PCAs to environmental variables
    if(!is.null(envVars)){
      # Format data
      pcas <- merge(pcas, omicsData$f_data[,which(colnames(omicsData$f_data) %in% c(attr(omicsData, "cnames")$fdata_cname, envVars))], by=attr(omicsData, "cnames")$fdata_cname)

      # Correlation test of environmental variables with module PCAs
      env.cor <- psych::corr.test(x=pcas[,grep("PC[12]_Module",colnames(pcas))], y=pcas[,which(colnames(pcas) %in% envVars)], method=cor.method, use=use, adjust=padjust)

      # Format correlation coefficient and p-value results
      env.r <- reshape2::melt(env.cor$r)
      colnames(env.r) <- c("Module","EnvVar","CorrCoeff")
      env.p <- reshape2::melt(env.cor$p)
      colnames(env.p) <- c("Module","EnvVar","p.value")
      env.cor <- merge(env.r, env.p, by=c("Module","EnvVar"))

      eigen_res <- list(pca=pcas, corr=env.cor)
    }else{
      eigen_res <- list(pca=pcas)
    }


  }

  attr(eigen_res, "group_DF") <- attr(omicsData, "group_DF")
  attr(eigen_res, "parameters") <- list(corrtest=cor.method, padjust=padjust, use=use)

  class(eigen_res) <- c("modEnv", class(eigen_res))

  return(eigen_res)


}
