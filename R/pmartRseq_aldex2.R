#' Modified version of ALDEx2 for pmartRseq
#'
#' This function calculates the centered-log ratio with draws from a Dirichlet distribution and runs a liner or linear mixed-effects model
#'
#' @param omicsData an object of the class 'seqData' reated by \code{\link{as.seqData}}.
#' @param mainEffects Optional, a character vector detailing which factors should be used as main effects in the model. The variable name must match a column name from \code{omicsData$f_data}. If NULL, will use the same main effects that were used in \code{\link{group_designation}}.
#' @param interactions Logical, indicating whether or not interactions should be considered in the model. If TRUE, all possible interactions will be used. If FALSE, none will be used.
#' @param randomEffect Optional, character specifying which factor to use as a random effect in the model. This will be used in the form (1|randomEffect). Must match a column name from \code{omicsData$f_data}
#' @param mc.samples The number of Monte Carlo instances to use. Default is 128.
#' @param denom Character vector indicating which features to use as the denominator for the geometric mean. Default is "all".
#' @param verbose Print diagnostic information while running. Default is FALSE.
#'
#' @details Modified version of ALDEx2, allowing for multiple factors and random effects
#'
#' @return An object of class paRes (also a list) containing to parts - clr, the clr transformed values for each Monte-Carlo Dirichlet instance, and res, a data frame with a p.value, statistic, degrees of freedom, and sum of squares values for each feature and each term in the mdoel.
#'
#' @references Fernandes, Unifying...
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @author Allison Thompson
#'
#' @export

pmartRseq_aldex2 <- function(omicsData, mainEffects=NULL, mc.samples=128, denom="all", verbose=FALSE,
                             interactions=FALSE, randomEffect=NULL, pval_thresh){

  library(ALDEx2)
  library(gtools)
  library(broom)
  library(car)

  aldex.clr.function <- function(omicsData, mc.samples=128, denom="all", verbose=FALSE) {

    # INPUT
    # The 'reads' data.frame MUST have row
    # and column names that are unique, and
    # looks like the following:
    #
    #              T1a T1b  T2  T3  N1  N2
    #   Gene_00001   0   0   2   0   0   1
    #   Gene_00002  20   8  12   5  19  26
    #   Gene_00003   3   0   2   0   0   0
    #       ... many more rows ...
    #
    # ---------------------------------------------------------------------

    # OUTPUT
    # The output returned is a list (x) that contains Monte-Carlo instances of
    # the centre log-ratio transformed values for each sample
    # Access to values
    # sample IDs: names(x)
    # number of features (genes, OTUs): length(x[[1]][,1])
    # number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
    # feature names: rownames(x[[1]])



    reads <- omicsData$e_data
    rownames(reads) <- reads[,which(colnames(reads) == attr(omicsData, "cnames")$edata_cname)]
    reads <- reads[,-which(colnames(reads) == attr(omicsData, "cnames")$edata_cname)]



    # make sure that mc.samples is an integer, despite it being a numeric type value
    as.numeric(as.integer(mc.samples))

    #  remove all rows with reads less than the minimum set by minsum
    minsum <- 0

    # get groups
    groups <- attr(omicsData, "group_DF")$Group

    # remove any row in which the sum of the row is 0
    z <- as.numeric(apply(reads, 1, sum))
    reads <- as.data.frame( reads[(which(z > minsum)),]  )

    if (verbose) print("removed rows with sums equal to zero")


    #  SANITY CHECKS ON THE DATA INPUT
    if ( any( round(reads) != reads ) ) stop("not all reads are integers")
    if ( any( reads < 0 ) )             stop("one or more reads are negative")

    for ( col in names(reads) ) {
      if ( any( ! is.finite( reads[[col]] ) ) )  stop("one or more reads are not finite")
    }

    if ( length(rownames(reads)) == 0 ) stop("rownames(reads) cannot be empty")
    if ( length(colnames(reads)) == 0 ) stop("colnames(reads) cannot be empty")

    if ( length(rownames(reads)) != length(unique(rownames(reads))) ) stop ("row names are not unique")
    if ( length(colnames(reads)) != length(unique(colnames(reads))) ) stop ("col names are not unique")
    if ( mc.samples < 128 ) warning("values are unreliable when estimated with so few MC smps")

    # add a prior expection to all remaining reads that are 0
    # this should be by a Count Zero Multiplicative approach, but in practice
    # this is not necessary because of the large number of features
    prior <- 0.5

    # This extracts the set of features to be used in the geometric mean computation
    feature.subset <- aldex.set.mode(reads, groups, denom)

    reads <- reads + prior

    if (verbose == TRUE) print("data format is OK")

    # ---------------------------------------------------------------------
    # Generate a Monte Carlo instance of the frequencies of each sample via the Dirichlet distribution,
    # returns frequencies for each feature in each sample that are consistent with the
    # feature count observed as a proportion of the total counts per sample given
    # technical variation (i.e. proportions consistent with error observed when resequencing the same library)

    nr <- nrow( reads )
    rn <- rownames( reads )

    #this returns a list of proportions that are consistent with the number of reads per feature and the
    #total number of reads per sample

    # environment test, runs in multicore if possible
    # if (has.BiocParallel){
    #     p <- bplapply( reads ,
    #         function(col) {
    #             q <- t( rdirichlet( mc.samples, col ) ) ;
    #             rownames(q) <- rn ;
    #             q })
    #     names(p) <- names(reads)
    # }
    #else{
    p <- lapply( reads ,
                 function(col) {
                   q <- t( rdirichlet( mc.samples, col ) ) ;
                   rownames(q) <- rn ; q } )
    #}

    # sanity check on the data, should never fail
    for ( i in 1:length(p) ) {
      if ( any( ! is.finite( p[[i]] ) ) ) stop("non-finite frequencies estimated")
    }

    if (verbose == TRUE) print("dirichlet samples complete")

    # ---------------------------------------------------------------------
    # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
    # i.e., do a centered logratio transformation as per Aitchison

    # apply the function over elements in a list, that contains an array

    # DEFAULT
    if(length(feature.subset) == nr)
    {
      # Default ALDEx2
      # if (has.BiocParallel){
      #     l2p <- bplapply( p, function(m) {
      #         apply( log2(m), 2, function(col) { col - mean(col) } )
      #     })
      #     names(l2p) <- names(p)
      # }
      #else{
      l2p <- lapply( p, function(m) {
        apply( log2(m), 2, function(col) { col - mean(col) } )
      })
      #}
    } else {
      ## IQLR or ZERO
      feat.result <- vector("list", length(unique(groups))) # Feature Gmeans
      condition.list <- vector("list", length(unique(groups)))    # list to store conditions

      for (i in 1:length(unique(groups)))
      {
        condition.list[[i]] <- which(groups == unique(groups)[i]) # Condition list
        feat.result[[i]] <- lapply( p[condition.list[[i]]], function(m) {
          apply(log2(m), 2, function(x){mean(x[feature.subset[[i]]])})
        })
      }
      set.rev <- unlist(feat.result, recursive=FALSE) # Unlist once to aggregate samples
      p.copy <- p
      for (i in 1:length(set.rev))
      {
        p.copy[[i]] <- as.data.frame(p.copy[[i]])
        p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
        p[[i]] <- t(p[[i]])
      }
      l2p <- p    # Save the set in order to generate the aldex.clr variable
    }


    # sanity check on data
    for ( i in 1:length(l2p) ) {
      if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
    }
    if (verbose == TRUE) print("clr transformation complete")

    return(new("aldex.clr",reads=reads,mc.samples=mc.samples,verbose=verbose,analysisData=l2p))
  }

  ###########################
  # a first run at a glm function for ALDEx2 output

  # April 14th, 2014

  # Initially, I wrote the function to extract model coefficients and 95%CI, but decided that
  # we would probably want to do the pairwise comparisons (with proper p-value correction)
  # for the features that have significant p-values by glm or Kruskal-Wallis. If people are interested
  # in the coefficients, I can reinsert them, but it ups computation time.

  # There is probably room for optimization to speed things up, but for the moment this works


  ##################
  # for glms
  ##################

  # Data structure returned by clr_test_function.r function
  # The output returned is a list (x) that contains Monte-Carlo instances of
  # the centre log-ratio transformed values for each sample
  # sample IDs: names(x)
  # number of features: length(x[[1]][,1])
  # number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
  # feature names: rownames(x[[1]])

  # INVOCATION
  # conditions is using selex dataset
  # conditions <- c(rep("N", 7), rep("S",7)
  # x.glm <- aldex.glm(x,conditions)

  #returns a dataframe of expected P and fdr statistics for each feature

  aldex.glm <- function(clr, mainEffects, interactions=FALSE, randomEffect=NULL){
    #library(nlme)
    #library(papeR)
    # make sure that the multicore package is in scope and return if available
    # is.multicore = FALSE
    #
    # if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    #     print("multicore environment is OK -- using the BiocParallel package")
    #     #require(BiocParallel)
    #     is.multicore = TRUE
    # }
    # else {
    #     print("operating in serial mode")
    # }

    # get dimensions, names, etc from the input data
    smpl.ids <- getSampleIDs(clr)
    feature.number <- numFeatures(clr)
    mc.instances <- numMCInstances(clr)
    feature.names <- getFeatureNames(clr)

    #conditions <- as.factor( conditions )
    #levels     <- levels( conditions )
    #conditions <- apply(conds, 2, as.factor)
    conditions <- omicsData$f_data[,mainEffects]
    #names(conditions) <- colnames(conds)
    #conditions <- apply(conditions, 2, as.factor)
    # site <- as.factor(site)
    # crop <- as.factor(crop)
    # agg <- as.factor(agg)

    levels <- apply(conditions, 2, levels)
    # levels.site <- levels(site)
    # levels.crop <- levels(crop)
    # levels.agg <- levels(agg)

    invisible(apply(conditions, 2, function(x) if(length(x) != numConditions(clr)) stop("mismatch between conds and names(clr)")))

    # if ( length( site ) !=  numConditions(clr) )  stop("mismatch btw 'length(site)' and 'length(names(clr))'")
    # if ( length( crop ) !=  numConditions(clr) )  stop("mismatch btw 'length(crop)' and 'length(names(clr))'")
    # if ( length( agg ) !=  numConditions(clr) )  stop("mismatch btw 'length(agg)' and 'length(names(clr))'")
    #
    #   levels.site <- vector( "list", length( levels.site ) )
    #   names( levels.site ) <- levels( site )
    #   sets.site <- names(levels.site)
    #   lsets.site <- length(sets.site)
    #
    #   levels.crop <- vector( "list", length( levels.crop ) )
    #   names( levels.crop ) <- levels( crop )
    #   sets.crop <- names(levels.crop)
    #   lsets.crop <- length(sets.crop)
    #
    #   levels.agg <- vector( "list", length( levels.agg ) )
    #   names( levels.agg ) <- levels( agg )
    #   sets.agg <- names(levels.agg)
    #   lsets.agg <- length(sets.agg)

    # set up the glm results containers

    n <- 2^(length(mainEffects)) - 1

    glm.matrices <- list()
    # p-values and BH p-values
    # glm.matrix.site <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
    # glm.matrix.crop <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
    # glm.matrix.agg <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
    # glm.matrix.inter <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
    #glm.matrix.pBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)

  #  glm.estes <- lapply(c(1:n), function(x) matrix(data = NA, nrow = feature.number, ncol = mc.instances))
    # glm.matrix.siteBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
    # glm.matrix.cropBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)

    # for Kruskal Wallis p
    # kw.p.matrix = matrix(data = NA, nrow = feature.number, ncol = mc.instances)
    # kw.pBH.matrix = matrix(data = NA, nrow = feature.number, ncol = mc.instances)

    #########
    # this is where optimization would probably be helpful to possibly avoid this for loop
    #########
    #mc.i is the monte carlo instance
    for(mc.i in 1:mc.instances){
      #print(mc.i)

      #generate a matrix of each Monte-Carlo instance, columns are samples, rows are features
      t.input <- sapply(getMonteCarloInstances(clr), function(y){y[,mc.i]})

      # do glms on each feature
      # make a list of glm outputs
      # new <- apply(t.input, 1, function(yy) {
      #   tmp <- data.frame(yy=yy, site=site, crop=crop, agg=agg, block=block)
      #   lmer(as.numeric(yy) ~ factor(site) + factor(crop) + factor(agg) + factor(site)*factor(crop) + (1|block), data=tmp)
      # })

      x <- apply(t.input, 1, function(yy) {

        if(!is.null(randomEffect)){
          tmp <- data.frame(Samples=names(yy), yy=yy)
          colnames(tmp)[1] <- attr(omicsData, "cnames")$fdata_cname
          tmp <- merge(tmp, omicsData$f_data[,which(colnames(attr(omicsData,"group_DF")) %in% c(mainEffects, randomEffect, attr(omicsData,"cnames")$fdata_cname))], by=attr(omicsData, "cnames")$fdata_cname)
          rownames(tmp) <- tmp[,which(colnames(tmp) == attr(omicsData, "cnames")$fdata_cname)]
          tmp <- tmp[,-which(colnames(tmp) == attr(omicsData, "cnames")$fdata_cname)]

          tmp$yy <- as.numeric(tmp$yy)

          if(interactions){
            form <- as.formula(paste("yy ~ ", paste(colnames(tmp)[-which(colnames(tmp) %in% c("yy",randomEffect))], collapse="*"), "+ (1|",randomEffect,")", sep=""))
            lmer(form, data=tmp, control=lmerControl(optimizer='Nelder_Mead'))
          }else{
            form <- as.formula(paste("yy ~ ", paste(colnames(tmp)[-which(colnames(tmp) %in% c("yy",randomEffect))], collapse="+"), "+ (1|",randomEffect,")", sep=""))
            lmer(form, data=tmp, control=lmerControl(optimizer='Nelder_Mead'))
          }
        }else{
          tmp <- data.frame(Samples=names(yy), yy=yy)
          colnames(tmp)[1] <- attr(omicsData, "cnames")$fdata_cname
          tmp <- merge(tmp, omicsData$f_data[,which(colnames(attr(omicsData,"group_DF")) %in% c(mainEffects, attr(omicsData,"cnames")$fdata_cname))], by=attr(omicsData, "cnames")$fdata_cname)
          rownames(tmp) <- tmp[,which(colnames(tmp) == attr(omicsData, "cnames")$fdata_cname)]
          tmp <- tmp[,-which(colnames(tmp) == attr(omicsData, "cnames")$fdata_cname)]

          tmp$yy <- as.numeric(tmp$yy)
          if(interactions){
            form <- as.formula(paste("yy ~ ", paste(colnames(tmp)[-which(colnames(tmp) %in% c("yy"))], collapse="*"), sep=""))
            lm(form, data=tmp)
          }else{
            form <- as.formula(paste("yy ~ ", paste(colnames(tmp)[-which(colnames(tmp) %in% c("yy"))], collapse="+"), sep=""))
            lm(form, data=tmp)
          }
        }


        # if(!is.null(random)){
        #   form <- as.numeric(yy) ~ .^(length(mainEffects))
        #   lme(form, random=~1|as.symbol(random), data=tmp)
        # }else{
        #   lme(as.numeric(yy) ~ as.symbol(paste(sapply(names(conditions), function(x) paste("factor(",x,")",sep="")),collapse="+")) +
        #         as.symbol(paste(sapply(names(conditions), function(x) paste("factor(",x,")",sep="")),collapse="*")),
        #       data=tmp)
        # }
      })

      # newx2 <- apply(t.input, 1, function(yy) {
      #   tmp <- data.frame(yy=yy, site=site, crop=crop, agg=agg, block=block)
      #   Anova(lme(as.numeric(yy) ~ factor(site) + factor(crop) + factor(agg) + factor(site)*factor(crop), random=~1|block, data=tmp),type="sequential")
      # })

      # x <- apply(t.input, 1, function(yy) {
      #   #tmp <- data.frame(yy=yy, site=site, crop=crop, agg=agg, block=block)
      #   lme(as.numeric(yy) ~ factor(site) + factor(crop) + factor(agg) + factor(site)*factor(crop), random=~1|block)
      # })
      # p-valuess
      # if (is.multicore == TRUE){
      #   pps <- bplapply(x, drop1, test = "Chis"  )
      # } else {
        #pps <- lapply(x, drop1, test = "Chis")
        # pps = lapply(x, function(x) summary(x)$tTable[2,5])
        # sitep  <- lapply(x, function(x) Anova(x)[2,4])
        # cropp <- lapply(x, function(x) Anova(x)[3,4])
        # aggp <- lapply(x, function(x) Anova(x)[4,4])
        # interp <- lapply(x, function(x) Anova(x)[5,4])
        glm.p <- lapply(c(1:length(x)), function(y){
          temp <- broom::tidy(Anova(x[[y]], type="III"))
          temp <- data.frame(Feature=rownames(t.input)[y], temp)
          colnames(temp)[3:5] <- paste(colnames(temp)[3:5], mc.i, sep=".")
          return(temp)
        })
        # pests = lapply(x, function(x) summary(x)$coeff$fixed[2])
        # site.ests <- lapply(x, function(x) summary(x)$coeff$fixed[2])
        # crop.ests <- lapply(x, function(x) summary(x)$coeff$fixed[3])
      #}
      # #glm.matrix.p[, mc.i] <- sapply(pps, function(x){x[[5]][2]})
      # glm.matrix.site[,mc.i] = unlist(sitep)
      # glm.matrix.crop[,mc.i] = unlist(cropp)
      # glm.matrix.agg[,mc.i] = unlist(aggp)
      # glm.matrix.inter[,mc.i] = unlist(interp)
      # # glm.matrix.pBH[, mc.i] <- unlist(pests)
      # glm.matrix.siteBH[,mc.i] = unlist(site.ests)
      # glm.matrix.cropBH[,mc.i] = unlist(crop.ests)

      glm.matrices[[mc.i]] <- do.call(rbind, glm.p)
      # cat(mc.i,"\n")
      # # Kruskal Wallis
      # kw.p.matrix[, mc.i] <- t(apply(t.input, 1, function(yy){
      #   kruskal.test(yy ~ factor(site) + factor(crop) + factor(agg) + factor(site)*factor(crop))[[3]]
      # }))
      # kw.pBH.matrix[, mc.i] <- as.numeric(p.adjust(kw.p.matrix[, mc.i], method = "BH"))

    }

    glm.res <- Reduce(function(x, y) merge(x, y, by=c("Feature", "term")), glm.matrices)

    glm.p.res <- apply(glm.res[,grep("p.value",colnames(glm.res))], 1, function(x) quantile(x, 0.5, na.rm=TRUE))
    glm.s.res <- apply(glm.res[,grep("statistic",colnames(glm.res))], 1, function(x) quantile(x, 0.5, na.rm=TRUE))
    glm.df.res <- apply(glm.res[,grep("df",colnames(glm.res))], 1, function(x) quantile(x, 0.5, na.rm=TRUE))
    glm.sumsq.res <- apply(glm.res[,grep("sumsq",colnames(glm.res))], 1, function(x) quantile(x, 0.5, na.rm=TRUE))

    glm.tot <- data.frame(glm.res[,c(1,2)], p.value=glm.p.res, statistic=glm.s.res, df=glm.df.res, sumsq=glm.sumsq.res)

    glm.tot <- glm.tot[-which(glm.tot$term == "Residuals"),]
    #glm.tot <- glm.tot[-which(glm.tot$term == "(Intercept)"),]
    rownames(glm.tot) <- c(1:nrow(glm.tot))

    return(glm.tot)

  }

  aldexclr <- aldex.clr.function(omicsData=omicsData, mc.samples=mc.samples, denom=denom, verbose=verbose)

  if(is.null(mainEffects)){
    mainEffects <- attr(attr(omicsData, "group_DF"), "main_effects")
  }

  if(is.null(mainEffects)){
    stop("Main effects were not given and cannot be found from the group data frame. Please run group designation function with desired main effects.")
  }

  aldexres <- aldex.glm(clr=aldexclr, mainEffects=mainEffects, interactions=interactions, randomEffect=randomEffect)

  results <- list(clr=aldexclr, results=aldexres)

  attr(results, "effects") <- list(mainEffects=mainEffects, interactions=ifelse(interactions, "All", "None"), randomEffect=randomEffect)

  attr(results, "pval_thresh") <- list(pval_thresh=pval_thresh)

  class(results) <- "paRes"

  return(results)



}
