#' Calculates Richness
#'
#' This function calculates the unique number of features seen in each sample.
#'
#' @param omicsData an object of the class 'seqData' created by \code{\link{as.seqData}}.
#' @param index a character vector stating which of the calculations to perform - "observed" for the observed richness, "chao1" the bias-corrected chao1 richness estimator, and/or "ace" for the abundance-based coverage richness estimator. Default is to perform all 3 calculations.
#'
#' @details Calculates richness of count data
#'
#' @return An object of class richRes (also a data.frame) containing the richness value for every sample in the data object.
#' @references Chao, Anne. Species Richness Estimation.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#' rRNA_richness <- richness_calc(omicsData = rRNA_data)
#' rRNA_richness
#' summary(rRNA_richness)
#' plot(rRNA_richness)
#' }
#'
#' @author Allison Thompson
#'
#' @export

richness_calc <- function(omicsData, index=c("observed","chao1","ace","break")){

  ## some initial checks ##

  # check that omicsData is of appropriate class #
  if(!class(omicsData) %in% c("seqData")) stop("omicsData must be of class 'seqData'")

  if(attr(omicsData, "data_info")$data_scale!='count'){
    warning("This function is meant for count data like 'rRNA', 'gDNA' or 'cDNA' data.")
  }

  if(attr(omicsData, "data_info")$data_norm){
    warning("This function works best on non-normalized data.")
  }

  if(any(!index %in% c("observed","chao1","ace","break"))){
    stop("index must only include 'observed', 'chao1', 'ace', or 'break'")
  }

  ## end initial checks ##

  # change 0 to NA, makes for easier calculation
  omicsData$e_data[omicsData$e_data == 0] <- NA

  edata_cname <- attr(omicsData, "cnames")$edata_cname

  res <- list()

  if("observed" %in% tolower(index)){
    # calculate observed richness
    observed <- apply(omicsData$e_data[,-which(colnames(omicsData$e_data)==edata_cname)], 2, function(y) length(which(!is.na(y))))
    res[["observed"]] <- observed
  }

  if("chao1" %in% tolower(index)){
    chao1 <- apply(omicsData$e_data[,-1], 2, function(x){

      # Chao1 needs un-normalized counts, uses the number of singletons and doubletons
      if(attr(omicsData, "data_info")$data_norm){
        warning("Chao1 calculation may be incorrect when computed on normalized data. Use raw counts for this calculation.")
      }

      # Remove any NA's and 0's
      if(any(is.na(x))){
        x <- x[-which(is.na(x))]
      }

      if(any(x==0)){
        x <- x[-which(x==0)]
      }

      # Calculate observed richness
      s <- length(x)

      # Calculate number of singletons
      one <- length(which(x==1))
      # Calculate number of doubletons
      two <- length(which(x==2))
      #s + (one ^ 2)/(2 * two)
      # Chao1's bias corrected estimator
      s + ((one * (one - 1)) / (2 * (two + 1)))
    })
    res[["chao1"]] <- chao1
    #attr(results, "indexValueRange")$chao1 <- data.frame(Min=min(res[["chao1"]], na.rm=TRUE), Max=max(res[["chao1"]], na.rm=TRUE))
  }

  if("ace" %in% tolower(index)){
    # Change NA to 0 and extract e_data - change edata_cname column to rownames and remove from data
    temp <- omicsData$e_data
    temp[temp == 0] <- NA
    e_data <- temp[,-which(colnames(omicsData$e_data) == attr(omicsData, "cnames")$edata_cname)]
    rownames(e_data) <- omicsData$e_data[,which(colnames(omicsData$e_data) == attr(omicsData, "cnames")$edata_cname)]

    ace <- apply(e_data, 2, function(x){
      data <- x
      data[is.na(data)] <- 0
      # Only interested in the values that are above 0
      data <- data[data > 0]

      # Calculate total observed richness
      S_obs <- length(data)
      # Calculate number of rare species (counts less than or equal to 10)
      S_rare <- length(which(data <= 10))
      # Calculate number of abundant species (counts more than 10)
      S_abund <- length(which(data > 10))
      # Total sum of all species with counts less than 11
      N_rare <- sum(data[data < 11])

      # Calculate how many species are seen once, twice, etc
      a <- unlist(lapply(c(1:10), function(x) length(which(data==x))))
      # Uses singletons
      C_ace <- 1 - a[1]/N_rare

      # Calculate gamma^2
      gamma2 <- max(S_rare * sum(unlist(lapply(c(1:10), function(x) x*(x-1)*a[x]))) /
                      (C_ace * (sum(unlist(lapply(c(1:10), function(y) y*a[y])))^2)) - 1, 0)

      S_P <- S_abund + (S_rare + a[1]*gamma2)/C_ace
      return(S_P)
    })
    #S_P = S_abund + S_rare/C_ace + a1/C_ace * gamma^2
    #C_ace = 1-a1/N_rare
    #gamma^2 = max(S_rare/C_ace(sum[i=1:10]i*(i-1)*a_i)/N_rare/(N_rare-1)-1,0)
    #a_i refers to number of species with abundance i
    #S_rare is the number of rare species
    #S_abund is the number of abundant species
    #N_rare is the number of individuals in rare species

    res[["ace"]] <- ace
  }

  if("break" %in% tolower(index)){
    library(breakaway)
    # Format data
    e_data <- omicsData$e_data[,-1]

    # Run breakaway function and only pull out the estimated richness value
    chats <- apply(e_data, 2, function(x) breakaway(as.data.frame(table(x)), print=FALSE, plot=FALSE, answers=TRUE, force=FALSE)$est)
    res[["break"]] <- chats
  }

  richness <- as.data.frame(do.call(rbind, res))

  # make a richness object
  attr(richness, "group_DF") <- attr(omicsData, "group_DF")
  attr(richness, "cnames") <- attr(omicsData, "cnames")
  attr(richness, "index") <- index
  class(richness) <- c("richRes", class(richness))

  return(richness)

}
