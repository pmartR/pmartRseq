#' Apply a filter to a mintR S3 object
#'
#' This function takes a filter object of class "countFilter" or "sampleFilter", and applies the filter to a dataset of class \code{cDNAdata}, \code{gDNAdata}, or \code{rRNAdata}.
#'
#' @param omicsData an object of the class \code{cDNAdata}, \code{gDNAdata}, or \code{rRNAdata} usually created by \code{\link{as.cDNAdata}}, \code{\link{as.gDNAdata}}, or \code{\link{as.rRNAdata}}, respectively.
#' @param filter_object an object of the class "countFilter" or "sampleFilter", created by \code{count_based_filter} \code{sample_based_filter}.
#'
#' @return An object of the class \code{cDNAdata}, \code{gDNAdata}, or \code{rRNAdata}, with specified cname_ids, edata_cnames, and emeta_cnames filtered out of the appropriate datasets.
#'
#' @examples
#' library(mintJansson)
#' data(rRNA_data)
#'
#' to_filter <- count_based_filter(omicsData = rRNA_data, fn="mean")
#' rRNAdata2 <- applyFilt(filter_object = to_filter, omicsData = rRNA_data, upper_lim = 2)
#' print(str(attributes(rRNAdata2)$filters))
#'
#' to_filter2 <- count_based_filter(omicsData = rRNAdata, fn="max")
#' rRNAdata3 <- applyFilt(filter_object = to_filter2, omicsData = rRNAdata, upper_lim = 2)
#' print(str(attributes(rRNAdata3)$filters))
#'
#' @seealso \code{\link{count_based_filter}} \code{\link{sample_based_filter}}
#'
#' @author Allison Thompson
#'
#' @export

applyFilt <- function(filter_object, omicsData, ...) {

  # check that omicsData is of mintR S3 class#
  if (!(class(omicsData) %in% c("cDNAdata", "gDNAdata", "rRNAdata"))) stop("omicsData must be of class 'cDNAdata', 'gDNAdata', or 'rRNAdata'")

  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  fdata_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname

  # generate warnings if data type is not present but user asks to filter #
  if (!is.null(filter_object$filter_edata) & is.null(edata_cname)) warning("e_data identifier column specified in filter_object is not present in omicsData$e_data. Specified e_data features cannot be filtered from data.")

  if (!is.null(filter_object$filter_emeta) & is.null(emeta_cname)) warning("e_meta identifier column specified in filter_object is not present in omicsData$e_meta. Specified e_meta features cannot be filtered from data.")

  if (!is.null(filter_object$filter_samples) & is.null(fdata_cname)) warning("Sample identifier column specified in filter_object is not present in omicsData. Specified samples cannot be filtered from data.")


  UseMethod("applyFilt")
}



# function for countFilter
#' @export
#' @name applyFilt
#' @rdname applyFilt
#' @param num_samps for k over a filtering, the minimum number of samples that need to have at least \code{upper_lim} features
#' @param  upper_lim OTUs must have a max/mean/percent/nonmiss/sum number of counts above this threshold. OTUs with a count less than or equal to this number will be removed.
applyFilt.countFilter <- function(filter_object, omicsData, upper_lim=2, num_samps=NULL) {

  # Check if k/a fn
  if (attr(filter_object, "function") == "ka") {
    if (is.null(num_samps)) {
      num_samps = 2
      warning("Minimum number of samples wasn't given, defaulting to 2.")
    }
    if (!(class(num_samps) %in% c("numeric","integer")) | num_samps <= 0 | num_samps > length(attr(filter_object, "sample_names"))) {
      stop("num_samps must be an integer greater than 0 and less than or equal to the total number of samples in the dataset")
    }
    if (length(num_samps) != 1) stop("num_samps must be of length 1")
  }else{
    if (!is.null(num_samps)) {
      warning("num_samps provided but not used for this filter function")
    }
  }

  # check that upper_lim is numeric and >=1 #
  if (!(class(upper_lim) %in% c("numeric","integer")) | upper_lim < 0) stop("upper_lim must be an integer greater than or equal to 0")
  # check that upper_lim is of length 1 #
  if (length(upper_lim) != 1) stop("upper_lim must be of length 1")

  edata_cname <- attributes(omicsData)$cnames$edata_cname
  fn <- attr(filter_object, "function")

  if (fn == "percent") {
    upper_lim = upper_lim / 100
  }

  if (fn == "ka") {
    num_obs <- filter_object[,paste("NumSamples",num_samps,sep="_")]
  }else{
    num_obs <- filter_object[,paste(fn,"OTUs",sep="")]
  }

  # get indices for which ones don't meet the min requirement #
  inds <- which(num_obs <= upper_lim)

  if (length(inds) < 1) {
    filter.edata <- NULL
  }else{
    filter.edata <- omicsData$e_data[, which(names(omicsData$e_data) == edata_cname)][inds]
  }

  filter_object_new = list(edata_filt = filter.edata, emeta_filt = NULL, samples_filt = NULL)

  # call the function that does the filter application
  results_pieces <- applyFilt_worker(omicsData = omicsData, filter_object = filter_object_new)

  # return filtered data object #
  results <- omicsData
  results$e_data <- results_pieces$temp.pep2
  results$f_data <- results_pieces$temp.samp2
  results$e_meta <- results_pieces$temp.meta1

  # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure
  if (!is.null(attr(results, "group_DF"))) {
    results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
    attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
    attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_frac_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

    if (!is.null(results$e_meta)) {
      # number of unique proteins that map to a peptide in e_data #
      if (!is.null(emeta_cname)) {
        num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
      }else{num_emeta = NULL}
    }else{
      num_emeta = NULL
    }
    attr(results, "data_info")$num_emeta = num_emeta
    ## end of update attributes (7/11/2016 by KS)
  }else{
    # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
    attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
    attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_frac_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

    if (!is.null(results$e_meta)) {
      # number of unique proteins that map to a peptide in e_data #
      if (!is.null(emeta_cname)) {
        num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
      }else{num_emeta = NULL}
    }else{
      num_emeta = NULL
    }
    attr(results, "data_info")$num_emeta = num_emeta
    ## end of update attributes (7/11/2016 by KS)
  }

  # set attributes for which filters were run
  attr(results, "filters")$countFilter <- list(report_text = "", threshold = c(), filtered = c())
  attr(results, "filters")$countFilter$report_text <- paste("A ", fn,"-based filter was applied to the data, removing molecules that have a ", fn," count less than ", upper_lim, ". A total of ", length(filter.edata)," molecules were filtered out of the dataset by this filter.", sep="")
  attr(results, "filters")$countFilter$threshold <- upper_lim
  attr(results, "filters")$countFilter$filtered <- filter.edata

  return(results)
}



# function for countFilter
#' @export
#' @name applyFilt
#' @rdname applyFilt
#' @param  upper_lim Samples must have a sum number of OTU reads above this threshold. Samples with a sum less than or equal to this number will be removed.
applyFilt.sampleFilter <- function(filter_object, omicsData, upper_lim=2) {

  # check that upper_lim is numeric and >=1 #
  if (!(class(upper_lim) %in% c("numeric","integer")) | upper_lim < 0) stop("upper_lim must be an integer greater than or equal to 0")
  # check that upper_lim is of length 1 #
  if (length(upper_lim) != 1) stop("upper_lim must be of length 1")

  edata_cname <- attr(omicsData, "cnames")$edata_cname
  fn <- attr(filter_object, "function")

  num_obs <- filter_object[,paste(fn,"Samps",sep="")]

  # get indices for which ones don't meet the min requirement #
  inds <- filter_object[which(num_obs <= upper_lim), "Sample"]

  if (length(inds) < 1) {
    filter.samples <- NULL
  }else{
    filter.samples <- inds
    #filter.samples <- omicsData$e_data[,c(which(colnames(omicsData$e_data) == edata_cname), which(colnames(omicsData$e_data) %in% as.character(inds)))]
    #filter.samples <- omicsData$e_data[, which(names(omicsData$f_data) == fdata_cname)][inds]
  }

  filter_object_new = list(edata_filt = NULL, emeta_filt = NULL, samples_filt = filter.samples)

  # call the function that does the filter application
  results_pieces <- applyFilt_worker(omicsData = omicsData, filter_object = filter_object_new)

  # return filtered data object #
  results <- omicsData
  results$e_data <- results_pieces$temp.pep2
  results$f_data <- results_pieces$temp.samp2
  results$e_meta <- results_pieces$temp.meta1

  # if group attribute is present, re-run group_designation in case filtering any items impacted the group structure
  if (!is.null(attr(results, "group_DF"))) {
    results <- group_designation(omicsData = results, main_effects = attr(attr(omicsData, "group_DF"), "main_effects"), covariates = attr(attr(omicsData, "group_DF"), "covariates"), time_course = attr(attr(omicsData, "group_DF"), "time_course"))
    # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
    attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
    attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_frac_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

    if (!is.null(results$e_meta)) {
      # number of unique proteins that map to a peptide in e_data #
      if (!is.null(emeta_cname)) {
        num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
      }else{num_emeta = NULL}
    }else{
      num_emeta = NULL
    }
    attr(results, "data_info")$num_emeta = num_emeta
    ## end of update attributes (7/11/2016 by KS)
  }else{
    # Update attributes (7/11/2016 by KS) - this is being done already in group_designation
    attributes(results)$data_info$num_edata = length(unique(results$e_data[, edata_cname]))
    attributes(results)$data_info$num_miss_obs = sum(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_frac_missing = mean(is.na(results$e_data[,-which(names(results$e_data)==edata_cname)]))
    attributes(results)$data_info$num_samps = ncol(results$e_data) - 1

    if (!is.null(results$e_meta)) {
      # number of unique proteins that map to a peptide in e_data #
      if (!is.null(emeta_cname)) {
        num_emeta = length(unique(results$e_meta[which(as.character(results$e_meta[, edata_cname]) %in% as.character(results$e_data[, edata_cname])), emeta_cname]))
      }else{num_emeta = NULL}
    }else{
      num_emeta = NULL
    }
    attr(results, "data_info")$num_emeta = num_emeta
    ## end of update attributes (7/11/2016 by KS)
  }

  # set attributes for which filters were run
  attr(results, "filters")$sampleFilter <- list(report_text = "", threshold = c(), filtered = c())
  attr(results, "filters")$sampleFilter$report_text <- paste("A ", fn,"-based filter was applied to the data, removing samples that have a ", fn," count less than ", upper_lim, ". A total of ", length(filter.samples)," samples were filtered out of the dataset by this filter.", sep="")
  attr(results, "filters")$sampleFilter$threshold <- upper_lim
  attr(results, "filters")$sampleFilter$filtered <- filter.samples

  return(results)
}


#' Remove items that need to be filtered out
#'
#' This function removes filtered objects
#'
#' @param omicsData an object of the class \code{cDNAdata}, \code{gDNAdata}, or \code{rRNAdata} usually created by \code{\link{as.cDNAdata}}, \code{\link{as.gDNAdata}}, \code{\link{as.rRNAdata}}, respectively.
#' @param filter_object a list created by the functions above
#' @return list similar to as.mintRData object
#' @author Kelly Stratton, Lisa Bramer, Allison Thompson
#'
applyFilt_worker <- function(filter_object, omicsData) {
  # pull column names from omicR_data attributes #
  col_nms = attr(omicsData, "cnames")
  fdata_cname = col_nms$fdata_cname
  edata_cname = col_nms$edata_cname
  emeta_cname = col_nms$emeta_cname

  # pull group_DF attribute #
  group_DF = attr(omicsData, "group_DF")

  ## check to see if e_meta is provided ##

  # if not provided we only need to worry about e_data and f_data #

  if (attr(omicsData, "meta_info") == FALSE) {

    ## remove entries in edata ##
    if (!is.null(filter_object$edata_filt) & !is.null(edata_cname)) {

      temp.pep = omicsData$e_data

      # have to check that at least one of the items is present in the data #
      edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_filt)

      if (length(edat_ids) > 0) {
        # identify which peptides in e_data match filter list and remove #
        temp.pep1 = temp.pep[-which(temp.pep[,edata_cname] %in% filter_object$edata_filt),]
      }else{temp.pep1 = temp.pep}

    }else{ # no entries in edata need to be removed
      temp.pep1 = omicsData$e_data
    }

    ## remove samples ##
    if (!is.null(filter_object$samples_filt) & !is.null(fdata_cname)) {
      # identify which samples in f_data match filter list #
      temp.samp = omicsData$f_data

      # check that at least one sample is in f_data and e_data #
      fdat_ids = which(temp.samp[,fdata_cname] %in% filter_object$samples_filt)
      edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_filt)

      if (length(fdat_ids) > 0) {
        temp.samp2 = temp.samp[-which(temp.samp[,fdata_cname] %in% filter_object$samples_filt),]
      }else{temp.samp2 = temp.samp}

      # identify which samples in e_data match filter list and remove #
      if (length(edat_ids2) > 0) {
        temp.pep2 = temp.pep1[, -which(names(temp.pep1) %in% filter_object$samples_filt)]
      }else{temp.pep2 = temp.pep1}

    }else{ # no entries in f_data need to be removed
      temp.samp2 = omicsData$f_data
      temp.pep2 = temp.pep1
    }

    temp.meta2 = NULL




  }else{ # e_meta is present, so we need to work with it as well
    ## remove entries in edata ##
    if (!is.null(filter_object$edata_filt) & !is.null(edata_cname)) {

      temp.pep = omicsData$e_data

      # have to check that at least one of the items is present in the data #
      edat_ids = which(temp.pep[,edata_cname] %in% filter_object$edata_filt)

      if (length(edat_ids) > 0) {
        # identify which peptides in e_data and e_meta match filter list and remove#
        temp.pep1 = temp.pep[-which(temp.pep[,edata_cname] %in% filter_object$edata_filt),]
      }else{temp.pep1 = temp.pep}

      temp.meta = omicsData$e_meta

      # check that at least one of the peptides is present in e_meta #
      emeta_ids = which(temp.meta[,edata_cname] %in% filter_object$edata_filt)

      if (length(emeta_ids) > 0) {
        temp.meta1 = temp.meta[-which(temp.meta[,edata_cname] %in% filter_object$edata_filt),]
      }else{temp.meta1 = temp.meta}

    }else{
      temp.pep1 = omicsData$e_data
      temp.meta1 = omicsData$e_meta
    }

    ## remove samples ##
    if (!is.null(filter_object$samples_filt) & !is.null(fdata_cname)) {
      # identify which samples in f_data match filter list #
      temp.samp = omicsData$f_data

      # check that at least one sample is in f_data and e_data #
      fdat_ids = which(temp.samp[,fdata_cname] %in% filter_object$samples_filt)
      edat_ids2 = which(names(temp.pep1) %in% filter_object$samples_filt)

      if (length(fdat_ids) > 0) {
        temp.samp2 = temp.samp[-which(temp.samp[,fdata_cname] %in% filter_object$samples_filt),]
      }else{temp.samp2 = temp.samp}

      # identify which samples in e_data match filter list and remove #
      if (length(edat_ids2) > 0) {
        temp.pep2 = temp.pep1[, -which(names(temp.pep1) %in% filter_object$samples_filt)]
      }else{temp.pep2 = temp.pep1}

    }else{
      temp.samp2 = omicsData$f_data
      temp.pep2 = temp.pep1
    }

    ## remove entries in emeta ##
    if (!is.null(filter_object$emeta_filt) & !is.null(emeta_cname)) {
      # identify which proteins in data match filter list and remove from e_meta #
      temp.meta = temp.meta1

      # check that at least one of the proteins is in e_meta #
      emeta_ids2 = which(temp.meta[,emeta_cname] %in% filter_object$emeta_filt)

      if (length(emeta_ids2) > 0) {
        temp.meta2 = temp.meta[-which(temp.meta[,emeta_cname] %in% filter_object$emeta_filt),]
      }else{temp.meta2 = temp.meta}
    }else{
      temp.meta2 = temp.meta1
    }


    # check for rogue entries in edata #
    edat_ids2 = which(!(temp.pep2[,edata_cname] %in% temp.meta1[,edata_cname]))

    # filter out edata entries which no longer have mappings to emeta entries #
    if (length(edat_ids2) > 0) {
      temp.pep2 = temp.pep2[-which(!(temp.pep2[,edata_cname] %in% temp.meta1[,edata_cname])),]
    }



  }

  output <- list(temp.pep2 = temp.pep2, temp.samp2 = temp.samp2, temp.meta1 = temp.meta2, edata_cname = edata_cname, emeta_cname = emeta_cname, fdata_cname = fdata_cname)

  # return the pieces needed to assemble a cDNAdata/gDNAdata/rRNAdata object
  return(output)
}
