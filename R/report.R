#' Creates a report for pmartRseq objects
#'
#' This function generates a .docx report for a pmartRseq object
#'
#' @param omicsData a list containing pmartRseq objects, including at least one data object (of the class 'seqData' created by \code{\link{as.seqData}} and any other pmartRseq objects to include in the report.
#' @param output_file an optional character string specifying the name of the file to create.
#'
#' @details This function generates a .docx report for a mintR data object and a mintR results object. The report includes information about the data, groups, filtering, normalization, and statistical tests, if relevant.
#'
#' @return A .docx report of the analyses performed on the data
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(cDNA_hiseq_data)
#' mycdnadata <- group_designation(omicsData = cDNA_hiseq_data, main_effects = c("treatment"), time_course=NULL)
#' mycdnadata_norm <- normalize_data(omicsData = mycdnadata, norm_fn = "percentile")
#' mycdnadata_results <- countSTAT(omicsData = mycdnadata_norm, comparisons = "all", control = NULL, test = c("dw", "eq", "el", "ef"), pval_adjust = "none", pval_thresh = 0.05)
#' report(omicsData = list(Norm=mycdnadata_norm, Statistics = mycdnadata_results), output_file = "cDNAdata_Report.docx")
#' }
#'
#' @author Allison Thompson
#'
#' @export
report <- function(omicsData, output_file=NULL){
  library(rmarkdown)

  if(!is.list(omicsData)){
    stop("omicsData must be a list of pmartRseq objects")
  }

  data <- omicsData
  classes <- unlist(lapply(data, class))

  params <- list(data=data, classes=classes)

  render("R/seqData_Report.Rmd", output_file=output_file, params=params)

}

