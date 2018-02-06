#' Import rRNA (16S/ITS/18S), metatranscript, or metagenomic data from a .biom,
#' .csv, or .txt file. Converts import to a list object to be passed to
#' \link{as.seqData}.
#'
#' @param e_data_filepath character string specifying the file path to a .csv,
#'   .txt, or .biom expression data file. The expression data file should
#'   contain a feature identifier/name.
#' @param f_data_filepath character string specifying the file path to a .csv or
#'   .txt sample information file. Each row must correspond to a sample with one
#'   column giving the unique sample identifiers found in e_data column names
#'   and other columns providing qualitative and/or quantitative traits of each
#'   sample. The sample information should also contain a sample identifier
#'   column.
#' @param e_meta_filepath optional character string specifying the file path of
#'   expression meta information. The meta information must contain an
#'   identifier column that matches the identifier column in the expression data
#'   file.
#'
#' @author Sarah Reehl
#' @seealso \code{\link{as.seqData}}
#'
#' @export
import_seqData <- function(e_data_filepath, f_data_filepath, e_meta_filepath = NULL){
  # initial checks
  if (!(grepl(pattern = "\\.biom$", x = e_data_filepath) | grepl(pattern = "\\.txt$", x = e_data_filepath) | grepl(pattern = "\\.csv$", x = e_data_filepath))){
    stop("Unsupported biom format. Must end in .biom, .csv, or .txt")
  }
  if (!(grepl(pattern = "\\.txt$", x = f_data_filepath) | grepl(pattern = "\\.csv$", x = f_data_filepath))){
    stop("Unsupported sample information format. Must end in .csv or .txt")
  }
  # check class structures for paths and attempt read #
  if (inherits(e_data_filepath, "character")) {
    biom_read_attempt <- try({
      # check if biom is a .biom file
      if (grepl(pattern = "\\.biom$", x = e_data_filepath)) {
        biom_read <- biomformat::read_biom(biom_file = e_data_filepath)
        otu_table <- as.matrix(biomformat::biom_data(biom_read))
        # create e_data with OTU identifier
        e_data = data.frame(OTU = rownames(otu_table), otu_table, stringsAsFactors = FALSE, check.names = FALSE)
        edata_cname <- colnames(e_data)[1]
        # create e_meta with OTU identifier
        otu_meta = biomformat::observation_metadata(biom_read)
        e_meta = data.frame(OTU = row.names(otu_meta), otu_meta, stringsAsFactors = FALSE, check.names = FALSE)
        taxa_cname <- colnames(e_meta)[2]
      }
      # check if e_data is a csv or txt
      if (grepl(pattern = "\\.csv$", x = e_data_filepath) | grepl(pattern = "\\.txt$", x = e_data_filepath)) {
        # create e_data
        e_data <- data.table::fread(e_data_filepath, data.table = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
        edata_cname <- colnames(e_data)[1]
        # create e_meta
        if (!is.null(e_meta_filepath) & inherits(e_meta_filepath, "character")) {
          e_meta <- data.table::fread(e_meta_filepath, data.table = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
          taxa_cname <- colnames(e_meta)[2]
        }
        if (is.null(e_meta_filepath)) e_meta <- NULL
      }
    }, silent = TRUE)

    if (inherits(biom_read_attempt, "try-error")) {
      stop("Biom import failed. Incorrect file path or unsupported biom format.")
    }
  }
  # try to read f_data
  if (inherits(f_data_filepath, "character")) {
    meta_read_attempt <- try({
      # check if f_data is a .txt
      if (grepl(pattern = "\\.txt$", x = f_data_filepath)){
        f = readLines(f_data_filepath)
        #find the last commented line and assume header info
        skipLines = which.max(grepl("#", x = f[1:length(f)]))
        f_data <- data.table::fread(input=paste0(f, collapse = "\n"), sep = "\t", header = TRUE, skip = skipLines - 1, data.table = FALSE, check.names = FALSE)
        # choose the sample identifier column (first for now)
        #fdata_cname = colnames(f_data)[1]
        # choose the sample identifier column by matching column names in e_data with a row in f_data
        ids <- table(do.call(rbind, lapply(colnames(e_data)[-which(colnames(e_data)==edata_cname)], function(x) which(f_data == x, arr.ind=TRUE)))[,2])
        ids <- ids[which(ids == max(ids))]
        ids <- names(ids)[1]
        fdata_cname <- colnames(f_data)[as.numeric(ids)]
      }
      # check if f_data is a .csv
      if (grepl(pattern = "\\.csv$", x = f_data_filepath)){
        # normal read assume no comments
        f_data <- data.table::fread(f_data_filepath, data.table = FALSE, check.names = FALSE)
        # choose the sample identifier column (first for now)
        fdata_cname = colnames(f_data)[1]
      }
    }, silent = TRUE)
  }
  if (inherits(meta_read_attempt, "try-error")) {
    stop("Sample metadata import failed. Incorrect file path or unsupported file format.")
  }

  imported_list <- list(e_data = e_data, e_meta = e_meta, f_data = f_data, guessed_edata_cname = edata_cname, guessed_fdata_cname = fdata_cname, guessed_taxa_cname = taxa_cname)
  return(imported_list)
}
