#' Roll up data to a specified taxonomic level
#'
#' This function rolls up data to the specified taxonomic level so that statistics may be calculated at various taxonomic levels.
#'
#' @param omicsData an object of the class 'seqData' created by \code{\link{as.seqData}}.
#' @param level taxonomic level to roll up to.
#' @param taxa_levels The levels of taxonomy (or other e_meta object) which might be used in the roll up. If NULL, will use c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), in that order. Default is NULL.
#'
#' @details Data will be rolled (summed) up to a specified taxonomic level. For example, data at the OTU level could be rolled (summed) up to the Genus level before statistics are computed.
#'
#' @return A seqData object of the same class as the input, where e_data and e_meta are rolled up to the specified level.
#'
#' @examples
#' \dontrun{
#' library(mintJansson)
#' data(rRNA_data)
#'
#' rRNA_split <- split_emeta(rRNA_data)
#'
#' rRNA_rollup <- taxa_rollup(omicsData = rRNA_split, level = "Phylum")
#'
#' dim(rRNA_rollup$e_data)
#' attributes(rRNA_rollup)
#' }
#'
#' @author Allison Thompson
#'
#' @export
taxa_rollup <- function(omicsData, level, taxa_levels=NULL){

  # Must have emeta to map to
  if(is.null(omicsData$e_meta)){
    stop("e_meta must not be NULL")
  }

  # Level to rollup to must be in emeta
  if(!(level %in% colnames(omicsData$e_meta))){
    stop("level for rollup must be in e_meta")
  }

  # Check that this is count data
  if(!(class(omicsData) %in% c("seqData"))){
    stop("omicsData is not an object of appropriate class")
  }

  if(is.null(taxa_levels)){
    taxa_levels <- c("kingdom","phylum","class","order","family","genus","species")
  }else{
    if(!is.character(taxa_levels)){
      stop("taxa_levels must be a character vector")
    }
  }
  numlevel <- which(taxa_levels %in% tolower(level))

  # Make new taxonomy column
  #idx <- which(colnames(omicsData$e_meta) == attr(omicsData, "cnames")$edata_cname)
  e_meta <- omicsData$e_meta
  e_meta <- e_meta[,c(1:which(colnames(e_meta)==level))]
  e_meta <- apply(e_meta, 1:2, function(x) as.character(x))

  newTaxa_data <- e_meta[,which(tolower(colnames(e_meta)) %in% taxa_levels)]
  newTaxa <- data.frame(ID=omicsData$e_meta[,attr(omicsData,"cnames")$edata_cname],
                        newTaxa=unlist(lapply(c(1:nrow(newTaxa_data)), function(y) paste(newTaxa_data[y,1:numlevel], collapse=";"))))

  e_meta <- as.data.frame(e_meta)
  e_meta <- merge(e_meta, newTaxa, by.x=attr(omicsData,"cnames")$edata_cname, by.y="ID")

  # Roll up data to appropriate level
  rollup <- lapply(1:length(levels(e_meta$newTaxa)), function(x){
    features <- e_meta[which(e_meta$newTaxa == levels(e_meta$newTaxa)[x]),attr(omicsData, "cnames")$edata_cname]
    sub <- omicsData$e_data[which(omicsData$e_data[,attr(omicsData,"cnames")$edata_cname] %in% features),]

    rollup <- t(as.data.frame(colSums(sub[,-which(colnames(sub)==attr(omicsData,"cnames")$edata_cname)], na.rm=TRUE)))
    rollup <- data.frame(Level=paste(level,x,sep=""), rollup)
    rownames(rollup) <- NULL

    return(rollup)
  })


  # Format new edata
  e_data <- do.call(rbind, rollup)
  names(e_data)[1] <- paste(level,"Id",sep="")
  e_data[e_data==0] <- NA

  # Format new emeta
  e_meta <- e_meta[!duplicated(e_meta$newTaxa), ]
  e_meta <- e_meta[match(levels(e_meta$newTaxa),e_meta$newTaxa),]
  # Rename ID column
  e_meta[,attr(omicsData, "cnames")$edata_cname] <- paste(level,c(1:nrow(e_meta)), sep="")
  names(e_meta)[which(names(e_meta)==attr(omicsData,"cnames")$edata_cname)] <- paste(level,"Id",sep="")
  rownames(e_meta) <- NULL
  e_meta <- e_meta[,-which(colnames(e_meta)=="newTaxa")]

  # Run function to make pmartRseq object again
  if(class(omicsData) == "seqData"){
    data <- as.seqData(e_data=e_data, f_data=omicsData$f_data, e_meta=e_meta, edata_cname=paste(level,"Id",sep=""),
                        taxa_cname=level, ec_cname=attr(omicsData,"cnames")$ec_cname, gene_cname=attr(omicsData,"cnames")$gene_cname,
                        fdata_cname=attr(omicsData,"cnames")$fdata_cname, data_type=attr(omicsData,"data_info")$data_type,
                        data_scale=attr(omicsData,"data_info")$data_scale, data_types=attr(omicsData,"data_info")$data_types,
                        db=attr(omicsData,"database")$db, db_version=attr(omicsData,"database")$db_version,
                        data_norm=attr(omicsData,"data_info")$data_norm, norm_method=attr(omicsData,"data_info")$norm_method,
                        location_param=attr(omicsData,"data_info")$location_param, scale_param=attr(omicsData,"data_info")$scale_param)
    attr(data,"group_DF") <- attr(omicsData, "group_DF")
  }

  attr(data, "rollup")$text = paste("Data was rolled up to the ", level, " level", sep="")
  attr(data, "rollup")$level = level

  return(data)

}
