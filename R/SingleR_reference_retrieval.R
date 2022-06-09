# Sandra Lilja
# R version 4.0.4
#'
#' Extract the reference into the right format for cell type identification
#' 
#' @param out_file full path to output file
#' 
#' @export
#'

SingleR_reference_retrieval <- function(out_file){
  library(SingleR)
  library(dplyr)
  mpca.se <- MouseRNAseqData()
  ReferenceDataMouse <- as.data.frame(t(mpca.se@assays@data@listData$logcounts))
  ReferenceDataMouse$'IDs' <- rownames(ReferenceDataMouse) 
  
  lbls <- data.frame(IDs = mpca.se@colData@rownames,
                     MainCellType = mpca.se$label.main,
                     FineCellType = mpca.se$label.fine)
  
  ReferenceDataMouse <- full_join(lbls, ReferenceDataMouse)
  ReferenceDataMouse$IDs <- NULL
  
  write.csv(ReferenceDataMouse, out_file, row.names = F)
}
