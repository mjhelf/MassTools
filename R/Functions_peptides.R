#' PeptideMF
#'
#' convert peptide sequences into an MFobjects
#' 
#' @param seq one or more peptide sequences in a character vector 
#' @param AAs molecular formulas for each amino acid letter, defaults to proteinogenic amino acids
#' 
#' @return a list of MFobjects
#' 
#' @export
PeptideMF <- function(seq, AAs = c("A" = "C3H7NO2",  "R" = "C6H14N4O2",  "N" = "C4H8N2O3",
                                   "D" = "C4H7NO4",  "C" = "C3H7NO2S",  "Q" = "C5H10N2O3",
                                   "E" = "C5H9NO4",  "G" = "C2H5NO2",  "H" = "C6H9N3O2",
                                   "I" = "C6H13NO2",  "L" = "C6H13NO2",  "K" = "C6H14N2O2",
                                   "M" = "C5H11NO2S",  "F" = "C9H11NO2",  "P" = "C5H9NO2",
                                   "S" = "C3H7NO3",  "T" = "C4H9NO3",  "W" = "C11H12N2O2", 
                                   "Y" = "C9H11NO3",  "V" = "C5H11NO2")){
  seq <- gsub("[[:space:]]","",seq)
  
  pepAAs <- strsplit(seq, split = "")
  Ecounts <- lapply(pepAAs, function(p){makeMF(paste0(paste(AAs[p], collapse = ""),paste(rep("H-2O-1",length(p)-1), collapse = "") ))})
  
  return(Ecounts)
}

#' countAAs
#'
#' count amino acids in a sequence
#' 
#' @param seq one or more peptide sequences in a character vector 
#' 
#' @return a list of named numeric vectors
#' 
#' @export
countAAs <- function(seq){
  
  seq <- gsub("[[:space:]]","",seq)
  
  pepAAs <- strsplit(seq, split = "")
  
  uniques <- lapply(pepAAs, unique)
  
  res <- mapply(function(u,p){counts <- sapply(u,function(aa){sum(p==aa)}); names(counts) <- u; return(counts)},uniques,pepAAs, SIMPLIFY = F)
  
  return(res)
}
