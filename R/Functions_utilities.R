#' matchNumbers
#'
#' match numbers between two numeric vectors with tolerances. 
#' A match is defined by complying wit either the ppm OR the 
#' abs tolerance values.
#' 
#' @return an indexing matrix of matching values, 
#' with matching pairs in each row, column 1 is the v1 match,
#'  row 2 is the v2 match.
#' 
#' @param v1,v2 numeric vectors
#' @param ppm matching tolerance in ppm
#' @param abs matching tolerance in absolute numbers 
#'
#' @export
matchNumbers <- function(v1, v2, ppm = 5, abs = 0.001){
  
  
  mx <- abs(outer(v1, v2, "-"))
  
  return(which(mx <= abs | mx <= v1*ppm*1e-6, arr.ind = T))
  
}


#' permutateMass
#' 
#' Use a data.frame with molecular formulas and combine them
#' 
#' @param baseMass numeric mass value to permutate
#' @param modifications data.frame with columns name, formula, min, max. 
#' @param sep separator between individual modifications in output;
#'  has to be escaped with \\ if it is a special character for regex
#' @param unmodifiedTag which tag to use for non-modified items. If NULL,
#'  will preceed modification tag with 0 (e.g. 0Me)
#'
#' @details More information
#' \subsection{ modifications data.frame}{
#' Columns that need to be in the modifications data.frame:
#' \itemize{
#' \item{min} minimum number of times this modification can be applied
#' \item{max} maximum number of times this modification can be applied
#' \item{mass} mass of this modification
#' \item{tag} label for this modification (will automatically be preceded by 
#' number of these modifications applied)
#' }
#' columns min and max should be the minimum and maximum times a modification 
#' is allowed in the combinations, can be positive, negative, or 0.
#' }
#'  
#' @return a data.fram, the mass in the resulting data.frame is the  monoisotopic mass.
#' 
#' @seealso \code{\link{permutatePeptideMass}}
#' 
#' @export
permutateMass <- function(baseMass, modifications,
                          sep = "\\|", unmodifiedTag = ""){
  
  if(!length(baseMass) 
     || (is.data.frame(baseMass) && !nrow(baseMass))){
    return( data.frame(mass = numeric(),
                       modifications = character()))}
  
  if(!is.data.frame(baseMass)){baseMass <- data.frame(mass = baseMass,
                                                      modifications = "")}
  
  if((missing(modifications) 
     || is.null(modifications) 
     ||!nrow(modifications))
     && !"modifications" %in% colnames(baseMass)){
    baseMass$modifications <- ""
    }
  
  puresep <- gsub("\\\\","",sep)
  
  #base case for recursion
  if(missing(modifications) 
     || is.null(modifications) 
     || nrow(modifications) < 1){
    
    #remove separators that follow each other 
    baseMass$modifications <- gsub(paste0(sep,"{2,}"),puresep,
                                   baseMass$modifications) 
    
    #remove separators at beginning and end
    baseMass$modifications <- gsub(paste0("^", sep,"*|", sep,"$"),
                                   "",baseMass$modifications)
    
    
    
    return(baseMass[!duplicated(baseMass),])
    
  }
  
  
  res <- lapply(seq(nrow(baseMass)),function(i){
    
    thismod <- paste0(seq(modifications[1,"min"],
                          modifications[1,"max"]),
                      modifications[1,"tag"])
    
    
    if(!is.null(unmodifiedTag)
       && any(seq(modifications[1,"min"],modifications[1,"max"]) == 0)){
      thismod[seq(modifications[1,"min"],
                  modifications[1,"max"]) == 0] <- unmodifiedTag
    }
    
    
    out <- data.frame(mass = seq(modifications[1,"min"],
                                 modifications[1,"max"])*modifications[1,"mass"] + baseMass[i,"mass"],
                      modifications = paste(baseMass[i,"modifications"],
                                            thismod, sep = puresep),
                      stringsAsFactors = F)
    
    
    return(out)
    
  })
  
  
  #now have to combine the modifications
  
  res <- do.call(rbind,res)
  
  return(permutateMass(res,modifications[-1,], sep = sep))
}