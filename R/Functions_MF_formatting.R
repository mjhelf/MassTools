#''reformatSF
#'
#' Reformats a molecula formula string so that it has a count for each element, even if the count is 1.
#' 
#' @param s Molecula formula string
#'
#' @export
reformatSF <- function(s){
  
  foundLL <- regexpr("[A-z][A-Z]|[A-z]$",s)
  
  if(foundLL == -1){
    return(s)
  }
  else{
    
    return( reformatSF(paste0(substr(s,1,foundLL),"1", substr(s,foundLL+1,nchar(s)))))
  }
  
}


#''consolidateSF
#'
#' Merges duplicates in an element count vector
#' 
#' @param x element count vector
#'
#' @export
consolidateSF <- function(x){
  
  dups <- duplicated(names(x))
  
  if(!any(dups)){
    class(counts) <- "MFobject"
    return(x)
  }
  
  eles <- unique(names(x))
  counts <- integer(length(eles))
  names(counts) <- eles
  
  for(i in names(counts)){
    
    counts[i] <- sum(x[names(x) == i])
    
  }
  
  class(counts) <- "MFobject"
  
  return(counts)
  
}


#''makeSF
#'
#' Reformat a molecula formula string into an element count vector
#' 
#' @param s Molecula formula string
#'
#' @export
makeSF <- function(s){
  
  s <- reformatSF(s)
  elnums <- strsplit(s, "[A-Z]|[A-Z][a-z]")[[1]]
  elnums <- as.integer(elnums[-1])
  names(elnums) <- strsplit(s, "-[0-9]+|[0-9]+")[[1]]
  return(consolidateSF(elnums))  
}

#''remakeSF
#'
#' Converts an element count vector back into a molecula formula string
#' 
#' @param s Molecula formula string
#'
#' @export
remakeSF <- function(x){
  
  return(paste(names(x[x!=0]), x[x!=0], collapse = "", sep = ""))
  
}