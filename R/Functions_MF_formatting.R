#''reformatMFstring
#'
#' Reformats a molecula formula string so that it has a count for each element, even if the count is 1.
#' 
#' @param s Molecula formula string
#'
#' @export
reformatMFstring <- function(s){
  
  s <- gsub(" ", "",s)
  
  #finds adjacent element names
  foundLL <- regexpr("[A-z][A-Z]|[A-z]$",s)
  
  sel <- foundLL != -1
  
  while(any(sel)){
   
    s[sel] <- paste0(substr(s[sel],1,foundLL[sel]),"1", substr(s[sel],foundLL[sel]+1,nchar(s[sel])))
    
    foundLL[sel] <- regexpr("[A-z][A-Z]|[A-z]$",s[sel])
    
    
    sel <- foundLL != -1
    
  }
  
  return(s)
  
}


#' consolidateMF
#'
#' Merges duplicates in an element count vector
#' 
#' @param x element count vector
#'
#' @export
consolidateMF <- function(x){
  
  target <- getOption("MassTools.elements")
  
 # x <- unclass(x)
  
  matched <- match(names(x),names(target))
  
  if(!any(duplicated(names(x)))){
    
    target[matched] <- x
    
  }else{
  
  for(i in unique(matched)){
    
    target[i] <- sum(x[matched == i])
    
  }
  }
    
  class(target) <- "MFobject"
  
  return(target)
  
}


#' makeMF
#'
#' construct an MFobject from a molecula formula string
#' 
#' @param s Molecula formula string
#' @param forcelist if \code{TRUE}, will always return a list, even when s is of length 1.
#'
#' @export
makeMF <- function(s, forcelist = F){
  
  s <- reformatMFstring(s)
  elnums <- strsplit(s, "[A-Z]|[A-Z][a-z]")
  elnums <- lapply(elnums,function(el){as.integer(el[-1])})
  elnames <- strsplit(s, "-[0-9]+|[0-9]+")
  
  elnums <- mapply(function(enu,ena){names(enu) <- ena
  return(enu)}, enu = elnums, ena = elnames, SIMPLIFY = F)
  
  consolidated <- lapply(elnums,consolidateMF)
  
  if(!forcelist && length(consolidated) ==1){
  return(consolidated[[1]])  
  }else{
    return(consolidated)  
    
  }
    
}

#' remakeMF
#'
#' Converts an element count vector back into a molecula formula string
#' 
#' @param x Molecula formula string
#'
#' @export
remakeMF <- function(x){

  return(paste(names(x[x!=0]), gsub("^1$","",x[x!=0]), collapse = "", sep = ""))
  
}

