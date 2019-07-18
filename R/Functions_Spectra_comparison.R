#' pairCompare
#'
#' Compare two vectors of equal length pairwise and calculate a similarity 
#' score (e.g. intensity values of matched MS peaks).
#'
#' 
#' @param v1 vector 1
#' @param v2 vector 2
#' @param NAasZero replace NA values with 0 if the value is NA in one vector, but not the other
#' @param method "cosine" will use lsa::cosine(), "pearson" will use stats::cor()
#' 
#' @importFrom stats cor
#' @importFrom lsa cosine
#'
#' @return the result of the method call, a numeric similarity score
#'
pairCompare <- function(v1, v2, NAasZero = T,
                        method = c("cosine", "pearson",
                                   "kendall", "spearman")){
  
  #remove instances where values in both vectors are NA
  remfeats <- is.na(v1) & is.na(v2)
  
  if(any(remfeats)){
  v1 <- v1[!remfeats]
  v2 <- v2[!remfeats]
  }
  
  if(NAasZero){
    
    # if only one vector is NA, set that value to 0.
    v1[is.na(v1)] <- 0
    v2[is.na(v2)] <- 0
    
    # selfeats <- which(!is.na(v1) | !is.na(v2))
    # 
    # v1[selfeats[which(is.na(v1[selfeats]))]] <- 0
    # v2[selfeats[which(is.na(v2[selfeats]))]] <- 0
    # 
    
    
  }
  
  if(method[1] == "cosine"){
    cosine(matrix(c(v1,v2),ncol = 2, byrow = F))[1,2]
  }else{
    cor(v1,v2, method = method[1], use = "pairwise.complete.obs")
  }
}

#' network1
#'
#' Compare two spectra and calculate a similarity score 
#' (e.g. intensity values of matched MS peaks)
#' NOTE: Current limitations: (1) peaks from one spectrum can match multiple 
#' peaks from the other scan (workaround: call quickMergeMS with ppm = 0 and 
#' mzdiff same as mzdiff here on both spectra prior to calling network1)
#'
#' 
#' @param spec1 spectrum 1
#' @param spec2 spectrum 2
#' @param mztol allowed m/z difference between peaks to match matching
#' @param parentshift m/z difference between parentmasses (importantly, 
#' has to be Parent(spec2) - Parent(spec1) ) to find alternative matches /
#'  neutral loss matches; will only be calculated if \code{abs(parentshift) >  abs(mztol)}
#' @param method method passed on to \code{\link{pairCompare}()}; 
#' "cosine" will use lsa::cosine(), "pearson" will use stats::cor()
#' @param minpeaks minimum number 
#' of peaks that have to be matched, otherwise returns 0
#' @param nonmatched if TRUE, will add non-matching peaks to calculation, 
#' with 0 intensity in the spectrum missing the peak
#' 
#' @return the result of the method call, a numeric similarity score
#'
#' @export
network1 <- function(spec1, spec2, mztol = 0.005,
                     parentshift = 0,
                     method = "cosine",
                     minpeaks = 6,
                     nonmatched = T){
  
  posx <- which(abs(outer(spec1[,1], spec2[,1], "-")) < abs(mztol), arr.ind = T)
  
  if(abs(parentshift) > abs(mztol)){
    
    posx <- rbind(posx, which(abs(outer(spec1[,1]+parentshift, spec2[,1], "-")) < abs(mztol), arr.ind = T))
    
    posx <- posx[!(duplicated(posx[,1]) & duplicated(posx[,2])),  ]
    
    if(length(posx) ==2){
      posx <- matrix(posx,nrow =1)
    }
    
  }
  
  if(nrow(posx) >= minpeaks && nrow(posx) > 0){
    
    if(nonmatched){
      
      #v1 <- numeric(length(spec1[,1]) + length(spec2[,2]) - length(unique(posx )
      
      v1 <- c(spec1[posx[,1],2], spec1[-posx[,1],2], rep(0, length(spec2[-posx[,2],2])))
      v2 <- c(spec2[posx[,2],2], rep(0, length(spec1[-posx[,1],2])), spec2[-posx[,2],2] )
      return(pairCompare(v1, v2, method = method))
      
    }
    
    return(pairCompare(spec1[posx[,1],2], spec2[posx[,2],2], method = method))
    
  }else{
    return(0)
  }
  
  
}

#' makeEdges
#'
#' Make an edgelist (data.frame) using a vector of parentmasses 
#' and a list of MS spectra. Wrapper for \code{\link{network1}()}.
#'
#' 
#' @param speclist (non-nested) list of MS spectra
#' @param parentmasses vector of parent m/z values (same length as speclist).
#' @param mztol max difference between matched peaks in m/z
#' @param method "cosine" will use lsa::cosine(), "pearson" will use stats::cor()
#' @param minpeaks minimum number of peaks that have to be matched, otherwise returns 0
#' @param nonmatched if TRUE, will add non-matching peaks to calculation, 
#' with 0 intensity in the spectrum missing the peak
#' 
#' @return a data.frame with edge information, see \code{Details}
#' 
#' @details Columns of the returned data.frame:
#' \itemize{
#' \item \code{from} vertex id
#' \item \code{to} vertex id
#' \item \code{cosine} similarity measure
#' \item \code{deltamz} parent mass difference
#' }
#'
#' @export
makeEdges <- function(speclist,
                      parentmasses = NULL,
                      mztol = 0.005, 
                      method = "cosine",
                      minpeaks = 6, 
                      nonmatched = T){
  
  #remove NULLs from the speclist, but keep track of the indices of the non-NULL scans:
  selNonNulls <- which(!sapply(speclist,function(x){return(is.null(x) || nrow(x) == 0)}))
  speclist <- speclist[selNonNulls]
  
  selectlist <- list()
  for(i in seq(length(speclist)-1)){
    selectlist[[i]] <- i:length(speclist)
  }
  
  if(!is.null(parentmasses)){
    parentmasses <- parentmasses[selNonNulls]
    
    #inefficient but prob wont take long
    pmassShifts <- list()
    for(n in seq(length(speclist)-1)){
      #spec2 - spec1!!
      pmassShifts[[n]] <- parentmasses[selectlist[[n]][1]] - parentmasses[selectlist[[n]][-1]]
    }
    
    
    alledges <-  mapply(function(sel, specs, pmasses){
      mapply(network1, spec1 = specs[sel[-1]],  parentshift = pmasses,
             MoreArgs = list(spec2 = specs[[sel[1]]],
                             mztol = mztol,
                             minpeaks = minpeaks,
                             method = method,
                             nonmatched = nonmatched))
    }, 
    sel = selectlist, 
    pmasses = pmassShifts,
    MoreArgs = list(specs = speclist),
    SIMPLIFY = F)
  }else{
    
    
    alledges <-  lapply(selectlist, function(sel, specs, mzt, mp, nonm){
      lapply(specs[sel[-1]],network1,
             spec2 = specs[[sel[1]]],
             method = method,
             mztol = mzt, 
             minpeaks = mp, 
             nonmatched = nonm)
    }, specs = speclist, mzt = mztol, mp = minpeaks, nonm = nonmatched)
    
  }
  
  
  
  
  
  sz <- sum(sapply(alledges, length))
  
  dt <- data.frame(from = integer(sz),
                   to = integer(sz),
                   cosine = numeric(sz),
                   stringsAsFactors = F)
  
  seq(6,1)
  dt$from = unlist(mapply(rep, 1:length(alledges), sapply(alledges, length)))
  dt$to = unlist(sapply(selectlist, "[", -1))
  
  if(!is.null(parentmasses)){
    
    dt$deltamz = unlist(pmassShifts)
    
  }
  
  dt$cosine = unlist(alledges)
  
  #put the node indices back together correctly
  dt$from <- selNonNulls[dt$from]
  dt$to <- selNonNulls[dt$to]
  
  return(dt)
  
}