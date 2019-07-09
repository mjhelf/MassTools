#' mergeMS
#'
#' Merge MS spectra, best used with noise-free spectra. Will merge all peaks that are EITHER within ppm OR within mzdiff range (mzdiff typically used to allow larger ppm difference in low mz range)
#'
#' Replacement for mergeMS(..., mergeOnly = T)
#' 
#' @param speclist list of spectra (matrices, all with two columns: \code{mz} and \code{intensity})
#' @param ppm min difference between peaks in ppm
#' @param mzdiff min difference between peaks in m/z
#' @param count if TRUE, will count how many peaks were combined into a peak, and add a \code{count} column in the resulting spectrum
#' @param iterative if TRUE, will iteratively merge two adjacent peaks within tolerance, 
#' and then check if there are peaks within tolerance of merged peak
#' @param removeZeros remove all entries with 0 intensity 
#' @param noiselevel all peaks below this intensity (relative to 
#' highest peak in merged spectrum at relative intensity 1) will be removed 
#' @param maxpeaks if not NULL, maximum number of peaks in merged spectrum.
#' Lowest intensity peaks will be removed from merged spectrum.
#'
#' @export
mergeMS <- function(speclist, ppm =5, mzdiff = 0.0005,
                    count = FALSE, iterative = T,
                    removeZeros = T,
                    noiselevel = 0,
                    maxpeaks = NULL){
  
  if(length(speclist) == 0){return(matrix(numeric(0),ncol = 2, dimnames = list(NULL,c("mz", "intensity"))))}
  
  if(is.list(speclist) || is(speclist,"Spectra")){
    
    if(is(speclist[[1]],"Spectrum")){
      SpectrumOut <- speclist[[1]]
      
      aspec <- do.call(rbind,lapply(speclist, function(x){matrix(c(x@mz, x@intensity),ncol = 2)}))
      

      }else{
       SpectrumOut <- NULL 
       
           aspec <- do.call(rbind,speclist)

      }
      

    
  }else{
    SpectrumOut <- NULL 
    
    aspec <- speclist
    
  }
  
  #safeguards
  if(is.null(aspec)){return(matrix(numeric(0),ncol = 2, dimnames = list(NULL,c("mz", "intensity"))))}
  
  
  #remove 0 intensity peaks because they will produce NaN mz values downstream
  
  if(removeZeros){
    aspec <- aspec[aspec[,2] != 0,, drop = FALSE]
  }
  
  #if(length(aspec) == 0){return(NULL)}
  
  
  if(nrow(aspec) > 1){
    aspec <- aspec[order(aspec[,1]),]
  }else{
   return(aspec) 
  }
  
    if(count){
    aspec <- cbind(aspec, matrix(rep(1,nrow(aspec)),ncol = 1, dimnames = list(NULL,"count")))
    }
    
    
  margins <- diff(aspec[,1])
  
  margins_ppm <- margins/aspec[-nrow(aspec),1]/(1e-6)
  
  belowmargin <- (margins <= mzdiff | margins_ppm <= ppm)
  
  
  #note: conveniently evaluates to FALSE if only one peak left (and belowmargin is logical(0) as result of )
  while(any(belowmargin)){
  
    #index of first peak in a group that is belowmargin, peaks following a first peak belowmargin should not count
    #should dorectly translate to index in aspec..
  firsts <- which( !c(FALSE,belowmargin[-length(belowmargin)]) & belowmargin)
  
  if(iterative){
  seconds <- firsts + 1
  }else{
    #finding the last entry in the group, iterative aggregation will merge all
    #into one without the "mz walking" iterative effect
  seconds <- which(belowmargin & !c(belowmargin[-1], FALSE))
    
  }
  
 
  
  sumints <- (aspec[firsts,2] + aspec[seconds,2]) #intensity sum
  
  
   if(!removeZeros && any(sumints == 0)){
     rem <- which(sumints == 0)
     firsts <- firsts[-rem]
     seconds <- seconds[-rem]
     }
  
  aspec[firsts,1] <-   (aspec[firsts,1]*aspec[firsts,2] + aspec[seconds,1]*aspec[seconds,2]) / sumints #mz weighted average
  
  aspec[firsts,2] <- sumints   
  
  if(count){
    aspec[firsts,3] <-  aspec[firsts,3] + 1
  }
  
  #remove the seconds             
  aspec <- aspec[-seconds,, drop = FALSE]  
  
  
  #reassess loop condition:
  margins <- diff(aspec[,1])
  
  margins_ppm <- margins/aspec[-nrow(aspec),1]/(1e-6)
  
  belowmargin <- (margins <= mzdiff | margins_ppm <= ppm)
  
  
  }
  
  #averaging is done at the end so that weighted averages are calculated with sum instead of mean intensities
  #mean intensities would not behave desirably e.g. if a large peak is followed by one or more very small peaks - 
  #the weight of the large peak gets "diluted" and the next large peak that it may be combined with is overrepresented
  if(is.list(speclist)){
  aspec[,2] <- aspec[,2]/length(speclist)
  }
  
  if(noiselevel > 0){
   
    aspec <- aspec[aspec[,2] >= max(aspec[,2])*noiselevel,, drop = F]
     
  }
  
  if(!is.null(maxpeaks) && maxpeaks < nrow(aspec)){
    if(maxpeaks < 1){
     aspec <- aspec[FALSE,,drop = F] 
    }else{
    sel <- seq(nrow(aspec)) %in% na.omit(order(aspec[,2], decreasing = T)[seq(maxpeaks)])
    aspec <- aspec[sel,, drop = F]
    }
  }
  
  if(exists("SpectrumOut") && !is.null(SpectrumOut)){
    
    SpectrumOut@mz <- aspec[,1]
    SpectrumOut@intensity <- aspec[,2]
    
    SpectrumOut@peaksCount <- nrow(aspec)
  
    return(SpectrumOut)
  }
    
    
  return(aspec)
  
}


#' annotateSpectrum
#'
#' Match peaks in a spectrum with those in a labeling data.frame 
#' and prepare them for plotting functions. Matches have to be within abs OR ppm 
#'
#' @param labels a data.frame as returned by permutatePeptideMass
#' @param spectrum a mass spectrum as data.frame or matrix
#' @param mzlabel if TRUE, will add mz values to the label
#' @param unmodifiedTag what tag to use for unmodified fragment ions
#' @param ppm relative matching tolerance in ppm
#' @param abs matching tolerance in absolute numbers 
#'
#' @return plots an annotated peptide sequence in the current plotting device
#'
#'@export
annotateSpectrum <- function(labels, spectrum,
                             mzlabel = F,
                             unmodifiedTag = "",
                             ppm = 5,
                             abs = 0
                             ){
  
  
  foundpeaks <- matchNumbers(labels$mz, spectrum[,1], ppm = ppm, abs = abs)
  
  foundpeaks_annotated <- cbind(as.data.frame(spectrum)[foundpeaks[,2],],
                                labels[foundpeaks[,1],])
  
  colnames(foundpeaks_annotated)[1:2] <- c("mzSpec","intensitySpec")
  
  foundpeaks_annotated$color <- "blue3"
  
  foundpeaks_annotated$label <- gsub("\\[\\]"," ",paste0(foundpeaks_annotated$type,
                                       "[",foundpeaks_annotated$pos,
                                       "]","^{+{}}*", ' "[',
                                       gsub("^$",unmodifiedTag,
                                            foundpeaks_annotated$modifications),
                                       ']"',
                                       if(mzlabel){paste0("(",format(round(foundpeaks_annotated$mzSpec,4),nsmall = 4, scientific = F),")")}else{""}
  ))
  
  return(foundpeaks_annotated)
}