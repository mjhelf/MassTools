#' mergeMS
#'
#' Merge MS spectra, best used with noise-free spectra. Will merge all peaks that are EITHER within ppm OR within mzdiff range (mzdiff typically used to allow larger ppm difference in low mz range)
#'
#' Replacement for mergeMS(..., mergeOnly = T)
#' 
#' @param speclist list of spectra (matrix with mz and intensity values)
#' @param ppm min difference between peaks in ppm
#' @param mzdiff min difference between peaks in m/z
#'
#' @export
mergeMS <- function(speclist, ppm =5, mzdiff = 0.0005){
  
  if(length(speclist) == 0){return(NULL)}
  
  if(is.list(speclist)){
    aspec <- do.call(rbind,speclist)
  }else{
    aspec <- speclist
    
  }
  
  #lots of safeguards
  if(is.null(aspec)){return(NULL)}
  
  
  #remove 0 intensity peaks because they will produce NaN mz values downstream
  
    aspec <- aspec[aspec[,2] != 0,, drop = F]
  
  
  if(length(aspec) == 0){return(NULL)}
  
  
  if(nrow(aspec) != 1){
    aspec <- aspec[order(aspec[,1]),]
  }else{
   return(aspec) 
  }
  
  
  margins <- diff(aspec[,1])
  
  margins_ppm <- margins/aspec[-nrow(aspec),1]/(1e-6)
  
  belowmargin <- (margins <= mzdiff | margins_ppm <= ppm)
  
  
  #note: conveniently evaluates to FALSE if only one peak left (and belowmargin is logical(0) as result of )
  while(any(belowmargin)){
  
    #index of first peak in a group that is belowmargin, peaks following a first peak belowmargin should not count
    #should dorectly translate to index in aspec..
  firsts <- which( !c(F,belowmargin[-length(belowmargin)]) & belowmargin)
  seconds <- firsts + 1
  
  sumints <- (aspec[firsts,2] + aspec[seconds,2]) #intensity mean
  
  aspec[firsts,1] <-   (aspec[firsts,1]*aspec[firsts,2] + aspec[seconds,1]*aspec[seconds,2]) / sumints #mz weighted average
  
  aspec[firsts,2] <- sumints   
  #remove the seconds             
  aspec <- aspec[-seconds,, drop = F]  
  
  
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
  
  return(aspec)
  
}