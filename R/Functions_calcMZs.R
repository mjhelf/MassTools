#' calcIons
#'
#' Calculate charge-to-mass ratios for charge states of ions
#' 
#' @param MF character vector with molecular formulas (MFs)
#' @param charges vector of positive or negative integers indicating the charge states to calculate
#' @param carrier charge carrier - Molecular formula will be added to or removed from MFs for each charge
#' @param monoisotopic if true, will calculate the monoisotopic mass m/z values, 
#' will calculate m/z of the most abundant isotope peak otherwise (preferrable e.g. for large organic molecules)
#' @param mf_column name of column in result data.frame that will contain the molecular formulas
#' @param adduct an adduct, e.g. "Na" or "Cl". If the charge is positive, each adduct replaces one charge carrier.
#' 
#' @return a data.frame with columns \code{mz}, \code{charge}, \code{ion} and the \code{mf_column} as specified
#'
#' @export
calcIons <- function(MF, charges = c(1), carrier = "H",
                    monoisotopic = T,
                    mf_column = "formula", adduct = NULL){
  
  
  if(charges[1] >= 0 ){
    
    inpforms <- paste0(MF, paste(rep(carrier, abs(charges[1])-length(adduct)), sep = "", collapse = ""), adduct)
    
    
  }else{
    
    inpforms <- paste0(MF, paste(rep(paste0(carrier,"-1"), abs(charges[1])), sep = "", collapse = ""), adduct)
    
    
  }
  

  if(monoisotopic){
    
    mzs <- (getExactMass(inpforms) - 5.48579909070e-4*charges[1])/max(c(abs(charges[1]),1))
    
  }else{
    mzs <- (getTopIsotope(inpforms) - 5.48579909070e-4*charges[1])/max(c(abs(charges[1]),1))
     }
  
  
  
  newdf <- data.frame(mz = mzs,
                      stringsAsFactors = F)
  
  newdf[[mf_column]] <- MF
  newdf$charge <- charges[1]
  
  hnumber <- if(charges[1] > 0){charges[1] - length(adduct)}else{charges[1] + length(adduct)}
  newdf$ion <- paste0("[M",
                      if(length(adduct)>0){paste0("+", adduct)}else{""},
                      if(hnumber > 0){paste0("+",if(hnumber!=1){hnumber}else{""}, carrier,"]")}else if (hnumber < 0){paste0(if(hnumber!=-1){hnumber}else{"-"}, carrier,"]")}else{"]"},
                      if(charges[1] > 0 ){paste0(charges[1],"+")}else if(charges[1] < 0){paste0(abs(charges[1]),"-")}else{""}  )
  
  #newdf$ion <- paste0("[M",  if(charges[1] > 0 ){"+"}else{"-"}, if(abs(charges[1]) > 1){abs(charges[1])}else{""}, if(abs(charges[1]) > 0){paste0(carrier, "]", abs(charges[1]), )}else{"]"})
  
  
  if(length(charges) > 1){
    
    return(do.call(rbind, list(a= newdf, b = calcIons(MF, charges = charges[-1],
                                                     carrier, monoisotopic, mf_column)   )))
    
  }
  return(newdf)
  
}