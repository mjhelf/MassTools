#' @title MassTools
#' @name MassTools
#' @import utils methods graphics stats
#'
NULL


.onLoad <- function(libname,pkgname){
  
  #if(!"MassTools.elements" %in% names(options())){
  
  
  elementNames <- c("C","H","N","O","P","S", 
                               "F", "Cl", "Br", "I", "Si",
                               "B", "As", "Se", 
                               "Li", "Na", "K", "Mg", "Ca",
                               "Fe", "Co", "Mn", "Cu", "Zn")
  
  maxValences <-c("C" = 4,
                  "H" = 1,
                  "N" = 4,
                  "O" = 2,
                  "P" = 5,
                  "S" = 6, 
                  "F" = 1, "Cl" = 1, "Br" = 1, "I" = 1, "Si" = 4,
                  "B" = 3, "As" = 5, "Se" = 6,
                  "Li" = 1, "Na" = 1, "K" = 1, "Mg" = 2, "Ca" = 2,
                  "Fe" = 3, "Co" = 4, "Mn" = 4, "Cu" = 2, "Zn" =2)
  
  op.MassTools <- list(
    MassTools.elements = integer(length(elementNames)),
    MassTools.maxValences = maxValences
  )
  
  names(op.MassTools$MassTools.elements) <- elementNames
  
  
  options(op.MassTools[!names(op.MassTools) %in% names(options())])
  
  
  
}