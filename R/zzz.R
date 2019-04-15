.onLoad <- function(libname,pkgname){
  
  if(!"MassTools.elements" %in% names(options())){
  
  
  elementNames <- c("C","H","N","O","P","S", 
                               "F", "Cl", "Br", "I", "Si", "B", "As",
                               "Li", "Na", "K", "Mg", "Ca",
                               "Fe", "Co", "Mn", "Cu", "Zn")
  op.MassTools <- list(
    MassTools.elements = integer(length(elementNames))
  )
  
  names(op.MassTools$MassTools.elements) <- elementNames
  
  options(op.MassTools)
  }
  
  
}