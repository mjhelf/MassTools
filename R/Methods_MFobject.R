#' Ops.MFobject
#'
#' Check if an object is a correctly formed MFobject
#' 
#' @param object object to check
#'
#' @export
is.MFobject <- function(object){
  
  is(object, "MFobject") && is.numeric(object) && !is.null(names(object)) && !any(duplicated(names(object)))
}


#' Ops.MFobject
#'
#' Operator behavior for MFobject S3 class
#' 
#' @param e1 numeric or MFobject
#' @param e2 numeric or MFobject
#'
#' @export
Ops.MFobject <- function(e1,e2){
  
  if (nargs() == 1L){
    switch(.Generic,
           `+` = ,
           `-` = NextMethod()
    )
    
  }else{
    
    
    
    
    if(is.MFobject(e1)
       && is.MFobject(e2)){
      switch(.Generic,
             `+` = NextMethod(),
             `-` = NextMethod(),
             `*` = ,
             `/` = ,
             `!=` = any(NextMethod()),
             `==` = all(NextMethod()),
             `<` = {e1 != e2 && all(e1-e2 <= 0)},
             `>` = {e1 != e2 && all(e1-e2 >= 0)},
             `<=` = {e1 == e2 || all(e1-e2 <= 0)},
             `>=` = {e1 == e2 || all(e1-e2 >= 0)}
             
      )
      
      
    }else{
      
      NextMethod() 
      
    }}
}

#' print.MFobject
#'
#' Print method for MFobject S3 class
#' 
#' @param x MFobject
#' @param ... other arguments to \code{base::print()}
#'
#' @export
print.MFobject <- function(x, ...){
  
  print(remakeMF(x), ...)
  
}

#' getSenior3
#'
#'  calculate the sum of valences minus twice (the number of atoms minus 1) (to check for Senior's third theorem, "Golden Rule #2")[1][2][3]
#'
#'
#' @references 
#' \enumerate{
#' \item Kind T, Fiehn O (2007) Seven Golden Rules for heuristic filtering of molecular formulas obtained by accurate mass spectrometry. BMC Bioinformatics 8:105. doi: 10.1186/1471-2105-8-105
#' \item Senior JK (1951) Partitions and Their Representative Graphs. Am J Math 73:663. doi: 10.2307/2372318
#' \item Morikawa T, Newbold BT (2003) Analogous odd-even parities in mathematics and chemistry. Chemistry 12:445–450
#'}
#' @param MF MFobject as constructed by \code{\link{makeMF}}
#'
getSenior3.MFobject <- function(MF){
  
  #if result is >= 0, formula is valid as per Senior's third theorem 
  sum(MF*getOption("MassTools.maxValences")) - 2*(sum(MF)-1)
  
}


#' getExactMass
#' 
#' get the exact mass of a molecule
#' 
#' @param x molecular formula as a character() or MFobject
#' @param ... additional arguments to \code{\link[Rdisop]{getMolecule}}
#' 
#' @export
getExactMass <- function(x, ...){
  
  UseMethod("getExactMass",x)
  
}


#' getExactMass.character
#' 
#' @param x molecular formula as a character()
#' @param ... additional arguments to \code{\link[Rdisop]{getMolecule}}
#' 
#' @importFrom Rdisop getMolecule
#' @export
getExactMass.character <- function(x, ...){
  
  x <- gsub("[[:space:]]","",x)
  
  
   res <- sapply(x,
           function(x){
             tryCatch({
               if(x == ""){return(NA)}
               if(length(grep("-",x))){
                 getExactMass(makeMF(x))
                 }else{
             getMolecule(x, ...)$isotopes[[1]][1]
                   }
             },
             error = function(e){return(NA)})
           })
   names(res) <- x
   
 return(res)
  
}


#' getExactMass.MFobject
#' 
#' @param x molecular formula as an MFobject
#' @param ... additional arguments to \code{\link[Rdisop]{getMolecule}}
#' 
#' @export
getExactMass.MFobject <- function(x, ...){
  
  plus <- if(any(x>0)){getExactMass(remakeMF(x[x>0], ...))}else{0}
  minus <- if(any(x<0)){getExactMass(remakeMF(-x[x<0], ...))}else{0}

  return(plus-minus)
  
}

#' getTopIsotope
#' 
#' get the mass of the most abundant isotope for a molecular formula
#' 
#' @details 
#' Warning: if mixing positive and negative values in one molecular formula,
#'  this may not give the correct result
#'  Warning: Default settings will not work correctly for very large molecules 
#'  (specify \code{maxisotopes} to a value >200 if you expect the largest isotope peak (by intensity)
#'   to be beyond the 200th smallest isotope peaks by mz) 
#' 
#' @param x molecular formula as a character() or MFobject
#' @param ... additional arguments to \code{\link[Rdisop]{getMolecule}}
#' 
#' @importFrom Rdisop getMolecule
#' @export
getTopIsotope <- function(x, ...){
  
  UseMethod("getTopIsotope",x)
  
}


#' @export
getTopIsotope.character <- function(x, ...){
  
  x <- gsub("[[:space:]]","",x)
  
  
  res <- sapply(x,
                function(x){
                  tryCatch({
                    if(x == ""){return(NA)}
                    if(length(grep("-",x))){
                      getTopIsotope(makeMF(x))
                    }else{
                      
                      if(length(list(...)) > 0 && "maxisotopes" %in% names(list(...))){
                      temp <- getMolecule(x, ...)$isotopes[[1]]
                      }else{
                        #200 max isotopes is a compromise between runtime and accuracy for large molecules
                      temp <- getMolecule(x, maxisotopes = 200, ...)$isotopes[[1]]
                        
                        }
                      
                      return(temp[1,which.max(temp[2,])])
                      
                    }
                  },
                  error = function(e){
                    warning(e)
                    return(NA)
                    })
                })
  names(res) <- x
  
  return(res)
  
}



#' @export
getTopIsotope.MFobject <- function(x, ...){
  
  plus <- if(any(x>0)){getTopIsotope(remakeMF(x[x>0], ...))}else{0}
  minus <- if(any(x<0)){getTopIsotope(remakeMF(-x[x<0], ...))}else{0}

  return(plus-minus)
  
}
