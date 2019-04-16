#' Ops.MFobject
#'
#' Check if an object is a correctly formed MFobject
#' 
#' @param a object
#'
#' @export
is.MFobject <- function(a){
  
 class(a) == "MFobject" && is.numeric(a) && !is.null(names(a)) && !any(duplicated(names(a)))
}


#' Ops.MFobject
#'
#' Operator behavior for MFobject S3 class
#' 
#' @param a numeric or MFobject
#' @param b numeric or MFobject
#'
#' @export
Ops.MFobject <- function(a,b){
  
  if (nargs() == 1L){
    switch(.Generic,
           `+` = ,
           `-` = NextMethod()
    )
    
  }else{
  
    
    
    
  if(class(a) == "MFobject"
    && class(b) == "MFobject"){
 switch(.Generic,
         `+` = NextMethod(),
         `-` = NextMethod(),
         `*` = ,
         `/` = ,
         `!=` = any(NextMethod()),
         `==` = all(NextMethod()),
         `<` = {a != b && all(a-b <= 0)},
         `>` = {a != b && all(a-b >= 0)},
         `<=` = {a == b || all(a-b <= 0)},
         `>=` = {a == b || all(a-b >= 0)}
         
  )
 
 
  }else if(class(a) != class(b)){
  
   NextMethod() 
    
  }}
}

#' print.MFobject
#'
#' Print method for MFobject S3 class
#' 
#' @param s MFobject
#'
#' @export
print.MFobject <- function(s){

   print(remakeMF(s))
   
}
