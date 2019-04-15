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
         `<` = {a != b && all(abs(a)-abs(b) <= 0)},
         `>` = {a != b && all(abs(a)-abs(b) >= 0)},
         `<=` = {a == b || all(abs(a)-abs(b) <= 0)},
         `>=` = {a == b || all(abs(a)-abs(b) >= 0)}
         
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
