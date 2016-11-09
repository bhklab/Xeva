#if(1==2){


FirstQuadrant <- setClass(
  # Set the name for the class
  "FirstQuadrant",

  # Define the slots
  slots = c(
    x = "numeric",
    y = "numeric",
    z = "list"
  ),

  # Set the default values for the slots. (optional)
  prototype=list(
    x = 0.0,
    y = 0.0
  )#,

  # Make a function that can test to see if the data is consistent.
  # This is not called if you have an initialize function defined!
  #validity=function(object)
  #{
  #  if((object@x < 0) || (object@y < 0)) {
  #    return("A negative number for one of the coordinates was given.")
  #  }
  #  return(TRUE)
  #}
)


x <- FirstQuadrant()

y <- FirstQuadrant(x=5,y=7, z=list(a=1:5))

#' @export
foox <- function(){print( "3" )}

#




































##
#}