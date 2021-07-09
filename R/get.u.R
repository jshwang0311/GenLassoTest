#' This function gets the dual solution of the genlasso problem for a fixed lambda.
#'
#' @param object ESPgenlasso object
#' @param lambda
#' @return the dual solution
#' @export
#'
get.u <- function(objet, lambda)
{
  object.index <- which(object$lambda == lambda)
  if(length(object.index) == 0){
    object.index <- which(object$lambda < lambda)[1] - 1
    if(object.index == 0){
      object.u <- object$u[, 1]
    } else{
      object.u <- -(object$u[, object.index]-object$u[, (object.index+1)])/
        (object$lambda[object.index] - object$lambda[object.index+1]) * (object$lambda[object.index]-lambda) + object$u[, object.index]
    }
  } else{
    object.u <- object$u[, object.index[1]]
  }
  return(object.u)
}
