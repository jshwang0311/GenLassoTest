#' This function gets the glasso solution for a fixed lambda.
#'
#' @param object ESPgenlasso object
#' @param lambda
#' @return the primal solution
#' @export
#'
get.beta <- function(objet, lambda, tol = 1e-10)
{
  object.index <- which(object$lambda == lambda)
  if(length(object.index) == 0){
    object.index <- which(object$lambda < lambda)[1] - 1
    if(object.index == 0){
      object.beta <- object$beta[, 1]
      object.beta[which(abs(object.beta) <= tol)] <- 0
    } else{
      object.beta <- -(object$beta[, object.index]-object$beta[, (object.index+1)])/
        (object$lambda[object.index] - object$lambda[object.index+1]) * (object$lambda[object.index]-lambda) + object$beta[, object.index]
      object.beta[which(abs(object.beta) <= tol)] <- 0
    }
  } else{
    object.beta <- object$beta[, object.index[1]]
    object.beta[which(abs(object.beta) <= tol)] <- 0
  }
  return(object.beta)
}
