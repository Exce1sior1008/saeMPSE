#' Title
#'
#' @param beta coefficient
#' @param A unknown sampling variance
#' @param D known sampling variance
#' @param m sampling number
#' @param p dimention of sample
#'
#' @return simulation data for Fay Herriot model
#' @export
#'
#' @examples 1
sim_data_FH = function(beta, A, D, m, p){
  X = matrix(runif(m * p), m, p)
  X[,1] = rep(1, m)
  theta = X %*% beta + rnorm(m, 0, sqrt(A))
  Y = theta + rnorm(m, 0, sqrt(D))
  return(list(X = X, Y = Y, D = D, theta = theta))
}
