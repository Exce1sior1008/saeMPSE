#' Title
#'
#' @param response a numeric vector. It represents the response or the observed value in the Fay Herriot Model
#' @param designmatrix a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.
#' @param sampling.var a numeric vector consisting of the known sampling variances of each of the small area levels.
#'
#' @return estimate
#' @export
#'
#' @examples 1
EmpiricalBayes = function( response, designmatrix , sampling.var ){
  p = ncol( designmatrix )
  m = length(response)
  invtXX =  solve(t(designmatrix) %*% designmatrix)
  beta_ols = invtXX %*% t(designmatrix) %*% response
  hii = numeric(m)
  for(i in 1:m){
    hii[i] = t(designmatrix[i,]) %*% invtXX %*% designmatrix[i,]
  }
  A_tilde = (sum((response - designmatrix %*% beta_ols)^2) - sum((1 - hii) * sampling.var))/(m-p)
  Ahat = max(0,A_tilde)
  return(Ahat)
}
