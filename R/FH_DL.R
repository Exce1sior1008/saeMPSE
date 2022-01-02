#' Title
#'
#' @param response a numeric vector. It represents the response or the observed value in the Fay Herriot Model
#' @param designmatrix a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.
#' @param sampling.var a numeric vector consisting of the known sampling variances of each of the small area levels.
#' @param method Ahat estimator method
#'
#' @return mspe
#' @export
#'
#' @examples 1
mspe_DL=function( response, designmatrix , sampling.var , method = "REML"){
  if(method == "REML"){
    A = resimaxilikelihood( response, designmatrix , sampling.var , 100)$estimate
  }
  if(method == "ML"){
    A = maximlikelihood( response, designmatrix , sampling.var )$estimate
  }
  m = length( response )
  p = ncol( designmatrix )
  g1iA = c()
  g2iA = c()
  g3iA = c()
  for(i in 1:m){
    g1iA[i] = A * sampling.var[i]/( A + sampling.var[i] )
    temp = matrix(0,p,p)
    for (u in 1:m) {
      temp = temp +  designmatrix [u,] %*% t( designmatrix [u,])/(A + sampling.var[u])
    }
    g2iA[i] = (( sampling.var[i] / (A + sampling.var[i] ))^2) * t( designmatrix [i,] ) %*% solve( temp ) %*%designmatrix [i,]
    g3iA[i] = ( sampling.var[i]^2)/((A + sampling.var[i])^3)/(sum((A + sampling.var)^-2))
  }
  mspe = g1iA + g2iA + 2 * g3iA
  return(mspe)
}


