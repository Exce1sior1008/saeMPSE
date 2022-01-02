#' Title
#'
#' @param response a numeric vector. It represents the response or the observed value in the Fay Herriot Model
#' @param designmatrix a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.
#' @param sampling.var a numeric vector consisting of the known sampling variances of each of the small area levels.
#'
#' @return 1
#' @export
#'
#' @examples 1
mspe_PR = function( response, designmatrix , sampling.var ){
  A = prasadraoest( response, designmatrix , sampling.var )$estimate
  m = length( response )
  p = ncol( designmatrix )
  g1iA = c()
  g2iA = c()
  g3iA = c()
  for(i in 1:m){
    g1iA[i] = A * sampling.var[i] / ( A+sampling.var[i] )
    g2iA[i] = ( ( sampling.var[i] / ( A + sampling.var[i] ) ) ^ 2) * t( designmatrix [i,]) %*% solve(t( designmatrix ) %*% diag((A + sampling.var)^(-1)) %*%  designmatrix  ) %*%  designmatrix [i,]
    varA    = 2 * m^(-2)*(A^2 + 2 * A * sum(sampling.var) / m + sum( sampling.var^2 ) / m )
    g3iA[i] = ( sampling.var[i]^2 ) / ( ( A + sampling.var[i] )^3 ) * varA
  }
  mspe = g1iA + g2iA + g3iA
  return(mspe)
}
