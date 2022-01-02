#' Title
#'
#' @param response a numeric vector. It represents the response or the observed value in the Fay Herriot Model
#' @param designmatrix a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.
#' @param sampling.var a numeric vector consisting of the known sampling variances of each of the small area levels.
#' @param method Ahat estimator method
#'
#' @return 1
#' @export
#'
#' @examples 1
mspe_Jackknife = function( response , designmatrix , sampling.var , method = "REML"){
  m <- length( response )
  p <- ncol( designmatrix )
  bhatFun = function(A, response , designmatrix , sampling.var ){
    V.inv = diag(1/(A + sampling.var ))
    bhat = solve(t( designmatrix ) %*% V.inv %*% designmatrix ) %*% t( designmatrix ) %*% V.inv %*% response
    return(bhat)
  }
  if(method == "REML"){
    A = resimaxilikelihood( response, designmatrix , sampling.var , 100)$estimate
  }
  if(method == "ML"){
    A = maximlikelihood( response, designmatrix , sampling.var )$estimate
  }
  if(method == "MOM"){
    A = prasadraoest( response, designmatrix , sampling.var )$estimate
  }
  bHat = bhatFun( A , response , designmatrix , sampling.var )
  mseU = function(sampling.var, index, designmatrix , response, m){
    bi.u = thetaHat.u = numeric(m)
    for(u in 1:m){
      dataX = designmatrix[-u, ]
      dataY = response[-u, ]

      if(method == "REML"){
        A.u = resimaxilikelihood( response, designmatrix , sampling.var , 100)$estimate
      }
      if(method == "ML"){
        A.u = maximlikelihood( response, designmatrix , sampling.var )$estimate
      }
      if(method == "MOM"){
        A.u = prasadraoest( response, designmatrix , sampling.var )$estimate
      }

      bHat.u = bhatFun( A.u, dataY , dataX ,  sampling.var[-u] )
      Bi.u = sampling.var[index]/(A.u +  sampling.var[index])
      bi.u[u] = A.u * Bi.u
      thetaHat.u[u] = designmatrix[index, ] %*% bHat.u + (1 - Bi.u)*(response[index] - designmatrix[index, ] %*% bHat.u)
    }
    return(list(bi.u, thetaHat.u))
  }

  mspe = numeric(m)
  Bi =  sampling.var/( sampling.var + A)
  bi = A * Bi
  thetaHat.i = designmatrix %*% bHat + (1 - Bi)*(response - designmatrix %*% bHat)
  for(i in 1:m){
    uHat = mseU(sampling.var, i, designmatrix, response ,m)
    bihatU = uHat[[1]]
    thetahatU = uHat[[2]]
    mspe[i] = bi[i] - ((m-1)/m) * sum(bihatU - bi[i]) + ((m-1)/m) * sum((thetahatU - thetaHat.i[i])^2)
  }
  return(mspe)
}



