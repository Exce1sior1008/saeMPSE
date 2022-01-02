#' Title
#'
#' @param ni sampled units observed in the ith small area
#' @param X auxiliary values matrix
#' @param Y response variable
#' @param X.mean mean auxiliary value of population
#' @param x.bar mean auxiliary value of small area
#' @param y.bar mean response value of small area
#' @param method estimation method of sig_e2 and sig_v2
#'
#' @return mspe of NER model using DL method
#' @export
#'
#' @examples 1
mspe_NER_DL = function(ni, X, Y, X.mean, x.bar, y.bar,method = "REML"){
  n = sum(ni)
  m = length(ni)
  p = ncol(X)

  # if(method == "ML"){
  #   phat = NER_ML(p, ni, XY, xbar, ybarMat)
  #   sig_e2 = phat$sigehat2
  #   sig_v2 = phat$sigvhat2
  # }
  if(method == "REML"){
    phat = NER_REML(p, cbind(X,Y,rep(1:m,ni)), x.bar, y.bar)
    sig_e2 = c(phat$sigehat2)
    sig_v2 = c(phat$sigvhat2)
  }
  Vi.inv = list()

  for (i in 1:m) {
    vtemp = 1/sig_e2 * diag(ni[i]) - sig_v2/((ni[i]*sig_v2+sig_e2)*sig_e2)*rep(1,ni[i])%*%t(rep(1,ni[i]))
    Vi.inv[[i]]<-vtemp
  }
  V.inv = as.matrix(bdiag(Vi.inv))
  gama = ni * sig_v2 / (ni*sig_v2 + sig_e2)
  g1temp = (1 - gama) * sig_v2

  w = sig_e2 + ni * sig_v2
  a = sum(ni^2 * w^(-2)) * sum((ni -1)* sig_e2^(-2) + w^(-2))-(sum(ni^2 * w^(-2)))^2
  Ivv = 2 * a^(-1) * sum((ni-1) * sig_e2^(-2) + w^(-2))
  Iee = 2 * a^(-1) * sum(ni^2 * w^(-2))
  Ive = -2 * a^(-1) * sum(ni * w^(-2))
  g1 = c();g2 = c();g3 = c();

  for (i in 1:m) {
    g1[i] = g1temp[i]
    g2[i] = t(X.mean[i, ] - gama[i] * x.bar[i, ])%*%solve(t(X) %*% V.inv %*% X)%*%(X.mean[i, ] - gama[i] * x.bar[i, ])
    g3[i] = (1/(ni[i]^2 * (sig_v2 + sig_e2 /ni[i])^3)) * (sig_e2^2 * Ivv+sig_v2^2 * Iee-2 * sig_e2 * sig_v2 * Ive)
  }
  return(g1+g2+2*g3)
}


