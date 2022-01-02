#' Title
#'
#' @param ni sampled units observed in the ith small area
#' @param X auxiliary values matrix
#' @param Y response variable
#' @param X.mean mean auxiliary value of population
#' @param x.bar mean auxiliary value of small area
#' @param y.bar mean response value of small area
#'
#' @return mspe of NER model using PR method
#' @export
#'
#' @examples 1
mspe_NER_PR = function(ni, X, Y, X.mean, x.bar, y.bar){
  n = sum(ni)
  m = length(ni)
  p = ncol(X)
  y.star = numeric(n)
  x.star = matrix(NA, n, p)
  xmat = 0
  for(j in 1:m){
    xmat = xmat + ni[j]^2 * x.bar[j, ] %*% t(x.bar[j, ])
    y.star[(sum(ni[1:j]) - ni[j] + 1) : sum(ni[1:j])] = Y[(sum(ni[1:j]) - ni[j] + 1) : sum(ni[1:j])] - y.bar[j]
    x.star[(sum(ni[1:j]) - ni[j] + 1) : sum(ni[1:j]), ] = X[((sum(ni[1:j]) - ni[j] + 1) : sum(ni[1:j])), ] -
      matrix(rep(x.bar[j, ], ni[j]), ni[j], 2, byrow = TRUE)
  }
  estar = residuals(lm(y.star ~ x.star[,2:p]))
  u     = residuals(lm(Y ~ X))
  ntemp = solve(t(X) %*% X)
  nstar = n - sum(diag(ntemp %*% xmat))
  sige2 = (1/(n - m - 2 + 1)) * sum(estar^2)
  sigv2 = max((1/nstar) * (sum(u^2) - (n - 2) * sige2),0)
  gama = ni * sigv2/(ni * sigv2 + sige2)
  Vi.inv = list()
  Zi = list()
  g1 = c()
  g2 = c()
  g3 = c()

  for(t in 1:m){
    vtemp = 1/sige2 * (diag(ni[t]) - (gama[t]/ni[t]) * rep(1,ni[t]) %*% t(rep(1,ni[t])))
    Vi.inv[[t]] = vtemp
    Zi[[t]] = rep(1,ni[t])
  }
  V.inv = as.matrix(bdiag(Vi.inv))
  Zmat = as.matrix(bdiag(Zi))
  M = diag(n) - X %*% ntemp %*% t(X)
  nstar2 = (sum(diag(M %*% Zmat %*% t(Zmat))))^2
  vare = (2/(n - m - 2 + 1)) * sige2^2
  varv = (2/nstar^2) * ((1/(n - m - 2 + 1)) * (m - 1) * (n - p) * sige2^2 +
                          2 * nstar * sige2 * sigv2 + nstar2 * sigv2^2)
  covev = -(m - 1) * (1/nstar) * vare
  msePR = numeric(m)
  for(i in 1:m){
    g1[i] = (1 - gama[i]) * sigv2
    g2[i] = t(X.mean[i, ] - gama[i] * x.bar[i, ])%*%solve(t(X) %*% V.inv %*% X)%*%(X.mean[i, ] - gama[i] * x.bar[i, ])
    g3[i] = (1/(ni[i]^2 * (sigv2 + sige2 /ni[i])^3)) *
            (sige2^2 * varv + sigv2^2 *
            vare - 2 * sige2 * sigv2 * covev)
  }
  msePR = g1 + g2 + 2 * g3
  return(msePR)
}




