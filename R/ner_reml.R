#' Title
#'
#' @param p dimension of auxiliary matrix X
#' @param XY data matrix containing X,Y and its area number
#' @param x.bar mean auxiliary value of small area
#' @param y.bar mean response value of small area
#' @param K searching number
#'
#' @return resimaxiliklihood estimation of NER model and its blup beta
#' @export
#'
#' @examples 1
NER_REML<-function(p, XY ,x.bar, y.bar, K = 1000){
  X = XY[ ,1:p]
  Y = XY[ ,p+1] 
  ni = as.data.frame(table(XY[,p+2]))[,2]
  m = dim(as.data.frame(table(XY[,p+2])))[1]
  xdot = ni * x.bar
  ydot = ni * y.bar
  n <- sum(ni)
  tXX.inv = solve(t(X) %*% X)
  pmat = diag(n) - X %*% tXX.inv %*% t(X)
  temp1 = t(Y) %*% pmat %*% Y
  temp2 = t(ydot) - t(Y) %*% X %*% tXX.inv  %*% t(xdot)
  ner.gamma<-function(gamma){
    temp3 = solve(diag((ni + 1/gamma), m,m) - xdot %*% tXX.inv %*% t(xdot))
    temp4 = (temp1 - temp2 %*% temp3 %*% t(temp2))/(n - p)
    temp5 = det(diag((1 + ni * gamma), m, m) - gamma * xdot %*% tXX.inv %*% t(xdot))
    return((n-p) * log(temp4) + log(temp5))
  }
  aa<-2*seq(1,K)/K
  bb<-aa
  for(i in 1:K) bb[i]<-ner.gamma(aa[i])
  minbb <- min(bb)
  ghat <- aa[seq(1,K)[bb==minbb]]
  temp6 = solve(diag((ni + 1/ghat), m,m) - xdot %*% tXX.inv %*% t(xdot))
  sehat2 = (temp1 - temp2 %*% temp6 %*% t(temp2))/(n - p) # variance of errors
  svhat = ghat * sehat2  #variance of random effects
  a <- diag(ghat/(1 + ni * ghat))
  temp7 <- t(X) %*% X - t(xdot) %*% a %*% xdot
  temp8 <- t(X) %*% Y - t(xdot) %*%  a %*% ydot
  bhat <- solve(temp7) %*% temp8  #estimators of betahat
  return(list(bhat = c(as.vector(bhat)), sigehat2 = sehat2, sigvhat2 = svhat))
}
