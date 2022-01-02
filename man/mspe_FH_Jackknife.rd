\name{mspe_FH_Jackknife}
\alias{mspe_FH_Jackknife}
\title{
Jackknife estimator of mspe
}
\description{
estimate mspe of FH model using jackknife method.
}
\usage{mspe_FH_Jackknife(p,m,X,Y,D,bhat,A_REML,maxiter)}
\arguments{
\item{p}{dimension of X.}

\item{m}{small area number.}

\item{X}{a numeric vector. It represents the response or the observed value in the Fay Herriot Model}

\item{Y}{a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.}

\item{D}{a numeric vector consisting of the known sampling variances of each of the small area levels.}

\item{bhat}{estimate of b.}

\item{A_REML}{estimate of A.}

\item{maxiter}{max iteration during REML estimation.}

}
\value{
mspe
}
\references{
Jiang, J., Lahiri, P., and Wan, S. M.(2002). A unified jackknife theory for empirical best prediction
with M-estimation. The Annals of Statistics, 30(6), 1782-1810.}
\examples{
1
## maybe str(mspe_FH_DB) ; plot(mspe_FH_DB) ...
}

