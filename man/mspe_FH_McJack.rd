\name{mspe_FH_McJack}
\alias{mspe_FH_McJack}
\title{
Mcjack estimator of mspe
}
\description{
estimate mspe of FH model using Monte-Carlo jackknife method.
}
\usage{mspe_FH_McJack()}
\arguments{
\item{m}{small area number.}

\item{X}{a numeric vector. It represents the response or the observed value in the Fay Herriot Model}

\item{Y}{a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.}

\item{D}{a numeric vector consisting of the known sampling variances of each of the small area levels.}

\item{b}{estimate of b.}

\item{Ahat}{estimate of A.}

\item{K}{first step bootstrap sample number.}

\item{maxiter}{max iteration during REML estimation.}

}
\value{
mspe
}
\references{
Jiang, J., Lahiri, P., and Nguyen, T.(2016). A unifed Monte-Carlo jackknife for small area
estimation after model selection. arXiv preprint arXiv:1602.05238.}
\examples{
1
## maybe str(mspe_FH_DB) ; plot(mspe_FH_DB) ...
}

