\name{mspe_FH_Sucma}
\alias{mspe_FH_Sucma}
\title{
Sucma estimator of mspe
}
\description{
estimate mspe of FH model using Sucma method.
}
\usage{mspe_FH_Sucma()}
\arguments{
\item{p}{dimension of X.}

\item{m}{small area number.}

\item{X}{a numeric vector. It represents the response or the observed value in the Fay Herriot Model}

\item{Y}{a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.}

\item{D}{a numeric vector consisting of the known sampling variances of each of the small area levels.}

\item{b}{estimate of b.}

\item{Ahat}{estimate of A.}

\item{K}{Monte-Carlo sample number.}

}
\value{
mspe
}
\references{
Jiang, J., and Torabi, M.(2020). Sumca: simple, unified, Monte-Carlo-assisted approach to
second-order unbiased mean-squared prediction error estimation. Journal of the Royal Statistical
Society: Series B (Statistical Methodology), 82(2), 467-485.}
\examples{
1
## maybe str(mspe_FH_DB) ; plot(mspe_FH_DB) ...
}