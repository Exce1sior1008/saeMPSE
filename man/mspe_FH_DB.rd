\name{mspe_FH_DB}
\alias{mspe_FH_DB}
\title{
Double bootstrap estimator of mspe
}
\description{
estimate mspe of FH model using double bootstrap sample method.
}
\usage{mspe_FH_DB()}
\arguments{
\item{X}{a numeric vector. It represents the response or the observed value in the Fay Herriot Model}

\item{Y}{a numeric matrix. The first column is a column of ones(also called the intercept). The other columns consist of observations of each of the covariates or the explanatory variable in Fay Herriot Model.}

\item{D}{a numeric vector consisting of the known sampling variances of each of the small area levels.}

\item{K}{first step bootstrap sample number.}

\item{C}{second step bootstrap sample number.}

\item{b}{estimate of b.}

\item{Ahat}{estimate of A.}

\item{maxiter}{max iteration during REML estimation.}

}
\value{
mspe
}
\references{
Hall, P., and Maiti, T.(2006a). On parametric bootstrap methods for small area prediction.
Hall, P., and Maiti, T.(2006b). Nonparametric estimation of mean-squared prediction error in
nested-error regression models. The Annals of Statistics, 34(4), 1733-1750.
}
\examples{
1
## maybe str(mspe_FH_DB) ; plot(mspe_FH_DB) ...
}

