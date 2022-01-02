\name{mspe_FH_PB}
\alias{mspe_FH_PB}
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

\item{K}{bootstrap sample number.}

}
\value{
mspe
}
\references{
Butar, F. B., and Lahiri, P.(2003). On measures of uncertainty of empirical Bayes small-area
estimators. Journal of Statistical Planning and Inference, 112, 63-76.
}
\examples{
1
## maybe str(mspe_FH_DB) ; plot(mspe_FH_DB) ...
}

