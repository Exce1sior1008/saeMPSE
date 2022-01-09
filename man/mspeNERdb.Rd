\name{mspeNERdb}
\alias{mspeNERdb}
\title{
Compute MSPE through double bootstrap(DB) method for Nested error regression model
}
\description{
This function returns MSPE estimator with double bootstrap method for Nested error regression model
}
\usage{
mspeNERdb(ni, X, Y, Xmean, K = 50, C = 50, method)
}
\arguments{
  \item{ni}{
a numeric vector. It represents the sample number for every small area.
}
  \item{X}{
a numeric matrix.  It represents the small area response.
}
  \item{Y}{
a numeric vector. It represents the design matrix.
}
  \item{Xmean}{
a numeric matrix. Stands for the population mean of auxiliary values.
}
  \item{method}{
The MSPE estimation method to be used. See "Details".
}
  \item{K}{
It represents the first bootstrap sample number.
}
  \item{C}{
It represents the second bootstrap sample number.
}
}
\details{
This method was proposed by P. Hall and T. Maiti. Double bootstrap method uses boostrap tool twice for NER model to avoid the unattractivitive bias correction: one is to estimate the estimator bias, the other is to correct for bias.

Default \code{method} is REML method and available variance component estimation method for each method is list as follows:
\code{method = "ML"}
\code{method = "MOM"}

}

\references{
F. B. Butar and P. Lahiri. On measures of uncertainty of empirical bayes small area estimators. Journal of Statistical Planning and Inference, 112(1-2):63-76, 2003.

N. G. N. Prasad and J. N. K. Rao. The estimation of the mean squared error of small-area estimators. Journal of the American Statistical Association, 85(409):163-171, 1990.

Peter Hall and T. Maiti. On parametric bootstrap methods for small area prediction. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 2006a.

H. T. Maiti and T. Maiti. Nonparametric estimation of mean squared prediction error in nested error regression models. Annals of Statistics, 34(4):1733-1750, 2006b.
}
\author{
Peiwen Xiao, Xiaohui Liu, Yuzi Liu, Jiming Jiang, and Shaochu Liu
}


\examples{
1
}

