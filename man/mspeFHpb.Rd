\name{mspeFHpb}
\alias{mspeFHpb}
\title{
Compute MSPE through parameter bootstrap method for Fay Herriot model
}
\description{
This function returns MSPE estimator with parameter bootstrap method for Fay Herriot model
}
\usage{
mspeFHpb(Y, X, D, K = 50)
}

\arguments{
  \item{Y}{
a numeric vector. It represents the response value for Fay Herriot model.
}
  \item{X}{
a numeric matrix. Stands for the available auxiliary values.
}
  \item{D}{
a numeric vector. It represents the knowing sampling variance for Fay Herriot model.
}
  \item{K}{
It represents the Monte-Carlo sample number.
}
}
\details{
This method was proposed by Peter Hall and T. Maiti. Parametric bootstrap (pb) method uses bootstrap-based method to measure the accuracy of EB estimator of theta.i.

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

