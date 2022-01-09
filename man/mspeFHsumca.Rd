\name{mspeFHsumca}
\alias{mspeFHsumca}
\title{
Compute MSPE through Sumca method for Fay Herriot model
}
\description{
This function returns MSPE estimator with the combination of linearization and resampling appoximation method called "Sumca", for Fay Herriot model
}
\usage{
mspeFHsumca(Y, X, D, K = 50, method = 2)
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
a numeric vector. It represents the Monte-Carlo sample number for "sumca".
}
  \item{method}{
It represents the variance component estimation method.
}
}

\details{
This method was proposed by J. Jiang, P. Lahiri, and T. Nguyen, sumca method combines the advantages of linearization and resampling methods and obtains unified, positive, low-computation burden and second-order unbiased MSPE estimators.

Default value for \code{method} is 2, \code{method = 2} represents the REML method and available variance component estimation method for each method is list as follows:
\code{method = 1} represents MOM method.

}

\references{
J. Jiang and M. Torabi. Sumca: simple; unified; monte carlo assisted approach to second order unbiased mean squared prediction error estimation. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 82(2):467-485, 2020.
}

\author{
Peiwen Xiao, Xiaohui Liu, Yuzi Liu, Jiming Jiang, and Shaochu Liu
}


\examples{
1
}

