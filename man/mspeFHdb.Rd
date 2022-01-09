\name{mspeFHdb}
\alias{mspeFHdb}
\title{
Compute MSPE through double bootstrap method for Fay Herriot model
}
\description{
This function returns MSPE estimator with double bootstrap appoximation method for Fay Herriot model
}
\usage{
mspeFHdb(Y, X, D, K = 50, C = 50, method)
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
It represents the first bootstrap sample number.
}
  \item{C}{
It represents the second bootstrap sample number.
}
  \item{method}{
It represents the variance component estimation method.
}
}
\details{
This method was proposed by P. Hall and T. Maiti. Double bootstrap method uses boostrap tool twice for Fay Herriot model to avoid the unattractivitive bias correction: one is to estimate the estimator bias, the other is to correct for bias.

Default value for \code{method} is 2, \code{method = 2} represents the REML method and available variance component estimation method for each method is list as follows:
\code{method = 1} represents MOM method,

}

\references{
P. Hall and T. Maiti. On parametric bootstrap methods for small area prediction. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 2006.
}
\author{
Peiwen Xiao, Xiaohui Liu, Yuzi Liu, Jiming Jiang, and Shaochu Liu
}


\examples{
1
}

