\name{mspeNERsumca}
\alias{mspeNERsumca}
\title{
Compute MSPE through Sumca method for Nested error regression model
}
\description{
This function returns MSPE estimator with the combination of linearization and resampling appoximation method for Nested error regression model
}
\usage{
mspeNERsumca(ni, X, Y, Xmean, K = 50, method)
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
  \item{K}{
a numeric vector. It represents the Monte-Carlo sample number for "sumca".
}
  \item{method}{
The MSPE estimation method to be used. See "Details".
}
}

\details{
This method was proposed by J. Jiang, P. Lahiri, and T. Nguyen, sumca method combines the advantages of linearization and resampling methods and obtains unified, positive, low-computation burden and second-order unbiased MSPE estimators.

Default \code{var.method} is "REML" and the available variance component estimation method for each method is list as follows:
\code{var.method = "ML"}
\code{var.method = "MOM"}

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

