\name{mspeNERjack}
\alias{mspeNERjack}
\title{
Compute MSPE through resampling method for Nested error regression model
}
\description{
This function returns MSPE estimator with jackknife appoximation method for Nested error regression model
}
\usage{
mspeNERjack(ni, X, Y, Xmean, method)
}

\arguments{
  \item{ni}{
a numeric vector. It represents the sample number for every small area.
}
  \item{X}{
a numeric matrix. Stands for the available auxiliary values.
}
  \item{Y}{
a numeric vector. It represents the response value for Nested error regression model.
}
  \item{Xmean}{
a numeric matrix. Stands for the population mean of auxiliary values.
}
  \item{method}{
The MSPE estimation method to be used. See "Details".
}

}

\details{
This method was proposed by J. Jiang and L. S. M. Wan, jackknife method is used to obtain the bias and variation of estimators.

Default \code{method} is REML method and available variance component estimation method for each method is list as follows:
\code{method = "ML"}
\code{method = "MOM"}
}

\references{
M. H. Quenouille. Approximate tests of correlation in time series. Journal of the Royal Statistical Society. Series B (Methodological), 11(1):68-84, 1949.

J. W. Tukey. Bias and confidence in not quite large samples. Annals of Mathematical Statistics, 29(2):614, 1958.

J. Jiang and L. S. M. Wan. A unified jackknife theory for empirical best prediction with m estimation. Annals of Statistics, 30(6):1782-1810, 2002.

J. Jiang, P. Lahiri, and T. Nguyen. A unified monte carlo jackknife for small area estimation after model selection, 2016.
}
\author{
Peiwen Xiao, Xiaohui Liu, Yuzi Liu, Jiming Jiang, and Shaochu Liu
}


\examples{
1
}

