\name{mspeNERlnr}
\alias{mspeNERlnr}
\title{
Compute MSPE through linearization method for Nested error regression model
}
\description{
This function returns MSPE estimator with linearization appoximation method for Nested error regression model
}
\usage{
mspeNERlnr(ni, X, Y, X.mean, method = "PR", var.method = "default")
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
  \item{X.mean}{
a numeric matrix. Stands for the population mean of auxiliary values.
}
  \item{method}{
The MSPE estimation method to be used. See "Details".
}
  \item{var.method}{
The variance component estimation method to be used. See "Details".
}
}
\details{
Default \code{method} for \code{mspeNERlnr} is "PR" ,proposed by N. G. N. Prasad and J. N. K. Rao, Prasad-Rao (PR) method uses Taylor series expansion to obtain a second-order approximation to
the MSPE. Function \code{mspeNERlnr} also provide the following method:

Method "DL" advanced PR method to cover the cases when the variance components are estimated by ML and REML estimator. Set \code{method = "DL"}.
}

\references{
N. G. N. Prasad and J. N. K. Rao. The estimation of the mean squared error of small-area estimators. Journal of the American Statistical Association, 85(409):163-171, 1990.

G. S. Datta and P. Lahiri. A unified measure of uncertainty of estimated best linear unbiased predictors in small area estimation problems. Statistica Sinica, 10(2):613-627, 2000.
}
\author{
Peiwen Xiao, Xiaohui Liu, Yuzi Liu, Jiming Jiang, and Shaochu Liu
}


\examples{
1
}

