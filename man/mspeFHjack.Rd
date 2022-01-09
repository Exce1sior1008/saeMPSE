\name{mspeFHjack}
\alias{mspeFHjack}
\title{
Compute MSPE through resampling method for Fay Herriot model
}
\description{
This function returns MSPE estimator with jackknife method for Fay Herriot model
}
\usage{
mspeFHjack(Y, X, D, method = 2)
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
  \item{method}{
It represents the variance component estimation method.
}

}

\details{
This method was proposed by J. Jiang and L. S. M. Wan, jackknife method is used to obtain the bias and variation of estimators.

Default value for \code{method} is 2, \code{method = 2} represents the REML method and available variance component estimation method for each method is list as follows:
\code{method = 1} represents MOM method.

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

