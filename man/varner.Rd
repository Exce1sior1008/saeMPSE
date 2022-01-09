\name{varner}
\alias{varner}
\title{
Estimate of the variance component using several methods for Nested error regression model.
}
\description{
This function returns the estimate of variance component with several existing method for Nested error regression model. This function does not accept missing values.
}
\usage{
varner(ni, X, Y, method)
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
  \item{method}{
The variance component estimation method to be used. See "Details".
}
}
\details{
Default \code{method} for \code{varner} is MOM, also called ANOVA estimator, set \code{method = 1}  the available variance component estimaion method are list as follows: \code{method = 1} represents the moment estimator; \code{method = 2} represents the restricted maximum likelihood(REML) estimator; \code{method = 3} represents the maximum likelihood(ML) estimator;  
}

\references{
J. Jiang. Linear and Generalized Linear Mixed Models and Their Applications. Linear and Generalized Linear Mixed Models and Their Applications, 2007.
}
\author{
Peiwen Xiao, Xiaohui Liu, Yuzi Liu, Jiming Jiang, and Shaochu Liu
}


\examples{
1
}

