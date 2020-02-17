% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/YAverage.R
\name{YAverage}
\alias{YAverage}
\title{YAverage}
\usage{
YAverage(Y, u, h)
}
\arguments{
\item{Y}{the observation, p * n matrix}

\item{u}{the condtion, 1 * n matrix}

\item{h}{the bandwidth, scalar}
}
\value{
This function will return list object which contains 3 variables:

\item{aver}{each column is the weighted average of Y given the condition u(i) and h. It's a p*n size matrix.}
\item{weight}{the n * n weight matrix, the ith row represents the weight for the u(i), it is also a sparse matrix.}
\item{alpha}{the coefficient adjusted after the i-th observation deleted.}
}
\description{
This function is used for estimating the local constant
kernel weighted average at each center u(i) given bandwidth h.
Here we use gaussian kerenel.
}
\examples{
\dontrun{

data(Y)
data(u)
data(LowerBoundary)
# h is fixed
h <- 0.1
result <- YAverage(Y, u, h)
}
}
\seealso{
\code{\link{CVLL}}
}