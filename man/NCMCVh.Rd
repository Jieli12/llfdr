% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NCMCVh.R
\name{NCMCVh}
\alias{NCMCVh}
\title{NCMCVh}
\usage{
NCMCVh(Y, u, h, ktype = "gaussian")
}
\arguments{
\item{Y}{the observation, p * n matrix}

\item{u}{the condition, it is a vector, one-dimensional array or
one-dimensional row(column) matrix}

\item{h}{the bandwidth, scalar}

\item{ktype}{the kernel type, can be "gaussian", "epanech", "triweight",
"biweight", "tricube", "triangular" and "cosine",
the default of ktype is "gaussian".}
}
\value{
the value of cross validation function
}
\description{
This routine calculate the cross validation value using NCM method.
}
\examples{
\dontrun{

data(u)
data(LowerBoundary)
data(Ystd_LC)
upper <- 1
h_grid <- matrix(seq(0.05, upper, length.out = 100), nrow = 100)
cv <- apply(h_grid, 1, NCMCVh, Y = Ystd_LC, u = u)
plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
# select the optimal bandwidth for diagonal entries of covariance
hncm <- optimise(NCMCVh, c(LowerBoundary, upper), tol = 1e-6,
                   Y = Ystd_LC, u = u)
abline(v = hncm$minimum, col="red")
}
}
\seealso{
\code{\link{CVLL}}
}
