% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_weight.R
\name{kernel_weight}
\alias{kernel_weight}
\title{Generate the kernel weight matrix}
\usage{
kernel_weight(u, u0 = NULL, ktype = "gaussian", bw = 0.1, poly_order = 0)
}
\arguments{
\item{u}{the condition, it is a vector, one-dimensional array or
one-dimensional row(column) matrix}

\item{u0}{the given condition to be estimated. The default is NULL.
If u0 is NULL, then u0 = u.}

\item{ktype}{the kernel type, can be "gaussian", "epanech", "triweight",
"biweight", "tricube", "triangular" and "cosine",
the default of ktype is "gaussian".}

\item{bw}{bandwidth, default is 0.1}

\item{poly_order}{0 or 1, 0 represents the local constant smoother
while 1 represents the local linear smoother, default is 0.}
}
\value{
the weight matrix
}
\description{
This function generated the kernel weight matrix depending on
different parameter settings.
}
\examples{
\dontrun{

n <- 100
n0 <- 10
u <- seq(-1, 1, length.out = n)
u0 <- seq(-0.9, 0.9, length.out = n0)
# Cosine kernel with bandwidth 0.4, local constant smoother
weight <- kernel_weight(u, u0, ktype = 'cosine', bw = 0.4, poly_order = 0)
}
}
\seealso{
\code{\link{computeUdiff}}, \code{\link{kernelCompute}}
}
