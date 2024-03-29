% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShrinkageCoefLC.R
\name{ShrinkageCoefLC}
\alias{ShrinkageCoefLC}
\title{ShrinkageCoefLC}
\usage{
ShrinkageCoefLC(
  YY,
  Cov_vec,
  u,
  u_test,
  h,
  ysquare,
  Sigmau,
  Cov_LL,
  ktype = "gaussian"
)
}
\arguments{
\item{YY}{the data, num_non_zero * n matrix}

\item{Cov_vec}{the weighted average of YY}

\item{u}{the condition}

\item{u_test}{the test center}

\item{h}{the bandwidth}

\item{ysquare}{the square of standadization of y}

\item{Sigmau}{the diagonal entries estimator at u_test}

\item{Cov_LL}{the covariance matrix}

\item{ktype}{the kernel type, can be "gaussian", "epanech", "triweight",
"biweight", "tricube", "triangular" and "cosine",
the default of ktype is "gaussian".}
}
\value{
the shrinkage covariance
}
\description{
This function computes the shrinkage covariance after threshold .
}
