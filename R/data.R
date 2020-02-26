#' @title The data
#'
#' @description  This R package is translated from our initial MATLAB code.
#' In order to keep functions performimg consistently. We load only one small dataset saved by
#' MATLAB. Here, Y is a data matrix with 50 rows and 100 columns.
#'
#' @usage data(Y)
#'
#' @format  A data \code{matrix} with 50 rows and 100 columns.
"Y"


#' @title The condition
#' @description A condition vector with 100 length.
#' @usage data(u)
"u"


#' @title The Covariance Matrix
#' @description The underlying covariance matrix.
#' @usage data(Sigma)
"Sigma"

#' @title The Noise Matrix
#' @description The noise matrix.
#' @usage data(Noise)
"Noise"

#' @title Lower Boundary
#' @description The lower boundary of bandwidth
#' @usage data(LowerBoundary)
"LowerBoundary"

#' @title Index of subset-y
#' @description The subset of y variables, size: 4 * 50. For more details, see
#'  Chen, Z. and Leng, C., 2016. Dynamic covariance models.
#'  Journal of the American Statistical Association, 111(515), pp.1196-1207.
#' @usage data(Index.y)
"Index.y"


#' @title Index of observations
#' @description The subset of observations, size: 78 * 100, is used for thresholding.
#'  For more details, see Bickel, P. J. and Levina, E. (2008).
#'  Covariance regularization by thresholding. The Annals of Statistics, 36(6):2577–2604.
#' @usage data(Index.obs)
"Index.obs"

#' @title Bandwidth
#' @description The bandwidth was computed by \code{\link{CVLC}}
#' @usage data(h1_lc)
"h1_lc"

#' @title Bandwidth
#' @description The bandwidth was computed by \code{\link{CVLL}}
#' @usage data(h1_ll)
"h1_ll"

#' @title The residual
#' @description The residual was computed by \code{\link{YResidualH1LC}}
#' @usage data(Yresid_lc)
"Yresid_lc"

#' @title The residual
#' @description The residual was computed by \code{\link{YResidualH1LL}}
#' @usage data(Yresid_ll)
"Yresid_ll"

#' @title The bandwidth
#' @description The bandwidth was computed by \code{\link{CVHstdLC}}
#' @usage data(hstd_lc)
"hstd_lc"

#' @title The Estimator of Diagonal Entries
#' @description The estimator of diagonal entries was computed by \code{\link{SigmaHstdLC}}
#' @usage data(Sigma_hstd_lc)
"Sigma_hstd_lc"

#' @title The standardization of Y
#' @description The standardization of Y
#' @usage data(Ystd_LC)
"Ystd_LC"


#' @title The bandwidth
#' @description The bandwidth was computed by \code{\link{CVHstdLL_LS}}
#' @usage data(hstd_LL_LS)
"hstd_LL_LS"

#' @title The Estimator of Diagonal Entries
#' @description The estimator of diagonal entries was computed by \code{\link{SigmaHstdLL_LS}}
#' @usage data(Sigma_LL_LS)
"Sigma_LL_LS"

#' @title The standardization of Y
#' @description The standardization of Y
#' @usage data(Ystd_LL_LS)
"Ystd_LL_LS"


#' @title The bandwidth
#' @description The bandwidth was computed by \code{\link{CVHstdLL_LH}}
#' @usage data(hstd_LL_LH)
"hstd_LL_LH"

#' @title The Estimator of Diagonal Entries
#' @description The estimator of diagonal entries was computed by \code{\link{SigmaHstdLL_LH}}
#' @usage data(Sigma_LL_LH)
"Sigma_LL_LH"

#' @title The standardization of Y
#' @description The standardization of Y
#' @usage data(Ystd_LL_LH)
"Ystd_LL_LH"

#' @title The bandwidth
#' @description The bandwidth was computed by \code{\link{DCMCVh}}
#' @usage data(hdcm1)
"hdcm1"

#' @title The bandwidth
#' @description The bandwidth was computed by \code{\link{DCMCVhMu}}
#' @usage data(hdcm2)
"hdcm2"

