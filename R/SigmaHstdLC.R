#' @title SigmaHstdLC
#'
#' @description  this routine computes the corresponding sigma matrix
#' with standard bandwidth at each u_i.
#'
#' @param Y the observation matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @param h the bandwidth, scalar
#'
#' @return the estimator of  diagonal entries at each given u_i.
#' @export
#' @seealso \code{\link{CVLC}},\code{\link{CVHstdLC}}
#' @examples \dontrun{
#'
#' data(Yresid_lc)
#' data(u)
#' data(hstd_lc)
#' # compute the estimator of diagonal entries
#' Sigma_hstd_lc <- SigmaHstdLC(Y = Yresid_lc, u = u, h = hstd_lc$minimum)
#' Ystd_LC <- Yresid_lc / sqrt(Sigma_hstd_lc)
#'}
SigmaHstdLC  <- function(Y, u, h, ktype = 'gaussian') {
    # p <- nrow(Y)
    # n <- ncol(Y)
    # u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    # x <- (t(u_mat) - u_mat) / h
    # kernel <- dnorm(x) / h
    # sumk <- rowSums(kernel)
    # weight <- kernel / replicate(n, sumk)
    weight <- kernel_weight(u, bw = h, ktype = ktype)
    Y2 <- Y^2
    # sigma <- apply(weight, 1, .Product, Y = Y2)
    sigma <- apply(weight, 1, tcrossprod, y = Y2)
    return(sigma)
}
