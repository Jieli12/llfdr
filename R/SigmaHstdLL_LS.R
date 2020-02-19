#' @title SigmaHstdLL_LS
#'
#' @description  this routine computes the corresponding sigma matrix
#' with standard bandwidth at each u_i using local linear smoother.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#'
#' @return the estimator of  diagonal entries at each given u_i.
#' @export
#' @seealso \code{\link{CVLL}},\code{\link{CVHstdLL}}
#' @examples \dontrun{
#'
#' data(YResidual_LL)
#' data(u)
#' data(h_std_ll_ls)
#' # compute the estimator of diagonal entries
#' Sigma_hstd_ll_ls <- SigmaHstdLL_LS(Y = YResidual_LL, u = u, h = h_std_ll$minimum)
#'}
SigmaHstdLL_LS  <- function(Y, u, h) {
    # p <- nrow(Y)
    # n <- ncol(Y)
    # u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    # x <- (t(u_mat) - u_mat) / h
    # kernel <- dnorm(x) / h
    # sumk <- rowSums(kernel)
    # weight <- kernel / replicate(n, sumk)
    weight = kernel_weight(u, bw = h, poly_order = 1)
    Y2 = Y^2
    # sigma <- apply(weight, 1, .Product, Y = Y2)
    sigma <- apply(weight, 1, tcrossprod, y = Y2)
    sigma <- apply(sigma, 2, interpolate)
    return(sigma)
}
