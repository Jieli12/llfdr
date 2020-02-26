#' @title NCMCVh
#' @description This routine calculate the cross validation value using NCM method.
#' @param Y the observation, p * n matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth, scalar
#' @export
#' @return the value of cross validation function
#' @seealso \code{\link{CVLL}}
#' @examples \dontrun{
#'
#' data(u)
#' data(LowerBoundary)
#' data(Ystd_LC)
#' upper <- 1
#' h_grid <- matrix(seq(0.05, upper, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, NCMCVh, Y = Ystd_LC, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' hncm <- optimise(NCMCVh, c(LowerBoundary, upper), tol = 1e-6,
#'                    Y = Ystd_LC, u = u)
#' abline(v = hncm$minimum, col="red")
#'}
NCMCVh <- function(Y, u, h, ktype = 'gaussian') {
    p <- nrow(Y)
    n <- ncol(Y)
    weight <- kernel_weight(u, bw = h, poly_order = 0, ktype = ktype)
    cv <- 0
    for(i in 1:n) {
        covi <- CovDelete(Y[,-i], weight[i,-i])
        diff <- covi - tcrossprod(Y[,i])
        cv <- cv + norm(diff, "F")^2 / p
    }
    return(cv/n)
}
