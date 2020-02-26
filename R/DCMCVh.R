#' @title DCMCVh
#' @description This routine calculate the cross validation value without
#' threshold andshrinkage. We use subset-y-variables method.
#'
#' @param Y the observation matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth, scalar
#' @param Index.y the index of subset-y-variables
#'
#' @return the value of cross validation
#' @export
#' @seealso \code{\link{CVLC}}
#' @examples \dontrun{
#'
#' data(u)
#' data(LowerBoundary)
#' data(Yresid_lc)
#' data(Index.y)
#' h_grid <- matrix(seq(0.2, 1, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, DCMCVh, Y = Yresid_lc, u = u, Index.y = Index.y)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' hdcm1 <- optimise(DCMCVh, c(LowerBoundary, 2), tol = 1e-6,
#'                    Y = Yresid_lc, u = u , Index.y = Index.y)
#' abline(v = hdcm1$minimum, col="red")
#'}
DCMCVh <- function(Y, u, h, Index.y, ktype = 'gaussian') {
    n <- ncol(Y)
    N <- ncol(Index.y)
    cv <- 0
    weight <- kernel_weight(u, bw = h, poly_order = 0, ktype = ktype)
    for (i in 1:n) {
        Cov <- CovDelete(Y[,-i], weight[i,-i])
        for (j in 1:N) {
            Indj <- Index.y[, j]
            Cov_indj <- Cov[Indj, Indj]
            Y_indj <- Y[Indj, i]
            cv <- cv + Y_indj %*% solve(Cov_indj, Y_indj) + log(det(Cov_indj))
        }
    }
    return(cv / (N * n))
}

