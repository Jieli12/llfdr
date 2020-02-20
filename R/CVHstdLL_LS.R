#' @title CVHstdLL_LS
#'
#' @description  This function calculates the cross validation values using
#' local linear smoother when bandwidth is given.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#'
#' @return the cross validation value
#' @export
#' @seealso \code{\link{CVLL}}
#' @examples \dontrun{
#'
#' data(Yresid_ll)
#' data(u)
#' data(hstd_ll)
#' h_grid <- matrix(seq(1.5, 2, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, CVHstdLL_LS, Y = Yresid_ll, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' hstd_LL_LS <- optimise(CVHstdLL_LS, c(LowerBoundary, 2), tol = 1e-6,
#'                    Y = Yresid_ll, u = u)
#' abline(v = hstd_LL_LS$minimum, col="red")
#'}
CVHstdLL_LS  <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    weight = kernel_weight(u, bw = h, poly_order = 1)
    Y_sq = Y^2
    cv <- 0
    for (i in 1:n) {
        w <- weight[i, -i]
        w <- w / sum(w)
        sigma_i <- Y_sq[,-i] %*% w
        cv <- cv + sum((Y_sq[, i] - sigma_i)^2)
    }
    cv <- cv / n / p
    return(cv)
}
