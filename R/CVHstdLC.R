#' @title CVHstdLC
#'
#' @description  This function calculates the cross validation values using local constant
#' method when bandwidth is given.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#'
#' @return the cross validation value
#' @export
#' @seealso \code{\link{CVLC}}
#' @examples \dontrun{
#'
#' data(u)
#' data(LowerBoundary)
#' data(Yresid_lc)
#' h_grid <- matrix(seq(0.2, 1, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, CVHstdLC, Y = Yresid_lc, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' hstd_lc <- optimise(CVHstdLC, c(LowerBoundary, 2), tol = 1e-6,
#'                    Y = Yresid_lc, u = u)
#' abline(v = hstd_lc$minimum, col="red")
#'}
CVHstdLC  <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    # u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    # x <- (t(u_mat) - u_mat) / h
    # kernel <- dnorm(x) / h
    # sumk <- rowSums(kernel)
    # weight <- kernel / replicate(n, sumk)
    weight = kernel_weight(u, bw = h, poly_order = 0)
    cv <- 0
    for (i in 1:n) {
        w <- weight[i, -i]
        w <- w / sum(w)
        sigma_i <- Y[,-i]^2 %*% w
        cv <- cv + sum(Y[,i]^2 / sigma_i + log(sigma_i))
    }
    cv <- cv / n / p
    return(cv)
}
