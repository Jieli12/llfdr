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
#' data(Y)
#' data(u)
#' data(LowerBoundary)
#' # find the optimal bandwidth for mean estimation
#' h_mean_lc <- optimise(CVLC, c(LowerBoundary, 2), tol = 1e-6, Y = Y, u = u)
#' # compute the residual
#' YResidual_LC <- YResidualLC(Y = Y, u = u, h = h_mean_lc$minimum)
#' h_grid <- matrix(seq(0.2, 1, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, CVHstdLC, Y = YResidual_LC, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' h_std_lc <- optimise(CVHstdLC, c(LowerBoundary, 2), tol = 1e-6,
#'                    Y = YResidual_LC, u = u)
#' abline(v = h_std_lc$minimum, col="red")
#'}
CVHstdLC  <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    x <- (t(u_mat) - u_mat) / h
    kernel <- dnorm(x) / h
    sumk <- rowSums(kernel)
    weight <- kernel / replicate(n, sumk)
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
