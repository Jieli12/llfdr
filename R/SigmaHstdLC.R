#' @title SigmaHstdLC
#'
#' @description  this routine computes the corresponding sigma matrix
#' with standard bandwidth at each u_i.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#'
#' @return the estimator of  diagonal entries at each given u_i.
#' @export
#' @seealso \code{\link{CVLC}},\code{\link{CVHstdLC}}
#' @examples \dontrun{
#'
#' data(Y)
#' data(u)
#' data(LowerBoundary)
#' # find the optimal bandwidth for mean estimation
#' h_mean <- optimise(CVLC, c(LowerBoundary, 2), tol = 1e-6, Y = Y, u = u)
#' # compute the residual
#' YResidual_LC <- YResidualLC(Y = Y, u = u, h = h_mean$minimum)
#' h_grid <- matrix(seq(0.2, 1, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, CVHstdLC, Y = YResidual_LC, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' h_std <- optimise(CVHstdLC, c(LowerBoundary, 2), tol = 1e-6,
#'                    Y = YResidual_LC, u = u)
#' abline(v = h_std$minimum, col="red")
#'}
SigmaHstdLC  <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    x <- (t(u_mat) - u_mat) / h
    kernel <- dnorm(x) / h
    sumk <- rowSums(kernel)
    weight <- kernel / replicate(n, sumk)
    Y2 = Y^2
    sigma <- apply(weight, 1, .Product, Y = Y2)
    return(sigma)
}
