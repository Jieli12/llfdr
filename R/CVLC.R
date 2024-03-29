#' CVLC
#'
#' This function compute the cross validation value of local constant estimator
#' using the Gaussian kernel. We trimmed  5% percent of upper and lower
#' observations to avoid boundary effects.
#'
#' @param Y the observation, p * n matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @param h the bandwidth, scalar
#'
#' @return the cross validation value
#' @export
#'
#' @seealso \code{\link{CVLL}}
#' @examples \dontrun{
#'
#' data(Y)
#' data(u)
#' data(LowerBoundary)
#' # h is fixed
#' h_grid <- matrix(seq(c(LowerBoundary), 0.3, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, CVLC, Y = Y, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # find the optimal bandwidth
#' h1_lc <- optimise(CVLC, c(LowerBoundary, 2), tol = 1e-6, Y = Y, u = u)
#' abline(v = h1_lc$minimum, col="red")
#' }
CVLC <- function(Y, u, h, ktype = 'gaussian') {
    p <- nrow(Y)
    n <- ncol(Y)
    Start <- trunc(0.05 * n)
    End <- trunc(0.95 * n)
    Res <- matrix(0, nrow = p, ncol = End - Start + 1)
    u_sort <- sort(u, index.return = TRUE)
    u <- u_sort$x
    Y <- Y[,u_sort$ix]
    if (is.vector(Y)) {
        Y <- matrix(Y, nrow = 1)
    }
    # u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    # x <- (t(u_mat) - u_mat) / h
    # kernel <- dnorm(x) / h
    # sumk <- rowSums(kernel)
    # weight <- kernel / replicate(n, sumk)
    # alpha <- sumk / (sumk - diag(kernel))
    weight <- kernel_weight(u, bw = h, ktype = ktype)
    alpha <- 1 / (1 - diag(weight))
    for (i in Start:End) {
        Res[,i-Start+1] <- alpha[i] * (Y[,i] - Y %*% weight[i,])
    }
    cv = sum(sum(Res^2)) / (n * p)
    return(cv)
}
