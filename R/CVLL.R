#' CVLL
#'
#' This function compute the cross validation value of local linear estimator
#' using the Gaussian kernel. We trimmed  5% percent of upper and lower
#' observations to avoid boundary effects.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#'
#' @return the cross validation value
#' @export
#'
#' @seealso \code{\link{CVLC}}
#' @examples \dontrun{
#'
#' data(Y)
#' data(u)
#' data(LowerBoundary)
#' # h is fixed
#' h_grid <- matrix(seq(LowerBoundary, 0.5, length.out = 100),nrow = 100)
#' cv <- apply(h_grid, 1, CVLL, Y = Y, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' }
CVLL <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    Start <- trunc(0.05 * n)
    End <- trunc(0.95 * n)
    Res <- matrix(0, nrow = p, ncol = End - Start + 1)
    u_sort <- sort(u, index.return = TRUE)
    u <- u_sort$x
    Y <- Y[,u_sort$ix]

    u_mat <- replicate(n, u)
    u_diff <- t(u_mat) - u_mat
    x <- u_diff / h
    kernel <- dnorm(x) / h
    u_diff_2 <- u_diff^2
    sn1 <- rowSums(kernel * u_diff)
    sn2 <- rowSums(kernel * u_diff_2)
    Sn1 <- replicate(n, sn1)
    Sn2 <- replicate(n, sn2)
    Weight <- kernel * (Sn2 - u_diff * Sn1)
    Weight_sum = rowSums(Weight)
    Weight = Weight / replicate(n, Weight_sum)
    for (i in Start:End) {
        Res[,i-Start+1] <- Y[,i] - Y[,-i] %*% Weight[i,-i]
    }
    cv = sum(sum(Res^2)) / (n * p)
    return(cv)
}
