#' @title YAverage
#'
#' @description This function is used for estimating the local constant
#' kernel weighted average at each center u(i) given bandwidth h.
#' Here we use gaussian kerenel.
#'
#' @param Y the observation, p * n matrix
#' @param u the condtion, 1 * n matrix
#' @param h the bandwidth, scalar
#'
#' @return This function will return list object which contains 3 variables:
#'
#' \item{aver}{each column is the weighted average of Y given the condition u(i) and h. It's a p*n size matrix.}
#' \item{weight}{the n * n weight matrix, the ith row represents the weight for the u(i), it is also a sparse matrix.}
#' \item{alpha}{the coefficient adjusted after the i-th observation deleted.}
#'
#' @export
#'
#' @seealso \code{\link{CVLL}}
#' @examples \dontrun{
#'
#' data(Y)
#' data(u)
#' data(LowerBoundary)
#' # h is fixed
#' h <- 0.1
#' result <- YAverage(Y, u, h)
#' }
YAverage <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    u_mat <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    x <- (t(u_mat) - u_mat) / h
    kernel <- dnorm(x) / h
    sumk <- rowSums(kernel)
    weight <- kernel / replicate(n, sumk)
    alpha <- sumk / (sumk - diag(kernel))
    # Product <- function(Y = Y, w) {
    #     return(Y %*% w)
    # }
    # aver <- apply(weight, 1, Product, Y = Y)
    aver <- apply(weight, 1, tcrossprod, y = Y)
    result <- list(aver = aver, weight = weight, alpha = alpha)
    return(result)
}
