#' CVHstdLC
#'
#' This function calculates the cross validation values using local constant
#' method when bandwidth is given.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#'
#' @return the cross validation value
#' @export
#'
#' @examples
CVHstdLC  <- function(Y, u, h) {
    p <- nrow(Y)
    n <- ncol(Y)
    u_mat <- replicate(n, u)
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
