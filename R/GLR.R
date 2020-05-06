#' GLR
#' @description this routine computes the p value using generalized likelihood
#' ratio statistics
#'
#' @param res_linear the residual of linear model
#' @param res_nonlinear the residual of nonlinear model
#' @param h  the bandwidth
#' @param len_interval the range of u
#' @param n sample size
#'
#' @return the p value of chi square test
#' @export
#'
GLR <- function(res_linear, res_nonlinear, h, len_interval, n){
    lambda <- log(sum(res_linear^2) / sum(res_nonlinear^2)) * n / 2
    rk <- 2.53750769446292
    ck <- 0.25789488451449
    df <- rk * ck * len_interval / h
    pvalue <- pchisq(rk * lambda, df, lower.tail = FALSE)
    return(pvalue)
}
