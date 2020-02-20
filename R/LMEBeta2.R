#' @title LMEBeta2
#' @description This function is used to compute the delete-i beta goal
#' function f = (sbuk / sbk - suk / sk)^2;
#'
#' @param v the beta when the i-th observation deleted
#' @param y2 the square of y
#' @param expu the exp(u_j - u_i)
#' @param k the kernel
#' @param uk the (u_j - u_i) * kernel
#' @param SK summation of kernel
#' @param SUK summation of kernel * u
#'
#' @return goal function value

LMEBeta2 <- function(v, y2, expu, k, uk, SK, SUK) {
    B <- y2 / (expu^v)
    SBK <- B %*% k
    SBUK <- B %*% uk
    f <- (SBUK / SBK - SUK / SK)^2
    return(f)
}
