#' @title LMEBeta2
#' @description This function is used to compute the delete-i beta goal
#' function f = sbuk - sbk * r0
#'
#' @param v the beta when the i-th observation deleted
#' @param y2 the square of y
#' @param expu the exp(u_j - u_i)
#' @param mat_uk the matrix, the first row is k, the secound row is u * k
#' @param r0 the ratio of SUK to SK
#' @export
#' @return goal function value

LMEBeta2 <- function(v, y2, expu, mat_uk, r0) {
    B <- y2 / (expu^v)
    C <- mat_uk %*% B
    return(C[2] - r0 * C[1])
}
