#' @title CovDelete
#' @description This routine computes the covariance when the ith observation
#' is deleted atcenter u(i).
#'
#' @param Y_del the matrix when the i-th column is deleted
#' @param weight_del the vector when the i-th entry is deleted
#'
#' @return covariance estimator
#' @export
CovDelete <- function(Y_del, weight_del) {
    tem <- sum(weight_del)
    if (tem == 0) {
        return(matrix(0, nrow = nrow(Y_del), ncol = nrow(Y_del)))
    } else {
        w <- weight_del / tem
        Cov <- tcrossprod(Y_del %*% diag(w), Y_del)
        return(Cov)
    }
}
