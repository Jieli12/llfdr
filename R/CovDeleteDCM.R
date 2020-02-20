#' @title CovDeleteDCM
#' @description This routine computes the covariance when the ith observation is
#'  deleted at center u(i). We use some algebra techniques to calculate the
#'  covariance estimation atcenter u(i) when the ith observation is deleted.
#'  This can avoid the loop of deleting ith observation, which can reduce
#'  the computation complexity. Actually, we can use the former computation
#'  results to estimate the covariance avoiding the computation of new dataset
#'  without ith observation. Note that, we use the covariance format in DCM.
#'  See details in Dynamiccovariance model.
#'
#' @param Y the p*n size matrix
#' @param Aver_i the Average of Y_i at center u(i), it comes from the results
#'  of function Average, p * 1 matrix
#' @param Residual_i the residual of Y_i - Aver_i, this can be calculted
#' outside the loop.
#' @param Weight_i the 1 * n weight matrix, it represents the weight for
#' the u(i), it is a sparse matrix.
#' @param Alpha_i the coefficient of alpha(i,u(i)), see my technical reports
#' for details.
#'
#' @return the covariance when deleting the ith observation

CovDeleteDCM <- function(Y, Aver_i, Residual_i, Weight_i, Alpha_i) {
    Cov1 <- tcrossprod(Y %*% diag(Weight_i), Y) - tcrossprod(Aver_i)

     # The residuals
    Cov2 <- tcrossprod(Residual_i)

    # Compute the covariance when the i-th observation is deleted
    Cov <- Alpha_i * Cov1 + Alpha_i * (1 - Alpha_i) * Cov2
    return(Cov)
}
