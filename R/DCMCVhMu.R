#' @title DCMCVhMu
#' @description This routine calculate the cross validation value without
#' threshold andshrinkage. We use subset-y-variables method.
#'
#' @param Y the observation matrix
#' @param u the condtion
#' @param h the bandwidth, scalar
#' @param Index.y the index of subset-y-variables
#' @param Y_res_h1 the residual
#'
#' @return the value of cross validation
#' @export
#' @seealso \code{\link{CVLC}}
#' @examples \dontrun{
#'
#' data(u)
#' data(LowerBoundary)
#' data(Yresid_lc)
#' data(Index.y)
#' data(Y)
#' h_grid <- matrix(seq(0.1, 1, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, DCMCVhMu, Y = Y, u = u,
#'             Y_res_h1 = Yresid_lc, Index.y = Index.y)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' hdcm2 <- optimise(DCMCVhMu, c(LowerBoundary, 2), tol = 1e-6,
#'                    Y = Y, u = u, Y_res_h1 = Yresid_lc, Index.y = Index.y)
#' abline(v = hdcm2$minimum, col="red")
#'}
DCMCVhMu <- function(Y, u, h, Y_res_h1, Index.y) {
    n <- ncol(Y)
    N <- ncol(Index.y)
    result <- YAverage(Y, u, h)
    Y_residual <- Y - result$aver
    cv <- 0
    for (i in 1:n) {
        Cov <- CovDeleteDCM(Y, result$aver[,i], Y_residual[,i],
                            result$weight[i,], result$alpha[i])
        for (j in 1:N) {
            Indj <- Index.y[, j]
            Cov_indj <- Cov[Indj, Indj]
            Y_indj <- Y_res_h1[Indj, i]
            cv <- cv + Y_indj %*% solve(Cov_indj, Y_indj) + log(det(Cov_indj))
        }
    }
    return(cv / (N * n))
}

