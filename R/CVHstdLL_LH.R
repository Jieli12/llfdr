#' @title  CVHstdLL_LH
#' @description This function is used to compute the delete-i cross validation value
#' @param Y the observation, p * n matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param h the bandwidth, scalar
#'
#' @return the value of cross validation
#' @export
#' @seealso \code{\link{CVLL}}
#' @examples \dontrun{
#'
#' data(Yresid_ll)
#' data(u)
#' data(LowerBoundary)
#' h_grid <- matrix(seq(0.5, 2, length.out = 100), nrow = 100)
#' cv <- apply(h_grid, 1, CVHstdLL_LH, Y = Yresid_ll, u = u)
#' plot(h_grid,cv, type = 'l', xlab = "Bandwidth", ylab = "CV Values", col = "blue")
#' # select the optimal bandwidth for diagonal entries of covariance
#' hstd_LL_LH <- optimise(CVHstdLL_LH, c(LowerBoundary, 2), tol = 1e-6,
#'                    Y = Yresid_ll, u = u)
#' abline(v = hstd_LL_LH$minimum, col="red")
#'}
CVHstdLL_LH <- function(Y, u, h, ktype = 'gaussian') {
    p <- nrow(Y)
    n <- ncol(Y)

    U_diff <- computeUdiff(u)
    K <- kernelCompute(U_diff, ktype = ktype, bw = h)
    UK_mat <- U_diff * K
    Y2 <- Y^2
    ExpU_mat <- exp(U_diff)
    cv <- 0
    for (i in 1:n) {
        expu_del_i <- ExpU_mat[i, -i]
        k_del_i <- K[i, -i]
        uk_del_i <- UK_mat[i, -i]
        mat_uk <- rbind(k_del_i,uk_del_i)
        SK <- sum(k_del_i)
        r0 <- sum(uk_del_i) / SK
        Y2_i <- Y2[, -i]
        # v <- rep(0, p)
        # for (j in 1:p) {
        #     y2_del_i <- Y2_i[j, ]
        #     v[j] <- Optimise_LMEBeta2(y2 = y2_del_i, expu = expu_del_i,
        #                            mat_uk = mat_uk, r0 = r0)
        # }
        v <- apply(Y2_i, 1, Optimise_LMEBeta2, expu = expu_del_i,
              mat_uk = mat_uk, r0 = r0)
        B <- Y2_i / t(outer(expu_del_i, v, '^'))
        SBK <- B %*% k_del_i
        Ai <- SBK / SK
        cv <- cv + sum(Y2[,i] / Ai + log(Ai))
    }
    return(cv)
}




