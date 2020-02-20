#' @title  CVHstdLL_LH
#' @description This function is used to compute the delete-i cross validation value
#' @inheritParams kernel_weight
#' @inheritParams CVll
#'
#' @return the value of cross validation

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
    U <- matrix(rep(t(u), n), ncol = ncol(u), byrow = TRUE)
    U_diff <- U - t(U)
    U_diff_H <- U_diff / h
    abs_diffh <- abs(U_diff_H)
    K <- switch(ktype,
                gaussian = dnorm(U_diff_H) / h,
                epanech = {
                    ifelse(abs_diffh <= 1, 3 / 4 * ( 1 - U_diff_H^2) / h, 0)
                },
                triweight = {
                    ifelse(abs_diffh <= 1, 35 / 32 * ( 1 - U_diff_H^2)^3 / h, 0)
                },
                biweight = {
                    ifelse(abs_diffh <= 1, 15 / 16 * ( 1 - U_diff_H^2)^2 / h, 0)
                },
                tricube = {
                    ifelse(abs_diffh <= 1, 70 / 81 * ( 1 - U_diff_H^3)^3 / h, 0)
                },
                triangular = {
                    ifelse(abs_diffh <= 1, ( 1 - U_diff_H) / h, 0)
                },
                cosine = {
                    ifelse(abs_diffh <= 1, pi / 4 * cos( pi / 2 * U_diff_H) / h , 0)
                }
    )
    UK_mat <- U_diff * K
    Y2 <- Y^2
    ExpU_mat <- exp(U_diff)
    cv <- 0
    for (i in 1:n) {
        expu_del_i <- ExpU_mat[i, -i]
        k_del_i <- K[i, -i]
        uk_del_i <- UK_mat[i, -i]
        SK <- sum(k_del_i)
        SUK <- sum(uk_del_i)
        Y2_i <- Y2[, -i]
        v <- rep(0, p)
        for (j in 1:p) {
            y2_del_i <- Y2_i[j, ]
            v[j] <- Optimise_LMEBeta2(y2 = y2_del_i, expu = expu_del_i,
                                   k = k_del_i, uk = uk_del_i,
                                   SK = SK, SUK = SUK)
        }
        B <- Y2_i / t(outer(expu_del_i, v, '^'))
        SBK <- B %*% k_del_i
        Ai <- SBK / SK
        cv <- cv + sum(Y2[,i] / Ai + log(Ai))
    }
    return(cv)
}




