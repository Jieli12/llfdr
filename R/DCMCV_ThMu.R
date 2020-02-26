#' @title DCMCV_ThMu
#' @description This routine computes the goal function for optimal threshold value
#'
#' @param Yi the p * n matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param Index.obs the index of observation
#' @param h the  bandwidth
#' @param r the threshold value
#' @export
#' @return the value of cross validation function
DCMCV_ThMu <- function(Yi, u, Index.obs, u0, h, r, ktype = 'gaussian') {
    N1 <- ncol(Index.obs)
    cv  <- 0
    p <- nrow(Yi)
    for (i in 1:N1) {
        Indi <- Index.obs[, i]
        Ys1 <- Yi[, Indi]
        Ys2 <- Yi[, -Indi]
        u1 <- u[Indi]
        u2 <- u[-Indi]
        Covs1 <- Cov_EstMu(Ys1, u1, u0, h, ktype)
        Covs2 <- Cov_EstMu(Ys2, u2, u0, h, ktype)
        Covs1 <- ifelse(abs(Covs1) >= r, Covs1, 0)
        cv <- cv + norm(Covs1 - Covs2, "F")^2 / p
    }
    return(cv / N1)
}
