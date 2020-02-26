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
#' @seealso \code{\link{DCMCVhMu}}, \code{\link{Cov_EstMu}}
#' @examples \dontrun{
#'
#' data(u)
#' data(Y)
#' data(Index.obs)
#' data(hdcm2)
#' u <- as.vector(u)
#'
#' u0_len <- 25
#' u0 <- seq(-0.9, 0.9, length.out = u0_len)
#' a <- 0
#' bcoef <- 0.6
#' thdcm2  <- rep(0, u0_len)
#' for (i in 1:u0_len) {
#'     b2 <- bcoef * max(Cov_EstMu(Y = Y, u = u, u0 = u0[i], h = hdcm2$minimum))
#'     thdcm2_temp <- optimise(DCMCV_ThMu, c(a, b2), tol = 1e-6, Yi = Y, u = u,
#'                   Index.obs = Index.obs, u0 = u0[i], h = hdcm2$minimum)
#'     thdcm2[i] <-thdcm2_temp$minimum
#' }
#' }
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
