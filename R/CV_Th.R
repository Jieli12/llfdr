#' @title CV_Th
#' @description This routine computes the goal function for optimal threshold value
#'
#' @param Yi the p * n matrix
#' @inheritParams kernelCompute
#' @inheritParams computeUdiff
#' @inheritParams kernel_weight
#' @param Index.obs the index of observation
#' @param h the  bandwidth
#' @param r the threshold value
#'
#' @return the value of cross validation function
#' @export
#' @seealso \code{\link{DCMCVh}}, \code{\link{Cov_Est}}
#' @examples \dontrun{
#'
#' data(u)
#' data(Yresid_lc)
#' data(Index.obs)
#' data(hdcm1)
#' u <- as.vector(u)
#'
#' u0_len <- 25
#' u0 <- seq(-0.9, 0.9, length.out = u0_len)
#' a <- 0
#' bcoef <- 0.6
#' thdcm1  <- rep(0, u0_len)
#' for (i in 1:u0_len) {
#'     b1 <- bcoef * max(Cov_Est(Y = Yresid_lc, u = u, u0 = u0[i], h = hdcm1$minimum))
#'     thdcm1_temp <- optimise(CV_Th, c(a, b1), tol = 1e-6, Yi = Yresid_lc, u = u,
#'                   Index.obs = Index.obs, u0 = u0[i], h = hdcm1$minimum)
#'     thdcm1[i] <-thdcm1_temp$minimum
#' }
#'
#' # NCM, threshold value selection
#' data(Ystd_LC)
#' data(hncm)
#' thncm  <- rep(0, u0_len)
#' for (i in 1:u0_len) {
#'     b1 <- bcoef * max(Cov_Est(Y = Ystd_LC, u = u, u0 = u0[i], h = hncm$minimum))
#'     thncm_temp <- optimise(CV_Th, c(a, b1), tol = 1e-8, Yi = Ystd_LC, u = u,
#'                   Index.obs = Index.obs, u0 = u0[i], h = hncm$minimum)
#'     thncm[i] <-thncm_temp$minimum
#' }
#' }
CV_Th <- function(Yi, u, Index.obs, u0, h, r, ktype = 'gaussian') {
    N1 <- ncol(Index.obs)
    cv  <- 0
    p <- nrow(Yi)
    for (i in 1:N1) {
        Indi <- Index.obs[, i]
        Ys1 <- Yi[, Indi]
        Ys2 <- Yi[, -Indi]
        u1 <- u[Indi]
        u2 <- u[-Indi]
        Covs1 <- Cov_Est(Ys1, u1, u0, h, ktype)
        Covs2 <- Cov_Est(Ys2, u2, u0, h, ktype)
        Covs1 <- ifelse(abs(Covs1) >= r, Covs1, 0)
        cv <- cv + norm(Covs1 - Covs2, "F")^2 / p
    }
    return(cv / N1)
}
