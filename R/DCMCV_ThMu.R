#' @title DCMCV_ThMu
#' @description This routine computes the goal function for optimal threshold value
#'
#' @param Yi the p * n matrix
#' @param Index.obs the index of observation
#' @param weight the weight at u0
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
#' Weight <- kernel_weight(u, u0 = u0, ktype = 'gaussian', bw = hdcm2$minimum)
#' for (i in 1:u0_len) {
#'     b2 <- bcoef * max(Cov_EstMu(Y = Y, weight = Weight[i,]))
#'     thdcm2_temp <- optimise(DCMCV_ThMu, c(a, b2), tol = 1e-6, Yi = Y,
#'                   weight = Weight[i,], Index.obs = Index.obs)
#'     thdcm2[i] <-thdcm2_temp$minimum
#' }
#' }
DCMCV_ThMu <- function(Yi, weight, Index.obs, r) {
    N1 <- ncol(Index.obs)
    cv  <- 0
    p <- nrow(Yi)
    for (i in 1:N1) {
        Indi <- Index.obs[, i]
        Ys1 <- Yi[, Indi]
        Ys2 <- Yi[, -Indi]
        weight1 <- weight[Indi]
        weight1 <- weight1 / sum(weight1)
        weight2 <- weight[-Indi]
        weight2 <- weight2 / sum(weight2)
        Covs1 <- Cov_EstMu(Ys1, weight1)
        Covs2 <- Cov_EstMu(Ys2, weight2)
        # Covs1 <- ifelse(abs(Covs1) >= r, Covs1, 0)
        Covs1[abs(Covs1) < r] <- 0
        cv <- cv + norm(Covs1 - Covs2, "F")^2 / p
    }
    return(cv / N1)
}
