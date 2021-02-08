#' @title CV_Th
#' @description This routine computes the goal function for optimal threshold value
#'
#' @param Yi the p * n matrix
#' @param Index.obs the index of observation
#' @param weight the weight at u0
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
#' Weight <- kernel_weight(u, u0 = u0, ktype = 'gaussian', bw = hdcm1$minimum)
#' for (i in 1:u0_len) {
#'     b1 <- bcoef * max(Cov_Est(Y = Yresid_lc, Weight[i,]))
#'     thdcm1_temp <- optimise(CV_Th, c(a, b1), tol = 1e-6, Yi = Yresid_lc,
#'                   weight = Weight[i,], Index.obs = Index.obs)
#'     thdcm1[i] <-thdcm1_temp$minimum
#' }
#'
#' # NCM, threshold value selection
#' data(Ystd_LC)
#' data(hncm)
#' thncm  <- rep(0, u0_len)
#' Weight <- kernel_weight(u, u0 = u0, ktype = 'gaussian', bw = hncm$minimum)
#' for (i in 1:u0_len) {
#'     b1 <- bcoef * max(Cov_Est(Y = Ystd_LC, Weight[i,]))
#'     thncm_temp <- optimise(CV_Th, c(a, b1), tol = 1e-8, Yi = Ystd_LC,
#'                   weight = Weight[i,], Index.obs = Index.obs)
#'     thncm[i] <-thncm_temp$minimum
#' }
#' }
# CV_Th <- function(Yi, u, Index.obs, u0, h, r, ktype = 'gaussian') {
CV_Th <- function(Yi, weight, Index.obs, r) {
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
        Covs1 <- Cov_Est(Ys1, weight1)
        Covs2 <- Cov_Est(Ys2, weight2)
        # Covs1 <- ifelse(abs(Covs1) >= r, Covs1, 0)
        Covs1[abs(Covs1) < r] <- 0
        cv <- cv + norm(Covs1 - Covs2, "F")^2 / p
    }
    return(cv / N1)
}
