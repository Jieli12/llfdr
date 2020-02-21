#' @title CV_Th
#' @description This routine computes the goal function for optimal threshold value
#'
#' @param Yi the p * n matrix
#' @param u the condition u
#' @param Index.obs the index of observation
#' @param u_test the given u
#' @param h the optimal bandwidth
#' @param r the threshold value
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#'
#' @return the value of cross validation function
CV_Th <- function(Yi, u, Index.obs, u_test, h, r, ktype = 'gaussian') {
    N1 <- ncol(Index.obs)
    cv  <- 0
    p <- nrow(Yi)
    for (i in 1:N1) {
        Indi <- Index.obs[, i]
        Ys1 <- Yi[, Indi]
        Ys2 <- Yi[, -Indi]
        u1 <- u[Indi]
        u2 <- u[-Indi]
        Covs1 <- Cov_Est(Ys1, u1, u_test, h, ktype)
        Covs2 <- Cov_Est(Ys2, u2, u_test, h, ktype)
        Covs1 <- ifelse(abs(Covs1) >= r, Covs1, 0)
        cv <- cv + norm(Covs1 - Covs2, "F")^2 / p
    }
    return(cv / N1)
}
