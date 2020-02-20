#' @title Optimise_LMEBeta2
#' @description  optimise function LMEBeta2
#'
#' @inheritParams LMEBeta2
#'
#' @return the optimal bandwidth

Optimise_LMEBeta2 <- function(y2, expu, k, uk, SK, SUK) {
    v <- optimise(LMEBeta2, c(-1000, 1000), tol = 1e-6,
             y2 = y2, expu = expu, k = k, uk = uk, SK = SK, SUK = SUK)
    return(v$minimum)
}