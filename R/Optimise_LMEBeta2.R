#' @title Optimise_LMEBeta2
#' @description  optimise function LMEBeta2
#'
#' @inheritParams LMEBeta2
#'
#' @return the optimal bandwidth

#' @export

Optimise_LMEBeta2 <- function(y2, expu, mat_uk, r0) {
    v <- uniroot(LMEBeta2, c(-100, 100), tol = 1e-6,
             y2 = y2, expu = expu, mat_uk = mat_uk, r0 = r0)
    return(v$root)
}
