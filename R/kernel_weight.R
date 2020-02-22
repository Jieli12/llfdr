#' @title Generate the kernel weight matrix
#'
#' @description This function generated the kernel weight matrix depending on
#' different parameter settings.
#' @inheritParams computeUdiff
#' @inheritParams kernelCompute
#' @param poly_order 0 or 1, 0 represents the local constant smoother
#'                  while 1 represents the local linear smoother, default is 0.
#'
#' @return the weight matrix
#' @seealso \code{\link{computeUdiff}}, \code{\link{kernelCompute}}
#' @examples \dontrun{
#'
#' n <- 100
#' n0 <- 10
#' u <- seq(-1, 1, length.out = n)
#' u0 <- seq(-0.9, 0.9, length.out = n0)
#' # Cosine kernel with bandwidth 0.4, local constant smoother
#' weight <- kernel_weight(u, u0, ktype = 'cosine', bw = 0.4, poly_order = 0)
#' }
kernel_weight <- function(u, u0 = NULL, ktype = 'gaussian',
                          bw = 0.1, poly_order = 0) {
    u_diff <- computeUdiff(u, u0)
    kernel <- kernelCompute(u_diff, ktype = ktype, bw = bw)
    if (poly_order == 0) {
        sumk <- rowSums(kernel)
        weight <- kernel / replicate(n1, sumk)
    } else if (poly_order == 1) {
        u_diff_sq <- u_diff^2
        sn1 <- rowSums(kernel * u_diff)
        sn2 <- rowSums(kernel * u_diff_sq)
        Sn1 <- replicate(n1, sn1)
        Sn2 <- replicate(n1, sn2)
        weight <- kernel * (Sn2 - u_diff * Sn1)
        weight_sum = rowSums(weight)
        weight = weight / replicate(n1, weight_sum)
    } else {
        stop("Argument 'poly_order' must be 0 or 1.")
    }
    return(weight)
}








