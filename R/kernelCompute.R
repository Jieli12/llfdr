#' @title kernelCompute
#'
#' @param u_diff a matrix of difference of u and u0
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#' @param bw bandwidth, default is 0.1
#'
#' @return kernel matrix
#'
#' @seealso \code{\link{computeUdiff}}, \code{\link{kernel_weight}}
#' @examples \dontrun{
#'
#' n <- 20
#' n0 <- 10
#' u <- seq(-1, 1, length.out = n)
#' u0 <- seq(-0.9, 0.9, length.out = n0)
#' u_diff <- computeUdiff(u, u0)
#' # Gaussian kernel with bandwidth 0.1, default
#' kernel <- kernelCompute(u_diff)
#' # Cosine kernel with bandwidth 0.4
#' kernel_cos <- kernelCompute(u_diff, ktype = 'cosine', bw = 0.4)
#' }
kernelCompute <- function(u_diff,  ktype = 'gaussian', bw = 0.1) {
    diffh <- u_diff / bw
    abs_diffh <- abs(diffh)
    kernel <- switch(ktype,
                     gaussian = dnorm(diffh) / bw,
                     epanech = {
                         ifelse(abs_diffh <= 1, 3 / 4 * ( 1 - diffh^2) / bw, 0)
                     },
                     triweight = {
                         ifelse(abs_diffh <= 1, 35 / 32 * ( 1 - diffh^2)^3 / bw, 0)
                     },
                     biweight = {
                         ifelse(abs_diffh <= 1, 15 / 16 * ( 1 - diffh^2)^2 / bw, 0)
                     },
                     tricube = {
                         ifelse(abs_diffh <= 1, 70 / 81 * ( 1 - abs_diffh^3)^3 / bw, 0)
                     },
                     triangular = {
                         ifelse(abs_diffh <= 1, ( 1 - abs_diffh) / bw, 0)
                     },
                     cosine = {
                         ifelse(abs_diffh <= 1, pi / 4 * cos( pi / 2 * diffh) / bw , 0)
                     }
    )
    return(kernel)
}
