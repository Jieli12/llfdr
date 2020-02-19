#' @title Generate the kernel weight matrix
#'
#' @description This function generated the kernel weight matrix depending on
#' differe parameter settings.
#' @param u the condition, row matrix, size 1 * n1.
#' @param u0 the center condition, row matrix, size 1 * n2, the default is u.
#' @param ktype the kernel type, can be "gaussian", "epanech", "triweight",
#'              "biweight", "tricube", "triangular" and "cosine",
#'              the default of ktype is "gaussian".
#' @param bw bandwidth, default is 0.1
#' @param poly_order 0 or 1, 0 represents the local constant smoother
#'                  while 1 represents the local linear smoother, default is 0.
#'
#' @return the weight matrix
#' @export
#'
#' @examples
kernel_weight <- function(u, u0 = u, ktype = 'gaussian',
                          bw = 0.1, poly_order = 0) {
    if (!is.matrix(u)){
        stop("argument 'u' must be row matrix")
    }
    if (!is.matrix(u0)){
        if (length(u0) == 1){
            u0 = matrix(u0)
        } else{
            stop("argument 'u0' must be row matrix")
        }
    }
    n1 <- length(u)
    n2 <- length(u0)
    u_rep <- rep(t(u), times = n2)
    u0_rep <- rep(t(u0), each = n1)
    u_diff <- matrix(u_rep - u0_rep, ncol = ncol(u), byrow = TRUE)
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
    if (poly_order == 0){
        sumk <- rowSums(kernel)
        weight <- kernel / replicate(n1, sumk)
    } else {
        u_diff_sq <- u_diff^2
        sn1 <- rowSums(kernel * u_diff)
        sn2 <- rowSums(kernel * u_diff_sq)
        Sn1 <- replicate(n1, sn1)
        Sn2 <- replicate(n1, sn2)
        weight <- kernel * (Sn2 - u_diff * Sn1)
        weight_sum = rowSums(weight)
        weight = weight / replicate(n1, weight_sum)
    }
    return(weight)
}








