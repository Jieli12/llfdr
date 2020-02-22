#' @title computeUdiff
#'
#' @param u the condition, it is a vector, one-dimensional array or
#' one-dimensional row(column) matrix
#' @param u0 the given condition to be estimated. The default is NULL.
#' If u0 is NULL, then u0 = u.
#' @details Supposing u and u0 are vectors with length n and n0 respectively,
#' this function will return a n0*n-size matrix, the i-th row of it is equal to
#' u - u0_i
#' @return a matrix of difference of u and u0
#' @seealso \code{\link{kernel_weight}}
#' @examples \dontrun{
#'
#' u <- 1:5
#' # Default u0
#' computeUdiff(u)
#' # Given u0
#' u0 <- seq(1.5, 4.5, by = 1)
#' computeUdiff(u, u0)

#'}
computeUdiff <- function(u, u0 = NULL) {
    if (is.matrix(u)) {
        u <- as.vector(u)
    }
    if (is.array(u)) {
        u <- as.vector(u)
    }
    if (!is.vector(u)) {
        stop("argument 'u' must be a vector")
    }
    if (is.null(u0)) {
        u0 <- u
    }
    if (is.matrix(u0)) {
        u0 <- as.vector(u0)
    }
    if (is.array(u0)) {
        u0 <- as.vector(u0)
    }
    if (!is.vector(u0)) {
        stop("argument 'u0' must be a vector")
    }
    n1 <- length(u)
    n2 <- length(u0)
    u_rep <- rep(u, each = n2)
    u0_rep <- rep(u0, times = n1)
    u_diff <- matrix(u_rep - u0_rep, nrow = n2)
    return(u_diff)
}
