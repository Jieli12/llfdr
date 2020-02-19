#' @title interpolate
#' @description this function inserts the positive value into the entries where
#' the value is negative. If all of the vector are negative, then replaces each
#' entry with 1e-4.
#'
#' @param x a numeric vector
#'
#' @return a vector which each entry is non-negative
#' @export
#'
#' @examples
interpolate <- function(x) {
    ind <- x < 0
    if (all(ind)) {
        x[ind] <- 1e-4
    } else if (any(ind)) {
        x[ind] <- mean(x[!ind])
    }
    return(x)
}
