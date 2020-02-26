#' @title Cubic
#'
#' @param a the
#' @param b the
#'
#' @return the real solution
#' @export

Cubic <- function(a, b) {
    p <- 3 * (a - 1) - b^2
    q <- 36 * b - 9 * a * b + 2 * b^3
    b3 <- b / 3
    sqrt_abs_p <- sqrt(abs(p))
    C <- q / sqrt_abs_p^(3) / 2
    x <- ifelse(p == 0, q^(1/3) / 3 + b3, 0)

    Ind1 <- p > 0
    if(any(Ind1)) {
        x[Ind1] <- 2/3 * sqrt_abs_p[Ind1] * sinh(1/3 * asinh(C[Ind1])) + b3[Ind1]
    }
    Ind2 <- p < 0 & C >= 1
    if(any(Ind2)) {
        x[Ind2] <- 2/3 * sqrt_abs_p[Ind2] * cosh(1/3 * acosh(C[Ind2])) + b3[Ind2]
    }
    Ind3 <- p < 0 & C <= -1
    if(any(Ind3)) {
        x[Ind3] <- -2/3 * sqrt_abs_p[Ind3] * cosh(1/3 * acosh(-C[Ind3])) + b3[Ind3]
    }
    Ind4 <- p < 0 & abs(C) < 1
    if(any(Ind4)){
        sqrt_abs_ind4 <- sqrt_abs_p[Ind4]
        acosC_ind4 <- 1/3 * acos(C[Ind4])
        b3_ind4 <- b3[Ind4]
        xreal1 <- 2/3 * sqrt_abs_ind4 * cos(acosC_ind4) + b3_ind4
        xreal2 <- 2/3 * sqrt_abs_ind4 * cos(2/3 * pi + acosC_ind4) + b3_ind4
        xreal3 <- 2/3 * sqrt_abs_ind4 * cos(4/3 * pi + acosC_ind4) + b3_ind4
        xreal <- rbind(xreal1, xreal2, xreal3)
        Ind_real <- abs(xreal) < 1
        SInd1 <- SInd2 <- SInd3  <- rep(FALSE, length(a))
        SInd <- colSums(Ind_real)
        SInd1temp <- SInd == 1
        SInd2temp <- SInd == 2
        SInd3temp <- SInd == 3
        if(any(SInd1temp)) {
            SInd1[Ind4] <- SInd1temp
            SInd1temp_mat <- matrix(rep(SInd1temp, each = 3), nrow = 3)
            x[SInd1] <- xreal[SInd1temp_mat & Ind_real]
        }
        if(any(SInd2temp)) {
            SInd2[Ind4] <- SInd2temp
            SInd2temp_mat <- matrix(rep(SInd2temp, each = 3), nrow = 3)
            xx <- matrix(xreal[SInd2temp_mat & Ind_real], nrow = 2)
            a_mat <- matrix(rep(a[SInd2], each = 2), nrow = 2)
            b_mat <- matrix(rep(b[SInd2], each = 2), nrow = 2)
            xx_mat <- (1 - xx^2)
            LH <- -a_mat / xx_mat + 2 * xx * b_mat / xx_mat - log(xx_mat)
            maxLH <- apply(LH, 2, max)
            x[SInd2] <- xx[match(maxLH, LH)]
        }
        if(any(SInd3temp)) {
            SInd3[Ind4] <- SInd3temp
            SInd3temp_mat <- matrix(rep(SInd3temp, each = 3), nrow = 3)
            xx <- matrix(xreal[SInd3temp_mat & Ind_real], nrow = 3)
            a_mat <- matrix(rep(a[SInd3], each = 3), nrow = 3)
            b_mat <- matrix(rep(b[SInd3], each = 3), nrow = 3)
            xx_mat <- (1 - xx^2)
            LH <- -a_mat / xx_mat + 2 * xx * b_mat / xx_mat - log(xx_mat)
            maxLH <- apply(LH, 2, max)
            x[SInd3] <- xx[match(maxLH, LH)]
        }
    }
    return(x)
}

