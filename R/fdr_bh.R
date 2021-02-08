#' @title fdr_bh
#' @description Executes the \insertCite{benjaminiControllingFalseDiscovery1995;textual}{llfdr}
#'           and the \insertCite{benjaminiControlFalseDiscovery2001;textual}{llfdr} procedure for controlling the false discovery
#'             rate (FDR) of a family of hypothesis tests. FDR is the expected
#'             proportion of rejected hypotheses that are mistakenly rejected
#'             (i.e., the null hypothesis is actually true for those tests).
#'             FDR is a somewhat less conservative/more powerful method for
#'             correcting for multiple comparisons than procedures like Bonferroni
#'             correction that provide strong control of the family-wise
#'             error rate (i.e., the probability that one or more null
#'              hypotheses are mistakenly rejected).
#' @param pvals A vector  containing the p-value of each individual test
#'              in a family of tests.
#' @param q The desired false discovery rate. The default is 0.05.
#' @param method 'pdep' or 'dep'. If 'pdep,' the original Bejnamini & Hochberg
#'               FDR procedure is used, which is guaranteed to be accurate if
#'               the individual tests are independent or positively dependent
#'               (e.g., Gaussian variables that are positively correlated or
#'               independent).  If 'dep,' the FDR procedure described in
#'               Benjamini & Yekutieli (2001) that is guaranteed to be accurate
#'               for any test dependency structure (e.g., Gaussian variables
#'               with any covariance matrix) is used. 'dep' is always
#'               appropriate to use but is less powerful than 'pdep.'
#'               The default is 'pdep'.
#'
#' @return A binary vector of the same size as the input "pvals." If the
#'         ith element of h is 1, then the test that produced the ith p-value
#'         in pvals is significant
#'         (i.e., the null hypothesisof the test is rejected).
#' @export
#' @references
#' \insertRef{benjaminiControllingFalseDiscovery1995}{llfdr}
#'
#' \insertRef{benjaminiControlFalseDiscovery2001}{llfdr}
#'
fdr_bh <- function(pvals, q = 0.05, method = 'pdep') {
    # sort_result <- sort(pvals, index.return = TRUE)
    # p_sorted <- sort_result$x
    # sort_ids <- sort_result$ix
    # sort_result1 <- sort(sort_ids, index.return = TRUE)
    # unsort_ids <- sort_result1$ix
    p_sorted <- sort(pvals)
    m <- length(p_sorted)
    seq_m <- 1:m
    if (method == 'pdep') {
        # BH procedure for independence or positive dependence
        thresh <- seq_m * q / m
        wtd_p <- m * p_sorted / seq_m
    } else if(method == 'dep'){
        # BH procedure for any dependency structure
        denom <- m * sum(1 / seq_m)
        thresh <- seq_m * q / denom
        wtd_p = denom * p_sorted / seq_m
    } else {
        stop('Argument method needs to be "pdep" or "dep".')
    }
    rej <- p_sorted <= thresh
    max_id <- max(which(rej, TRUE)) # find greatest significant pvalue
    if (is.infinite(max_id)) {
        crit_p <- 0
        h <- pvals * 0
    } else {
        crit_p <- p_sorted[max_id];
        h <- pvals <= crit_p
    }
    return(h)
}





