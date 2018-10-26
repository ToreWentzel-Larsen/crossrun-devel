#' Joint Probabilities for Crossings and Runs, Symmetric Case
#' @description Joint probability distribution for the number of crossings
#' C and the longest run L in a sequence of n independent Bernoulli observations
#' with success probability p. To enhance precision, results are stored
#' in mpfr arrays and the probabilities are multiplied by \eqn{m^{n-1}}
#' for a multiplier m. This is for the symmetric case with success
#' probability 0.5, in which the multiplied probabilities are
#' integers for the default value 2 of the multiplier.
#'
#' @param nmax ; max sequence length.
#' @param mult ; multiplier for joint probabilities. Default 2.
#' @param prec ; mpft precision.
#' @param printn ; logical for including progress output.
#' @return pt, list of joint probabilities, multiplied with \eqn{m^{n-1}}.
#' In addition cumulative probabilities qt within each row are also included.
#' @examples
#' crs10 <- crossrunsymm(nmax=10,printn=TRUE)
#' @export
crossrunsymm <- function(nmax = 100, mult = 2, prec = 120,
    printn = FALSE) {
    nill <- Rmpfr::mpfr(0, prec)
    one <- Rmpfr::mpfr(1, prec)
    two <- Rmpfr::mpfr(2, prec)
    multm <- Rmpfr::mpfr(mult, prec)
    pmultm <- multm/two
    pt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    qt <- pt
    for (nn in 2:nmax) {
        pt[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, nn))
        rownames(pt[[nn]]) <- c(0:(nn - 1))
        colnames(pt[[nn]]) <- c(1:nn)
        pt[[nn]][1, nn] <- pmultm^(nn - 1)
        for (ff in 2:nn) {
            if (nn - ff + 1 <= ff - 1)
                {
                  f1 <- ff
                  pt[[nn]][(1 + 1):(nn - f1 + 2), f1 - 1] <- pt[[nn]][2:(nn -
                    f1 + 2), f1 - 1] + (pmultm^(f1 - 2)) * qt[[nn -
                    f1 + 1]][1:(nn - f1 + 1), nn - f1 + 1]
                }  # end if last part shortest
            if (nn - ff + 1 > ff - 1)
                {
                  f2 <- ff
                  pt[[nn]][2:(nn - f2 + 2), f2 - 1] <- pt[[nn]][2:(nn -
                    f2 + 2), f2 - 1] + (pmultm^(f1 - 2)) * qt[[nn -
                    f2 + 1]][1:(nn - f2 + 1), f2 - 1]
                  pt[[nn]][2:(nn - f2 + 2), f2:(nn - f2 + 1)] <- pt[[nn]][2:(nn -
                    f2 + 2), f2:(nn - f2 + 1)] + (pmultm^(f1 - 2)) *
                    pt[[nn - f2 + 1]][1:(nn - f2 + 1), f2:(nn - f2 +
                      1)]
                }  # end if last part longest
        }  # end for ff
        qt[[nn]] <- cumsumm(pt[[nn]])
        rownames(qt[[nn]]) <- c(0:(nn - 1))
        colnames(qt[[nn]]) <- c(1:nn)
        if (printn)
            {
                print(nn)
                print(Sys.time())
            }  # end optional timing information
    }  # end for nn
    names(pt) <- paste("pt", 1:nmax, sep = "")
    names(qt) <- paste("qt", 1:nmax, sep = "")
    return(list(pt = pt, qt = qt))
}  # end function crossrunsymm
