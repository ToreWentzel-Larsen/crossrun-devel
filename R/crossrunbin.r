#' Joint Distribution for Crossings and Runs
#'
#' @description Joint probability distribution for the number of crossings
#' C and the longest run L in a sequence of n independent Bernoulli observations
#' with success probability p. To enhance precision, results are stored
#' in mpfr arrays and the probabilities are multiplied by \eqn{m^{n-1}}
#' for a multiplier m.
#' @param nmax max sequence length.
#' @param prob success probability.
#' @param mult multiplier for joint probabilities.
#' @param prec mpft precision.
#' @param printn logical for progress output.
#' @return list of joint probabilities.
#' @examples
#' crb10.6 <- crossrunbin(nmax=10, prob=.6, printn=TRUE)
#' print(crb10.6$pt[[10]])
#' @export
crossrunbin <- function(nmax = 100, prob = 0.5, mult = 2, prec = 120,
    printn = FALSE) {
    nill <- Rmpfr::mpfr(0, prec)
    one <- Rmpfr::mpfr(1, prec)
    multm <- Rmpfr::mpfr(mult, prec)
    pm <- Rmpfr::mpfr(prob, prec)
    qm <- one - pm
    pmultm <- pm * multm
    qmultm <- qm * multm
    # conditioning of S= first value, pat: above 0, pbt: below 0 suffix
    # t: probabilities times multm^(n-1).  n=1:
    pat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    pbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    pt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    qat <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    qbt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    qt <- list(pt1 = Rmpfr::mpfr2array(one, dim = c(1, 1)))
    for (nn in 2:nmax) {
        pat[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, nn))
        pbt[[nn]] <- Rmpfr::mpfr2array(rep(nill, nn * nn), dim = c(nn, nn))
        rownames(pat[[nn]]) <- c(0:(nn - 1))
        rownames(pbt[[nn]]) <- c(0:(nn - 1))
        colnames(pat[[nn]]) <- c(1:nn)
        colnames(pbt[[nn]]) <- c(1:nn)
        pat[[nn]][1, nn] <- (pmultm^(nn - 1))  # from cond on no crossing
        pbt[[nn]][1, nn] <- (qmultm^(nn - 1))  # from cond on no crossing
        for (ff in 2:nn) {
            # from cond on first crossing at ff if last part shortest:
            if (nn - ff + 1 <= ff - 1)
                {
                  f1 <- ff  # unnecessary, but makes code checking easier
                  pat[[nn]][2:(nn - f1 + 2), f1 - 1] <- pat[[nn]][2:(nn -
                    f1 + 2), f1 - 1] + (pmultm^(f1 - 2)) * qmultm *
                    qbt[[nn - f1 + 1]][1:(nn - f1 + 1), nn - f1 + 1]
                  pbt[[nn]][2:(nn - f1 + 2), f1 - 1] <- pbt[[nn]][2:(nn -
                    f1 + 2), f1 - 1] + (qmultm^(f1 - 2)) * pmultm *
                    qat[[nn - f1 + 1]][1:(nn - f1 + 1), nn - f1 + 1]
                }  # end if last part shortest
            if (nn - ff + 1 > ff - 1)
                {
                  # if last part longest
                  f2 <- ff  # unnecessary, but makes code checking easier
                  pat[[nn]][2:(nn - f2 + 2), f2 - 1] <- pat[[nn]][2:(nn -
                    f2 + 2), f2 - 1] + (pmultm^(f2 - 2)) * qmultm *
                    qbt[[nn - f2 + 1]][1:(nn - f2 + 1), f2 - 1]
                  pat[[nn]][2:(nn - f2 + 2), f2:(nn - f2 + 1)] <- pat[[nn]][2:(nn -
                    f2 + 2), f2:(nn - f2 + 1)] + (pmultm^(f2 - 2)) *
                    qmultm * pbt[[nn - f2 + 1]][1:(nn - f2 + 1), f2:(nn -
                    f2 + 1)]
                  pbt[[nn]][2:(nn - f2 + 2), f2 - 1] <- pbt[[nn]][2:(nn -
                    f2 + 2), f2 - 1] + (qmultm^(f2 - 2)) * pmultm *
                    qat[[nn - f2 + 1]][1:(nn - f2 + 1), f2 - 1]
                  pbt[[nn]][2:(nn - f2 + 2), f2:(nn - f2 + 1)] <- pbt[[nn]][2:(nn -
                    f2 + 2), f2:(nn - f2 + 1)] + (qmultm^(f2 - 2)) *
                    pmultm * pat[[nn - f2 + 1]][1:(nn - f2 + 1), f2:(nn -
                    f2 + 1)]
                }  # end if last part longest
        }  # end for ff
        pt[[nn]] <- pm * pat[[nn]] + qm * pbt[[nn]]
        qat[[nn]] <- cumsumm(pat[[nn]])
        qbt[[nn]] <- cumsumm(pbt[[nn]])
        qt[[nn]] <- pm * qat[[nn]] + qm * qbt[[nn]]
        rownames(pt[[nn]]) <- c(0:(nn - 1))
        colnames(pt[[nn]]) <- c(1:nn)
        rownames(qat[[nn]]) <- c(0:(nn - 1))
        colnames(qat[[nn]]) <- c(1:nn)
        rownames(qbt[[nn]]) <- c(0:(nn - 1))
        rownames(qat[[nn]]) <- c(0:(nn - 1))
        colnames(qt[[nn]]) <- c(1:nn)
        colnames(qt[[nn]]) <- c(1:nn)
        if (printn)
            {
                print(nn)
                print(Sys.time())
            }  # end optional timing information
    }  # end for nn
    names(pat) <- paste("pat", 1:nmax, sep = "")
    names(pbt) <- paste("pbt", 1:nmax, sep = "")
    names(pt) <- paste("pt", 1:nmax, sep = "")
    names(qat) <- paste("qat", 1:nmax, sep = "")
    names(qbt) <- paste("qbt", 1:nmax, sep = "")
    names(qt) <- paste("qt", 1:nmax, sep = "")
    return(list(pat = pat, pbt = pbt, pt = pt, qat = qat, qbt = qbt,
        qt = qt))
}

#' wrapper for crossrunbin, succes probability=pnorm(shift).
#'
#' @param nmax max sequence length.
#' @param shift mean of normal distribution.
#' @param mult multiplier for joint probabilities.
#' @param prec mpft precision.
#' @param printn logical for progress output.
#' @return list pt of joint probabilities. Cumulative probabilities
#' qt within each row are also included. Further, mostly for code
#' checking, lists pat and qat conditional on starting with a success,
#' and pbt and qbt conditional of starting with a failure, are
#' included.
#' @examples
#' crs14 <- crossrunshift(nmax=14,printn=TRUE)
#' print(crs14$pt[[14]])
#' @export
crossrunshift <- function(nmax = 100, shift = 0, mult = 2, prec = 120,
    printn = FALSE) {
    prob <- stats::pnorm(shift)
    return(crossrunbin(nmax = nmax, prob = prob, mult = mult, prec = prec,
        printn = printn))
}

