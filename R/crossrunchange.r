#' Joint Distribution for Crossings and Runs, Varying Success Probability.
#'
#' @description Joint probability distribution for the number of crossings
#' C and the longest run L in a sequence of n independent Bernoulli observations
#' with p ossibly varying success probability. To enhance precision, results are stored
#' in mpfr arrays and the probabilities are multiplied by \eqn{m^{n-1}}
#' for a multiplier m.

#' @param nmax max sequence length.
#' @param prob success probabilities.
#' @param mult multiplier for joint probabilities.
#' @param prec mpft precision.
#' @param printn logical for progress output.
#' @return list pt of joint probabilities. Cumulative probabilities
#' qt within each row are also included. Further, mostly for code
#' checking, lists pat and qat conditional on starting with a success,
#' and pbt and qbt conditional of starting with a failure, are
#' included.
#' @examples
#' prob10 <- c(rep(.5,5),rep(.7,5))
#' crchange10 <- crossrunchange(nmax=10, prob=prob10,printn=TRUE)
#' print(crchange10$pt[[10]])
#' @export
crossrunchange <- function(nmax = 100, prob = rep(0.5, 100), mult = 2,
    prec = 120, printn = FALSE) {
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
        pat[[nn]][1, nn] <- prod(pmultm[(nmax + 2 - nn):nmax])  # from cond on no crossing
        pbt[[nn]][1, nn] <- prod(qmultm[(nmax + 2 - nn):nmax])  # from cond on no crossing
        for (ff in 2:nn) {
            # from cond on first crossing at ff if last part shortest:
            if (nn - ff + 1 <= ff - 1)
                {
                  f1 <- ff  # unnecessary, but makes code checking easier
                  if (f1 == 2) {
                    prodmulta <- one
                    prodmultb <- one
                  } else {
                    prodmulta <- prod(pmultm[(nmax + 2 - nn):(nmax -
                      nn + f1 - 1)])
                    prodmultb <- prod(qmultm[(nmax + 2 - nn):(nmax -
                      nn + f1 - 1)])
                  }
                  pat[[nn]][2:(nn - f1 + 2), f1 - 1] <- pat[[nn]][2:(nn -
                    f1 + 2), f1 - 1] + prodmulta * qmultm[nmax - nn +
                    f1] * qbt[[nn - f1 + 1]][1:(nn - f1 + 1), nn - f1 +
                    1]
                  pbt[[nn]][2:(nn - f1 + 2), f1 - 1] <- pbt[[nn]][2:(nn -
                    f1 + 2), f1 - 1] + prodmultb * pmultm[nmax - nn +
                    f1] * qat[[nn - f1 + 1]][1:(nn - f1 + 1), nn - f1 +
                    1]
                }  # end if last part shortest
            if (nn - ff + 1 > ff - 1)
                {
                  # if last part longest
                  f2 <- ff  # unnecessary, but makes code checking easier
                  if (f2 == 2) {
                    prodmulta <- one
                    prodmultb <- one
                  } else {
                    prodmulta <- prod(pmultm[(nmax + 2 - nn):(nmax -
                      nn + f2 - 1)])
                    prodmultb <- prod(qmultm[(nmax + 2 - nn):(nmax -
                      nn + f2 - 1)])
                  }
                  pat[[nn]][2:(nn - f2 + 2), f2 - 1] <- pat[[nn]][2:(nn -
                    f2 + 2), f2 - 1] + prodmulta * qmultm[nmax - nn +
                    f2] * qbt[[nn - f2 + 1]][1:(nn - f2 + 1), f2 - 1]
                  pat[[nn]][2:(nn - f2 + 2), f2:(nn - f2 + 1)] <- pat[[nn]][2:(nn -
                    f2 + 2), f2:(nn - f2 + 1)] + prodmulta * qmultm[nmax -
                    nn + f2] * pbt[[nn - f2 + 1]][1:(nn - f2 + 1), f2:(nn -
                    f2 + 1)]
                  pbt[[nn]][2:(nn - f2 + 2), f2 - 1] <- pbt[[nn]][2:(nn -
                    f2 + 2), f2 - 1] + prodmultb * pmultm[nmax - nn +
                    f2] * qat[[nn - f2 + 1]][1:(nn - f2 + 1), f2 - 1]
                  pbt[[nn]][2:(nn - f2 + 2), f2:(nn - f2 + 1)] <- pbt[[nn]][2:(nn -
                    f2 + 2), f2:(nn - f2 + 1)] + prodmultb * pmultm[nmax -
                    nn + f2] * pat[[nn - f2 + 1]][1:(nn - f2 + 1), f2:(nn -
                    f2 + 1)]
                }  # end if last part longest
        }  # end for ff
        pt[[nn]] <- pm[nmax - nn + 1] * pat[[nn]] + qm[nmax - nn + 1] *
            pbt[[nn]]
        qat[[nn]] <- cumsumm(pat[[nn]])
        qbt[[nn]] <- cumsumm(pbt[[nn]])
        qt[[nn]] <- pm[nmax - nn + 1] * qat[[nn]] + qm[nmax - nn + 1] *
            qbt[[nn]]
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


