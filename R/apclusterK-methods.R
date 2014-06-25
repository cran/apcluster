apclusterK.matrix <- function(s, x, K, prc=10, bimaxit=20, exact=FALSE,
                              maxits=1000, convits=100, lam=0.9,
                              includeSim=FALSE, details=FALSE, nonoise=FALSE,
                              seed=NA, verbose=TRUE)
{
    if (!is.na(seed)) set.seed(seed)

    #
    # check input data
    #
    if (length(dim(s)) != 2 || ncol(s) != nrow(s))
        stop("'s' must be a square matrix")

    N <- nrow(s)

    if (K < 2 || K >= N)
        stop("number of data samples is ", N, ".\n",
             "\tmeaningful range for K: 2 to ", N - 1)

    pminmax <- preferenceRange(s, exact)

    lopref <- pminmax[1]
    hipref <- pminmax[2]
    lok    <- 1
    hik    <- N

    if (is.na(lopref))
        stop("could not find numeric entries in matrix")

    # In case user did not remove degeneracies from the input similarities,
    # avoid degenerate solutions by adding a small amount of noise to the
    # input similarities; we do this here before running apcluster() in order
    # to have deterministic behavior during bisection
    if (!nonoise)
    {
        randomMat <- matrix(rnorm(N*N),N)

        s <- s + (.Machine$double.eps * s + .Machine$double.xmin * 100) *
                 randomMat
    }

    # try to guess better lower bound  before starting with bisection
    ex <- -3
    dn <- FALSE

    while (!dn)
    {
        tmppref <- hipref - 10^ex * (hipref - lopref)

        if (verbose)
            cat("Trying p =", tmppref, "\n")

        apresultObj <- apcluster(s, p=tmppref, nonoise=TRUE)

        tmpk <- length(apresultObj@exemplars)

        if (verbose)
            cat("   Number of clusters:", tmpk, "\n");

        if (tmpk < K)
        {
            lok <- tmpk
            lopref <- tmppref
            dn <- TRUE
        }
        else if (ex == -1)
        {
            dn <- TRUE
        }
        else
        {
            ex <- ex + 1
        }
    }

    # now do bisection (if still necessary)

    ntries <- 0

    while ((abs(tmpk - K) * 100 / K) > prc && ntries < bimaxit)
    {
        ntries <- ntries + 1

        tmppref <- (lopref + hipref) / 2

        if (verbose)
            cat("Trying p =", tmppref, "(bisection step no.", ntries, ")\n")

        apresultObj <- apcluster(s, p=tmppref, nonoise=TRUE,
                                 maxits=maxits, convits=convits, lam=lam,
                                 details=details)

        tmpk <- length(apresultObj@exemplars)

        if (verbose)
            cat("   Number of clusters:", tmpk, "\n");

        if (K > tmpk)
        {
            lopref <- tmppref
            lok    <- tmpk
        }
        else
        {
            hipref <- tmppref
            hik    <- tmpk
        }
    }

    if (verbose)
        cat("\nNumber of clusters:", tmpk, "for p =", tmppref, "\n")

    if ((abs(tmpk - K) * 100 / K) > prc)
        warning("number of clusters not in desired range; Increase 'bimaxit'",
                " to improve accuracy of bisection.")

    apresultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        apresultObj@sim <- s

    apresultObj
}

setMethod("apclusterK", signature(s="matrix", x="missing"), apclusterK.matrix)


apclusterK.function <- function(s, x, K, prc=10, bimaxit=20, exact=FALSE,
                                maxits=1000, convits=100, lam=0.9,
                                includeSim=TRUE, details=FALSE,
                                nonoise=FALSE, seed=NA, verbose=TRUE, ...)
{
    if (!is.na(seed)) set.seed(seed)

    if (is.data.frame(x))
        x <- as.matrix(x[, sapply(x, is.numeric)])

    if (is.matrix(x))
        N <- nrow(x)
    else
        N <- length(x)

    if (N < 2)
        stop("cannot cluster less than 2 samples")

    if (!is.function(s))
    {
        if (!is.character(s) || !exists(s, mode="function"))
            stop("invalid distance function")

        s <- match.fun(s)
    }

    sim <- s(x=x, ...)

    if (!is.matrix(sim) || (nrow(sim) != N) || ncol(sim) != N)
        stop("computation of similarity matrix failed")

    apres <- apclusterK(s=sim, K=K, prc=prc, bimaxit=bimaxit, exact=exact,
                        maxits=maxits, convits=convits, lam=lam,
                        includeSim=FALSE,
                        details=details, nonoise=nonoise, verbose=verbose)

    apres@call <- deparse(sys.call(-1))

    if (includeSim)
        apres@sim <- sim

    apres
}

setMethod("apclusterK", signature(s="function" , x="ANY"), apclusterK.function)
setMethod("apclusterK", signature(s="character", x="ANY"), apclusterK.function)
