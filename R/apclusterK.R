apclusterK <- function(s, K, prc=10, bimaxit=20, exact=FALSE,
                       nonoise=FALSE, seed=NA, ...)
{
    if (!is.na(seed)) set.seed(seed)

    #
    # check input data
    #

    if (length(dim(s)) != 2 || ncol(s) != nrow(s))
    {
        stop("s must be a square matrix")
    }

    N <- nrow(s)

    if (K < 2 || K >= N)
    {
        stop("Number of data samples is ", N, ".\n",
             "       Meaningful range for K: 2 to ", N - 1, "\n", call.=FALSE)
    }

    pminmax <- preferenceRange(s, exact)

    lopref <- pminmax[1]
    hipref <- pminmax[2]
    lok    <- 1
    hik    <- N

    if (is.na(lopref)) stop("Could not find numeric entries in matrix\n")

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

        cat("Trying p =", tmppref, "\n")

        apresultObj <- apcluster(s, p=tmppref, nonoise=TRUE, ...)

        tmpk <- length(apresultObj@exemplars)

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
        
        message("Trying p = ", tmppref, " (bisection step no. ", ntries, ")\n")

        apresultObj <- apcluster(s, p=tmppref, nonoise=TRUE, ...)

        tmpk <- length(apresultObj@exemplars)

        message("   Number of clusters:", tmpk, "\n");

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

    message("\nNumber of clusters: ", tmpk, " for p = ", tmppref, "\n")

    if ((abs(tmpk - K) * 100 / K) > prc)
    {
        warning("Number of clusters not in desired range. Increase bimaxit",
                " to improve accuracy of bisection.")
    }

    apresultObj
}
