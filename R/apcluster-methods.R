apcluster.matrix <- function(s, x, p=NA, q=NA, maxits=1000, convits=100,
                             lam=0.9, includeSim=FALSE, details=FALSE,
                             nonoise=FALSE, seed=NA)
{
    if (!is.na(seed)) set.seed(seed)

    apresultObj <- new("APResult") # create the result object to be returned

    apresultObj@call <- deparse(sys.call(-1))

    #
    # check input data
    #
    if (!is.na(p[1]) && (!is.numeric(p) || !is.vector(p)))
        stop("p must be a number or vector")

    if (length(dim(s)) != 2 || ncol(s) != nrow(s))
        stop("s must be a square matrix")

    N <- nrow(s)

    if (length(p) > 1)
    {
        if (length(p) < N)
            stop("vector p is shorter than number of samples")
        else if (length(p) > N)
            p <- p[1:N] # truncate unnecessarily long p
    }

    if (any(is.na(p)) && !is.na(q) && !is.numeric(q))
        stop("q must be a number")

    if (lam > 0.9)
        warning("Large damping factor in use. Turn on details\n",
                "and call plot() to monitor net similarity. The\n",
                "algorithm will change decisions slowly, so consider using\n",
                "a larger value of convits.")

    # If argument p is not given, p is set to median of s
    if (any(is.na(p)))
    {
        if (is.na(q))
            p <- median(s[setdiff(which(s > -Inf), 0:(N - 1) * N + 1:N)])
        else
            p <- quantile(s[setdiff(which(s > -Inf), 0:(N - 1) * N + 1:N)], q)
    }

    apresultObj@l <- N

    # In case user did not remove degeneracies from the input similarities,
    # avoid degenerate solutions by adding a small amount of noise to the
    # input similarities
    if (!nonoise)
    {
        randomMat <- matrix(rnorm(N * N),N)

        s <- s + (.Machine$double.eps * s + .Machine$double.xmin * 100) *
                 randomMat
    }

    attributes(p) <- NULL

    # Place preferences on the diagonal of s (recycled if p is scalar)
    diag(s) <- p

    # store p into result object for future reference
    apresultObj@p <- p

    # replace -Inf (for numerical stability) and NA with -realmax
    infelem <- which(s < -.Machine$double.xmax | is.na(s))

    if (length(infelem) > 0)
        s[infelem] <- -.Machine$double.xmax

    infelem <- which(s > .Machine$double.xmax)

    if (length(infelem) > 0)
        stop("+Inf similarities detected: change to a large positive value,",
             " but smaller than realmax")

    res <- .Call("apclusterC", s, as.integer(maxits),
                 as.integer(convits), as.double(lam), as.logical(details))

    K <- res$K
    I <- res$I[1:K] + 1
    i <- res$it

    if (details)
    {
        apresultObj@idxAll    <- res$idxAll[,1:i] + 1
        apresultObj@netsimAll <- res$netsimAll[1:i]
        apresultObj@dpsimAll  <- res$dpsimAll[1:i]
        apresultObj@exprefAll <- res$exprefAll[1:i]
    }

    if (K > 0)
    {
        i <- i + 1

        c <- max.col(s[, I])

        c[I] <- 1:K # Identify clusters
        c[is.na(c)] <- 0 # R inserts NAs by default, so replace them with 0s
                         # to get the same result as the Matlab code

        # Refine the final set of exemplars and clusters and return results
        for (k in 1:K)
        {
            ii <- which(c == k)
            I[k] <- ii[which.max(colSums(s[ii, ii, drop=FALSE]))]
        }

        names(I) <- colnames(s)[I]

        notI <- matrix(sort(setdiff(1:N, I)), ncol=1)
        c <- max.col(s[, I])
        c[I] <- 1:K

        tmpidx <- I[c]


        tmpdpsim <- sum(s[sub2ind(N, notI, tmpidx[notI])])
        tmpexpref <- sum(diag(s)[I])
        tmpnetsim <- tmpdpsim + tmpexpref

        apresultObj@exemplars <- as.numeric(levels(factor(tmpidx)))

        apresultObj@clusters <- list()

        for (c in 1:length(apresultObj@exemplars))
            apresultObj@clusters[[c]] <- which(tmpidx ==
                                               apresultObj@exemplars[c])

        if (length(colnames(s)) == N)
        {
            names(apresultObj@exemplars) <- colnames(s)[apresultObj@exemplars]

            for (c in 1:length(apresultObj@exemplars))
                 names(apresultObj@clusters[[c]]) <-
                    colnames(s)[apresultObj@clusters[[c]]]
        }
    }
    else
    {
        tmpidx    <- rep(NaN, N)
        tmpnetsim <- NaN
        tmpdpsim <- NaN
        tmpexpref <- NaN

        apresultObj@exemplars <- numeric(0)
        apresultObj@clusters  <- list()
    }

    apresultObj@netsim <- tmpnetsim
    apresultObj@dpsim  <- tmpdpsim
    apresultObj@expref <- tmpexpref
    apresultObj@idx    <- tmpidx
    apresultObj@it     <- i

    if (details)
    {
        apresultObj@netsimAll <- c(apresultObj@netsimAll, tmpnetsim)
        apresultObj@dpsimAll  <- c(apresultObj@dpsimAll, tmpdpsim)
        apresultObj@exprefAll <- c(apresultObj@exprefAll, tmpexpref)
        apresultObj@idxAll    <- cbind(apresultObj@idxAll, tmpidx)
    }

    if (res$unconv)
        warning("Algorithm did not converge. Turn on details\n",
                "and call plot() to monitor net similarity. Consider\n",
                "increasing maxits and convits, and, if oscillations occur\n",
                "also increasing damping factor lam.")

    if (includeSim)
        apresultObj@sim <- s

    apresultObj
}

setMethod("apcluster", signature(s="matrix", x="missing"), apcluster.matrix)


apcluster.function <- function(s, x, p=NA, q=NA, maxits=1000, convits=100,
                               lam=0.9, includeSim=TRUE, details=FALSE,
                               nonoise=FALSE, seed=NA, ...)
{
    if (!is.na(seed)) set.seed(seed)

    if (is.data.frame(x))
        x <- as.matrix(x[, sapply(x, is.numeric)])

    if (is.matrix(x))
        N <- nrow(x)
    else
        N <- length(x)

    if (N < 2) stop("cannot cluster less than 2 samples")

    if (!is.function(s))
    {
        if (!is.character(s) || !exists(s, mode="function"))
            stop("Invalid distance function")

        s <- match.fun(s)
    }

    sim <- s(x=x, ...)

    if (!is.matrix(sim) || (nrow(sim) != N) || ncol(sim) != N)
        stop("Computation of similarity matrix failed")

    apres <- apcluster(s=sim, p=p, q=q, maxits=maxits, convits=convits, lam=lam,
                       details=details, nonoise=nonoise)

    apres@call <- deparse(sys.call(-1))

    if (includeSim)
        apres@sim <- sim

    apres
}

setMethod("apcluster", signature(s="function", x="ANY"), apcluster.function)
setMethod("apcluster", signature(s="character", x="ANY"), apcluster.function)


# Linear index from multiple subscripts.
#   sub2ind is used to determine the equivalent single index
#   corresponding to a given set of subscript values.
sub2ind <- function(N, I, J) (I + (N * (J - 1)))
