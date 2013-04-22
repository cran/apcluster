apclusterL.matrix <- function(s, x, sel, p=NA, q=NA, maxits=1000, convits=100,
                              lam=0.9, includeSim=FALSE, nonoise=FALSE, seed=NA)
{
    if (!is.na(seed)) set.seed(seed)

    #
    # check input data
    #
    if (length(dim(s)) != 2)
        stop("s must be a matrix")

    M <- ncol(s)
    N <- nrow(s)

    if (M > N)
        stop("no. of columns of s may not be larger than number of rows")

    if (!is.vector(sel) || !is.numeric(sel))
        stop("sel must be a numeric vector")

    if (length(sel) != M)
        stop("vector sel is shorter or longer than number of selected samples")

    if (max(sel) > N || min(sel) < 1)
        stop("sample index in sel must be between one and number of samples")

    if (!is.na(p) && (!is.numeric(p) || !is.vector(p)))
        stop("p must be a number or vector")

    if (length(p) > 1)
    {
        if (length(p) < N)
            stop("vector p is shorter than number of samples")
        else if (length(p) > N)
            p <- p[1:N] # truncate unnecessarily long p
    }

    if (is.na(p) && !is.na(q) && !is.numeric(q))
        stop("q must be a number")

    if (lam > 0.9)
        warning("Large damping factor in use. The algorithm\n",
                "will change decisions slowly, so consider using\n",
                "a larger value of convits.")

    # If argument p is not given, p is set to median of s
    if (is.na(p))
    {
        if (is.na(q))
            p <- median(s[setdiff(which(s > -Inf), (1:M-1) * N + sel)])
        else
            p <- quantile(s[setdiff(which(s > -Inf), (1:M-1) * N + sel)], q)
    }

    attributes(p) <- NULL

    apresultObj <- new("APResult") # create the result object to be returned

    apresultObj@call <- deparse(sys.call(-1))

    # store p into result object for future reference
    apresultObj@p <- p

    if (length(p) == 1)
        p <- rep(p, N)

    apresultObj@l <- N
    apresultObj@sel <- sel

    # In case user did not remove degeneracies from the input similarities,
    # avoid degenerate solutions by adding a small amount of noise to the
    # input similarities
    if (!nonoise)
    {
        randomMat <- matrix(rnorm(N*M),N,M)

        s <- s + (.Machine$double.eps * s + .Machine$double.xmin * 100) *
                 randomMat
    }

    # Append preferences as additional column to s
    s <- cbind(s, p)

    # Numerical stability -- replace -Inf with -realmax
    infelem <- which(s < -.Machine$double.xmax)

    if (length(infelem) > 0)
    {
        warning("-Inf similarities detected: changing to -realmax ",
                " to ensure numerical stability")

        s[infelem] <- -.Machine$double.xmax
    }

    infelem <- which(s > .Machine$double.xmax)

    if (length(infelem) > 0)
        stop("+Inf similarities detected: change to a large positive value,",
             " but smaller than realmax")

    # convert to C sample indices
    selC <- sel - 1

    res <- .Call("apclusterLeveragedC", s, selC, as.integer(maxits),
                 as.integer(convits), as.double(lam))

    K <- res$K

    # convert cluster center indices to R
    I <- res$I[1:K] + 1
    i <- res$it

    if (K > 0)
    {
        i <- i + 1

        ee <- which(sel %in% I)
        if (length(ee) < 1)
            stop("Internal error: no exemplars in selected samples")

        # select exemplar with the highest similarity
        # diagonal elements are not relevant, they are
        # overwritten with 1:K
        c <- rep(NA, N)
        for (j in 1:N)
        {
            if (j %in% sel[ee])
            {
                sc <- s[j,ee]
                sc[which(j %in% sel[ee])] <- s[j, M + 1]
                c[j] <- which(I %in% sel[ee[which.max(sc)]])
            }
            else
            {
                c[j] <- which(I %in% sel[ee[which.max(s[j, ee])]])
            }
        }

        c[I] <- 1:K # Identify clusters
        c[is.na(c)] <- 0 # R inserts NAs by default, so replace them with 0s
                         # to get the same result as the Matlab code

        # Refine the final set of exemplars and clusters and return results
        for (k in 1:K)
        {
            jj <- which(c == k)
            ii <- which(sel %in% jj)
            ns <- s[jj, M + 1]

            if (length(ii) > 0)
            {
                ns[match(sel[ii],jj)] <- colSums(s[jj, ii, drop=FALSE]) +
                                         s[sel[ii], M + 1] -
                                             diag(s[sel[ii], ii, drop=FALSE])
                I[k] <- jj[match(sel[ii],jj)[which.max(ns[match(sel[ii], jj)])]]
            }
            else
            {
                if (length(jj) > 1)
                    stop("Internal error: wrong number of exemplar candidates")
            }
        }

        names(I) <- rownames(s)[I]
        notI <- matrix(sort(setdiff(1:N, I)), ncol=1)

        ee <- which(sel %in% I)

        for (j in 1:N)
        {
            if (j %in% sel[ee])
            {
                sc <- s[j, ee]
                sc[which(j %in% sel[ee])] <- s[j, M + 1]
                c[j] <- which(I %in% sel[ee[which.max(sc)]])
            }
            else
                 c[j] <- which(I %in% sel[ ee[which.max(s[j, ee])]])
        }

        c[I] <- 1:K
        tmpidx <- I[c]

        # Self similarities not relevant
        tmpdpsim <- sum(s[sub2ind(N, notI, match(tmpidx[notI], sel))])
        tmpexpref <- sum(s[I, M + 1])
        tmpnetsim <- tmpdpsim + tmpexpref

        apresultObj@exemplars <- as.numeric(levels(factor(tmpidx)))

        apresultObj@clusters <- list()

        for (c in 1:length(apresultObj@exemplars))
            apresultObj@clusters[[c]] <- which(tmpidx ==
                                               apresultObj@exemplars[c])

        if (length(rownames(s)) == N)
        {
            names(apresultObj@exemplars) <- rownames(s)[apresultObj@exemplars]

            for (c in 1:length(apresultObj@exemplars))
                names(apresultObj@clusters[[c]]) <-
                    rownames(s)[apresultObj@clusters[[c]]]
        }
    }
    else
    {
        tmpidx    <- rep(NaN, N)
        tmpnetsim <- NaN
        tmpdpsim  <- NaN
        tmpexpref <- NaN

        apresultObj@exemplars <- numeric(0)
        apresultObj@clusters  <- list()
    }

    apresultObj@netsim <- tmpnetsim
    apresultObj@dpsim  <- tmpdpsim
    apresultObj@expref <- tmpexpref
    apresultObj@idx    <- tmpidx
    apresultObj@it     <- i

    if (res$unconv)
        warning("Algorithm did not converge. Turn on details\n",
                "and call plot() to monitor net similarity. Consider\n",
                "increasing maxits and convits, and, if oscillations occur\n",
                "also increasing damping factor lam.")

    apresultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        apresultObj@sim <- s

    apresultObj
}

setMethod("apclusterL", signature(s="matrix", x="missing"), apclusterL.matrix)


apclusterL.function <- function(s, x, frac, sweeps, p=NA, q=NA,
                                maxits=1000, convits=100, lam=0.9,
                                includeSim=TRUE, nonoise=FALSE, seed=NA, ...)
{
    if (frac <=0 || frac > 1)
        stop("Invalid fraction of samples specified")

    if (!is.na(seed)) set.seed(seed)

    if (is.data.frame(x))
        x <- as.matrix(x[, sapply(x, is.numeric)])

    if (is.matrix(x))
        N <- nrow(x)
    else
        N <- length(x)

    if (N < 2) stop("cannot cluster less than 2 samples")

    nsel <- max(ceiling(N * frac), 2)
    sel <- sort(sample(1:N, nsel))

    if (!is.function(s))
    {
        if (!is.character(s) || !exists(s, mode="function"))
            stop("Invalid distance function")

        s <- match.fun(s)
    }

    apresultObj <- new("APResult") # create the result object to be returned
    apresultObj@netsim <- -Inf
    netsimL <- rep(-Inf, sweeps)

    for (i in 1:sweeps)
    {
        sim <- s(x=x, sel=sel, ...)

        if (!is.matrix(sim) || nrow(sim) != N || ncol(sim) != length(sel))
            stop("Computation of similarity matrix failed")

        apres <- apclusterL(s=sim, sel=sel, p=p, q=q, maxits=maxits,
                            convits=convits, lam=lam, nonoise=nonoise)

        netsimL[i] <- apres@netsim

        if (apres@netsim > apresultObj@netsim || apresultObj@netsim == -Inf)
        {
            apresultObj <- apres

            if (includeSim)
                apresultObj@sim <- sim
            else
                apresultObj@sim <- matrix(NA, 1, 1)
        }

        sel <- sort(unique(apresultObj@idx))

        if (nsel - length(sel) > 0)
        {
            otherSamples <- setdiff(1:N, sel)
            sel <- sort(c(sel, sample(otherSamples, nsel - length(sel))))

            if (length(rownames(sim)) > 0)
                names(sel) <- rownames(sim)[sel]
        }
        else
            break
    }

    apresultObj@call <- deparse(sys.call(-1))
    apresultObj@sweeps <- sweeps
    apresultObj@netsimLev <- netsimL
    apresultObj
}

setMethod("apclusterL", signature(s="function" , x="ANY"), apclusterL.function)
setMethod("apclusterL", signature(s="character", x="ANY"), apclusterL.function)
