apclusterLM <- function(s, p=NA, q=NA, maxits=1000, convits=100,
                      lam=0.9, details=FALSE, nonoise=FALSE, seed=NA)
{
    if (!is.na(seed)) set.seed(seed)

    apresultObj <- new("APResult") # create the result object to be returned

    #
    # check input data
    #

    if (!is.na(p) && (!is.numeric(p) || !is.vector(p)))
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

    if (lam > 0.9)
        warning("Large damping factor in use. Turn on details\n",
                "and call plot() to monitor net similarity. The\n",
                "algorithm will change decisions slowly, so consider using\n",
                "a larger value of convits.")

    # If argument p is not given, p is set to median of s
    if (is.na(p))
    {
        if (is.na(q))
        {
            p <- median(s[setdiff(which(s > -Inf), 0:(N-1) * N + 1:N)])
        }
        else
        {
            p <- quantile(s[setdiff(which(s > -Inf), 0:(N-1) * N + 1:N)], q)
        }
    }

    apresultObj@l <- N

    # In case user did not remove degeneracies from the input similarities,
    # avoid degenerate solutions by adding a small amount of noise to the
    # input similarities
    if (!nonoise)
    {
        randomMat <- matrix(rnorm(N*N),N)

        s <- s + (.Machine$double.eps * s + .Machine$double.xmin * 100) *
                 randomMat
    }

    attributes(p) <- NULL

    # Place preferences on the diagonal of s (recycled if p is scalar)
    diag(s) <- p

    # store p into result object for future reference
    apresultObj@p <- p

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

    # create temporary storage
    e   <- matrix(0, N, convits)
    dn  <- FALSE
    i   <- 0
    dS  <- diag(s)
    A   <- matrix(0, N, N)
    R   <- matrix(0, N, N)
    t   <- 1

    if (details)
    {
        apresultObj@idxAll    <- matrix(0, N, maxits + 1)
        apresultObj@netsimAll <- seq(0, 0, length.out=(maxits + 1))
        apresultObj@dpsimAll  <- seq(0, 0, length.out=(maxits + 1))
        apresultObj@exprefAll <- seq(0, 0, length.out=(maxits + 1))
    }

    # main loop

    while (!dn)
    {
        i <- i + 1;

        # Compute responsibilities

        for (ii in 1:N)
        {
            old <- R[ii,]

            AS <- A[ii,] + s[ii,]
            I <- which.max(AS) # get greatest element
            Y <- AS[I]
            AS[I] <- -Inf

            Y2 <- max(AS); # get second-greatest element

            R[ii,] <- s[ii,] - Y # make update
            R[ii,I] <- s[ii, I] - Y2

            R[ii,] <- (1 - lam) * R[ii,] + lam * old # damping

            # truncate too large values
            R[ii,R[ii,] > .Machine$double.xmax] <- .Machine$double.xmax
        }

        # Compute availabilities

        for (jj in 1:N)
        {
            old <- A[,jj]
            Rp <- (R[,jj] + abs(R[,jj])) / 2
            Rp[jj] <- R[jj,jj]

            auxsum <- sum(Rp)
            A[,jj] <- auxsum - Rp
            dA <- A[jj,jj]
            A[,jj] <- (A[,jj] - abs(A[,jj])) / 2
            A[jj,jj] <- dA
            A[,jj] <- (1-lam) * A[,jj] + lam * old # Damping
        }

        # determine clusters and check for convergence
        E <- as.numeric((diag(A) + diag(R)) > 0)
        e[,((i - 1) %% convits) + 1] <- E
        K <- sum(E)

        if (i >= convits || i >= maxits)
        {
            se <- rowSums(e)

            unconverged <- (sum((se == convits) + (se == 0)) != N)

            if ((!unconverged && (K > 0)) || (i == maxits))
            {
                dn <- TRUE
            }
        }

        if (K==0)
        {
            tmpnetsim <- NaN
            tmpdpsim <- NaN
            tmpexpref <- NaN
            tmpidx <- NaN
        }
        else
        {
            I <- which(E != 0)
            notI <- which(E == 0)

            c         <- max.col(s[,I])
            c[I]      <- 1:K
            tmpidx    <- I[c]
            tmpdpsim  <- sum(s[sub2ind(N, notI, tmpidx[notI])])
            tmpexpref <- sum(dS[I])
            tmpnetsim <- tmpdpsim + tmpexpref
        }

        if (details)
        {
            apresultObj@netsimAll[i] <- tmpnetsim
            apresultObj@dpsimAll[i]  <- tmpdpsim
            apresultObj@exprefAll[i] <- tmpexpref
            apresultObj@idxAll[,i]   <- tmpidx
        }
    } # iterations

    I <- which((diag(A) + diag(R)) > 0)
    K <- length(I) # Identify exemplars

    if (K > 0)
    {
        c <- max.col(s[,I])

        c[I] <- 1:K # Identify clusters
        c[is.na(c)] <- 0 # R inserts NAs by default, so replace them with 0s
                         # to get same result as the Matlab code

        # Refine the final set of exemplars and clusters and return results
        for (k in 1:K)
        {
            ii <- which(c == k)
            j <- which.max(colSums(s[ii,ii,drop=FALSE]))
            I[k] <- ii[j[1]]
        }

        notI <- matrix(sort(setdiff(1:N, I)), ncol=1)
        c <- max.col(s[,I])
        c[I] <- 1:K
        tmpidx <- I[c]
        tmpdpsim <- sum(s[sub2ind(N, notI, tmpidx[notI])])
        tmpexpref <- sum(dS[I])
        tmpnetsim <- tmpdpsim + tmpexpref

        apresultObj@exemplars <- as.numeric(levels(factor(tmpidx)))

        apresultObj@clusters <- list()

        for (c in 1:length(apresultObj@exemplars))
        {
            apresultObj@clusters[[c]] <- which(tmpidx ==
                                               apresultObj@exemplars[c])
        }

        if (length(colnames(s)) == N)
        {
            names(apresultObj@exemplars) <- colnames(s)[apresultObj@exemplars]

            for (c in 1:length(apresultObj@exemplars))
            {
                names(apresultObj@clusters[[c]]) <-
                    colnames(s)[apresultObj@clusters[[c]]]
            }
        }
    }
    else
    {
        tmpidx    <- matrix(NaN, N, 1)
        tmpnetsim <- NaN
        tmpexpref <- NaN

        apresultObj@exemplars <- c()
        apresultObj@clusters  <- list()
    }

    apresultObj@netsim <- tmpnetsim
    apresultObj@dpsim  <- tmpdpsim
    apresultObj@expref <- tmpexpref
    apresultObj@idx    <- tmpidx
    apresultObj@it     <- i

    if (details)
    {
        apresultObj@netsimAll[i+1] <- tmpnetsim
        apresultObj@dpsimAll[i+1]  <- tmpdpsim
        apresultObj@exprefAll[i+1] <- tmpexpref
        apresultObj@idxAll[,i+1]   <- tmpidx

        apresultObj@idxAll    <- apresultObj@idxAll[,1:i+1]
        apresultObj@netsimAll <- apresultObj@netsimAll[1:i+1]
        apresultObj@dpsimAll  <- apresultObj@dpsimAll[1:i+1]
        apresultObj@exprefAll <- apresultObj@exprefAll[1:i+1]
    }

    if (unconverged)
        warning("Algorithm did not converge. Turn on details\n",
                "and call plot() to monitor net similarity. Consider\n",
                "increasing maxits and convits, and, if oscillations occur\n",
                "also increasing damping factor lam.")

    apresultObj
}
