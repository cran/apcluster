preferenceRange <- function(s, exact=FALSE)
{
    if (length(dim(s)) != 2 || ncol(s) != nrow(s))
    {
        stop("s must be a square matrix")
    }

    N <- nrow(s)

    diag(s) <- 0

    dpsim1 <- max(colSums(s))

    if (dpsim1 == -Inf)
    {
        pmin <- NaN
    }
    else if (exact)
    {
        dpsim2 <- -Inf

        for (j21 in 1:(N - 1))
        {
            for (j22 in (j21 + 1):N)
            {
                tmp <- sum(pmax(s[,j21], s[,j22]))

                if (tmp > dpsim2) dpsim2 <- tmp
            }
        }

        pmin <- dpsim1 - dpsim2
    
        diag(s) <- -Inf
    }
    else
    {
        diag(s) <- -Inf

        m <- apply(s, 1, max)

        ii <- which.min(m)
        yy <- m[ii]

        pmin <- dpsim1 - sum(m) + yy + min(m[-ii])
    }

    as.vector(c(pmin, max(s)))
}
