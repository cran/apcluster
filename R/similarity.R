negDistMat <- function(x, r=1, ...)
{
    dm <- as.matrix(dist(x, diag=TRUE, upper=TRUE, ...))

    if (r <= 0)
    {
         warning("ignoring r<=0; using default r=1")
         r <- 1
    }

    if (r != 1)
        dm <- apply(dm, c(1, 2), function(x){x^r})

    -dm
}

expSimMat <- function(x, r=2, w=1, ...)
{
    dm <- as.matrix(dist(x, diag=TRUE, upper=TRUE, ...))

    if (r <= 0)
    {
        warning("ignoring r<=0; using default r=1")
        r <- 1
    }

    if (w <= 0)
    {
        warning("ignoring w<=0; using default w=1")
        w <- 1
    }

    apply(dm, c(1,2), function(x){exp(-(x / w)^r)})
}

linSimMat <- function(x, w=1, ...)
{
    dm <- as.matrix(dist(x, diag=TRUE, upper=TRUE, ...))

    if (w <= 0)
    {
        warning("ignoring w<=0; using default w=1")
        w <- 1
    }

    apply(dm, c(1,2), function(x){max(0, 1 - x / w)})
}

linKernel <- function(x, normalize=FALSE)
{
    mat <- x %*% t(x)

    if (normalize)
    {
        di <- 1 / sqrt(diag(mat))
        di[which(is.infinite(di))] <- 0

        mat <- mat * (di %o% di)
    }

    if (length(dim(x)) == 0) rownames(mat) <- colnames(mat)

    mat
}
