sparseToFull <- function(s, fill=-Inf)
{
    if (length(dim(s)) != 2 || ncol(s) != 3)
        stop("'s' must be a matrix with 3 columns")

    if (min(s[, 1:2]) <= 0)
        stop("indices in 's' must be >= 1")

    N <- max(s[, 1:2])

    S <- matrix(fill, N, N)

    S[s[, 1] + N * (s[, 2] - 1)] <- s[, 3]

    S
}
