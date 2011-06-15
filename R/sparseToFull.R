sparseToFull <- function(s, fill=-Inf)
{
    if (length(dim(s)) != 2 || ncol(s) != 3)
        stop("s must be a 2D matrix with 3 columns")

    if (min(min(s[,1]), min(s[,2])) <= 0)
        stop("data point indices in s must be >= 1")

    N <- max(max(s[,1]),max(s[,2]))

    S <- matrix(fill, N, N)

    for (j in 1:nrow(s))
        S[s[j,1],s[j,2]] <- s[j,3]

    S
}
