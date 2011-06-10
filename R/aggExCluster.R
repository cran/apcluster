aggExCluster <- function(s, cl)
{
    if (length(dim(s)) != 2 || ncol(s) != nrow(s))
        stop("s must be a square matrix")

    AggResultObj <- new("AggExResult")

    K <- ncol(s)
    AggResultObj@l <- K

    preserveNames <- (length(colnames(s)) == ncol(s))

    if (missing(cl) || is.null(cl)) ## no prior clustering
    {
        AggResultObj@maxNoClusters <- K
        AggResultObj@clusters[[K]] <- as.list(1:K)
        AggResultObj@exemplars[[K]] <- 1:K

        if (preserveNames)
        {
            AggResultObj@labels <- colnames(s)
            names(AggResultObj@exemplars[[K]]) <- colnames(s)

            for (i in 1:K)
                names(AggResultObj@clusters[[K]][[i]]) <- colnames(s)[i]
        }
        else
            AggResultObj@labels <- as.character(1:K)
    }
    else ## prior clustering
    {
        if (class(cl) != "APResult" && class(cl) != "ExClust")
            stop("cl must be missing or object of classes 'APResult'",
                 " or 'ExClust'")
        else if (cl@l != ncol(s))
            stop("data set sizes of s and cl do not match")

        K <- length(cl@exemplars)

        if (K < 1)
            stop("cl empty or corrupted")

        AggResultObj@maxNoClusters <- K
        AggResultObj@clusters[[K]] <- cl@clusters
        AggResultObj@exemplars[[K]] <- cl@exemplars
        AggResultObj@labels <- paste("Cluster", 1:K)
    }

    if (K < 2)
    {
        warning("there is nothing to cluster")
        return(invisible(AggResultObj))
    }

    objMat <- matrix(NA, K, K) ## matrix of objective values for pairs
    exeMat <- matrix(NA, K, K) ## matrix of joint exemplars
    ## note: only the upper triangle of these matrices is non-NA

    actClust <- AggResultObj@clusters[[K]]
    actExem <- AggResultObj@exemplars[[K]]
    actLabels <- -(1:K)

    AggResultObj@merge <- matrix(NA, K - 1, 2)
    AggResultObj@height <- rep(0, K - 1)

    ## compute complete matrices before starting joining
    for (i in 1:(K - 1))
    {
        for (j in (i + 1):K)
        {
             joint <- c(actClust[[i]], actClust[[j]])

             cM <- colMeans(s[joint, joint, drop=FALSE])
             ex <- joint[which.max(cM)]
             exeMat[i, j] <- ex

             objMat[i, j] <- (mean(s[ex, actClust[[i]]]) +
                              mean(s[ex, actClust[[j]]])) / 2
        }
    }

    ## agglomeration loop
    for (k in (K - 1):1)
    {
        tojoin <- which.max(objMat) - 1 ## determine pair to join
        I <- tojoin %% K + 1
        J <- floor(tojoin / K) + 1

        newClust <- c(actClust[[I]], actClust[[J]]) ## join 'em
        actClust[c(I, J)] <- NULL
        actClust[[k]] <- newClust
        actExem <- c(actExem[c(-I, -J)], exeMat[I, J])

        AggResultObj@clusters[[k]] <- actClust
        AggResultObj@merge[K - k, ] <- c(actLabels[I], actLabels[J])
        actLabels <- c(actLabels[c(-I, -J)], K - k)
        AggResultObj@height[K - k] <- objMat[I, J]
        AggResultObj@exemplars[[k]] <- actExem

        if (preserveNames)
            names(AggResultObj@exemplars[[k]]) <- colnames(s)[actExem]

        if (k == 1) break

        ## rearrange matrices objMat and exeMat
        ## put values for unchanged clusters in the first k-1 rows/columns
        indexVec <- 1:(k + 1)
        indexVec <- indexVec[c(-I, -J)]

        exeMat[1:(k - 1), 1:(k - 1)] <- exeMat[indexVec, indexVec, drop=FALSE]
        objMat[1:(k - 1), 1:(k - 1)] <- objMat[indexVec, indexVec, drop=FALSE]

        ## wipe out k+1-st column
        exeMat[, k + 1] <- NA
        objMat[, k + 1] <- NA

        ## update k-th column with objective values and joint exemplars of
        ## unchanged clusters and the newly joined cluster
        for (i in 1:(k - 1))
        {
             joint <- c(actClust[[i]], actClust[[k]])

             cM <- colMeans(s[joint, joint, drop=FALSE])
             ex <- joint[which.max(cM)]
             exeMat[i, k] <- ex

             objMat[i, k] <- (mean(s[ex, actClust[[i]]]) +
                              mean(s[ex, actClust[[k]]])) / 2
        }
    }

    ## finally, determine reordering for dendrogram plotting
    AggResultObj@order <- determineOrder(AggResultObj@merge,
                                         AggResultObj@height, K - 1)

    AggResultObj
}

## auxiliary function for determining the order for dendrogram plotting
## fills up order recursively starting from the last merge
determineOrder <- function(merge, height, k)
{
    I <- merge[k, 1] ## I and J are the clusters merged in the k-th step
    J <- merge[k, 2]

    if (I < 0 && J < 0) ## if both are singletons, list I first
        return(c(-I, -J))
    else if (I < 0) ## if I is a singleton and J is not, list it first
        return(c(-I, determineOrder(merge, height, J)))
    else if (J < 0) ## if J is a singleton and I is not, list it first
        return(c(-J, determineOrder(merge, height, I)))
    else ## if both are non-singleton clusters, list the "tighter" cluster
    {    ## on the left-hand side (see ?hclust)
        if (height[I] > height[J])
            return(c(determineOrder(merge, height, I),
                     determineOrder(merge, height, J)))
        else
            return(c(determineOrder(merge, height, J),
                     determineOrder(merge, height, I)))
    }
}
