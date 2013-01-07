heatmap.APResult <- function(x, y, ...)
{
    if(!identical(dim(x@sim), as.integer(c(1, 1))))
    {
        aggres <- heatmap(as(x, "ExClust"), x@sim, ...)
        return(invisible(aggres))
    }
    else
        stop("Similarity matrix is missing for heatmap plotting.\n",
             "       Please use the apcluster or apclusterL argument\n",
             "       includeSim=TRUE in order to add similarity data to\n",
             "       the APResult object or pass externally computed\n",
             "       similarity matrix as second parameter.")
}

setMethod("heatmap", signature(x="APResult", y="missing"), heatmap.APResult)


heatmap.APResult.matrix <- function(x, y, ...)
{
    aggres <- heatmap(as(x, "ExClust"), y, ...)
    return(invisible(aggres))
}

setMethod("heatmap", signature(x="APResult", y="matrix"),
          heatmap.APResult.matrix)


heatmap.ExClust <- function(x, y, ...)
{
    if(!identical(dim(x@sim), as.integer(c(1, 1))))
    {
        aggres <- heatmap(x, x@sim, ...)
        return(invisible(aggres))
    }
    else
        stop("Similarity matrix is missing for heatmap plotting.\n",
             "       Make sure that similarity data is part of\n",
             "       the ExClust object or pass externally computed\n",
             "       similarity matrix as second parameter.")
}

setMethod("heatmap", signature(x="ExClust", y="missing"), heatmap.ExClust)


heatmap.ExClust.matrix <- function(x, y, ...)
{
    if (identical(dim(y), as.integer(c(1, 1))))
        stop("y must be a nonempty similarity matrix")
    else if (x@l != nrow(y))
        stop("size of clustering result does not fit to size of data set")
    else if (nrow(y) != ncol(y) && length(x@sel) == 0)
        stop("similarity matrix must be quadratic")

    aggres <- aggExCluster(y, x)

    heatmap(aggres, y, ...)

    return(invisible(aggres))
}

setMethod("heatmap", signature(x="ExClust", y="matrix"), heatmap.ExClust.matrix)


heatmap.AggExResult <- function(x, y, ...)
{
    if(!identical(dim(x@sim), as.integer(c(1, 1))))
        heatmap(x, x@sim, ...)
    else
        stop("Similarity matrix is missing for heatmap plotting.\n",
             "       Use aggExCluster() argument includeSim=TRUE or pass\n",
             "       externally computed similarity matrix as second\n",
             "       parameter.")
}

setMethod("heatmap", signature(x="AggExResult", y="missing"),
          heatmap.AggExResult)


heatmap.AggExResult.matrix <- function(x, y, ...)
{
    if (identical(dim(y), as.integer(c(1, 1))))
        stop("y must be a nonempty matrix")
    else if (x@l != nrow(y))
        stop("size of clustering result does not fit to size of data set")
    else if (length(x@sel)==0 && ncol(y) != nrow(y))
        stop("y must be quadratic")
    else if (x@maxNoClusters < 3 || x@maxNoClusters < x@l)
    {
        if (length(x@sel) > 0)
        {
            colInd <- c(rep(0, nrow(y)))
            for (i in 1:length(x@sel))
            colInd[x@sel[i]] <- i
        }

        order <- unlist(x@clusters[[x@maxNoClusters]][x@order])
        colVec <- rainbow(x@maxNoClusters)
        colors <- unlist(lapply(1:x@maxNoClusters,
                                function(i)
                                rep(colVec[x@order[i]],
                                    length(x@clusters[[x@maxNoClusters]]
                                           [[x@order[i]]]))))

        ver <- getRversion()
        ver <- as.integer(c(ver[[c(1, 1)]], ver[[c(1, 2)]]))

        rowColors <- colors
        if (!(ver[1] > 2 || (ver[1] == 2 && ver[2] >= 15)))
            rowColors <- rev(colors)

        if (length(x@sel) > 0)
        {
            colColors <-
            unlist(lapply(1:x@maxNoClusters,
                          function(i)
                          rep(colVec[x@order[i]],
                              length(intersect(x@clusters[[x@maxNoClusters]]
                                               [[x@order[i]]],
                                               x@sel)))))

            stats::heatmap(y[order, colInd[intersect(order, x@sel)]], Rowv=NA,
                           Colv=NA, revC=TRUE, ColSideColors=colColors,
                           RowSideColors=rowColors, ...)
        }
        else
            stats::heatmap(y[order, order], symm=TRUE, Rowv=NA, Colv=NA,
                           revC=TRUE, ColSideColors=colors,
                           RowSideColors=rowColors, ...)
    }
    else
    {
        mini <- min(x@height)
        maxi <- max(x@height)
        auxH <- x@height <- 0.05 + 0.95 * (-x@height + maxi) / (maxi - mini)

        hCl <- list(merge=x@merge, height=auxH, labels=x@labels,
                    order=x@order)
        class(hCl) <- "hclust"
        dend <- as.dendrogram(hCl)

        stats::heatmap(y, symm=TRUE, Rowv=dend, Colv=dend, revC=TRUE, ...)

        return(invisible(dend))
    }
}

setMethod("heatmap", signature(x="AggExResult", y="matrix"),
          heatmap.AggExResult.matrix)

setMethod("heatmap", signature(x="matrix", y="missing"),
          function(x, y, ...) stats::heatmap(x=x, ...))

setMethod("heatmap", signature(x="missing", y="matrix"),
          function(x, y, ...) stats::heatmap(x=y, ...))
