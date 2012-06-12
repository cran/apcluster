# Plot graph(s) with objective values (works only if details were switched on)
setMethod("plot", signature(x="APResult", y="missing"),
    function(x, type=c("netsim", "dpsim", "expref"),
             xlab="# Iterations",
             ylab="Similarity", ...)
    {
        plotnetsim <- FALSE
        plotexpref <- FALSE
        plotdpsim  <- FALSE

        legtxt <- c()
        legcol <- c()

        ymin <- .Machine$double.xmax
        ymax <- -.Machine$double.xmax

        if (is.element("netsim", type))
        {
            tmp = x@netsimAll[which(!is.nan(x@netsimAll))]

            if (length(tmp) > 0)
            {
                ymin <- min(tmp, ymin, na.rm = TRUE)
                if (ymin == -Inf) ymin <- -.Machine$double.xmax

                ymax <- max(tmp, ymax, na.rm = TRUE)
                if (ymax == Inf) ymax <- .Machine$double.xmax

                plotnetsim <- TRUE

                legtxt <- c(legtxt, "Fitness (overall net similarity)")
                legcol <- c(legcol, "red")
            }
        }

        if (is.element("expref", type))
        {
            tmp = x@exprefAll[which(!is.nan(x@exprefAll))]

            if (length(tmp) > 0)
            {
                ymin <- min(tmp, ymin, na.rm = TRUE)
                if (ymin == -Inf) ymin <- -.Machine$double.xmax

                ymax <- max(tmp, ymax, na.rm = TRUE)
                if (ymax == Inf) ymax <- .Machine$double.xmax

                plotexpref <- TRUE

                legtxt <- c(legtxt, "Sum of exemplar preferences")
                legcol <- c(legcol, "green")
            }
        }

        if (is.element("dpsim", type))
        {
            tmp = x@dpsimAll[which(!is.nan(x@dpsimAll))]

            if (length(tmp) > 0)
            {
                ymin <- min(tmp, ymin, na.rm = TRUE)
                if (ymin == -Inf) ymin <- -.Machine$double.xmax

                ymax <- max(tmp, ymax, na.rm = TRUE)
                if (ymax == Inf) ymax <- .Machine$double.xmax

                plotdpsim <- TRUE

                legtxt <- c(legtxt, "Sum of similarities to exemplars")
                legcol <- c(legcol, "blue")
            }
        }

        if (length(legtxt) > 0)
        {
            plot(x=NULL, y=NULL,
                 xlim=c(0, x@it + 1), ylim=c(ymin, ymax),
                 xlab=xlab, ylab=ylab, ...)

            if (plotnetsim) lines(x@netsimAll, col="red")
            if (plotexpref) lines(x@exprefAll, col="green")
            if (plotdpsim)  lines(x@dpsimAll, col="blue")

            legend(x="bottomright", legend=legtxt, col=legcol, lwd=1)
        }
        else
        {
            stop("No valid data was found for plotting.\n",
                 "       Please use the apcluster argument details=TRUE in\n",
                 "       order to add plotting data to the APResult object.")
        }
    }
)

# Plot clustering result along with data set
setMethod("plot", signature(x="APResult", y="matrix"),
    function(x, y, connect=TRUE, xlab="", ylab="", ...)
    {
        plot(as(x, "ExClust"), y, connect, xlab, ylab, ...)
    }
)

setMethod("plot", signature(x="ExClust", y="matrix"),
    function(x, y, connect=TRUE, xlab="", ylab="", ...)
    {
        if (length(dim(y)) != 2)
            stop("y must be a matrix")
        else if (x@l != nrow(y))
            stop("size of clustering result does not fit to size of data set")
        else if (ncol(y) == 2)
        {
            xlim <- c(min(y[,1]), max(y[,1]))
            ylim <- c(min(y[,2]), max(y[,2]))

            plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,
                 ...)

            num <- length(x@exemplars)

            if (num <= 0)
            {
                warning("No exemplars defined in clustering result. Plotting ",
                        "data set as it is.")

                points(y[,1], y[,2], col="black", pch=19, cex=0.8)
            }
            else
            {
                cols <- rainbow(num)

                for (i in 1:num)
                {
                    data <- y[x@clusters[[i]],,drop=FALSE]
                    exem <- y[x@exemplars[i],]

                    points(data[,1], data[,2], col=cols[i],
                           pch=19, cex=0.8)

                    if (connect)
                    {
                        for (j in 1:nrow(data))
                        {
                            lines(c(data[j,1], exem[1]), c(data[j,2], exem[2]),
                                  col=cols[i])
                        }
                    }

                    points(exem[1], exem[2], col="black",
                           type="p", pch=22, cex=1.5)
                }
            }
        }
        else if (nrow(y) == ncol(y))
        {
            aggres <- aggExCluster(y, x)

            plot(aggres, y, ...)

            return(invisible(aggres))
        }
        else
            stop("y must be quadratic or two-column")
    }
)


# Plot clustering result
setMethod("plot", signature(x="AggExResult", y="missing"),
    function(x, main="Cluster dendrogram", xlab="",
             ylab="Balanced avg. similarity to exemplar", ticks=4,
             digits=2, ...)
    {
        if (x@maxNoClusters < 3)
            stop("cannot plot dendrogram with less than 3 clusters")

        mini <- min(x@height)
        maxi <- max(x@height)
        auxH <- x@height <- 0.05 + 0.95 * (-x@height + maxi) / (maxi - mini)

        hCl <- list(merge=x@merge, height=auxH, labels=x@labels, order=x@order)
        class(hCl) <- "hclust"
        dend <- as.dendrogram(hCl)

        plot(dend, axes=FALSE, xlab=xlab, ylab=ylab, main=main, ...)

        axis(side=2, at=seq(0.05, 1, length=ticks), tick=TRUE,
             labels=as.character(format(seq(maxi, mini, length=ticks),
                                        digits=digits)))

        return(invisible(dend))
    }
)


# Plot clustering result along with data set
setMethod("plot", signature(x="AggExResult", y="matrix"),
    function(x, y, k=NA, h=NA, ...)
    {
        if (length(dim(y)) != 2)
            stop("y must be a matrix")
        else if (x@l != nrow(y))
            stop("size of clustering result does not fit to size of data set")
        else if (ncol(y) == 2)
        {
            if (is.na(k) || !is.numeric(k) || k > x@maxNoClusters)
                k <- x@maxNoClusters

            if (k< 1)
                k <- 1

            excl <- cutree(x, k, h)

            plot(excl, y, ...)

            return(invisible(excl))
        }
        else if (ncol(y) != nrow(y))
            stop("y must be quadratic or two-column")
        else if (x@maxNoClusters < 3 || x@maxNoClusters < x@l)
        {
            order <- unlist(x@clusters[[x@maxNoClusters]][x@order])
            colVec <- rainbow(x@maxNoClusters)
            colors <- unlist(lapply(1:x@maxNoClusters,
                function(i)
                    rep(colVec[x@order[i]],
                        length(x@clusters[[x@maxNoClusters]][[x@order[i]]]))))

            ver <- getRversion()
            ver <- as.integer(c(ver[[c(1, 1)]], ver[[c(1, 2)]]))

            if (ver[1] > 2 || (ver[1] == 2 && ver[2] >= 15))
                heatmap(y[order, order], symm=TRUE, Rowv=NA, Colv=NA, revC=TRUE,
                        ColSideColors=colors, RowSideColors=colors, ...)
            else
                heatmap(y[order, order], symm=TRUE, Rowv=NA, Colv=NA, revC=TRUE,
                        ColSideColors=colors, RowSideColors=rev(colors), ...)
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

            heatmap(y, symm=TRUE, Rowv=dend,
                    Colv=dend, revC=TRUE, ...)

            return(invisible(dend))
        }
    }
)

