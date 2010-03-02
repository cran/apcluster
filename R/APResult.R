# S4 class definition for the result object
setClass("APResult",
    representation = representation
    (
        l         = "numeric",
        it        = "numeric",
        p         = "numeric",
        netsim    = "numeric",
        dpsim     = "numeric",
        expref    = "numeric",
        exemplars = "numeric",
        clusters  = "list",
        idx       = "numeric",
        netsimAll = "numeric",
        dpsimAll  = "numeric",
        exprefAll = "numeric",
        idxAll    = "matrix"
    ),
    prototype = prototype
    (
        l         = 0,
        it        = 0,
        p         = 0,
        netsim    = NaN,
        dpsim     = NaN,
        expref    = NaN,
        exemplars = c(),
        clusters  = list(),
        idx       = c(),
        netsimAll = NaN,
        dpsimAll  = NaN,
        exprefAll = NaN,
        idxAll    = matrix(NA,1,1)
    )
)

# Display clustering results
setMethod("show", signature(object="APResult"),
    function(object)
    {
        cat("\nAPResult object\n")

        if (!is.finite(object@l) || !is.finite(object@it))
        {
            stop("Object is not result of an affinity propagation run.",
                 "It is pointless to create APResult objects yourself.")
        }

        cat("\nNumber of samples    = ", object@l, "\n")
        cat("Number of iterations = ", object@it, "\n")
        cat("Input preference     = ", object@p, "\n")
        cat("Sum of similarities  = ", object@dpsim, "\n")
        cat("Sum of preferences   = ", object@expref, "\n")
        cat("Net similarity       = ", object@netsim, "\n")
        cat("Number of exemplars  = ", length(object@exemplars), "\n\n")

        if (length(object@exemplars) > 0)
        {
            cat("Exemplars:\n")
            cat(object@exemplars, fill=TRUE, labels="  ")
            cat("Clusters:\n")

            for (i in 1:length(object@exemplars))
            {
                cat("   Cluster ", i, ", exemplar ", object@exemplars[i], ":\n",
                    sep="")
                cat(object@clusters[[i]], fill=TRUE, labels="     ")
            }
        }
        else
        {
            cat("No clusters identified.\n")
        }
    }
)

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
                 "       order to add plotting data to the APResult object.",
                 call.=FALSE);
        }
    }
)

# Plot clustering result along with data set (works only for 2D data)
setMethod("plot", signature(x="APResult", y="matrix"),
    function(x, y, connect=TRUE, xlab="", ylab="", ...)
    {
        if (length(dim(y)) != 2 || ncol(y) != 2)
        {
            stop("Second argument y must be a two-column matrix\n")
        }
        else if (x@l != nrow(y))
        {
            stop("Size of clustering result does not fit to size of data set\n")
        }

        xlim <- c(min(y[,1]), max(y[,1]))
        ylim <- c(min(y[,2]), max(y[,2]))

        plot(x=NULL, y=NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)

        num <- length(x@exemplars)

        if (num <= 0)
        {
            warning("No exemplars defined in clustering result. Plotting ",
                    "data set as it is.\n", call.=FALSE)

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
)
