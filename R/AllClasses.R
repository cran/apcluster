# S4 class definition for the result object of affinity propagation clustering
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
        idxAll    = matrix(NA, 1, 1)
    )
)

# S4 class definition for the result object of the aggExCluster algorithm
setClass("AggExResult",
    representation = representation
    (
        l             = "numeric",
        maxNoClusters = "numeric",
        clusters      = "list",
        exemplars     = "list",
        merge         = "matrix",
        height        = "numeric",
        order         = "numeric",
        labels        = "character"
    ),
    prototype = prototype
    (
        l             = 0,
        maxNoClusters = 0,
        clusters      = list(),
        exemplars     = list(),
        merge         = matrix(NA, 1, 1),
        height        = c(),
        order         = c(),
        labels        = c()
    )
)


# S4 class definition for exemplar-based clustering
setClass("ExClust",
    representation = representation
    (
        l         = "numeric",
        exemplars = "numeric",
        clusters  = "list",
        idx       = "numeric"
    ),
    prototype = prototype
    (
        l         = 0,
        exemplars = c(),
        clusters  = list(),
        idx       = c()
    )
)


