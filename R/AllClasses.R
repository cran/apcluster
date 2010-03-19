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

