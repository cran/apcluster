setMethod("[", signature(x="APResult"),
    function(x, i)
    {
        x@clusters[i]
    }
)

setMethod("[[", signature(x="APResult"),
    function(x, i)
    {
        x@clusters[[i]]
    }
)

setMethod("[", signature(x="ExClust"),
    function(x, i)
    {
        x@clusters[i]
    }
)

setMethod("[[", signature(x="ExClust"),
    function(x, i)
    {
        x@clusters[[i]]
    }
)
