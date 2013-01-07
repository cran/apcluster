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

setMethod("[", signature(x="AggExResult"),
    function(x, i)
    {
        lapply(i, function(index) cutree(x, k=index))
    }
)

setMethod("[[", signature(x="AggExResult"),
    function(x, i)
    {
        cutree(x, k=i)
    }
)

setMethod("similarity", signature(x="APResult"), function(x) x@sim)

setMethod("similarity", signature(x="AggExResult"), function(x) x@sim)

setMethod("similarity", signature(x="ExClust"), function(x) x@sim)