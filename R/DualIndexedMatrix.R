#' Dual-indexed matrix representation
#'
#' Implements a dual-indexed \linkS4class{DelayedMatrix} representation that enables efficient row and column access at the cost of doubling the data.
#'
#'
#' @author Aaron Lun
#'
#' @examples
#' library(Matrix)
#' by.column <- rsparsematrix(100, 200, 0.1)
#' class(by.column) # dgCMatrices enable efficient column access.
#'
#' # When we transpose, we obtain a new dgCMatrix where the columns
#' # represent the original rows (and thus enable fast row access).
#' by.row <- t(by.column) 
#' class(by.row)
#'
#' # We store both in the DualIndexedMatrix, which will determine
#' # which to use for best performance for a given access pattern.
#' mat <- DualIndexedMatrix(by.row, by.column, row.transposed=TRUE) 
#'
#' @name DualIndexedMatrix-class
#' @aliases
#' DualIndexedMatrix
#' DualIndexedMatrix-class
#' DualIndexedMatrixSeed
#' DualIndexedMatrixSeed-class
#' dim,DualIndexedMatrixSeed-method
#' dimnames,DualIndexedMatrixSeed-method
#' extract_array,DualIndexedMatrixSeed-method
#' extract_sparse_array,DualIndexedMatrixSeed-method
#' type,DualIndexedMatrixSeed-method
#' is_sparse,DualIndexedMatrixSeed-method
NULL

#' @export
setClass("DualIndexedMatrixSeed", slots=c(row="ANY", column="ANY", row.transposed="logical", column.transposed="logical"))

#' @export
DualIndexedMatrixSeed <- function(row, column, row.transposed=FALSE, column.transposed=FALSE) {
    new("DualIndexedMatrixSeed", row=row, column=column, row.transposed=row.transposed, column.transposed=column.transposed)
}

setValidity2("DualIndexedMatrixSeed", function(object) {
    rdim <- dim(object@row)
    if (object@row.transposed) {
        rdim <- rev(rdim)
    }

    cdim <- dim(object@column)
    if (object@column.transposed) {
        cdim <- rev(cdim)
    }

    if (length(rdim)!=2 || length(cdim)!=2) {
        return("row and column seeds should refer to 2-dimensional arrays")
    }

    if (!identical(rdim, cdim)) {
        return("inconsistent dimensions for row and column seeds")
    }

    if (type(rdim)!=type(cdim)) {
        return("inconsistent types for row and column seeds")
    }

    if (is_sparse(rdim)!=is_sparse(cdim)) {
        return("inconsistent sparsity for row and column seeds")
    }

    TRUE
})

#' @export
setMethod("type", "DualIndexedMatrixSeed", function(x) type(x@row))

#' @export
setMethod("is_sparse", "DualIndexedMatrixSeed", function(x) is_sparse(x@row))

#' @export
setMethod("dim", "DualIndexedMatrixSeed", function(x) {
    d <- dim(x@row)
    if (x@row.transposed) {
        d <- rev(d)
    }
    d
})

#' @export
setMethod("dimnames", "DualIndexedMatrixSeed", function(x) {
    d <- dimnames(x@row)
    if (x@row.transposed) {
        d <- rev(d)
    }
    d
})

#' @export
setMethod("extract_sparse_array", "DualIndexedMatrixSeed", function(x, index) .dual_extract(x, index, FUN=extract_sparse_array))

#' @export
setMethod("extract_array", "DualIndexedMatrixSeed", function(x, index) .dual_extract(x, index, FUN=extract_array))

.dual_extract <- function(x, index, FUN) {
    if (is.null(index[[1]])) {
        .extract_transposed(x@row, index, x@row.transposed, FUN=FUN)
    } else if (is.null(index[[2]])) {
        .extract_transposed(x@column, index, x@column.transposed, FUN=FUN)
    } else {
        costs <- integer(length(index))
        sanitized <- index

        for (i in seq_along(sanitized)) {
            # Guaranteed to be non-NULL at this point, no need to test.
            current <- sanitized[[i]]
            if (!all(diff(current)==1L)) {
                sanitized[[i]] <- sort(unique(current))
                costs[i] <- sum(diff(index[[i]])!=1L) + 1L # number of non-consecutive runs.
            } else {
                costs[i] <- 1L
            }
        }

        if (costs[1] <= costs[2]) {
            out <- .extract_transposed(x@row, sanitized, x@row.transposed, FUN=FUN)
        } else {
            out <- .extract_transposed(x@column, sanitized, x@column.transposed, FUN=FUN)
        }

        # Quick and dirty reorganization for the reordering that we had to do.
        reindex <- vector("list", length(index))
        for (i in seq_along(index)) {
            if (costs[i] != 1L && !identical(index[[i]], sanitized[[i]])) {
                reindex[[i]] <- match(index[[i]], sanitized[[i]])
            }
        }

        if (is.null(reindex[[1]]) && is.null(reindex[[2]])) {
            out
        } else {
            FUN(out, reindex)
        }
    }
}

.extract_transposed <- function(seed, index, transposed, FUN) {
    if (transposed) {
        index <- rev(index)
    }

    out <- FUN(seed, index)

    if (transposed) {
        out <- t(out)
    }
    out
}

#' @export
setClass("DualIndexedMatrix",
    contains="DelayedMatrix",
    representation(seed="DualIndexedMatrixSeed")
)

#' @export
setMethod("DelayedArray", "DualIndexedMatrixSeed", function(seed) new_DelayedArray(seed, Class="DualIndexedMatrix"))

#' @export
DualIndexedMatrix <- function(row, column, row.transposed=FALSE, column.transposed=FALSE) {
    DelayedArray(DualIndexedMatrixSeed(row, column, row.transposed, column.transposed))
}

