#' Dual-indexed matrix representation
#'
#' Implements a dual-indexed \linkS4class{DelayedMatrix} representation that enables efficient row and column access at the cost of doubling the data.
#' This is done by storing two seeds where each has been indexed to enable fast random access to the values of a single dimension.
#'
#' @param row A matrix-like object or a DelayedMatrix seed that enables fast access of row data.
#' @param column A matrix-like object or a DelayedMatrix seed that enables fast access of column data.
#' This should contain the same dataset as \code{row}.
#' @param row.transposed Logical scalar indicating whether \code{row} contains the data in a transposed form.
#' @param column.transposed Logical scalar indicating whether \code{column} contains the data in a transposed form.
#' @param seed A DualIndexedMatrixSeed object.
#'
#' @return
#' The \code{DualIndexedMatrixSeed} constructor will return an instance of a DualIndexedMatrixSeed object.
#'
#' The \code{DualIndexedMatrix} and \code{DelayedArray} constructors will return an instance of a DualIndexedMatrix object.
#'
#' @details
#' The idea behind this class is to store two copies of the same dataset, where each copy is indexed differently to enable access on different dimensions.
#' For a given access request, the class will then dynamically choose the more efficient representation for data extraction.
#' In this manner, we can achieve fast row and column access at the cost of doubling the dataset.
#' 
#' This class is mostly intended to be used with file-backed representations, where we can achieve major performance improvements by minimizing costly I/O.
#' The assumption is that disk space is cheap such that the doubling of the size/number of data files is not a major concern.
#' However, this class can also be used with in-memory representations that are indexed for optimal access in one dimension.
#'
#' Setting \code{row.transposed=TRUE} indicates that the rows of the original dataset are stored in the columns of \code{row}.
#' Similarly, setting \code{column.transposed=TRUE} indicates that the columns of the original dataset are stored in the rows of \code{column}.
#' These are occasionally necessary to adapt to representations that have inherent preferred access patterns, e.g., \linkS4class{dgCMatrix}es in the Example.
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
#' @docType class
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
#' @rdname DualIndexedMatrix-class
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
        # Extracting an entire column.
        .extract_transposed(x@column, index, x@column.transposed, FUN=FUN)
    } else if (is.null(index[[2]])) {
        # Extracting an entire row.
        .extract_transposed(x@row, index, x@row.transposed, FUN=FUN)
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

        # Computing the number of reads in either orientation. We multiply the number of reads
        # per row (i.e., the number of separate column runs) by the number of rows to extract,
        # and vice versa for the number of columns reads.
        nreads.row <- length(sanitized[[1]]) * costs[2]
        nreads.col <- costs[1] * length(sanitized[[2]])

        if (nreads.col <= nreads.row) {
            out <- .extract_transposed(x@column, sanitized, x@column.transposed, FUN=FUN)
        } else {
            out <- .extract_transposed(x@row, sanitized, x@row.transposed, FUN=FUN)
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
#' @rdname DualIndexedMatrix-class
setMethod("DelayedArray", "DualIndexedMatrixSeed", function(seed) new_DelayedArray(seed, Class="DualIndexedMatrix"))

#' @export
#' @rdname DualIndexedMatrix-class
DualIndexedMatrix <- function(row, column, row.transposed=FALSE, column.transposed=FALSE) {
    DelayedArray(DualIndexedMatrixSeed(row, column, row.transposed, column.transposed))
}

