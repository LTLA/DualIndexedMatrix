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
    requested <- lengths(index)
    if (is.null(index[[1]])) {
        requested[1] <- nrow(x)
    }
    if (is.null(index[[2]])) {
        requested[2] <- ncol(x)
    }

    if (requested[1] > requested[2]) {
        # Fewer columns to extract than rows, so let's use a column format,
        # as this will presumably involve fewer random accesses; the required
        # rows can then be fished out cheaply from the extracted columns. 
        .extract_transposed(x@column, index, x@column.transposed, FUN=FUN)
    } else {
        # Conversely, fewer random row accesses required here.
        .extract_transposed(x@row, index, x@row.transposed, FUN=FUN)
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

