% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DualIndexedMatrix.R
\docType{class}
\name{DualIndexedMatrix-class}
\alias{DualIndexedMatrix-class}
\alias{DualIndexedMatrix}
\alias{DualIndexedMatrixSeed}
\alias{DualIndexedMatrixSeed-class}
\alias{dim,DualIndexedMatrixSeed-method}
\alias{dimnames,DualIndexedMatrixSeed-method}
\alias{extract_array,DualIndexedMatrixSeed-method}
\alias{extract_sparse_array,DualIndexedMatrixSeed-method}
\alias{type,DualIndexedMatrixSeed-method}
\alias{is_sparse,DualIndexedMatrixSeed-method}
\alias{DelayedArray,DualIndexedMatrixSeed-method}
\title{Dual-indexed matrix representation}
\usage{
DualIndexedMatrixSeed(
  row,
  column,
  row.transposed = FALSE,
  column.transposed = FALSE
)

\S4method{DelayedArray}{DualIndexedMatrixSeed}(seed)

DualIndexedMatrix(
  row,
  column,
  row.transposed = FALSE,
  column.transposed = FALSE
)
}
\arguments{
\item{row}{A matrix-like object or a DelayedMatrix seed that enables fast access of row data.}

\item{column}{A matrix-like object or a DelayedMatrix seed that enables fast access of column data.
This should contain the same dataset as \code{row}.}

\item{row.transposed}{Logical scalar indicating whether \code{row} contains the data in a transposed form.}

\item{column.transposed}{Logical scalar indicating whether \code{column} contains the data in a transposed form.}

\item{seed}{A DualIndexedMatrixSeed object.}
}
\value{
The \code{DualIndexedMatrixSeed} constructor will return an instance of a DualIndexedMatrixSeed object.

The \code{DualIndexedMatrix} and \code{DelayedArray} constructors will return an instance of a DualIndexedMatrix object.
}
\description{
Implements a dual-indexed \linkS4class{DelayedMatrix} representation that enables efficient row and column access at the cost of doubling the data.
This is done by storing two seeds where each has been indexed to enable fast random access to the values of a single dimension.
}
\details{
The idea behind this class is to store two copies of the same dataset, where each copy is indexed differently to enable access on different dimensions.
For a given access request, the class will then dynamically choose the more efficient representation for data extraction.
In this manner, we can achieve fast row and column access at the cost of doubling the dataset.

This class is mostly intended to be used with file-backed representations, where we can achieve major performance improvements by minimizing costly I/O.
The assumption is that disk space is cheap such that the doubling of the size/number of data files is not a major concern.
However, this class can also be used with in-memory representations that are indexed for optimal access in one dimension.

Setting \code{row.transposed=TRUE} indicates that the rows of the original dataset are stored in the columns of \code{row}.
Similarly, setting \code{column.transposed=TRUE} indicates that the columns of the original dataset are stored in the rows of \code{column}.
These are occasionally necessary to adapt to representations that inherently favor access on one dimension, e.g., compressed sparse column matrices.
}
\examples{
library(HDF5Array)
common.mat <- matrix(runif(2000000), ncol=1000)

rowpath <- tempfile(fileext=".hdf5")
invisible(writeHDF5Array(common.mat, rowpath, name="data", chunkdim=c(1, ncol(common.mat))))

colpath <- tempfile(fileext=".hdf5")
invisible(writeHDF5Array(common.mat, colpath, name="data", chunkdim=c(nrow(common.mat), 1)))

# We store both in the DualIndexedMatrix, which will determine
# which to use for best performance for a given access pattern.
mat <- DualIndexedMatrix(HDF5ArraySeed(rowpath, "data"), HDF5ArraySeed(colpath, "data"))

# Another example with sparse matrices. In this case, we transpose the data and set
# row.transposed=TRUE so that we can keep using a dgCMatrix for rapid row access.
library(Matrix)
by.column <- rsparsematrix(100, 200, 0.1)
class(by.column) # dgCMatrices enable efficient column access.

by.row <- t(by.column) 
class(by.row)

mat <- DualIndexedMatrix(by.row, by.column, row.transposed=TRUE) 

}
\author{
Aaron Lun
}
