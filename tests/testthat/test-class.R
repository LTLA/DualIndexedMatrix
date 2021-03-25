# This checks the various aspects of the DIM class.
# library(testthat); library(DualIndexedMatrix); source("test-class.R")

library(Matrix)
by.column <- rsparsematrix(100, 200, 0.1)
by.row <- t(by.column)
mat <- DualIndexedMatrix(by.row, by.column, row.transposed=TRUE)

test_that("construction and validity methods work as expected", {
    expect_s4_class(mat, "DualIndexedMatrix")
    expect_identical(type(mat), "double")
    expect_true(is_sparse(mat))

    expect_error(DualIndexedMatrix(by.row, by.column), "inconsistent dimensions")

    expect_error(DualIndexedMatrix(by.row!=0, by.column, row.transposed=TRUE), "inconsistent type")

    expect_error(DualIndexedMatrix(as.matrix(by.row), by.column, row.transposed=TRUE), "inconsistent sparsity")

    mat2 <- DualIndexedMatrix(as.matrix(by.row), as.matrix(by.column), row.transposed=TRUE)
    expect_false(is_sparse(mat2))
})

test_that("extract_array methods work as expected", {
    COMPARE <- function(x, y, index) {
        expect_identical(extract_array(x, index), extract_array(y, index))
    }

    COMPARE(mat, by.column, list(1:10, NULL))
    COMPARE(mat, by.column, list(NULL, 1:10))
    COMPARE(mat, by.column, list(10:1, NULL))
    COMPARE(mat, by.column, list(NULL, 10:1))

    chosen <- c(2,6,4,10,8)
    COMPARE(mat, by.column, list(chosen, NULL))
    COMPARE(mat, by.column, list(NULL, chosen))

    chosen2 <- c(chosen, chosen)
    COMPARE(mat, by.column, list(chosen2, NULL))
    COMPARE(mat, by.column, list(NULL, chosen2))

    COMPARE(mat, by.column, list(5, NULL))
    COMPARE(mat, by.column, list(NULL, 5))

    COMPARE(mat, by.column, list(chosen, 1:10))
    COMPARE(mat, by.column, list(10:1, chosen))
})

test_that("extract_array methods work as expected", {
    COMPARE <- function(x, y, index) {
        X <- extract_sparse_array(x, index)
        Y <- extract_sparse_array(y, index)
        expect_identical(as(X, "dgCMatrix"), as(Y, "dgCMatrix")) # enforce ordering.
    }

    COMPARE(mat, by.column, list(1:10, NULL))
    COMPARE(mat, by.column, list(NULL, 1:10))
    COMPARE(mat, by.column, list(10:1, NULL))
    COMPARE(mat, by.column, list(NULL, 10:1))

    chosen <- c(2,6,4,10,8)
    COMPARE(mat, by.column, list(chosen, NULL))
    COMPARE(mat, by.column, list(NULL, chosen))

    chosen2 <- c(chosen, chosen)
    COMPARE(mat, by.column, list(chosen2, NULL))
    COMPARE(mat, by.column, list(NULL, chosen2))

    COMPARE(mat, by.column, list(5, NULL))
    COMPARE(mat, by.column, list(NULL, 5))

    COMPARE(mat, by.column, list(chosen, 1:10))
    COMPARE(mat, by.column, list(10:1, chosen))
})


