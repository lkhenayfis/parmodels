setup_test_data <- function() {
    data("dummyseries", envir = environment())
    serie <- ts(dummyseries[, 2], start = c(1931, 1), frequency = 12)
    return(serie)
}

test_that("build_RHO", {
    f <- build_RHO
    expect_true(is.function(f))
    test_that("build_RHO returns a square matrix with dimensions equal to lag_max", {
        serie <- setup_test_data()
        lag_max <- 5
        result <- f(serie, m = 1, lag_max = lag_max, est = "n")
        expect_true(is.matrix(result))
        expect_equal(nrow(result), lag_max)
        expect_equal(ncol(result), lag_max)
    })
    test_that("build_RHO diagonal elements are all 1", {
        serie <- setup_test_data()
        lag_max <- 5
        result <- f(serie, m = 1, lag_max = lag_max, est = "n")
        expect_equal(diag(result), rep(1, lag_max))
    })
    test_that("build_RHO matrix is symmetric", {
        serie <- setup_test_data()
        lag_max <- 5
        result <- f(serie, m = 1, lag_max = lag_max, est = "n")
        expect_equal(result, t(result))
    })
    test_that("build_RHO with est='n' vs est='n-1' produces different results", {
        serie <- setup_test_data()
        lag_max <- 5
        result_n <- f(serie, m = 1, lag_max = lag_max, est = "n")
        result_n1 <- f(serie, m = 1, lag_max = lag_max, est = "n-1")
        expect_false(identical(result_n, result_n1))
    })
    test_that("build_RHO works with different seasons", {
        serie <- setup_test_data()
        lag_max <- 5
        result_m1 <- f(serie, m = 1, lag_max = lag_max, est = "n")
        result_m6 <- f(serie, m = 6, lag_max = lag_max, est = "n")
        result_m12 <- f(serie, m = 12, lag_max = lag_max, est = "n")
        expect_true(is.matrix(result_m1))
        expect_true(is.matrix(result_m6))
        expect_true(is.matrix(result_m12))
        expect_false(identical(result_m1, result_m6))
        expect_false(identical(result_m6, result_m12))
    })
    test_that("build_RHO returns finite numeric values", {
        serie <- setup_test_data()
        lag_max <- 5
        result <- f(serie, m = 1, lag_max = lag_max, est = "n")
        expect_true(is.numeric(result))
        expect_true(all(is.finite(result)))
    })
})

test_that("build_rho", {
    f <- build_rho
    expect_true(is.function(f))
    test_that("build_rho returns a vector with length equal to lag_max", {
        serie <- setup_test_data()
        lag_max <- 5
        result <- f(serie, m = 1, lag_max = lag_max, est = "n")
        expect_true(is.numeric(result))
        expect_equal(length(result), lag_max)
    })
    test_that("build_rho returns finite numeric values", {
        serie <- setup_test_data()
        lag_max <- 5
        result <- f(serie, m = 1, lag_max = lag_max, est = "n")
        expect_true(is.numeric(result))
        expect_true(all(is.finite(result)))
    })
    test_that("build_rho with est='n' vs est='n-1' produces different results", {
        serie <- setup_test_data()
        lag_max <- 5
        result_n <- f(serie, m = 1, lag_max = lag_max, est = "n")
        result_n1 <- f(serie, m = 1, lag_max = lag_max, est = "n-1")
        expect_false(identical(result_n, result_n1))
    })
    test_that("build_rho works with different seasons", {
        serie <- setup_test_data()
        lag_max <- 5
        result_m1 <- f(serie, m = 1, lag_max = lag_max, est = "n")
        result_m6 <- f(serie, m = 6, lag_max = lag_max, est = "n")
        result_m12 <- f(serie, m = 12, lag_max = lag_max, est = "n")
        expect_true(is.numeric(result_m1))
        expect_true(is.numeric(result_m6))
        expect_true(is.numeric(result_m12))
        expect_false(identical(result_m1, result_m6))
        expect_false(identical(result_m6, result_m12))
    })
    test_that("build_rho works with different lag_max values", {
        serie <- setup_test_data()
        result_3 <- f(serie, m = 1, lag_max = 3, est = "n")
        result_7 <- f(serie, m = 1, lag_max = 7, est = "n")
        expect_equal(length(result_3), 3)
        expect_equal(length(result_7), 7)
    })
})

test_that("perpacf", {
    f <- perpacf
    expect_true(is.function(f))
    test_that("perpacf returns an object of class periodic_pacf", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, lag_max = 5, est = "n", plot = FALSE)
        expect_true(inherits(result, "periodic_pacf"))
        expect_true(is.list(result))
    })
    test_that("perpacf returned list has all required elements", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, lag_max = 5, est = "n", plot = FALSE)
        expect_true("phi" %in% names(result))
        expect_true("n_used" %in% names(result))
        expect_true("m" %in% names(result))
        expect_true("lag_max" %in% names(result))
        expect_true("rho" %in% names(result))
        expect_true("RHO" %in% names(result))
    })
    test_that("perpacf phi vector has length equal to lag_max", {
        serie <- setup_test_data()
        lag_max <- 5
        result <- f(serie, m = 1, lag_max = lag_max, est = "n", plot = FALSE)
        expect_equal(length(result$phi), lag_max)
    })
    test_that("perpacf n_used equals floor(length(serie) / frequency(serie))", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, lag_max = 5, est = "n", plot = FALSE)
        expected_n_used <- floor(length(serie) / frequency(serie))
        expect_equal(result$n_used, expected_n_used)
    })
    test_that("perpacf m parameter is preserved in output", {
        serie <- setup_test_data()
        result <- f(serie, m = 3, lag_max = 5, est = "n", plot = FALSE)
        expect_equal(result$m, 3)
    })
    test_that("perpacf with default lag_max", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, est = "n", plot = FALSE)
        expect_equal(result$lag_max, frequency(serie) - 1)
        expect_equal(length(result$phi), frequency(serie) - 1)
    })
    test_that("perpacf with custom lag_max value", {
        serie <- setup_test_data()
        lag_max <- 6
        result <- f(serie, m = 1, lag_max = lag_max, est = "n", plot = FALSE)
        expect_equal(result$lag_max, lag_max)
    })
    test_that("perpacf with est='n' via match.arg", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, lag_max = 5, est = "n", plot = FALSE)
        expect_true(is.numeric(result$phi))
        expect_true(all(is.finite(result$phi)))
    })
    test_that("perpacf with est='n-1'", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, lag_max = 5, est = "n-1", plot = FALSE)
        expect_true(is.numeric(result$phi))
        expect_true(all(is.finite(result$phi)))
    })
    test_that("perpacf with plot=FALSE doesn't error", {
        serie <- setup_test_data()
        expect_silent(f(serie, m = 1, lag_max = 5, est = "n", plot = FALSE))
    })
    test_that("perpacf phi values are finite", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, lag_max = 5, est = "n", plot = FALSE)
        expect_true(all(is.finite(result$phi)))
    })
    test_that("perpacf RHO and rho match helper functions", {
        serie <- setup_test_data()
        lag_max <- 5
        m <- 1
        result <- f(serie, m = m, lag_max = lag_max, est = "n", plot = FALSE)
        expected_RHO <- build_RHO(serie, m = m, lag_max = lag_max, est = "n")
        expected_rho <- build_rho(serie, m = m, lag_max = lag_max, est = "n")
        expect_equal(result$RHO, expected_RHO)
        expect_equal(result$rho, expected_rho)
    })
})
