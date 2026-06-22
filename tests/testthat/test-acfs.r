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
        detrended_serie <- scale_by_season(serie)[[1]]
        lag_max <- 5
        m <- 1
        result <- f(serie, m = m, lag_max = lag_max, est = "n", plot = FALSE)
        expected_RHO <- build_RHO(detrended_serie, m = m, lag_max = lag_max, est = "n")
        expected_rho <- build_rho(detrended_serie, m = m, lag_max = lag_max, est = "n")
        expect_equal(result$RHO, expected_RHO)
        expect_equal(result$rho, expected_rho)
    })
})

test_that("build_RHO_A augments build_RHO with the annual border (A last)", {
    serie <- setup_test_data()
    ss <- scale_by_season(serie, est = "n")[[1]]
    sa <- scale_by_season(calcula_medias_anuais(serie), est = "n")[[1]]
    lag_max <- 4
    m <- 3
    RA <- build_RHO_A(ss, sa, m, lag_max, "n")
    expect_equal(dim(RA), c(lag_max + 1, lag_max + 1))
    expect_equal(RA[seq_len(lag_max), seq_len(lag_max)], build_RHO(ss, m, lag_max, "n"))
    expect_equal(RA, t(RA))
    expect_equal(RA[lag_max + 1, lag_max + 1], 1)
    expect_true(all(is.finite(RA)))
})

test_that("build_rho_A appends rho_A0 to build_rho (A last)", {
    serie <- setup_test_data()
    ss <- scale_by_season(serie, est = "n")[[1]]
    sa <- scale_by_season(calcula_medias_anuais(serie), est = "n")[[1]]
    lag_max <- 4
    m <- 3
    rA <- build_rho_A(ss, sa, m, lag_max, "n")
    expect_equal(length(rA), lag_max + 1)
    expect_equal(rA[seq_len(lag_max)], build_rho(ss, m, lag_max, "n"))
    expect_true(all(is.finite(rA)))
})

test_that("percacf returns a full conditional FACP (first lag filled)", {
    serie <- setup_test_data()
    lag_max <- 6
    m <- 3
    cc <- percacf(serie, m, lag_max = lag_max, est = "n")
    expect_equal(length(cc$phi), lag_max)
    expect_true(all(is.finite(cc$phi)))
    expect_true(cc$phi[1] != 0)
})

test_that("percacf conditional FACP differs from the ordinary perpacf", {
    serie <- setup_test_data()
    lag_max <- 6
    m <- 3
    cc <- percacf(serie, m, lag_max = lag_max, est = "n")
    pp <- perpacf(scale_by_season(serie, est = "n")[[1]], m, lag_max)
    expect_false(isTRUE(all.equal(cc$phi, pp$phi)))
})

test_that("percacf reproduces CEPEL/GEVAZP CORRELOGRAMA PARCIAL for the Furnas series", {
    # Golden master: dummyseries[, 2] is posto FURNAS, fitted by GEVAZP (see exemplo_furnas). The
    # conditional FACP is the PARTIAL CORRELATION (normalized, in [-1, 1]) conditioned on the annual
    # component -- NOT the augmented-system regression coefficient. Guards both the cross-correlation
    # divisor convention and the partial-correlation read-out. Rows are order 1..11, cols JAN..DEZ.
    pacf_cep <- matrix(c(
        0.47769, 0.43127, 0.49418, 0.59960, 0.63464, 0.58379, 0.74643, 0.77672, 0.58863, 0.51179, 0.55649, 0.49778,
       -0.06997,-0.09241, 0.13286, 0.30904, 0.29508, 0.13376, 0.51280,-0.16600, 0.32834, 0.07341, 0.18269, 0.21787,
       -0.10078,-0.02237,-0.03768,-0.00821, 0.29356, 0.01200, 0.10347, 0.07870, 0.24329, 0.27765,-0.10032, 0.02612,
       -0.13525,-0.08909,-0.03760, 0.00657, 0.00207,-0.01420, 0.03176,-0.15397,-0.15368, 0.23969,-0.05275, 0.16634,
        0.05361,-0.14526, 0.07487,-0.03822, 0.09650, 0.16497, 0.03152, 0.07305, 0.11883,-0.07372,-0.13487, 0.02752,
        0.04041, 0.12842,-0.19397,-0.08391,-0.12012, 0.11100,-0.00869,-0.20050,-0.01531,-0.05186,-0.06375, 0.01933,
       -0.14645, 0.00976, 0.20850, 0.12197, 0.05768,-0.05545, 0.05416,-0.18588,-0.05683,-0.09575,-0.17798,-0.08879,
        0.14455,-0.04416, 0.05902, 0.29880, 0.04368,-0.20706, 0.08258, 0.01955,-0.13300, 0.09851, 0.05191, 0.00174,
        0.01074, 0.31252,-0.09573,-0.05090, 0.03106,-0.02804,-0.09686, 0.19126,-0.05571,-0.00835, 0.02769,-0.10209,
        0.12798,-0.03317, 0.07369,-0.16472, 0.18027,-0.11515,-0.25510,-0.05296, 0.06218, 0.10293, 0.05216,-0.08906,
       -0.04272, 0.09933, 0.02784,-0.24392, 0.07984, 0.12038,-0.14512,-0.04072, 0.19310,-0.07111,-0.01978, 0.02094),
        nrow = 11, byrow = TRUE)

    serie <- setup_test_data()
    for (m in seq_len(12)) {
        phi <- percacf(serie, m, lag_max = 11, est = "n")$phi
        expect_lt(max(abs(phi - pacf_cep[, m])), 1e-3)
    }
})
