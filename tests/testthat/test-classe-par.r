setup_test_data <- function() {
    data("dummyseries", envir = environment())
    serie <- ts(dummyseries[, 2], start = c(1931, 1), frequency = 12)
    return(serie)
}

test_that("new_par", {
    f <- new_par
    expect_true(is.function(f))
    test_that("new_par creates an object of class par", {
        serie <- setup_test_data()
        phis <- list(c(0.5), c(0.3))
        result <- f(x = serie, phis = phis, sigma2 = 1.5, sigma2_norm = 1.2, residuals = NULL, call = quote(par()))
        expect_true(inherits(result, "par"))
        expect_true(inherits(result, "list"))
        expect_equal(class(result), c("par", "list"))
    })
    test_that("new_par has all required elements", {
        serie <- setup_test_data()
        phis <- list(c(0.5), c(0.3))
        result <- f(x = serie, phis = phis, sigma2 = 1.5, sigma2_norm = 1.2, residuals = NULL, call = quote(par()))
        expect_true("phis" %in% names(result))
        expect_true("sigma2" %in% names(result))
        expect_true("sigma2_norm" %in% names(result))
        expect_true("x" %in% names(result))
        expect_true("residuals" %in% names(result))
        expect_true("call" %in% names(result))
    })
    test_that("new_par elements are assigned correctly", {
        serie <- setup_test_data()
        phis <- list(c(0.5, 0.2), c(0.3))
        sigma2_val <- 1.5
        sigma2_norm_val <- 1.2
        call_val <- quote(par())
        result <- f(x = serie, phis = phis, sigma2 = sigma2_val, sigma2_norm = sigma2_norm_val,
            residuals = NULL, call = call_val)
        expect_equal(result$phis, phis)
        expect_equal(result$sigma2, sigma2_val)
        expect_equal(result$sigma2_norm, sigma2_norm_val)
        expect_equal(result$x, serie)
        expect_null(result$residuals)
        expect_equal(result$call, call_val)
    })
    test_that("new_par with NULL residuals", {
        serie <- setup_test_data()
        phis <- list(c(0.5))
        result <- f(x = serie, phis = phis, sigma2 = 1.5, sigma2_norm = 1.2, residuals = NULL, call = quote(par()))
        expect_null(result$residuals)
    })
    test_that("new_par with actual residuals", {
        serie <- setup_test_data()
        phis <- list(c(0.5))
        resid <- ts(rnorm(length(serie)), start = start(serie), frequency = frequency(serie))
        result <- f(x = serie, phis = phis, sigma2 = 1.5, sigma2_norm = 1.2, residuals = resid, call = quote(par()))
        expect_equal(result$residuals, resid)
    })
})

test_that("idordem_yulewalker", {
    f <- idordem_yulewalker
    expect_true(is.function(f))
    test_that("idordem_yulewalker returns a single integer value", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, max_p = 5)
        expect_true(is.numeric(result))
        expect_equal(length(result), 1)
        expect_true(result == floor(result))
    })
    test_that("idordem_yulewalker returned order is >= 1", {
        serie <- setup_test_data()
        result <- f(serie, m = 1, max_p = 5)
        expect_true(result >= 1)
    })
    test_that("idordem_yulewalker returned order is <= max_p", {
        serie <- setup_test_data()
        max_p <- 5
        result <- f(serie, m = 1, max_p = max_p)
        expect_true(result <= max_p)
    })
    test_that("idordem_yulewalker works with different seasons", {
        serie <- setup_test_data()
        result_m1 <- f(serie, m = 1, max_p = 5)
        result_m6 <- f(serie, m = 6, max_p = 5)
        result_m12 <- f(serie, m = 12, max_p = 5)
        expect_true(is.numeric(result_m1))
        expect_true(is.numeric(result_m6))
        expect_true(is.numeric(result_m12))
    })
    test_that("idordem_yulewalker works with different max_p values", {
        serie <- setup_test_data()
        result_3 <- f(serie, m = 1, max_p = 3)
        result_8 <- f(serie, m = 1, max_p = 8)
        expect_true(result_3 <= 3)
        expect_true(result_8 <= 8)
    })
    test_that("idordem_yulewalker is deterministic", {
        serie <- setup_test_data()
        result1 <- f(serie, m = 1, max_p = 5)
        result2 <- f(serie, m = 1, max_p = 5)
        expect_equal(result1, result2)
    })
})

test_that("par", {
    f <- par
    expect_true(is.function(f))
    test_that("par returns an object of class par", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(inherits(result, "par"))
    })
    test_that("par with ps='auto' automatic order selection", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1980, 1), end = c(2000, 12))
        result <- f(serie_subset, ps = "auto", max_ps = 3, metodo = "YuleWalker")
        expect_true(inherits(result, "par"))
        expect_true(is.list(result$phis))
    })
    test_that("par with ps as scalar gets repeated for all seasons", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_equal(length(result$phis), frequency(serie_subset))
        expect_true(all(sapply(result$phis, length) == 2))
    })
    test_that("par with ps as vector", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        ps_vec <- rep(2, 12)
        result <- f(serie_subset, ps = ps_vec, metodo = "YuleWalker")
        expect_equal(length(result$phis), 12)
    })
    test_that("par throws error for non-seasonal series", {
        serie_no_season <- ts(rnorm(100), frequency = 1)
        expect_error(f(serie_no_season, ps = 2), "'serie' nao possui sazonalidade")
    })
    test_that("par phis is a list with length equal to frequency", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(is.list(result$phis))
        expect_equal(length(result$phis), frequency(serie_subset))
    })
    test_that("par each element of phis is a numeric vector", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(all(sapply(result$phis, is.numeric)))
        expect_true(all(sapply(result$phis, function(x) is.vector(x) || is.null(x))))
    })
    test_that("par sigma2 is a single numeric value > 0", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(is.numeric(result$sigma2))
        expect_equal(length(result$sigma2), 1)
        expect_true(result$sigma2 > 0)
    })
    test_that("par residuals are calculated and stored", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(!is.null(result$residuals))
        expect_true(is.ts(result$residuals))
    })
    test_that("par x original series is stored", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_equal(result$x, serie_subset)
    })
    test_that("par call is stored", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(!is.null(result$call))
    })
    test_that("par length of residuals equals length of x", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        result <- f(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_equal(length(result$residuals), length(result$x))
    })
})

test_that("fitted.par", {
    f <- fitted.par
    expect_true(is.function(f))
    test_that("fitted.par returns a time series", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result <- f(model)
        expect_true(is.ts(result))
    })
    test_that("fitted.par length equals length of original series", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result <- f(model)
        expect_equal(length(result), length(model$x))
    })
    test_that("fitted.par fitted plus residuals approximately equals original series", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        fitted_vals <- f(model)
        resids <- residuals(model)
        reconstructed <- fitted_vals + resids
        non_na_idx <- !is.na(reconstructed)
        expect_equal(as.numeric(reconstructed[non_na_idx]), as.numeric(model$x[non_na_idx]), tolerance = 1e-10)
    })
    test_that("fitted.par uses cached residuals when available", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(!is.null(model$residuals))
        result <- f(model)
        expect_true(is.ts(result))
    })
})

test_that("fitted.values.par", {
    f <- fitted.values.par
    expect_true(is.function(f))
    test_that("fitted.values.par is an alias for fitted.par", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result1 <- fitted.par(model)
        result2 <- f(model)
        expect_equal(result1, result2)
    })
})

test_that("residuals.par", {
    f <- residuals.par
    expect_true(is.function(f))
    test_that("residuals.par returns a time series", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result <- f(model)
        expect_true(is.ts(result))
    })
    test_that("residuals.par length equals length of original series", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result <- f(model)
        expect_equal(length(result), length(model$x))
    })
    test_that("residuals.par uses cached residuals when available", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        expect_true(!is.null(model$residuals))
        result <- f(model)
        expect_equal(result, model$residuals)
    })
})

test_that("predict.par", {
    f <- predict.par
    expect_true(is.function(f))
    test_that("predict.par returns a time series", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result <- f(model, n.ahead = 12)
        expect_true(is.ts(result))
    })
    test_that("predict.par length equals n.ahead", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        n_ahead <- 12
        result <- f(model, n.ahead = n_ahead)
        expect_equal(length(result), n_ahead)
    })
    test_that("predict.par start time is correct", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result <- f(model, n.ahead = 12)
        expected_start <- c(1996, 1)
        expect_equal(start(result), expected_start)
    })
    test_that("predict.par frequency is preserved", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result <- f(model, n.ahead = 12)
        expect_equal(frequency(result), frequency(model$x))
    })
    test_that("predict.par works with different n.ahead values", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1990, 1), end = c(1995, 12))
        model <- par(serie_subset, ps = 2, metodo = "YuleWalker")
        result_1 <- f(model, n.ahead = 1)
        result_6 <- f(model, n.ahead = 6)
        result_12 <- f(model, n.ahead = 12)
        result_24 <- f(model, n.ahead = 24)
        expect_equal(length(result_1), 1)
        expect_equal(length(result_6), 6)
        expect_equal(length(result_12), 12)
        expect_equal(length(result_24), 24)
    })
})
