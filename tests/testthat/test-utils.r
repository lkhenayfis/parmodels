setup_test_data <- function() {
    data("dummyseries", envir = environment())
    serie <- ts(dummyseries[, 2], start = c(1931, 1), frequency = 12)
    return(serie)
}

test_that("sd2", {
    f <- sd2
    expect_true(is.function(f))
    test_that("sd2 calculates standard deviation with n denominator", {
        x <- c(1, 2, 3, 4, 5)
        n <- length(x)
        result <- f(x)
        expected <- sqrt((n - 1) / n) * sd(x)
        expect_equal(result, expected)
    })
    test_that("sd2 handles NA values with na.rm = TRUE", {
        x <- c(1, 2, NA, 4, 5)
        n_valid <- sum(!is.na(x))
        result <- f(x, na.rm = TRUE)
        expected <- sqrt((n_valid - 1) / n_valid) * sd(x, na.rm = TRUE)
        expect_equal(result, expected)
    })
    test_that("sd2 returns NA when vector contains NA and na.rm = FALSE", {
        x <- c(1, 2, NA, 4, 5)
        result <- f(x, na.rm = FALSE)
        expect_true(is.na(result))
    })
    test_that("sd2 works with dummyseries data", {
        serie <- setup_test_data()
        result <- f(serie)
        expect_true(is.numeric(result))
        expect_true(length(result) == 1)
        expect_true(result > 0)
    })
})

test_that("cov2", {
    f <- cov2
    expect_true(is.function(f))
    test_that("cov2 calculates covariance with n denominator", {
        x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), ncol = 2)
        n <- nrow(x)
        result <- f(x)
        expected <- (n - 1) / n * cov(x)
        expect_equal(result, expected)
    })
    test_that("cov2 handles NA values with na.rm = TRUE", {
        x <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8, 9), ncol = 3)
        n <- nrow(x)
        result <- f(x, na.rm = TRUE)
        expected <- (n - 1) / n * cov(x, use = "pairwise.complete.obs")
        expect_equal(result, expected)
    })
    test_that("cov2 handles NA values with na.rm = FALSE", {
        x <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8, 9), ncol = 3)
        result <- f(x, na.rm = FALSE)
        expect_true(any(is.na(result)))
    })
    test_that("cov2 works with time series matrix", {
        serie <- setup_test_data()
        x <- matrix(serie[1:24], ncol = 2)
        result <- f(x)
        expect_true(is.matrix(result))
        expect_equal(dim(result), c(2, 2))
    })
})

test_that("split_by_season", {
    f <- split_by_season
    expect_true(is.function(f))
    test_that("split_by_season splits time series by season correctly", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        expect_true(is.list(result))
        expect_equal(length(result), 12)
        expect_true(!is.null(attr(result, "orig_order")))
    })
    test_that("split_by_season preserves all observations", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        total_obs <- sum(sapply(result, length))
        expect_equal(total_obs, length(serie_subset))
    })
    test_that("split_by_season orig_order attribute reconstructs original order", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        orig_order <- attr(result, "orig_order")
        reconstructed <- unlist(result)[orig_order]
        expect_equal(as.numeric(reconstructed), as.numeric(serie_subset))
    })
    test_that("split_by_season groups same seasons together", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        expect_equal(length(result[[1]]), 3)
    })
})

test_that("seasonal_mean", {
    f <- seasonal_mean
    expect_true(is.function(f))
    test_that("seasonal_mean returns vector with length equal to frequency", {
        serie <- setup_test_data()
        result <- f(serie)
        expect_equal(length(result), frequency(serie))
    })
    test_that("seasonal_mean calculates correct means for each season", {
        serie <- setup_test_data()
        result <- f(serie)
        jan_values <- serie[cycle(serie) == 1]
        expected_jan_mean <- mean(jan_values, na.rm = TRUE)
        expect_equal(result[1], expected_jan_mean)
    })
    test_that("seasonal_mean works with partial years", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 3), end = c(1933, 10))
        result <- f(serie_subset)
        expect_equal(length(result), 12)
        expect_true(all(is.finite(result)))
    })
    test_that("seasonal_mean handles NA values", {
        serie <- setup_test_data()
        serie_with_na <- serie
        serie_with_na[1:5] <- NA
        result <- f(serie_with_na)
        expect_equal(length(result), 12)
        expect_true(all(is.finite(result)))
    })
})

test_that("seasonal_sd", {
    f <- seasonal_sd
    expect_true(is.function(f))
    test_that("seasonal_sd returns vector with length equal to frequency", {
        serie <- setup_test_data()
        result <- f(serie, est = "n")
        expect_equal(length(result), frequency(serie))
    })
    test_that("seasonal_sd with est='n' uses sd2", {
        serie <- setup_test_data()
        result_n <- f(serie, est = "n")
        jan_values <- serie[cycle(serie) == 1]
        expected_jan_sd <- sd2(jan_values, na.rm = TRUE)
        expect_equal(result_n[1], expected_jan_sd)
    })
    test_that("seasonal_sd with est='n-1' uses standard sd", {
        serie <- setup_test_data()
        result_n1 <- f(serie, est = "n-1")
        jan_values <- serie[cycle(serie) == 1]
        expected_jan_sd <- sd(jan_values, na.rm = TRUE)
        expect_equal(result_n1[1], expected_jan_sd)
    })
    test_that("seasonal_sd results differ between est='n' and est='n-1'", {
        serie <- setup_test_data()
        result_n <- f(serie, est = "n")
        result_n1 <- f(serie, est = "n-1")
        expect_false(identical(result_n, result_n1))
        expect_true(all(result_n < result_n1))
    })
})

test_that("scale_by_season", {
    f <- scale_by_season
    expect_true(is.function(f))
    test_that("scale_by_season returns list with 2 elements", {
        serie <- setup_test_data()
        result <- f(serie)
        expect_true(is.list(result))
        expect_equal(length(result), 2)
    })
    test_that("scale_by_season first element is scaled series with same attributes", {
        serie <- setup_test_data()
        result <- f(serie)
        scaled_serie <- result[[1]]
        expect_true(is.ts(scaled_serie))
        expect_equal(length(scaled_serie), length(serie))
        expect_equal(tsp(scaled_serie), tsp(serie))
        expect_equal(frequency(scaled_serie), frequency(serie))
    })
    test_that("scale_by_season second element contains means and sds", {
        serie <- setup_test_data()
        result <- f(serie)
        scales <- result[[2]]
        expect_true(is.list(scales))
        expect_equal(length(scales), 2)
        expect_equal(length(scales[[1]]), 12)
        expect_equal(length(scales[[2]]), 12)
    })
    test_that("scale_by_season scaling is correct", {
        serie <- setup_test_data()
        result <- f(serie)
        scaled <- result[[1]]
        means <- result[[2]][[1]]
        sds <- result[[2]][[2]]
        first_season_idx <- which(cycle(serie) == 1)[1]
        expected_scaled <- (serie[first_season_idx] - means[1]) / sds[1]
        expect_equal(as.numeric(scaled[first_season_idx]), as.numeric(expected_scaled))
    })
    test_that("scale_by_season accepts provided means and sds", {
        serie <- setup_test_data()
        means <- seasonal_mean(serie)
        sds <- seasonal_sd(serie, est = "n")
        result1 <- f(serie)
        result2 <- f(serie, means = means, sds = sds)
        expect_equal(result1[[1]], result2[[1]])
    })
    test_that("scale_by_season respects est parameter", {
        serie <- setup_test_data()
        result_n <- f(serie, est = "n")
        result_n1 <- f(serie, est = "n-1")
        expect_false(identical(result_n[[1]], result_n1[[1]]))
    })
    test_that("scale_by_season is reversible", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1935, 12))
        result <- f(serie_subset)
        scaled <- result[[1]]
        means <- result[[2]][[1]]
        sds <- result[[2]][[2]]
        unscaled <- scaled
        for (i in seq_along(scaled)) {
            m <- cycle(scaled)[i]
            unscaled[i] <- scaled[i] * sds[m] + means[m]
        }
        expect_equal(as.numeric(unscaled), as.numeric(serie_subset), tolerance = 1e-10)
    })
})

test_that("pad_series", {
    f <- pad_series
    expect_true(is.function(f))
    test_that("pad_series adds padding to complete seasonal cycles", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 5), end = c(1933, 8))
        result <- f(serie_subset)
        expect_equal(length(result), length(serie_subset) + 4 + 4)
    })
    test_that("pad_series pads with NA by default", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 5), end = c(1933, 8))
        result <- f(serie_subset)
        expect_true(all(is.na(result[1:4])))
        expect_true(all(is.na(tail(result, 4))))
    })
    test_that("pad_series accepts custom pad value", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 5), end = c(1933, 8))
        result <- f(serie_subset, pad = 0)
        expect_equal(result[1:4], rep(0, 4))
        expect_equal(tail(result, 4), rep(0, 4))
    })
    test_that("pad_series preserves frequency", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 5), end = c(1933, 8))
        result <- f(serie_subset)
        expect_equal(frequency(result), frequency(serie_subset))
    })
    test_that("pad_series handles series already starting at season 1", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        expect_equal(length(result), length(serie_subset))
    })
    test_that("pad_series adjusts start time correctly", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 5), end = c(1933, 8))
        result <- f(serie_subset)
        expect_equal(start(result), c(1931, 1))
    })
})

test_that("ts2matrix", {
    f <- ts2matrix
    expect_true(is.function(f))
    test_that("ts2matrix returns a matrix", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        expect_true(is.matrix(result))
    })
    test_that("ts2matrix number of columns equals frequency", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        expect_equal(ncol(result), frequency(serie))
    })
    test_that("ts2matrix pads series to complete cycles", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 5), end = c(1933, 8))
        result <- f(serie_subset)
        expect_true(any(is.na(result)))
    })
    test_that("ts2matrix organizes by seasonal cycles", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1933, 12))
        result <- f(serie_subset)
        expect_equal(nrow(result), 3)
        expect_equal(result[1, ], as.numeric(window(serie_subset, start = c(1931, 1), end = c(1931, 12))))
    })
    test_that("ts2matrix preserves values correctly", {
        serie <- setup_test_data()
        serie_subset <- window(serie, start = c(1931, 1), end = c(1932, 12))
        result <- f(serie_subset)
        result_vec <- as.vector(t(result))
        result_vec <- result_vec[!is.na(result_vec)]
        expect_equal(result_vec, as.numeric(serie_subset))
    })
})