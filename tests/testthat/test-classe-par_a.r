setup_test_data <- function() {
    data("dummyseries", envir = environment())
    serie <- ts(dummyseries[, 2], start = c(1931, 1), frequency = 12)
    return(serie)
}

test_that("par_a returns a par-inheriting object with frequency-wide coefficients", {
    serie <- setup_test_data()
    mod <- par_a(serie)
    expect_equal(class(mod), c("par_a", "par", "list"))
    expect_true(inherits(mod, "par"))
    expect_equal(length(mod$phis), frequency(serie))
    expect_true(all(vapply(mod$phis, length, integer(1)) == frequency(serie)))
})

test_that("par_a stores psi, ps and the annual component", {
    serie <- setup_test_data()
    mod <- par_a(serie)
    s <- frequency(serie)
    expect_equal(length(mod$psi), s)
    expect_equal(length(mod$ps), s)
    expect_true(is.ts(mod$annual))
    expect_equal(length(mod$annual), length(serie))
    expect_equal(mod$annual, calcula_medias_anuais(serie))
})

test_that("par_a methods dispatch by inheritance and reconstruct the series for t > s", {
    serie <- setup_test_data()
    mod <- par_a(serie)
    s <- frequency(serie)
    recon <- fitted(mod) + residuals(mod)
    idx <- (s + 1):length(serie)
    expect_equal(as.numeric(recon)[idx], as.numeric(serie)[idx])
    expect_equal(length(predict(mod, n.ahead = s)), s)
})

test_that("fit_par_a_yulewalker returns phis of the requested order plus a scalar psi", {
    serie <- setup_test_data()
    serie_A <- calcula_medias_anuais(serie)
    fit <- fit_par_a_yulewalker(serie, serie_A, m = 3, p = 2, max_p = 11)
    expect_equal(length(fit$phis), 2)
    expect_true(is.numeric(fit$psi) && length(fit$psi) == 1)
})

test_that("idordem_cacf identifies order from the conditional FACP within bounds", {
    serie <- setup_test_data()
    serie_A <- calcula_medias_anuais(serie)
    ordem <- idordem_cacf(serie, serie_A, m = 3, max_p = 11)
    expect_true(ordem >= 1 && ordem <= 11)
    expect_equal(ordem, as.integer(ordem))
})

test_that("unfold_par_a matches the closed-form annual stride", {
    sds  <- rep(1, 12)
    sdsA <- rep(1, 12)
    phis <- c(0.5, 0.2)
    psi  <- 0.3
    out  <- unfold_par_a(phis, psi, sds, sdsA, m = 4, s = 12)
    expect_equal(length(out), 12)
    # l = 1 (<= p): phi_1 + psi / (12 * sdsA[m-1]) * sds[m-1]
    expect_equal(out[1], 0.5 + 0.3 / 12)
    # l = 3 (> p): 0 + psi / (12 * sdsA[m-1]) * sds[m-3]
    expect_equal(out[3], 0.3 / 12)
})
