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

test_that("fit_par_a_yulewalker reproduces CEPEL/GEVAZP coefficients for the Furnas series", {
    # Golden master: dummyseries[, 2] is the incremental inflow of posto FURNAS, the exact series
    # fitted by CEPEL's GEVAZP (see exemplo_furnas). With CEPEL's per-month AR orders, the augmented
    # Yule-Walker fit must reproduce the published phis and the "PARTE A12" psi to printed precision.
    # Guards the cross-correlation divisor convention (sum / N, not mean over valid pairs).
    serie   <- setup_test_data()
    serie_A <- calcula_medias_anuais(serie)

    cep_p   <- c(1, 1, 1, 8, 3, 8, 2, 1, 3, 4, 1, 2)
    cep_psi <- c(-0.0472, 0.0279, 0.155, -0.318, 0.214, 0.290,
                  0.0363, 0.118, -0.203, -0.0343, 0.111, 0.140)
    cep_phi <- list(
        0.593, 0.490, 0.520,
        c(0.557, 0.361, 0.132, 0.0633, 0.0576, -0.0794, 0.0104, 0.375),
        c(0.395, 0.200, 0.202),
        c(0.562, 0.107, -0.00468, -0.0775, 0.0344, 0.106, -0.0174, -0.216),
        c(0.523, 0.426), 0.823, c(0.433, 0.255, 0.342),
        c(0.348, -0.0998, 0.245, 0.372), 0.639, c(0.404, 0.257))

    # CEPEL prints coefficients to ~3 significant figures, so compare on absolute error (relative
    # tolerance is meaningless for near-zero coefficients such as -0.00468).
    for (m in seq_len(12)) {
        fit <- fit_par_a_yulewalker(serie, serie_A, m, cep_p[m], max_p = 11)
        err <- max(abs(c(fit$phis - cep_phi[[m]], fit$psi - cep_psi[m])))
        expect_lt(err, 1e-3)
    }
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
