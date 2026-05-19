#' PACF Periodica De Series Temporais
#' 
#' Calcula a funcao de autocorrelacao parcial periodica (PACF) de uma serie temporal sazonal
#' 
#' @param serie serie temporal sazonal
#' @param m estacao do ano para a qual calcular a PACF
#' @param lag_max numero maximo de lags a considerar. Default é frequencia da serie - 1
#' @param est um de `c("n", "n-1")` indicando qual denominador utilizar para desvio padrao
#' @param plot se TRUE, plota a PACF periodica
#'
#' @return objeto da classe `periodic_pacf`, uma lista contendo os seguintes elementos:
#'     * `phi`: vetor com os valores da PACF periodica
#'     * `n.used`: numero de observacoes utilizadas no calculo
#'     * `m`: estacao do ano para a qual a PACF foi calculada
#'     * `lag_max`: numero maximo de lags considerados
#'     * `rho`: vetor com os valores da funcao de autocorrelacao (ACF) periodica
#'     * `RHO`: matriz com os valores da matriz de autocorrelacao (ACF) periodica
#' 
#' @export

perpacf <- function(serie, m, lag_max = frequency(serie) - 1, est = c("n", "n-1"), plot = FALSE) {
    est <- match.arg(est)

    RHO <- build_RHO(serie, m, lag_max, est)
    rho <- build_rho(serie, m, lag_max, est)

    phi <- double(lag_max)
    for (i in seq_len(lag_max)) {
        phi[i] <- solve(RHO[1:i, 1:i], rho[1:i])[i]
    }

    perpacf <- list(phi = phi, n_used = floor(length(serie) / frequency(serie)),
        m = m, lag_max = lag_max, rho = rho, RHO = RHO)
    class(perpacf) <- c("periodic_pacf", "list")

    if (plot) print(plot(perpacf))

    return(perpacf)
}

#' Funcoes Auxiliares De `perpacf`
#' 
#' Funcoes internas para calculo da funcao de autocorrelacao periodica
#' 
#' Estas funcoes calculam o vetor `rho` e a matriz `RHO` necessarios para execucao do metodo de
#' Yule-Walker para calculo de autocorrelacoes parciais periodicas.
#' 
#' @param serie serie temporal sazonal
#' @param m estacao do ano para a qual calcular a PACF
#' @param lag_max numero maximo de lags a considerar. Default é frequencia da serie - 1
#' @param est um de `c("n", "n-1")` indicando qual denominador utilizar para desvio padrao
#' 
#' @rdname perpacf_helpers

build_RHO <- function(serie, m, lag_max, est) {
    serie <- ts2matrix(serie)
    N     <- nrow(serie)
    s     <- ncol(serie)
    fsd   <- ifelse(est == "n-1", sd, sd2)

    RHO <- diag(1, lag_max, lag_max)

    for (i in seq_len(lag_max)) {

        # Identifica a coluna correspondente ao lag i do mes m
        col1 <- wrap_season(m - i, s)

        for (j in seq_len(lag_max)[-seq_len(i)]) {

            # Identifica a serie lag
            col2 <- wrap_season(col1 - (j - i), s)

            if (col1 < col2) {
                vec1 <- serie[2:N, col1]
                vec2 <- serie[1:(N - 1), col2]
                RHO[i, j] <- 1 / N * sum(vec1 * vec2, na.rm = TRUE) /
                    (fsd(serie[, col1], na.rm = TRUE) * fsd(serie[, col2], na.rm = TRUE))
                RHO[j, i] <- RHO[i, j]
            } else {
                vec1 <- serie[, col1]
                vec2 <- serie[, col2]
                RHO[i, j] <- mean(vec1 * vec2, na.rm = TRUE) /
                    (fsd(vec1, na.rm = TRUE) * fsd(vec2, na.rm = TRUE))
                RHO[j, i] <- RHO[i, j]
            }
        }
    }

    return(RHO)
}

#' @rdname perpacf_helpers

build_rho <- function(serie, m, lag_max, est) {
    serie <- ts2matrix(serie)
    N     <- nrow(serie)
    s     <- ncol(serie)
    fsd   <- ifelse(est == "n-1", sd, sd2)

    rho <- rep(NA_real_, lag_max)

    for (i in seq_len(lag_max)) {
        col2 <- wrap_season(m - i, s)

        if (m < col2) {
            vec1 <- serie[2:N, m]
            vec2 <- serie[1:(N - 1), col2]
            rho[i] <- 1 / N * sum(vec1 * vec2, na.rm = TRUE) /
                (fsd(serie[, m], na.rm = TRUE) * fsd(serie[, col2], na.rm = TRUE))
        } else {
            vec1 <- serie[, m]
            vec2 <- serie[, col2]
            rho[i] <- mean(vec1 * vec2, na.rm = TRUE) /
                (fsd(vec1, na.rm = TRUE) * fsd(vec2, na.rm = TRUE))
        }
    }

    return(rho)
}