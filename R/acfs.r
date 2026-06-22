# PERIODIC PARTIAL ACF -----------------------------------------------------------------------------

#' PACF Periodica De Series Temporais
#' 
#' Calcula a funcao de autocorrelacao parcial periodica (PACF) de uma serie temporal sazonal
#' 
#' @param serie serie temporal sazonal. Sera padronizada internamente
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

    if (!inherits(serie, "detrended_ts")) serie <- scale_by_season(serie)[[1]]

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

# CONDITIONAL PERIODIC ACF -------------------------------------------------------------------------

percacf <- function(serie, m, serie_anual = NULL, lag_max = frequency(serie) - 1,
    est = c("n", "n-1"), plot = FALSE) {

    est <- match.arg(est)

    if (is.null(serie_anual)) serie_anual <- calcula_medias_anuais(serie)
    serie_anual <- scale_by_season(serie_anual)[[1]]

    if (!inherits(serie, "detrended_ts")) serie <- scale_by_season(serie)[[1]]

    RHO <- build_RHO_A(serie, serie_anual, m, lag_max, est)
    rho <- build_rho_A(serie, serie_anual, m, lag_max, est)

    # FACP condicional: em cada lag k, ajusta o sistema com os lags 1..k MAIS o termo anual (ultima
    # linha/coluna) e le o coeficiente do lag k. phi[k] e a autocorrelacao parcial do lag k
    # condicionada nos lags anteriores e na componente anual.
    a   <- lag_max + 1
    phi <- double(lag_max)
    for (k in seq_len(lag_max)) {
        idx    <- c(seq_len(k), a)
        sol    <- solve(RHO[idx, idx], rho[idx])
        phi[k] <- sol[k]
    }

    percacf <- list(phi = phi, n_used = floor(length(serie) / frequency(serie)),
        m = m, lag_max = lag_max, rho = rho, RHO = RHO)
    class(percacf) <- c("periodic_cacf", "list")

    if (plot) print(plot(percacf))

    return(percacf)
}

#' Funcoes Auxiliares De `percacf`
#'
#' Funcoes internas para montagem do sistema de Yule-Walker aumentado (com componente anual)
#'
#' Reaproveitam `build_RHO`/`build_rho` para o bloco Z-Z e apenas acrescentam a borda anual: a
#' componente anual entra como ultima linha/coluna do sistema (layout CEPEL / `methodology_report` 4).
#'
#' @param serie serie temporal sazonal (ja padronizada)
#' @param serie_anual componente anual (ja padronizada) alinhada a `serie`
#' @param m estacao do ano para a qual calcular
#' @param lag_max numero maximo de lags a considerar
#' @param est um de `c("n", "n-1")` indicando qual denominador utilizar para desvio padrao
#'
#' @rdname percacf_helpers

build_RHO_A <- function(serie, serie_anual, m, lag_max, est) {
    RHO_zz      <- build_RHO(serie, m, lag_max, est)
    serie       <- ts2matrix(serie)
    serie_anual <- ts2matrix(serie_anual)
    N     <- nrow(serie)
    s     <- ncol(serie)
    fsd   <- ifelse(est == "n-1", sd, sd2)

    a   <- lag_max + 1
    RHO <- diag(1, a, a)
    RHO[seq_len(lag_max), seq_len(lag_max)] <- RHO_zz

    # Borda anual (ultima linha/coluna): rho_{A,k} = corr(A no mes m, Z no lag k), k = 1..lag_max.
    # Mesmo pareamento de build_rho, trocando vec1 pela serie anual no mes m.
    for (k in seq_len(lag_max)) {
        col2 <- wrap_season(m - k, s)

        if (m < col2) {
            vec1 <- serie_anual[2:N, m]
            vec2 <- serie[1:(N - 1), col2]
            RHO[a, k] <- 1 / N * sum(vec1 * vec2, na.rm = TRUE) /
                (fsd(serie_anual[, m], na.rm = TRUE) * fsd(serie[, col2], na.rm = TRUE))
        } else {
            vec1 <- serie_anual[, m]
            vec2 <- serie[, col2]
            RHO[a, k] <- mean(vec1 * vec2, na.rm = TRUE) /
                (fsd(vec1, na.rm = TRUE) * fsd(vec2, na.rm = TRUE))
        }
        RHO[k, a] <- RHO[a, k]
    }
    # canto rho_AA = 1 ja vem de diag()

    return(RHO)
}

#' @rdname percacf_helpers

build_rho_A <- function(serie, serie_anual, m, lag_max, est) {
    rho_zz      <- build_rho(serie, m, lag_max, est)
    serie       <- ts2matrix(serie)
    serie_anual <- ts2matrix(serie_anual)
    fsd         <- ifelse(est == "n-1", sd, sd2)

    # rho_{A,0} = corr(A no mes m, Z no mes m) na mesma posicao (lag 0); ultima entrada do RHS.
    rho_A0 <- mean(serie_anual[, m] * serie[, m], na.rm = TRUE) /
        (fsd(serie_anual[, m], na.rm = TRUE) * fsd(serie[, m], na.rm = TRUE))

    return(c(rho_zz, rho_A0))
}
