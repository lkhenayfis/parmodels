# UTILS GERAIS -------------------------------------------------------------------------------------

#' Calculo do desvio padrao normalizado por 1/n
#'
#' Funcao alternativa para calcular desvio padrao atraves do estimador "n"
#' 
#' @param x vetor numerico do qual se deseja calcular o desvio padrao
#' @param na.rm logico indicando se valores NA devem ser removidos antes do calculo
#' 
#' @return desvio padrao calculado

sd2 <- function(x, na.rm = FALSE) {
    n <- length(x)
    sqrt((n - 1) / n) * sd(x, na.rm)
}

#' Calculo da covariancia normalizada por 1/n
#'
#' Funcao alternativa para calcular covariancia atraves do estimador "n"
#' 
#' @param x matriz numerica na qual se deseja calcular a covariancia
#' @param na.rm logico indicando se valores NA devem ser removidos antes do calculo
#' 
#' @return covariancia calculada

cov2 <- function(x, na.rm = FALSE) {

    n <- nrow(x)
    use <- ifelse(na.rm, "pairwise.complete.obs", "everything")
    (n - 1) / n * cov(x = x, use = use)
}

split_by_season <- function(serie) {
    seasons <- c(cycle(serie))
    orig_order <- order(unlist(split(seq_along(serie), seasons)))
    splitted   <- split(serie, seasons)
    attr(splitted, "orig_order") <- orig_order
    return(splitted)
}

#' Normaliza Serie Sazonalmente
#' 
#' Normaliza uma serie temporal com sazonalidade utilizando medias e desvios padrao por estacao
#' 
#' Caso `means` e `sds` sejam fornecidos, devem ser vetores com comprimento igual a a frequencia
#' sazonal da serie, em que a primeira posicao corresponde à primeira estacao e assim por diante.
#' 
#' @param serie serie temporal com sazonalidade para normalizar
#' @param est um de `c("n", "n-1")` indicando qual denominador utilizar para desvio padrao
#' @param means,sds medias e desvios padrao sazonais ja calculadas (opcional)
#' 
#' @return lista de dois elementos: serie temporal padronizada sazonalmente e lista contendo medias
#'     e desvios padrao sazonais

scale_by_season <- function(serie, est = "n", means = NULL, sds = NULL) {

    attr0 <- attributes(serie)

    if (is.null(means)) means <- seasonal_mean(serie)
    if (is.null(sds))   sds   <- seasonal_sd(serie, est)

    serie <- split_by_season(serie)
    order0 <- attr(serie, "orig_order")
    serie <- mapply(function(x, m, s) (x - m) / s, serie, means, sds, SIMPLIFY = FALSE)
    serie <- unname(unlist(serie)[order0])

    attributes(serie) <- attr0

    out <- list(serie, list(means, sds))
    return(out)
}

seasonal_mean <- function(serie) {
    splitted <- split_by_season(serie)
    means <- sapply(splitted, function(x) mean(x, na.rm = TRUE))
    return(means)
}

seasonal_sd <- function(serie, est = "n") {
    fsd <- ifelse(est == "n-1", sd, sd2)
    splitted <- split_by_season(serie)
    sds <- sapply(splitted, function(x) fsd(x, na.rm = TRUE))
    return(sds)
}

#' Padding De Serie Temporal
#' 
#' Preenche uma serie temporal com NA no inicio e no fim para completar ciclos sazonais
#' 
#' @param serie serie temporal a ser preenchida
#' @param pad valor a ser utilizado para preencher as pontas (default: NA)
#' 
#' @return serie temporal preenchida com NA nas pontas

pad_series <- function(serie, pad = NA) {
    tsp   <- tsp(serie)
    freq  <- frequency(serie)
    season <- cycle(serie)

    deltaini <- head(season, 1) - 1
    deltafim <- freq - tail(season, 1)

    padini <- rep(pad, deltaini)
    padfim <- rep(pad, deltafim)

    out <- ts(c(padini, serie, padfim), start = tsp[1] - (deltaini) / tsp[3], freq = tsp[3])
    return(out)
}

#' Converte Serie Temporal em Matriz
#' 
#' Converte uma serie temporal em matriz, preenchendo com NA nas pontas para completar ciclos sazon
#' 
#' @param serie serie temporal a ser convertida
#' 
#' @return matriz com as observacoes organizadas por linhas de ciclos sazonais

ts2matrix <- function(serie) {
    dat <- pad_series(serie)
    out <- matrix(dat, ncol = frequency(serie), byrow = TRUE)
    return(out)
}

# AUXILIARES DE MODELO -----------------------------------------------------------------------------

#' Filtra Serie Temporal Com Modelo PAR
#' 
#' Filtra uma serie temporal utilizando um modelo PAR ajustado
#' 
#' @param object objeto da classe `par` contendo o modelo ajustado
#' 
#' @return serie temporal filtrada

filter_series <- function(object) {
    x <- scale_by_season(object$x)
    scales <- x[[2]]
    x <- x[[1]]

    s <- frequency(x)
    n <- length(x)
    fitted <- rep(NA_real_, n)

    for (t in seq_len(n)) {
        m <- cycle(x)[t]
        p <- length(object$phis[[m]])
        if (t > p) {
            fitted[t] <- sum(object$phis[[m]] * rev(x[(t - p):(t - 1)]), na.rm = TRUE) *
                scales[[2]][m] + scales[[1]][m]
        }
    }

    fitted <- ts(fitted, start = start(x), frequency = s)

    return(fitted)
}

#' Previsao Com Modelo PAR
#' 
#' Realiza previsao de uma serie temporal utilizando um modelo PAR ajustado
#' 
#' @param object objeto da classe `par` contendo o modelo ajustado
#' @param n.ahead numero de periodos a serem previstos
#' 
#' @return serie temporal contendo as previsoes

predict_series <- function(object, n.ahead) {
    x <- scale_by_season(object$x)
    scales <- x[[2]]
    x <- x[[1]]

    s <- frequency(x)
    n <- length(x)
    preds <- ts(rep(NA_real_, n.ahead), start = tail(time(x), 1) + 1 / s, frequency = s)

    for (h in seq_len(n.ahead)) {
        t <- n + h
        m <- cycle(preds)[h]
        p <- length(object$phis[[m]])
        if (h > p) {
            vals <- preds[(h - p):(h - 1)]
        } else {
            vals <- c(x, head(preds, h - 1))[(t - p):(t - 1)]
        }
        preds[h] <- sum(object$phis[[m]] * rev(vals), na.rm = TRUE)
    }

    for (h in seq_len(n.ahead)) {
        m <- cycle(preds)[h]
        preds[h] <- preds[h] * scales[[2]][m] + scales[[1]][m]
    }

    preds <- ts(preds, start = end(x) + c(0, 1), frequency = s)

    return(preds)
}