
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

    out <- list(serie, c(means, sds))
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
