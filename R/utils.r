
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
