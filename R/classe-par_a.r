#' Estimacao De `par_a`s
#'
#' Estimacao de modelos autorregressivos periodicos com componente anual (PAR(p)-A)
#'
#' Um `par_a` e ajustado resolvendo, por estacao, o sistema de Yule-Walker *aumentado* pela componente
#' anual `A` (media dos `frequency(serie)` valores imediatamente anteriores). O coeficiente anual `psi`
#' resultante e entao "desdobrado" nos coeficientes de cada um dos 12 lags do ultimo ano, de modo que o
#' objeto retornado e, mecanicamente, um `par` de ordem `frequency(serie)`: todos os metodos de `par`
#' (`fitted`/`residuals`/`predict`/`simulate`) funcionam por heranca, sem recursao especializada.
#'
#' O argumento `ps` segue a mesma convencao de [par()]: escalar (replicado) ou vetor de comprimento
#' igual a frequencia sazonal; `"auto"` identifica a ordem automaticamente.
#'
#' @param serie serie temporal com sazonalidade
#' @param ps ordens autorregressivas de cada modelo periodico. Veja [par()]
#' @param max_ps valor maximo a considerar na identificacao automatica de ordem
#' @param metodo metodo de estimacao. Atualmente, apenas "Yule-Walker" e suportado
#'
#' @return objeto de classe `c("par_a", "par", "list")`, um `par` cujos `phis` ja incorporam o termo
#'     anual desdobrado, acrescido de:
#'     * `psi`: vetor com o coeficiente anual estimado de cada estacao
#'     * `ps`: vetor com a ordem autorregressiva identificada de cada estacao
#'     * `annual`: serie temporal das medias anuais utilizadas no ajuste (para validacao)
#'
#' @examples
#' serie <- ts(dummyseries[, 2], start = c(1931, 1), frequency = 12)
#' mod <- par_a(serie)
#'
#' @export

par_a <- function(serie, ps = "auto", max_ps = frequency(serie) - 1, metodo = "YuleWalker") {

    s <- frequency(serie)
    if (s == 1) stop("'serie' nao possui sazonalidade")

    if (length(ps) < s) ps <- rep(ps, length.out = s)
    if (length(max_ps) < s) max_ps <- rep(max_ps, length.out = s)

    serie_A <- calcula_medias_anuais(serie)

    fits <- lapply(seq_len(s),
        function(m) fit_par_a(serie, serie_A, m, ps[m], metodo, max_p = max_ps[m]))

    sds  <- seasonal_sd(serie, est = "n")
    sdsA <- seasonal_sd(serie_A, est = "n")

    phis <- lapply(seq_len(s),
        function(m) unfold_par_a(fits[[m]]$phis, fits[[m]]$psi, sds, sdsA, m, s))

    sigma2_norm <- vapply(fits, function(f) f$sigma2_norm, numeric(1))
    sigma2      <- sigma2_norm * sds^2

    cc  <- match.call()
    new <- new_par(x = serie, phis = phis, sigma2 = sigma2, sigma2_norm = sigma2_norm,
        residuals = NULL, call = cc)
    class(new) <- c("par_a", "par", "list")
    new$psi    <- vapply(fits, function(f) f$psi, numeric(1))
    new$ps     <- vapply(fits, function(f) length(f$phis), integer(1))
    new$annual <- serie_A

    new$residuals <- residuals(new)

    return(new)
}

#' Ajuste Por Estacao De `par_a`
#'
#' Funcoes internas para ajuste do modelo periodico anual de uma unica estacao
#'
#' @param serie serie temporal sazonal (bruta)
#' @param serie_A componente anual ja calculada, repassada para evitar recomputo a cada estacao
#' @param m estacao do ano a ajustar
#' @param p ordem autorregressiva, ou `"auto"` para identificacao automatica
#' @param metodo metodo de estimacao
#' @param max_p ordem maxima na identificacao automatica
#' @param ... demais argumentos repassados ao metodo
#'
#' @return lista com `phis` (coeficientes autorregressivos), `psi` (coeficiente anual) e
#'     `sigma2_norm` (variancia da inovacao na escala normalizada sazonalmente)
#'
#' @rdname fit_par_a

fit_par_a <- function(serie, serie_A, m, p, metodo, ...) {
    fitfun  <- paste0("fit_par_a_", tolower(metodo))
    cc      <- match.call()
    cc[[1]] <- as.name(fitfun)
    fit     <- eval(cc, parent.frame(), parent.frame())
    return(fit)
}

#' @rdname fit_par_a

fit_par_a_yulewalker <- function(serie, serie_A, m, p, max_p, ...) {
    if (p == "auto") p <- idordem_cacf(serie, serie_A, m, max_p)
    cacf <- percacf(serie, m, serie_A, lag_max = p, est = "n")
    sol  <- as.numeric(solve(cacf$RHO, cacf$rho))
    sigma2_norm <- 1 - sum(sol * cacf$rho)
    list(phis = sol[seq_len(p)], psi = sol[[p + 1]], sigma2_norm = sigma2_norm)
}

#' Identificacao Automatica De Ordem De Modelos `par_a`
#'
#' Funcao interna que identifica a ordem autorregressiva via FACP condicional (`percacf`), seguindo o
#' mesmo procedimento de [idordem_yulewalker()]: mantem os lags cujo coeficiente excede a banda de 95%
#' (`+/- z_0.975 / sqrt(n)`), ate `max_p`.
#'
#' @param serie serie temporal sazonal (bruta)
#' @param serie_A componente anual ja calculada
#' @param m estacao do ano
#' @param max_p numero maximo de lags a considerar
#'
#' @return ordem identificada
#'
#' @rdname idordem_cacf

idordem_cacf <- function(serie, serie_A, m, max_p) {
    cacf  <- percacf(serie, m, serie_A, lag_max = max_p, est = "n")
    conf  <- qnorm((1 + .95) / 2) / sqrt(cacf$n_used)
    over  <- abs(cacf$phi) >= conf
    ordem <- ifelse(any(over), max(which(over)), 1)
    return(ordem)
}

#' Desdobramento Do Termo Anual Em Coeficientes De Lag
#'
#' Converte `(phis, psi)` na estacao `m` no vetor de `s` coeficientes equivalentes de um `par`, somando
#' a contribuicao do termo anual padronizado a cada lag do ultimo ano:
#' `c_l = phi_l + psi / (s * sdsA[m-1]) * sds[m-l]`, com `phi_l = 0` para `l > p`.
#'
#' @param phis coeficientes autorregressivos da estacao
#' @param psi coeficiente do termo anual da estacao
#' @param sds desvios padrao sazonais de `Z` (denominador "n")
#' @param sdsA desvios padrao sazonais da componente anual (denominador "n")
#' @param m estacao do ano
#' @param s frequencia sazonal
#'
#' @return vetor de comprimento `s` com os coeficientes desdobrados
#'
#' @rdname unfold_par_a

unfold_par_a <- function(phis, psi, sds, sdsA, m, s) {
    p   <- length(phis)
    c_l <- numeric(s)
    for (l in seq_len(s)) {
        phil   <- if (l <= p) phis[l] else 0
        c_l[l] <- phil + psi / (s * sdsA[wrap_season(m - 1, s)]) * sds[wrap_season(m - l, s)]
    }
    return(c_l)
}
