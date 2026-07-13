#' Estimacao De `par`s
#' 
#' Estimacao de modelos autorregressivos periodicos
#' 
#' O argumento `ps` pode ser ou um escalar ou um vetor. No caso de ser escalar, sera repetido para N
#' vezes, onde N é a frequencia sazonal da serie. Caso seja um vetor, deve ter comprimento igual a
#' frequencia sazonal da serie. Em ambos os casos, o valor em cada posicao indica a ordem do modelo
#' autorregressivo periodico para a estacao correspondente. Sendo escalar ou vetor, o valor "auto"
#' pode ser utilizado para indicar que a ordem deve ser identificada automaticamente utilizando o
#' criterio correspondente a `metodo`: se Yule-Walker, analise das pacfs periodicas.
#' 
#' @param serie serie temporal com sazonalidade
#' @param ps ordens de cada modelo periodico. Veja Detalhes
#' @param max_ps se nao forem fornecidas ordens, valor maximo a ser considerado para identificacao
#'     automatica
#' @param metodo metodo de estimacao. Atualmente, apenas "Yule-Walker" é suportado
#' 
#' @return objeto da classe `par`, uma lista contendo os seguintes elementos:
#'     * `phis`: lista contendo os coeficientes de cada modelo autorregressivo periodico
#'     * `sigma2`: vetor de comprimento igual a frequencia sazonal com a variancia (escala bruta) da
#'       inovacao de cada estacao, calculada em forma fechada via Yule-Walker
#'     * `x`: serie temporal utilizada na estimacao
#'     * `call`: chamada da funcao
#'
#' @export

par <- function(serie, ps = "auto", max_ps = frequency(serie) - 1, metodo = "YuleWalker") {

    s <- frequency(serie)
    if (s == 1) stop("'serie' nao possui sazonalidade")

    if (length(ps) < s) ps <- rep(ps, length.out = s)
    if (length(max_ps) < s) max_ps <- rep(max_ps, length.out = s)

    fits <- lapply(seq_len(s), function(m) fit_par(serie, m, ps[m], metodo, max_p = max_ps[m]))
    phis <- lapply(fits, function(f) f$phis)
    sigma2_norm <- vapply(fits, function(f) f$sigma2_norm, numeric(1))

    sds <- seasonal_sd(serie, est = "n")
    sigma2 <- sigma2_norm * sds^2

    cc <- match.call()
    new <- new_par(x = serie, phis = phis, sigma2 = sigma2, sigma2_norm = sigma2_norm,
        residuals = NULL, call = cc)

    new$residuals <- residuals(new)

    return(new)
}

fit_par <- function(serie, m, p, metodo, ...) {
    fitfun <- paste0("fit_par_", tolower(metodo))
    cc <- match.call()
    cc[[1]] <- as.name(fitfun)
    fit <- eval(cc, parent.frame(), parent.frame())
    return(fit)
}

#' Ajuste Por Estacao De `par`
#'
#' Funcao interna de ajuste de Yule-Walker para uma unica estacao
#'
#' Alem dos coeficientes autorregressivos, calcula em forma fechada a variancia da inovacao na
#' escala normalizada sazonalmente: `sigma2_norm = 1 - sum(phis * rho)`, onde `rho` e o vetor de
#' correlacoes do lado direito do sistema de Yule-Walker (ver [perpacf()]).
#'
#' @param serie serie temporal sazonal (bruta)
#' @param m estacao do ano a ajustar
#' @param p ordem autorregressiva, ou `"auto"` para identificacao automatica
#' @param max_p ordem maxima na identificacao automatica
#' @param ... demais argumentos repassados ao metodo
#'
#' @return lista com `phis` (coeficientes autorregressivos) e `sigma2_norm` (variancia da inovacao
#'     na escala normalizada sazonalmente)

fit_par_yulewalker <- function(serie, m, p, max_p, ...) {
    if (p == "auto") p <- idordem_yulewalker(serie, m, max_p)
    pacf_res <- perpacf(serie, m, p)
    phis <- as.numeric(solve(pacf_res$RHO, pacf_res$rho))
    sigma2_norm <- 1 - sum(phis * pacf_res$rho)
    list(phis = phis, sigma2_norm = sigma2_norm)
}

#' Construtor Interno De `par`
#' 
#' Funcao interna para construcao de objetos da classe `par`
#' 
#' @param x serie temporal utilizada na estimacao
#' @param phis lista contendo os coeficientes de cada modelo autorregressivo periodico
#' @param sigma2 vetor de comprimento igual a frequencia sazonal com a variancia (escala bruta) da
#'     inovacao de cada estacao
#' @param sigma2_norm vetor de comprimento igual a frequencia sazonal com a variancia da inovacao
#'     de cada estacao na escala normalizada sazonalmente
#' @param residuals residuos do modelo
#' @param call chamada da funcao
#'
#' @return objeto da classe `par`, uma lista contendo os seguintes elementos:
#'     * `phis`: lista contendo os coeficientes de cada modelo autorregressivo periodico
#'     * `sigma2`: vetor de comprimento igual a frequencia sazonal com a variancia (escala bruta)
#'       da inovacao de cada estacao
#'     * `sigma2_norm`: idem, na escala normalizada sazonalmente
#'     * `x`: serie temporal utilizada na estimacao
#'     * `residuals`: residuos do modelo
#'     * `call`: chamada da funcao

new_par <- function(x, phis, sigma2, sigma2_norm, residuals, call) {
    new <- list(
        phis = phis,
        sigma2 = sigma2,
        sigma2_norm = sigma2_norm,
        x = x,
        residuals = residuals,
        call = call
    )
    class(new) <- c("par", "list")
    return(new)
}

# METODOS ------------------------------------------------------------------------------------------

#' Metodos De Modelos `par`
#' 
#' Metodos S3 para objetos da classe `par`
#' 
#' @param object objeto da classe `par` contendo o modelo ajustado
#' @param n.ahead numero de passos a frente a serem previstos
#' @param nsim inteiro, numero de simulacoes a realizar
#' @param seed semente para simulacao
#' @param ... argumentos adicionais
#' 
#' @return varios, dependendo do metodo
#' 
#' @rdname par_methods
#' 
#' @export

fitted.par <- function(object, ...) {
    if (!is.null(residuals(object))) {
        res <- residuals(object)
        fitted <- object$x - res
    } else {
        fitted <- filter_series(object)
    }

    return(fitted)
}

#' @rdname par_methods
#' 
#' @export

fitted.values.par <- function(object, ...) fitted.par(object, ...)

#' @rdname par_methods
#' 
#' @export

residuals.par <- function(object, ...) {
    if (is.null(object$residuals)) {
        fitted <- filter_series(object)
        res <- object$x - fitted
    } else {
        res <- object$residuals
    }

    return(res)
}

#' @rdname par_methods
#' 
#' @export

predict.par <- function(object, n.ahead, ...) {
    run_model_recursion(object, n.ahead)
}

#' @rdname par_methods
#' 
#' @export 

simulate.par <- function(object, nsim, seed = 1234, n.ahead = 1, ...) {
    set.seed(seed)
    out <- lapply(seq_len(nsim), function(i) run_model_recursion(object, n.ahead, "simulation"))
    list2mts(out)
}

list2mts <- function(list) {
    names(list) <- paste0("simulation_", seq_along(list))
    tsp <- tsp(list[[1]])
    mts <- do.call(cbind, list)
    mts <- ts(mts)
    tsp(mts) <- tsp
    mts
}

# AUXILIARES ---------------------------------------------------------------------------------------

#' Identificacao Automatica De Ordem De Modelos Autorregressivos Periodicos
#' 
#' Funcao interna para identificar automaticamente a ordem de modelos autorregressivos periodicos
#' 
#' @param serie serie temporal sazonal
#' @param m estacao do ano para a qual identificar a ordem
#' @param max_p numero maximo de lags a considerar
#' 
#' @return ordem identificada

idordem_yulewalker <- function(serie, m, max_p) {
    phis <- perpacf(serie, m, max_p)
    conf <- qnorm((1 + .95) / 2) / sqrt(phis$n_used)
    over <- abs(phis$phi) >= conf
    ordem <- ifelse(any(over), max(which(over)), 1)
    return(ordem)
}
