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
#'     * `sigma2`: variancia dos residuos do modelo
#'     * `x`: serie temporal utilizada na estimacao
#'     * `call`: chamada da funcao
#' 
#' @export

par <- function(serie, ps = "auto", max_ps = frequency(serie) - 1, metodo = "YuleWalker") {

    s <- frequency(serie)
    if (s == 1) stop("'serie' nao possui sazonalidade")

    if (length(ps) < s) ps <- rep(ps, length.out = s)
    if (length(max_ps) < s) max_ps <- rep(max_ps, length.out = s)

    phis <- lapply(seq_len(s), function(m) fit_par(serie, m, ps[m], metodo, max_p = max_ps[m]))

    cc <- match.call()
    new <- new_par(x = serie, phis = phis, sigma2 = NA_real_, sigma2_norm = NA_real_,
        residuals = NULL, call = cc)

    res <- residuals(new)
    sigma2 <- var(res, na.rm = TRUE)
    res_norm <- scale_by_season(res)[[1]]
    sigma2_norm <- var(res_norm, na.rm = TRUE)

    new$residuals <- res
    new$sigma2 <- sigma2
    new$sigma2_norm <- sigma2_norm

    return(new)
}

fit_par <- function(serie, m, p, metodo, ...) {
    fitfun <- paste0("fit_par_", tolower(metodo))
    cc <- match.call()
    cc[[1]] <- as.name(fitfun)
    fit <- eval(cc, parent.frame(), parent.frame())
    return(fit)
}

fit_par_yulewalker <- function(serie, m, p, max_p, ...) {
    serie <- scale_by_season(serie, est = "n")[[1]]
    if (p == "auto") p <- idordem_yulewalker(serie, m, max_p)
    pacf_res <- perpacf(serie, m, p)
    phis <- as.numeric(solve(pacf_res$RHO, pacf_res$rho))
    return(phis)
}

#' Construtor Interno De `par`
#' 
#' Funcao interna para construcao de objetos da classe `par`
#' 
#' @param x serie temporal utilizada na estimacao
#' @param phis lista contendo os coeficientes de cada modelo autorregressivo periodico
#' @param sigma2 variancia dos residuos do modelo
#' @param residuals residuos do modelo
#' @param call chamada da funcao
#' 
#' @return objeto da classe `par`, uma lista contendo os seguintes elementos:
#'     * `phis`: lista contendo os coeficientes de cada modelo autorregressivo periodico
#'     * `sigma2`: variancia dos residuos do modelo
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

simulate.par <- function(object, n.ahead, nsims, seed = 1234, ...) {
    set.seed(seed)
    out <- lapply(seq_len(nsims), function(i) run_model_recursion(object, n.ahead, "simulation"))
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
    serie <- scale_by_season(serie, est = "n")[[1]]
    phis <- perpacf(serie, m, max_p)
    conf <- qnorm((1 + .95) / 2) / sqrt(phis$n_used)
    over <- abs(phis$phi) >= conf
    ordem <- ifelse(any(over), max(which(over)), 1)
    return(ordem)
}
