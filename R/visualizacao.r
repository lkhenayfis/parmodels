#' Plota a PACF periodica
#' 
#' Plota a funcao de autocorrelacao parcial periodica calculada pela funcao `perpacf`.
#' 
#' @param x objeto da classe `periodic_pacf`
#' @param ... argumentos adicionais para a funcao `plot`
#' 
#' @return nenhum, apenas plota o grafico
#' 
#' @export

plot.periodic_pacf <- function(x, ...) {

    lags <- seq_along(x$phi)
    acfs <- x$phi
    conf <- qnorm((1 + .95) / 2) / sqrt(x$n_used) # intervalo 95%
    ylim <- range(-conf, conf, range(acfs))

    plot(lags, acfs, type = "h", ylim = ylim, ylab = "PACF Periodica", ...)
    title(paste0("Periodo m = ", x$m))
    abline(h = 0)
    abline(h = c(-conf, conf), lty = 2, col = "blue")
}

#' @export

plot_simulation <- function(simul, predicted = NULL, ref = NULL) {
    plot(simul[, 1], ylim = range(simul), col = 4)
    for (i in seq_len(ncol(simul))[-1]) lines(simul[, i], col = 4)

    if (!is.null(ref)) lines(ref, lty = 1, col = 2)
    if (!is.null(predicted)) lines(predicted, lty = 2, col = 1)
}