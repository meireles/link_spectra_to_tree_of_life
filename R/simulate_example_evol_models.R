#' Simulate a single realization of the OU process
#'
#' @param sig2 rate of the Wiener process
#' @param alpha strength of the mean reversion
#' @param theta mean longrun
#' @param init initial value
#' @param time over which the process unfolds
#' @param n number of time slices. integer. NULL defaults to 100
ou_sim = function(sig2, alpha, theta, init, time = 1, n = NULL){
    n     = ifelse(is.null(n),max(time * 100), n)
    dt    = time / n
    dw    = rnorm(n, 0, sqrt(dt))
    pull  = function(theta, x, dt){ alpha * (theta - x) *  dt}
    y     = c(init, rep(NA, n - 1))
    for (i in 2 : n){
        y[i] = y[i - 1] + pull(theta, y[i - 1], dt) + sig2 * dw[i - 1]
    }
    setNames(y, nm = seq(0, time, length.out = n))
}

#set.seed(08282018)

n = 3

brownian = sapply(1:n, function(x){
    ou_sim(sig2   = 0.5,
           alpha  = 0,
           theta  = 1,
           init   = 0,
           time   = 2,
           n      = 100)
})

ou_weak = sapply(1:n, function(x){
    ou_sim(sig2   = 0.5,
           alpha  = log(2) / 2,
           theta  = 2,
           init   = 0,
           time   = 2,
           n      = 100)
})

ou_strong = sapply(1:n, function(x){
    ou_sim(sig2   = 0.5,
           alpha  = log(2) * 3,
           theta  = 2,
           init   = 0,
           time   = 2,
           n      = 100)
})


png("figures/example_evol_model.png", width = 9, height = 3, units = "in", res = 800)
par(mfrow = c(1, 3))

yrange = range(brownian, ou_weak, ou_strong)

matplot(brownian, type = "l", lty = 1, lwd = 2,
        ylab = "trait", xlab = "time",
        main = "BM (random wiggle)", ylim = yrange,
        col = c("black", "blue", "orange") )

legend("topleft", legend = "rate = 0.2",
       fill = NULL, border = NA, bty = "n", adj = 0.25)

mtext("a", side = 1, adj = 0.01, outer = T, font = 2, line = -1.5)


matplot(ou_weak, type = "l", lty = 1, lwd = 2,
        ylab = "trait", xlab = "time",
        main = "OU (wiggle + weak pull)", ylim = yrange,
        col = c("black", "blue", "orange"))
arrows(116, 2, 106, 2, xpd = TRUE, lwd = 2.5, cex = 1, length = 0.1, col = "red")

legend("topleft", legend = c("rate = 0.2", "pull = 0.35"),
       fill = NULL, border = NA, bty = "n", adj = 0.25)
mtext("b", side = 1, adj = 0.35, outer = T, font = 2, line = -1.5)


matplot(ou_strong, type = "l", lty = 1, lwd = 2,
        ylab = "trait", xlab = "time",
        main = "OU (wiggle + strong pull)", ylim = yrange,
        col = c("black", "blue", "orange"))
arrows(116, 2, 106, 2, xpd = TRUE, lwd = 2.5, cex = 1, length = 0.1, col = "red")

legend("topleft", legend = c("rate = 0.2", "pull = 2.1"),
       fill = NULL, border = NA, bty = "n", adj = 0.25)
mtext("c", side = 1, adj = 0.70, outer = T, font = 2, line = -1.5)

dev.off()
