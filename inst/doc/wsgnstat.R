## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6, fig.height = 4,
                      dpi = 72)
has_circular <- requireNamespace("circular", quietly = TRUE)
if (has_circular) {
  wind_env <- new.env()
  utils::data("wind", package = "circular", envir = wind_env)
}

## ----setup--------------------------------------------------------------------
library(wsgnstat)

## ----family, fig.cap = "Wrapped densities for three tail parameters, symmetric (solid) and skew (dashed)."----
th <- seq(0, 2 * pi, length.out = 721)
op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (b in c(0.8, 1.5, 3)) {
  plot(th, dwsgn(th, mu = pi, beta = b, lam = 0.5), type = "l", ylim = c(0, 1.2),
       xlab = "angle (radians)", ylab = "density", main = paste("beta =", b))
  lines(th, dwsgn(th, mu = pi, beta = b, lam = 0.25), lty = 2)
}
par(op)

## ----dpqr---------------------------------------------------------------------
mu <- 1
beta <- 1.5
lam <- 0.3
dwsgn(c(0.5, 1, 2), mu, beta, lam)
pwsgn(c(0.5, 1, 2), mu, beta, lam)
qwsgn(c(0.1, 0.5, 0.9), mu, beta, lam)

## ----dpqr-check---------------------------------------------------------------
c(pwsgn = pwsgn(2, mu, beta, lam),
  integrate = integrate(dwsgn, 0, 2, mu = mu, beta = beta, lam = lam)$value,
  total = integrate(dwsgn, 0, 2 * pi, mu = mu, beta = beta, lam = lam)$value)

## ----rwsgn, fig.cap = "A sample of 5000 draws with the wrapped density overlaid."----
set.seed(1)
x <- rwsgn(5000, mu, beta, lam)
hist(x, breaks = seq(0, 2 * pi, length.out = 49), freq = FALSE, main = "",
     xlab = "angle (radians)")
curve(dwsgn(x, mu, beta, lam), 0, 2 * pi, add = TRUE, lwd = 2)

## ----ks-----------------------------------------------------------------------
u <- sort(pwsgn(x, mu, beta, lam))
n <- length(u)
c(KS = max(pmax((1:n) / n - u, u - (0:(n - 1)) / n)), five_per_cent_critical = 1.358 / sqrt(n))

## ----wn-----------------------------------------------------------------------
tt <- c(0.5, 2, 4.5)
k <- -60:60
wrapped_normal <- vapply(tt, function(t) sum(dnorm(t + 2 * pi * k, 0, sqrt(1 / 2))), 0)
rbind(wsgnstat = dwsgn(tt, mu = 0, beta = 2, lam = 0.5), wrapped_normal = wrapped_normal)

## wrapped split normal: standard deviation sigma / (sqrt(2) lam) on the right of mu and
## sigma / (sqrt(2) (1 - lam)) on the left, with left mass lam
l <- 0.3
s <- sqrt(l * (1 - l))
split_normal <- function(x) {
  sr <- s / (sqrt(2) * l)
  sl <- s / (sqrt(2) * (1 - l))
  ifelse(x >= 0, 2 * (1 - l) * dnorm(x, 0, sr), 2 * l * dnorm(x, 0, sl))
}
wrapped_split <- vapply(tt, function(t) sum(split_normal(t + 2 * pi * k)), 0)
rbind(wsgnstat = dwsgn(tt, mu = 0, beta = 2, lam = l), wrapped_split_normal = wrapped_split)

## ----tm-----------------------------------------------------------------------
tm <- trig_moments(4, mu = 0, beta = 2, lam = 0.5)
rbind(engine = tm$A, closed_form = exp(-(1:4)^2 / 4))

## ----wsn----------------------------------------------------------------------
dwsn <- function(theta, xi, omega, alpha, k = -60:60) {
  vapply(theta, function(t) {
    z <- (t + 2 * pi * k - xi) / omega
    sum(2 / omega * dnorm(z) * pnorm(alpha * z))
  }, 0)
}
rbind(wrapped_skew_normal_alpha0 = dwsn(tt, 0, sqrt(1 / 2), 0),
      wsgnstat = dwsgn(tt, mu = 0, beta = 2, lam = 0.5))

## ----sim----------------------------------------------------------------------
truth <- c(mu = 4, beta = 0.8, lam = 0.6)
set.seed(17)
y <- rwsgn(2000, truth["mu"], truth["beta"], truth["lam"])
ft <- fit_wsgn_trig(y)
fm <- fit_wsgn_mle(y)
round(rbind(truth = truth, moments = ft$par, moments_se = ft$se, likelihood = fm$par), 3)

## ----boot---------------------------------------------------------------------
truth2 <- c(mu = 1, beta = 1.5, lam = 0.3)
set.seed(18)
y2 <- rwsgn(500, truth2["mu"], truth2["beta"], truth2["lam"])
fb <- fit_wsgn_mle(y2, B = 30)
round(rbind(truth = truth2, estimate = fb$par, bootstrap_se = fb$se), 3)

## ----wind, eval = has_circular------------------------------------------------
w <- as.numeric(wind_env$wind) %% (2 * pi)
fw_t <- fit_wsgn_trig(w)
fw <- fit_wsgn_mle(w)
round(rbind(moments = fw_t$par, moments_se = fw_t$se, likelihood = fw$par), 3)

## ----wind-compare, eval = has_circular----------------------------------------
wc <- circular::circular(w, type = "angles", units = "radians", template = "none",
                         modulo = "asis", zero = 0, rotation = "counter")
f_wn <- circular::mle.wrappednormal(wc)
f_wc <- circular::mle.wrappedcauchy(wc)
ll <- c(wsgn = fw$logLik,
        wrapped_normal = sum(log(circular::dwrappednormal(wc, f_wn$mu, f_wn$rho))),
        wrapped_cauchy = sum(log(circular::dwrappedcauchy(wc, f_wc$mu, f_wc$rho))))
npar <- c(3, 2, 2)
round(cbind(logLik = ll, parameters = npar, AIC = -2 * ll + 2 * npar), 1)

## ----wind-plot, eval = has_circular, fig.cap = "Wind directions with the fitted wrapped skew generalised normal density."----
hist(w, breaks = seq(0, 2 * pi, length.out = 37), freq = FALSE, main = "",
     xlab = "direction (radians)")
curve(dwsgn(x, fw$par["mu"], fw$par["beta"], fw$par["lam"]), 0, 2 * pi, add = TRUE, lwd = 2)

## ----wind-gof, eval = has_circular, fig.cap = "Probability-probability plot of the fitted model."----
u <- sort(pwsgn(w, fw$par["mu"], fw$par["beta"], fw$par["lam"]))
n <- length(u)
plot((1:n) / (n + 1), u, xlab = "uniform plotting position", ylab = "fitted distribution function",
     pch = 20, cex = 0.6, asp = 1)
abline(0, 1, lty = 2)
V <- max((1:n) / n - u) + max(u - (0:(n - 1)) / n)
Vstar <- V * (sqrt(n) + 0.155 + 0.24 / sqrt(n))
kuiper_p <- function(v) {
  j <- 1:100
  min(1, sum(2 * (4 * j^2 * v^2 - 1) * exp(-2 * j^2 * v^2)))
}
c(V = V, V_star = Vstar, approximate_p = kuiper_p(Vstar))

