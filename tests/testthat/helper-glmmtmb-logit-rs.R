# Independent adaptive integration and analytic gradients for the Stata fixture.
glmmtmb_logit_rs_native_margins <- function(b, V, data) {
  X <- cbind(1, data$x, data$z); x <- data$x; eta <- drop(X %*% b[1:3])
  A <- exp(2*b[4]); B <- exp(2*b[5]); rho <- tanh(b[6]); C <- sqrt(A*B)*rho
  E <- sqrt(A*B)*(1-rho^2); sigma <- sqrt(A+2*C*x+B*x^2)
  m <- vapply(0:4, function(k) vapply(seq_along(x), function(i) integrate(function(u) {
    p <- plogis(eta[i]+sigma[i]*u); p1 <- p*(1-p)
    switch(k+1L, p, p1, p1*(1-2*p), p1*(1-6*p+6*p^2),
      p1*(1-14*p+36*p^2-24*p^3))*dnorm(u)
  }, -Inf, Inf, rel.tol = 1e-12, abs.tol = 1e-13)$value, numeric(1)), numeric(length(x)))
  XX <- cbind(X, 0, 0, 0)
  Vg <- cbind(0, 0, 0, 2*(A+C*x), 2*(B*x^2+C*x), 2*E*x)
  G <- m[, 2]*XX+.5*m[, 3]*Vg
  H <- b[2]*(m[, 3]*XX+.5*m[, 4]*Vg)+(C+B*x)*(m[, 4]*XX+.5*m[, 5]*Vg)+
    m[, 2]*matrix(rep(c(0, 1, 0, 0, 0, 0), each = length(x)), ncol = 6)+
    m[, 3]*cbind(0, 0, 0, rep(C, length(x)), C+2*B*x, rep(E, length(x)))
  gradient <- rbind(colMeans(G), colMeans(H))
  list(estimates = c(mean = mean(m[, 1]), slope_x = mean(b[2]*m[, 2]+(C+B*x)*m[, 3])),
    gradient = gradient, vcov = gradient %*% V %*% t(gradient))
}
