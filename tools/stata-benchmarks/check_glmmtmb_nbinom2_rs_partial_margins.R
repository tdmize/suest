# Run here after parsing the partial return; independently check native margins.
f <- readRDS('../../tests/testthat/fixtures/glmmtmb-nbinom2-rs-stata-partial-v2.rds')
for (y in names(f$components)) {
 s <- f$components[[y]]; b <- s$coefficients
 d <- f$data[!is.na(f$data[[y]]), ]; x <- d$x; X <- model.matrix(~x+z,d)
 g <- function(b) { A <- exp(2*b[5]); B <- exp(2*b[6]); C <- sqrt(A*B)*tanh(b[7]); mu <- exp(drop(X%*%b[1:3])+(A+2*C*x+B*x^2)/2); c(mean(mu),mean(mu*(b[2]+C+B*x))) }
 A <- exp(2*b[5]); B <- exp(2*b[6]); rho <- tanh(b[7]); C <- sqrt(A*B)*rho; E <- sqrt(A*B)*(1-rho^2)
 mu <- exp(drop(X%*%b[1:3])+(A+2*C*x+B*x^2)/2)
 G <- mu*cbind(X,0,A+C*x,B*x^2+C*x,E*x)
 H <- G*(b[2]+C+B*x)+mu*cbind(0,1,0,0,C,C+2*B*x,E)
 J <- rbind(colMeans(G),colMeans(H)); Jn <- numDeriv::jacobian(g,b)
 vs <- vapply(s$native_margins,`[[`,numeric(1),'variance'); va <- diag(J%*%s$vcov%*%t(J))
 print(data.frame(outcome=y,margin=names(vs),stata_variance=vs,analytic_variance=va,difference=vs-va,relative_SE=sqrt(vs/va)-1))
 cat('analytic vs Richardson gradient:',max(abs(J-Jn)),'\n')
 stopifnot(max(abs(J-Jn)) < 1e-7, max(abs(sqrt(vs/va)-1)) < 1e-5)
}
cat("NB2_RS_PARTIAL_NATIVE_MARGINS_COMPLETE=1\n")
