# Run here after the independent audit; diagnose native versus Fisher curvature.
pkgload::load_all('../..', quiet=TRUE)
f <- readRDS('../../tests/testthat/fixtures/glmmtmb-nbinom2-rs-stata-disjoint-v3.rds')
a <- readRDS('../../tests/testthat/fixtures/glmmtmb-nbinom2-rs-disjoint-audit.rds'); d <- f$data; d$id <- factor(d$id)
r <- readRDS('glmmtmb_nbinom2_rs_r_reference.rds')
for (i in 1:2) {
 y <- c('y1_left','y2_right')[i]
 m <- glmmTMB::glmmTMB(reformulate(c('x','z','(1+x|id)'),y),data=d[!is.na(d[[y]]),],family=glmmTMB::nbinom2('log'))
 b <- m$fit$par; J <- diag(7); J[7,7] <- 1/sqrt(1+b[7]^2)
 B <- a$R$components[[i]]$bread; gradient <- a$R$components[[i]]$gradient
 T <- matrix(0,7,7); T[7,7] <- gradient[7]*tanh(asinh(b[7]))
 Bnative <- solve(solve(B)+T)
 V <- J%*%m$sdr$cov.fixed%*%J
 cat('\n',y,'coef reference delta',max(abs(suest:::.suest_model_components(m,"glmm_nbinom2_rs","glmmTMB::glmmTMB")$parameters-r$components[[y]]$coefficients)),'\n')
 cat('gradient:', gradient,'\n')
 cat('native transform effect:',max(abs(Bnative-B)),'\n')
 stopifnot(max(abs(V-Bnative)) < 1e-6)
 cat('original native vs audit:',max(abs(V-B)), 'native coordinate audit:',max(abs(V-Bnative)),'\n')
 for (h in c(.001,.0005,.00025,.0001)) {
  H <- optimHess(b, m$obj$fn, m$obj$gr, control=list(ndeps=rep(h,7)))
  W <- J%*%solve(H)%*%J
  cat('step=',h,'vs supplied=',max(abs(W-V)),'vs independent native=',max(abs(W-Bnative)),'vs independent Fisher=',max(abs(W-B)),'\n')
 }
}
cat("NB2_RS_DISJOINT_NATIVE_CURVATURE_COMPLETE=1\n")
