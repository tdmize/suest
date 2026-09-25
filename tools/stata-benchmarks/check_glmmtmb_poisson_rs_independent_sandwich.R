# Run from tools/stata-benchmarks after building the v3 fixture.
repo <- '../..'
f <- readRDS(file.path(repo,'tests/testthat/fixtures/glmmtmb-poisson-rs-stata-components.rds'))
r <- readRDS(file.path(repo,'tools/stata-benchmarks/glmmtmb_poisson_rs_r_reference.rds'))$systems$balanced
s <- readRDS(file.path(repo,'tests/testthat/fixtures/glmmtmb-poisson-rs-stata-balanced-v3.rds'))
d <- f$data; groups <- split(seq_len(nrow(d)),d$id); X <- cbind(1,d$x,d$z); Z <- cbind(1,d$x)
# Independent conditional-mode solver and bivariate Laplace group likelihood.
group_ll <- function(b,y) {
  sd <- exp(b[4:5]); rho <- tanh(b[6]); S <- outer(sd,sd)*matrix(c(1,rho,rho,1),2)
  P <- solve(S); eta <- drop(X%*%b[1:3]); logdet <- as.numeric(determinant(S,logarithm=TRUE)$modulus)
  vapply(groups,function(idx) {
    zi <- Z[idx,,drop=FALSE]; yi <- y[idx]; ei <- eta[idx]
    objective <- function(u) sum(exp(ei+drop(zi%*%u))-yi*(ei+drop(zi%*%u)))+sum(u*(P%*%u))/2
    u <- c(0,0)
    for(k in 1:100) {
      mu <- exp(ei+drop(zi%*%u)); g <- drop(crossprod(zi,mu-yi)+P%*%u)
      if(max(abs(g)) < 1e-10) break
      H <- crossprod(zi,zi*mu)+P; step <- drop(solve(H,g)); scale <- 1
      # Permit rounding-level equality near the mode; gradient verifies accuracy.
      while(objective(u-scale*step) > objective(u)+1e-12 && scale>1e-8) scale <- scale/2
      u <- u-scale*step
    }
    mu <- exp(ei+drop(zi%*%u)); g <- drop(crossprod(zi,mu-yi)+P%*%u)
    stopifnot(max(abs(g))<1e-8)
    H <- crossprod(zi,zi*mu)+P
    sum(dpois(yi,mu,log=TRUE))-sum(u*(P%*%u))/2-logdet/2-as.numeric(determinant(H,logarithm=TRUE)$modulus)/2
  },numeric(1))
}
results <- lapply(c('y1','y2'),function(y) {
  b <- f$components[[y]]$coefficients; V <- f$components[[y]]$vcov
  scores <- lapply(c(1e-5,5e-6),function(h) sapply(seq_along(b),function(j) {
    plus <- minus <- b; plus[j] <- plus[j]+h; minus[j] <- minus[j]-h
    (group_ll(plus,d[[y]])-group_ll(minus,d[[y]]))/(2*h)
  }))
  cat(y,'LL=',sum(group_ll(b,d[[y]])),'StataLL=',f$components[[y]]$logLik,'\n')
  cat(y,'max score change when halving finite-difference step=',max(abs(scores[[1]]-scores[[2]])),'\n')
  g <- colSums(scores[[2]]); A <- exp(2*b[4]); B <- exp(2*b[5]); rho <- tanh(b[6])
  var_slope_gradient <- g[5]/(2*B)-g[6]*rho/(2*B*(1-rho^2))
  cat(y,'sum gradient in native var_slope coordinate=',var_slope_gradient,'\n')
  list(influence=scores[[2]]%*%V, gradient=g)
})
V <- crossprod(cbind(results[[1]]$influence,results[[2]]$influence))*100/99
cat('Independent sandwich versus R maxabsV=',max(abs(V-r$vcov)),', maxrelativeSE=',max(abs(sqrt(diag(V)/diag(r$vcov))-1)),'\n')
cat('Independent sandwich versus GSEM maxabsV=',max(abs(V-s$vcov)),', maxrelativeSE=',max(abs(sqrt(diag(s$vcov)/diag(V))-1)),'\n')
stopifnot(max(abs(V-r$vcov))<1e-6,max(abs(sqrt(diag(V)/diag(r$vcov))-1))<1e-4)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) saveRDS(list(vcov=V,components=results),args[1])
cat('INDEPENDENT_POISSON_RS_SANDWICH_COMPLETE=1\n')
