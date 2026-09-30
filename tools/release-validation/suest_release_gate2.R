# Independent extended-family release gate. Run from the package source root.
# source("tools/release-validation/suest_release_gate2.R")
# References use distribution primitives and numDeriv, never suest/sandwich
# scores, bread, covariance, coefficient setters, or prediction helpers.
run_suest_release_gate2 <- function(root = getwd(), output = root,
                                    selected = NULL, patterns_selected = NULL,
                                    truncreg_method = "NR") {
  root <- normalizePath(root, mustWork = TRUE)
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  output <- normalizePath(output, mustWork = TRUE)
  con <- file(file.path(output, "suest_release_gate2.log"), "wt")
  old_message <- sink.number(type = "message")
  sink(con, split = TRUE); sink(con, type = "message")
  on.exit({
    sink(if (old_message == 2L) NULL else getConnection(old_message), type = "message")
    sink(); close(con)
  }, add = TRUE)
  old <- options(warn = 1, digits = 15, width = 160)
  on.exit(options(old), add = TRUE)
  cat("SUEST EXTENDED-FAMILY RELEASE GATE 2\nStarted:", format(Sys.time(), tz = "UTC"), "UTC\n")
  required <- c("pkgload", "numDeriv", "marginaleffects", "betareg", "survival",
                "censReg", "truncreg", "Rchoice", "pscl")
  stopifnot(all(vapply(required, requireNamespace, logical(1), quietly = TRUE)))
  pkgload::load_all(root, quiet = TRUE)
  cat("Source:", root, "\n")
  print(tools::md5sum(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE)))
  print(sapply(required, function(p) as.character(utils::packageVersion(p))))
  bounds <- c(loglik_gap = 2e-6, stationary = 2e-3, derivative_stability = 5e-5, score_stability = 5e-5,
              covariance = 3e-4, estimate = 2e-7, contrast_se = 3e-4)
  cat("Preset bounds:\n"); print(bounds)
  cat("Truncated-regression optimizer:", truncreg_method, "\n")
  results <- list()
  run_case <- function(name, fun) {
    cat("\nCASE:", name, "\n")
    warnings <- character()
    ans <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
    }), error = function(e) list(status = "ERROR", error = conditionMessage(e)))
    if (length(warnings) && identical(ans$status, "PASS")) ans$status <- "REVIEW_WARNING"
    ans$warnings <- unique(warnings)
    results[[name]] <<- ans
    print(ans)
    saveRDS(results, file.path(output, "suest_release_gate2_results.rds"))
  }
  set.seed(20260929)
  n <- 1600L
  d <- data.frame(id = seq_len(n), cluster = rep(seq_len(n/4L), each = 4L),
                  x = runif(n, -1, 1), z = runif(n, -1, 1))
  shared <- rep(rnorm(n/4L), each = 4L)
  e1 <- (.35*shared + rnorm(n))/sqrt(1+.35^2)
  e2 <- (.25*shared + .35*e1 + rnorm(n))/sqrt(1+.25^2+.35^2+2*.25*.35*.35/sqrt(1+.35^2))
  for (j in 1:2) {
    e <- if (j == 1) e1 else e2
    eta <- if (j == 1) .2+.55*d$x-.25*d$z else -.1+.35*d$x+.2*d$z
    u <- pnorm(e)
    mu <- plogis(eta); phi <- exp(2.7+.3*d$z)
    d[[paste0("beta",j)]] <- qbeta(u, mu*phi, (1-mu)*phi)
    d[[paste0("latent",j)]] <- eta+.8*e
    d[[paste0("binary",j)]] <- as.integer(eta+exp(.3*d$z)*e>0)
    zi <- plogis(-1+.25*d$x+.4*d$z)
    structural <- runif(n)<zi
    d[[paste0("zip",j)]] <- ifelse(structural, 0, qpois(u,exp(eta+.7)))
    d[[paste0("zinb",j)]] <- ifelse(structural, 0, qnbinom(u,size=2.5,mu=exp(eta+.7)))
    d[[paste0("lognormal",j)]] <- exp(.3+.25*eta+.55*e)
    d[[paste0("weibull",j)]] <- exp(.3+.25*eta)*(-log1p(-u))^.7
    d[[paste0("loglogistic",j)]] <- exp(.3+.25*eta+.5*qlogis(u))
  }
  linkfun <- function(eta, link) switch(link,
    logit = plogis(eta), probit = pnorm(eta), cloglog = -expm1(-exp(eta)),
    loglog = exp(-exp(-eta)), stop("Unknown reference link"))
  linkderiv <- function(eta, link) switch(link,
    logit = plogis(eta)*plogis(-eta), probit = dnorm(eta),
    cloglog = exp(eta-exp(eta)), loglog = exp(-eta-exp(-eta)))
  specs <- list()
  for (l in c("logit", "probit", "cloglog", "loglog"))
    specs[[paste0("beta_",l)]] <- list(kind="beta", link=l, phi="log", variable=TRUE, info="expected")
  specs$beta_constant_identity <- list(kind="beta",link="logit",phi="identity",variable=FALSE,info="expected")
  specs$beta_observed <- list(kind="beta",link="logit",phi="log",variable=TRUE,info="observed")
  for (s in c("left","right","interval","fixed"))
    specs[[paste0("survreg_gaussian_",s)]] <- list(kind="survreg",dist="gaussian",censor=s)
  for (s in c("weibull","lognormal","loglogistic")) {
    specs[[paste0("survreg_",s)]] <- list(kind="survreg",dist=s,censor="right")
    specs[[paste0("survreg_",s,"_interval")]] <- list(kind="survreg",dist=s,censor="interval")
  }
  specs$censreg_left <- list(kind="censreg",left=0,right=Inf)
  specs$censreg_right <- list(kind="censreg",left=-Inf,right=1)
  specs$censreg_twosided <- list(kind="censreg",left=-.5,right=1.3)
  specs$truncreg_left <- list(kind="truncreg",point=0,direction="left")
  specs$truncreg_right <- list(kind="truncreg",point=1,direction="right")
  for (l in c("probit","logit")) specs[[paste0("heteroskedastic_",l)]] <- list(kind="hetprob",link=l)
  for (f in c("zip","zinb")) for (l in c("logit","probit"))
    specs[[paste0(f,"_",l)]] <- list(kind=f,link=l)
  stopifnot(length(specs)==27L)
  patterns <- list(identical=list(rep(TRUE,n),rep(TRUE,n)),
    partial=list(d$id%%4!=0,d$id%%4!=1),
    disjoint_shared_clusters=list(d$id%%2==0,d$id%%2==1),
    disjoint_clusters=list(d$cluster<=n/8,d$cluster>n/8))
  if (!is.null(selected)) specs <- specs[selected]
  if (!is.null(patterns_selected)) patterns <- patterns[patterns_selected]
  make_component <- function(spec, keep, j) {
    dat <- d[keep, ]
    X <- model.matrix(~x+z, dat)
    Z <- model.matrix(~z,dat)
    kind <- spec$kind
    fitting_warnings <- character()
    fit_model <- function(expr) withCallingHandlers(expr, warning = function(w) {
      fitting_warnings <<- c(fitting_warnings, paste(paste(deparse(conditionCall(w)), collapse=" "), conditionMessage(w), sep=": "))
    })
    if (kind == "beta") {
      dat$y <- dat[[paste0("beta",j)]]
      if (!spec$variable) Z <- matrix(1,nrow(dat),1)
      form <- if (spec$variable) y~x+z|z else y~x+z|1
      m <- fit_model(betareg::betareg(form, data=dat, link=spec$link,link.phi=spec$phi,
        type="ML",model=TRUE,x=TRUE,y=TRUE,
        control=betareg::betareg.control(hessian=spec$info=="observed")))
      par <- coef(m)
      ll <- function(b) {
        mu <- linkfun(drop(X%*%b[1:3]),spec$link)
        phi_eta <- drop(Z%*%b[-(1:3)])
        phi <- if(spec$phi=="log") exp(phi_eta) else phi_eta
        dbeta(dat$y,mu*phi,(1-mu)*phi,log=TRUE)
      }
      info <- if(spec$info=="expected") function(b) {
        eta <- drop(X%*%b[1:3]); mu <- linkfun(eta,spec$link)
        pe <- drop(Z%*%b[-(1:3)])
        phi <- if(spec$phi=="log") exp(pe) else pe
        dp <- if(spec$phi=="log") phi else rep(1,length(phi))
        dm <- linkderiv(eta,spec$link)
        a <- trigamma(mu*phi); bta <- trigamma((1-mu)*phi)
        aa <- crossprod(X,X*(phi^2*(a+bta)*dm^2))
        ab <- crossprod(X,Z*(phi*(mu*a-(1-mu)*bta)*dm*dp))
        bb <- crossprod(Z,Z*((mu^2*a+(1-mu)^2*bta-trigamma(phi))*dp^2))
        rbind(cbind(aa,ab),cbind(t(ab),bb))
      } else NULL
      pred <- function(b,nd) linkfun(drop(model.matrix(~x+z,nd)%*%b[1:3]),spec$link)
    } else if (kind %in% c("zip","zinb")) {
      dat$y <- dat[[paste0(kind,j)]]
      m <- fit_model(pscl::zeroinfl(y~x+z|z,data=dat,dist=if(kind=="zip")"poisson" else "negbin",
        link=spec$link,model=TRUE,x=TRUE,y=TRUE,control=pscl::zeroinfl.control(maxit=1000,reltol=1e-12)))
      par <- if(kind=="zip")coef(m) else c(coef(m),ln_theta=log(m$theta))
      ll <- function(b) {
        mu <- exp(drop(X%*%b[1:3])); pi <- linkfun(drop(Z%*%b[4:5]),spec$link)
        logcount <- if(kind=="zip")dpois(dat$y,mu,log=TRUE) else dnbinom(dat$y,mu=mu,size=exp(b[6]),log=TRUE)
        ifelse(dat$y==0,log(pi+(1-pi)*exp(logcount)),log1p(-pi)+logcount)
      }
      info <- NULL
      pred <- function(b,nd) exp(drop(model.matrix(~x+z,nd)%*%b[1:3]))*
        (1-linkfun(drop(model.matrix(~z,nd)%*%b[4:5]),spec$link))
    } else if (kind == "hetprob") {
      dat$y <- dat[[paste0("binary",j)]]
      m <- fit_model(eval(substitute(Rchoice::hetprob(y~x+z|z,data=dat,link=L),list(L=spec$link))))
      par <- coef(m)
      ll <- function(b) {
        eta <- drop(X%*%b[1:3])/exp(dat$z*b[4])
        if(spec$link=="probit")pnorm((2*dat$y-1)*eta,log.p=TRUE) else plogis((2*dat$y-1)*eta,log.p=TRUE)
      }
      info <- NULL
      pred <- function(b,nd) linkfun(drop(model.matrix(~x+z,nd)%*%b[1:3])/exp(nd$z*b[4]),spec$link)
    } else if (kind == "truncreg") {
      dat$y <- dat[[paste0("latent",j)]]
      dat <- dat[if(spec$direction=="left")dat$y>spec$point else dat$y<spec$point, ]
      X <- model.matrix(~x+z,dat)
      m <- fit_model(truncreg::truncreg(y~x+z,data=dat,point=spec$point,direction=spec$direction,
                             model=TRUE,x=TRUE,y=TRUE,iterlim=200,method=truncreg_method))
      par <- coef(m)
      ll <- function(b) dnorm(dat$y,drop(X%*%b[1:3]),b[4],log=TRUE)-
        pnorm((spec$point-drop(X%*%b[1:3]))/b[4],lower.tail=spec$direction=="right",log.p=TRUE)
      info <- NULL
      pred <- function(b,nd) drop(model.matrix(~x+z,nd)%*%b[1:3])
    } else if (kind == "censreg") {
      dat$latent <- dat[[paste0("latent",j)]]
      dat$y <- pmin(spec$right,pmax(spec$left,dat$latent))
      m <- fit_model(censReg::censReg(y~x+z,data=dat,left=spec$left,right=spec$right))
      par <- coef(m)
      ll <- function(b) {
        mu <- drop(X%*%b[1:3]); sd <- exp(b[4])
        ifelse(dat$y<=spec$left,pnorm((spec$left-mu)/sd,log.p=TRUE),
          ifelse(dat$y>=spec$right,pnorm((spec$right-mu)/sd,lower.tail=FALSE,log.p=TRUE),
                 dnorm(dat$y,mu,sd,log=TRUE)))
      }
      info <- NULL
      pred <- function(b,nd) drop(model.matrix(~x+z,nd)%*%b[1:3])
    } else if (kind == "survreg") {
      fixed <- identical(spec$censor,"fixed")
      if(spec$dist=="gaussian") {
        dat$latent <- dat[[paste0("latent",j)]]
        if(spec$censor=="interval") {
          dat$lower <- floor(dat$latent/.4)*.4; dat$upper <- dat$lower+.4
          form <- survival::Surv(lower,upper,type="interval2")~x+z
        } else {
          left <- spec$censor %in% c("left","fixed")
          dat$y <- if(left)pmax(0,dat$latent) else pmin(1,dat$latent)
          dat$event <- if(left)dat$latent>0 else dat$latent<1
          form <- if(left)survival::Surv(y,event,type="left")~x+z else survival::Surv(y,event)~x+z
        }
      } else {
        dat$time <- dat[[paste0(spec$dist,j)]]
        if (spec$censor == "interval") {
          dat$lower <- exp(floor(log(dat$time)/.4)*.4); dat$upper <- dat$lower*exp(.4)
          form <- survival::Surv(lower,upper,type="interval2")~x+z
        } else {
          dat$y <- pmin(2,dat$time); dat$event <- dat$time<=2
          form <- survival::Surv(y,event)~x+z
        }
      }
      m <- fit_model(survival::survreg(form,data=dat,dist=spec$dist,scale=if(fixed).8 else 0,model=TRUE,x=TRUE,y=TRUE))
      par <- if(fixed)coef(m) else c(coef(m),`Log(scale)`=log(m$scale))
      ll <- function(b) {
        mu <- drop(X%*%b[1:3]); s <- if(fixed).8 else exp(b[4])
        if(spec$dist=="gaussian") {
          if(spec$censor=="interval")return(log(pnorm((dat$upper-mu)/s)-pnorm((dat$lower-mu)/s)))
          return(ifelse(dat$event,dnorm(dat$y,mu,s,log=TRUE),
            pnorm((dat$y-mu)/s,lower.tail=spec$censor %in% c("left","fixed"),log.p=TRUE)))
        }
        if (spec$censor == "interval") {
          lo <- (log(dat$lower)-mu)/s; hi <- (log(dat$upper)-mu)/s
          cdf <- switch(spec$dist, lognormal=pnorm, loglogistic=plogis,
            weibull=function(z)-expm1(-exp(z)))
          return(log(cdf(hi)-cdf(lo)))
        }
        if(spec$dist=="weibull")return(ifelse(dat$event,dweibull(dat$y,shape=1/s,scale=exp(mu),log=TRUE),
          pweibull(dat$y,shape=1/s,scale=exp(mu),lower.tail=FALSE,log.p=TRUE)))
        if(spec$dist=="lognormal")return(ifelse(dat$event,dlnorm(dat$y,mu,s,log=TRUE),
          plnorm(dat$y,mu,s,lower.tail=FALSE,log.p=TRUE)))
        zz <- (log(dat$y)-mu)/s
        ifelse(dat$event,dlogis(zz,log=TRUE)-log(s)-log(dat$y),plogis(zz,lower.tail=FALSE,log.p=TRUE))
      }
      info <- NULL
      pred <- function(b,nd) {
        eta <- drop(model.matrix(~x+z,nd)%*%b[1:3])
        if(spec$dist=="gaussian")eta else exp(eta)
      }
    } else stop("Unknown specification")
    stopifnot(length(par)==length(coef(m))+(kind=="zinb")+(kind=="survreg" && !fixed))
    U4 <- numDeriv::jacobian(ll,par,method.args=list(r=4))
    U <- numDeriv::jacobian(ll,par,method.args=list(r=6))
    score_gap <- max(abs(U4-U)/pmax(1,abs(U)))
    total <- function(b)sum(ll(b))
    I4 <- -numDeriv::hessian(total,par,method.args=list(r=4))
    I6 <- -numDeriv::hessian(total,par,method.args=list(r=6))
    derivative_gap <- max(abs(I4-I6)/sqrt(outer(diag(I6),diag(I6))))
    I <- if(is.null(info))I6 else info(par)
    stopifnot(all(is.finite(U)),all(is.finite(I)),min(eigen(I,symmetric=TRUE,only.values=TRUE)$values)>0)
    list(model=m,data=dat,par=par,score=U,information=I,pred=pred,
      loglik_gap=abs(total(par)-as.numeric(logLik(m))),
      stationary=max(abs(colSums(U))/sqrt(diag(I))),derivative_gap=derivative_gap,score_gap=score_gap,
      fitting=list(converged=m$converged,code=m$code,optim_code=m$optim$convergence,
        iterations=if(kind=="truncreg")m$est.stat$nb.iter else m$iter,
        message=if(kind=="truncreg")m$est.stat$message else m$message,
        fail=m$fail,warnings=unique(fitting_warnings)),
      information_kind=if(is.null(info))"independent observed Hessian" else "independent expected beta information")
  }
  nd <- d[,c("x","z")]
  nd0 <- nd1 <- nd; nd0$x <- -.5; nd1$x <- .5
  for (name in names(specs)) for (pattern in names(patterns)) {
    run_case(paste(name,pattern),function() {
      a <- make_component(specs[[name]],patterns[[pattern]][[1]],1)
      b <- make_component(specs[[name]],patterns[[pattern]][[2]],2)
      joint <- suest::suest(a$model,b$model,model_names=c("A","B"),
        observation_id=list(a$data$id,b$data$id),cluster=list(a$data$cluster,b$data$cluster))
      stopifnot(isTRUE(all.equal(unname(coef(joint)),unname(c(a$par,b$par)),tolerance=1e-12)))
      clusters <- sort(unique(c(a$data$cluster,b$data$cluster)))
      aggregate_score <- function(comp) {
        out <- matrix(0,length(clusters),length(comp$par))
        rows <- rowsum(comp$score,comp$data$cluster,reorder=FALSE)
        out[match(as.integer(rownames(rows)),clusters),] <- rows
        out
      }
      q1 <- length(a$par); q2 <- length(b$par)
      inverse <- matrix(0,q1+q2,q1+q2)
      inverse[seq_len(q1),seq_len(q1)] <- solve(a$information)
      inverse[q1+seq_len(q2),q1+seq_len(q2)] <- solve(b$information)
      U <- cbind(aggregate_score(a),aggregate_score(b))
      V <- inverse%*%crossprod(U)%*%t(inverse)*length(clusters)/(length(clusters)-1)
      W <- unname(vcov(joint))
      stopifnot(all(is.finite(V)),all(diag(V)>0))
      out <- list(loglik_gap=max(a$loglik_gap,b$loglik_gap),
        stationary=max(a$stationary,b$stationary),derivative_stability=max(a$derivative_gap,b$derivative_gap),
        score_stability=max(a$score_gap,b$score_gap),fitting=list(A=a$fitting,B=b$fitting),
        covariance=max(abs(W-V)/sqrt(outer(diag(V),diag(V)))),
        reference_crossblock=max(abs(V[seq_len(q1),q1+seq_len(q2)])),
        information=c(a$information_kind,b$information_kind),n=c(nrow(a$data),nrow(b$data)),clusters=length(clusters))
      metrics <- unlist(out[names(bounds)[1:5]])
      stopifnot(identical(names(metrics), names(bounds)[1:5]), all(is.finite(metrics)))
      passed <- all(metrics < bounds[names(metrics)])
      for(kind in c("prediction","comparison")) {
        target <- function(par,comp) if(kind=="prediction")mean(comp$pred(par,nd)) else
          mean(comp$pred(par,nd1)-comp$pred(par,nd0))
        expected <- target(a$par,a)-target(b$par,b)
        g <- c(numDeriv::grad(function(p)target(p,a),a$par),-numDeriv::grad(function(p)target(p,b),b$par))
        se <- sqrt(drop(g%*%V%*%g))
        tab <- if(kind=="prediction")marginaleffects::avg_predictions(joint,newdata=nd,
          numderiv=list("fdcenter",eps=1e-5)) else marginaleffects::avg_comparisons(joint,
          newdata=nd,variables=list(x=c(-.5,.5)),numderiv=list("fdcenter",eps=1e-5))
        stopifnot(nrow(tab)==2,setequal(as.character(tab$group),c("A","B")))
        signs <- ifelse(as.character(tab$group)=="A",1,-1)
        contrast <- marginaleffects::hypotheses(tab,hypothesis=matrix(signs,ncol=1))
        gap <- abs(contrast$estimate-expected); se_gap <- abs(contrast$std.error/se-1)
        out[[kind]] <- c(reference_estimate=expected,reference_se=se,estimate_gap=gap,relative_se_gap=se_gap)
        passed <- passed && gap<bounds["estimate"] && se_gap<bounds["contrast_se"]
      }
      out$status <- if(passed)"PASS" else "FAIL"
      out
    })
  }
  cat("\nSUMMARY\n"); print(vapply(results,`[[`,character(1),"status"))
  cat("\nSESSION INFO\n"); print(sessionInfo())
  cat("\nGATE2_COMPLETED=1\n")
  invisible(results)
}

if (!isTRUE(getOption("suest.gate2.define_only")))
  suest_release_gate2_results <- run_suest_release_gate2()
