# Independent panel/integrated/survey contrast gate. Run from package root.
run_suest_release_gate3 <- function(root=getwd(),output=root,selected=NULL,patterns_selected=NULL) {
  root <- normalizePath(root,mustWork=TRUE)
  dir.create(output,recursive=TRUE,showWarnings=FALSE)
  output <- normalizePath(output,mustWork=TRUE)
  con <- file(file.path(output,"suest_release_gate3.log"),"wt")
  sink(con,split=TRUE); sink(con,type="message")
  on.exit({sink(type="message");sink();close(con)},add=TRUE)
  old <- options(warn=1,digits=15,width=160,survey.lonely.psu="fail",survey.adjust.domain.lonely=FALSE)
  on.exit(options(old),add=TRUE)
  required <- c("pkgload","numDeriv","marginaleffects","plm","nlme","geepack","pglm","glmmTMB","survey","maxLik")
  stopifnot(all(vapply(required,requireNamespace,logical(1),quietly=TRUE)))
  pkgload::load_all(root,quiet=TRUE)
  cat("SUEST RELEASE GATE 3\nStarted:",format(Sys.time(),tz="UTC"),"UTC\n")
  print(sapply(required,function(p)as.character(utils::packageVersion(p))))
  print(tools::md5sum(list.files(file.path(root,"R"),"\\.R$",full.names=TRUE)))
  bounds <- c(covariance=2e-6,gradient=2e-5,estimate=2e-6,se=3e-4,
              continuous_estimate=1e-5,continuous_se=1e-3,
              loglik=2e-6,stationary=2e-3,information_stability=5e-5)
  cat("Preset bounds:\n");print(bounds)
  set.seed(20260929)
  groups <- 160L; periods <- 5L; n <- groups*periods
  d <- data.frame(obs=seq_len(n),id=rep(seq_len(groups),each=periods),time=rep(seq_len(periods),groups))
  d$higher <- ceiling(d$id/4); d$strata <- ceiling(d$higher/10)
  d$x <- runif(n,-1,1); d$z <- runif(n,-1,1)
  d$w <- 1+(d$obs%%11)/4; d$pop <- 20+5*d$strata
  shared <- rnorm(40,sd=.25)[d$higher]
  u1 <- rnorm(groups,sd=.65)[d$id]+shared
  u2 <- rnorm(groups,sd=.5)[d$id]+.35*u1+shared
  e1 <- rnorm(n);e2 <- .35*e1+rnorm(n)
  for(j in 1:2) {
    eta <- if(j==1).15+.6*d$x-.3*d$z else -.2+.35*d$x+.25*d$z
    u <- if(j==1)u1 else u2; e <- if(j==1)e1 else e2
    d[[paste0("gaussian",j)]] <- eta+u+.7*e
    d[[paste0("binary",j)]] <- rbinom(n,1,plogis(eta+u))
    d[[paste0("poisson",j)]] <- rpois(n,exp(eta+.4+u))
  }
  patterns <- list(common=list(rep(TRUE,n),rep(TRUE,n)),
    partial=list(!(d$time==5 & d$id%%3==0),
                 !(d$time==1 & d$id%%4==0 | d$time==4 & d$id%%5==0)),
    disjoint_shared_clusters=list(d$id%%2==0,d$id%%2==1),
    disjoint_clusters=list(d$higher%%2==0,d$higher%%2==1))
  if(!is.null(patterns_selected))patterns<-patterns[patterns_selected]
  specs <- list(fe=list(engine="plm",kind="within",link="identity"),
    be=list(engine="plm",kind="between",link="identity"),
    re=list(engine="plm",kind="random",link="identity"),
    ml=list(engine="lme",link="identity"))
  for(cor in c("independence","exchangeable"))for(link in c("identity","logit","probit","cloglog","log"))
    specs[[paste("gee",link,cor,sep="_")]] <- list(engine="gee",link=link,cor=cor)
  for(link in c("logit","probit","log"))specs[[paste0("pglm_",link)]]<-list(engine="pglm",link=link)
  for(link in c("logit","log"))specs[[paste0("glmm_",link)]]<-list(engine="glmm",link=link)
  for(fpc in c(FALSE,TRUE))for(link in c("identity","logit","probit"))
    specs[[paste("survey",link,if(fpc)"fpc" else "wr",sep="_")]]<-list(engine="survey",link=link,fpc=fpc)
  stopifnot(length(specs)==25L)
  if(!is.null(selected))specs<-specs[selected]
  linkparts <- function(t,link) switch(link,
    identity=list(mu=t,d1=rep(1,length(t)),d2=rep(0,length(t))),
    log=list(mu=exp(t),d1=exp(t),d2=exp(t)),
    logit={p<-plogis(t);a<-p*(1-p);list(mu=p,d1=a,d2=a*(1-2*p))},
    probit=list(mu=pnorm(t),d1=dnorm(t),d2=-t*dnorm(t)),
    cloglog={a<-exp(t-exp(t));list(mu=-expm1(-exp(t)),d1=a,d2=a*(1-exp(t)))})
  family_for <- function(link) if(link=="identity")gaussian() else if(link=="log")poisson() else binomial(link)
  # Golub-Welsch rule for a standard normal: independent of package quadrature.
  normal_rule <- function(k) {
    J <- matrix(0,k,k); J[cbind(1:(k-1),2:k)]<-sqrt(1:(k-1));J<-J+t(J)
    eig<-eigen(J,symmetric=TRUE)
    list(nodes=eig$values,weights=eig$vectors[1,]^2)
  }
  reference <- function(b,spec,nd,weights,kind="prediction",continuous=FALSE) {
    weights<-weights/sum(weights)
    if(kind=="comparison") {
      hi<-lo<-nd;hi$x<-.5;lo$x<--.5
      a<-reference(b,spec,hi,weights,"prediction",continuous)
      z<-reference(b,spec,lo,weights,"prediction",continuous)
      return(list(value=a$value-z$value,gradient=a$gradient-z$gradient))
    }
    X<-model.matrix(~x+z,nd);eta<-drop(X%*%b[1:3]);q<-length(b)
    integrated<-spec$engine%in%c("pglm","glmm")&&spec$link%in%c("logit","probit")
    sigma<-if(integrated)if(spec$engine=="glmm")exp(b[4]) else b[4] else 0
    ds<-if(spec$engine=="glmm")sigma else 1
    at_node <- function(z) {
      part<-linkparts(eta+sigma*z,spec$link)
      if(kind=="slope") {
        val<-b[2]*sum(weights*part$d1)
        grad<-c(b[2]*colSums(X*(weights*part$d2)),rep(0,q-3))
        grad[2]<-grad[2]+sum(weights*part$d1)
        if(integrated)grad[4]<-b[2]*sum(weights*part$d2*z)*ds
      } else {
        val<-sum(weights*part$mu)
        grad<-c(colSums(X*(weights*part$d1)),rep(0,q-3))
        if(integrated)grad[4]<-sum(weights*part$d1*z)*ds
      }
      c(val,grad)
    }
    if(integrated) {
      if(continuous) {
        if(spec$link=="probit") {
          scale<-sqrt(1+sigma^2);v<-eta/scale;phi<-dnorm(v)
          if(kind=="prediction") {
            ans<-c(sum(weights*pnorm(v)),colSums(X*(weights*phi/scale)),
              sum(weights*(-eta*sigma*phi/scale^3))*ds)
          } else {
            slope<-phi/scale
            gb<-b[2]*colSums(X*(weights*(-eta*phi/scale^3)))
            gb[2]<-gb[2]+sum(weights*slope)
            gs<-b[2]*sum(weights*sigma*phi*(eta^2/scale^5-1/scale^3))*ds
            ans<-c(b[2]*sum(weights*slope),gb,gs)
          }
        } else ans<-vapply(seq_len(q+1),function(k)integrate(function(z)
          vapply(z,function(zz)at_node(zz)[k]*dnorm(zz),numeric(1)),
          -Inf,Inf,rel.tol=1e-10,abs.tol=1e-11)$value,numeric(1))
      } else {
        rule<-normal_rule(if(spec$engine=="pglm")12L else 20L)
        ans<-drop(vapply(rule$nodes,at_node,numeric(q+1))%*%rule$weights)
      }
    } else if(spec$engine=="glmm" && spec$link=="log") {
      variance<-exp(2*b[4]);mu<-exp(eta+variance/2)
      grad<-c(colSums(X*(weights*mu)),sum(weights*mu)*variance)
      val<-sum(weights*mu)
      if(kind=="slope") {grad<-b[2]*grad;grad[2]<-grad[2]+val;val<-b[2]*val}
      ans<-c(val,grad)
    } else ans<-at_node(0)
    list(value=unname(ans[1]),gradient=unname(ans[-1]))
  }
  fit_component <- function(spec,keep,j,design=NULL) {
    dat<-d[keep,];dat$y<-dat[[paste0(if(spec$link=="identity")"gaussian" else if(spec$link=="log")"poisson" else "binary",j)]]
    warnings<-character();warning_trials<-list()
    m<-withCallingHandlers({
      if(spec$engine=="plm")plm::plm(y~x+z,data=dat,index=c("id","time"),model=spec$kind,random.method="swar")
      else if(spec$engine=="lme")nlme::lme(y~x+z,random=~1|id,data=dat,method="ML")
      else if(spec$engine=="gee")geepack::geeglm(y~x+z,id=id,data=dat,family=family_for(spec$link),corstr=spec$cor,
        control=geepack::geese.control(epsilon=1e-10,maxit=200))
      else if(spec$engine=="glmm")glmmTMB::glmmTMB(y~x+z+(1|id),data=dat,family=family_for(spec$link))
      else if(spec$engine=="pglm") {
        maxLik<-maxLik::maxLik
        start<-c(coef(glm(y~x+z,data=dat,family=family_for(spec$link))),if(spec$link=="log")c(alpha=.6) else c(sigma=.7))
        pglm::pglm(y~x+z,data=dat,family=family_for(spec$link),model="random",effect="individual",index=c("id","time"),
          R=12,method="bfgs",print.level=0,reltol=1e-12,start=start,other=if(spec$link=="log")"sd" else NULL)
      } else {
        form<-reformulate(c("x","z"),paste0(if(spec$link=="identity")"gaussian" else "binary",j))
        survey::svyglm(form,design=design,subset=keep,
          family=if(spec$link=="identity")gaussian() else quasibinomial(spec$link),
          control=list(epsilon=1e-12,maxit=100),rescale=(j==1))
      }
    },warning=function(w) {
      warnings<<-c(warnings,paste(deparse(conditionCall(w)),conditionMessage(w),collapse=" "))
      # Capture the offending native likelihood trial, before its stack unwinds.
      frames<-Filter(function(e)all(vapply(c("param","X","y","id","other"),
        exists,logical(1),envir=e,inherits=FALSE)),sys.frames())
      trial<-NA_real_
      if(length(frames)) {
        frame<-frames[[length(frames)]]
        if(identical(get("other",frame),"sd"))trial<-tail(get("param",frame),1)
      }
      warning_trials[[length(warning_trials)+1L]]<<-list(message=conditionMessage(w),alpha=unname(trial))
    })
    X<-model.matrix(~x+z,dat);y<-dat$y
    b<-if(spec$engine=="lme")c(nlme::fixef(m),sigma_u=sqrt(as.numeric(nlme::getVarCov(m)[1,1])),sigma_e=m$sigma)
       else if(spec$engine=="glmm")c(glmmTMB::fixef(m)$cond,log_sigma=log(attr(glmmTMB::VarCorr(m)$cond[[1]],"stddev")[1]))
       else if(spec$engine=="plm"&&spec$kind=="within")c(`(Intercept)`=mean(y)-sum(colMeans(X[,2:3])*coef(m)),coef(m))
       else coef(m)
    if(spec$engine=="pglm" && spec$link!="log")b[4]<-abs(b[4])
    likelihood_check<-NULL
    if(spec$engine=="pglm" && spec$link=="log") {
      panels<-split(seq_len(nrow(dat)),dat$id)
      ll<-function(par) {
        mu<-exp(drop(X%*%par[1:3]));shape<-1/par[4]
        sum(vapply(panels,function(i) {
          Y<-sum(y[i]);M<-sum(mu[i])
          sum(y[i]*log(mu[i])-lgamma(y[i]+1))+lgamma(shape+Y)-lgamma(shape)+
            shape*log(shape)-(shape+Y)*log(shape+M)
        },numeric(1)))
      }
      I4<- -numDeriv::hessian(ll,b,method.args=list(r=4))
      I6<- -numDeriv::hessian(ll,b,method.args=list(r=6))
      stopifnot(all(is.finite(I6)),min(eigen(I6,symmetric=TRUE,only.values=TRUE)$values)>0)
      native_warnings<-character()
      native<-withCallingHandlers(m$objectiveFn(coef(m)),warning=function(w)
        native_warnings<<-c(native_warnings,conditionMessage(w)))
      likelihood_check<-list(loglik_gap=abs(ll(b)-as.numeric(logLik(m))),
        stationary=max(abs(numDeriv::grad(ll,b))/sqrt(diag(I6))),
        information_stability=max(abs(I4-I6)/sqrt(outer(diag(I6),diag(I6)))),
        native_final_finite=all(is.finite(c(native,attr(native,"gradient"),attr(native,"hessian")))),
        native_final_warnings=native_warnings,
        all_warning_trials_negative_alpha=length(warning_trials)>0&&all(vapply(warning_trials,
          function(x)is.finite(x$alpha)&&x$alpha<0&&identical(x$message,"NaNs produced"),logical(1))))
      stopifnot(likelihood_check$loglik_gap<bounds["loglik"],likelihood_check$stationary<bounds["stationary"],
        likelihood_check$information_stability<bounds["information_stability"],likelihood_check$native_final_finite,
        !length(native_warnings))
    }
    influence<-NULL
    if(spec$engine=="plm") {
      idx<-split(seq_len(nrow(dat)),dat$id)
      means<-function(A) {if(is.null(dim(A)))A<-matrix(A,ncol=1);t(vapply(idx,function(i)colMeans(A[i,,drop=FALSE]),numeric(ncol(A))))}
      Xbar<-means(X);ybar<-drop(means(y));repeated<-match(as.character(dat$id),names(idx))
      if(spec$kind=="within") {
        T<-X[,2:3]-Xbar[repeated,2:3];r<-y-ybar[repeated]-drop(T%*%b[2:3])
        influence<-(T*r)%*%solve(crossprod(T))%*%t(rbind(-colMeans(X[,2:3]),diag(2)))
        df_n<-nrow(dat);df_k<-3
      } else if(spec$kind=="between") {
        r<-ybar-drop(Xbar%*%b);Ip<-(Xbar*r)%*%solve(crossprod(Xbar))
        influence<-matrix(0,nrow(dat),3);influence[vapply(idx,`[`,integer(1),1L),]<-Ip
        df_n<-nrow(Xbar);df_k<-3
      } else {
        theta<-m$ercomp$theta
        if(length(theta)>1)theta<-theta[match(as.character(dat$id),names(theta))]
        T<-X-Xbar[repeated,,drop=FALSE]*as.numeric(theta)
        r<-y-ybar[repeated]*as.numeric(theta)-drop(T%*%b)
        influence<-(T*r)%*%solve(crossprod(T));df_n<-nrow(dat);df_k<-3
      }
      G<-length(unique(dat$higher));influence<-influence*sqrt(G/(G-1)*(df_n-1)/(df_n-df_k))
    } else if(spec$engine=="gee") {
      part<-linkparts(drop(X%*%b),spec$link);v<-if(spec$link=="identity")rep(1,nrow(dat)) else if(spec$link=="log")part$mu else part$mu*(1-part$mu)
      rho<-if(spec$cor=="independence")0 else unname(m$geese$alpha)
      phi<-unname(m$geese$gamma);U<-matrix(0,nrow(dat),3);I<-matrix(0,3,3)
      for(rows in split(seq_len(nrow(dat)),dat$id)) {
        R<-matrix(rho,length(rows),length(rows));diag(R)<-1
        W<-solve(phi*R*sqrt(outer(v[rows],v[rows])))
        D<-X[rows,,drop=FALSE]*part$d1[rows]
        U[rows,]<-D*drop(W%*%(y[rows]-part$mu[rows]))
        I<-I+crossprod(D,W%*%D)
      }
      G<-length(unique(dat$higher));influence<-U%*%solve(I)*sqrt(G/(G-1))
    } else if(spec$engine=="survey") {
      part<-linkparts(drop(X%*%b),spec$link)
      v<-if(spec$link=="identity")rep(1,nrow(dat)) else part$mu*(1-part$mu)
      U<-X*(dat$w*(y-part$mu)*part$d1/v)
      I<-crossprod(X,X*(dat$w*part$d1^2/v))
      influence<-U%*%solve(I)
    }
    list(model=m,data=dat,parameters=b,influence=influence,
      fitting=list(warnings=unique(warnings),warning_trials=warning_trials,likelihood_check=likelihood_check,converged=m$converged,
        code=if(spec$engine=="glmm")m$fit$convergence else if(spec$engine=="gee")m$geese$error else m$code))
  }
  results<-list()
  for(name in names(specs))for(pattern in names(patterns)) {
    label<-paste(name,pattern);cat("\nCASE:",label,"\n")
    warnings<-character()
    ans<-tryCatch(withCallingHandlers({
      spec<-specs[[name]];keep<-patterns[[pattern]]
      design<-if(spec$engine=="survey") {
        if(spec$fpc)survey::svydesign(~higher,strata=~strata,weights=~w,fpc=~pop,data=d) else
          survey::svydesign(~higher,strata=~strata,weights=~w,data=d)
      } else NULL
      a<-fit_component(spec,keep[[1]],1,design);b<-fit_component(spec,keep[[2]],2,design)
      joint<-if(spec$engine=="survey")suest::suest(a$model,b$model,model_names=c("A","B"),observation_id="obs",survey_design=design) else
        suest::suest(a$model,b$model,model_names=c("A","B"),observation_id=list(a$data$obs,b$data$obs),cluster=list(a$data$higher,b$data$higher))
      pa<-a$parameters;pb<-b$parameters;pars<-c(pa,pb);q<-length(pa)
      stopifnot(length(pb)==q,isTRUE(all.equal(unname(coef(joint)),unname(pars),tolerance=1e-10)))
      V<-unname(vcov(joint));covgap<-NA_real_
      if(!is.null(a$influence)) {
        align<-function(comp) {out<-matrix(0,n,ncol(comp$influence));out[comp$data$obs,]<-comp$influence;out}
        U<-cbind(align(a),align(b))
        if(spec$engine=="survey") {
          Vref<-matrix(0,2*q,2*q)
          for(h in unique(d$strata)) {
            rows<-d$strata==h;S<-rowsum(U[rows,,drop=FALSE],d$higher[rows]);G<-nrow(S)
            S<-sweep(S,2,colMeans(S));fpc<-if(spec$fpc)1-G/d$pop[which(rows)[1]] else 1
            Vref<-Vref+crossprod(S)*G/(G-1)*fpc
          }
        } else Vref<-crossprod(rowsum(U,d$higher))
        covgap<-max(abs(V-Vref)/sqrt(outer(diag(Vref),diag(Vref))))
        V<-Vref
      }
      stopifnot(all(is.finite(V)),all(diag(V)>0))
      out<-list(covariance_reference=if(is.na(covgap))"audited suest covariance; new check is conditional propagation" else
          if(spec$engine=="survey")"independent expected-information influence and PSU algebra" else "independent transformed/mean equations; native working covariance parameters",
        covariance_gap=covgap,fitting=list(A=a$fitting,B=b$fitting),n=c(nrow(a$data),nrow(b$data)),
        reference_crossblock=max(abs(V[1:q,q+1:q])))
      passed<-is.na(covgap)||covgap<bounds["covariance"]
      w<-if(spec$engine=="survey")d$w else rep(1,n)
      # A fixed evaluation population away from the FE training means avoids
      # the structurally zero conditional variance of its average fitted mean.
      nd<-d[,c("x","z")];nd$x<-nd$x+.25;nd$z<-nd$z+.15
      is_integrated<-spec$engine%in%c("pglm","glmm")&&spec$link%in%c("logit","probit")
      for(kind in c("prediction","comparison","slope")) {
        ra<-reference(pa,spec,nd,w,kind);rb<-reference(pb,spec,nd,w,kind)
        expected<-ra$value-rb$value;g<-c(ra$gradient,-rb$gradient)
        numerical<-numDeriv::grad(function(p)reference(p[1:q],spec,nd,w,kind)$value-
          reference(p[q+1:q],spec,nd,w,kind)$value,pars)
        gradient_gap<-max(abs(g-numerical)/pmax(1,abs(g)))
        se<-sqrt(drop(g%*%V%*%g));stopifnot(is.finite(se),se>0)
        args<-list(model=joint,newdata=nd,wts=w,numderiv=list("fdcenter",eps=1e-5))
        tab<-if(kind=="prediction")do.call(marginaleffects::avg_predictions,args) else if(kind=="comparison")
          do.call(marginaleffects::avg_comparisons,c(args,list(variables=list(x=c(-.5,.5))))) else
          do.call(marginaleffects::avg_slopes,c(args,list(variables="x",eps=1e-5)))
        stopifnot(nrow(tab)==2,setequal(as.character(tab$group),c("A","B")))
        signs<-ifelse(as.character(tab$group)=="A",1,-1)
        contrast<-marginaleffects::hypotheses(tab,hypothesis=matrix(signs,ncol=1))
        gap<-abs(contrast$estimate-expected);se_gap<-abs(contrast$std.error/se-1)
        gn<-g;if(q>3)gn[c(4:q,q+4:q)]<-0
        Vn<-V;Vn[1:q,q+1:q]<-Vn[q+1:q,1:q]<-0
        metric<-c(reference_estimate=expected,reference_se=se,gradient_gap=gradient_gap,
          estimate_gap=gap,relative_se_gap=se_gap,
          se_change_omit_nuisance=sqrt(drop(gn%*%V%*%gn))/se-1,
          se_change_omit_crosscov=sqrt(drop(g%*%Vn%*%g))/se-1)
        passed<-passed&&gradient_gap<bounds["gradient"]&&gap<bounds["estimate"]&&se_gap<bounds["se"]
        if(is_integrated) {
          ca<-reference(pa,spec,nd,w,kind,TRUE);cb<-reference(pb,spec,nd,w,kind,TRUE)
          cg<-c(ca$gradient,-cb$gradient);cs<-sqrt(drop(cg%*%V%*%cg))
          ce<-abs(expected-(ca$value-cb$value));cr<-abs(se/cs-1)
          metric<-c(metric,continuous_estimate_gap=ce,continuous_relative_se_gap=cr)
          passed<-passed&&ce<bounds["continuous_estimate"]&&cr<bounds["continuous_se"]
        }
        out[[kind]]<-metric
      }
      out$reviewed_warning_count<-if(spec$engine=="pglm"&&spec$link=="log"&&
        a$fitting$likelihood_check$all_warning_trials_negative_alpha&&
        b$fitting$likelihood_check$all_warning_trials_negative_alpha)
          length(a$fitting$warning_trials)+length(b$fitting$warning_trials) else 0L
      out$status<-if(passed)"PASS" else "FAIL";out
    },warning=function(w)warnings<<-c(warnings,conditionMessage(w))),error=function(e)list(status="ERROR",error=conditionMessage(e)))
    ans$warnings<-unique(warnings)
    if(length(warnings)&&identical(ans$status,"PASS"))ans$status<-
      if(identical(length(warnings),ans$reviewed_warning_count))"PASS_REVIEWED_FITTING_WARNING" else "REVIEW_WARNING"
    results[[label]]<-ans;print(ans);saveRDS(results,file.path(output,"suest_release_gate3_results.rds"))
  }
  cat("\nSUMMARY\n");print(vapply(results,`[[`,character(1),"status"))
  print(sessionInfo());cat("\nGATE3_COMPLETED=1\n")
  invisible(results)
}
if(!isTRUE(getOption("suest.gate3.define_only")))suest_release_gate3_results<-run_suest_release_gate3()
