# Remaining release checks. Run from the candidate package root.
run_suest_release_gate4 <- function(root=getwd(),output=root,selected=NULL,patterns_selected=NULL) {
 dir.create(output,recursive=TRUE,showWarnings=FALSE)
 con<-file(file.path(output,'suest_release_gate4.log'),'wt');sink(con,split=TRUE);sink(con,type='message')
 on.exit({sink(type='message');sink();close(con)},add=TRUE)
 pkgload::load_all(root,quiet=TRUE);options(digits=15,warn=1)
 cat('GATE 4; fixed bounds: engine b/est 2e-5, V 2e-4, SE 5e-4; independent V 2e-5, est 2e-6, SE 3e-4\n')
 set.seed(29264);n<-480;d<-data.frame(obs=1:n,cluster=rep(1:60,each=8),x=runif(n,-1,1),z=rnorm(n),iv=rnorm(n))
 u<-rnorm(60,sd=.3)[d$cluster];e<-rnorm(n)
 d$endog<-.6*d$iv+.3*d$x+.4*e+rnorm(n)
 for(j in 1:2) {
  eta<-(if(j==1).2 else -.15)+(if(j==1).6 else .35)*d$x+.2*d$z+u
  d[[paste0('bin',j)]]<-rbinom(n,1,plogis(eta))
  d[[paste0('count',j)]]<-rpois(n,exp(eta))
  d[[paste0('ord',j)]]<-ordered(cut(eta+rlogis(n),c(-Inf,-.5,.6,Inf),labels=c('low','mid','high')))
  ee<-mvtnorm::rmvnorm(n,sigma=matrix(c(1,.35,.35,1),2))
  d[[paste0('a',j)]]<-as.integer(eta+ee[,1]>0)
  d[[paste0('b',j)]]<-as.integer(-.1+.4*d$x-.15*d$z+u+ee[,2]>0)
  d[[paste0('linear',j)]]<-eta+(if(j==1).8 else .5)*d$endog+e+rnorm(n)
 }
 patterns<-list(common=list(rep(TRUE,n),rep(TRUE,n)),partial=list(d$obs%%5!=0,d$obs%%7!=0),
   disjoint_shared=list(d$obs%%2==0,d$obs%%2==1),disjoint_clusters=list(d$cluster%%2==0,d$cluster%%2==1))
 if(!is.null(patterns_selected))patterns<-patterns[patterns_selected]
 specs<-c('glm_logit','glm_probit','glm_poisson','ordinal_logit','ordinal_probit','biprobit','iv')
 if(!is.null(selected))specs<-intersect(specs,selected)
 nd<-d[1:40,c('x','z','endog','iv')];nd$x<-nd$x+.2
 bv<-function(a,b,r) vapply(seq_along(a),function(i)as.numeric(mvtnorm::pmvnorm(upper=c(a[i],b[i]),corr=matrix(c(1,r,r,1),2),algorithm=mvtnorm::TVPACK())),numeric(1))
 target<-function(b,nd,kind,engine) {
  if(kind=='comparison') {hi<-lo<-nd;hi$x<-.5;lo$x<--.5;return(target(b,hi,'prediction',engine)-target(b,lo,'prediction',engine))}
  if(engine=='iv') {X<-cbind(1,nd$endog,nd$x,nd$z);return(if(kind=='slope')unname(b[3]) else mean(X%*%b))}
  X<-model.matrix(~x+z,nd);a<-drop(X%*%b[1:3]);c<-drop(X%*%b[4:6]);r<-tanh(b[7]);s<-sqrt(1-r*r)
  if(kind=='prediction')mean(bv(a,c,r)) else mean(b[2]*dnorm(a)*pnorm((c-r*a)/s)+b[5]*dnorm(c)*pnorm((a-r*c)/s))
 }
 predgrad<-function(b,nd) {X<-model.matrix(~x+z,nd);a<-drop(X%*%b[1:3]);c<-drop(X%*%b[4:6]);r<-tanh(b[7]);s<-sqrt(1-r*r)
  c(colMeans(X*(dnorm(a)*pnorm((c-r*a)/s))),colMeans(X*(dnorm(c)*pnorm((a-r*c)/s))),mean(exp(-(a*a-2*r*a*c+c*c)/(2*s*s))/(2*pi*s))*s*s)
 }
 fitone<-function(spec,dat,j,alt=FALSE) {
  f<-function(y)reformulate(c('x','z'),y)
  if(grepl('^glm_',spec)) {link<-sub('glm_','',spec);fam<-if(link=='poisson')poisson() else binomial(link);fn<-if(alt)glm2::glm2 else stats::glm
   return(fn(f(paste0(if(link=='poisson')'count' else 'bin',j)),data=dat,family=fam,control=glm.control(epsilon=1e-11,maxit=100)))}
  if(grepl('^ordinal_',spec)) {link<-sub('ordinal_','',spec)
   if(alt)return(ordinal::clm(f(paste0('ord',j)),data=dat,link=link,control=ordinal::clm.control(gradTol=1e-9)))
   return(MASS::polr(f(paste0('ord',j)),data=dat,method=if(link=='logit')'logistic' else link,Hess=TRUE,control=list(reltol=1e-12,maxit=1000)))}
  if(spec=='biprobit') {form<-as.formula(paste0('cbind(a',j,',b',j,')~x+z'));trials<-list();m<-withCallingHandlers(mvProbit::mvProbit(form,data=dat,algorithm=mvtnorm::TVPACK(),method='BFGS',finalHessian=TRUE,iterlim=500,reltol=1e-11,printLevel=0),warning=function(w) {
      frames<-sys.frames();vals<-lapply(frames,function(fr) {if(exists('coef',fr,inherits=FALSE)){v<-get('coef',fr,inherits=FALSE);if(is.list(v)&&!is.null(v$sigma))return(list(beta=v$beta,rho=v$sigma[1,2]))};NULL})
      trials[[length(trials)+1]]<<-list(message=conditionMessage(w),call=deparse(conditionCall(w)),parameters=Filter(Negate(is.null),vals))
    });m$call$formula<-form;attr(m,'gate_warning_trials')<-trials;return(m)}
  fixest::feols(as.formula(paste0('linear',j,'~x+z|endog~iv')),data=dat)
 }
 tabs<-function(joint,kind) {
  args<-list(model=joint,newdata=nd,numderiv=list('fdcenter',eps=1e-5))
  if(kind=='prediction')do.call(marginaleffects::avg_predictions,args) else if(kind=='comparison')do.call(marginaleffects::avg_comparisons,c(args,list(variables=list(x=c(-.5,.5))))) else do.call(marginaleffects::avg_slopes,c(args,list(variables='x',eps=1e-5)))
 }
 results<-list()
 for(spec in specs)for(pat in names(patterns)) {
  label<-paste(spec,pat);cat('\nCASE',label,'\n');warnings<-character()
  ans<-tryCatch(withCallingHandlers({
   ds<-lapply(patterns[[pat]],function(k)d[k,]);mods<-lapply(1:2,function(j)fitone(spec,ds[[j]],j))
   combine<-function(mm)suest::suest(mm[[1]],mm[[2]],model_names=c('A','B'),observation_id=lapply(ds,`[[`,'obs'),cluster=lapply(ds,`[[`,'cluster'))
   if(spec=='biprobit')saveRDS(mods,file.path(output,paste0('fits-',pat,'.rds')))
   fit<-combine(mods);V<-unname(vcov(fit));p<-unname(coef(fit));q<-length(p)/2
   out<-list();passed<-TRUE
   if(spec=='biprobit')out$fitting<-lapply(mods,function(m)list(code=m$code,maximum=m$maximum,warning_trials=attr(m,'gate_warning_trials')))
   if(grepl('^(glm|ordinal)_',spec)) {
    alt<-combine(lapply(1:2,function(j)fitone(spec,ds[[j]],j,TRUE)))
    stopifnot(setequal(names(coef(fit)),names(coef(alt))))
    order_alt<-match(names(coef(fit)),names(coef(alt)))
    out$coefficient_gap<-max(abs(coef(fit)-coef(alt)[order_alt]));out$covariance_gap<-max(abs(V-vcov(alt)[order_alt,order_alt])/sqrt(outer(diag(V),diag(V))))
    passed<-out$coefficient_gap<2e-5&&out$covariance_gap<2e-4
    for(kind in c('prediction','comparison','slope')) {
     t1<-tabs(fit,kind);t2<-tabs(alt,kind)
     # Use group/term labels from output: one contrast per response category.
     groups<-as.character(t1$group);stopifnot(identical(groups,as.character(t2$group)))
     ka<-which(grepl('^A($|:)',groups));kb<-which(grepl('^B($|:)',groups))
     if(!length(ka)){ka<-which(groups=='A');kb<-which(groups=='B')}
     stopifnot(length(ka)>0,length(ka)==length(kb))
     H<-matrix(0,nrow(t1),length(ka));H[cbind(ka,seq_along(ka))]<-1;H[cbind(kb,seq_along(kb))]<--1
     a<-marginaleffects::hypotheses(t1,hypothesis=H);b<-marginaleffects::hypotheses(t2,hypothesis=H)
     met<-c(estimate_gap=max(abs(a$estimate-b$estimate)),relative_se_gap=max(abs(a$std.error/b$std.error-1)))
     out[[kind]]<-met;passed<-passed&&all(is.finite(met))&&met[1]<2e-5&&met[2]<5e-4
    }
   } else {
    U<-matrix(0,n,2*q)
    for(j in 1:2) {
     dat<-ds[[j]];b<-p[(j-1)*q+seq_len(q)]
     if(spec=='iv') {X<-cbind(1,dat$endog,dat$x,dat$z);Z<-cbind(1,dat$iv,dat$x,dat$z);Xh<-qr.fitted(qr(Z),X)
      influence<-(Xh*drop(dat[[paste0('linear',j)]]-X%*%b))%*%solve(crossprod(Xh,X))
     } else {
      X<-model.matrix(~x+z,dat);s1<-2*dat[[paste0('a',j)]]-1;s2<-2*dat[[paste0('b',j)]]-1
      ll<-function(bb) {a<-drop(X%*%bb[1:3]);c<-drop(X%*%bb[4:6]);r<-tanh(bb[7]);vapply(seq_along(a),function(i)log(bv(s1[i]*a[i],s2[i]*c[i],s1[i]*s2[i]*r)),numeric(1))}
      scores<-numDeriv::jacobian(ll,b);info<--numDeriv::hessian(function(bb)sum(ll(bb)),b)
      stopifnot(min(eigen(info,symmetric=TRUE,only.values=TRUE)$values)>0)
      influence<-scores%*%solve(info)*sqrt(length(unique(d$cluster[patterns[[pat]][[1]]|patterns[[pat]][[2]]]))/(length(unique(d$cluster[patterns[[pat]][[1]]|patterns[[pat]][[2]]]))-1))
      out[[paste0('score_sum_',j)]]<-max(abs(colSums(scores)))
     }
     U[dat$obs,(j-1)*q+seq_len(q)]<-influence
    }
    Vr<-crossprod(rowsum(U,d$cluster));out$covariance_gap<-max(abs(V-Vr)/sqrt(outer(diag(Vr),diag(Vr))))
    passed<-out$covariance_gap<2e-5
    for(kind in c('prediction','comparison','slope')) {
     value<-function(pp)target(pp[1:q],nd,kind,spec)-target(pp[q+1:q],nd,kind,spec)
     g<-numDeriv::grad(value,p)
     if(spec=='iv') {ga<-if(kind=='prediction')colMeans(cbind(1,nd$endog,nd$x,nd$z)) else c(0,0,1,0);ga<-c(ga,-ga);stopifnot(max(abs(g-ga))<2e-6);g<-ga}
     if(spec=='biprobit'&&kind!='slope') {
      grad<-function(b)if(kind=='prediction')predgrad(b,nd) else {hi<-lo<-nd;hi$x<-.5;lo$x<--.5;predgrad(b,hi)-predgrad(b,lo)}
      ga<-c(grad(p[1:q]),-grad(p[q+1:q]));out[[paste0(kind,'_gradient_gap')]]<-max(abs(g-ga));stopifnot(max(abs(g-ga))<2e-6);g<-ga
     }
     expected<-value(p);se<-sqrt(drop(g%*%Vr%*%g));tab<-tabs(fit,kind);stopifnot(nrow(tab)==2)
     a<-marginaleffects::hypotheses(tab,hypothesis=matrix(ifelse(tab$group=='A',1,-1),ncol=1))
     gn<-g;if(spec=='biprobit')gn[c(q,2*q)]<-0
     vz<-Vr;vz[1:q,q+1:q]<-vz[q+1:q,1:q]<-0
     met<-c(estimate_gap=abs(a$estimate-expected),relative_se_gap=abs(a$std.error/se-1),reference_se=se,
       omitted_rho_se_change=sqrt(drop(gn%*%Vr%*%gn))/se-1,omitted_cross_se_change=sqrt(drop(g%*%vz%*%g))/se-1)
     out[[kind]]<-met;passed<-passed&&all(is.finite(met))&&met[1]<2e-6&&met[2]<3e-4
    }
   }
   out$status<-if(passed)'PASS' else 'FAIL';out
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e)list(status='ERROR',error=conditionMessage(e)))
  ans$warnings<-warnings;if(length(warnings)&&ans$status=='PASS')ans$status<-'REVIEW_WARNING'
  results[[label]]<-ans;print(ans);saveRDS(results,file.path(output,'suest_release_gate4_results.rds'))
 }
 print(vapply(results,`[[`,character(1),'status'));print(sessionInfo());cat('GATE4_COMPLETED=1\n');invisible(results)
}
if(!isTRUE(getOption('suest.gate4.define_only')))run_suest_release_gate4()
