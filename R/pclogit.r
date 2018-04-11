pclogit=function(x,y,stra=NULL,alpha=1.0,nlam=100,lambda=NULL,eps=NULL,maxit=100000,maxp=p,thre=1e-6,group=NULL,type=c("ring","fcon"),sgnc=NULL){

  type <- match.arg(type)
  x<-as.matrix(x)
  n=as.integer(nrow(x))
  p=as.integer(ncol(x))
  if (length(y)!=n) stop("x and y have different number of rows.")
  y[y!=0]<-1
  if (is.null(stra)) {
    sy<-sum(y)/length(y)
    if (sy<1e-5) stop("the proportion of cases is too small.")
    else if (sy>(1-1e-5)) stop("the proportion of controls is too small.")
    stra<-rep(0,n)
    }
  else {
    if (length(stra)!=n) stop("strata must match the rows of x.")
    if (!identical(as.integer(sort(unique(stra))),as.integer(seq(max(stra))))) stop("strata must be consecutively numbered 1,2,3,...")
    ytab<-table(as.integer(y),as.integer(stra))
    if (sum(as.numeric(ytab[y="1",])!=1)>0) stop("each stratum must have exact one case.")
    if (min(as.numeric(ytab[y="0",]))<1) stop("each stratum must have at least one control.")
    }
  if (is.null(group)) {
    loc<-matrix(0,2,p)
    idg<-rep(0,p)
    }
  else {
    if (!is.matrix(group)) {
      group<-as.integer(group)
      if (any(group<=0)) stop("group must be postive integers.")
      if (sum(group)!=p) stop("sum of group should be equal to p.")    
      }
    if (!is.null(sgnc)) {
      sgnc=sign(sgnc)
      if (length(sgnc)!=p) stop("the length of sgnc should be equal to p.")
    }     
    gr<-group.matrix(group,type,sgnc)
    loc<-gr$loc
    idg<-gr$idg  
    }
  if (is.null(colnames(x))) colnames(x) <- paste("V",1:p,sep="")
  if (is.null(lambda)){
    if (is.null(eps)) flmin=as.double(ifelse(n<p,5e-2,1e-4))
    else flmin=as.double(max(eps,1e-6))
    ulam=double(1)
    nlam=as.integer(nlam)
  }
  else {
    flmin=as.double(1)
    if (any(lambda<0)) stop("lambdas must be non-negative")
    ulam=as.double(rev(sort(lambda)))
    nlam=as.integer(length(lambda))
  }
  nx=as.integer(min(maxp,p))
  y=as.factor(y)
  ntab=table(y)
  nc=as.integer(length(ntab))
  y=diag(nc)[as.numeric(y),]
  stab<-as.numeric(table(stra))
  q0<-rep(1/stab,stab)
  lq<-as.integer(nrow(loc))
  fit=.Fortran("pelogit",as.double(alpha),n,p,as.double(x),y,nx,nlam,flmin,ulam,as.double(thre),as.integer(maxit),
               lmu=integer(1),a0=double(nlam),ca=double(nx*nlam),ia=integer(nx),nin=integer(nlam),nulldev=double(1),
               dev=double(nlam),alm=double(nlam),nlp=integer(1),jerr=integer(1),as.integer(loc),as.double(idg),lq,
               as.integer(stra),as.double(q0))
  jerr<-fit$jerr
  if(jerr>0) {
        if (jerr<7777) msg="Memory allocation error"
        else if (jerr==7777) msg="All used predictors have zero variance"
        else msg="Unknown error"
        stop(paste("from pclogit Fortran code - ",msg),call.=FALSE)
        }
  if (jerr<0) {
        if (jerr> -10000) msg=paste("Convergence for ",-jerr,"th lambda value not reached after maxit=",maxit," iterations; solutions for larger lambdas returned",sep="")
        if (jerr< -10000) msg=paste("Number of nonzero coefficients along the path exceeds maxp=",nx, " at ",-jerr-10000,"th lambda value; solutions for larger lambdas returned",sep="")
        warning(paste("from pclogit Fortran code - ",msg),call.=FALSE)
        }
  lmu=fit$lmu
  nin=fit$nin[seq(lmu)]
  ninmax=max(nin)
  lam=fit$alm[seq(lmu)]
  dd=c(p,lmu)
  cnames=list(colnames(x),paste("s",seq(lmu)-1,sep=""))
  if (ninmax>0) {
        ca=matrix(fit$ca[seq(nx*lmu)],nx,lmu)[seq(ninmax),,drop=FALSE]
        df=apply(abs(ca)>0,2,sum)
        ja=fit$ia[seq(ninmax)]
        oja=order(ja)
        ja=rep(ja[oja],lmu)
        ia=cumsum(c(1,rep(ninmax,lmu)))
        ca=ca[oja,]
       }
  else {
        ca=rep(0,lmu)
        df=rep(0,lmu)
        ja=rep(1,lmu)
        ia=seq(lmu+1)
       }
  beta=new("dgCMatrix",Dim=dd,Dimnames=cnames,x=as.vector(ca),p=as.integer(ia-1),i=as.integer(ja-1))
  if (is.null(lambda)) {
        llam=log(lam)
        lam[1]=exp(2*llam[2]-llam[3])
      }
  dev=fit$dev[seq(fit$lmu)]    
  results<-list(beta=-beta,lambda=lam,df=df,nobs=n,alpha=alpha,iterations=fit$nlp,jerr=jerr)
  if (sum(stra)<1) results<-c(list(b0=as.numeric(-fit$a0)),results)
  else results<-c(list(strata=as.integer(stra)),results)
  class(results)<-"pclogit"
  return(results)
}
