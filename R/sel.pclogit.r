sel.pclogit<-function(x,y,stra=NULL,...,psub=0.5,N.lam=5,K=100) {
  
   g<-pclogit(x,y,stra,...)
   if (psub<0.5 || psub>=1) stop("the proportion of subsamples should be between 0.5 and 1.")
   if (is.null(stra)) {
      wc<-which(y==1)
      wt<-which(y==0)
      nc<-floor(length(wc)*psub)
      nt<-floor(length(wt)*psub)
      N<-min(K,choose(length(wc),nc)*choose(length(wt),nt))
      }
   else {
      ns<-floor(max(stra)*psub)
      N<-min(K,choose(max(stra),ns))
      }
   u0<-as.integer(rnorm(1,-1000,1000))   
   for (i in 1:N) {
      set.seed(u0+i)
      if (is.null(stra)) {
        ss<-c(sample(wc,nc),sample(wt,nt))
        if (i==1) {
          nlam<-pclogit(x[sort(ss),],y[sort(ss)],stra=NULL,...)$lambda
          N.lam<-min(N.lam,length(nlam))
          lam<-seq(nlam[1],min(nlam),length.out=N.lam)
          out<-matrix(0,ncol(x),length(lam))
          }
        g<-pclogit(x[sort(ss),],y[sort(ss)],stra=NULL,lambda=lam,...)
      }
      else {
        ws<-sort(which(stra%in%sample(unique(stra),ns)))
        ng<-rep(seq(ns),as.integer(table(stra[ws])))
        nst<-as.integer(ng[rank(stra[ws])])
        if (i==1) {
          nlam<-pclogit(x[ws,],y[ws],stra=nst,...)$lambda
          N.lam<-min(N.lam,length(nlam))
          lam<-seq(nlam[1],min(nlam),length.out=N.lam)
          out<-matrix(0,ncol(x),length(lam))
          }
        g<-pclogit(x[ws,],y[ws],stra=nst,lambda=lam,...)
      }
      out0<-matrix(0,ncol(x),length(lam))
      beta<-as.matrix(g$beta)
      beta[abs(beta)>0]<-1
      out0[,1:ncol(beta)]<-beta  
      out<-out+out0
   }   
   if (is.null(colnames(x))) rownames(out) <- paste("V",1:ncol(x),sep="")
   else rownames(out) <- colnames(x)
   colnames(out) <- paste("s", 1:length(lam), sep="")
   maxs <- apply(out, 1, max)/N
   u <- order(maxs, decreasing=TRUE)
   mat <- cbind(u, maxs[u])
   colnames(mat)<-c("variable","sel.prob")
return(list(beta=beta,maxsel=mat,lambda=lam,K=N))
}
