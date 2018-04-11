group.matrix<-function(group,type,sgnc) {
   
   if (is.matrix(group)) {
      adjm<-group
      n<-nrow(adjm)
      p<-ncol(adjm)
      if (n!=p) stop("adjacency matrix should be a square matrix.")
      diag(adjm)<-0
      uniq<-sort(unique(as.numeric(adjm)))
      if (!identical(as.integer(uniq),as.integer(c(0,1)))) stop("adjacency matrix should have only 0 and 1.")
      tadjm<-t(adjm)
      low<-t(t(adjm[lower.tri(adjm)]))
      upp<-t(t(tadjm[lower.tri(tadjm)]))
      if (!identical(as.integer(low),as.integer(upp))) stop("adjacency matrix should be symmetric.")
      dg<-apply(adjm,1,sum)
      loc<-matrix(0,max(2,max(dg)),p)
      ws<-which(adjm==1,arr.ind=TRUE)
      k<-unique(ws[,2])
      for (j in 1:length(k)) {
        ki<-ws[ws[,2]==k[j],1]
        loc[1:length(ki),k[j]]<-ki
      }
   }  
   else {
      if (type=="ring") {
        p<-sum(group)
        loc<-matrix(0,p,2)
        cc<-cumsum(group)+1
        cc<-c(1,cc)
        for (i in 1:length(group)) {
          gi<-group[i]
          sq<-cc[i]:(cc[i+1]-1)
          if (gi==2) {
              loc[sq[1],1]<-sq[2]
              loc[sq[2],1]<-sq[1]
              }
          else if (gi>2) {
              loc[sq[1],1]<-sq[gi]
              loc[sq[1],2]<-sq[2]
              loc[sq[gi],1]<-sq[gi-1]
              loc[sq[gi],2]<-sq[1]
              for (j in 2:(gi-1)) {
                loc[sq[j],1]<-sq[j-1]
                loc[sq[j],2]<-sq[j+1]
              }
           }
        }
      loc<-t(loc)
      }
      if (type=="fcon") {
        p<-sum(group)
        loc<-matrix(0,p,(max(group)))
        cc<-cumsum(group)+1
        cc<-c(1,cc)
        for (i in 1:length(group)) {
          gi<-group[i]
          sq<-cc[i]:(cc[i+1]-1)
          if (gi>=2) {
             for (j in 1:gi) loc[sq[j],1:(gi-1)]<-sq[-j]
             }
        }
      loc<-t(loc[,-ncol(loc)])    
      }
   }
   idg<-rep(0,ncol(loc))
   dg<-apply(loc>0,2,sum)
   w<-which(dg!=0)
   if (length(w)>=1) idg[w]<-1/sqrt(dg[w])
   if (!is.null(sgnc)) idg<-idg*sgnc
     
return(list(loc=loc,idg=idg))
}
