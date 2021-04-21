      subroutine pelogit(parm,no,ni,x,y,nx,nlam,flmin,ulam,thr,maxit
     *,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,loc,idg,lq,st,q0)
      real x(no,ni),y(no,2),ulam(nlam),idg(ni),q0(no)
      real ca(nx,nlam),a0(nlam),dev(nlam),alm(nlam)
      integer ia(nx),nin(nlam),loc(lq,ni),st(no)
      real, dimension (:), allocatable :: xm,xs,ww
      integer, dimension (:), allocatable :: ju
      allocate(ww(1:no),stat=jerr)
      allocate(ju(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(xs(1:ni),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      call chkvars(no,ni,x,ju)
      if(maxval(ju) .gt. 0)goto 11871
      jerr=7777
      return
11871 continue
11880 do 11881 i=1,no
      ww(i)=sum(y(i,:))
      y(i,:)=y(i,:)/ww(i)
11881 continue
      sw=sum(ww)
      ww=ww/sw
      call standard(no,ni,x,ww,ju,xm,xs)
      if (sum(st).lt. 1.0) goto 11891
      call pclogit(parm,no,ni,x,y(:,1),ww,ju,nx,nlam,flmin,ulam,thr
     *,maxit,lmu,ca,ia,nin,dev0,dev,alm,nlp,jerr,loc,idg,lq,st,q0)
      goto 11893
11891 continue
      call pulogit(parm,no,ni,x,y(:,1),ww,ju,nx,nlam,flmin,ulam,thr
     *,maxit,lmu,a0,ca,ia,nin,dev0,dev,alm,nlp,jerr,loc,idg,lq)
11893 continue
      if(jerr.gt.0) return
      dev0=2.0*sw*dev0
11920 do 11921 k=1,lmu
      nk=nin(k)
11960 do 11961 l=1,nk
      ca(l,k)=ca(l,k)/xs(ia(l))
11961 continue
      a0(k)=a0(k)-dot_product(ca(1:nk,k),xm(ia(1:nk)))
11921 continue
      deallocate(ww,ju,xm,xs)
      return
      end

      subroutine pclogit(parm,no,ni,x,y,w,ju,nx,nlam,flmin,ulam,shri                        
     *,maxit,lmu,a,m,kin,dev0,dev,alm,nlp,jerr,loc,idg,lq,st,q0)
      parameter(sml=1.0e-5, pmin=1.0e-5,  big=9.9e30,mnlam=5
     *,devmax=0.999)
      real x(no,ni),y(no),w(no),ulam(nlam),q0(no)
      real a(nx,nlam),dev(nlam),alm(nlam),idg(ni)
      integer ju(ni),m(nx),kin(nlam),loc(lq,ni),st(no)
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga,fi
      integer, dimension (:), allocatable :: mm,ixx
      allocate(b(1:ni),stat=jerr)
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ga(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(bs(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ixx(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(r(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(v(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(fi(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(q(1:no),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      fmax=log(1.0/pmin-1.0)
      fmin=-fmax
      bta=parm
      omb=1.0-bta
      ixx=0
      al=0.0
      v=w*q0*(1-q0)
      r=w*(y-q0)
      q=q0
      dev0=deviance(no,w,y,q,pmin)
      xv=0.25
      if (flmin .ge. 1.0) goto 12151
      alf=flmin**(1.0/(nlam-1))
12151 continue      
      m=0
      mm=0
      nlp=0
      nin=nlp
      mnl=min(mnlam,nlam)
      bs=0.0
      b(1:ni)=0.0
      shr=shri*dev0
12160 do 12161 j=1,ni
      if(ju(j).eq.0)goto 12161
      ga(j)=abs(dot_product(r,x(:,j)))
12161 continue
12170 do 12171 ilm=1,nlam
      al0=al
      if(flmin .lt. 1.0)goto 12191
      al=ulam(ilm)
      goto 12181
12191 if(ilm .le. 2)goto 12201
      al=al*alf
      goto 12181
12201 if(ilm .ne. 1)goto 12211
      al=big
      goto 12221
12211 continue
      al0=0.0
12230 do 12231 j=1,ni
      if(ju(j).eq.0)goto 12231
      al0=max(al0,ga(j))
12231 continue
      al0=al0/max(bta,1.0e-3)
      al=alf*al0
12221 continue
12181 continue
      al2=al*omb
      al1=al*bta
      tlam=bta*(2.0*al-al0)
12240 do 12241 k=1,ni
      if(ixx(k).eq.1)goto 12241
      if(ju(k).eq.0)goto 12241
      if(ga(k).gt.tlam) ixx(k)=1
12241 continue
10680 continue
12251 continue
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))
12280 do 12281 j=1,ni
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)
12281 continue
12291 continue
      nlp=nlp+1
      dlx=0.0
12300 do 12301 k=1,ni
      if(ixx(k).eq.0)goto 12301
      bk=b(k)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k)*b(k)
      if(sum(abs(idg)) .eq. 0)goto 12327 
      u3=0.0
      do 12326 kk=1,lq
      if (loc(kk,k).eq.0) goto 12325                                                                                         
      u2=b(loc(kk,k))*idg(loc(kk,k))
      u3=u3+u2
12325 continue
12326 continue
      u=u+al2*u3*idg(k)
12327 continue                                                                  
      au=abs(u)-al1
      if(au .gt. 0.0)goto 12321
      b(k)=0.0
      goto 12331
12321 continue
      b(k)=sign(au,u)/(xv(k)+al2)
12331 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 12301
      dlx=max(dlx,xv(k)*d**2)
      r=r-d*v*x(:,k)
      if(mm(k) .ne. 0)goto 12351
      nin=nin+1
      if(nin.gt.nx)goto 12302
      mm(k)=nin
      m(nin)=k
12351 continue
12301 continue
12302 continue
      if(nin.gt.nx)goto 12292
      if(dlx.lt.shr)goto 12292
      if(nlp .le. maxit)goto 12391
      jerr=-ilm
      return
12391 continue
12401 continue
      nlp=nlp+1
      dlx=0.0
12410 do 12411 l=1,nin
      k=m(l)
      bk=b(k)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k)*b(k)
      if(sum(abs(idg)) .eq. 0)goto 12427 
      u3=0.0
      do 12426 kk=1,lq
      if (loc(kk,k).eq.0) goto 12425                                                                                         
      u2=b(loc(kk,k))*idg(loc(kk,k))
      u3=u3+u2
12425 continue
12426 continue
      u=u+al2*u3*idg(k)
12427 continue                                                                  
      au=abs(u)-al1
      if(au .gt. 0.0)goto 12431
      b(k)=0.0
      goto 12441
12431 continue
      b(k)=sign(au,u)/(xv(k)+al2)
12441 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 12411
      dlx=max(dlx,xv(k)*d**2)
      r=r-d*v*x(:,k)
12411 continue
      if(dlx.lt.shr)goto 12402
      if(nlp .le. maxit)goto 12481
      jerr=-ilm
      return
12481 continue
      goto 12401
12402 continue
      goto 12291
12292 continue
      if(nin.gt.nx)goto 12252
12490 do 12491 i=1,no
      fi(i)=0.0
      if(nin.gt.0) fi(i)=dot_product(b(m(1:nin)),x(i,m(1:nin)))
12491 continue
      call stratum(fi,st,no,q,pmin)
      v=w*q*(1.0-q)
      r=w*(y-q)
      ix=0
12560 do 12561 j=1,nin
      k=m(j)
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 12561
      ix=1
      goto 12562
12561 continue
12562 continue
      if(ix .ne. 0)goto 12581
12590 do 12591 k=1,ni
      if(ixx(k).eq.1)goto 12591
      if(ju(k).eq.0)goto 12591
      ga(k)=abs(dot_product(r,x(:,k)))
      if(ga(k) .le. al1)goto 12611
      ixx(k)=1
      ix=1
12611 continue
12591 continue
      if(ix.eq.1) go to 10680
      goto 12252
12581 continue
      goto 12251
12252 continue
      if(nin .le. nx)goto 12631
      jerr=-10000-ilm
      goto 12172
12631 continue
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))
      kin(ilm)=nin
      alm(ilm)=al
      lmu=ilm
      devi=deviance(no,w,y,q,pmin)
      dev(ilm)=(dev0-devi)/dev0
      if(ilm.lt.mnl)goto 12171
      if(flmin.ge.1.0)goto 12171
      if(dev(ilm).gt.devmax)goto 12172
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 12172
12171 continue
12172 continue
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx,fi)
      return
      end

      subroutine pulogit(parm,no,ni,x,y,w,ju,nx,nlam,flmin,ulam,shri                        
     *,maxit,lmu,a0,a,m,kin,dev0,dev,alm,nlp,jerr,loc,idg,lq)
      parameter(sml=1.0e-5, pmin=1.0e-5,  big=9.9e30,mnlam=5
     *,devmax=0.999)
      real x(no,ni),y(no),w(no),ulam(nlam)
      real a(nx,nlam),a0(nlam),dev(nlam),alm(nlam),idg(ni)
      integer ju(ni),m(nx),kin(nlam),loc(lq,ni)
      real, dimension (:), allocatable :: b,bs,v,r,xv,q,ga
      integer, dimension (:), allocatable :: mm,ixx
      allocate(b(0:ni),stat=jerr)
      allocate(xv(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ga(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(bs(0:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(mm(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(ixx(1:ni),stat=ierr)
      jerr=jerr+ierr
      allocate(r(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(v(1:no),stat=ierr)
      jerr=jerr+ierr
      allocate(q(1:no),stat=ierr)
      jerr=jerr+ierr
      if(jerr.ne.0) return
      fmax=log(1.0/pmin-1.0)
      fmin=-fmax
      vmin=(1.0+pmin)*pmin*(1.0-pmin)
      bta=parm
      omb=1.0-bta
      q0=dot_product(w,y)
      ixx=0
      al=0.0
      b(0)=log(q0/(1.0-q0))
      v=w*q0*(1-q0)
      r=w*(y-q0)
      q=q0
      xmz=sum(v)
      dev0=dev2(no,w,y,q,pmin)
      xv=0.25
      if (flmin .ge. 1.0) goto 22151
      alf=flmin**(1.0/(nlam-1))
22151 continue      
      m=0
      mm=0
      nlp=0
      nin=nlp
      mnl=min(mnlam,nlam)
      bs=0.0
      b(1:ni)=0.0
      shr=shri*dev0
22160 do 22161 j=1,ni
      if(ju(j).eq.0)goto 22161
      ga(j)=abs(dot_product(r,x(:,j)))
22161 continue
22170 do 22171 ilm=1,nlam
      al0=al
      if(flmin .lt. 1.0)goto 22191
      al=ulam(ilm)
      goto 22181
22191 if(ilm .le. 2)goto 22201
      al=al*alf
      goto 22181
22201 if(ilm .ne. 1)goto 22211
      al=big
      goto 22221
22211 continue
      al0=0.0
22230 do 22231 j=1,ni
      if(ju(j).eq.0)goto 22231
      al0=max(al0,ga(j))
22231 continue
      al0=al0/max(bta,1.0e-3)
      al=alf*al0
22221 continue
22181 continue
      al2=al*omb
      al1=al*bta
      tlam=bta*(2.0*al-al0)
22240 do 22241 k=1,ni
      if(ixx(k).eq.1)goto 22241
      if(ju(k).eq.0)goto 22241
      if(ga(k).gt.tlam) ixx(k)=1
22241 continue
20680 continue
22251 continue
      bs(0)=b(0)
      if(nin.gt.0) bs(m(1:nin))=b(m(1:nin))
22280 do 22281 j=1,ni
      if(ixx(j).gt.0) xv(j)=dot_product(v,x(:,j)**2)
22281 continue
22291 continue
      nlp=nlp+1
      dlx=0.0
22300 do 22301 k=1,ni
      if(ixx(k).eq.0)goto 22301
      bk=b(k)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k)*b(k)
      if(sum(abs(idg)) .eq. 0)goto 22327 
      u3=0.0
      do 22326 kk=1,lq
      if (loc(kk,k).eq.0) goto 22325                                                                                         
      u2=b(loc(kk,k))*idg(loc(kk,k))
      u3=u3+u2
22325 continue
22326 continue
      u=u+al2*u3*idg(k)
22327 continue                                                                  
      au=abs(u)-al1
      if(au .gt. 0.0)goto 22321
      b(k)=0.0
      goto 22331
22321 continue
      b(k)=sign(au,u)/(xv(k)+al2)
22331 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 22301
      dlx=max(dlx,xv(k)*d**2)
      r=r-d*v*x(:,k)
      if(mm(k) .ne. 0)goto 22351
      nin=nin+1
      if(nin.gt.nx)goto 22302
      mm(k)=nin
      m(nin)=k
22351 continue
22301 continue
22302 continue
      if(nin.gt.nx)goto 22292
      d=sum(r)/xmz
      if(d .eq. 0.0)goto 22371
      b(0)=b(0)+d
      dlx=max(dlx,xmz*d**2)
      r=r-d*v
22371 continue
      if(dlx.lt.shr)goto 22292
      if(nlp .le. maxit)goto 22391
      jerr=-ilm
      return
22391 continue
22401 continue
      nlp=nlp+1
      dlx=0.0
22410 do 22411 l=1,nin
      k=m(l)
      bk=b(k)
      gk=dot_product(r,x(:,k))
      u=gk+xv(k)*b(k)
      if(sum(abs(idg)) .eq. 0)goto 22427 
      u3=0.0
      do 22426 kk=1,lq
      if (loc(kk,k).eq.0) goto 22425                                                                                         
      u2=b(loc(kk,k))*idg(loc(kk,k))
      u3=u3+u2
22425 continue
22426 continue
      u=u+al2*u3*idg(k)
22427 continue                                                                  
      au=abs(u)-al1
      if(au .gt. 0.0)goto 22431
      b(k)=0.0
      goto 22441
22431 continue
      b(k)=sign(au,u)/(xv(k)+al2)
22441 continue
      d=b(k)-bk
      if(abs(d).le.0.0)goto 22411
      dlx=max(dlx,xv(k)*d**2)
      r=r-d*v*x(:,k)
22411 continue
      d=sum(r)/xmz
      if(d .eq. 0.0)goto 22461
      b(0)=b(0)+d
      dlx=max(dlx,xmz*d**2)
      r=r-d*v
22461 continue
      if(dlx.lt.shr)goto 22402
      if(nlp .le. maxit)goto 22481
      jerr=-ilm
      return
22481 continue
      goto 22401
22402 continue
      goto 22291
22292 continue
      if(nin.gt.nx)goto 22252
22490 do 22491 i=1,no
      fi=b(0)                                                        
      if(nin.gt.0) fi=fi+dot_product(b(m(1:nin)),x(i,m(1:nin)))           
      if(fi .ge. fmin)goto 22511                                          
      q(i)=0.0                                                            
      goto 22501                                                          
22511 if(fi .le. fmax)goto 22521                                          
      q(i)=1.0                                                            
      goto 22531                                                          
22521 continue                                                            
      q(i)=1.0/(1.0+exp(-fi))                                             
22531 continue                                                             
22501 continue
22491 continue 
      v=w*q*(1.0-q)
      xmz=sum(v)
      if(xmz.le.vmin)goto 22252
      r=w*(y-q)
      if(xmz*(b(0)-bs(0))**2 .ge. shr)goto 22551      
      ix=0
22560 do 22561 j=1,nin
      k=m(j)
      if(xv(k)*(b(k)-bs(k))**2.lt.shr)goto 22561
      ix=1
      goto 22562
22561 continue
22562 continue
      if(ix .ne. 0)goto 22581
22590 do 22591 k=1,ni
      if(ixx(k).eq.1)goto 22591
      if(ju(k).eq.0)goto 22591
      ga(k)=abs(dot_product(r,x(:,k)))
      if(ga(k) .le. al1)goto 22611
      ixx(k)=1
      ix=1
22611 continue
22591 continue
      if(ix.eq.1) go to 20680
      goto 22252
22581 continue
22551 continue
      goto 22251
22252 continue
      if(nin .le. nx)goto 22631
      jerr=-10000-ilm
      goto 22172
22631 continue
      if(nin.gt.0) a(1:nin,ilm)=b(m(1:nin))
      kin(ilm)=nin
      a0(ilm)=b(0)
      alm(ilm)=al
      lmu=ilm
      devi=dev2(no,w,y,q,pmin)
      dev(ilm)=(dev0-devi)/dev0
      if(xmz.le.vmin)goto 22172
      if(ilm.lt.mnl)goto 22171
      if(flmin.ge.1.0)goto 22171
      if(dev(ilm).gt.devmax)goto 22172
      if(dev(ilm)-dev(ilm-1).lt.sml)goto 22172
22171 continue
22172 continue
      deallocate(b,bs,v,r,xv,q,mm,ga,ixx)
      return
      end

      subroutine standard(no,ni,x,w,ju,xm,xs)
      real x(no,ni),w(no),xm(ni),xs(ni)
      integer ju(ni)
11970 do 11971 j=1,ni
      if(ju(j).eq.0)goto 11971
      xm(j)=dot_product(w,x(:,j))
      x(:,j)=x(:,j)-xm(j)
      xs(j)=sqrt(dot_product(w,x(:,j)**2))
      x(:,j)=x(:,j)/xs(j)
11971 continue
      return
      end

      subroutine chkvars(no,ni,x,ju)
      real x(no,ni)
      integer ju(ni)
10860 do 10861 j=1,ni
      ju(j)=0
      t=x(1,j)
10870 do 10871 i=2,no
      if(x(i,j).eq.t)goto 10871
      ju(j)=1
      goto 10872
10871 continue
10872 continue
10861 continue
      return
      end

      function dev2(n,w,y,p,pmin)                                         
      real w(n),y(n),p(n)                                                 
      pmax=1.0-pmin                                                       
      s=0.0                                                               
14650 do 14651 i=1,n                                                      
      pi=min(max(pmin,p(i)),pmax)                                         
      s=s-w(i)*(y(i)*log(pi)+(1.0-y(i))*log(1.0-pi))                      
14651 continue                                                            
      dev2=s                                                              
      return                                                              
      end                                                                 

      function deviance(n,w,y,p,pmin)
      real w(n),y(n),p(n)
      pmax=1.0-pmin
      s=0.0
12650 do 12651 i=1,n
      pi=min(max(pmin,p(i)),pmax)
      s=s-w(i)*y(i)*log(pi)
12651 continue
      deviance=s
      return
      end

      subroutine stratum(x,st,n,out,pmin)
      real x(n),out(n)
      integer st(n)
      fmax=log(1.0/pmin-1.0)
      fmin=-fmax
      do 18010 i=1,n
      tti=0.0
      do 18070 j=1,n
      if (st(j).ne.st(i)) goto 18053
      xx=x(j)-x(i)
      if (xx .ge. fmin) goto 18050
      goto 18053
18050 if (xx .le. fmax) goto 18052      
      goto 18055
18052 continue      
      tti=tti+exp(xx)
18053 continue
18070 continue
      out(i)=1/tti
      goto 18059
18055 out(i)=0.0
18059 continue
18010 continue
      return
      end
