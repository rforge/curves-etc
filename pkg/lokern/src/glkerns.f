      subroutine glkerns(t,x,n,tt,m,nue,kord,ihom,irnd,
     .     ismo,m1,tl,tu,s,sig,wn,w1,b,y)
c----------------------------------------------------------------------*
c-----------------------------------------------------------------------
c       short-version: oct 1996
c
c       purpose:
c
c       general subroutine for kernel smoothing:
c       computation of iterative plug-in algorithm for global bandwidth
c       selection for kernels with (nue,kord) = (0,2),(0,4),(1,3) or 
c       (2,4).
c-----------------------------------------------------------------------
c  used subroutines: coff, resest, kernel with further subroutines
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      dimension x(n),t(n),tt(m),y(m),wn(0:n,5),s(0:n)
      dimension w1(m1,3)
      dimension bias(2,0:2),vark(2,0:2),fak2(2:4)
c-
c-------- 1. initialisations 
      data bias/.2,.04762,.4286,.1515,1.33,.6293/
      data vark/.6,1.250,2.143,11.93,35.0,381.6/
      data fak2/4.,36.,576./

      nyg=0
      inputs=0
c
      if(nue.gt.4.or.nue.lt.0) stop
      if(nue.gt.2.and.ismo.eq.0) stop
      if(n.le.2) stop
      if(m.lt.1) stop
      if(m1.lt.3) stop
      kk=(kord-nue)/2 
      if(2*kk+nue.ne.kord) kord=nue+2
      if(kord.gt.4.and.ismo.eq.0) kord=nue+2
      if(kord.gt.6.or.kord.le.nue) kord=nue+2
      if(ismo.ne.0.and.b.le.0) ismo=0
      rvar=sig
c-
c-------- 2. computation of s-sequence 
      s0=1.5*t(1)-0.5*t(2)
      sn=1.5*t(n)-0.5*t(n-1)
      if(s(n).le.s(0)) then
         inputs=1
         do 20 i=1,n-1
            s(i)=.5*(t(i)+t(i+1))
 20      continue
         s(0)=s0
         s(n)=sn
         if(ismo.ne.0.and.irnd.ne.0) goto 160
      else
         if(ismo.ne.0) goto 160
      end if
c-
c-------- 3. computation of minimal, maximal allowed bandwidth
      bmax=(sn-s0)*.5
      bmin=(sn-s0)/dble(n)*dble(kord-1)*.6
c-
c-------- 4. compute tl,tu 
      itt=0
40    if (tu.le.tl) then
        tl=.933*s0+.067*sn
        tu=.067*s0+.933*sn
        itt=itt+1
      end if
      tl=max(s0,tl)
      tu=min(sn,tu)
c-
c-------- 5. compute indices 
      il=1
      iu=n
      wn(1,1)=0.0
      wn(n,1)=0.0
      do 50 i=1,n
        if(t(i).le.tl.or.t(i).ge.tu) wn(i,1)=0.0
        if(t(i).gt.tl.and.t(i).lt.tu) wn(i,1)=1.0
        if(t(i).lt.tl) il=i+1
        if(t(i).le.tu) iu=i
 50   continue
      nn=iu-il+1
      if(nn.eq.0.and.itt.eq.0) then
        tu=tl-1.0
        goto 40
      end if
      if(nn.eq.0.and.itt.eq.1) then
        tu=sn
        tl=s0
        goto 40
      end if
c-
c-------- 6. compute t-grid for integral approximation 
      do 60 i=1,m1
         w1(i,2)=1.0
         w1(i,1)=tl+(tu-tl)*dble(i-1)/dble(m1-1)
 60   continue
c-
c-------- 7. calculation of weight function
      alpha=1.d0/dble(13)
      do 70 i=il,iu
         xi=(t(i) - tl)/alpha/(tu-tl)
         if(xi.gt.1) goto 71
         wn(i,1)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
 70   continue
 71   do 72 i=iu,il,-1
        xi=(tu-t(i))/alpha/(tu-tl)
        if(xi.gt.1) goto 73
        wn(i,1)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
 72   continue
 73   do 74 i=1,m1
         xi=(w1(i,1)-tl)/alpha/(tu-tl)
         if(xi.gt.1) goto 75
         w1(i,2)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
 74   continue 
 75   do 76 i=m1,1,-1
         xi=(tu-w1(i,1))/alpha/(tu-tl)
         if(xi.gt.1) goto 77
         w1(i,2)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
 76   continue 
 77   continue 
c-
c-------- 8. compute constants for iteration
      ex=1./dble(kord+kord+1)
      kk2=(kord-nue)
      kk=kk2/2
c-
c-------- 9. estimating variance and smoothed pseudoresiduals
      rvar=sig
      if((sig.le..0).and.(ihom.eq.0)) 
     .     call resest(t(il),x(il),nn,wn(il,2),r2,sig)
      if(ihom.ne.0) then 
         call resest(t,x,n,wn(1,2),snr,sig)
         bres=max(bmin,.2*nn**(-.2)*(s(iu)-s(il-1)))
         do 91 i=1,n
            wn(i,3)=t(i)
            wn(i,2)=wn(i,2)*wn(i,2)
 91      continue 
         call kernel(t,wn(1,2),n,bres,0,kk2,nyg,s,
     .        wn(il,3),nn,wn(il,4))
      else
         call coff(wn(1,4),n,sig)
      end if
c-
c-------- 10. estimate/compute integral constant
100   vi=0.
      do 101 i=il,iu
         vi=vi+wn(i,1)*n*(s(i)-s(i-1))**2*wn(i,4)       
 101  continue
c- 
c-------- 11. refinement of s-sequence for random design
      if(inputs.eq.1.and.irnd.eq.0) then
        do 110 i=0,n
          wn(i,5)=dble(i)/dble(n+1)
          wn(i,2)=(dble(i)+.5)/dble(n+1)
110       wn(i,3)=wn(i,2)
        exs=-dble(3*kord+1)/dble(6*kord+3)
        exsvi=dble(kord)/dble(6*kord+3)
        bs=0.1*(vi/(sn-s0)**2)**exsvi*n**exs
        call kernel(wn(1,5),t,n,bs,0,2,nyg,wn(0,3),wn(0,2),n+1,s(0))
        isort=0
        vi=0.0
111     do 112 i=1,n
        vi=vi+wn(i,1)*n*(s(i)-s(i-1))**2*wn(i,4)       
          if(s(i).lt.s(i-1)) then
            ssi=s(i-1)
            s(i-1)=s(i)
            s(i)=ssi
            isort=1
          end if
112     continue
        if(isort.eq.1) goto 111
        if(ismo.ne.0) goto 160
      end if
      b=bmin*2.
c-
c-------- 12. compute inflation constant and exponent and loop of iterations
      const=dble(2*nue+1)*fak2(kord)*vark(kk,nue)*vi
     .       /(dble(2*kord-2*nue)*bias(kk,nue)**2*dble(n))
      fac=1.1*(1.+(nue/10.)+0.05*(kord-nue-2.))
     .       *n**(2./dble((2*kord+1)*(2*kord+3)))
      itende=1+2*kord+kord*(2*kord+1)

      do 120 it=1,itende
c-
c-------- 13. estimate derivative of order kord in iterations
        b2=b*fac
        b2=max(b2,bmin/dble(kord-1)*dble(kord+1))
        b2=min(b2,bmax)
        call kernel(t,x,n,b2,kord,kord+2,nyg,s,w1(1,1),m1,w1(1,3))
c-
c-------- 14. estimate integralfunctional in iterations
        xmy2=.75*(w1(1,2)*w1(1,3)*w1(1,3)+w1(m1,2)*w1(m1,3)*w1(m1,3))
        do 140 i=2,m1-1
140       xmy2=xmy2+w1(i,2)*w1(i,3)*w1(i,3)
        xmy2=xmy2*(tu-tl)/dble(m1)
c-
c-------- 15. finish of iterations
        b=(const/xmy2)**ex
        b=max(bmin,b)
        b=min(bmax,b)
120     continue
c-
c-------- 16  compute smoothed function with plug-in bandwidth
160   call kernel(t,x,n,b,nue,kord,nyg,s,tt,m,y)
c-
c-------- 17. variance check 
      irnd=1-irnd
      if(ihom.ne.0) sig=rvar
      if(rvar.eq.sig.or.r2.lt..88.or.ihom.ne.0.or.nue.gt.0) return
      ii=0
      iil=0
      j=2
      tll=max(tl,tt(1))
      tuu=min(tu,tt(m))
      do 170 i=il,iu
         if(t(i).lt.tll.or.t(i).gt.tuu) goto 170
         ii=ii+1
         if(iil.eq.0) iil=i
171       if(tt(j).lt.t(i)) then
            j=j+1
            if(j.le.m) goto 171
            end if
       wn(ii,3)=x(i)-y(j)+(y(j)-y(j-1))*(tt(j)-t(i))/(tt(j)-tt(j-1))
170    continue
      if(iil.eq.0.or.ii-iil.lt.10) then
        call resest(t(il),wn(1,3),nn,wn(1,4),snr,rvar)
       else
        call resest(t(iil),wn(1,3),ii,wn(1,4),snr,rvar)
      end if
      q=sig/rvar
      if(q.le.2.) return
      if(q.gt.5..and.r2.gt..95) rvar=rvar*.5
      sig=rvar
      call coff(wn(1,4),n,sig)
      goto 100
      end

