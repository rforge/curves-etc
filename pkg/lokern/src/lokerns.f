      subroutine lokern_s(t,x, tt,y, n,m,nue,kord, hetero,isrand,
     .     inputb,m1,tl,tu,s,sig,wn,w1, wm,ban)
c----------------------------------------------------------------------*
c-----------------------------------------------------------------------
c       Short-version: January 1997
c
c       Purpose:
c
c       General subroutine for kernel smoothing:
c       Computation of iterative plug-in algorithm for local bandwidth
c       selection for kernels with
c       (nue,kord) = (0,2), (0,4), (1,3) or (2,4).
c-----------------------------------------------------------------------
c  used subroutines: constV, resest, kernel with further subroutines
c-----------------------------------------------------------------------
c Args
      integer n, m, nue,kord
      double precision t(n),x(n), tt(m), tl,tu, s(0:n), sig
      logical hetero, isrand, inputb
c	 inputb (was "smo", now same as in R):
c	 if TRUE, do not compute bandwidths but use ban(.)
      integer m1
      double precision wn(0:n,5), w1(m1,3), wm(m),ban(m), y(m)
c Var
      logical inputs, needsrt
      integer nyg, i,ii,iil,itt,il,iu,itende,it, j,
     1     kk,kk2, kordv, nn, nuev,nyl
      double precision bias(2,0:2),vark(2,0:2),fak2(2:4),
     1     rvar, s0,sn, b,b2,bmin,bmax,bres,bvar,bs,alpha,ex,exs,exsvi,
     2     r2,snr,vi,ssi,const,fac, g1,g2,dist, q,tll,tuu, wstep,
     3     xi,xh,xxh,xmy2
c-
c-------- 1. initialisations
      data bias/.2, .04762, .4286, .1515, 1.33, .6293/
      data vark/.6,  1.250, 2.143, 11.93, 35.0, 381.6/
      data fak2/4.,36.,576./
      nyg=0
      inputs = .false.
      nyl=1
c Stop for invalid inputs (impossible when called from R's lokerns())

c     0 <= nue <= 4;  nue <= 2 if(! smo)
      if(nue.gt.4 .or. nue.lt.0) call rexit("nue must be in 0..4")
      if(nue.gt.2 .and. .not. inputb)
     +     call rexit("nue must be in 0..2 if not 'inputb'")
      if(n .le. 2) call rexit("n <= 2")
      if(m .lt. 1) call rexit("m < 1")
      if(m1.lt. 3) call rexit("m1 < 3")

c     kord - nue must be even :
      kk=(kord-nue)/2
      if(2*kk + nue .ne. kord)        kord=nue+2
      if(kord.gt.4 .and. .not.inputb) kord=nue+2
      if(kord.gt.6 .or. kord.le.nue)  kord=nue+2
      rvar=sig

c-------- 2. computation of s-sequence
      s0=1.5*t(1)-0.5*t(2)
      sn=1.5*t(n)-0.5*t(n-1)
      if(s(n).le.s(0)) then
         inputs= .true.
         do i=1,n-1
            s(i)=.5*(t(i)+t(i+1))
         end do
         s(0)=s0
         s(n)=sn
         if(inputb .and. .not.isrand) goto 230
      else
         if(inputb) goto 230
      end if
c-
c-------- 3. computation of minimal, maximal allowed global bandwidth
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
      do i=1,n
        if(t(i).le.tl .or.  t(i).ge.tu) wn(i,1)=0.0
        if(t(i).gt.tl .and. t(i).lt.tu) wn(i,1)=1.0
        if(t(i).lt.tl) il=i+1
        if(t(i).le.tu) iu=i
      end do
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
      do i=1,m1
         w1(i,2)=1.0
         w1(i,1)=tl+(tu-tl)*dble(i-1)/dble(m1-1)
      end do
c-
c-------- 7. calculation of weight function
      alpha=1.d0/dble(13)
      do i=il,iu
         xi=(t(i) - tl)/alpha/(tu-tl)
         if(xi.gt.1) goto 71
         wn(i,1)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
      end do
 71   do i=iu,il,-1
        xi=(tu-t(i))/alpha/(tu-tl)
        if(xi.gt.1) goto 73
        wn(i,1)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
      end do
 73   do i=1,m1
         xi=(w1(i,1)-tl)/alpha/(tu-tl)
         if(xi.gt.1) goto 75
         w1(i,2)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
      end do
 75   do i=m1,1,-1
         xi=(tu-w1(i,1))/alpha/(tu-tl)
         if(xi.gt.1) goto 77
         w1(i,2)=(10.0-15*xi+6*xi*xi)*xi*xi*xi
      end do
 77   continue
c-
c-------- 8. compute constants for iteration
      ex=1./dble(kord+kord+1)
      kk2=(kord-nue)
      kk=kk2/2
c-
c-------- 9. estimating variance and smoothed pseudoresiduals
      if(sig .le. 0. .and. .not.hetero) then
        call resest(t(il),x(il),nn,wn(il,2),r2,sig)
      endif
      if(hetero) then
        call resest(t,x,n,wn(1,2),snr,sig)
        bres=max(bmin,.2*nn**(-.2)*(s(iu)-s(il-1)))
        do i=1,n
           wn(i,3)=t(i)
           wn(i,2)=wn(i,2)*wn(i,2)
        end do
        call kernel(t,wn(1,2),n,bres,0,kk2,nyg,s, wn(1,3),n,wn(1,4))
      else
c       not hetero
        call constV(wn(1,4),n,sig)
      end if
c-
c-------- 10. [LOOP:] estimate/compute integral constant
100   vi=0.
      do i=il,iu
         vi=vi+ wn(i,1)*n*(s(i)-s(i-1))**2 * wn(i,4)
      end do
c-
c-------- 11. refinement of s-sequence for random design
      if(inputs .and. isrand) then
        do i=0,n
          wn(i,5)=dble(i)/dble(n+1)
          wn(i,2)=(dble(i)+.5)/dble(n+1)
          wn(i,3)=wn(i,2)
        end do
        exs= -dble(3*kord+1) / dble(6*kord+3)
        exsvi=dble(kord)     / dble(6*kord+3)
        bs=0.1*(vi/(sn-s0)**2)**exsvi * n**exs
        call kernel(wn(1,5),t,n,bs,0,2,nyg,wn(0,3),wn(0,2),n+1,s(0))

        vi=0.0
111     needsrt=.false.
        do i=1,n
           vi=vi+ wn(i,1)*n*(s(i)-s(i-1))**2 * wn(i,4)
           if(s(i).lt.s(i-1)) then
              ssi=s(i-1)
              s(i-1)=s(i)
              s(i)=ssi
              needsrt=.true.
           end if
        end do
        if(needsrt) goto 111
        if(inputb) goto 230
      end if
      b=bmin*2.
c-
c-------- 12. compute inflation constant and exponent and loop of iterations
      const=dble(2*nue+1)*fak2(kord)*vark(kk,nue)*vi
     .       /(dble(2*kord-2*nue)*bias(kk,nue)**2*dble(n))
      fac=1.1*(1.+(nue/10.)+0.05*(kord-nue-2.))
     .       *n**(2./dble((2*kord+1)*(2*kord+3)))

c     itende=1+2*kord+kord*(2*kord+1)
      itende = (1 + 2*kord) * (1 + kord)
c     ^^^^^^  *fixed* number of iterations ( <== theory !)

      do it=1,itende
c-
c-------- 13. estimate derivative of order kord in iterations
        b2=b*fac
        b2=max(b2,bmin/dble(kord-1)*dble(kord+1))
        b2=min(b2,bmax)
        call kernel(t,x,n,b2,kord,kord+2,nyg,s,w1(1,1),m1,w1(1,3))
c-
c-------- 14. estimate integralfunctional in iterations
        xmy2= .75*(w1(1,2)*w1(1,3)*w1(1,3) + w1(m1,2)*w1(m1,3)*w1(m1,3))
        do i=2,m1-1
           xmy2=xmy2+w1(i,2)*w1(i,3)*w1(i,3)
        end do
        xmy2=xmy2*(tu-tl)/dble(m1)
c-
c-------- 15. finish of iterations
        b=(const/xmy2)**ex
        b=max(bmin,b)
        b=min(bmax,b)
      end do

c-------- 16  compute smoothed function with global plug-in bandwidth
      call kernel(t,x,n,b,nue,kord,nyg,s,tt,m,y)
c-
c-------- 17. variance check
      if(hetero) sig=rvar
      if(hetero.or. r2.lt.0.88 .or. nue.gt.0) goto 180
      ii=0
      iil=0
      j=2
      tll=max(tl,tt(1))
      tuu=min(tu,tt(m))
      do 170 i=il,iu
         if(t(i).lt.tll.or.t(i).gt.tuu) goto 170
         ii=ii+1
         if(iil.eq.0) iil=i
 171     if(tt(j).lt.t(i)) then
            j=j+1
            if(j.le.m) goto 171
         end if
         wn(ii,3)=x(i)-y(j)+(y(j)-y(j-1))*(tt(j)-t(i))/(tt(j)-tt(j-1))
 170  continue
      if(iil.eq.0.or.ii-iil.lt.10) then
         call resest(t(il),wn(1,3),nn,wn(1,4),snr,rvar)
      else
         call resest(t(iil),wn(1,3),ii,wn(1,4),snr,rvar)
      end if
      q=sig/rvar
      if(q.le.2.) then
         call constV(wn(1,4),n,sig)
         goto 180
      end if
      if(q.gt.5. .and. r2.gt..95) rvar=rvar*.5
      sig=rvar
      call constV(wn(1,4),n,sig)
      goto 100
c     -------- end loop
c
c-------- 18. local initializations
 180  bvar=b
      nuev=0
      kordv=2
c-
c-------- 19. compute inner bandwidths
      g1=0.86*(1.+dble(kord-nue-2)*.05)*b
      g1=g1*dble(n)**(4./dble(2*kord+1)/(2*kord+5))
      g1=max(g1,bmin/dble(kord-1)*dble(kord+1))
      g1=min(g1,bmax)

      g2=1.4*(1.+dble(kord-nue-2)*0.05)*b
      g2=g2*dble(n)**(2./dble(2*kord+1)/dble(2*kord+3))
      g2=max(g2,bmin)
      g2=min(g2,bmax)
c-
c-------- 20. estimate/compute integral constant vi locally
      do i=1,n
         wn(i,4)=dble(n)*wn(i,4)*(s(i)-s(i-1))
      end do
      do j=1,m
         ban(j)=bvar
         wm(j)=tt(j)
         if(tt(j).lt.s(0)+g1) then
            dist=((tt(j)-g1-s(0))/g1)**2
            ban(j)=bvar*(1.0+1.0*dist)
            ban(j)=min(ban(j),bmax)
            wm(j)=tt(j)+.5*dist*g1
         else if(tt(j).gt.s(n)-g1) then
            dist=((tt(j)-s(n)+g1)/g1)**2
            ban(j)=bvar*(1.0+1.0*dist)
            ban(j)=min(ban(j),bmax)
            wm(j)=tt(j)-.5*dist*g1
         end if
      end do
      call kernel(t,wn(1,4),n,bvar,nuev,kordv,nyl,s,wm,m,ban)
c-
c-------- 21. estimation of kord-th derivative locally
      wstep=(tt(m)-tt(1))/dble(m1-2)
      do j=2,m1
         w1(j,2)=tt(1)+dble(j-2)*wstep
         w1(j,1)=tt(1)+dble(j-1.5)*wstep
      end do
      w1(1,1)=tt(1)+.5*wstep

      call kernel(t,x,n,g1,kord,kord+2,nyg,s,w1(2,2),m1-1,w1(2,3))

      do j=2,m1
        w1(j,3)=w1(j,3)*w1(j,3)
      end do
      do j=1,m
         y(j)=g2
         if(tt(j).lt.s(0)+g1) then
            y(j)=g2*(1.0+1.0*((tt(j)-g1-s(0))/g1)**2)
            y(j)=min(y(j),bmax)
         else if(tt(j).gt.s(n)-g1) then
            y(j)=g2*(1.0+1.0*((tt(j)-s(n)+g1)/g1)**2)
            y(j)=min(y(j),bmax)
         end if
      end do

      call kernp(w1(2,2),w1(2,3),m1-1,g2,nuev,kordv,nyl, w1(1,1),wm,m,y)
c-
c-------- 22. finish
c-what for? irnd=1-irnd
      do j=1,m
         xh=bmin**(2*kord+1)*abs(y(j))*vi/const
         xxh=const*abs(ban(j))/vi/bmax**(2*kord+1)
         if(ban(j).lt.xh) then
            ban(j)=bmin
         else if(y(j).lt.xxh) then
            ban(j)=bmax
         else
            ban(j)=(const*ban(j)/y(j)/vi)**ex
         end if
      end do
c-
c-------- 23. compute smoothed function with local plug-in bandwidth
 230  do j=1,m
        y(j)=ban(j)
      end do
      call kernel(t,x,n,b,nue,kord,nyl,s,tt,m,y)
c-  return #{iter} (iff aplicable)
      m1=itende
      return
      end
c     --- of lokern_s()
