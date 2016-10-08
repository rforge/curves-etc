      subroutine glkern_s(t,x, tt,y, n,m,nue,kord, hetero,isrand,
     .     inputb,m1,tl,tu,s,sig,wn,w1,b, trace)
c----------------------------------------------------------------------*
c-----------------------------------------------------------------------
c       Short-version: Oct 1996
c
c       Purpose:
c
c       General subroutine for kernel smoothing:
c       Computation of iterative plug-in algorithm for global bandwidth
c       selection for kernels with
c       (nue,kord) = (0,2), (0,4), (1,3) or (2,4).
c-----------------------------------------------------------------------
c  used subroutines: constV, resest, kernel with further subroutines
c-----------------------------------------------------------------------
c Args
      integer n, m, nue,kord
      double precision t(n),x(n), tt(m), tl,tu, s(0:n), sig
      logical hetero, isrand, inputb
c			      inputb (was "smo", now same as in R):
c	 if TRUE, do not compute bandwidth, but *use* 'b' given as input
      integer m1, trace
      double precision wn(0:n,5),w1(m1,3), b, y(m)
c Var
      logical inputs, needsrt
      integer nyg, i,ii,iil,itt,il,iu,itende,it, j, kk,kk2, nn
      double precision bias(2,0:2),vark(2,0:2),fak2(2:4),
     1     rvar, s0,sn, b2,bmin,bmax,bres,bs,alpha,ex,exs,exsvi,
     2     r2,snr,vi,ssi,const,fac, q,tll,tuu, xi,xmy2
c-
c-------- 1. initialisations
      data bias/.2, .04762, .4286, .1515, 1.33, .6293/
      data vark/.6,  1.250, 2.143, 11.93, 35.0, 381.6/
      data fak2/4.,36.,576./
      nyg=0
      inputs = .false.

c Stop for invalid inputs (impossible when called from R's glkerns())

c     0 <= nue <= 4;  nue <= 2 if(! inputb)
      if(nue.gt.4 .or. nue.lt.0) call rexit("nue must be in 0..4")
      if(nue.gt.2 .and. .not. inputb)
     +     call rexit("nue must be in 0..2 if not 'inputb'")
      if(n .le. 2) call rexit("n <= 2")
      if(m .lt. 1) call rexit("m < 1")
      if(m1.lt. 3) call rexit("m1 < 3")

c     kord - nue must be even  <<--- MM: *UN*desirable (want to fix kord, vary 'nue') !
c                        ----  <<---     but a short glimpse at coff*() ./auxkerns.f
c                                        reveals how much work would be needed to change this
      kk=(kord-nue)/2
      if(2*kk + nue .ne. kord)       kord=nue+2
      if(kord.gt.4 .and. .not.inputb)   kord=nue+2
      if(kord.gt.6 .or. kord.le.nue) kord=nue+2
      if(inputb.and.b.le.0) inputb=.false.
      rvar=sig
c- -Wall (erronously warning if not)
      bmin=1
      bmax=1
      ex=1
      itende= -1

      il=1
      iu=n

      if(trace .gt. 0) call monit0(0, n, m, nue, kord,
     +     inputb, isrand, b, trace)

c-------- 2. computation of s-sequence
      if(trace .gt. 0) call monit1(2, trace)
      s0=1.5*t(1)-0.5*t(2)
      sn=1.5*t(n)-0.5*t(n-1)
      if(s(n).le.s(0)) then
         inputs= .true.
         do i=1,n-1
            s(i)=.5*(t(i)+t(i+1))
         end do
         s(0)=s0
         s(n)=sn
         if(inputb.and. .not.isrand) goto 160
      else
         if(inputb) goto 160
      end if
c-
c-------- 3. computation of minimal, maximal allowed global bandwidth
      if(trace .gt. 0) call monit1(3, trace)
      bmax=(sn-s0)*.5
      bmin=(sn-s0)/dble(n)*dble(kord-1)*.6
c-
c-------- 4. compute tl,tu
      if(trace .gt. 0) call monit1(4, trace)
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
      if(trace .gt. 0) call monit1(5, trace)
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
      if(trace .gt. 0) call monit1(6, trace)
      do i=1,m1
         w1(i,2)=1.0
         w1(i,1)=tl+(tu-tl)*dble(i-1)/dble(m1-1)
      end do
c-
c-------- 7. calculation of weight function
      if(trace .gt. 0) call monit1(7, trace)
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
      if(trace .gt. 0) call monit1(8, trace)
      ex=1./dble(kord+kord+1)
      kk2=(kord-nue)
      kk=kk2/2
c-
c-------- 9. estimating variance and smoothed pseudoresiduals
      if(trace .gt. 0) call monit1(9, trace)
      rvar=sig
      if(hetero) then
         call resest(t,x,n,wn(1,2),snr,sig)
         bres=max(bmin, .2* dble(nn)**(-.2) * (s(iu)-s(il-1)))
         do i=1,n
            wn(i,3)=t(i)
            wn(i,2)=wn(i,2)*wn(i,2)
         end do
         call kernel(t,wn(1,2),n,bres,0,kk2,nyg,s,
     .        wn(il,3),nn,wn(il,4), trace)
      else !-- not hetero
         if(sig .le. 0.) then
            call resest(t(il),x(il),nn,wn(il,2),r2,sig)
         end if
         call constV(wn(1,4),n,sig)
      end if
c-
c-------- 10. [LOOP:] estimate/compute integral constant
100   vi=0.
      if(trace .gt. 0) call monit1(10, trace)
      do i=il,iu
         vi=vi+ wn(i,1)*n*(s(i)-s(i-1))**2 * wn(i,4)
      end do
c-
c-------- 11. refinement of s-sequence for random design
      if(trace .ge. 2) call monit1(11, trace)
      if(inputs .and. isrand) then
        do i=0,n
          wn(i,5)=dble(i)/dble(n+1)
          wn(i,2)=(dble(i)+.5)/dble(n+1)
          wn(i,3)=wn(i,2)
        end do
        exs= -dble(3*kord+1) / dble(6*kord+3)
        exsvi=dble(kord)     / dble(6*kord+3)
        bs=0.1*(vi/(sn-s0)**2)**exsvi * n**exs
        call kernel(wn(1,5),t,n,bs,0,2,nyg,wn(0,3),wn(0,2),n+1,s(0),
     .       trace)
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
        if(inputb) goto 160
      end if
      b=bmin*2.
c-
c-------- 12. compute inflation constant and exponent and loop of iterations
      if(trace .ge. 2) call monit1(12, trace)
      const=dble(2*nue+1)*fak2(kord)*vark(kk,nue)*vi
     .       /(dble(2*kord-2*nue)*bias(kk,nue)**2*dble(n))
      fac=1.1*(1.+(nue/10.)+0.05*(kord-nue-2.))
     .       * dble(n)**(2./dble((2*kord+1)*(2*kord+3)))

c     itende=1+2*kord+kord*(2*kord+1)
      itende = (1 + 2*kord) * (1 + kord)
c     ^^^^^^  *fixed* number of iterations ( <== theory !)

      do it=1,itende
c-
c-------- 13. estimate derivative of order kord in iterations
        if(trace .ge. 3) call monit1(13, trace)
        b2 = min(bmax, max(b*fac, bmin/dble(kord-1)*dble(kord+1)))
        call kernel(t,x,n,b2,kord,kord+2,nyg,s,w1(1,1),m1,w1(1,3),trace)
c-
c-------- 14. estimate integralfunctional in iterations
        if(trace .ge. 3) call monit1(14, trace)
        xmy2= .75*(w1(1,2)*w1(1,3)*w1(1,3) + w1(m1,2)*w1(m1,3)*w1(m1,3))
        do i=2,m1-1
           xmy2=xmy2+w1(i,2)*w1(i,3)*w1(i,3)
        end do
        xmy2=xmy2*(tu-tl)/dble(m1)
c-
c-------- 15. finish of iterations
        if(trace .ge. 3) call monit1(15, trace)
        b = min(bmax, max(bmin, (const/xmy2)**ex))
      end do

c-------- 16  compute smoothed function with global plug-in bandwidth
 160  if(trace .ge. 2) call monit1(16, trace)
      call kernel(t,x,n,b,nue,kord,nyg,s,tt,m,y, trace)
c-  return #{iter} (iff aplicable)
      m1=itende
c-------- 17. variance check
      if(trace .ge. 2) call monit1(17, trace)
      if(hetero) sig=rvar
      if(hetero.or. r2.lt.0.88 .or. nue.gt.0) return
      ii=0
      iil=0
      j=2
      tll=max(tl,tt(1))
      tuu=min(tu,tt(m))
      do i=il,iu
         if(t(i).lt.tll .or. t(i).gt.tuu) goto 170
         ii=ii+1
         if(iil.eq.0) iil=i
 171     if(tt(j).lt.t(i)) then
            j=j+1
            if(j.le.m) goto 171
         end if
         wn(ii,3)=x(i)-y(j)+(y(j)-y(j-1))*(tt(j)-t(i))/(tt(j)-tt(j-1))
      end do
 170  continue
      if(iil.eq.0.or.ii-iil.lt.10) then
         call resest(t(il),wn(1,3),nn,wn(1,4),snr,rvar)
      else
         call resest(t(iil),wn(1,3),ii,wn(1,4),snr,rvar)
      end if
      q=sig/rvar
      if(q.le.2.) return
      if(q.gt.5. .and. r2.gt..95) rvar=rvar*.5
      sig=rvar
      if(trace .ge. 2) call monit_s(r2, q, sig, trace)
      call constV(wn(1,4),n,sig)
      goto 100
c     -------- end loop
      end

