      subroutine lpridge_s(t,x,n,b,nue,p,kord,wk,tt,m,mnew,imoms,moms,y,
     .              leng,nmoms,nvar,var,ridge,nsins)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     local polynomials and local polynomial ridge regression with
c     polynomial weights for regression functions and derivatives
c
c     22. 02. 94    B. Seifert                    last edit: 10. 04. 95
c
c     Calculates in O(n) steps the fit and its variance
c
c     input    t(n)     inputgrid
c     input    x(n)     data
c     input    n        number of data
c     in/out   b(m)     bandwidth
c     input    nue      derivative to be estimated (0 <= nue <= p)
c     input    p        polynomial approximation order (0 <= p <= 10)
c     input    kord     order of kernel weights used (0 <= kord <= 12-p)
c     input    wk(0:kord) coefficients of kernel weights
c     input    tt(m)    outputgrid
c     input    m        size of outputgrid
c     input    mnew     restart parameter: restart is forced after mnew
c                       updating steps. (mnew >= 0)
c                       mnew = 0 --> no updating
c     in/out   imoms(nmoms) 0,1 - block (not) computed, used in lpnew.
c                       Set imoms = 0 at the beginning and after
c                       changes other than b, nue
c     in/out   moms(nmoms,4*(max(kord,2)+p+nvar))
c                       moments for blocks, used in lpnew.
c                       Do not change or set imoms = 0.
c     output   y(m)     estimated function
c     input    leng     length of blocks stored. Proposal: leng = 10
c     input    nmoms    dimension of arrays. nmoms >= n/leng
c     input    nvar     computation of variance 0 = no, 1 = yes
c     output   var(m)   variance of estimator / sigma^2
c     input    ridge    ridge parameter (ridge >= 0)
c                       ridge = 0 --> local polynomial regression
c     output   nsins    number of singularity exceptions
c
c     calls:            lpadd, lpnew, lpslv, lpsub, lpsv
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer n,nue,p,kord,m,mnew,imoms(*),leng,nmoms,nvar,nsins
      double precision t(*),x(*),b(*),wk(0:*),tt(*),moms(nmoms,2,0:*),
     .   y(*),var(*),ridge

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     workarrays and local parameters
c
c     input    nmini    minimal number of points in smoothing interval
c     input    pmax     dimension of arrays, where
c                       pmax >= 2*(1+max(kord-2,0)+p+nvar)
c                       p <= 10 & p+kord <= 12 --> pmax = 24
c     work     w(0:pmax,5) - (*,1)=s_ (*,2)=t_ (*,4)=S_ (*,5)=T_
c     work     w1((p+1)*(p+1)) matrix for lpslv
c     work     bin(0:pmax,0:pmax) binomial coefficients
c     work     work(2*pmax) work array for lpnew, lpslv and lpsv
c     work     bb(0:powmax) powers of 1/b(i)
c     work     tti(0:powmax) powers of (tt(i)-tbar)
c     work     wk2(0:2*kord) coefficients of squared kernel weights
c     input    sin(2)   relative bounds for singularity.
c                       (1) = update, (2) = restart. (0 < sin < 1)
c                       Proposal: sin = (.99,.01)
c     input    zer      absolute lower bound for Cholesky-factors
c                       (zer > 0)
c                       Proposal: zer = 1d-10*n*((p+1)/(2*n))**(p+p-1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision w(0:24,5),w1(121),bin(0:24,0:24),work(48),
     .   bb(0:24),tti(0:24),wk2(0:24),sin(2),zer
      integer nmini,pmax

      integer again,dif,i,iaux,io,ioold,irec,iu,iuold,iup,j,k,kk,l,m2,
     .   na,nmin,nsin,nsub,nzer,pow,powmax
      double precision chol(20,2),nuefak,sino,sins(2,2,2),
     .   tbar,tk,tkt,tleft,tright,tttj,tttl,xbar,xsin,xx,ysin
c - ideal Cholesky-factors for uniform and Epanechnikov weights
      data chol /1.d0,          1.d0,  5.33333333d-01,2.28571429d-01,
     .   8.70748299d-02,3.07840308d-02,1.03331012d-02,3.33838656d-03,
     .   1.04733696d-03,3.21010399d-04,9.65444809d-05,2.85835627d-05,
     .   8.35137136d-06,2.41261839d-06,6.90199908d-07,1.95774166d-07,
     .   5.51153184d-08,1.54132174d-08,4.28476026d-09,1.18479489d-09,
     .           1.d0,          1.d0,  4.44444444d-01,1.60000000d-01,
     .   5.22448980d-02,1.61249685d-02,4.79751129d-03,1.39099440d-03,
     .   3.95660629d-04,1.10894501d-04,3.07186985d-05,8.42848643d-06,
     .   2.29433279d-06,6.20387584d-07,1.66798306d-07,4.46249887d-08,
     .   1.18876192d-08,3.15473249d-09,8.34421263d-10,2.19982450d-10/

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c - internal and parameters for numerical stability and speed
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       nmini=1
       pmax=24
       sin(1)=.99d0
       sin(2)=1d-2
       zer=1d-10
       if (p.gt.0) zer=1d-10*dble(n)*(dble(p+1)/dble(n+n))**(p+p-1)

c - limits for singularity: sins(i,*,*)=lower and upper,
c      (*,i,*)=update and restart, (*,*,i)=weighted and unweighted
       do 10 j=1,2
          sins(2,1,j)=chol(p+1,j)/sin(1)
          sins(2,2,j)=chol(p+1,j)/zer
          do 10 i=1,2
10           sins(1,i,j)=chol(p+1,j)*sin(i)
c - maximum order of moments
      pow=p+max(kord-2,0)
      if(nvar.gt.0) pow=pow+1
      powmax=2+2*pow
      na=p+1
c - nuefak = faktorial of nue
      nuefak=1d0
      do 20 i=1,nue
20       nuefak=nuefak*i
c - binomial coefficients
      do 40 kk=0,powmax
         bin(kk,0)=1
         do 30 k=1,kk-1
30          bin(kk,k)=bin(kk-1,k-1)+bin(kk-1,k)
40       bin(kk,kk)=1
c - squared kernel weights
      if (nvar.le.0) goto 100
      do 50 j=0,kord+kord
50       wk2(j)=0d0
      do 60 k=0,kord
         wk2(k+k)=wk2(k+k)-wk(k)*wk(k)
         do 60 l=0,k
60          wk2(k+l)=wk2(k+l)+2d0*wk(k)*wk(l)
c**********************************************************************
c Smoothing loop over outputgrid
c
c - loop-variables: i - index of output grid
c                   k - number of points in smoothing interval
c**********************************************************************
100   nsins=0
      iup=mnew
      nsub=1
      iuold=0
      ioold=0
      iu=1
      io=2
      m2=(m+1)/2
      do 990 iaux=1,m
         if (iaux.lt.m2) then
c - left half of smoothing interval is done from the left to the right
            i=iaux
         else
            if (iaux.gt.m2) then
c - right half of smoothing interval is done from the right to the left
               i=m+m2-iaux
            else
c - if center of smoothing interval is reached,
c      switch to the right boundary and restart
               i=m
               iup=mnew
               iu=n
               io=n+1
            endif
         endif
c - minimal number of points in the smoothing interval
         nmin=max(nmini,p+1)
c
c - determine smoothing interval
c
110      tleft=tt(i)-b(i)
         tright=tt(i)+b(i)
c - determine indices corresponding to the smoothing interval
c   - first left boundary (iu)
c     - check whether we have to move left
120      if ((iu.gt.1).and.(t(iu-1).gt.tleft)) then
            iu=iu-1
            goto 120
         endif
c     - check whether we have to move right
130      if ((iu.le.n).and.(t(iu).lt.tleft)) then
            iu=iu+1
            goto 130
         endif
c   - now the right boundary (io)
c     - check whether we have to move left
140      if ((io.gt.1).and.(t(io-1).gt.tright)) then
            io=io-1
            goto 140
         endif
c     - check whether we have to move right
150      if ((io.le.n).and.(t(io).lt.tright)) then
            io=io+1
            goto 150
         endif
c - Now iu and io have their new values
c - check if there are enough points in smoothing interval
160      if (io-iu.lt.nmin) then
            again=0
170         if (iu.gt.1) then
               if (io.le.n) then
                  if (abs(tt(i)-t(io)).lt.abs(tt(i)-t(iu-1))) then
                     tkt=t(io)
                     io=io+1
                  else
                     iu=iu-1
                     tkt=t(iu)
                  endif
               else
                  iu=iu-1
                  tkt=t(iu)
               endif
            else
               if (io.le.n) then
                  tkt=t(io)
                  io=io+1
               else
                  nsins=nsins+1
                  y(i)=0d0
                  goto 990
               endif
            endif
            if (io-iu.lt.nmin) goto 170
            if (again.eq.0) then
               tk=tkt
               again=1
               goto 170
            endif
            b(i)=(abs(tt(i)-tk)+abs(tt(i)-tkt))/2
            goto 110
         endif
c - Now we have enough points
c
c - compute sums
c
200      iup=iup+1
c - Restart Sum ?
         if(iup.gt.mnew.and.(nsub.gt.0.or.iuold.lt.iu.or.io.lt.ioold))
     .      then
c - restart sum!
            iup=1
            nsub=0
c   - sum over the smoothing interval [tt(i)-b(i),tt(i)+b(i)]
            call lpnew(w(0,2),w(0,1),t,x,tbar,xbar,leng,nmoms,
     .              imoms,moms,work,pow,pmax,k,bin,iu,io-1)
c - Now the sum is (re)calculated.
c      t(iu) is the leftmost point in the smoothing intervall
c      and t(io) is the leftmost point above the smoothing interval
         else
c - update sum
c   - interval has moved too far -> restart sum
            if(2*iuold+ioold-1.lt.3*iu.or.3*io-1.lt.iuold+2*ioold) then
               iup=mnew
               goto 200
            endif
c   - subtract terms no longer needed
c     - at left boundary
            if (iuold.lt.iu) then
               call lpsub(w(0,2),w(0,1),t,x,tbar,xbar,pow,pmax,k,bin,
     .                 iuold,iu-1)
               nsub=nsub+1
            endif
c     - at right boundary
            if (io.lt.ioold) then
               call lpsub(w(0,2),w(0,1),t,x,tbar,xbar,pow,pmax,k,bin,
     .                 io,ioold-1)
               nsub=nsub+1
            endif
c   - add new terms
c     - at left boundary
            if (iu.lt.iuold)
     .         call lpadd(w(0,2),w(0,1),t,x,tbar,xbar,pow,pmax,k,bin,
     .                 iu,iuold-1)
c     - at right boundary
            if (ioold.lt.io)
     .         call lpadd(w(0,2),w(0,1),t,x,tbar,xbar,pow,pmax,k,bin,
     .                 ioold,io-1)
         endif
c - Now the sum has its new value
c
c compute and solve LSE
c
c - calculate powers of b and (tbar-tt(i))
         bb(0)=1d0
         tti(0)=1d0
         do 310 kk=1,powmax
            bb(kk)=bb(kk-1)/b(i)
310         tti(kk)=tti(kk-1)*(tbar-tt(i))
c - for p+1 points or kord=0 use unweighted LSE
         if (io-iu.le.p+1) iup=mnew
         if (io-iu.le.p+1.or.kord.eq.0) goto 400
c - weighted LSE
c    - calculate Sn,j
         do 330 j=0,2*p
            w(j,4)=wk(0)*w(j,1)
            do 330 kk=1,kord
               xx=tti(kk)*w(j,1)
               do 320 l=1,kk
320               xx=xx+bin(kk,l)*tti(kk-l)*w(j+l,1)
330         w(j,4)=w(j,4)+bb(kk)*wk(kk)*xx
c    - calculate Tn,j
         do 350 j=0,p
            w(j,5)=wk(0)*w(j,2)+w(j,4)*xbar
            do 350 kk=1,kord
               xx=tti(kk)*w(j,2)
               do 340 l=1,kk
340               xx=xx+bin(kk,l)*tti(kk-l)*w(j+l,2)
350            w(j,5)=w(j,5)+bb(kk)*wk(kk)*xx
c    - construct matrix
         do 370 j=0,p
            do 370 l=1,j+1
370            w1(j*na+l)=w(l+j-1,4)
         irec=1
         goto 500
c - unweighted LSE
c    - store Sn,j for ridging
400      do 410 j=0,1
410         w(j,4)=w(j,1)
c    - store T_n,l
         do 420 l=0,p
420         w(l,5)=w(l,2)+w(l,1)*xbar
c    - construct matrix
         do 430 j=0,p
            do 430 l=1,j+1
430            w1(j*na+l)=w(l+j-1,1)
         irec=2
c - add ridge-matrix
500      if (ridge.le.0d0) goto 550
         xx=tt(i)-tbar-w(1,4)/w(0,4)
         tttj=xx
         do 540 j=nue+1,p
            tttl=xx
            do 530 l=nue+1,j
               w1(j*na+l+1)=w1(j*na+l+1)
     .            +ridge*bin(j,nue)*tttj*bin(l,nue)*tttl
530            tttl=tttl*xx
540         tttj=tttj*xx
c    - solve linear equations
550      if (nsub.eq.0) then
            xsin=sins(1,2,irec)
            ysin=sins(2,2,irec)
         else
            xsin=sins(1,1,irec)
            ysin=sins(2,1,irec)
         endif
         dif=p+1-nue
         call lpslv(w1,work,w(0,5),na,nsin,nzer,sino,xsin,zer,dif)
         if (sino.gt.ysin) nsin=nsin+1
         if (nsin+nzer.le.0) goto 600
c the matrix is singular!
c - 1. try restart
         if (nsub.gt.0) then
            iup=mnew
            goto 200
         endif
c - 1.b singularity caused by ridging
         if (sino.gt.ysin) then
            nsins=nsins+1
            y(i)=0d0
            goto 990
         endif
c - 2. try unweighted LSE
         if (irec.eq.1) goto 400
c - 3. try more points
         ioold=io
         iuold=iu
         nmin=nmin+1
         goto 160
c - limits for singularity: sins(i,*,*)=lower and upper,
c      (*,i,*)=update and restart, (*,*,i)=weighted and unweighted
600      if (nsin.eq.0.and.nsub.eq.0) then
            sins(2,1,irec)=sino/sin(1)
            sins(1,1,irec)=sino*sin(1)
         endif
c - compute estimate y(i)
         y(i)=0d0
         kk=1
         do 640 j=nue,p
            y(i)=y(i)+kk*bin(j,nue)*tti(j-nue)*w(j,5)
640         kk=-kk
         y(i)=y(i)*nuefak
c - store old values
         ioold=io
         iuold=iu
c - the estimator is numerically stable
c
c - compute variance
c
         if (nvar.le.0) goto 990
c     - w(,3) = u, u from y = u' * beta
         do 710 j=0,nue
710         w(j,3)=0d0
         kk=1
         do 740 j=nue,p
            w(j,3)=kk*bin(j,nue)*tti(j-nue)*nuefak
740         kk=-kk
c    - w(,3) = (H + X'W X)**-1 * u
         call lpsv(w1,work,w(0,3),na,nzer,sino,xsin,zer,na)
         if (irec.eq.2) goto 850
c    - w(,4) = X'W**2 X
c    - weighted LSE
         do 830 j=0,2*p
            w(j,4)=wk2(0)*w(j,1)
            do 830 kk=1,kord+kord
               xx=tti(kk)*w(j,1)
               do 820 l=1,kk
820               xx=xx+bin(kk,l)*tti(kk-l)*w(j+l,1)
830            w(j,4)=w(j,4)+bb(kk)*wk2(kk)*xx
         goto 880
c    - unweighted LSE
850      do 860 l=0,2*p
860         w(l,4)=w(l,1)
c    - variance (sigma = 1)
c    - var = u'* (H + X'W X)**-1 X'W**2 X (H + X'W X)**-1 * u
880      var(i)=0d0
         do 890 j=0,p
            var(i)=var(i)-w(j+j,4)*w(j,3)*w(j,3)
            do 890 l=0,j
890            var(i)=var(i)+2d0*w(l+j,4)*w(j,3)*w(l,3)

c**********************************************************************
c - end of smoothing loop over outputgrid
c**********************************************************************
990   continue
c
      return
      end
