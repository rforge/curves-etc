      subroutine lpepa(t,x,n,b,nue,p,tt,m,mnew,imoms,moms,y,
     .              leng,nmoms,nvar,var)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     local polynomials with Epanechnikov(!) weights for regression
c     functions and derivatives
c
c     04. 01. 93    B. Seifert                    last edit: 19. 08. 93
c
c                                (hopefully) cosmetic edits: 23. 03. 95
c
c     Calculates in O(n) steps the fit and its variance
c
c     input    t(n)     inputgrid
c     input    x(n)     data
c     input    n        number of data
c     in/out   b(m)     bandwidth
c     input    nue      derivative to be estimated (0 <= nue <= p)
c     input    p        polynomial approximation order (0 <= p <= 10)
c     input    tt(m)    outputgrid
c     input    m        size of outputgrid
c     input    mnew     restart parameter: restart is forced after mnew
c                       updating steps. (mnew >= 0)
c                       mnew = 0 --> no updating
c     in/out   imoms(nmoms) 0,1 - block (not) computed, used in lpnew.
c                       Set imoms = 0 at the beginning and after
c                       changes other than b, nue
c     in/out   moms(nmoms,4*(2+p+nvar))
c                       moments for blocks, used in lpnew.
c                       Do not change or set imoms = 0.
c     output   y(m)     estimated function
c     input    leng     length of blocks stored. Proposal: leng = 10
c     input    nmoms    dimension of arrays. nmoms >= n/leng
c     input    nvar     computation of variance 0 = no, 1 = yes
c     output   var(m)   variance of estimator / sigma^2
c
c     calls:            lpadd, lpnew, lpslv, lpsub, lpsv
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer n,nue,p,m,mnew,imoms(*),leng,nmoms,nvar
      double precision t(*),x(*),b(*),tt(*),moms(nmoms,0:*),y(*),var(*)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     workarrays and local parameters
c
c     input    pmax     dimension of arrays,
c                       where pmax >= 2+2*p+2*nvar --> pmax = 24
c     work     w(0:pmax,5) - (*,1)=s_ (*,2)=t_ (*,4)=S_ (*,5)=T_
c     work     w1((p+1)*(p+1)) matrix for lpslv
c     work     bin(0:pmax,0:pmax) binomial coefficients
c     work     work(2*pmax) work array for lpslv and lpnew
c     input    sin(2)   relative bounds for singularity - (1) = update,
c                       (2) = restart. Proposal: sin = (.99,.01)
c     input    zer      absolute lower bound for Cholesky-factors.
c                       Proposal: zer = 1d-10*n*((p+1)/(2*n))**(p+p-1)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision w(0:24,5),w1(121),bin(0:24,0:24),work(48),
     .   sin(2),zer
      integer pmax

      integer again,dif,i,iaux,io,ioold,irec,iu,iuold,iup,j,k,kk,l,m2,
     .   na,nmin,nsin,nsub,nzer,pow,powmax
      double precision bb,chol(20,2),nuefak,sino,sins(2,2,2),
     .   tbar,tk,tkt,tleft,tright,tti,ttti,xbar,xsin,ysin
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

       pmax=24
       sin(1)=.99d0
       sin(2)=1d-2
       zer=1d-10
       if (p.gt.0) zer=1d-10*dble(n)*(dble(p+1)/dble(n+n))**(p+p-1)

c - limits for singularity: sins(i,*,*)=lower and upper,
c      (*,i,*)=update and restart, (*,*,i)=weighted and unweighted
       do 10 i=1,2
          do 10 j=1,2
             sins(2,i,j)=2d0
             if (sin(i).gt.0d0) sins(2,i,j)=chol(p+1,j)/sin(i)
10           sins(1,i,j)=chol(p+1,j)*sin(i)
c - maximum order of moments
      powmax=2+2*p
      if(nvar.gt.0) powmax=powmax+2
      pow=p
      if(nvar.gt.0) pow=pow+1
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
c**********************************************************************
c Smoothing loop over outputgrid
c
c - loop-variables: i - index of output grid
c                   k - number of points in smoothing interval
c**********************************************************************
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
         nmin=p+1
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
               tkt=t(io)
               io=io+1
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
c - for p+1 points use unweighted LSE
         if (io-iu.le.p+1) goto 400
c - weighted LSE
c    - calculate Sn,l
         bb=1d0/(b(i)*b(i))
         tti=tt(i)-tbar
         do 310 l=0,2*p
310         w(l,4)=w(l,1)-bb*(tti*tti*w(l,1)-2d0*tti*w(l+1,1)+w(l+2,1))
c    - calculate Tn,l
         do 320 l=0,p
320         w(l,5)=w(l,2)-bb*(tti*tti*w(l,2)-2d0*tti*w(l+1,2)+w(l+2,2))
c    - construct matrix
         do 360 j=0,p
            do 360 l=1,j+1
360            w1(j*na+l)=w(l+j-1,4)
         irec=1
         goto 500
c - unweighted LSE
c    - store t_n,l
400      do 420 l=0,p
420         w(l,5)=w(l,2)
c    - construct matrix
         do 430 j=0,p
            do 430 l=1,j+1
430            w1(j*na+l)=w(l+j-1,1)
         iup=mnew
         irec=2
c    - solve linear equations
500      if (nsub.eq.0) then
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
            sins(2,1,irec)=2d0
            if (sin(1).gt.0d0) sins(2,1,irec)=sino/sin(1)
            sins(1,1,irec)=sino*sin(1)
         endif
c - compute estimate y(i)
         y(i)=0d0
         tti=tt(i)-tbar
         ttti=1d0
         do 640 j=nue,p
            y(i)=y(i)+bin(j,nue)*ttti*w(j,5)
640         ttti=ttti*tti
         if (nue.eq.0) y(i)=y(i)+xbar
         y(i)=y(i)*nuefak
c - store old values
         ioold=io
         iuold=iu
c - the estimator is numerically stable
c
c - compute variance
c
         if (nvar.le.0) goto 990
         do 710 j=0,nue
710         w(j,3)=0d0
         ttti=1d0
         do 740 j=nue,p
            w(j,3)=bin(j,nue)*ttti*nuefak
740         ttti=ttti*tti
         call lpsv(w1,work,w(0,3),na,nsin,nzer,sino,xsin,zer,na)
         if (irec.eq.2) goto 850
c    - weighted LSE
         bb=1d0/(b(i)*b(i))
         do 810 l=0,2*p
810         w(l,4)=w(l,1)
     .         -2d0*bb*(tti*tti*w(l,1)-2d0*tti*w(l+1,1)+w(l+2,1))
     .         +bb*bb*(tti*tti*tti*tti*w(l,1)-4d0*tti*tti*tti*w(l+1,1)
     .         +6d0*tti*tti*w(l+2,1)-4d0*tti*w(l+3,1)+w(l+4,1))
         goto 880
c    - unweighted LSE
850      do 860 l=0,2*p
860         w(l,4)=w(l,1)
c    - variance (sigma = 1)
c    - var = u'* (X'W X)**-1 X'W**2 X (X'W X)**-1 * u
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
      subroutine lpadd(t,s,to,x,tbar,xbar,p,pmax,n,bin,iu,io)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   21. 12. 92 B. Seifert                         last edit: 22. 12. 92
c
c     Calculates new centered moments when one observation added
c
c     in/out   t(0:p)	vector of centered moments multiplied with
c                       centered y-values, i.e. t_{j,n}
c     in/out   s(0:p)	vector of centered moments, i.e. s_{j,n}
c     input    to(*)    inputgrid
c     input    x(*)     data
c     in/out   tbar     mean of inputgrid-points
c     in/out   xbar     mean of data
c     input    p	number of moments
c     input    pmax     true dimension of arrays
c     input    n        number of data used
c     input    bin(0:p,0:p)  matrix of binomial coefficients
c     input    iu,io    intervall of data used
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer n,p,ps,j,k,pmax,iu,io,i
      double precision s(0:pmax),t(0:pmax),to(*),x(*),tbar,xbar,
     .         e,d,dp,sx,tx,bin(0:pmax,0:pmax),dnp,nn
c
      ps=2+2*p
      do 500 i=iu,io
         n=n+1
         d=tbar
         e=xbar
         tbar=tbar+(to(i)-tbar)/n
         xbar=xbar+(x(i)-xbar)/n
         d=d-tbar
         e=e-xbar
         nn=1-n
         do 200 j=p+2,1,-1
		dp=1d0
		dnp=1d0
		tx=0d0
		do 100 k=j,1,-1
			tx=tx+bin(j,k)*(t(k)+e*s(k))*dp
			dp=dp*d
			dnp=dnp*nn
100             continue
		t(j)=tx-e*dp*nn*(1d0-dnp)
200      continue
C
         do 400 j=ps,2,-1
		dnp=nn
		dp=1d0
		sx=0d0
		do 300 k=j,2,-1
			sx=sx+bin(j,k)*s(k)*dp
			dp=dp*d
			dnp=dnp*nn
300		continue
		dp=dp*d
		s(j)=sx+dp*(dnp-nn)
400	continue
500   continue
      s(0)=dble(n)
      return
      end
      subroutine lpnew(t,s,to,x,tbar,xbar,leng,nmoms,imoms,moms,mom,
     .              p,pmax,n,bin,iu,io)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     04. 01. 93    B. Seifert                    last edit: 14. 01. 93
c
c     Calculates new centered moments at restart
c
c     output   t(0:p+2) vector of centered moments multiplied with
c                       centered y-values, i.e. t_{j,n}
c     output   s(0:2+2p) vector of centered moments, i.e. s_{j,n}
c     input    to(*)    inputgrid
c     input    x(*)     data
c     output   tbar     mean of inputgrid-points
c     output   xbar     mean of data
c     input    leng     length of blocks stored
c     in/out   imoms(*) 0,1 - block (not) computed
c     in/out   moms(nmoms,2,0:2+2p) moments for blocks - (*,1,*)=t
c                       (*,2,*)=s (*,*,0)=means
c     work     mom(*)   local moms
c     input    p        number of moments
c     input    pmax     true dimension of arrays
c     output   n        number of data used
c     input    bin(0:p,0:p)  matrix of binomial coefficients
c     input    iu,io    intervall of data used
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer leng,nmoms,imoms(*),p,pmax,n,iu,io
      double precision t(0:pmax),s(0:pmax),to(*),x(*),tbar,xbar,
     .   moms(nmoms,2,0:pmax),mom(2,*),bin(0:pmax,0:pmax)

      integer ps,j,k,im1,im2,i1,i2,ii,ii1,ii2,i,n2
      double precision e,d,dp,sx,tx,dnp,xy,xx,tt,n1,nn
c
      n=0
      ps=2+2*p
      im1=(iu+leng-2)/leng+1
      im2=io/leng
      i1=(im1-1)*leng
      if(i1.gt.io) i1=io
      i2=im2*leng+1
      if(i2.le.i1) i2=i1+1
c         write(*,'(''restart'',7i4)')iu,io,ps,im1,im2,i1,i2
      tbar=0d0
      xbar=0d0
      do 10 j=0,ps
         t(j)=0d0
10       s(j)=0d0
c
c compute first incomplete block without storing
c
      if(iu.le.i1) then
         n=i1-iu+1
         do 110 i=iu,i1
             tbar=tbar+to(i)
110          xbar=xbar+x(i)
         tbar=tbar/n
         xbar=xbar/n
         do 140 i=iu,i1
            xy=x(i)-xbar
            do 130 j=1,p+2
               xy=xy*(to(i)-tbar)
130            t(j)=t(j)+xy
            xy=to(i)-tbar
            do 140 j=2,ps
               xy=xy*(to(i)-tbar)
140            s(j)=s(j)+xy
      endif
c
c compute central blocks
c
      if(im2.ge.im1) then
         do 480 i=im1,im2
            if(imoms(i).le.0) then
c - store block i
               imoms(i)=1
               ii1=(i-1)*leng+1
               ii2=i*leng
               do 210 j=1,ps
                  moms(i,1,j)=0d0
210               moms(i,2,j)=0d0
               tt=0d0
               xx=0d0
               do 220 ii=ii1,ii2
                  tt=tt+to(ii)
220               xx=xx+x(ii)
               tt=tt/leng
               xx=xx/leng
               moms(i,1,0)=tt
               moms(i,2,0)=xx
               do 240 ii=ii1,ii2
                  xy=x(ii)-xx
                  do 230 j=1,p+2
                     xy=xy*(to(ii)-tt)
230                  moms(i,1,j)=moms(i,1,j)+xy
                  xy=to(ii)-tt
                  do 240 j=2,ps
                     xy=xy*(to(ii)-tt)
240                  moms(i,2,j)=moms(i,2,j)+xy
            endif
c
            n1=n
            n=n+leng
            d=tbar
            e=xbar
            tbar=tbar+(moms(i,1,0)-tbar)*leng/n
            xbar=xbar+(moms(i,2,0)-xbar)*leng/n
            d=d-tbar
            e=e-xbar
            nn=-n1/leng
            do 320 j=p+2,1,-1
               dp=1d0
               dnp=1d0
               tx=0d0
               do 310 k=j,1,-1
                  tx=tx+bin(j,k)*(t(k)+e*s(k)
     .               +dnp*(moms(i,1,k)+nn*e*moms(i,2,k)))*dp
                  dp=dp*d
310               dnp=dnp*nn
320               t(j)=tx+e*dp*n1*(1d0-dnp)
C
            do 340 j=ps,2,-1
               dnp=1d0
               dp=1d0
               sx=0d0
               do 330 k=j,2,-1
                  sx=sx+bin(j,k)*(s(k)+dnp*moms(i,2,k))*dp
                  dp=dp*d
330               dnp=dnp*nn
               dp=dp*d
340            s(j)=sx+dp*n1*(1d0-dnp)
480      continue
      endif
c
c compute last incomplete block without storing
c
      if(io.ge.i2) then
         tt=0d0
         xx=0d0
         n2=io-i2+1
         do 510 i=i2,io
             tt=tt+to(i)
510          xx=xx+x(i)
         tt=tt/n2
         xx=xx/n2
         do 520 j=1,ps
            mom(1,j)=0d0
520         mom(2,j)=0d0
         do 540 i=i2,io
            xy=x(i)-xx
            do 530 j=1,p+2
               xy=xy*(to(i)-tt)
530            mom(1,j)=mom(1,j)+xy
            xy=to(i)-tt
            do 540 j=2,ps
               xy=xy*(to(i)-tt)
540            mom(2,j)=mom(2,j)+xy
c
         n1=n
         n=n+n2
         d=tbar
         e=xbar
         tbar=tbar+(tt-tbar)*n2/n
         xbar=xbar+(xx-xbar)*n2/n
         d=d-tbar
         e=e-xbar
         nn=-n1/n2
         do 620 j=p+2,1,-1
            dp=1d0
            dnp=1d0
            tx=0d0
            do 610 k=j,1,-1
               tx=tx+bin(j,k)*(t(k)+e*s(k)
     .            +dnp*(mom(1,k)+nn*e*mom(2,k)))*dp
               dp=dp*d
610            dnp=dnp*nn
620            t(j)=tx+e*dp*n1*(1d0-dnp)
C
         do 640 j=ps,2,-1
            dnp=1d0
            dp=1d0
            sx=0d0
            do 630 k=j,2,-1
               sx=sx+bin(j,k)*(s(k)+dnp*mom(2,k))*dp
               dp=dp*d
630            dnp=dnp*nn
            dp=dp*d
640         s(j)=sx+dp*n1*(1d0-dnp)
      endif
c
      s(0)=dble(n)
c         write(*,'(5g12.4)')(t(i),i=0,ps)
c         write(*,'(5g12.4)')(s(i),i=0,ps)
c
      return
      end
        SUBROUTINE lpslv(A,D,y,NA,NSIN,nzer,SINout,sin,ZER,dif)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       28. 12. 92    B. Seifert                  last edit: 19. 05. 93
c
C       CHOLESKY-DECOMPOSITION and SOLUTION OF LINEAR EQUATION
C       FORMULA: R(T)*D*R=A
C       UPPER TRIANGULAR OF NONNEGATIVE A IS USED
C       R TRANSPOSED (R(T)) STORED AS LOWER TRIANGULAR
c          WITHOUT DIAG (=1) ON A
C       then solve: R(T)*D*R*X=Y
C       RESULTING X OVERWRITES Y
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none
        double precision A(NA,NA),D(NA),y(na),SINout,sin,ZER
        INTEGER dif,NA,NSIn,nzer
        double precision XX,YY
        INTEGER II,JJ,KK,NN1,NN2
C
C CHOLESKY-DECOMPOSITION
c
        NSIN=0
        DO 130 II=1,NA
           D(II)=A(II,II)
           NN1=II-1
           IF (NN1) 100,100,30
30         DO 90 JJ=1,NN1
              XX=A(JJ,II)
              NN2=JJ-1
              IF (NN2) 60,60,40
40            DO 50 KK=1,NN2
                 XX=XX-A(II,KK)*A(JJ,KK)*D(KK)
50            CONTINUE
60            YY=0D0
              IF (D(JJ)-ZER) 80,80,70
70            YY=XX/D(JJ)
80            D(II)=D(II)-XX*YY
              A(II,JJ)=YY
90         CONTINUE
100        IF (D(II)-A(II,II)*SIN) 110,110,130
110        D(II)=0D0
           NSIN=NSIN+1
130     CONTINUE
        sinout=d(na)/a(na,na)
C
C SOLUTION OF LINEAR EQUATION
C
        NZER=0
        IF (NA-1) 300,300,210
210     DO 230 II=2,NA
           NN1=II-1
           DO 220 KK=1,NN1
              Y(II)=Y(II)-A(II,KK)*Y(KK)
220        CONTINUE
230     CONTINUE
c - compute only last dif elements
300     DO 390 II=1,dif
           JJ=NA+1-II
           XX=0D0
           IF (D(JJ)-ZER) 320,320,310
310        XX=Y(JJ)/D(JJ)
320        IF (D(JJ)-ZER) 330,330,350
330        IF (D(JJ)) 340,350,340
340        NZER=NZER+1
350        NN1=JJ+1
           IF (NA-NN1) 380,360,360
360        DO 370 KK=NN1,NA
              XX=XX-A(KK,JJ)*Y(KK)
370        CONTINUE
380        Y(JJ)=XX
390     CONTINUE
C
        RETURN
        END
      subroutine lpsub(t,s,to,x,tbar,xbar,p,pmax,n,bin,iu,io)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   21. 12. 92   B. Seifert                       last edit: 22. 12. 92
c
c     Calculates new centered moments when one observation dropped
c		 (corresonds to t{j,n-1}, s_{j,n-1}, 0<=j<=p)
c
c     in/out   t(0:p)	vector of centered moments multiplied with
c				centered y-values, i.e. t_{j,n-1}
c     in/out   s(0:p)	vector of centered moments, i.e. s_{j,n-1}
c     input    to(*)    inputgrid
c     input    x(*)     data
c     in/out   tbar     mean of inputgrid-point
c     in/out   xbar     mean of data
c     input    p	number of moments
c     input    pmax     true dimension of arrays
c     in/out   n        number of data
c     input    bin(0:p,0:p)   matrix of binomialcoefficients
c     input    iu,io    interval of data used
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer n,p,ps,j,k,pmax,iu,io,i
      double precision s(0:pmax),t(0:pmax),to(*),x(*),tbar,xbar,
     .              e,d,dp,sx,tx,bin(0:pmax,0:pmax),dnp,nn
c
      ps=2+2*p
      do 500 i=iu,io
         d=tbar
         e=xbar
         tbar=tbar-(to(i)-tbar)/(n-1)
         xbar=xbar-(x(i)-xbar)/(n-1)
         d=d-tbar
         e=e-xbar
         nn=n
         do 100 j=p+2,1,-1
		dp=1d0
		dnp=1d0
		tx=0d0
		do 200 k=j,1,-1
			tx=tx+bin(j,k)*(t(k)+e*s(k))*dp
			dp=dp*d
			dnp=dnp*nn
200		continue
		t(j)=tx+e*dp*nn*(1d0-dnp)
100       continue
c
          do 400 j=ps,2,-1
		dnp=nn
		dp=1d0
		sx=0d0
		do 300 k=j,2,-1
			sx=sx+bin(j,k)*s(k)*dp
			dp=dp*d
			dnp=dnp*nn
300		continue
		dp=dp*d
		s(j)=sx+dp*(nn-dnp)
400        continue
        n=n-1
500   continue
      s(0)=dble(n)
      return
      end

        SUBROUTINE lpsv(A,D,y,NA,NSIN,nzer,SINout,sin,ZER,dif)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       19. 05. 93    B. Seifert                  last edit: 21. 05. 93
c
C       SOLUTION OF LINEAR EQUATION without CHOLESKY-DECOMPOSITION
C       FORMULA: R(T)*D*R=A
C       UPPER TRIANGULAR OF NONNEGATIVE A IS USED
C       R TRANSPOSED (R(T)) STORED AS LOWER TRIANGULAR
c          WITHOUT DIAG (=1) ON A
C       solve: R(T)*D*R*X=Y
C       RESULTING X OVERWRITES Y
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        implicit none
        double precision A(NA,NA),D(NA),y(na),SINout,sin,ZER
        INTEGER dif,NA,NSIn,nzer
        double precision XX
        INTEGER II,JJ,KK,NN1
C
C SOLUTION OF LINEAR EQUATION
C
        IF (NA-1) 300,300,210
210     DO 230 II=2,NA
           NN1=II-1
           DO 220 KK=1,NN1
              Y(II)=Y(II)-A(II,KK)*Y(KK)
220        CONTINUE
230     CONTINUE
c - compute only last dif elements
300     DO 390 II=1,dif
           JJ=NA+1-II
           XX=0D0
           IF (D(JJ)-ZER) 350,350,310
310        XX=Y(JJ)/D(JJ)
350        NN1=JJ+1
           IF (NA-NN1) 380,360,360
360        DO 370 KK=NN1,NA
              XX=XX-A(KK,JJ)*Y(KK)
370        CONTINUE
380        Y(JJ)=XX
390     CONTINUE
C
        RETURN
        END
