C---------- LP auxiliary routines   used BOTH in   lpridge()  and  lpepa() ---------

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
        integer NA,NSIn,nzer,dif
        double precision A(NA,NA),D(NA),y(na),SINout,sin,ZER

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
        INTEGER NA,NSIn,nzer,dif
        double precision A(NA,NA),D(NA),y(na),SINout,sin,ZER

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
