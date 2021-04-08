ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     several kernel smoothing subroutines which are used by
cc     glkern.f and lokern.f, version oct 1996
cc
cc     glkerns() & lokerns() directly only call the first 3 and constV()
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     this file contains:
cc
cc     subroutine resest(t,x,n,res,snr,sigma2)
cc                for variance estimation
cc
cc     subroutine kernel(t,x,n,b,nue,kord,ny,s,tt,m,y, trace)
cc                driver subroutine for kernel regression estimation
cc                calls fast or convential kernel routine
cc
cc     subroutine kernp(t,x,n,b,nue,kord,ny,s,tt,m,y)
cc                driver subroutine for kernel regression estimation
cc                without use of boundary kernels
cc---------------------------------------------------------------------
cc     subroutine kernfa(t,x,n,b,nue,kord,ny,s,tt,m,y)
cc                fast algorithm for kernel estimation
cc     subroutine dreg(sw,a1,a2,iord,x,sl,sr,t,b,iflop)
cc                used by subroutine kernfa,kernfp
cc     subroutine lreg(sw,a3,iord,d,dold,q,c)
cc                used by subroutine kernfa,kernfp
cc     subroutine freg(sw,nue,kord,iboun,y,c,icall,a)
cc                used by subroutine kernfa,kernfp
cc     subroutine kernfp(t,x,n,b,nue,kord,ny,s,tt,m,y)
cc                fast algorithm for kernp estimation without boundary
cc---------------------------------------------------------------------
cc     subroutine kerncl(t,x,n,b,nue,kord,ny,s,tt,m,y)
cc                conventional algorithm for kernel estimation
cc     subroutine smo(s,x,n,tau,wid,nue,iord,iboun,ist,s1,c,y)
cc                single estimation step, used by kerncl
cc     subroutine kerncp(t,x,n,b,nue,kord,ny,s,tt,m,y)
cc                conventional algorithm without boundary kernels
cc     subroutine smop(s,x,n,tau,wid,nue,iord,iboun,ist,s1,c,y)
cc                single estimation step, used by kerncp
cc---------------------------------------------------------------------
cc     subroutine coffi(nue,kord,c)
cc                kernel coefficient of polynomial kernels used by
cc                kerncl,kerncp and kernfp
cc     subroutine coffb(nue,kord,q,iboun,c)
cc                kernel coefficient of polynomial boundary kernels
cc                used by kernfa and kerncl
cc---------------------------------------------------------------------
cc     subroutine constV(x,n,fa)
cc                simple subroutine for array initialization
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c==  called from other Fortran *and* from R varNPreg() in ../R/varNPreg.R

      subroutine resest(t,x, n, res, snr,sigma2)
c-----------------------------------------------------------------------
c       version: june, 1996
c
c       purpose:
c
c       computes one-leave-out residuals for nonparametric estimation
c       of residual variance (local linear approximation followed by
c       reweighting)
c
c   Parameters:
c
c   Input
c       t(n)      abscissae (ordered: t(i) <= t(i+1))
c       x(n)      data
c       n         length of data ( >2 )
c   Output
c       res(n)    residuals at t(1),...,t(n)
c       snr       explained variance "R^2" of the true curve
c       sigma2    estimation of sigma^2 (residual variance)
c
c-----------------------------------------------------------------------
      implicit none
c
c     Arguments
      integer n
      double precision x(n),t(n), res(n), snr, sigma2
c Var
      integer i
      double precision ex,ex2,tt,g1,g2,sx,dn

      sigma2=0.
      ex =x(1)*(t(2)-t(1))
      ex2=x(1)*ex
      do i=2,n-1
         tt=t(i+1)-t(i-1)
         if(tt.ne.0.) then
            g1=(t(i+1)-t(i))/tt
         else
            g1=.5
         endif
         g2=1.-g1
         res(i)=(x(i)-g1*x(i-1)-g2*x(i+1))/sqrt(1.+g1*g1+g2*g2)
         sigma2=sigma2+res(i)*res(i)
         sx = x(i)*tt
         ex =ex + sx
         ex2=ex2+x(i)*sx
      end do
c     first points (ex & ex2 done at beginning)
      tt=t(3)-t(2)
      if(tt.ne.0.) then
         g1=(t(1)-t(2))/tt
      else
         g1=.5
      endif
      g2=1.-g1
      res(1)=(x(1)-g1*x(3)-g2*x(2))/sqrt(1.+g1*g1+g2*g2)
c     last points
      tt=t(n-1)-t(n-2)
      if(tt.ne.0.) then
         g1=(t(n-1)-t(n))/tt
      else
         g1=.5
      endif
      g2=1.-g1
      res(n)=(x(n)-g1*x(n-2)-g2*x(n-1))/sqrt(1.+g1*g1+g2*g2)

      sigma2=(sigma2+ res(1)*res(1)+ res(n)*res(n))/n
c  snr := explained variance ( = R^2 )
      sx=x(n)*(t(n)-t(n-1))
      dn=2.*(t(n)-t(1))
      ex =(ex + sx    )/dn
      ex2=(ex2+x(n)*sx)/dn
      if(ex2.gt.0) then
         snr= 1 - sigma2/(ex2-ex*ex)
      else
         snr= 0.
      endif
      return
      end ! resest

      subroutine kernel(t,x,n, b, nue,kord,ny, s,tt,m, y, trace)
c-----------------------------------------------------------------------
c       short-version may, 1995
c
c       driver subroutine for kernel smoothing,
c       chooses between standard and O(n) algorithm
c
c  parameters :
c
c  input    t(n)         input grid (regression design)
c  input    x(n)         data, given on t(n)
c  input    n            length of x
c  input    b            one sided bandwidth (for ny=1 mean bandwidth)
c  input    nue          order of derivative (0-4)
c  input    kord         order of kernel (<=6); default is kord=nue+2
c  input    ny           0: global bandwidth (default)
c                        1: variable bandwidths, given in y as input
c  input    s(0:n)       interpolation sequence
c  input    tt(m)        output grid. must be part of input grid for ieq=0
c  input    m            number of points where function is estimated,
c                         or  length of tt. default is m=400
c  input    y(m)         bandwith sequence for ny=1, dummy for ny=0
c  input    trace        integer: > 0 means "print tracing info"
c  output   y(m)         estimated regression function
c
c-----------------------------------------------------------------------
      implicit none
c Args
      integer n, m
      double precision t(n), x(n), b
      integer nue, kord, ny
      double precision s(0:n), tt(m),y(m)
      integer trace
c Vars
      double precision chan, chR
      integer classic
c
c------  computing change point
      chan=(5.+kord) * max(1., sqrt(float(n)/float(m)))
c------
      chR = chan * (t(n)-t(1)) / (n-1)

      if(trace .gt. 0) then
         if(b .lt. chR) then
            classic = 1
         else
            classic = 0
         end if
         call monitk0(0, n, m, b, ny, chan, chR, classic)
      endif

      if(b .lt. chR) then ! small bandwidth ==> classical kernel
         call kerncl(t,x,n,b, nue,kord,ny, s,tt,m,y, trace)
      else ! fast kernel
         call kernfa(t,x,n,b, nue,kord,ny, s,tt,m,y, trace)
      end if
      return
      end ! kernel


      subroutine kernp(t,x,n,b,nue,kord,ny,s,tt,m,y, trace)
c-----------------------------------------------------------------------
c       short-version january, 1997
c
c       driver subroutine for *local bandwidth* kernel smoothing,
c       chooses between standard and O(n) algorithm
c       without using boundary kernels
c
c  parameters :
c
c  input    t(n)         input grid (regression design)
c  input    x(n)         data, given on t(n)
c  input    n            length of x
c  input    b            one sided bandwidth (for ny=1 mean bandwidth)
c  input    nue          order of derivative (0-4)
c  input    kord         order of kernel (<=6); default is kord=nue+2
c  input    ny           0: global bandwidth (default)
c                        1: variable bandwidths, given in y as input
c  input    s(0:n)       interpolation sequence
c  input    tt(m)        output grid. must be part of input grid for ieq=0
c  input    m            number of points where function is estimated,
c                         or  length of tt. default is m=400
c  input    y(m)         bandwith sequence for ny=1, dummy for ny=0
c  input    trace        integer: > 0 means "print tracing info"
c  output   y(m)         estimated regression function
c
c-----------------------------------------------------------------------
      integer n,nue, kord,ny,m, trace, classic
      double precision t(n),x(n), b,s(0:n),tt(m),y(m)
c
      double precision chan, chR
c
c------  computing change point
      chan=(5.+kord)*max(1.,sqrt(float(n)/float(m)))
c------
      chR = chan * (t(n)-t(1)) / (n-1)

      if(trace .gt. 0) then
         if(b .lt. chR) then
            classic = 1
         else
            classic = 0
         end if
         call monitk0(1, n, m, b, ny, chan, chR, classic)
      endif

      if(b .lt. chR) then
         call kerncp(t,x,n,b,nue,kord,ny,s,tt,m,y, trace)
      else
         call kernfp(t,x,n,b,nue,kord,ny,s,tt,m,y, trace)
      end if
      return
      end ! kernp

      subroutine kernfa(t,x,n,b,nue,kord,ny,s,tt,m,y, trace)
c-----------------------------------------------------------------------
c       short-version: may, 1995
c
c       purpose:
c
c       computation of kernel estimate using O(n) algorithm based on
c       legendre polynomials, general spaced design and local
c       bandwidth allowed. (new initialisations of the legendre sums
c       for numerical reasons)
c
c  parameters :
c
c  input    t(n)         input grid
c  input    x(n)         data, given on t(n)
c  input    n            length of x
c  input    b            one sided bandwidth
c  input    nue          order of derivative (0-4)
c  input    kord         order of kernel (<=6)
c  input    ny           0, global bandwidth; 1, local bandwidth in y
c  input    s(0:n)       half point interpolation sequence
c  input    tt(m)        output grid
c  input    m            number of points to estimate
c  input    y(m)         bandwith sequence for ny=1, dummy for ny=0
c  input    trace        integer: > 0 means "print tracing info"
c  output   y(m)         estimated function
c
c
c-----------------------------------------------------------------------
      implicit none

      integer n,nue,kord,ny,m, trace
      double precision t(n),x(n),s(0:n),tt(m),y(m),b
c Var
      integer j,k,iord,init,icall,i,iboun, jl,jr,jnr,jnl
      double precision c(7),sw(7),xf(7),dold
      double precision a(7,7),a1(7),a2(7),a3(7,7),cm(7,6)
      double precision s0,sn,bmin,bmax,bb,wwl,wwr,wid,wr,wido
c-
c- Shut up over zealous compiler warnings from -Wmaybe-uninitialized :
      wr = 0.
      jl = -1
      jr = -1

      if(trace .gt. 0) call monitfp(0, n, b, nue, kord, ny, m, trace)

c------ compute constants for later use
      s0=1.5*t(1)-0.5*t(2)
      sn=1.5*t(n)-0.5*t(n-1)
      bmin=(sn-s0)*.6d0/dble(n)*dble(kord-1)
      bmax=(s(n)-s(0))*.5
      if(kord.eq.2) bmin=bmin*0.1d0
      iord=kord+1
      do k=3,iord
        a1(k)=dble(2*k-1)/dble(k)
        a2(k)=dble(1-k)/dble(k)
      end do
c-
      init=0
      icall=0
      dold=0.d0
c-
c------ smoothing loop
      do 100 i=1,m
        bb=b
        if (ny .eq. 1) bb=y(i)
        if(bb.lt.bmin) bb=bmin
        if(bb.gt.bmax) bb=bmax
        iboun=0
c-
c------ compute left boundary kernel
        if(tt(i).lt.s(0)+bb) then
          wwl=s(0)
          wwr=s(0)+bb+bb
          wid=wwr-tt(i)
          iboun=1
          call coffb(nue,kord,(tt(i)-s(0))/wid,iboun,c)
        end if
c-
c------ compute right boundary kernel
        if(tt(i)+bb.gt.s(n)) then
          wwl=s(n)-(bb+bb)
          wwr=s(n)
          wid=tt(i)-wwl
          iboun=-1
          call coffb(nue,kord,(s(n)-tt(i))/wid,iboun,c)
        end if
c-
c------ no boundary
        if(iboun.eq.0) then
          wid=bb
          wwl=tt(i)-bb
          wwr=tt(i)+bb
        end if
c-
c------ initialisation for init=0
        if(init.eq.0) then
          do k=1,iord
            sw(k)=0.
          end do
          jl=1
          do j=1,n
            if(s(j-1).lt.wwl) then
               jl=j+1
            else
               if(s(j).gt.wwr) goto 488 ! break
               call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
            end if
          end do
488       jr=j-1
          wr=wwr
          init=1
          goto 6666
        else
          init=init+1
        end if
c-
c------ compare old sum with new smoothing intervall tt(i)-b,tt(i)+b
        if(s(jr-1).ge.wwl) then
          jnr=jr
          jnl=jl
          if(s(jr).gt.wwr) then
            do j=jr,jl,-1
               call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i-1),wido,-1)
               jnr=j-1
               if(s(jnr).le.wwr) goto 2011
            end do
          end if
 2011     continue
          if(s(jl-1).lt.wwl) then
            do j=jl,jr
               call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i-1),wido,-1)
               jnl=j+1
               if(s(j).ge.wwl) goto 3011
            end do
          end if
 3011     continue
c-
c------ updating of sw
          call lreg(sw,a3,iord,(tt(i)-tt(i-1))/wid,dold,wido/wid,cm)
          if(jnr.eq.jr) then
            do j=jr+1,n
              if(s(j).gt.wwr) goto 4011
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
              jnr=j
            end do
          end if
4011      continue
          jr=jnr
          if(jl.eq.jnl) then
            do j=jl-1,1,-1
              if(s(j-1).lt.wwl) goto 4022
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
              jnl=j
            end do
           end if
 4022      continue
          jl=jnl
        else
c-
c------ new initialisation of sw
          do k=1,iord
            sw(k)=0.
          end do
          do j=jr,n
            if(s(j-1).lt.wwl) then
              jl=j+1
            else
              if(s(j).gt.wwr) goto 2022 ! break
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
            end if
          end do
 2022     continue
          jr=j-1
          wr=wwr
        end if
6666    continue
c-
c------ if bandwidth is too small no smoothing
        if(s(jr-1).le.wwl .and. wwr.le.s(jr)) then
          y(i)=x(jr)
          if(nue.gt.0) y(i)=0.d0
        else
c-
c------ add first and last point of the smoothing interval
          do k=1,iord
            xf(k)=sw(k)
          end do
          if(jl.ne.1)
     .       call dreg(xf,a1,a2,iord,x(jl-1),wwl,s(jl-1),tt(i),wid,1)
          if(jr.ne.n)
     .       call dreg(xf,a1,a2,iord,x(jr+1),s(jr),wwr,tt(i),wid,1)
c-
c------ now the sums are built that are needed to compute the estimate
          call freg(xf,nue,kord,iboun,y(i),c,icall,a)
          if(nue.gt.0) y(i)=y(i)/(wid**nue)
        end if
c-
c------ new initialisation ?
        if(jl.gt.jr .or. wwl.gt.wr .or. init.gt.100) init=0
        wido=wid
c-
 100  continue
c-
      return
      end ! kernfa

      subroutine kernfp(t,x,n,b,nue,kord,ny,s,tt,m,y, trace)
c-----------------------------------------------------------------------
c       short-version: january, 1997
c
c       purpose:
c
c       computation of kernel estimate using o(n) algorithm based on
c       legendre polynomials, general spaced design and local
c       bandwidth allowed. (new initialisations of the legendre sums
c       for numerical reasons) without boundary kernels, just normalizing
c
c  parameters :
c
c  input    t(n)         input grid
c  input    x(n)         data, given on t(n)
c  input    n            length of x
c  input    b            one sided bandwidth
c  input    nue          order of derivative (0-4)
c  input    kord         order of kernel (<=6)
c  input    ny           0, global bandwidth; 1, local bandwidth in y
c  input    s(0:n)       half point interpolation sequence
c  input    tt(m)        output grid
c  input    m            number of points to estimate
c  input    y(m)         bandwith sequence for ny=1, dummy for ny=0
c  output   y(m)         estimated function
c
c
c-----------------------------------------------------------------------
      implicit none

      integer n,nue,kord,ny,m, trace
      double precision t(n),x(n),s(0:n),tt(m),y(m),b
c Var
      integer j,k,iord,init,icall,i,iboun, jl,jr,jnr,jnl
      double precision c(7),sw(7),xf(7),dold,qq,q,xnor
      double precision a(7,7),a1(7),a2(7),a3(7,7),cm(7,6)
      double precision s0,sn,bmin,bmax,bb,wwl,wwr,wid,wr,wido
c-
      if(trace .gt. 0) call monitfp(1, n, b, nue, kord, ny, m, trace)
c- Shut up over zealous compiler warnings from -Wmaybe-uninitialized
      jl = -1
      jr = -1
      wr = 0.

c------ compute constants for later use
      s0=1.5*t(1)-0.5*t(2)
      sn=1.5*t(n)-0.5*t(n-1)
      bmin=(sn-s0)*.6d0/dble(n)*dble(kord-1)
      bmax=(s(n)-s(0))*.5
      if(kord.eq.2) bmin=bmin*0.1d0
      iord=kord+1
      call coffi(nue,kord,c)
      do k=3,iord
        a1(k)=dble(2*k-1)/dble(k)
        a2(k)=dble(1-k)/dble(k)
      end do
c-
      init=0
      icall=0
      dold=0.d0
c-
c------ smoothing loop
      do i=1,m
        bb=b
        if (ny .eq. 1) bb=y(i)
        if(bb.lt.bmin) bb=bmin
        if(bb.gt.bmax) bb=bmax
        iboun=0
c-
c------ compute left boundary
        if(tt(i).lt.s(0)+bb) then
          wwl=s(0)
          wwr=s(0)+bb+bb
          wid=wwr-tt(i)
          iboun=1
        end if
c-
c------ compute right boundary
        if(tt(i)+bb.gt.s(n)) then
          wwl=s(n)-(bb+bb)
          wwr=s(n)
          wid=tt(i)-wwl
          iboun=-1
        end if
c-
c------ no boundary
        if(iboun.eq.0) then
          wid=bb
          wwl=tt(i)-bb
          wwr=tt(i)+bb
          xnor=1.d0
        end if
c-
c------ compute normalizing constant
        if(iboun.ne.0) then
          q=0. ! -Wall
          if(iboun.eq.1) q=(tt(i)-s(0))/wid
          if(iboun.eq.-1) q=(s(n)-tt(i))/wid
          qq=q*q
          xnor=c(1)*(1.d0+q)
          do k=3,iord,2
             q=q*qq
             xnor=xnor+c(k)*(1.d0+q)
          end do
          iboun=0
        end if
c-
c------ initialisation for init=0
        if(init.eq.0) then
          do k=1,iord
            sw(k)=0.
          end do
          jl=1
          do j=1,n
            if(s(j-1).lt.wwl) then
              jl=j+1
            else
              if(s(j).gt.wwr) goto 488
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
            end if
          end do
488       jr=j-1
          wr=wwr
          init=1
          goto 6666
        else
          init=init+1
        end if
c-
c------ compare old sum with new smoothing intervall tt(i)-b,tt(i)+b
        if(s(jr-1).ge.wwl) then
          jnr=jr
          jnl=jl
          if(s(jr).gt.wwr) then
            do j=jr,jl,-1
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i-1),wido,-1)
              jnr=j-1
              if(s(jnr).le.wwr) goto 2011
            end do
          end if
2011      continue
          if(s(jl-1).lt.wwl) then
            do j=jl,jr
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i-1),wido,-1)
              jnl=j+1
              if(s(j).ge.wwl) goto 3011
            end do
          end if
3011      continue
c-
c------ updating of sw
          call lreg(sw,a3,iord,(tt(i)-tt(i-1))/wid,dold,wido/wid,cm)
          if(jnr.eq.jr) then
            do j=jr+1,n
              if(s(j).gt.wwr) goto 4011
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
              jnr=j
            end do
          end if
4011      continue
          jr=jnr
          if(jl.eq.jnl) then
            do j=jl-1,1,-1
              if(s(j-1).lt.wwl) goto 4022
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
              jnl=j
            end do
          end if
4022      continue
          jl=jnl
        else
c-
c------ new initialisation of sw
          do k=1,iord
            sw(k)=0.
          end do
          do j=jr,n
            if(s(j-1).lt.wwl) then
              jl=j+1
            else
              if(s(j).gt.wwr) goto 2022
              call dreg(sw,a1,a2,iord,x(j),s(j-1),s(j),tt(i),wid,1)
            end if
          end do
2022      jr=j-1
          wr=wwr
        end if
6666    continue
c-
c------ if bandwidth is too small no smoothing
        if(s(jr-1).le.wwl .and. wwr.le.s(jr)) then
          y(i)=x(jr)
          if(nue.gt.0) y(i)=0.d0
        else
c-
c------ add first and last point of the smoothing interval
          do k=1,iord
            xf(k)=sw(k)
          end do
          if(jl.ne.1)
     .       call dreg(xf,a1,a2,iord,x(jl-1),wwl,s(jl-1),tt(i),wid,1)
          if(jr.ne.n)
     .       call dreg(xf,a1,a2,iord,x(jr+1),s(jr),wwr,tt(i),wid,1)
c-
c------ now the sums are built that are needed to compute the estimate
          call freg(xf,nue,kord,iboun,y(i),c,icall,a)
          if(nue.gt.0) y(i)=y(i)/(wid**nue)
          if(xnor.ne.1.d0) y(i)=y(i)/xnor
        end if
c-
c------ new initialisation ?
        if(jl.gt.jr .or. wwl.gt.wr .or. init.gt.100) init=0
        wido=wid
c-
      end do
c-    ------ end{smoothing loop}
      return
      end ! kernfp

      subroutine dreg(sw,a1,a2, iord, x,sl,sr, t,b, iflop)
c-----------------------------------------------------------------------
c       version: may, 1995
c
c       purpose:
c
c       computes new legendre sums (for regression)
c
c       parameters:
c                    **************   input   *******************
c
c        sw(iord)  :   old sum of data weights for legendre polynom.
c        a1(7)     :   constants of recursive formula for legendre pol.
c        a2(7)     :                             "
c        iord      :   order of kernel polynomial
c        x         :   data point
c        sl        :   left s-value
c        sr        :   right s-value
c        t         :   point where the smoothed value is to be estimated
c        b         :   bandwidth
c        iflop     :   1: addition, else subtraction
c
c
c                    **************   output   *******************
c
c        sw(iord)  :   new sum of data weights for legendre polynom.
c
      implicit none
c-----------------------------------------------------------------------

      double precision sw(7), a1(7), a2(7)
      integer iord, iflop
      double precision x,sl,sr, t,b
c Var
      integer k
      double precision p(7,2)

c------  compute legendre polynomials
      p(1,1)=(t-sl)/b
      p(1,2)=(t-sr)/b
      p(2,1)=1.5d0*p(1,1)*p(1,1)-.5d0
      p(2,2)=1.5d0*p(1,2)*p(1,2)-.5d0
      do k=3,iord
         p(k,1)=a1(k)*p(k-1,1)*p(1,1)+ a2(k)*p(k-2,1)
         p(k,2)=a1(k)*p(k-1,2)*p(1,2)+ a2(k)*p(k-2,2)
      end do

c------compute new legendre sums
      if(iflop.eq.1) then
        do k=1,iord
          sw(k)=sw(k)+(p(k,1)-p(k,2))*x
        end do
      else
        do k=1,iord
          sw(k)=sw(k)+(p(k,2)-p(k,1))*x
        end do
      end if
      return
      end

      subroutine lreg(sw,a3, iord, d,dold, q,c)
c------------------------------------------------------------------
c       version: may, 1995
c
c       purpose:
c
c       update of sw-sequence according to new bandwidth and new data
c       (version for regression)
c
c       parameters :
c                          **************   input   *******************
c
c               sw(iord)  :  sum of data weights for legendre polynom.
c               iord      :  order of kernel polynomial
c               d         :  dist. to the next point divided by bandw.
c               dold      :  d       previous step
c               q         :  new bandwidth divided by old bandwidth
c                          **************    work   *******************
c
c               a3(7,7)   :  matrix (p*q*p)**(-1)
c               c(7,6)    :  matrix of coefficients
c                          **************   output   ******************
c
c               sw(iord)  :  updated version of sw
c
      implicit none
c---------------------------------------------------------------------
      double precision sw(7), a3(7,7)
      integer iord
      double precision d,dold, q, c(7,6)
c Var
      integer k,i,l
      double precision dd,ww,qq,xx

c- build up matrix
      if(dold.ne.d .or. dold.eq.0) then
        dold=d
        dd=d*d
c-
        if(iord.eq.7) then
          c(7,6)=13.d0*d
          c(7,5)=71.5d0*dd
          c(7,4)=(214.5d0*dd+9.d0)*d
          c(7,3)=(375.375d0*dd+77.d0)*dd
          c(7,2)=((375.375d0*dd+247.5d0)*dd+5.d0)*d
          c(7,1)=((187.6875d0*dd+346.5d0)*dd+40.5d0)*dd
        end if
c-
        if(iord.ge.6) then
          c(6,5)=11.d0*d
          c(6,4)=49.5d0*dd
          c(6,3)=(115.5d0*dd+7.d0)*d
          c(6,2)=(144.375d0*dd+45.d0)*dd
          c(6,1)=((86.625d0*dd+94.5d0)*dd+3.d0)*d
        end if
c-
        if(iord.ge.5) then
          c(5,4)=9.d0*d
          c(5,3)=31.5d0*dd
          c(5,2)=(52.5d0*dd+5.d0)*d
          c(5,1)=(39.375d0*dd+21.d0)*dd
        end if
c-
        if(iord.ge.4) then
          c(4,3)=7.d0*d
          c(4,2)=17.5d0*dd
          c(4,1)=(17.5d0*dd+3.d0)*d
        end if
c-
        if(iord.ge.3) then
          c(3,2)=5.d0*d
          c(3,1)=7.5d0*dd
        end if
c-
        c(2,1)=3.d0*d
      end if
      if(q.lt..9999.or.q.gt.1.0001) then
c-
c------- built up matrix a3=p*q*p**-1
        a3(1,1)=q
        do k=2,iord
          a3(k,k)=a3(k-1,k-1)*q
        end do
        ww=q*q-1.d0
        do k=1,iord-2
          ww=ww*q
          a3(k+2,k)=(k+.5d0)*ww
        end do
c-
        qq=0. ! -Wall (stupid compiler)
        if(iord.ge.5) then
          qq=a3(2,2)
          a3(5,1)=q*(1.875d0+qq*(-5.25d0+qq*3.375d0))
        end if
        if(iord.ge.6) a3(6,2)=qq*(4.375d0+qq*(-11.25d0+qq*6.875d0))
        if(iord.eq.7) then
          a3(7,1)=q*(-2.1875d0+qq*(11.8125d0+qq*(-18.5625d0+qq*
     .          8.9375d0)))
          a3(7,3)=q*qq*(7.875d0+qq*(-19.25d0+qq*11.375d0))
        end if
c-
c------- compute a*c and new legendre sums
        do i=iord,2,-1
          xx=0.
          do k=1,i
            ww=0.
            do l=k,i-1,2
               ww=ww+a3(l,k)*c(i,l)
            end do
            if(mod(i-k,2).eq.0) ww=ww+a3(i,k)
            xx=xx+ww*sw(k)
          end do
          sw(i)=xx
        end do
        sw(1)=a3(1,1)*sw(1)
      else
        do i=iord,2,-1
           do k=1,i-1
              sw(i)=sw(i)+c(i,k)*sw(k)
           end do
        end do
      end if
      return
      end

      subroutine freg(sw, nue,kord,iboun, y,c, icall, a)
c------------------------------------------------------------------
c       short-version: may, 1995
c
c       purpose:
c
c       final computation of a smoothed value via legendre polynomials
c
c       parameters :
c                          **************   input   *******************
c
c               sw(kord+1):  sum of data weights for legendre polynom.
c               nue       :  order of derivative
c               kord      :  order of kernel
c               iboun     :  0: interior kernel, else boundary kernel
c               c(kord+1) :  sequence of polyn. coeff. for bound. kernel
c               icall     :  parameter used to initialise computation
c                         :   of a matrix
c                          **************    work    ******************
c
c               a(7,7)    :  matrix of coefficients
c                          **************   output   ******************
c
c               y          :  computed estimate
c
      implicit none
c--------------------------------------------------------------------
      double precision sw(7)
      integer nue, kord, iboun
      double precision y, c(7)
      integer icall
      double precision a(7,7)
c Var
      integer i,j
      double precision ww

c------- definition of legendre coefficients for boundary
      if(icall.eq.0.and.iboun.ne.0) then
             a(2,2)=2./3.
         a(1,3)=.6
                 a(3,3)=.4
             a(2,4)=4./7.
                     a(4,4)=8./35.
         a(1,5)=27./63.
                 a(3,5)=28./63.
                         a(5,5)=8./63.
             a(2,6)=110./231.
                     a(4,6)=72./231.
                             a(6,6)=16./231.
         a(1,7)=143./429.
                 a(3,7)=182./429.
                         a(5,7)=88./429.
                                 a(7,7)=16./429.
         icall=1
      end if
      if(iboun.ne.0) then
c-
c------- computation of the smoothed value at boundary
        y=c(1)*sw(1)+c(2)*a(2,2)*sw(2)
        do j=3,kord+1
          ww=a(j,j)*sw(j)
          do i=j-2,1,-2
            ww=ww+a(i,j)*sw(i)
          end do
          y=y+c(j)*ww
        end do
c-
      else
c------- computation of the smoothed value at interior
        if(nue.eq.0) then
          if(kord.eq.2) y=-.1*sw(3)+.6*sw(1)
          if(kord.eq.4) y=(sw(5)-4.*sw(3)+9.*sw(1))/12.
          if(kord.eq.6) y=-7.2115379e-02*sw(7)+.25961537*sw(5)
     .                    -.4375*sw(3)+.75*sw(1)
        end if
        if(nue.eq.1) then
          if(kord.eq.3) y=(3.*sw(4)-10.*sw(2))/14.
          if(kord.eq.5) y=(-15.*sw(6)+48.*sw(4)-55.*sw(2))/44.
        end if
        if(nue.eq.2) then
          if(kord.eq.4) y=(-5.*sw(5)+14.*sw(3)-9.*sw(1))/6.
          if(kord.eq.6) y=2.01923*sw(7)-5.76923*sw(5)+5.25*sw(3)
     .                   -1.5*sw(1)
        end if
        if(nue.eq.3) y=4.772727*sw(6)-12.272727*sw(4)+7.5*sw(2)
        if(nue.eq.4) y=-36.34615*sw(7)+88.84615*sw(5)-52.5*sw(3)
      end if
c-
      return
      end

      subroutine kerncl(t,x, n, b, nue,kord,ny, s,tt,m,y, trace)
c-----------------------------------------------------------------------
c       short-version january 1995
c
c       kernel smoothing, conventional algorithm,general
c       design, local bandwidth allowed
c
c  parameters :
c
c  input    t(n)         input grid
c  input    x(n)         data, given on t(n)
c  input    n            length of x
c  input    b            one sided bandwidth (for ny=1 mean bandwidth)
c  input    nue          order of derivative (0-4)
c  input    kord         order of kernel (<=6)
c  input    s(0:n)       half point interpolation sequence
c  input    ny           0: global, 1: local bandwidth inputed in y
c  input    tt(m)        output grid
c  input    m            number of points to estimate
c  input    y(m)         bandwith sequence for ny=1, dummy for ny=0
c  output   y(m)         estimated regression function
c
      implicit none
c-----------------------------------------------------------------------

      integer n, m
      double precision t(n), x(n), b
      integer nue, kord, ny
      double precision s(0:n), tt(m), y(m)
      integer trace
c Var
      double precision c(7),c1(7)
      integer ist,i,iboun,iord
      double precision bb, s0,s1,sn,bmin,bmax, wid
c-
c      if(trace .gt. 0) call intpr('  kerncl()',-1, 0,0)

c------  compute kernel coefficients for interior and some constants
      call coffi(nue,kord,c)
      iord=kord+1
      bb=b
      s0=1.5*t(1)-0.5*t(2)
      sn=1.5*t(n)-0.5*t(n-1)
      bmin=(sn-s0)*.6d0/dble(n)*dble(kord-1)
      bmax=(s(n)-s(0))*.5
      if(kord.eq.2) bmin=0.1d0*bmin
      ist=1
c-
c------- Loop over output grid ------------------------------
      do i=1,m
        if(ny.ne.0) bb=y(i)
        if(bb.gt.bmax) bb=bmax
        if(bb.lt.bmin) bb=bmin
        wid=bb
        s1=tt(i)-bb
        iboun=0
c-
c-------  compute left boundary kernel
        if(s1.lt.s(0)) then
          s1=s(0)
          wid=bb+bb+s(0)-tt(i)
          call coffb(nue,kord,(tt(i)-s(0))/wid,1,c1)
          iboun=1
        end if
c-
c-------  compute right boundary kernel
        if(tt(i)+bb.gt.s(n)) then
          s1=s(n)-(bb+bb)
          wid=tt(i)-s1
          call coffb(nue,kord,(s(n)-tt(i))/wid,-1,c1)
          iboun=-1
        end if
c-
c------  search first s-point of smoothing interval
c
c MM: "internal logic error":  ist is used as index for x(.) and s(.)
c     but s = s(0:n)  whereas  x = x(1:n)
c     --> below, we cannot allow ist < 1
2       if(s(ist).le.s1) then
          ist=ist+1
          goto 2
        end if
3       if(ist .gt. 1) then
           if(s(ist-1) .gt. s1) then
              ist=ist-1
              goto 3
           end if
        end if
c-
        if(s(ist).ge.tt(i)+wid .or. ist.eq.n) then
c         if bandwidth is too small no smoothing
           y(i)=x(ist)
           if(nue.gt.0) y(i)=0.
        else ! compute smoothed data at tt(i)
           if (iboun.ne.0) then ! boundary kernel
              call smo(s,x,n,tt(i),wid,nue,iord,iboun,ist,s1,
     .             c1,y(i),trace)
           else
              call smo(s,x,n,tt(i),wid,nue,iord,iboun,ist,s1,
     .             c, y(i),trace)
           end if
        end if
      end do
c-
      return
      end ! kerncl()

      subroutine smo(s,x,n, tau,wid, nue,iord,iboun,ist, s1,c,y, trace)
c-----------------------------------------------------------------------
c       short-version january 1995
c
c       performs one smoothing step
c
c  parameters :
c
c  input    s(0:n)       half point interpolation sequence
c  input    x(n)         data
c  input    n            length of x
c  input    tau          point where function is estimated
c  input    wid          one sided bandwidth
c  input    nue          order of derivative (0-4)
c  input    iord          order of kernel polynomial
c  input    iboun        type of boundary (-1: right, +1, left; 0: *no* bndry)
c  input    ist          index of first point of smoothing interval
c  input    s1           left boundary of smoothing interval
c  input    c(7)         kernel coefficients
c  OUTPUT   y            smoothed value at tau
c  work     wo(7)        work array
c
      implicit none
c-----------------------------------------------------------------------
      integer n
      double precision s(0:n), x(n)
      double precision tau, wid
      integer nue, iord, iboun, ist
      double precision s1, c(7), y
      integer trace
c Var
      double precision wo(7), yy,yyy,w
      integer jend,ibeg,incr,i,j
      logical nu_odd
c-
      nu_odd = (nue.eq.1 .or. nue.eq.3)
      y=0.
      jend=0
      if(iboun.ne.0 .or. .not.nu_odd) then
         ibeg=1
      else
         ibeg=2
      end if

      if(iboun.ne.0) then
         incr=1
      else
         incr=2
      end if
      if(trace .ge. 2) call monits(tau, ist, n, iboun)
c-
c------  compute initial kernel values
      if(iboun.gt.0) then
        yy=(tau-s1)/wid
        wo(ibeg)=yy
        do i=ibeg,iord-incr,incr
          wo(i+incr)=wo(i)*yy
        end do
      else
        do i=ibeg,iord,incr
          wo(i)=1.
        end do
      end if
c-
c------  loop over smoothing interval
      do j=ist,n
        yy=(tau-s(j))/wid
        if(yy.lt.-1.) then
          yy=-1.
          jend=1
        end if
        yyy=yy
        if(iboun.eq.0) then
          yy=yy*yy
          if(nu_odd) yyy=yy
        end if

c       loop for computing weights
        w=0.
        do i=ibeg,iord,incr
          w=w+c(i)*(wo(i)-yyy)
          wo(i)=yyy
          yyy=yyy*yy
        end do
        y=y+w*x(j)
        if(jend.eq.1) goto 110 ! break
      end do
 110  continue

c     -- normalizing for nue > 0
      if(nue.gt.0) y= y/(wid**nue)
c-
      return
      end ! smo

      subroutine kerncp(t,x,n, b, nue,kord,ny, s,tt,m, y, trace)
c-----------------------------------------------------------------------
c       short-version january 1997
c
c       kernel smoothing, conventional algorithm,general
c       design, local bandwidth allowed, without boundary kernels,
c       just normalizing
c
c  parameters :
c
c  input    t(n)         input grid
c  input    x(n)         data, given on t(n)
c  input    n            length of x
c  input    b            one sided bandwidth (for ny=1 mean bandwidth)
c  input    nue          order of derivative (0-4)
c  input    kord         order of kernel (<=6)
c  input    s(0:n)       half point interpolation sequence
c  input    ny           0: global, 1: local bandwidth inputed in y
c  input    tt(m)        output grid
c  input    m            number of points to estimate
c  input    y(m)         bandwith sequence for ny=1, dummy for ny=0
c  output   y(m)         estimated regression function
c
      implicit none
c-----------------------------------------------------------------------
      integer n, m
      double precision t(n),x(n), b
      integer nue, kord, ny
      double precision s(0:n), tt(m), y(m)
      integer trace
c Var
      integer ist,i,iboun,iord
      double precision c(7), c1(7), bb, bmax, wid, s0,s1,sn, bmin

c------  compute kernel coefficients for interior and some constants
      call coffi(nue,kord,c)
      iord=kord+1
      bb=b
      s0=1.5*t(1)-0.5*t(2)
      sn=1.5*t(n)-0.5*t(n-1)
      bmin=(sn-s0)*.6d0/dble(n)*dble(kord-1)
      bmax=(s(n)-s(0))*.5
      if(kord.eq.2) bmin=0.1d0*bmin
      ist=1
c-
c-------  loop over output grid
      do i=1,m
        if(ny.ne.0) bb=y(i)
        if(bb.gt.bmax) bb=bmax
        if(bb.lt.bmin) bb=bmin
        wid=bb
        s1=tt(i)-bb
        iboun=0
c-
c-------  compute left boundary kernel
        if(s1.lt.s(0)) then
          s1=s(0)
          wid=bb+bb+s(0)-tt(i)
          call coffb(nue,kord,(tt(i)-s(0))/wid,1,c1)
          iboun=1
        end if
c-
c-------  compute right boundary kernel
        if(tt(i)+bb.gt.s(n)) then
          s1=s(n)-(bb+bb)
          wid=tt(i)-s1
          iboun=-1
        end if
c-
c------  search first s-point of smoothing interval
2       if(s(ist).le.s1) then
          ist=ist+1
          goto 2
        end if
3       if(s(ist-1).gt.s1) then
          ist=ist-1
          goto 3
        end if
c-
c-------  if bandwidth is too small no smoothing
        if(s(ist).ge.tt(i)+wid.or.ist.eq.n) then
         y(i)=x(ist)
         if(nue.gt.0) y(i)=0.
        else
c-
c-----  compute smoothed data at tt(i)
          call smop(s,x,n,tt(i),wid,nue,iord,iboun,ist,s1,c,y(i),trace)
        end if
      end do
c-
      return
      end ! kerncp()

      subroutine smop(s,x,n, tau,wid, nue,iord,iboun,ist, s1,c,y, trace)
c-----------------------------------------------------------------------
c       short-version january 1995
c
c       performs one smoothing step
c
c  parameters :
c
c  input    s(0:n)       half point interpolation sequence
c  input    x(n)         data
c  input    n            length of x
c  input    tau          point where function is estimated
c  input    wid          one sided bandwidth
c  input    nue          order of derivative (0-4)
c  input    iord          order of kernel polynomial
c  input    iboun        type of boundary (-1: right, +1, left; 0: *no* bndry)
c  input    ist          index of first point of smoothing interval
c  input    s1           left boundary of smoothing interval
c  input    c(7)         kernel coefficients
c  OUTPUT   y            smoothed value at tau
c  work     wo(7)        work array
c
      implicit none
c-----------------------------------------------------------------------

      integer n
      double precision s(0:n), x(n),  tau, wid
      integer nue, iord, iboun, ist
      double precision s1, c(7), y
      integer trace
c Var
      double precision wo(7), yy,yyy,w,ww
      integer jend,ibeg,incr,i,j
      logical nu_odd
c-
      nu_odd = (nue.eq.1 .or. nue.eq.3)
      y=0.
      ww=0.
      jend=0
      if(nu_odd) then
         ibeg=2
      else
         ibeg=1
      end if
      incr=2
      if(trace .ge. 2) call monits(tau, ist, n, iboun)
c-
c------  compute initial kernel values
      if(iboun.gt.0) then
        yy=(tau-s1)/wid
        wo(ibeg)=yy
        yy=yy*yy
        if(nu_odd) wo(ibeg)=yy
        do i=ibeg,iord-incr,incr
          wo(i+incr)=wo(i)*yy
        end do
      else
        do i=ibeg,iord,incr
          wo(i)=1.
        end do
      end if
c-
c------  loop over smoothing interval
      do j=ist,n
        yy=(tau-s(j))/wid
        if(yy.lt.-1.) then
          yy=-1.
          jend=1
        end if
        yyy=yy
        yy=yy*yy
        if(nu_odd) yyy=yy

c       loop for computing weights
        w=0.
        do i=ibeg,iord,incr
          w=w+c(i)*(wo(i)-yyy)
          wo(i)=yyy
          yyy=yyy*yy
        end do
        y=y+w*x(j)
        ww=ww+w
        if(jend.eq.1) goto 110 ! break
      end do
 110  continue

      if(ww.ne.0) y= y/ww
c     -- normalizing for nue > 0
      if(nue.gt.0) y= y/(wid**nue)
c-
      return
      end ! smop

      subroutine coffi(nue,kord,c)
c-----------------------------------------------------------------------
c       short-version: january 1995
c
c       purpose:
c
c       defines polynomial kernel coefficients for interior.
c
c  parameters:
c
c  input  nue        order of derivative (0-4)
c  input  kord       order of kernel (nue+i, i=2,4,6;  kord<=6)
c  output c(7)       polynomial kernel coefficients
c
      implicit none
c-----------------------------------------------------------------------
      integer nue,kord
      double precision c(7)
c
      integer i

      do i=1,7
        c(i)=0.
      end do
      if(nue.eq.0 .and. kord.eq.2) then
          c(1)=0.75d0
          c(3)=-0.25d0
      end if
c
      if(nue.eq.0.and.kord.eq.4) then
        c(1)=1.40625d0
        c(3)=-1.5625d0
        c(5)=0.65625d0
      end if
c
      if(nue.eq.0.and.kord.eq.6) then
        c(1)=2.05078125d0
        c(3)=-4.78515625d0
        c(5)=5.16796875d0
        c(7)=-1.93359375d0
      end if
c
      if(nue.eq.1.and.kord.eq.3) then
        c(2)=-1.875d0
        c(4)=0.9375d0
      end if
c
      if(nue.eq.1.and.kord.eq.5) then
        c(2)=-8.203125d0
        c(4)=11.484375d0
        c(6)=-4.921875d0
      end if
c
      if(nue.eq.2.and.kord.eq.4) then
        c(1)=-6.5625d0
        c(3)=13.125d0
        c(5)=-6.5625d0
      end if
c
      if(nue.eq.2.and.kord.eq.6) then
        c(1)=-24.609375d0
        c(3)=103.359375d0
        c(5)=-132.890625d0
        c(7)=54.140625d0
      end if
c
      if(nue.eq.3.and.kord.eq.5) then
        c(2)=88.59375d0
        c(4)=-147.65625d0
        c(6)=68.90625d0
      end if
c
      if(nue.eq.4.and.kord.eq.6) then
        c(1)=324.84375d0
        c(3)=-1624.21875d0
        c(5)=2273.90625d0
        c(7)=-974.53125d0
      end if
c
      return
      end ! coffi

      subroutine coffb(nue,kord,q,iboun,c)
c-----------------------------------------------------------------------
c       short-version: january 1995
c
c       purpose:
c
c       computes coefficients of polynomial boundary kernels,
c       following Gasser + Mueller preprint 38 SFB 123, Heidelberg
c       and unpublished results
c
c  parameters:
c
c  input  nue        order of derivative (0-4)
c  input  kord       order of kernel (nue+i, i=2,4,6;  kord<=6)
c  input  q          percentage of wid at boundary
c  input  iboun      < 0 right boundary of data
c                    > 0 left boundary of data
c  output c(7)       polynomial kernel coefficients
c
c-----------------------------------------------------------------------
c Arguments
      integer nue,kord,iboun
      double precision q, c(7)
c Var
      integer i,j
      double precision p,p1,p3,p12,p6,d
c
      do i=1,7
         c(i)=0.
      end do
      p=-q
      p1=1.+q
      p3=p1*p1*p1
c
c Compute c(j) depending on derivative and kernel order  (nue, kord) :
c
      if(nue.eq.0.and.kord.eq.2) then
        d=1./(p3*p1)
        c(1)=(6.+p*(12.+p*18.))*d
        c(2)=9.*(1.+p)*(1.+p)*d
        c(3)=(4.+p*8.)*d
      end if
c
      if(nue.eq.0.and.kord.eq.4) then
        d=p1/(p3*p3*p3)
        p12=(1.+p)*(1.+p)*d
        c(1)=20.*(1.+p*(12.+p*(78.+p*(164.+p*(165.+p*(60.+p*10.))))))*d
        c(2)=100.*(1.+p*(5.+p))**2*p12
        c(3)=200.*(1.+p*(12.+p*(33.+p*(36.+p*(14.+p+p)))))*d
        c(4)=175.*(1.+p*(10.+p*3.))*p12
        c(5)=56.*(1.+p*(12.+p*(18.+p*4.)))*d
      end if
c
      if(nue.eq.0.and.kord.eq.6) then
        p6=p3*p3
        d=1./(p6*p6)
        p12=(1.+p)*(1.+p)*d
        c(1)=42.*(1.+p*(30.+p*(465.+p*(3000.+p*(10050.+p*(17772.+p
     .       *(17430.+p*(9240.+p*(2625.+p*(350.+p*21.))))))))))*d
        c(2)=441.*(1.+p*(14.+p*(36.+p*(14.+p))))**2*p12
        c(3)=1960.*(1.+p*(30.+p*(255.+p*(984.+p*(1902.+p*(1956.+p
     .       *(1065.+p*(300.+p*(39.+p+p)))))))))*d
        c(4)=4410.*(1.+p*(28.+p*(156.+p*(308.+p*(188.+p*(42.+p*3.))))))
     .       *p12
        c(5)=5292.*(1.+p*(30.+p*(185.+p*(440.+p*(485.+p*(250.+p*(57.
     .       +p*4.)))))))*d
        c(6)=3234.*(1.+p*(28.+p*(108.+p*(56.+p*5.))))*p12
        c(7)=792.*(1.+p*(30.+p*(150.+p*(200.+p*(75.+p*6.)))))*d
       end if
c
      if(nue.eq.1.and.kord.eq.3) then
        d=-1./(p3*p3)
        p12=(1.+p)*(1.+p)*d
        c(1)=(60.+p*240.)*p12
        c(2)=120.*(2.+p*(6.+p*(6.+p)))*d
        c(3)=300.*p12
        c(4)=(120.+p*180.)*d
      end if
c
      if(nue.eq.1.and.kord.eq.5) then
        d=-1./(p3*p3*p3*p1)
        p12=(1.+p)*(1.+p)*d
        c(1)=420.*(1.+p*(18.+p*(98.+p*(176.+p*(75.+p*10.)))))*p12
        c(2)=2100.*(2.+p*(25.+p*(120.+p*(245.+p*(238.+p*(105.+p*(20.
     .       +p)))))))*d
        c(3)=14700.*(1.+p*(4.+p))**2*p12
        c(4)=5880.*(4.+p*(35.+p*(90.+p*(95.+p*(40.+p*6.)))))*d
        c(5)=17640.*(1.+p*(6.+p+p))*p12
        c(6)=2520.*(2.+p*(15.+p*(20.+p*5.)))*d
      end if
c
      if(nue.eq.2.and.kord.eq.4) then
        d=p1/(p3*p3*p3)
        p12=(1.+p)*(1.+p)*d
        c(1)=840.*(1.+p*(12.+p*(28.+p*(24.+p*5.))))*d
        c(2)=2100.*(3.+p*(10.+p))*p12
        c(3)=1680.*(9.+p*(28.+p*(27.+p*6.)))*d
        c(4)=14700.*p12
        c(5)=(5040.+p*6720.)*d
      end if
c
      if(nue.eq.2.and.kord.eq.6) then
        p6=p3*p3
        d=1./(p6*p6)
        p12=(1.+p)*(1.+p)*d
        c(1)=5040.*(2.+p*(60.+p*(489.+p*(1786.+p*(3195.+p*(2952.+p
     .       *(1365.+p*(294.+p*21.))))))))*d
        c(2)=52920.*(3.+p*(42.+p*(188.+p*(308.+p*(156.+p*(28.
     .       +p))))))*p12
        c(3)=141120.*(6.+p*(68.+p*(291.+p*(570.+p*(555.+p*(264.+p
     .       *(57.+p*4.)))))))*d
        c(4)=529200.*(2.+p*(7.+p+p))**2*p12
        c(5)=90720.*(30.+p*(228.+p*(559.+p*(582.+p*(255.
     .       +p*40.)))))*d
        c(6)=582120.*(3.+p*(14.+p*5.))*p12
        c(7)=221760.*(2.+p*(12.+p*(15.+p*4.)))*d
       end if
c
      if(nue.eq.3.and.kord.eq.5) then
        d=-1./(p3*p3*p3*p1)
        p12=(1.+p)*(1.+p)*d
        c(1)=15120.*(1.+p*(18.+p*(38.+p*6.)))*p12
        c(2)=45360.*(4.+p*(35.+p*(80.+p*(70.+p*(20.+p)))))*d
        c(3)=352800.*(2.+p*(6.+p))*p12
        c(4)=151200.*(8.+p*(25.+p*(24.+p*6.)))*d
        c(5)=952560.*p12
        c(6)=70560.*(4.+p*5.)*d
      end if
c
      if(nue.eq.4.and.kord.eq.6) then
        p6=p3*p3
        d=1./(p6*p6)
        p12=(1.+p)*(1.+p)*d
        c(1)=332640.*(1.+p*(30.+p*(171.+p*(340.+p*(285.+p*(90.+p*7.
     .      ))))))*d
        c(2)=1164240.*(5.+p*(56.+p*(108.+p*(28.+p))))*p12
        c(3)=6652800.*(5.+p*(38.+p*(85.+p*(76.+p*(25.+p+p)))))*d
        c(4)=17463600.*(5.+p*(14.+p*3.))*p12
        c(5)=4656960.*(25.+p*(78.+p*(75.+p*20.)))*d
        c(6)=76839840.*p12
        c(7)=3991680.*(5.+p*6.)*d
      end if
c
      if(iboun.gt.0) return
      j=2
      if(nue.eq.1.or.nue.eq.3) j=1
      do i=j,kord,2
        c(i)=-c(i)
      end do
      return
      end ! coffb

      subroutine constV(x,n,fa)
      integer n
      double precision  x(n),fa

      integer i
      do i=1,n
        x(i)=fa
      end do
      return
      end
