       SUBROUTINE LOKERNS(T,X,N,TT,M,NUE,KORD,IHOM,IRND,
     .        ISMO,M1,TL,TU,S,SIG,WN,W1,WM,BAN,Y)
C----------------------------------------------------------------------*
C-----------------------------------------------------------------------
C       SHORT-VERSION: JANUARY 1997
C
C       PURPOSE:
C
C       GENERAL SUBROUTINE FOR KERNEL SMOOTHING:
C       COMPUTATION OF ITERATIVE PLUG-IN ALGORITHM FOR LOCAL BANDWIDTH
C       SELECTION FOR KERNELS WITH (NUE,KORD) = (0,2),(0,4),(1,3) OR 
C       (2,4).
C-----------------------------------------------------------------------
C  USED SUBROUTINES: COFF, RESEST, KERNEL WITH FURTHER SUBROUTINES
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),T(N),TT(M),Y(M),WN(0:N,5),S(0:N)
      DIMENSION W1(M1,3),BAN(M),WM(M)
      DIMENSION BIAS(2,0:2),VARK(2,0:2),FAK2(2:4)
C-
C-------- 1. INITIALISATIONS 
      DATA BIAS/.2,.04762,.4286,.1515,1.33,.6293/
      DATA VARK/.6,1.250,2.143,11.93,35.0,381.6/
      DATA FAK2/4.,36.,576./
      NYG=0
      INPUTS=0

      IF(NUE.GT.4.OR.NUE.LT.0) STOP
      IF(NUE.GT.2.AND.ISMO.EQ.0) STOP
      IF(N.LE.2) STOP
      IF(M.LT.1) STOP
      IF(M1.LT.3) STOP
      KK=(KORD-NUE)/2
      IF(2*KK+NUE.NE.KORD) KORD=NUE+2
      IF(KORD.GT.4.AND.ISMO.EQ.0) KORD=NUE+2
      IF(KORD.GT.6.OR.KORD.LE.NUE) KORD=NUE+2
      RVAR=SIG
C-
C-------- 2. COMPUTATION OF S-SEQUENCE 
      S0=1.5*T(1)-0.5*T(2)
      SN=1.5*T(N)-0.5*T(N-1)
      IF(S(N).LE.S(0)) THEN 
        INPUTS=1
        DO 20 I=1,N-1
20        S(I)=.5*(T(I)+T(I+1))
        S(0)=S0
        S(N)=SN
        IF(ISMO.NE.0.AND.IRND.NE.0) GOTO 230
      ELSE
        IF(ISMO.NE.0) GOTO 230
      END IF
C-
C-------- 3. COMPUTATION OF MINIMAL, MAXIMAL ALLOWED GLOBAL BANDWIDTH
      BMAX=(SN-S0)*.5
      BMIN=(SN-S0)/DBLE(N)*DBLE(KORD-1)*.6
C-
C-------- 4. COMPUTE TL,TU
      ITT=0
40    IF (TU.LE.TL) THEN
        TL=.933*S0+.067*SN
        TU=.067*S0+.933*SN
        ITT=ITT+1
      END IF
      TL=MAX(S0,TL)
      TU=MIN(SN,TU)
C-
C-------- 5. COMPUTE INDICES
      IL=1
      IU=N
      WN(1,1)=0.0
      WN(N,1)=0.0
      DO 50 I=1,N
        IF(T(I).LE.TL.OR.T(I).GE.TU) WN(I,1)=0.0
        IF(T(I).GT.TL.AND.T(I).LT.TU) WN(I,1)=1.0
        IF(T(I).LT.TL) IL=I+1
50      IF(T(I).LE.TU) IU=I
      NN=IU-IL+1
      IF(NN.EQ.0.AND.ITT.EQ.0) THEN
        TU=TL-1.0
        GOTO 40
      END IF
      IF(NN.EQ.0.AND.ITT.EQ.1) THEN
        TU=SN
        TL=S0
        GOTO 40
      END IF
C-
C-------- 6. COMPUTE T-GRID FOR INTEGRAL APPROXIMATION 
      DO 60 I=1,M1
         W1(I,2)=1.0
60       W1(I,1)=TL+(TU-TL)*DBLE(I-1)/DBLE(M1-1)
C-
C-------- 7. CALCULATION OF WEIGHT FUNCTION
      ALPHA=1.D0/DBLE(13)
      DO 70 I=IL,IU
        XI=(T(I) - TL)/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 71
70      WN(I,1)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
71    DO 72 I=IU,IL,-1
        XI=(TU-T(I))/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 73
72      WN(I,1)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
73    DO 74 I=1,M1
        XI=(W1(I,1)-TL)/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 75
74      W1(I,2)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
75    DO 76 I=M1,1,-1
        XI=(TU-W1(I,1))/ALPHA/(TU-TL)
        IF(XI.GT.1) GOTO 77
76      W1(I,2)=(10.0-15*XI+6*XI*XI)*XI*XI*XI
77    CONTINUE 
C-
C-------- 8. COMPUTE CONSTANTS FOR ITERATION
      EX=1./DBLE(KORD+KORD+1)
      KK2=(KORD-NUE)
      KK=KK2/2
C-
C-------- 9. ESTIMATING VARIANCE AND SMOOTHED PSEUDORESIDUALS
      IF((SIG.LE..0).AND.(IHOM.EQ.0)) 
     .     CALL RESEST(T(IL),X(IL),NN,WN(IL,2),R2,SIG)
      IF(IHOM.NE.0) THEN 
        CALL RESEST(T,X,N,WN(1,2),SNR,SIG)
        BRES=MAX(BMIN,.2*NN**(-.2)*(S(IU)-S(IL-1)))
        DO 91 I=1,N
          WN(I,3)=T(I)
91        WN(I,2)=WN(I,2)*WN(I,2)
        CALL KERNEL(T,WN(1,2),N,BRES,0,KK2,NYG,S,
     .            WN(1,3),N,WN(1,4))
      ELSE
        CALL COFF(WN(1,4),N,SIG)
      END IF
C-
C-------- 10. ESTIMATE/COMPUTE INTEGRAL CONSTANT
100   VI=0.
      DO 101 I=IL,IU
101      VI=VI+WN(I,1)*N*(S(I)-S(I-1))**2*WN(I,4)       
C- 
C-------- 11. REFINEMENT OF S-SEQUENCE FOR RANDOM DESIGN
      IF(INPUTS.EQ.1.AND.IRND.EQ.0) THEN
        DO 110 I=0,N
          WN(I,5)=DBLE(I)/DBLE(N+1)
          WN(I,2)=(DBLE(I)+.5)/DBLE(N+1)
110       WN(I,3)=WN(I,2)
        EXS=-DBLE(3*KORD+1)/DBLE(6*KORD+3)
        EXSVI=DBLE(KORD)/DBLE(6*KORD+3)
        BS=0.1*(VI/(SN-S0)**2)**EXSVI*N**EXS
        CALL KERNEL(WN(1,5),T,N,BS,0,2,NYG,WN(0,3),WN(0,2),N+1,S(0))
        ISORT=0
        VI=0.0
111     DO 112 I=1,N
        VI=VI+WN(I,1)*N*(S(I)-S(I-1))**2*WN(I,4)       
          IF(S(I).LT.S(I-1)) THEN
            SSI=S(I-1)
            S(I-1)=S(I)
            S(I)=SSI
            ISORT=1
          END IF
112     CONTINUE
        IF(ISORT.EQ.1) GOTO 111
        IF(ISMO.NE.0) GOTO 230
      END IF
      B=BMIN*2.
C-
C-------- 12. COMPUTE INFLATION CONSTANT AND EXPONENT AND LOOP OF ITERATIONS
      CONST=DBLE(2*NUE+1)*FAK2(KORD)*VARK(KK,NUE)*VI
     .       /(DBLE(2*KORD-2*NUE)*BIAS(KK,NUE)**2*DBLE(N))
      FAC=1.1*(1.+(NUE/10.)+0.05*(KORD-NUE-2.))
     .       *N**(2./DBLE((2*KORD+1)*(2*KORD+3)))
      ITENDE=1+2*KORD+KORD*(2*KORD+1)

      DO 120 IT=1,ITENDE
C-
C-------- 13. ESTIMATE DERIVATIVE OF ORDER KORD IN ITERATIONS
        B2=B*FAC
        B2=MAX(B2,BMIN/DBLE(KORD-1)*DBLE(KORD+1))
        B2=MIN(B2,BMAX)
        CALL KERNEL(T,X,N,B2,KORD,KORD+2,NYG,S,W1(1,1),M1,W1(1,3))
C-
C-------- 14. ESTIMATE INTEGRALFUNCTIONAL IN ITERATIONS
        XMY2=.75*(W1(1,2)*W1(1,3)*W1(1,3)+W1(M1,2)*W1(M1,3)*W1(M1,3))
        DO 140 I=2,M1-1
140       XMY2=XMY2+W1(I,2)*W1(I,3)*W1(I,3)
        XMY2=XMY2*(TU-TL)/DBLE(M1)
C-
C-------- 15. FINISH OF ITERATIONS
        B=(CONST/XMY2)**EX
        B=MAX(BMIN,B)
        B=MIN(BMAX,B)

120     CONTINUE
C-------- 16  COMPUTE SMOOTHED FUNCTION WITH GLOBAL PLUG-IN BANDWIDTH
160   CALL KERNEL(T,X,N,B,NUE,KORD,NYG,S,TT,M,Y)
C-
C-------- 17. VARIANCE CHECK 
      IF(IHOM.NE.0) SIG=RVAR
      IF(RVAR.EQ.SIG.OR.R2.LT..88.OR.IHOM.NE.0.OR.NUE.GT.0) GOTO 180
      II=0
      IIL=0
      J=2
      TLL=MAX(TL,TT(1))
      TUU=MIN(TU,TT(M))
      DO 170 I=IL,IU
         IF(T(I).LT.TLL.OR.T(I).GT.TUU) GOTO 170
         II=II+1
         IF(IIL.EQ.0) IIL=I
171       IF(TT(J).LT.T(I)) THEN
            J=J+1
            IF(J.LE.M) GOTO 171
            END IF
       WN(II,3)=X(I)-Y(J)+(Y(J)-Y(J-1))*(TT(J)-T(I))/(TT(J)-TT(J-1))
170    CONTINUE
      IF(IIL.EQ.0.OR.II-IIL.LT.10) THEN
        CALL RESEST(T(IL),WN(1,3),NN,WN(1,4),SNR,RVAR)
       ELSE
        CALL RESEST(T(IIL),WN(1,3),II,WN(1,4),SNR,RVAR)
      END IF
      Q=SIG/RVAR
      IF(Q.LE.2.) CALL COFF(WN(1,4),N,SIG) 
      IF(Q.LE.2.) GOTO 180
      IF(Q.GT.5..AND.R2.GT..95) RVAR=RVAR*.5
      SIG=RVAR
      CALL COFF(WN(1,4),N,SIG)
      GOTO 100
C-
C-------- 18. LOCAL INITIALIZATIONS
180	BVAR=B
 	NUEV=0
	KORDV=2
        NYL=1
C-
C-------- 19. COMPUTE INNER BANDWIDTHS
	G1=0.86*(1.+DBLE(KORD-NUE-2)*.05)*B
	G1=G1*DBLE(N)**(4./DBLE(2*KORD+1)/(2*KORD+5))
	G1=MAX(G1,BMIN/DBLE(KORD-1)*DBLE(KORD+1))
	G1=MIN(G1,BMAX)

	G2=1.4*(1.+DBLE(KORD-NUE-2)*0.05)*B
	G2=G2*DBLE(N)**(2./DBLE(2*KORD+1)/DBLE(2*KORD+3))
	G2=MAX(G2,BMIN)
	G2=MIN(G2,BMAX)
C-
C-------- 20. ESTIMATE/COMPUTE INTEGRAL CONSTANT VI LOCALLY
        DO 200 I=1,N
200	   WN(I,4)=DBLE(N)*WN(I,4)*(S(I)-S(I-1))
	DO 201 J=1,M
	   BAN(J)=BVAR
             WM(j)=tt(j)
	   IF(TT(J).LT.S(0)+G1) THEN
             DIST=((TT(J)-G1-S(0))/G1)**2
             BAN(J)=BVAR*(1.0+1.0*DIST)
             BAN(J)=MIN(BAN(J),BMAX)
             WM(J)=TT(J)+.5*DIST*G1
           ELSE IF(TT(J).GT.S(N)-G1) THEN
	     DIST=((TT(J)-S(N)+G1)/G1)**2
             BAN(J)=BVAR*(1.0+1.0*DIST)
             BAN(J)=MIN(BAN(J),BMAX)
             WM(J)=TT(J)-.5*DIST*G1
           END IF
201	CONTINUE
	CALL KERNEL(T,WN(1,4),N,BVAR,NUEV,KORDV,NYL,S,WM,M,BAN)
C-
C-------- 21. ESTIMATION OF KORDTH DERIVATIVE LOCALLY
        WSTEP=(TT(M)-TT(1))/DBLE(M1-2)
	DO 210 J=2,M1
           W1(J,2)=TT(1)+DBLE(J-2)*WSTEP
210	   W1(J,1)=TT(1)+DBLE(J-1.5)*WSTEP
        W1(1,1)=TT(1)+.5*WSTEP

	CALL KERNEL(T,X,N,G1,KORD,KORD+2,NYG,S,W1(2,2),M1-1,W1(2,3))

	DO 211 J=2,M1
211        W1(J,3)=W1(J,3)*W1(J,3)
	DO 212 J=1,M
	   Y(J)=G2
	   IF(TT(J).LT.S(0)+G1) THEN           
             Y(J)=G2*(1.0+1.0*((TT(J)-G1-S(0))/G1)**2)
             Y(J)=MIN(Y(J),BMAX)
           ELSE IF(TT(J).GT.S(N)-G1) THEN
             Y(J)=G2*(1.0+1.0*((TT(J)-S(N)+G1)/G1)**2)
             Y(J)=MIN(Y(J),BMAX)
           END IF
212	CONTINUE

	CALL KERNP(W1(2,2),W1(2,3),M1-1,G2,NUEV,KORDV,NYL,
     .	            W1(1,1),WM,M,Y)
C-
C-------- 22. FINISH
        IRND=1-IRND
	DO 220 J=1,M
          XH=BMIN**(2*KORD+1)*ABS(Y(J))*VI/CONST
          XXH=CONST*ABS(BAN(J))/VI/BMAX**(2*KORD+1)
          IF(BAN(J).LT.XH) THEN
     	     BAN(J)=BMIN
 	  ELSE
             IF(Y(J).LT.XXH) THEN
 	        BAN(J)=BMAX
             ELSE
                BAN(J)=(CONST*BAN(J)/Y(J)/VI)**EX
	     END IF
          END IF
220	CONTINUE
C-
C-------- 23. COMPUTE SMOOTHED FUNCTION WITH LOCAL PLUG-IN BANDWIDTH
230     DO 231 J=1,M
231	   Y(J)=BAN(J)
	CALL KERNEL(T,X,N,B,NUE,KORD,NYL,S,TT,M,Y)

	RETURN
	END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC     SEVERAL KERNEL SMOOTHING SUBROUTINES WHICH ARE USED BY GLKERN.F 
CC     AND LOKERN.F, VERSION OCT 1996
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC     THIS FILE CONTAINS:
CC
CC     SUBROUTINE RESEST(T,X,N,RES,SNR,SIGMA2)
CC                FOR VARIANCE ESTIMATION
CC
CC     SUBROUTINE KERNEL(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
CC                DRIVER SUBROUTINE FOR KERNEL REGRESSION ESTIMATION
CC                CALLS FAST OR CONVENTIAL KERNEL ROUTINE
CC
CC     SUBROUTINE KERNP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
CC                DRIVER SUBROUTINE FOR KERNEL REGRESSION ESTIMATION
CC                WITHOUT USE OF BOUNDARY KERNELS
CC---------------------------------------------------------------------
CC     SUBROUTINE KERNFA(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
CC                FAST ALGORITHM FOR KERNEL ESTIMATION
CC     SUBROUTINE DREG(SW,A1,A2,IORD,X,SL,SR,T,B,IFLOP)
CC                USED BY SUBROUTINE KERNFA,KERNFP
CC     SUBROUTINE LREG(SW,A3,IORD,D,DOLD,Q,C)
CC                USED BY SUBROUTINE KERNFA,KERNFP
CC     SUBROUTINE FREG(SW,NUE,KORD,IBOUN,Y,C,ICALL,A)
CC                USED BY SUBROUTINE KERNFA,KERNFP
CC     SUBROUTINE KERNFP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
CC                FAST ALGORITHM FOR KERNP ESTIMATION WITHOUT BOUNDARY
CC---------------------------------------------------------------------
CC     SUBROUTINE KERNCL(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
CC                CONVENTIONAL ALGORITHM FOR KERNEL ESTIMATION
CC     SUBROUTINE SMO(S,X,N,TAU,WID,NUE,IORD,IBOUN,IST,S1,C,Y)
CC                SINGLE ESTIMATION STEP, USED BY KERNCL
CC     SUBROUTINE KERNCP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
CC                CONVENTIONAL ALGORITHM WITHOUT BOUNDARY KERNELS
CC     SUBROUTINE SMOP(S,X,N,TAU,WID,NUE,IORD,IBOUN,IST,S1,C,Y)
CC                SINGLE ESTIMATION STEP, USED BY KERNCP
CC---------------------------------------------------------------------
CC     SUBROUTINE COFFI(NUE,KORD,C)
CC                KERNEL COEFFICIENT OF POLYNOMIAL KERNELS USED BY
CC                KERNCL,KERNCP AND KERNFP
CC     SUBROUTINE COFFB(NUE,KORD,Q,IBOUN,C)
CC                KERNEL COEFFICIENT OF POLYNOMIAL BOUNDARY KERNELS
CC                USED BY KERNFA AND KERNCL
CC---------------------------------------------------------------------
CC     SUBROUTINE COFF(X,N,FA)
CC                SIMPLE SUBROUTINE FOR ARRAY INITIALIZATION
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE RESEST(T,X,N,RES,SNR,SIGMA2)
C-----------------------------------------------------------------------
C       VERSION: JUNE, 1996
C
C       PURPOSE:
C
C       COMPUTES ONE-LEAVE-OUT RESIDUALS FOR NONPARAMETRIC ESTIMATION
C       OF RESIDUAL VARIANCE (LOCAL LINEAR APPROXIMATION FOLLOWED BY
C       REWEIGHTING) 
C
C     PARAMETERS:
C
C     INPUT   T(N)      ABSCISSAE (ORDERED: T(I)<=T(I+1))
C     INPUT   X(N)      DATA
C     INPUT   N         LENGTH OF DATA ( >2 )
C     OUTPUT  RES(N)    RESIDUALS AT T(1),...,T(N)
C     OUTPUT  SNR       EXPLAINED VARIANCE OF THE TRUE CURVE
C     OUTPUT  SIGMA2    ESTIMATION OF SIGMA**2 (RESIDUAL VARIANCE)
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(N),T(N),RES(N)
C-
      SIGMA2=0.
      EX=X(1)*(T(2)-T(1))
      EX2=X(1)*EX
      DO 1 I=2,N-1
         TT=T(I+1)-T(I-1)
         IF(TT.NE.0.) G1=(T(I+1)-T(I))/TT
         IF(TT.EQ.0.) G1=.5
         G2=1.-G1
         RES(I)=(X(I)-G1*X(I-1)-G2*X(I+1))/SQRT(1.+G1*G1+G2*G2)
         SIGMA2=SIGMA2+RES(I)*RES(I)
         SX=X(I)*TT
         EX=EX+SX
         EX2=EX2+X(I)*SX
1        CONTINUE
      TT=T(3)-T(2)
      IF(TT.NE.0.) G1=(T(1)-T(2))/TT
      IF(TT.EQ.0.) G1=.5
      G2=1.-G1
      RES(1)=(X(1)-G1*X(3)-G2*X(2))/SQRT(1.+G1*G1+G2*G2)
      TT=T(N-1)-T(N-2)
      IF(TT.NE.0.) G1=(T(N-1)-T(N))/TT
      IF(TT.EQ.0.) G1=.5
      G2=1.-G1
      RES(N)=(X(N)-G1*X(N-2)-G2*X(N-1))/SQRT(1.+G1*G1+G2*G2)
      SIGMA2=(SIGMA2+RES(1)*RES(1)+RES(N)*RES(N))/N
C-
      SX=X(N)*(T(N)-T(N-1))
      DN=2.*(T(N)-T(1))
      EX=(EX+SX)/DN
      EX2=(EX2+X(N)*SX)/DN
      IF(EX2.EQ.0) SNR=0.
      IF(EX2.GT.0) SNR=1-SIGMA2/(EX2-EX*EX)
      RETURN
      END

      SUBROUTINE KERNEL(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION MAY, 1995
C
C       DRIVER SUBROUTINE FOR KERNEL SMOOTHING, CHOOSES BETWEEN
C       STANDARD AND O(N) ALGORITHM 
C
C  PARAMETERS :
C
C  INPUT    T(N)         INPUT GRID (REGRESSION DESIGN)
C  INPUT    X(N)         DATA, GIVEN ON T(N)
C  INPUT    N            LENGTH OF X
C  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    KORD         ORDER OF KERNEL (<=6); DEFAULT IS KORD=NUE+2
C  INPUT    NY           0: GLOBAL BANDWIDTH (DEFAULT)
C                        1: VARIABLE BANDWIDTHS, GIVEN IN Y AS INPUT
C  INPUT    S(0:N)       INTERPOLATION SEQUENCE
C  INPUT    TT(M)        OUTPUT GRID. MUST BE PART OF INPUT GRID FOR
C                        IEQ=0
C  INPUT    M            NUMBER OF POINTS WHERE FUNCTION IS ESTIMATED,
C                         OR  LENGTH OF TT. DEFAULT IS M=400
C  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
C  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION T(N),X(N),B,S(0:N),TT(M),Y(M),CHAN
      INTEGER N,NUE,KORD,NY,M
C-
C------  COMPUTING CHANGE POINT
      CHAN=(5.+KORD)*MAX(1.,SQRT(FLOAT(N)/FLOAT(M)))
C------
      IF(B*(N-1)/(T(N)-T(1)).LT.CHAN) THEN
            CALL KERNCL(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
         ELSE
            CALL KERNFA(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
       END IF
C-
      RETURN
      END


      SUBROUTINE KERNP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION JANUARY, 1997
C
C       DRIVER SUBROUTINE FOR KERNEL SMOOTHING, CHOOSES BETWEEN
C       STANDARD AND O(N) ALGORITHM WITHOUT USING BOUNDARY KERNELS 
C
C  PARAMETERS :
C
C  INPUT    T(N)         INPUT GRID (REGRESSION DESIGN)
C  INPUT    X(N)         DATA, GIVEN ON T(N)
C  INPUT    N            LENGTH OF X
C  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    KORD         ORDER OF KERNEL (<=6); DEFAULT IS KORD=NUE+2
C  INPUT    NY           0: GLOBAL BANDWIDTH (DEFAULT)
C                        1: VARIABLE BANDWIDTHS, GIVEN IN Y AS INPUT
C  INPUT    S(0:N)       INTERPOLATION SEQUENCE
C  INPUT    TT(M)        OUTPUT GRID. MUST BE PART OF INPUT GRID FOR
C                        IEQ=0
C  INPUT    M            NUMBER OF POINTS WHERE FUNCTION IS ESTIMATED,
C                         OR  LENGTH OF TT. DEFAULT IS M=400
C  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
C  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION T(N),X(N),B,S(0:N),TT(M),Y(M),CHAN
      INTEGER N,NUE,KORD,NY,M
C-
C------  COMPUTING CHANGE POINT
      CHAN=(5.+KORD)*MAX(1.,SQRT(FLOAT(N)/FLOAT(M)))
C------
      IF(B*(N-1)/(T(N)-T(1)).LT.CHAN) THEN
            CALL KERNCP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
         ELSE
            CALL KERNFP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
       END IF
C-
      RETURN
      END

      SUBROUTINE KERNFA(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION: MAY, 1995
C
C       PURPOSE:
C
C       COMPUTATION OF KERNEL ESTIMATE USING O(N) ALGORITHM BASED ON
C       LEGENDRE POLYNOMIALS, GENERAL SPACED DESIGN AND LOCAL 
C       BANDWIDTH ALLOWED. (NEW INITIALISATIONS OF THE LEGENDRE SUMS
C       FOR NUMERICAL REASONS)
C
C  PARAMETERS :
C
C  INPUT    T(N)         INPUT GRID
C  INPUT    X(N)         DATA, GIVEN ON T(N)
C  INPUT    N            LENGTH OF X
C  INPUT    B            ONE SIDED BANDWIDTH
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    KORD         ORDER OF KERNEL (<=6)
C  INPUT    NY           0, GLOBAL BANDWIDTH; 1, LOCAL BANDWIDTH IN Y
C  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
C  INPUT    TT(M)        OUTPUT GRID
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE
C  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
C  OUTPUT   Y(M)         ESTIMATED FUNCTION
C
C
C-----------------------------------------------------------------------
      INTEGER N,NUE,KORD,NY,M
      INTEGER J,K,IORD,INIT,ICALL,I,IBOUN
      INTEGER JL,JR,JNR,JNL
      DOUBLE PRECISION X(N),T(N),S(0:N),TT(M),Y(M),B
      DOUBLE PRECISION C(7),SW(7),XF(7),DOLD
      DOUBLE PRECISION A(7,7),A1(7),A2(7),A3(7,7),CM(7,6)
      DOUBLE PRECISION BMIN,BMAX,BB,WWL,WWR,WID,WR,WIDO
C-
C------ COMPUTE CONSTANTS FOR LATER USE
      S0=1.5*T(1)-0.5*T(2)
      SN=1.5*T(N)-0.5*T(N-1)
      BMIN=(SN-S0)*.6D0/DBLE(N)*DBLE(KORD-1)
      BMAX=(S(N)-S(0))*.5
      IF(KORD.EQ.2) BMIN=BMIN*0.1D0
      IORD=KORD+1
      DO 2 K=3,IORD
        A1(K)=DBLE(2*K-1)/DBLE(K)
2       A2(K)=DBLE(1-K)/DBLE(K)
C-
      INIT=0
      ICALL=0
      DOLD=0.D0
C-
C------ SMOOTHING LOOP
      DO 100 I=1,M
        BB=B
        IF (NY .EQ. 1) BB=Y(I)
        IF(BB.LT.BMIN) BB=BMIN
        IF(BB.GT.BMAX) BB=BMAX
        IBOUN=0
C-
C------ COMPUTE LEFT BOUNDARY KERNEL
        IF(TT(I).LT.S(0)+BB) THEN
          WWL=S(0)
          WWR=S(0)+BB+BB
          WID=WWR-TT(I)
          IBOUN=1
          CALL COFFB(NUE,KORD,(TT(I)-S(0))/WID,IBOUN,C)
        END IF
C-
C------ COMPUTE RIGHT BOUNDARY KERNEL
        IF(TT(I)+BB.GT.S(N)) THEN
          WWL=S(N)-(BB+BB)
          WWR=S(N)
          WID=TT(I)-WWL
          IBOUN=-1
          CALL COFFB(NUE,KORD,(S(N)-TT(I))/WID,IBOUN,C)
        END IF
C-
C------ NO BOUNDARY
        IF(IBOUN.EQ.0) THEN
          WID=BB
          WWL=TT(I)-BB
          WWR=TT(I)+BB
        END IF
C-
C------ INITIALISATION FOR INIT=0
        IF(INIT.EQ.0) THEN
          DO 44 K=1,IORD
44          SW(K)=0.
          JL=1
          DO 48 J=1,N
            IF(S(J-1).LT.WWL) THEN
              JL=J+1
            ELSE
              IF(S(J).GT.WWR) GOTO 488
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
            END IF
48          CONTINUE
488       JR=J-1
          WR=WWR
          INIT=1
          GOTO 6666
        ELSE
          INIT=INIT+1
        END IF
C-
C------ COMPARE OLD SUM WITH NEW SMOOTHING INTERVALL TT(I)-B,TT(I)+B
        IF(S(JR-1).GE.WWL) THEN
          JNR=JR
          JNL=JL
          IF(S(JR).GT.WWR) THEN
            DO 201 J=JR,JL,-1
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I-1),WIDO,-1)
              JNR=J-1
              IF(S(JNR).LE.WWR) GOTO 2011
201           CONTINUE
2011      END IF
          IF(S(JL-1).LT.WWL) THEN
            DO 301 J=JL,JR      
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I-1),WIDO,-1)
              JNL=J+1
              IF(S(J).GE.WWL) GOTO 3011
301           CONTINUE
3011      END IF
C-
C------ UPDATING OF SW
          CALL LREG(SW,A3,IORD,(TT(I)-TT(I-1))/WID,DOLD,WIDO/WID,CM)
          IF(JNR.EQ.JR) THEN
            DO 401 J=JR+1,N
              IF(S(J).GT.WWR) GOTO 4011
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
              JNR=J
401           CONTINUE
4011      END IF
          JR=JNR
          IF(JL.EQ.JNL) THEN
            DO 402  J=JL-1,1,-1
              IF(S(J-1).LT.WWL) GOTO 4022
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
              JNL=J
402           CONTINUE
4022      END IF
          JL=JNL
        ELSE
C-
C------ NEW INITIALISATION OF SW
          DO 22 K=1,IORD
22          SW(K)=0.
          DO 202 J=JR,N
            IF(S(J-1).LT.WWL) THEN
              JL=J+1
            ELSE
              IF(S(J).GT.WWR) GOTO 2022
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
            END IF
202         CONTINUE
2022      JR=J-1
          WR=WWR
        END IF
6666    CONTINUE
C-
C------ IF BANDWIDTH IS TOO SMALL NO SMOOTHING
        IF(WWL.GE.S(JR-1).AND.WWR.LE.S(JR)) THEN
          Y(I)=X(JR)
          IF(NUE.GT.0) Y(I)=0.D0
        ELSE
C-
C------ ADD FIRST AND LAST POINT OF THE SMOOTHING INTERVAL
          DO 501 K=1,IORD
501         XF(K)=SW(K)
          IF(JL.NE.1)
     .       CALL DREG(XF,A1,A2,IORD,X(JL-1),WWL,S(JL-1),TT(I),WID,1)
          IF(JR.NE.N)
     .       CALL DREG(XF,A1,A2,IORD,X(JR+1),S(JR),WWR,TT(I),WID,1)
C-
C------ NOW THE SUMS ARE BUILT THAT ARE NEEDED TO COMPUTE THE ESTIMATE
          CALL FREG(XF,NUE,KORD,IBOUN,Y(I),C,ICALL,A)
          IF(NUE.GT.0) Y(I)=Y(I)/(WID**NUE)
        END IF
C-
C------ NEW INITIALISATION ?
        IF(JL.GT.JR.OR.WWL.GT.WR.OR.INIT.GT.100) INIT=0
        WIDO=WID
C-
100     CONTINUE
C-
      RETURN
      END

      SUBROUTINE KERNFP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION: JANUARY, 1997
C
C       PURPOSE:
C
C       COMPUTATION OF KERNEL ESTIMATE USING O(N) ALGORITHM BASED ON
C       LEGENDRE POLYNOMIALS, GENERAL SPACED DESIGN AND LOCAL 
C       BANDWIDTH ALLOWED. (NEW INITIALISATIONS OF THE LEGENDRE SUMS
C       FOR NUMERICAL REASONS) WITHOUT BOUNDARY KERNELS, JUST NORMALIZING
C
C  PARAMETERS :
C
C  INPUT    T(N)         INPUT GRID
C  INPUT    X(N)         DATA, GIVEN ON T(N)
C  INPUT    N            LENGTH OF X
C  INPUT    B            ONE SIDED BANDWIDTH
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    KORD         ORDER OF KERNEL (<=6)
C  INPUT    NY           0, GLOBAL BANDWIDTH; 1, LOCAL BANDWIDTH IN Y
C  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
C  INPUT    TT(M)        OUTPUT GRID
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE
C  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
C  OUTPUT   Y(M)         ESTIMATED FUNCTION
C
C
C-----------------------------------------------------------------------
      INTEGER N,NUE,KORD,NY,M
      INTEGER J,K,IORD,INIT,ICALL,I,IBOUN
      INTEGER JL,JR,JNR,JNL
      DOUBLE PRECISION X(N),T(N),S(0:N),TT(M),Y(M),B
      DOUBLE PRECISION C(7),SW(7),XF(7),DOLD,QQ,Q,XNOR
      DOUBLE PRECISION A(7,7),A1(7),A2(7),A3(7,7),CM(7,6)
      DOUBLE PRECISION BMIN,BMAX,BB,WWL,WWR,WID,WR,WIDO
C-
C------ COMPUTE CONSTANTS FOR LATER USE
      S0=1.5*T(1)-0.5*T(2)
      SN=1.5*T(N)-0.5*T(N-1)
      BMIN=(SN-S0)*.6D0/DBLE(N)*DBLE(KORD-1)
      BMAX=(S(N)-S(0))*.5
      IF(KORD.EQ.2) BMIN=BMIN*0.1D0
      IORD=KORD+1
      CALL COFFI(NUE,KORD,C)
      DO 2 K=3,IORD
        A1(K)=DBLE(2*K-1)/DBLE(K)
2       A2(K)=DBLE(1-K)/DBLE(K)
C-
      INIT=0
      ICALL=0
      DOLD=0.D0
C-
C------ SMOOTHING LOOP
      DO 100 I=1,M
        BB=B
        IF (NY .EQ. 1) BB=Y(I)
        IF(BB.LT.BMIN) BB=BMIN
        IF(BB.GT.BMAX) BB=BMAX
        IBOUN=0
C-
C------ COMPUTE LEFT BOUNDARY
        IF(TT(I).LT.S(0)+BB) THEN
          WWL=S(0)
          WWR=S(0)+BB+BB
          WID=WWR-TT(I)
          IBOUN=1
        END IF
C-
C------ COMPUTE RIGHT BOUNDARY
        IF(TT(I)+BB.GT.S(N)) THEN
          WWL=S(N)-(BB+BB)
          WWR=S(N)
          WID=TT(I)-WWL
          IBOUN=-1
        END IF
C-
C------ NO BOUNDARY
        IF(IBOUN.EQ.0) THEN
          WID=BB
          WWL=TT(I)-BB
          WWR=TT(I)+BB
          XNOR=1.D0
        END IF
C-
C------ COMPUTE NORMALIZING CONSTANT 
        IF(IBOUN.NE.0) THEN
          IF(IBOUN.EQ.1) Q=(TT(I)-S(0))/WID
          IF(IBOUN.EQ.-1) Q=(S(N)-TT(I))/WID
          QQ=Q*Q
          XNOR=C(1)*(1.D0+Q)
          DO 3 K=3,IORD,2
             Q=Q*QQ
3            XNOR=XNOR+C(K)*(1.D0+Q)
          IBOUN=0
        END IF 
C-
C------ INITIALISATION FOR INIT=0
        IF(INIT.EQ.0) THEN
          DO 44 K=1,IORD
44          SW(K)=0.
          JL=1
          DO 48 J=1,N
            IF(S(J-1).LT.WWL) THEN
              JL=J+1
            ELSE
              IF(S(J).GT.WWR) GOTO 488
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
            END IF
48          CONTINUE
488       JR=J-1
          WR=WWR
          INIT=1
          GOTO 6666
        ELSE
          INIT=INIT+1
        END IF
C-
C------ COMPARE OLD SUM WITH NEW SMOOTHING INTERVALL TT(I)-B,TT(I)+B
        IF(S(JR-1).GE.WWL) THEN
          JNR=JR
          JNL=JL
          IF(S(JR).GT.WWR) THEN
            DO 201 J=JR,JL,-1
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I-1),WIDO,-1)
              JNR=J-1
              IF(S(JNR).LE.WWR) GOTO 2011
201           CONTINUE
2011      END IF
          IF(S(JL-1).LT.WWL) THEN
            DO 301 J=JL,JR      
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I-1),WIDO,-1)
              JNL=J+1
              IF(S(J).GE.WWL) GOTO 3011
301           CONTINUE
3011      END IF
C-
C------ UPDATING OF SW
          CALL LREG(SW,A3,IORD,(TT(I)-TT(I-1))/WID,DOLD,WIDO/WID,CM)
          IF(JNR.EQ.JR) THEN
            DO 401 J=JR+1,N
              IF(S(J).GT.WWR) GOTO 4011
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
              JNR=J
401           CONTINUE
4011      END IF
          JR=JNR
          IF(JL.EQ.JNL) THEN
            DO 402  J=JL-1,1,-1
              IF(S(J-1).LT.WWL) GOTO 4022
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
              JNL=J
402           CONTINUE
4022      END IF
          JL=JNL
        ELSE
C-
C------ NEW INITIALISATION OF SW
          DO 22 K=1,IORD
22          SW(K)=0.
          DO 202 J=JR,N
            IF(S(J-1).LT.WWL) THEN
              JL=J+1
            ELSE
              IF(S(J).GT.WWR) GOTO 2022
              CALL DREG(SW,A1,A2,IORD,X(J),S(J-1),S(J),TT(I),WID,1)
            END IF
202         CONTINUE
2022      JR=J-1
          WR=WWR
        END IF
6666    CONTINUE
C-
C------ IF BANDWIDTH IS TOO SMALL NO SMOOTHING
        IF(WWL.GE.S(JR-1).AND.WWR.LE.S(JR)) THEN
          Y(I)=X(JR)
          IF(NUE.GT.0) Y(I)=0.D0
        ELSE
C-
C------ ADD FIRST AND LAST POINT OF THE SMOOTHING INTERVAL
          DO 501 K=1,IORD
501         XF(K)=SW(K)
          IF(JL.NE.1)
     .       CALL DREG(XF,A1,A2,IORD,X(JL-1),WWL,S(JL-1),TT(I),WID,1)
          IF(JR.NE.N)
     .       CALL DREG(XF,A1,A2,IORD,X(JR+1),S(JR),WWR,TT(I),WID,1)
C-
C------ NOW THE SUMS ARE BUILT THAT ARE NEEDED TO COMPUTE THE ESTIMATE
          CALL FREG(XF,NUE,KORD,IBOUN,Y(I),C,ICALL,A)
          IF(NUE.GT.0) Y(I)=Y(I)/(WID**NUE)
          IF(XNOR.NE.1.D0) Y(I)=Y(I)/XNOR
        END IF
C-
C------ NEW INITIALISATION ?
        IF(JL.GT.JR.OR.WWL.GT.WR.OR.INIT.GT.100) INIT=0
        WIDO=WID
C-
100     CONTINUE
C-
      RETURN
      END

      SUBROUTINE DREG(SW,A1,A2,IORD,X,SL,SR,T,B,IFLOP)
C-----------------------------------------------------------------------
C       VERSION: MAY, 1995
C
C       PURPOSE:
C
C       COMPUTES NEW LEGENDRE SUMS (FOR REGRESSION)
C
C       PARAMETERS:
C                    **************   INPUT   *******************
C
C        SW(IORD)  :   OLD SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.
C        A1(7)     :   CONSTANTS OF RECURSIVE FORMULA FOR LEGENDRE POL.
C        A2(7)     :                             "
C        IORD      :   ORDER OF KERNEL POLYNOMIAL
C        X         :   DATA POINT
C        SL        :   LEFT S-VALUE
C        SR        :   RIGHT S-VALUE
C        T         :   POINT WHERE THE SMOOTHED VALUE IS TO BE ESTIMATED
C        B         :   BANDWIDTH
C        IFLOP     :   1: ADDITION, ELSE SUBTRACTION
C
C
C                    **************   OUTPUT   *******************
C
C        SW(IORD)  :   NEW SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.
C
C-----------------------------------------------------------------------
      INTEGER IORD,IFLOP,K
      DOUBLE PRECISION SW(7),A1(7),A2(7),X,SL,SR,T,B
      DOUBLE PRECISION P(7,2)
C-
C------  COMPUTE LEGENDRE POLYNOMIALS
      P(1,1)=(T-SL)/B
      P(1,2)=(T-SR)/B
      P(2,1)=1.5D0*P(1,1)*P(1,1)-.5D0
      P(2,2)=1.5D0*P(1,2)*P(1,2)-.5D0
      DO 1 K=3,IORD
        P(K,1)=A1(K)*P(K-1,1)*P(1,1)+A2(K)*P(K-2,1)
1       P(K,2)=A1(K)*P(K-1,2)*P(1,2)+A2(K)*P(K-2,2)
C-
C------  COMPUTE NEW LEGENDRE SUMS
      IF(IFLOP.EQ.1) THEN
        DO 2 K=1,IORD
2         SW(K)=SW(K)+(P(K,1)-P(K,2))*X
      ELSE
        DO 3 K=1,IORD
3         SW(K)=SW(K)+(P(K,2)-P(K,1))*X
      END IF
      RETURN
      END

      SUBROUTINE LREG(SW,A3,IORD,D,DOLD,Q,C)
C------------------------------------------------------------------
C       VERSION: MAY, 1995
C
C       PURPOSE:
C
C       UPDATE OF SW-SEQUENCE ACCORDING TO NEW BANDWIDTH AND NEW DATA
C       (VERSION FOR REGRESSION)
C
C       PARAMETERS :
C                          **************   INPUT   *******************
C
C               SW(IORD)  :  SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.
C               IORD      :  ORDER OF KERNEL POLYNOMIAL
C               D         :  DIST. TO THE NEXT POINT DIVIDED BY BANDW.
C               DOLD      :  D       PREVIOUS STEP
C               Q         :  NEW BANDWIDTH DIVIDED BY OLD BANDWIDTH
C                          **************    WORK   *******************
C
C               A3(7,7)   :  MATRIX (P*Q*P)**(-1)        
C               C(7,6)    :  MATRIX OF COEFFICIENTS
C                          **************   OUTPUT   ******************
C
C               SW(IORD)  :  UPDATED VERSION OF SW
C
C---------------------------------------------------------------------
      INTEGER IORD,K,I,L
      DOUBLE PRECISION D,DOLD,Q,DD,WW,QQ,XX
      DOUBLE PRECISION A3(7,7),C(7,6),SW(7)
C-
C- BUILD UP MATRIX
      IF(DOLD.NE.D.OR.DOLD.EQ.0) THEN
        DOLD=D
        DD=D*D
C-
        IF(IORD.EQ.7) THEN
          C(7,6)=13.D0*D
          C(7,5)=71.5D0*DD
          C(7,4)=(214.5D0*DD+9.D0)*D
          C(7,3)=(375.375D0*DD+77.D0)*DD
          C(7,2)=((375.375D0*DD+247.5D0)*DD+5.D0)*D
          C(7,1)=((187.6875D0*DD+346.5D0)*DD+40.5D0)*DD
        END IF
C-
        IF(IORD.GE.6) THEN
          C(6,5)=11.D0*D
          C(6,4)=49.5D0*DD
          C(6,3)=(115.5D0*DD+7.D0)*D
          C(6,2)=(144.375D0*DD+45.D0)*DD
          C(6,1)=((86.625D0*DD+94.5D0)*DD+3.D0)*D
        END IF
C-
        IF(IORD.GE.5) THEN
          C(5,4)=9.D0*D
          C(5,3)=31.5D0*DD
          C(5,2)=(52.5D0*DD+5.D0)*D
          C(5,1)=(39.375D0*DD+21.D0)*DD
        END IF
C-
        IF(IORD.GE.4) THEN
          C(4,3)=7.D0*D
          C(4,2)=17.5D0*DD
          C(4,1)=(17.5D0*DD+3.D0)*D
        END IF
C-
        IF(IORD.GE.3) THEN
          C(3,2)=5.D0*D
          C(3,1)=7.5D0*DD
        END IF
C-
        C(2,1)=3.D0*D
      END IF
      IF(Q.LT..9999.OR.Q.GT.1.0001) THEN
C-
C------- BUILT UP MATRIX A3=P*Q*P**-1
        A3(1,1)=Q
        DO 1 K=2,IORD
1         A3(K,K)=A3(K-1,K-1)*Q
        WW=Q*Q-1.D0
        DO 2 K=1,IORD-2
          WW=WW*Q
2         A3(K+2,K)=(K+.5D0)*WW
C-
        IF(IORD.GE.5) THEN
          QQ=A3(2,2)
          A3(5,1)=Q*(1.875D0+QQ*(-5.25D0+QQ*3.375D0))
        END IF
        IF(IORD.GE.6) A3(6,2)=QQ*(4.375D0+QQ*(-11.25D0+QQ*6.875D0))
        IF(IORD.EQ.7) THEN
          A3(7,1)=Q*(-2.1875D0+QQ*(11.8125D0+QQ*(-18.5625D0+QQ*
     .          8.9375D0)))
          A3(7,3)=Q*QQ*(7.875D0+QQ*(-19.25D0+QQ*11.375D0))
        END IF
C-
C------- COMPUTE A*C AND NEW LEGENDRE SUMS
        DO 10 I=IORD,2,-1
          XX=0.
          DO 20 K=1,I
            WW=0.
            DO 30 L=K,I-1,2
30            WW=WW+A3(L,K)*C(I,L)
            IF(MOD(I-K,2).EQ.0) WW=WW+A3(I,K)
            XX=XX+WW*SW(K)
20          CONTINUE
          SW(I)=XX
10        CONTINUE
        SW(1)=A3(1,1)*SW(1)
      ELSE
        DO 111 I=IORD,2,-1
          DO 112 K=1,I-1
112         SW(I)=SW(I)+C(I,K)*SW(K)
111       CONTINUE
      END IF
      RETURN
      END

      SUBROUTINE FREG(SW,NUE,KORD,IBOUN,Y,C,ICALL,A)
C------------------------------------------------------------------
C       SHORT-VERSION: MAY, 1995
C
C       PURPOSE:
C
C       FINAL COMPUTATION OF A SMOOTHED VALUE VIA LEGENDRE POLYNOMIALS
C
C       PARAMETERS :
C                          **************   INPUT   *******************
C
C               SW(KORD+1):  SUM OF DATA WEIGHTS FOR LEGENDRE POLYNOM.
C               NUE       :  ORDER OF DERIVATIVE
C               KORD      :  ORDER OF KERNEL
C               IBOUN     :  0: INTERIOR KERNEL, ELSE BOUNDARY KERNEL
C               C(KORD+1) :  SEQUENCE OF POLYN. COEFF. FOR BOUND. KERNEL
C               ICALL     :  PARAMETER USED TO INITIALISE COMPUTATION
C                         :   OF A MATRIX
C                          **************    WORK    ******************
C
C               A(7,7)    :  MATRIX OF COEFFICIENTS
C                          **************   OUTPUT   ******************
C
C               Y          :  COMPUTED ESTIMATE
C
C--------------------------------------------------------------------
      INTEGER NUE,KORD,IBOUN,ICALL,I,J
      DOUBLE PRECISION SW(7),C(7),A(7,7),Y,WW
C-
C------- DEFINITION OF LEGENDRE COEFFICIENTS FOR BOUNDARY
      IF(ICALL.EQ.0.AND.IBOUN.NE.0) THEN
             A(2,2)=2./3.
         A(1,3)=.6
                 A(3,3)=.4
             A(2,4)=4./7.
                     A(4,4)=8./35.
         A(1,5)=27./63.
                 A(3,5)=28./63.
                         A(5,5)=8./63.
             A(2,6)=110./231.
                     A(4,6)=72./231.
                             A(6,6)=16./231.
         A(1,7)=143./429.
                 A(3,7)=182./429.
                         A(5,7)=88./429.
                                 A(7,7)=16./429.
         ICALL=1
      END IF
      IF(IBOUN.NE.0) THEN
C-
C------- COMPUTATION OF THE SMOOTHED VALUE AT BOUNDARY
        Y=C(1)*SW(1)+C(2)*A(2,2)*SW(2)
        DO 1 J=3,KORD+1
          WW=A(J,J)*SW(J)
          DO 2 I=J-2,1,-2
2           WW=WW+A(I,J)*SW(I)
          Y=Y+C(J)*WW
1         CONTINUE
      ELSE
C-
C------- COMPUTATION OF THE SMOOTHED VALUE AT INTERIOR
        IF(NUE.EQ.0) THEN
          IF(KORD.EQ.2) Y=-.1*SW(3)+.6*SW(1)
          IF(KORD.EQ.4) Y=(SW(5)-4.*SW(3)+9.*SW(1))/12.
          IF(KORD.EQ.6) Y=-7.2115379E-02*SW(7)+.25961537*SW(5)
     .                    -.4375*SW(3)+.75*SW(1)
        END IF
        IF(NUE.EQ.1) THEN
          IF(KORD.EQ.3) Y=(3.*SW(4)-10.*SW(2))/14.
          IF(KORD.EQ.5) Y=(-15.*SW(6)+48.*SW(4)-55.*SW(2))/44.
        END IF
        IF(NUE.EQ.2) THEN
          IF(KORD.EQ.4) Y=(-5.*SW(5)+14.*SW(3)-9.*SW(1))/6.
          IF(KORD.EQ.6) Y=2.01923*SW(7)-5.76923*SW(5)+5.25*SW(3)
     .                   -1.5*SW(1)
        END IF
        IF(NUE.EQ.3) Y=4.772727*SW(6)-12.272727*SW(4)+7.5*SW(2)
        IF(NUE.EQ.4) Y=-36.34615*SW(7)+88.84615*SW(5)-52.5*SW(3)
      END IF
C-
      RETURN
      END

      SUBROUTINE KERNCL(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION JANUARY 1995
C
C       KERNEL SMOOTHING, CONVENTIONAL ALGORITHM,GENERAL
C       DESIGN, LOCAL BANDWIDTH ALLOWED
C
C  PARAMETERS :
C
C  INPUT    T(N)         INPUT GRID
C  INPUT    X(N)         DATA, GIVEN ON T(N)
C  INPUT    N            LENGTH OF X
C  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    KORD         ORDER OF KERNEL (<=6)
C  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
C  INPUT    NY           0: GLOBAL, 1: LOCAL BANDWIDTH INPUTED IN Y
C  INPUT    TT(M)        OUTPUT GRID
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE
C  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
C  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION X(N),T(N),S(0:N),TT(M),Y(M)
      DOUBLE PRECISION C(7),C1(7)
      INTEGER N,NUE,KORD,NY,M,IST,I,IBOUN,IORD
      DOUBLE PRECISION B, BB, BMAX, WID, S1, BMIN
C-
C------  COMPUTE KERNEL COEFFICIENTS FOR INTERIOR AND SOME CONSTANTS
      CALL COFFI(NUE,KORD,C)
      IORD=KORD+1
      BB=B
      S0=1.5*T(1)-0.5*T(2)
      SN=1.5*T(N)-0.5*T(N-1)
      BMIN=(SN-S0)*.6D0/DBLE(N)*DBLE(KORD-1)
      BMAX=(S(N)-S(0))*.5
      IF(KORD.EQ.2) BMIN=0.1D0*BMIN
      IST=1
C-
C-------  LOOP OVER OUTPUT GRID
      DO 100 I=1,M
        IF(NY.NE.0) BB=Y(I)
        IF(BB.GT.BMAX) BB=BMAX
        IF(BB.LT.BMIN) BB=BMIN
        WID=BB
        S1=TT(I)-BB
        IBOUN=0
C-
C-------  COMPUTE LEFT BOUNDARY KERNEL
        IF(S1.LT.S(0)) THEN
          S1=S(0)
          WID=BB+BB+S(0)-TT(I)
          CALL COFFB(NUE,KORD,(TT(I)-S(0))/WID,1,C1)
          IBOUN=1
        END IF
C-
C-------  COMPUTE RIGHT BOUNDARY KERNEL
        IF(TT(I)+BB.GT.S(N)) THEN
          S1=S(N)-(BB+BB)
          WID=TT(I)-S1
          CALL COFFB(NUE,KORD,(S(N)-TT(I))/WID,-1,C1)
          IBOUN=-1
        END IF
C-
C------  SEARCH FIRST S-POINT OF SMOOTHING INTERVAL
2       IF(S(IST).LE.S1) THEN
          IST=IST+1
          GOTO 2
        END IF
3       IF(S(IST-1).GT.S1) THEN
          IST=IST-1
          GOTO 3
        END IF
C-
C-------  IF BANDWIDTH IS TOO SMALL NO SMOOTHING
        IF(S(IST).GE.TT(I)+WID.OR.IST.EQ.N) THEN
         Y(I)=X(IST)
         IF(NUE.GT.0) Y(I)=0.
        ELSE
C-
C-----  COMPUTE SMOOTHED DATA AT TT(I)
         IF (IBOUN.NE.0) THEN
          CALL SMO(S,X,N,TT(I),WID,NUE,IORD,IBOUN,IST,S1,C1,Y(I))
         ELSE
          CALL SMO(S,X,N,TT(I),WID,NUE,IORD,IBOUN,IST,S1,C,Y(I))
         END IF
        END IF
100   CONTINUE
C-
      RETURN
      END

      SUBROUTINE SMO(S,X,N,TAU,WID,NUE,IORD,IBOUN,IST,S1,C,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION JANUARY 1995
C
C       PERFORMS ONE SMOOTHING STEP 
C
C  PARAMETERS :
C
C  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
C  INPUT    X(N)         DATA
C  INPUT    N            LENGTH OF X
C  INPUT    TAU          POINT WHERE FUNCTION IS ESTIMATED
C  INPUT    WID          ONE SIDED BANDWIDTH
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    IORD          ORDER OF KERNEL POLYNOMIAL
C  INPUT    IBOUN         TYPE OF BOUNDARY
C  INPUT    IST          INDEX OF FIRST POINT OF SMOOTHING INTERVAL
C  INPUT    S1           LEFT BOUNDARY OF SMOOTHING INTERVAL
C  INPUT    C(7)         KERNEL COEFFICIENTS 
C  OUTPUT   Y            SMOOTHED VALUE AT TAU
C  WORK     WO(7)        WORK ARRAY
C
C-----------------------------------------------------------------------
      DOUBLE PRECISION X(N),S(0:N),WO(7)
      DOUBLE PRECISION C(7)
      INTEGER N,NUE,IORD,IBOUN,IST,JEND,IBEG,INCR,I,J
      DOUBLE PRECISION TAU,WID,S1,Y,YY,YYY,W,WIDNUE
C-
      Y=0.
      JEND=0
      IBEG=2
      IF(IBOUN.NE.0.OR.(NUE.NE.1.AND.NUE.NE.3)) IBEG=1
      INCR=2
      IF(IBOUN.NE.0) INCR=1
C-
C------  COMPUTE INITIAL KERNEL VALUES
      IF(IBOUN.GT.0) THEN
        YY=(TAU-S1)/WID
        WO(IBEG)=YY
        DO 1 I=IBEG,IORD-INCR,INCR
1         WO(I+INCR)=WO(I)*YY
      ELSE
        DO 2 I=IBEG,IORD,INCR
2         WO(I)=1.
      END IF
C-
C------  LOOP OVER SMOOTHING INTERVAL
      DO 100 J=IST,N
        YY=(TAU-S(J))/WID
        IF(YY.LT.-1.) THEN
          YY=-1.
          JEND=1
        END IF
        YYY=YY
        IF(IBOUN.EQ.0) THEN
          YY=YY*YY
          IF(NUE.EQ.1.OR.NUE.EQ.3) YYY=YY
        END IF
C-
C------  LOOP FOR COMPUTING WEIGHTS
        W=0.
        DO 3 I=IBEG,IORD,INCR
          W=W+C(I)*(WO(I)-YYY)
          WO(I)=YYY
          YYY=YYY*YY
3         CONTINUE
        Y=Y+W*X(J)
        IF(JEND.EQ.1) GOTO 110
100     CONTINUE
C-
C-------  NORMALIZING FOR NUE>0
110   IF(NUE.GT.0) THEN
        WIDNUE=WID**NUE
        Y=Y/WIDNUE
      END IF
C-
      RETURN
      END

      SUBROUTINE KERNCP(T,X,N,B,NUE,KORD,NY,S,TT,M,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION JANUARY 1997
C
C       KERNEL SMOOTHING, CONVENTIONAL ALGORITHM,GENERAL
C       DESIGN, LOCAL BANDWIDTH ALLOWED, WITHOUT BOUNDARY KERNELS,
C       JUST NORMALIZING
C
C  PARAMETERS :
C
C  INPUT    T(N)         INPUT GRID
C  INPUT    X(N)         DATA, GIVEN ON T(N)
C  INPUT    N            LENGTH OF X
C  INPUT    B            ONE SIDED BANDWIDTH (FOR NY=1 MEAN BANDWIDTH)
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    KORD         ORDER OF KERNEL (<=6)
C  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
C  INPUT    NY           0: GLOBAL, 1: LOCAL BANDWIDTH INPUTED IN Y
C  INPUT    TT(M)        OUTPUT GRID
C  INPUT    M            NUMBER OF POINTS TO ESTIMATE
C  INPUT    Y(M)         BANDWITH SEQUENCE FOR NY=1, DUMMY FOR NY=0
C  OUTPUT   Y(M)         ESTIMATED REGRESSION FUNCTION
C
C-----------------------------------------------------------------------

      DOUBLE PRECISION X(N),T(N),S(0:N),TT(M),Y(M)
      DOUBLE PRECISION C(7),C1(7)
      INTEGER N,NUE,KORD,NY,M,IST,I,IBOUN,IORD
      DOUBLE PRECISION B, BB, BMAX, WID, S1, BMIN
C-
C------  COMPUTE KERNEL COEFFICIENTS FOR INTERIOR AND SOME CONSTANTS
      CALL COFFI(NUE,KORD,C)
      IORD=KORD+1
      BB=B
      S0=1.5*T(1)-0.5*T(2)
      SN=1.5*T(N)-0.5*T(N-1)
      BMIN=(SN-S0)*.6D0/DBLE(N)*DBLE(KORD-1)
      BMAX=(S(N)-S(0))*.5
      IF(KORD.EQ.2) BMIN=0.1D0*BMIN
      IST=1
C-
C-------  LOOP OVER OUTPUT GRID
      DO 100 I=1,M
        IF(NY.NE.0) BB=Y(I)
        IF(BB.GT.BMAX) BB=BMAX
        IF(BB.LT.BMIN) BB=BMIN
        WID=BB
        S1=TT(I)-BB
        IBOUN=0
C-
C-------  COMPUTE LEFT BOUNDARY KERNEL
        IF(S1.LT.S(0)) THEN
          S1=S(0)
          WID=BB+BB+S(0)-TT(I)
          CALL COFFB(NUE,KORD,(TT(I)-S(0))/WID,1,C1)
          IBOUN=1
        END IF
C-
C-------  COMPUTE RIGHT BOUNDARY KERNEL
        IF(TT(I)+BB.GT.S(N)) THEN
          S1=S(N)-(BB+BB)
          WID=TT(I)-S1
          IBOUN=-1
        END IF
C-
C------  SEARCH FIRST S-POINT OF SMOOTHING INTERVAL
2       IF(S(IST).LE.S1) THEN
          IST=IST+1
          GOTO 2
        END IF
3       IF(S(IST-1).GT.S1) THEN
          IST=IST-1
          GOTO 3
        END IF
C-
C-------  IF BANDWIDTH IS TOO SMALL NO SMOOTHING
        IF(S(IST).GE.TT(I)+WID.OR.IST.EQ.N) THEN
         Y(I)=X(IST)
         IF(NUE.GT.0) Y(I)=0.
        ELSE
C-
C-----  COMPUTE SMOOTHED DATA AT TT(I)
          CALL SMOP(S,X,N,TT(I),WID,NUE,IORD,IBOUN,IST,S1,C,Y(I))
        END IF
100   CONTINUE
C-
      RETURN
      END

      SUBROUTINE SMOP(S,X,N,TAU,WID,NUE,IORD,IBOUN,IST,S1,C,Y)
C-----------------------------------------------------------------------
C       SHORT-VERSION JANUARY 1995
C
C       PERFORMS ONE SMOOTHING STEP 
C
C  PARAMETERS :
C
C  INPUT    S(0:N)       HALF POINT INTERPOLATION SEQUENCE
C  INPUT    X(N)         DATA
C  INPUT    N            LENGTH OF X
C  INPUT    TAU          POINT WHERE FUNCTION IS ESTIMATED
C  INPUT    WID          ONE SIDED BANDWIDTH
C  INPUT    NUE          ORDER OF DERIVATIVE (0-4)
C  INPUT    IORD          ORDER OF KERNEL POLYNOMIAL
C  INPUT    IBOUN         TYPE OF BOUNDARY
C  INPUT    IST          INDEX OF FIRST POINT OF SMOOTHING INTERVAL
C  INPUT    S1           LEFT BOUNDARY OF SMOOTHING INTERVAL
C  INPUT    C(7)         KERNEL COEFFICIENTS 
C  OUTPUT   Y            SMOOTHED VALUE AT TAU
C  WORK     WO(7)        WORK ARRAY
C
C-----------------------------------------------------------------------
      DOUBLE PRECISION X(N),S(0:N),WO(7)
      DOUBLE PRECISION C(7)
      INTEGER N,NUE,IORD,IBOUN,IST,JEND,IBEG,INCR,I,J
      DOUBLE PRECISION TAU,WID,S1,Y,YY,YYY,W,WIDNUE,WW
C-
      Y=0.
      WW=0.
      JEND=0
      IBEG=2
      IF(NUE.NE.1.AND.NUE.NE.3) IBEG=1
      INCR=2
C-
C------  COMPUTE INITIAL KERNEL VALUES
      IF(IBOUN.GT.0) THEN
        YY=(TAU-S1)/WID
        WO(IBEG)=YY
        YY=YY*YY
        IF(NUE.EQ.1.OR.NUE.EQ.3) WO(IBEG)=YY        
        DO 1 I=IBEG,IORD-INCR,INCR
1         WO(I+INCR)=WO(I)*YY
      ELSE
        DO 2 I=IBEG,IORD,INCR
2         WO(I)=1.
      END IF
C-
C------  LOOP OVER SMOOTHING INTERVAL
      DO 100 J=IST,N
        YY=(TAU-S(J))/WID
        IF(YY.LT.-1.) THEN
          YY=-1.
          JEND=1
        END IF
        YYY=YY
        YY=YY*YY
        IF(NUE.EQ.1.OR.NUE.EQ.3) YYY=YY
C-
C------  LOOP FOR COMPUTING WEIGHTS
        W=0.
        DO 3 I=IBEG,IORD,INCR
          W=W+C(I)*(WO(I)-YYY)
          WO(I)=YYY
          YYY=YYY*YY
3         CONTINUE
        Y=Y+W*X(J)
        WW=WW+W
        IF(JEND.EQ.1) GOTO 110
100     CONTINUE
C-
C-------  NORMALIZING FOR NUE>0
110   IF(WW.NE.0) Y=Y/WW
      IF(NUE.GT.0) THEN
        WIDNUE=WID**NUE
        Y=Y/WIDNUE
      END IF
C-
      RETURN
      END

      SUBROUTINE COFFI(NUE,KORD,C)
C-----------------------------------------------------------------------
C       SHORT-VERSION: JANUARY 1995
C
C       PURPOSE:
C
C       DEFINES POLYNOMIAL KERNEL COEFFICIENTS FOR INTERIOR.
C
C  PARAMETERS:
C
C  INPUT  NUE        ORDER OF DERIVATIVE (0-4)
C  INPUT  KORD       ORDER OF KERNEL (NUE+I, I=2,4,6;  KORD<=6)
C  OUTPUT C(7)       POLYNOMIAL KERNEL COEFFICIENTS
C
C-----------------------------------------------------------------------
      DOUBLE PRECISION C(7)
C-
      DO 10 I=1,7
10      C(I)=0.
      IF(NUE.EQ.0.AND.KORD.EQ.2) THEN
          C(1)=0.75D0
          C(3)=-0.25D0
      END IF
C
      IF(NUE.EQ.0.AND.KORD.EQ.4) THEN
        C(1)=1.40625D0
        C(3)=-1.5625D0
        C(5)=0.65625D0
      END IF
C
      IF(NUE.EQ.0.AND.KORD.EQ.6) THEN
        C(1)=2.05078125D0
        C(3)=-4.78515625D0
        C(5)=5.16796875D0
        C(7)=-1.93359375D0
      END IF
C
      IF(NUE.EQ.1.AND.KORD.EQ.3) THEN
        C(2)=-1.875D0
        C(4)=0.9375D0
      END IF
C
      IF(NUE.EQ.1.AND.KORD.EQ.5) THEN
        C(2)=-8.203125D0
        C(4)=11.484375D0
        C(6)=-4.921875D0
      END IF
C
      IF(NUE.EQ.2.AND.KORD.EQ.4) THEN
        C(1)=-6.5625D0
        C(3)=13.125D0
        C(5)=-6.5625D0
      END IF
C
      IF(NUE.EQ.2.AND.KORD.EQ.6) THEN
        C(1)=-24.609375D0
        C(3)=103.359375D0
        C(5)=-132.890625D0
        C(7)=54.140625D0
      END IF
C
      IF(NUE.EQ.3.AND.KORD.EQ.5) THEN
        C(2)=88.59375D0
        C(4)=-147.65625D0
        C(6)=68.90625D0
      END IF
C
      IF(NUE.EQ.4.AND.KORD.EQ.6) THEN
        C(1)=324.84375D0
        C(3)=-1624.21875D0
        C(5)=2273.90625D0
        C(7)=-974.53125D0
      END IF
C
      RETURN
      END

      SUBROUTINE COFFB(NUE,KORD,Q,IBOUN,C)
C-----------------------------------------------------------------------
C       SHORT-VERSION: JANUARY 1995
C
C       PURPOSE:
C
C       COMPUTES COEFFICIENTS OF POLYNOMIAL BOUNDARY KERNELS,
C       FOLLOWING GASSER + MUELLER PREPRINT 38 SFB 123 HEIDELBERG
C       AND UNPUBLISHED RESULTS
C
C  PARAMETERS:
C
C  INPUT  NUE        ORDER OF DERIVATIVE (0-4)
C  INPUT  KORD       ORDER OF KERNEL (NUE+I, I=2,4,6;  KORD<=6)
C  INPUT  Q          PERCENTAGE OF WID AT BOUNDARY
C  INPUT  IBOUN      < 0 RIGHT BOUNDARY OF DATA
C                    > 0 LEFT BOUNDARY OF DATA
C  OUTPUT C(7)       POLYNOMIAL KERNEL COEFFICIENTS
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION C(7)
C
      DO 10 I=1,7
10       C(I)=0.
      P=-Q
      P1=1.+Q
      P3=P1*P1*P1
C
      IF(NUE.EQ.0.AND.KORD.EQ.2) THEN
        D=1./(P3*P1)
        C(1)=(6.+P*(12.+P*18.))*D
        C(2)=9.*(1.+P)*(1.+P)*D
        C(3)=(4.+P*8.)*D
      END IF
C
      IF(NUE.EQ.0.AND.KORD.EQ.4) THEN
        D=P1/(P3*P3*P3)
        P12=(1.+P)*(1.+P)*D
        C(1)=20.*(1.+P*(12.+P*(78.+P*(164.+P*(165.+P*(60.+P*10.))))))*D
        C(2)=100.*(1.+P*(5.+P))**2*P12
        C(3)=200.*(1.+P*(12.+P*(33.+P*(36.+P*(14.+P+P)))))*D
        C(4)=175.*(1.+P*(10.+P*3.))*P12
        C(5)=56.*(1.+P*(12.+P*(18.+P*4.)))*D
      END IF
C
      IF(NUE.EQ.0.AND.KORD.EQ.6) THEN
        P6=P3*P3
        D=1./(P6*P6)
        P12=(1.+P)*(1.+P)*D
        C(1)=42.*(1.+P*(30.+P*(465.+P*(3000.+P*(10050.+P*(17772.+P
     .       *(17430.+P*(9240.+P*(2625.+P*(350.+P*21.))))))))))*D
        C(2)=441.*(1.+P*(14.+P*(36.+P*(14.+P))))**2*P12
        C(3)=1960.*(1.+P*(30.+P*(255.+P*(984.+P*(1902.+P*(1956.+P
     .       *(1065.+P*(300.+P*(39.+P+P)))))))))*D
        C(4)=4410.*(1.+P*(28.+P*(156.+P*(308.+P*(188.+P*(42.+P*3.))))))
     .       *P12
        C(5)=5292.*(1.+P*(30.+P*(185.+P*(440.+P*(485.+P*(250.+P*(57.
     .       +P*4.)))))))*D
        C(6)=3234.*(1.+P*(28.+P*(108.+P*(56.+P*5.))))*P12
        C(7)=792.*(1.+P*(30.+P*(150.+P*(200.+P*(75.+P*6.)))))*D
       END IF
C
      IF(NUE.EQ.1.AND.KORD.EQ.3) THEN
        D=-1./(P3*P3)
        P12=(1.+P)*(1.+P)*D
        C(1)=(60.+P*240.)*P12
        C(2)=120.*(2.+P*(6.+P*(6.+P)))*D
        C(3)=300.*P12
        C(4)=(120.+P*180.)*D
      END IF
C
      IF(NUE.EQ.1.AND.KORD.EQ.5) THEN
        D=-1./(P3*P3*P3*P1)
        P12=(1.+P)*(1.+P)*D
        C(1)=420.*(1.+P*(18.+P*(98.+P*(176.+P*(75.+P*10.)))))*P12
        C(2)=2100.*(2.+P*(25.+P*(120.+P*(245.+P*(238.+P*(105.+P*(20.
     .       +P)))))))*D
        C(3)=14700.*(1.+P*(4.+P))**2*P12
        C(4)=5880.*(4.+P*(35.+P*(90.+P*(95.+P*(40.+P*6.)))))*D
        C(5)=17640.*(1.+P*(6.+P+P))*P12
        C(6)=2520.*(2.+P*(15.+P*(20.+P*5.)))*D
      END IF
C
      IF(NUE.EQ.2.AND.KORD.EQ.4) THEN
        D=P1/(P3*P3*P3)
        P12=(1.+P)*(1.+P)*D
        C(1)=840.*(1.+P*(12.+P*(28.+P*(24.+P*5.))))*D
        C(2)=2100.*(3.+P*(10.+P))*P12
        C(3)=1680.*(9.+P*(28.+P*(27.+P*6.)))*D
        C(4)=14700.*P12
        C(5)=(5040.+P*6720.)*D
      END IF
C
      IF(NUE.EQ.2.AND.KORD.EQ.6) THEN
        P6=P3*P3
        D=1./(P6*P6)
        P12=(1.+P)*(1.+P)*D
        C(1)=5040.*(2.+P*(60.+P*(489.+P*(1786.+P*(3195.+P*(2952.+P
     .       *(1365.+P*(294.+P*21.))))))))*D
        C(2)=52920.*(3.+P*(42.+P*(188.+P*(308.+P*(156.+P*(28.
     .       +P))))))*P12
        C(3)=141120.*(6.+P*(68.+P*(291.+P*(570.+P*(555.+P*(264.+P
     .       *(57.+P*4.)))))))*D
        C(4)=529200.*(2.+P*(7.+P+P))**2*P12
        C(5)=90720.*(30.+P*(228.+P*(559.+P*(582.+P*(255.
     .       +P*40.)))))*D
        C(6)=582120.*(3.+P*(14.+P*5.))*P12
        C(7)=221760.*(2.+P*(12.+P*(15.+P*4.)))*D
       END IF
C
      IF(NUE.EQ.3.AND.KORD.EQ.5) THEN
        D=-1./(P3*P3*P3*P1)
        P12=(1.+P)*(1.+P)*D
        C(1)=15120.*(1.+P*(18.+P*(38.+P*6.)))*P12
        C(2)=45360.*(4.+P*(35.+P*(80.+P*(70.+P*(20.+P)))))*D
        C(3)=352800.*(2.+P*(6.+P))*P12
        C(4)=151200.*(8.+P*(25.+P*(24.+P*6.)))*D
        C(5)=952560.*P12
        C(6)=70560.*(4.+P*5.)*D
      END IF
C
      IF(NUE.EQ.4.AND.KORD.EQ.6) THEN
        P6=P3*P3
        D=1./(P6*P6)
        P12=(1.+P)*(1.+P)*D
        C(1)=332640.*(1.+P*(30.+P*(171.+P*(340.+P*(285.+P*(90.+P*7.
     .      ))))))*D
        C(2)=1164240.*(5.+P*(56.+P*(108.+P*(28.+P))))*P12
        C(3)=6652800.*(5.+P*(38.+P*(85.+P*(76.+P*(25.+P+P)))))*D
        C(4)=17463600.*(5.+P*(14.+P*3.))*P12
        C(5)=4656960.*(25.+P*(78.+P*(75.+P*20.)))*D
        C(6)=76839840.*P12
        C(7)=3991680.*(5.+P*6.)*D
      END IF
C
      IF(IBOUN.GT.0) RETURN
      J=2
      IF(NUE.EQ.1.OR.NUE.EQ.3) J=1
      DO 2 I=J,KORD,2
2       C(I)=-C(I)
      RETURN
      END

      SUBROUTINE COFF(X,N,FA)
      DOUBLE PRECISION  X(N),FA
      INTEGER I,N
      DO 1 I=1,N
1       X(I)=FA
      RETURN
      END



