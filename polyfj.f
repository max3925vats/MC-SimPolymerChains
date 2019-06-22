      PROGRAM POLYHC   
C     MONTE CARLO SIMULATION OF HARD CHAIN
C     NEAR A LARGER HARD SPHERE
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      INTEGER HEAD
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (NBMAX=10000,NMAX=100,MAXBIN=500)
      PARAMETER (MCMAX=25,NCM=MCMAX*MCMAX*MCMAX,MAPMAX=27*NCM)
      DIMENSION X1(NBMAX),Y1(NBMAX),Z1(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION N1X(MAXBIN),NX1(MAXBIN),NX2(MAXBIN)
      DIMENSION NR(MAXBIN),SR(MAXBIN)
      DIMENSION NC(MAXBIN),RGX(MAXBIN),RGY(MAXBIN),RGZ(MAXBIN)
      DIMENSION REX(MAXBIN),REY(MAXBIN),REZ(MAXBIN)
      DIMENSION RG(MAXBIN),RE(MAXBIN)
      DIMENSION ACIG(MAXBIN),BCIG(MAXBIN),CCIG(MAXBIN)
      DIMENSION AI(3,3),DI(3),VI(3,3)
      DIMENSION HEAD(NCM),MAP(MAPMAX),LIST(NBMAX),IDI(NBMAX)
      COMMON /POS1/ X1,Y1,Z1,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /RON/ DLR, DINT
      COMMON /BOX/ AL,RSPL2
      COMMON /INTVAR/ NMOL1,N,NVL
      COMMON /SIGS/ DEFF,CELLI
      COMMON /OLAP/ HEAD,MAP,IDI,LIST,MC
C
C     GENERATE RANDOM NUMBER SEED 
      NSEED = 2*INT(SECNDS(0.0)) + 2937
C
C     READ RUN PARAMETERS 
C 
      OPEN(UNIT=1,FILE='polyfj.inp',STATUS='OLD')
C     TOTAL MOVES 
      READ(1,*)NCON
C     MOVES PER ACCUMULATION 
      READ(1,*)NSKIP
C     BIN SIZE 
      READ(1,*)BSZ
      BSZCM = 0.10D0
C     CHAIN TRANSLATION 
      READ(1,*)DLR
C     INTERNAL DISPLACEMENT 
      READ(1,*)DINT
C     FRACTION OF DICKMAN MOVES
      READ(1,*)FDICK
C     FRACTION OF REPTATION 
      READ(1,*)FREPT
      CLOSE(UNIT=1,STATUS='KEEP')

      OPEN(UNIT=2,FILE='polyfj.ic',STATUS='OLD')
      READ(2,*)PFC
      READ(2,*)
      READ(2,*)NMOL1
      READ(2,*)N
      READ(2,*)
      READ(2,*)RSP
      READ(2,*)AL
      READ(2,*)
      NCTOT = NMOL1*N
      DO I = 1, NCTOT
       READ(2,*)ID,JD,X1(I),Y1(I),Z1(I)
       X1(I) = X1(I)/AL
       Y1(I) = Y1(I)/AL
       Z1(I) = Z1(I)/AL
      ENDDO
      CLOSE(UNIT=2,STATUS='KEEP')

      NBIN = INT(AL/2.0D0/BSZ)
      NBINC = INT(AL/2.0D0/BSZCM)

      IF(NBIN.GT.MAXBIN)PAUSE 'Bin dimension'
C
      IF(MOD(N,2).EQ.0)THEN
       NM1 = N/2
       NM2 = NM1 + 1
      ELSE
       NM1 = (N-1)/2 + 1
       NM2 = NM1
      ENDIF
C
C     INITIALIZE BINS
      DO I = 1, NBIN
       N1X(I) = 0
       NX1(I) = 0
       NX2(I) = 0
       NR(I) = 0
       SR(I) = 0.0
      ENDDO
      DO I = 1, NBINC
       NC(I) = 0
       RGX(I) = 0.0D0
       RGY(I) = 0.0D0
       RGZ(I) = 0.0D0
       REX(I) = 0.0D0
       REY(I) = 0.0D0
       REZ(I) = 0.0D0
       RG(I) = 0.0D0
       RE(I) = 0.0D0
       ACIG(I) = 0.0D0
       BCIG(I) = 0.0D0
       CCIG(I) = 0.0D0
      ENDDO

C     PARAMETERS
      BSZR = BSZ/AL
      BSZC = BSZCM/AL
      DLR = DLR/AL
      DINT = DINT/AL
      DEFF = (1.0D0/AL)**2
      RSPL2 = ((RSP+0.5D0)/AL)**2

C     LINKED LIST
      DCMAX = 1.0D0
      MC = MIN(MCMAX,INT(AL/DCMAX))
      NCELL = MC*MC*MC
      CELLI = DBLE(MC)
      MAPSIZ = 27*NCELL
      DO I = 1, MAPSIZ
       MAP(I) = 0
      ENDDO
      DO ICELL = 1, NCELL
       ICX = MOD(ICELL,MC)
       IF(ICX.EQ.0)ICX=MC
       ICY = MOD( (ICELL-ICX)/MC , MC ) + 1
       ICZ = (ICELL-ICX-MC*(ICY-1))/MC/MC+1
       KCELL = 1
       DO IX = ICX-1,ICX+1
        IXR=IX
        IF(IXR.EQ.0) IXR = MC
        IF(IXR.EQ.MC+1) IXR = 1
        DO IY = ICY-1,ICY+1
         IYR = IY
         IF(IYR.EQ.0)IYR=MC
         IF(IYR.EQ.MC+1)IYR=1
         DO IZ = ICZ-1,ICZ+1
          IZR = IZ
          IF(IZR.EQ.0)IZR=MC
          IF(IZR.EQ.MC+1)IZR=1
          ICT = MC*MC*(IZR-1)+MC*(IYR-1)+IXR
          KD = (ICELL-1)*27 + KCELL
          MAP(KD) = ICT
          KCELL = KCELL + 1
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C
      DO I = 1, NCELL
       HEAD(I) = 0
      ENDDO
      DO I = 1, NBMAX
       LIST(I) = 0
      ENDDO

C     MAKE A LIST
      DO I = 1, NMOL1
       DO J = 1, N
        IJ = (I-1)*N+J
        X11 = X1(IJ)
        Y11 = Y1(IJ)
        Z11 = Z1(IJ)
        X11 = X11 - DNINT(X11) + 0.50D0
        Y11 = Y11 - DNINT(Y11) + 0.50D0
        Z11 = Z11 - DNINT(Z11) + 0.50D0
        ICELL = 1 + INT( X11*CELLI )
     C            + INT( Y11*CELLI ) *MC
     C            + INT( Z11*CELLI ) *MC*MC
        LIST(IJ) = HEAD(ICELL)
        HEAD(ICELL) = IJ
       ENDDO
      ENDDO

C     INITIALIZE COUNTERS 
      NTD = 0
      NTJ = 0
      NTR = 0
      NSD = 0
      NSJ = 0
      NSR = 0
      NAVER = 0
C
C     BEGIN SIMULATION
      DO 1000 ICON = 1, NCON

C     MAKE A MOVE
      
      I = INT( NMOL1*RAN(NSEED) ) + 1
      IMOL = (I-1)*N
      XRAN = RAN(NSEED)
C      write(*,*)'xran',xran
C     PICK DICKMAN REPTATION OR CCB
      IF(XRAN.LT.FREPT)THEN 
       NTR = NTR + 1
       CALL REPT(I,ISUC)
C       write(*,*)'rept',isuc
       IF(ISUC.EQ.0)NSR = NSR + 1
C       IF(ISUC.EQ.0)write(*,*)I,ISUC,XRAN,'REPT'
      ELSEIF(XRAN.LT.FDICK+FREPT)THEN
       NTD = NTD + 1
       CALL DICK(I,ISUC)
C       write(*,*)'dick',isuc
       IF(ISUC.EQ.0)NSD = NSD + 1
C       IF(ISUC.EQ.0)write(*,*)I,ISUC,XRAN,'DICKMANN'
      ELSE
       NTJ = NTJ + 1
       CALL CCB(I,ISUC)
C       write(*,*)'ccb',isuc
       IF(ISUC.EQ.0)NSJ = NSJ + 1
C       IF(ISUC.EQ.0)write(*,*)I,ISUC,XRAN,'CCB'
      ENDIF
C
      IF(ISUC.EQ.0)THEN
C     UPDATE POSITIONS AND LIST ARRAYS
       DO 10 J = 1, N
        IJ = IMOL + J
        XO = X1(IJ) - DNINT(X1(IJ)) + 0.50D0
        YO = Y1(IJ) - DNINT(Y1(IJ)) + 0.50D0
        ZO = Z1(IJ) - DNINT(Z1(IJ)) + 0.50D0
        XN = XITR(J) - DNINT(XITR(J)) + 0.50D0
        YN = YITR(J) - DNINT(YITR(J)) + 0.50D0
        ZN = ZITR(J) - DNINT(ZITR(J)) + 0.50D0
        ICO = 1 + INT( XO*CELLI )
     C          + INT( YO*CELLI ) *MC
     C          + INT( ZO*CELLI ) *MC*MC
        ICN = 1 + INT( XN*CELLI )
     C          + INT( YN*CELLI ) *MC
     C          + INT( ZN*CELLI ) *MC*MC

C       REMOVE X1(IJ) FROM LIST
        JMOL = HEAD(ICO)
        IF(JMOL.NE.IJ)THEN
         DO WHILE (LIST(JMOL).NE.IJ)
          JMOL = LIST(JMOL)
         ENDDO
         LIST(JMOL) = LIST(IJ)
        ELSE
         HEAD(ICO) = LIST(JMOL)
        ENDIF
C       ADD XITR(J) TO THE LIST
        JMOL = HEAD(ICN)
        IF(JMOL.NE.0)THEN
         DO WHILE (LIST(JMOL).NE.0)
          JMOL = LIST(JMOL)
         ENDDO
         LIST(JMOL) = IJ
        ELSE
         HEAD(ICN) = IJ
        ENDIF
        LIST(IJ) = 0
        X1(IJ) = XITR(J)
        Y1(IJ) = YITR(J)
        Z1(IJ) = ZITR(J)
10     CONTINUE
      ENDIF
C
      IF(MOD(ICON,NSKIP).EQ.0)THEN
       NAVER = NAVER + 1
C      DENSITY PROFILES
       DO I = 1, NMOL1
        DO J = 1, N
         IJ = (I-1)*N + J
         X1I = X1(IJ) - DNINT(X1(IJ))
         Y1I = Y1(IJ) - DNINT(Y1(IJ))
         Z1I = Z1(IJ) - DNINT(Z1(IJ))
         DIST = DSQRT(X1I**2+Y1I**2+Z1I**2)
         K = 1 + INT(DIST/BSZR)
         IF(K.LE.NBIN) THEN
          N1X(K) = N1X(K) + 1
          IF(J.EQ.1 .OR. J.EQ.N) NX1(K) = NX1(K) + 1
          IF(J.EQ.NM1 .OR. J.EQ.NM2) NX2(K) = NX2(K) + 1
         ENDIF
C
         IF(J.GE.2)THEN
          IJ1 = (I-1)*N + J
          IJ2 = (I-1)*N + J-1
          X11 = DSQRT(X1(IJ1)**2+Y1(IJ1)**2+Z1(IJ1)**2)
          X22 = DSQRT(X1(IJ2)**2+Y1(IJ2)**2+Z1(IJ2)**2)
          XM = DSQRT((X1(IJ2)-X1(IJ1))**2+(Y1(IJ2)-Y1(IJ1))**2
     C         +(Z1(IJ2)-Z1(IJ1))**2)/2.D0
          COSO = (X1(IJ1)*X1(IJ2)+Y1(IJ1)*Y1(IJ2)
     C           +Z1(IJ1)*Z1(IJ2))/(X11*X22)
           K = 1 + INT(XM/BSZR)
           NR(K) = NR(K) + 1
           SR(K) = SR(K) + COSO*AL*AL
          ENDIF
C
        ENDDO
       ENDDO
C
C      RG PROFILE
       DO I = 1, NMOL1
        XCM = 0.0D0
        YCM = 0.0D0
        ZCM = 0.0D0
        X2B = 0.0D0
        Y2B = 0.0D0
        Z2B = 0.0D0
        DO J = 1, N
         IJ = (I-1)*N + J
C         X1I = X1(IJ) - DNINT(X1(IJ))
C         Y1I = Y1(IJ) - DNINT(Y1(IJ))
C         Z1I = Z1(IJ) - DNINT(Z1(IJ))
         X1I = X1(IJ) 
         Y1I = Y1(IJ) 
         Z1I = Z1(IJ) 
         XCM = XCM + X1I
         YCM = YCM + Y1I
         ZCM = ZCM + Z1I
         X2B = X2B + X1I**2
         Y2B = Y2B + Y1I**2
         Z2B = Z2B + Z1I**2
        ENDDO
        XCM = XCM/DBLE(N)
        YCM = YCM/DBLE(N)
        ZCM = ZCM/DBLE(N)
        X1CM = XCM - DNINT(XCM)
        Y1CM = YCM - DNINT(YCM)
        Z1CM = ZCM - DNINT(ZCM)
        X2B = X2B/DBLE(N)
        Y2B = Y2B/DBLE(N)
        Z2B = Z2B/DBLE(N)
        XCMB = DSQRT(X1CM**2+Y1CM**2+Z1CM**2)
C
        I1 = (I-1)*N + 1
        IN = (I-1)*N + N
        XE = X1(IN) - X1(I1)
        YE = Y1(IN) - Y1(I1)
        ZE = Z1(IN) - Z1(I1)
        RE1I = XE*XE
        RE2I = YE*YE
        RE3I = ZE*ZE
C
C       CALCULATE SEMI-AXIS LENGTHS
C       CALCULATE MOMENT OF INERTIA TENSOR
        DO IA = 1, 3
         DO IB = 1, 3
          AI(IA,IB) = 0.0D0
         ENDDO
        ENDDO
        DO J = 1, N
         IJ = (I-1)*N + J
         RJX = X1(IJ) - XCM
         RJY = Y1(IJ) - YCM
         RJZ = Z1(IJ) - ZCM
         AI(1,1) = AI(1,1) + RJY**2 + RJZ**2
         AI(2,2) = AI(2,2) + RJX**2 + RJZ**2
         AI(3,3) = AI(3,3) + RJX**2 + RJY**2
         AI(1,2) = AI(1,2) - RJX*RJY
         AI(1,3) = AI(1,3) - RJX*RJZ
         AI(2,3) = AI(2,3) - RJY*RJZ
        ENDDO
        AI(3,2) = AI(2,3)
        AI(3,1) = AI(1,3)
        AI(2,1) = AI(1,2)
C       DIAGONALIZE IT
        CALL JACOBI(AI,3,3,DI,VI,NROT)
C       FIND THE SMALLEST EIGENVALUE
        EIG = 1.D+10
        DO K = 1, 3
         IF(DI(K).LT.EIG)THEN
          IEIG = K
          EIG = DI(K)
         ENDIF
        ENDDO
        EIG1 = EIG
        IF(IEIG.EQ.1)EIG2 = MIN(DI(2),DI(3))
        IF(IEIG.EQ.1)EIG3 = MAX(DI(2),DI(3))
        IF(IEIG.EQ.2)EIG2 = MIN(DI(1),DI(3))
        IF(IEIG.EQ.2)EIG3 = MAX(DI(1),DI(3))
        IF(IEIG.EQ.3)EIG2 = MIN(DI(1),DI(2))
        IF(IEIG.EQ.3)EIG3 = MAX(DI(1),DI(2))
        CIGAI = DSQRT(2.5*(EIG2+EIG3-EIG1)/DFLOAT(N))
        CIGBI = DSQRT(2.5*(EIG1+EIG3-EIG2)/DFLOAT(N))
        CIGCI = DSQRT(2.5*(EIG2+EIG1-EIG3)/DFLOAT(N))
C
        K = 1 + INT(XCMB/BSZC)
        IF(K.LE.NBINC)THEN
         RG1I = X2B - X1CM*X1CM
         RG2I = Y2B - Y1CM*Y1CM
         RG3I = Z2B - Z1CM*Z1CM
         RGX(K) = RGX(K) + RG1I*AL*AL
         RGY(K) = RGY(K) + RG2I*AL*AL
         RGZ(K) = RGZ(K) + RG3I*AL*AL
         REX(K) = REX(K) + RE1I*AL*AL
         REY(K) = REY(K) + RE2I*AL*AL
         REZ(K) = REZ(K) + RE3I*AL*AL
         RG(K) = RG(K) + (RG1I+RG2I+RG3I)*AL*AL
         RE(K) = RE(K) + (RE1I+RE2I+RE3I)*AL*AL
         NC(K) = NC(K) + 1
         ACIG(K) = ACIG(K) + CIGAI*AL
         BCIG(K) = BCIG(K) + CIGBI*AL
         CCIG(K) = CCIG(K) + CIGCI*AL
        ENDIF
       ENDDO
      ENDIF
1000  CONTINUE
C     WRITE RESULTS
      IF(NTD.GT.0)FSD=DBLE(NSD)/DBLE(NTD)
      IF(NTJ.GT.0)FSJ=DBLE(NSJ)/DBLE(NTJ)
      IF(NTR.GT.0)FSR=DBLE(NSR)/DBLE(NTR)
      DO I = 1, NCTOT
       X1(I) = AL*X1(I)
       Y1(I) = AL*Y1(I)
       Z1(I) = AL*Z1(I)
      ENDDO

      DLR = DLR*AL
      DINT = DINT*AL
C
      OPEN(UNIT=2,FILE='polyfj.fc',STATUS='UNKNOWN')
       WRITE(2,111)PFC,NMOL1,N,RSP,AL
       DO I = 1, NMOL1
        DO J = 1, N
         IJ = (I-1)*N + J
         WRITE(2,112)I,J,X1(IJ),Y1(IJ),Z1(IJ)
        ENDDO
       ENDDO
       CLOSE(UNIT=2,STATUS='KEEP')
111    FORMAT(3X,F10.4,3X,'PACKING FRACTION' /
     C       /
     C      2X, I8, 6X,'NUMBER OF MOLECULES' /
     C      2X, I6, 6X,'CHAIN LENGTH' /
     C      /
     C      2X, F10.4,5X,'RADIUS OF SOLUTE'/
     C      4X, E18.10, 2X, 'PERIODIC LENGTH '/ )
112   FORMAT(1X,I6,2X,I6,2X,3E18.10)
C
      OPEN(UNIT=1,FILE='polyfj.out',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='polyden.out',STATUS='UNKNOWN')
C
      VOL = AL**3 - 4.D0*PI*RSP**3/3.D0
      PFCMC = (PI/6.)*NCTOT/VOL
C
      WRITE(1,*)'RESULTS OF CHAIN NEAR A SPHERE RUN'
      WRITE(1,*)
      WRITE(1,113)PFCMC,N,NMOL1,DLR,DINT,RSP
      WRITE (1,114) NCON, NSKIP
      WRITE (1,116) BSZ,NBIN
      WRITE(1,*)'Move statistics'
      WRITE (1,117) NTD, NTR, NTJ, NSD, NSR, NSJ,
     C              FSD, FSR, FSJ
      WRITE (1,*)'DENSITY PROFILES NEAR A SPHERE'
      CONST = 4.0D0*PI/3.0D0
      DO I=1, NBIN
       RLOWER = DBLE(I-1)*BSZ
       RUPPER = RLOWER + BSZ
       VOL = CONST*(RUPPER**3 - RLOWER**3)*DBLE(NAVER)
       DIST = RLOWER + BSZ/2.0
       IF(NAVER.GT.0)AN1X = DBLE(N1X(I))/VOL
       IF(NAVER.GT.0)ANX1 = DBLE(NX1(I))/VOL
       IF(NAVER.GT.0)ANX2 = DBLE(NX2(I))/VOL
       IF(NAVER.GT.0)WRITE (1,1140)I,DIST,AN1X,ANX1,ANX2
       IF(NAVER.GT.0)WRITE (4,1140)I,DIST,AN1X,ANX1,ANX2
      ENDDO
      CLOSE(UNIT=4,STATUS='KEEP')
C
C     SEGMENTAL ORDER PARAMETER
      WRITE(1,*)'SEGMENTAL ORDER PARAMETER'
      OPEN(UNIT=4,FILE='seg.out',STATUS='UNKNOWN')
      DO I=1, NBIN
       RLOWER = DBLE(I-1)*BSZ
       DIST = RLOWER + BSZ/2.0
       IF(NR(I).GT.0)SRI = SR(I)/DFLOAT(NR(I))
       SRI = (3.0D0*SRI-1.0D0)/2.0D0
       IF(NR(I).GT.0)WRITE (4,1140)I,DIST,SRI
       IF(NR(I).GT.0)WRITE (1,1140)I,DIST,SRI
      ENDDO
      CLOSE(UNIT=4,STATUS='KEEP')
C
      OPEN(UNIT=4,FILE='rgfj.out',STATUS='UNKNOWN')
      OPEN(UNIT=2,FILE='g2fj.out',STATUS='UNKNOWN')
      WRITE(1,*)'CENTRE OF MASS AND R2G PROFILE'
      WRITE(1,*)
      WRITE(1,113)PFC,N,NMOL1,DLR,DINT,RSP
      WRITE (1,114) NCON, NSKIP
      WRITE (1,118) BSZ,BSZCM,NBIN,NBINC
      WRITE(1,*)'Move statistics'
      WRITE (1,117) NTD, NTR, NTJ, NSD, NSR, NSJ,
     C              FSD, FSR, FSJ
      CONST = (NBINC*BSZCM)**3 
      DO I = 1, NBINC
       RLOWER = DBLE(I-1)*BSZCM
       RUPPER = RLOWER + BSZCM
       VOL = (RUPPER**3 - RLOWER**3)*DBLE(NAVER)/CONST
       DIST = RLOWER + BSZCM/2.0
       IF(NAVER.GT.0)ANC = DBLE(NC(I))/VOL/NMOL1
       IF(NC(I).GT.0)RGXI = RGX(I)/DBLE(NC(I))
       IF(NC(I).GT.0)RGYI = RGY(I)/DBLE(NC(I))
       IF(NC(I).GT.0)RGZI = RGZ(I)/DBLE(NC(I))
       IF(NC(I).GT.0)REXI = REX(I)/DBLE(NC(I))
       IF(NC(I).GT.0)REYI = REY(I)/DBLE(NC(I))
       IF(NC(I).GT.0)REZI = REZ(I)/DBLE(NC(I))
       IF(NC(I).GT.0)RGI = RG(I)/DBLE(NC(I))
       IF(NC(I).GT.0)REI = RE(I)/DBLE(NC(I))
       IF(NC(I).GT.0) ACIGI = DSQRT(ACIG(I)/DBLE(NC(I)))
       IF(NC(I).GT.0) BCIGI = DSQRT(BCIG(I)/DBLE(NC(I)))
       IF(NC(I).GT.0) CCIGI = DSQRT(CCIG(I)/DBLE(NC(I)))
       WRITE(1,1140)I,DIST,ANC,RGXI,RGYI,RGZI,REXI,REYI,REZI,RGI,REI
       WRITE(4,1140)I,DIST,ANC,RGXI,RGYI,RGZI,REXI,REYI,REZI,RGI,REI
       WRITE(1,1140)I,DIST,ACIGI,BCIGI,CCIGI
       WRITE(2,1140)I,DIST,ACIGI,BCIGI,CCIGI
      ENDDO
      CLOSE(UNIT=4,STATUS='KEEP')
      CLOSE(UNIT=2,STATUS='KEEP')
113   FORMAT(1X,'PACKING FRACTION ',4X,F10.4 /
     C       1X,' NUMBER OF BEADS ',6X, I6  /
     C       1X,'NUMBER OF CHAINS ',6X, I8  /
     C       1X,'     DELTA CHAIN ',5X, F10.4, /
     C       1X,'  DELTA INTERNAL ',5X, F10.4, /
     C       1X,'RADIUS OF LARGER SPHERE',5X, F10.4 )
114   FORMAT (1X / 1X, 'NO. TOTAL MOVES', 2X, I8 /
     C             1X, 'NUMBER OF SKIPS', 3X, I7 /)
116   FORMAT (1X,  '      BIN SIZE', 3X, F10.4/
     C        1X, 'NUMBER OF BINS', 4X, I6 / )
118   FORMAT (1X,  '      BIN SIZE', 3X, F10.4, 2X, F10.4 /
     C        1X, 'NUMBER OF BINS', 4X,I6 , 3X, I6 / )
117   FORMAT (1X / 13X, 'TRANSLATION DICKMAN  REPTATION       CCB' /
     C        1X, 'ATTEMPTED', 3X, I10, 1X, I10, 1X, I10 /
     C        1X, 'SUCCESSFUL', 2X, I10, 1X, I10, 1X, I10 /
     C 1X, 'FRACTION',5X, F10.5, 2X, F10.5, 2X, F10.5 /) 
1140  FORMAT (1X,I10,1X,F12.5,3X,9(E10.4,4X,E10.4,4X))
      STOP
      END
C
      SUBROUTINE DICK(I,ISUC)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      INTEGER HEAD
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (NBMAX=10000,NMAX=100)
      PARAMETER (MCMAX=25,NCM=MCMAX*MCMAX*MCMAX,MAPMAX=27*NCM)
      DIMENSION X1(NBMAX),Y1(NBMAX),Z1(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION HEAD(NCM),MAP(MAPMAX),LIST(NBMAX),IDI(NBMAX)
      COMMON /POS1/ X1,Y1,Z1,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /RON/ DLR, DINT
      COMMON /BOX/ AL,RSPL2
      COMMON /INTVAR/ NMOL1,N,NVL
      COMMON /SIGS/ DEFF,CELLI
      COMMON /OLAP/ HEAD,MAP,IDI,LIST,MC
      ISUC = 1
      ZS1 = RAN(NSEED) - .5
      ZS2 = RAN(NSEED) - .5
      ZS3 = RAN(NSEED) - .5
      IMOL = (I-1)*N
      IBD = IMOL + 1
C     DISPLACE FIRST SEGMENT
      XITR(1) = X1(IBD) + ZS1*DLR
      YITR(1) = Y1(IBD) + ZS2*DLR
      ZITR(1) = Z1(IBD) + ZS3*DLR
C     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
      XT = XITR(1) - DNINT(XITR(1))
      YT = YITR(1) - DNINT(YITR(1))
      ZT = ZITR(1) - DNINT(ZITR(1))
      RDIST = XT**2+YT**2+ZT**2
      IF(RDIST.LT.RSPL2)RETURN
C     CHECK FOR INTERMOLECULAR OVERLAP
      CALL OVER(I,1)
      IF(NVL.EQ.1)RETURN
C     DISPLACE SUBSEQUENT SEGMENTS
      DO 200 II = 2, N
       IBD = IMOL + II
       D1 = X1(IBD)-X1(IBD-1)
       D2 = Y1(IBD)-Y1(IBD-1)
       D3 = Z1(IBD)-Z1(IBD-1)
       IF(RAN(NSEED).GT.0.25)GOTO 211
       CALL RUV(ZS1,ZS2,ZS3)
       D1 = D1 + DINT*ZS1
       D2 = D2 + DINT*ZS2
       D3 = D3 + DINT*ZS3
       CNM = D1*D1 + D2*D2 + D3*D3
       CN = 1.0D0/DSQRT(CNM)/AL
       D1 = D1*CN
       D2 = D2*CN
       D3 = D3*CN
211    CONTINUE
       XITR(II) = XITR(II-1) + D1
       YITR(II) = YITR(II-1) + D2
       ZITR(II) = ZITR(II-1) + D3
C     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
      XT = XITR(II) - DNINT(XITR(II))
      YT = YITR(II) - DNINT(YITR(II))
      ZT = ZITR(II) - DNINT(ZITR(II))
      RDIST = XT**2+YT**2+ZT**2
      IF(RDIST.LT.RSPL2)RETURN
C      CHECK FOR INTERMOLECULAR OVERLAP
       CALL OVER(I,II)
       IF(NVL.EQ.1)RETURN
C      CHECK FOR INTRA OVERLAP
       IF(II.GT.2)THEN
        DO 215 JJ = 1, II - 2
         DX = XITR(JJ)-XITR(II)
         DY = YITR(JJ)-YITR(II)
         DZ = ZITR(JJ)-ZITR(II)
         DX = DX - DNINT(DX)
         DY = DY - DNINT(DY)
         DZ = DZ - DNINT(DZ)
         RDIS = DX*DX+DY*DY+DZ*DZ
         IF( RDIS .LT. DEFF ) RETURN
215     CONTINUE
       ENDIF
200   CONTINUE
C     MOVE ACCEPTED
      ISUC = 0
      RETURN
      END
C
      SUBROUTINE REPT(I,ISUC)
C     PERFORMS A REPTATION MOVE ON LINEAR CHAIN I
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      INTEGER HEAD
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (NBMAX=10000,NMAX=100)
      PARAMETER (MCMAX=25,NCM=MCMAX*MCMAX*MCMAX,MAPMAX=27*NCM)
      DIMENSION X1(NBMAX),Y1(NBMAX),Z1(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION HEAD(NCM),MAP(MAPMAX),LIST(NBMAX),IDI(NBMAX)
      COMMON /POS1/ X1,Y1,Z1,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /RON/ DLR, DINT
      COMMON /BOX/ AL,RSPL2
      COMMON /INTVAR/ NMOL1,N,NVL
      COMMON /SIGS/ DEFF,CELLI
      COMMON /OLAP/ HEAD,MAP,IDI,LIST,MC
C
      ISUC = 1
      IF(RAN(NSEED).GT.0.5)CALL CVERT(I)
      IMOL = (I-1)*N
C
C     CUT OFF END AND ATTACH TO BEAD 1
      CALL RUV(DELX,DELY,DELZ)
      XITR(1) = X1(IMOL+1) + DELX /AL
      YITR(1) = Y1(IMOL+1) + DELY /AL
      ZITR(1) = Z1(IMOL+1) + DELZ /AL
C     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
      XT = XITR(1) - DNINT(XITR(1))
      YT = YITR(1) - DNINT(YITR(1))
      ZT = ZITR(1) - DNINT(ZITR(1))
      RDIST = XT**2+YT**2+ZT**2
      IF(RDIST.LT.RSPL2)RETURN
C     CHECK FOR INTER CHAIN OVERLAP
      CALL OVER(I,1)
      IF(NVL.EQ.1)RETURN
C     CHECK FOR INTRA CHAIN OVERLAP
      DO 320 II = 2, N-1
       I1 = IMOL + II
       DX = XITR(1)-X1(I1)
       DY = YITR(1)-Y1(I1)
       DZ = ZITR(1)-Z1(I1)
       DX = DX - DNINT(DX)
       DY = DY - DNINT(DY)
       DZ = DZ - DNINT(DZ)
       RDIS = DX*DX+DY*DY+DZ*DZ
       IF( RDIS .LT. DEFF )RETURN
320   CONTINUE
C     MOVE ACCEPTED
C     RESET NEW CHAIN BEADS
      DO J = 2, N
       I1 = IMOL + J-1
       XITR(J) = X1(I1)
       YITR(J) = Y1(I1)
       ZITR(J) = Z1(I1)
      ENDDO
      ISUC = 0
      RETURN
      END
C
      SUBROUTINE CCB(I,ISUC)
C
C     PERFORMS A CCB MOVE ON CHAIN I
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      INTEGER HEAD
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER ( NSAMP=15 )
      PARAMETER (NBMAX=10000,NMAX=100)
      PARAMETER (MCMAX=25,NCM=MCMAX*MCMAX*MCMAX,MAPMAX=27*NCM)
      DIMENSION X1(NBMAX),Y1(NBMAX),Z1(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION XT(NSAMP),YT(NSAMP),ZT(NSAMP),ET(NSAMP),
     C  WS(NMAX),STX(NMAX),STY(NMAX),STZ(NMAX),
     C  XN(NMAX),YN(NMAX),ZN(NMAX)
      DIMENSION HEAD(NCM),MAP(MAPMAX),LIST(NBMAX),IDI(NBMAX)
      COMMON /POS1/ X1,Y1,Z1,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /RON/ DLR, DINT
      COMMON /BOX/ AL,RSPL2
      COMMON /INTVAR/ NMOL1,N,NVL
      COMMON /SIGS/ DEFF,CELLI
      COMMON /OLAP/ HEAD,MAP,IDI,LIST,MC
C
      ISUC = 1
      WN = 1.0D0
      WO = 1.0D0
      IF(RAN(NSEED).GT.0.5)CALL CVERT(I)
      IMOL = (I-1)*N
      DO J = 1, N
       I1 = IMOL + J
       XN(J) = X1(I1)
       YN(J) = Y1(I1)
       ZN(J) = Z1(I1)
      ENDDO
      ICUT = INT((N-1)*RAN(NSEED)) + 2
C      write(*,*)ICUT
      DO 50 J = ICUT, N
       SUM = 0.0D0
       DO 20 K = 1, NSAMP
        ET(K) = 1.0D0
        CALL RUV(DELX,DELY,DELZ)
        XT(K) = XN(J-1) + DELX/AL
        YT(K) = YN(J-1) + DELY/AL
        ZT(K) = ZN(J-1) + DELZ/AL
C     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
        XT1 = XT(K) - DNINT(XT(K))
        YT1 = YT(K) - DNINT(YT(K))
        ZT1 = ZT(K) - DNINT(ZT(K))
        RDIST = XT1**2+YT1**2+ZT1**2
        IF(RDIST.LT.RSPL2)THEN
         ET(K) = 0.0D0
         GOTO 19
        ENDIF
C      CHECK FOR INTER CHAIN OVERLAP
        XITR(J) = XT(K)
        YITR(J) = YT(K)
        ZITR(J) = ZT(K)
        CALL OVER(I,J)
        IF(NVL.EQ.1)THEN
         ET(K) = 0.0D0
         GOTO 19
        ENDIF
C     CHECK FOR INTRA CHAIN OVERLAP
        IF(J.GE.3)THEN
         DO 16 II = 1, J - 2
          DX = XT(K)-XN(II)
          DY = YT(K)-YN(II)
          DZ = ZT(K)-ZN(II)
          DX = DX - DNINT(DX)
          DY = DY - DNINT(DY)
          DZ = DZ - DNINT(DZ)
          RDIS = DX*DX+DY*DY+DZ*DZ
          IF( RDIS .LT. DEFF )THEN
           ET(K) = 0.0D0
           GOTO 19
          ENDIF
16       CONTINUE
        ENDIF
19      CONTINUE
        SUM = SUM + ET(K)
20     CONTINUE
       IF(SUM.LT.1.0D-10)RETURN
       SUM = 1.0D0/SUM
       DO 25 K = 1, NSAMP
        ET(K) = ET(K)*SUM
25     CONTINUE
       XRAN = RAN(NSEED)
       S = 0.0D0
       DO 30 K = 1, NSAMP
        S = S + ET(K)
        IF(XRAN.LT.S)GOTO 35
30     CONTINUE
35     CONTINUE
       WS(J) = ET(K)
       XN(J) = XT(K)
       YN(J) = YT(K)
       ZN(J) = ZT(K)
       WN = WN*WS(J)
50    CONTINUE
      DO 60 J = 1, N
       STX(J) = XN(J)
       STY(J) = YN(J)
       STZ(J) = ZN(J)
60    CONTINUE
C  COMPUTE WEIGHT OF EXISTING CHAIN
      DO 65 J = 1, N
       I1 = IMOL + J
       XN(J) = X1(I1)
       YN(J) = Y1(I1)
       ZN(J) = Z1(I1)
65    CONTINUE

      DO 100 J = ICUT, N
       SUM = 0.0D0
       DO 80 K = 1, NSAMP
        ET(K) = 1.0D0
        IF(K.EQ.1) GOTO 79
        CALL RUV(DELX,DELY,DELZ)
        XT(K) = XN(J-1) + DELX /AL
        YT(K) = YN(J-1) + DELY/ AL
        ZT(K) = ZN(J-1) + DELZ/ AL
C     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
        XT1 = XT(K) - DNINT(XT(K))
        YT1 = YT(K) - DNINT(YT(K))
        ZT1 = ZT(K) - DNINT(ZT(K))
        RDIST = XT1**2+YT1**2+ZT1**2
        IF(RDIST.LT.RSPL2)THEN
         ET(K) = 0.0D0
         GOTO 79
        ENDIF
C      CHECK FOR INTER CHAIN OVERLAP
        XITR(J) = XT(K)
        YITR(J) = YT(K)
        ZITR(J) = ZT(K)
        CALL OVER(I,J)
        IF(NVL.EQ.1)THEN
         ET(K) = 0.0D0
         GOTO 79
        ENDIF
C     CHECK FOR INTRA CHAIN OVERLAP
        IF(J.GE.3)THEN
         DO 75 II = 1, J - 2
          DX = XT(K)-XN(II)
          DY = YT(K)-YN(II)
          DZ = ZT(K)-ZN(II)
          DX = DX - DNINT(DX)
          DY = DY - DNINT(DY)
          DZ = DZ - DNINT(DZ)
          RDIS = DX*DX+DY*DY+DZ*DZ
          IF( RDIS .LT. DEFF )THEN
           ET(K) = 0.0D0
           GOTO 79
          ENDIF
75       CONTINUE
        ENDIF
79      CONTINUE
        SUM = SUM + ET(K)
80     CONTINUE
       ET(1) = ET(1)/SUM
       WO = WO*ET(1)
100   CONTINUE
      BOLTZ = WO/WN
      IF(BOLTZ.GT.RAN(NSEED))THEN
       DO J = 1, N
        XITR(J) = STX(J)
        YITR(J) = STY(J)
        ZITR(J) = STZ(J)
       ENDDO
       ISUC = 0
C
      ENDIF 
C 
      RETURN
      END
C
      SUBROUTINE OVER1(ITRY,JTRY)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER HEAD
      PARAMETER (NBMAX=10000,NMAX=100,MCMAX=25)
      PARAMETER (NCM=MCMAX*MCMAX*MCMAX,MAPMAX=27*NCM)
      DIMENSION X1(NBMAX),Y1(NBMAX),Z1(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION HEAD(NCM),MAP(MAPMAX),LIST(NBMAX),IDI(NBMAX)
      DIMENSION JNEAR(NBMAX)
      COMMON /POS1/ X1,Y1,Z1,XITR,YITR,ZITR
      COMMON /BOX/ AL,RSPL2
      COMMON /INTVAR/ NMOL1,N,NVL
      COMMON /SIGS/ DEFF,CELLI
      COMMON /OLAP/ HEAD,MAP,IDI,LIST,MC
      NVL = 1
      X11 = XITR(JTRY)
      Y11 = YITR(JTRY)
      Z11 = ZITR(JTRY)
      X11 = X11 - DNINT(X11) + 0.5D0
      Y11 = Y11 - DNINT(Y11) + 0.5D0
      Z11 = Z11 - DNINT(Z11) + 0.5D0
      ICELL = 1 + INT( X11*CELLI )
     C          + INT( Y11*CELLI ) *MC
     C          + INT( Z11*CELLI ) *MC*MC
      ID = (ICELL-1)*27
      NITER = 0
      DO NC = 1, 27
       JCELL = MAP(ID+NC)
       IF(JCELL.EQ.0)GOTO 20
       IM = HEAD(JCELL)
       DO WHILE (IM.NE.0)
        IF(IDI(IM).EQ.ITRY)GOTO 10
        NITER = NITER + 1
        JNEAR(NITER) = IM
10      CONTINUE
        IM = LIST(IM)
       ENDDO
      ENDDO
20    CONTINUE
      DO NN = 1, NITER
       IM = JNEAR(NN)
       XT = X1(IM) - XITR(JTRY)
       YT = Y1(IM) - YITR(JTRY)
       ZT = Z1(IM) - ZITR(JTRY)
       XT = XT - DNINT(XT)
       YT = YT - DNINT(YT)
       ZT = ZT - DNINT(ZT)
       RDIS = XT*XT+YT*YT+ZT*ZT
       IF(RDIS.LT.DEFF)RETURN
      ENDDO
C
      NVL = 0
      RETURN
      END
C
      SUBROUTINE OVER(ITRY,JTRY)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NBMAX=10000,NMAX=100)
      DIMENSION X1(NBMAX),Y1(NBMAX),Z1(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      COMMON /POS1/ X1,Y1,Z1,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /RON/ DLR, DINT
      COMMON /BOX/ AL,RSPL2
      COMMON /INTVAR/ NMOL1,N,NVL
      COMMON /SIGS/ DEFF,CELLI
      NVL = 1
      DO I = 1, NMOL1
       DO J = 1, N
        IJ = (I-1)*N + J
        IF (I .NE. ITRY) THEN
         XT = X1(IJ) - XITR(JTRY)
         YT = Y1(IJ) - YITR(JTRY)
         ZT = Z1(IJ) - ZITR(JTRY)
         XT = XT - DNINT(XT)
         YT = YT - DNINT(YT)
         ZT = ZT - DNINT(ZT)
         RDIS = XT*XT+YT*YT+ZT*ZT 
         IF(RDIS.LT.DEFF)RETURN
        ENDIF
       ENDDO
      ENDDO
      NVL = 0
      RETURN
      END
C
      SUBROUTINE RUV(A,B,C)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      COMMON /SEED/ NSEED

1     B1 = 1.0D0 - 2.0D0*RAN(NSEED)
      B2 = 1.0D0 - 2.0D0*RAN(NSEED)
      BSQ = B1*B1 + B2*B2
      IF(BSQ.GT.1.0D0) THEN
C      REJECT
       GOTO 1
      ELSE
       BH = DSQRT(1.D0 - BSQ)
       A = 2.D0 * B1 * BH
       B = 2.D0 * B2 * BH
       C = 1.D0 - 2.D0*BSQ
      ENDIF
      RETURN
      END
C
      SUBROUTINE CVERT(I)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      INTEGER HEAD
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (NBMAX=10000,NMAX=100)
      PARAMETER (MCMAX=25,NCM=MCMAX*MCMAX*MCMAX,MAPMAX=27*NCM)
      DIMENSION X1(NBMAX),Y1(NBMAX),Z1(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION HEAD(NCM),MAP(MAPMAX),LIST(NBMAX),IDI(NBMAX)
      COMMON /POS1/ X1,Y1,Z1,XITR,YITR,ZITR
      COMMON /INTVAR/ NMOL1,N,NVL
      COMMON /SIGS/ DEFF,CELLI
      COMMON /OLAP/ HEAD,MAP,IDI,LIST,MC
C
      IMOL = (I-1)*N
      DO J = 1, N
       K = N - J + 1
       IK = IMOL + K
       XITR(J) = X1(IK)
       YITR(J) = Y1(IK)
       ZITR(J) = Z1(IK)
      ENDDO
      DO J = 1, N
       IJ = IMOL + J
        XO = X1(IJ) - DNINT(X1(IJ)) + 0.50D0
        YO = Y1(IJ) - DNINT(Y1(IJ)) + 0.50D0
        ZO = Z1(IJ) - DNINT(Z1(IJ)) + 0.50D0
        XN = XITR(J) - DNINT(XITR(J)) + 0.50D0
        YN = YITR(J) - DNINT(YITR(J)) + 0.50D0
        ZN = ZITR(J) - DNINT(ZITR(J)) + 0.50D0
        ICO = 1 + INT( XO*CELLI )
     C          + INT( YO*CELLI ) *MC
     C          + INT( ZO*CELLI ) *MC*MC
        ICN = 1 + INT( XN*CELLI )
     C          + INT( YN*CELLI ) *MC
     C          + INT( ZN*CELLI ) *MC*MC
C       REMOVE X1(IJ) FROM LIST
        JMOL = HEAD(ICO)
        IF(JMOL.NE.IJ)THEN
         DO WHILE (LIST(JMOL).NE.IJ)
          JMOL = LIST(JMOL)
         ENDDO
         LIST(JMOL) = LIST(IJ)
        ELSE
         HEAD(ICO) = LIST(JMOL)
        ENDIF
C       ADD XITR(J) TO THE LIST
        JMOL = HEAD(ICN)
        IF(JMOL.NE.0)THEN
         DO WHILE (LIST(JMOL).NE.0)
          JMOL = LIST(JMOL)
         ENDDO
         LIST(JMOL) = IJ
        ELSE
         HEAD(ICN) = IJ
        ENDIF
        LIST(IJ) = 0
       X1(IJ) = XITR(J)
       Y1(IJ) = YITR(J)
       Z1(IJ) = ZITR(J)
      ENDDO
      RETURN
      END
C
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP = 1, N
       DO 11 IQ = 1, N
        V(IP,IQ) = 0.0D0
11     CONTINUE
       V(IP,IP) = 1.0D0
12    CONTINUE
      DO 13 IP = 1, N
       B(IP) = A(IP,IP)
       D(IP) = B(IP)
       Z(IP) = 0.0D0
13    CONTINUE
      NROT = 0
      DO 24 I = 1, 50
       SM = 0.0D0
       DO 15 IP = 1, N-1
        DO 14 IQ = IP + 1, N
         SM = SM + DABS(A(IP,IQ))
14      CONTINUE
15     CONTINUE
       IF(SM.EQ.0.0D0)RETURN
       IF(I.LT.4)THEN
        TRESH = 0.2D0*SM/N**2
       ELSE
        TRESH = 0.0D0
       ENDIF
       DO 22 IP = 1, N - 1
        DO 21 IQ = IP + 1, N
         G = 100.0D0*DABS(A(IP,IQ))
         IF((I.GT.4).AND.(DABS(D(IP))+G.EQ.DABS(D(IP)))
     C         .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ))))THEN
          A(IP,IQ) = 0.0D0
         ELSE IF(DABS(A(IP,IQ)).GT.TRESH)THEN
          H = D(IQ) - D(IP)
          IF(DABS(H)+G.EQ.DABS(H))THEN
           T = A(IP,IQ)/H
          ELSE
           THETA = 0.50D0*H/A(IP,IQ)
           T = 1.0D0/(DABS(THETA)+DSQRT(1.0D0+THETA**2))
           IF(THETA.LT.0.0D0)T = -T
          ENDIF
          C = 1./DSQRT(1.0D0+T**2)
          S = T*C
          TAU = S/(1.0D0+C)
          H = T*A(IP,IQ)
          Z(IP) = Z(IP) - H
          Z(IQ) = Z(IQ) + H
          D(IP) = D(IP) - H
          D(IQ) = D(IQ) + H
          A(IP,IQ)=0.0D0
          DO 16 J = 1, IP - 1
           G = A(J,IP)
           H = A(J,IQ)
           A(J,IP) = G - S*(H+G*TAU)
           A(J,IQ) = H + S*(G-H*TAU)
16        CONTINUE
          DO 17 J = IP + 1, IQ - 1
           G = A(IP,J)
           H = A(J,IQ)
           A(IP,J) = G - S*(H+G*TAU)
           A(J,IQ) = H + S*(G-H*TAU)
17        CONTINUE
          DO 18 J = IQ + 1, N
           G = A(IP,J)
           H = A(IQ,J)
           A(IP,J) = G - S*(H+G*TAU)
           A(IQ,J) = H + S*(G-H*TAU)
18        CONTINUE
          DO 19 J = 1, N
           G = V(J,IP)
           H = V(J,IQ)
           V(J,IP) = G - S*(H+G*TAU)
           V(J,IQ) = H + S*(G-H*TAU)
19        CONTINUE
          NROT = NROT + 1
         ENDIF
21      CONTINUE
22     CONTINUE
       DO 23 IP = 1, N
        B(IP) = B(IP) + Z(IP)
        D(IP) = B(IP)
        Z(IP) = 0.0D0
23     CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END
