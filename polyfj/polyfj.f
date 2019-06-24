!     MONTE CARLO SIMULATION OF HARD CHAIN
!     NEAR A LARGER HARD SPHERE
      PROGRAM POLYHC   
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
!
!     GENERATE RANDOM NUMBER SEED 
      NSEED = 2*INT(SECNDS(0.0)) + 2937

!     READ RUN PARAMETERS 
 
      OPEN(UNIT=1,FILE='polyfj.inp',STATUS='OLD')
!     TOTAL MOVES 
      READ(1,*)NCON
!     MOVES PER ACCUMULATION 
      READ(1,*)NSKIP
!     BIN SIZE 
      READ(1,*)BSZ
      BSZCM = 0.10D0
!     CHAIN TRANSLATION 
      READ(1,*)DLR
!     INTERNAL DISPLACEMENT 
      READ(1,*)DINT
!     FRACTION OF DICKMAN MOVES
      READ(1,*)FDICK
!     FRACTION OF REPTATION 
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
!
      IF(MOD(N,2).EQ.0)THEN
         NM1 = N/2
         NM2 = NM1 + 1
      ELSE
         NM1 = (N-1)/2 + 1
         NM2 = NM1
      ENDIF
!
!     INITIALIZE BINS
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

!     PARAMETERS
      BSZR = BSZ/AL
      BSZC = BSZCM/AL
      DLR = DLR/AL
      DINT = DINT/AL
      DEFF = (1.0D0/AL)**2
      RSPL2 = ((RSP+0.5D0)/AL)**2

!     LINKED LIST
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
!
      DO I = 1, NCELL
         HEAD(I) = 0
      ENDDO
      DO I = 1, NBMAX
         LIST(I) = 0
      ENDDO

!     MAKE A LIST
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

!     INITIALIZE COUNTERS 
      NTD = 0
      NTJ = 0
      NTR = 0
      NSD = 0
      NSJ = 0
      NSR = 0
      NAVER = 0
!
!     BEGIN SIMULATION
      DO 1000 ICON = 1, NCON

!     MAKE A MOVE
      
         I = INT( NMOL1*RAN(NSEED) ) + 1
         IMOL = (I-1)*N
         XRAN = RAN(NSEED)
!      write(*,*)'xran',xran
!     PICK DICKMAN REPTATION OR CCB
         IF(XRAN.LT.FREPT)THEN 
            NTR = NTR + 1
            CALL REPT(I,ISUC)
!       write(*,*)'rept',isuc
            IF(ISUC.EQ.0)NSR = NSR + 1
!       IF(ISUC.EQ.0)write(*,*)I,ISUC,XRAN,'REPT'
         ELSEIF(XRAN.LT.FDICK+FREPT)THEN
            NTD = NTD + 1
            CALL DICK(I,ISUC)
!       write(*,*)'dick',isuc
            IF(ISUC.EQ.0)NSD = NSD + 1
!       IF(ISUC.EQ.0)write(*,*)I,ISUC,XRAN,'DICKMANN'
         ELSE
            NTJ = NTJ + 1
            CALL CCB(I,ISUC)
!       write(*,*)'ccb',isuc
            IF(ISUC.EQ.0)NSJ = NSJ + 1
!       IF(ISUC.EQ.0)write(*,*)I,ISUC,XRAN,'CCB'
         ENDIF
!
         IF(ISUC.EQ.0)THEN
!     UPDATE POSITIONS AND LIST ARRAYS
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

!       REMOVE X1(IJ) FROM LIST
            JMOL = HEAD(ICO)
            IF(JMOL.NE.IJ)THEN
               DO WHILE (LIST(JMOL).NE.IJ)
                  JMOL = LIST(JMOL)
               ENDDO
               LIST(JMOL) = LIST(IJ)
            ELSE
               HEAD(ICO) = LIST(JMOL)
            ENDIF
!       ADD XITR(J) TO THE LIST
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
10       CONTINUE
      ENDIF
!
      IF(MOD(ICON,NSKIP).EQ.0)THEN
         NAVER = NAVER + 1
!      DENSITY PROFILES
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
!
               IF(J.GE.2)THEN
                  IJ1 = (I-1)*N + J
                  IJ2 = (I-1)*N + J-1
                  X11 = DSQRT(X1(IJ1)**2+Y1(IJ1)**2+Z1(IJ1)**2)
                  X22 = DSQRT(X1(IJ2)**2+Y1(IJ2)**2+Z1(IJ2)**2)
                  XM = DSQRT((X1(IJ2)-X1(IJ1))**2+(Y1(IJ2)-Y1(IJ1))**2
     C                      +(Z1(IJ2)-Z1(IJ1))**2)/2.D0
                  COSO = (X1(IJ1)*X1(IJ2)+Y1(IJ1)*Y1(IJ2)
     C                     +Z1(IJ1)*Z1(IJ2))/(X11*X22)
                  K = 1 + INT(XM/BSZR)
                  NR(K) = NR(K) + 1
                  SR(K) = SR(K) + COSO*AL*AL
               ENDIF
            ENDDO
         ENDDO
!
!      RG PROFILE
         DO I = 1, NMOL1
            XCM = 0.0D0
            YCM = 0.0D0
            ZCM = 0.0D0
            X2B = 0.0D0
            Y2B = 0.0D0
            Z2B = 0.0D0
            DO J = 1, N
               IJ = (I-1)*N + J
!         X1I = X1(IJ) - DNINT(X1(IJ))
!         Y1I = Y1(IJ) - DNINT(Y1(IJ))
!         Z1I = Z1(IJ) - DNINT(Z1(IJ))
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
!
            I1 = (I-1)*N + 1
            IN = (I-1)*N + N
            XE = X1(IN) - X1(I1)
            YE = Y1(IN) - Y1(I1)
            ZE = Z1(IN) - Z1(I1)
            RE1I = XE*XE
            RE2I = YE*YE
            RE3I = ZE*ZE
!
!       CALCULATE SEMI-AXIS LENGTHS
!       CALCULATE MOMENT OF INERTIA TENSOR
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
!       DIAGONALIZE IT
            CALL JACOBI(AI,3,3,DI,VI,NROT)
!       FIND THE SMALLEST EIGENVALUE
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
!
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

!     WRITE RESULTS
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
!
      OPEN(UNIT=2,FILE='polyfj.fc',STATUS='UNKNOWN')
      WRITE(2,111)PFC,NMOL1,N,RSP,AL
      DO I = 1, NMOL1
         DO J = 1, N
            IJ = (I-1)*N + J
            WRITE(2,112)I,J,X1(IJ),Y1(IJ),Z1(IJ)
         ENDDO
      ENDDO
      CLOSE(UNIT=2,STATUS='KEEP')
111   FORMAT(3X,F10.4,3X,'PACKING FRACTION' /
     C       /
     C      2X, I8, 6X,'NUMBER OF MOLECULES' /
     C      2X, I6, 6X,'CHAIN LENGTH' /
     C      /
     C      2X, F10.4,5X,'RADIUS OF SOLUTE'/
     C      4X, E18.10, 2X, 'PERIODIC LENGTH '/ )
112   FORMAT(1X,I6,2X,I6,2X,3E18.10)
!
      OPEN(UNIT=1,FILE='polyfj.out',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='polyden.out',STATUS='UNKNOWN')
!
      VOL = AL**3 - 4.D0*PI*RSP**3/3.D0
      PFCMC = (PI/6.)*NCTOT/VOL
!
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
!
!     SEGMENTAL ORDER PARAMETER
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
!
      OPEN(UNIT=4,FILE='rgfj.out',STATUS='UNKNOWN')
      OPEN(UNIT=2,FILE='g2fj.out',STATUS='UNKNOWN')
      WRITE(1,*)'CENTRE OF MASS AND R2G PROFILE'
      WRITE(1,*)
      WRITE(1,113)PFC,N,NMOL1,DLR,DINT,RSP
      WRITE(1,114) NCON, NSKIP
      WRITE(1,118) BSZ,BSZCM,NBIN,NBINC
      WRITE(1,*)'Move statistics'
      WRITE(1,117) NTD, NTR, NTJ, NSD, NSR, NSJ,
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

      SUBROUTINE RUV(A,B,C)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      COMMON /SEED/ NSEED

1     B1 = 1.0D0 - 2.0D0*RAN(NSEED)
      B2 = 1.0D0 - 2.0D0*RAN(NSEED)
      BSQ = B1*B1 + B2*B2
      IF(BSQ.GT.1.0D0) THEN
!      REJECT
         GOTO 1
      ELSE
         BH = DSQRT(1.D0 - BSQ)
         A = 2.D0 * B1 * BH
         B = 2.D0 * B2 * BH
         C = 1.D0 - 2.D0*BSQ
      ENDIF
      RETURN
      END
!
      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NMAX=100)
      DIMENSION A(NP,NP),D(NP),V(NP,NP),B(NMAX),Z(NMAX)
      DO 12 IP = 1, N
         DO 11 IQ = 1, N
            V(IP,IQ) = 0.0D0
11       CONTINUE
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
14          CONTINUE
15       CONTINUE
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
     C             .AND.(DABS(D(IQ))+G.EQ.DABS(D(IQ))))THEN
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
16                CONTINUE
                  DO 17 J = IP + 1, IQ - 1
                     G = A(IP,J)
                     H = A(J,IQ)
                     A(IP,J) = G - S*(H+G*TAU)
                     A(J,IQ) = H + S*(G-H*TAU)
17                CONTINUE
                  DO 18 J = IQ + 1, N
                     G = A(IP,J)
                     H = A(IQ,J)
                     A(IP,J) = G - S*(H+G*TAU)
                     A(IQ,J) = H + S*(G-H*TAU)
18                CONTINUE
                  DO 19 J = 1, N
                     G = V(J,IP)
                     H = V(J,IQ)
                     V(J,IP) = G - S*(H+G*TAU)
                     V(J,IQ) = H + S*(G-H*TAU)
19                CONTINUE
                  NROT = NROT + 1
               ENDIF
21          CONTINUE
22       CONTINUE
         DO 23 IP = 1, N
            B(IP) = B(IP) + Z(IP)
            D(IP) = B(IP)
            Z(IP) = 0.0D0
23       CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
      END
