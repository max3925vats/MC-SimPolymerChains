!---------------------------------------------------------!
!     MONTE CARLO SIMULATION OF HARD CORE + FULL YUKAWA
!     CHAINS AROUND A BIG BEAD
!---------------------------------------------------------!
      PROGRAM BBEADMC
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      INTEGER HEAD
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (NBMAX=10000,NMAX=40,MCMAX=25,MAXBIN=400)
      DIMENSION X(NBMAX),Y(NBMAX),Z(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION NR(MAXBIN),PVOL(MAXBIN)
      COMMON /POS/ X,Y,Z,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /RON/ DLR, DINT
      COMMON /BOX/ AL,BDIA
      COMMON /INTVAR/ NMOL,N,NVL
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /BIGBEAD/ BBEPS,RSP,RSPL2
      COMMON /SIGS/ DEFF,CEFF

!     GENERATE RANDOM NUMBER SEED
      NSEED = 2*INT(SECNDS(0.0)) + 2937

!     READ RUN PARAMETERS

      OPEN(UNIT=1,FILE='runt.inp',STATUS='OLD')
!     TOTAL MOVES
      READ(1,*)NCON
!     MOVES PER ACCUMULATION
      READ(1,*)NSKIP
!     BIN SIZE
      READ(1,*)BSZ
!     CHAIN TRANSLATION
      READ(1,*)DLR
!     INTERNAL DISPLACEMENT
      READ(1,*)DINT
!     FRACTION OF DICKMAN MOVES
      READ(1,*)FDICK
!     FRACTION OF REPTATION
      READ(1,*)FREPT
!     INVERSE TEMPERATURE
      READ(1,*)BEPS
!     BIG BEAD ATTRACTION
      READ(1,*)BBEPS
!     CENTRAL SPHERE DIAMETER
      READ(1,*)RSP
      CLOSE(UNIT=1,STATUS='KEEP')

      OPEN(UNIT=2,FILE='runt.ic',STATUS='OLD')
      READ(2,*)PFC
      READ(2,*)
      READ(2,*)NMOL
      READ(2,*)N
      READ(2,*)
      READ(2,*)
      READ(2,*)AL

      NBTOT = NMOL*N
      DO I = 1, NBTOT
         READ(2,*)ID,JD,X(I),Y(I),Z(I)
         X(I) = X(I)/AL
         Y(I) = Y(I)/AL
         Z(I) = Z(I)/AL
      ENDDO

      CLOSE(UNIT=2,STATUS='KEEP')

      BDIA = (AL*SQRT(3.0))
      BDIA02 = BDIA/2.0D0

      NBIN = INT(BDIA02/BSZ)
      IF(NBIN.GT.MAXBIN)WRITE(*,*) 'Bin dimension'

      RSPL2 = ((RSP+0.5D0)/AL)**2
      AL1 = AL/AL
      AL02 = AL1/2.0D0
      BSZR = BSZ/AL
      BDIA = BDIA/AL
      BDIA02 = BDIA02/AL
      BDIA2 = BDIA02*BDIA02
      DLR = DLR/AL
      DINT = DINT/AL
      DEFF = (1.0D0/AL)**2
      CEFF = 1.0D0**2 + (1.0D0/2.0D0)**2 + (1.0D0/2.0D0)**2

!     INITIALIZE BINS
      DO I = 1, NBIN
         NR(I) = 0
      ENDDO
!      DRSP = RSP
!      DO I = 1, NBIN
!         PVOL(I) = (4.0D0/3.0D0)*PI*((DRSP+BSZ)**3 - DRSP**3)
!         DRSP = DRSP + BSZ
!      ENDDO

!     INITIALIZE COUNTERS
      NTD = 0
      NTJ = 0
      NTR = 0
      NSD = 0
      NSJ = 0
      NSR = 0
      NAVER = 0

      AVENER = 0.0D0

!     CALCULATE INITIAL ENERGY
      CALL ENERCAL(ENERGY)
      WRITE(6,*)'INITIAL ENERGY',ENERGY

!     BEGIN SIMULATION
      DO 1000 ICON = 1, NCON
 
         IF (ICON.EQ.1)WRITE(6,*)'Simulation begun'
         IF (ICON.EQ.NCON)WRITE(6,*)'Ending simulation'
!      MAKE A MOVE

         I = INT( NMOL*RAN(NSEED) ) + 1
         IMOL = (I-1)*N

!      PICK DICKMAN, REPTATION OR CCB

         XRAN = RAN(NSEED)
         IF(XRAN.LT.FDICK)THEN
            NTD = NTD + 1
            CALL DICK(I,ISUC)
            IF(ISUC.EQ.0)NSD = NSD + 1
         ELSEIF(XRAN.LT.FDICK+FREPT)THEN
            NTR = NTR + 1
            CALL REPT(I,ISUC)
            IF(ISUC.EQ.0)NSR = NSR + 1
         ELSE
            NTJ = NTJ + 1
            CALL CCB(I,ISUC)
            IF(ISUC.EQ.0)NSJ = NSJ + 1
         ENDIF
!
         IF(ISUC.EQ.0)THEN
!       Update Positions And List Arrays
            ENERGY = ENERGY + ENEW - EOLD
            DO 100 J = 1, N
               IJ = IMOL + J
               X(IJ) = XITR(J)
               Y(IJ) = YITR(J)
               Z(IJ) = ZITR(J)
100         CONTINUE
         ENDIF

         IF(MOD(ICON,NSKIP).EQ.0)THEN
!       Average Properties
            NAVER = NAVER + 1
 
            AVENER = AVENER + ENERGY
!       Density Profiles
            DO I = 1, NMOL
               DO J = 1, N
                  IJ = (I-1)*N + J
                  X1 = X(IJ)
                  Y1 = Y(IJ)
                  Z1 = Z(IJ)
                  DSX = X1 - DNINT(X1)
                  DSY = Y1 - DNINT(Y1)
                  DSZ = Z1 - DNINT(Z1)
                  SQRCP = DSQRT(DSX**2 + DSY**2 + DSZ**2)
                  K = 1 + INT(SQRCP/BSZR)
                  IF(K.LE.NBIN)THEN
                     NR(K) = NR(K) + 1
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
1000  CONTINUE

      DLR = DLR*AL
      DINT = DINT*AL

      WRITE(6,*)'Writing results'
!     Write Results
      IF(NTD.GT.0)FSD=DBLE(NSD)/DBLE(NTD)
      IF(NTJ.GT.0)FSJ=DBLE(NSJ)/DBLE(NTJ)
      IF(NTR.GT.0)FSR=DBLE(NSR)/DBLE(NTR)
      IF(NAVER.GT.0)AVENER = AVENER/DBLE(NAVER)/DBLE(NBTOT)
      WRITE(6,*)AVENER

      OPEN(UNIT=3,FILE='runt.fc',STATUS='UNKNOWN')
      WRITE(3,111)PFC,NMOL,N,AL
      DO I = 1, NMOL
         DO J = 1, N
            IJ = (I-1)*N + J
            X(IJ) = X(IJ)*AL
            Y(IJ) = Y(IJ)*AL
            Z(IJ) = Z(IJ)*AL
            WRITE(3,112)I,J,X(IJ),Y(IJ),Z(IJ)
         ENDDO
      ENDDO
      CLOSE(UNIT=3,STATUS='KEEP')
111   FORMAT(3X,F6.4,3X,'PACKING FRACTION' /
     C      2X, I3, 6X,'NUMBER OF MOLECULES' /
     C      2X, I3, 6X,'CHAIN LENGTH' /
     C      4X, E18.10, 2X, 'PERIODIC LENGTH ' )
!
112   FORMAT(1X,I3,2X,I3,2X,3E18.10)
!
      OPEN(UNIT=5,FILE='runt.out',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='denrunt.out',STATUS='UNKNOWN')
!
      VOL = AL**3
      PFC = (PI/6.0D0)*NMOL*N/VOL
      WRITE(5,*)'RESULTS OF HARD-YUKAWA CHAIN BIG BEAD RUN'
      WRITE(5,*)
      WRITE(5,113)PFC,N,NMOL,DLR,DINT
      WRITE (5,114) NCON, NSKIP
      WRITE (5,115) BSZ,NBIN
      WRITE (5,116) BEPS,BBEPS
      WRITE(5,*)'Move statistics'
      WRITE (5,117) NTD, NTR, NTJ, NSD, NSR, NSJ,
     C              FSD, FSR, FSJ
      WRITE (5,118) ALH,AL
      WRITE (5,*)'DENSITY PROFILES NEAR THE BIG BEAD'
      CONST = 4.0D0*PI/3.0D0
      DO I=1, NBIN
         RLOWER = DBLE(I-1)*BSZ
         RUPPER = RLOWER + BSZ
         VOL = CONST*(RUPPER**3 - RLOWER**3)*DBLE(NAVER)
         DIST = RLOWER + BSZ/2.0
         IF(NAVER.GT.0)ANR = DBLE(NR(I))/VOL
         IF(NAVER.GT.0)WRITE (5,1140)I,DIST,ANR
         IF(NAVER.GT.0)WRITE (4,1140)I,DIST,ANR
      ENDDO
      CLOSE(UNIT=5,STATUS='KEEP')
      CLOSE(UNIT=4,STATUS='KEEP')
113   FORMAT(1X,'PACKING FRACTION ',4X,F6.4 /
     C       1X,' NUMBER OF BEADS ',6X, I4  /
     C       1X,'NUMBER OF CHAINS ',6X, I4  /
     C       1X,'     DELTA CHAIN ',5X, F5.3, /
     C       1X,'  DELTA INTERNAL ',5X, F5.3 )
114   FORMAT (1X / 1X, 'NO. TOTAL MOVES', 2X, I8 /
     C             1X, 'NUMBER OF SKIPS', 3X, I7 /)
115   FORMAT (1X  '      BIN SIZE', 3X, F5.3/
     C        1X, 'NUMBER OF BINS', 4X,I4 / )
116   FORMAT (1X  '      BEPS', 3X, F10.4/
     C        1X, '     BBEPS', 3X, F10.4 / )
118   FORMAT (1X  'BOX LENGTH', 3X, E16.8 /
     C        1X, 'PERIODIC LENGTH' 4X, E16.8 / )
117   FORMAT (1X / 13X, 'DICKMAN  REPTATION       CCB' /
     C        1X, 'ATTEMPTED', 3X, I7, 4X, I7, 3X, I7 /
     C        1X, 'SUCCESSFUL', 2X, I7, 4X, I7, 3X, I7 /
     C        1X, 'FRACTION',5X, F6.4, 5X, F6.4, 4X, F6.4 / )
1140  FORMAT (1X,I4,1X,F6.3,3X,3(F9.5,1X,F9.5,2X))
!
      WRITE(6,*)'files created'
      STOP
      END
!
!
!
      SUBROUTINE DICK(I,ISUC)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (NBMAX=10000,NMAX=40)
      DIMENSION X(NBMAX),Y(NBMAX),Z(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      COMMON /POS/ X,Y,Z,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /RON/ DLR, DINT
      COMMON /BOX/ AL,BDIA
      COMMON /INTVAR/ NMOL,N,NVL
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /SIGS/ DEFF,CEFF
      COMMON /BIGBEAD/ BBEPS,RSP,RSPL2
      ISUC = 1
      Z1 = RAN(NSEED) - 0.5
      Z2 = RAN(NSEED) - 0.5
      Z3 = RAN(NSEED) - 0.5
      ENEW = 0
      EOLD = 0
      IMOL = (I-1)*N
      IBD = IMOL + 1
!     DISPLACE FIRST SEGMENT
      XITR(1) = X(IBD) + Z1*DLR
      YITR(1) = Y(IBD) + Z2*DLR
      ZITR(1) = Z(IBD) + Z3*DLR
!     CHECK FOR OVERLAP WITH CENTRAL SPHERE
      XATT = XITR(1)
      YATT = YITR(1)
      ZATT = ZITR(1)
      DSX = XATT - DNINT(XATT)
      DSY = YATT - DNINT(YATT)
      DSZ = ZATT - DNINT(ZATT)
      RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
      IF (RCP.LT.RSPL2)RETURN

!     CHECK FOR INTERMOLECULAR OVERLAP
      CALL NEWEN1(I,1)
      IF(NVL.EQ.1)RETURN
      ENEW = ENEW + ENER + BBATT(RCP)

!     DISPLACE SUBSEQUENT SEGMENTS
      DO 200 II = 2, N
         IBD = IMOL + II
         D1 = X(IBD)-X(IBD-1)
         D2 = Y(IBD)-Y(IBD-1)
         D3 = Z(IBD)-Z(IBD-1)
         IF(RAN(NSEED).GT.0.25)GOTO 211
         CALL RUV(Z1,Z2,Z3)
         D1 = D1 + DINT*Z1
         D2 = D2 + DINT*Z2
         D3 = D3 + DINT*Z3
         CNM = D1*D1 + D2*D2 + D3*D3
         CN = 1.0D0/DSQRT(CNM)/AL
         D1 = D1*CN
         D2 = D2*CN
         D3 = D3*CN
211   CONTINUE
      XITR(II) = XITR(II-1) + D1
      YITR(II) = YITR(II-1) + D2
      ZITR(II) = ZITR(II-1) + D3
!    CHECK FOR OVERLAP WITH CENTRAL SPHERE
      XATT = XITR(II)
      YATT = YITR(II)
      ZATT = ZITR(II)
      DSX = XATT - DNINT(XATT)
      DSY = YATT - DNINT(YATT)
      DSZ = ZATT - DNINT(ZATT)
      RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
      IF (RCP.LT.RSPL2)RETURN
!    CHECK FOR INTERMOLECULAR OVERLAP
      CALL NEWEN1(I,II)
      IF(NVL.EQ.1)RETURN
      ENEW = ENEW + ENER + BBATT(RCP)
!     CHECK FOR INTRA OVERLAP
      IF(II.GT.2)THEN
         DO 215 JJ = 1, II - 2
            DX = XITR(JJ)-XITR(II)
            DY = YITR(JJ)-YITR(II)
            DZ = ZITR(JJ)-ZITR(II)
            DX = DX - DNINT(DX)
            DY = DY - DNINT(DY)
            DZ = DZ - DNINT(DZ)
            RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
            IF( RDIS .LT. DEFF ) RETURN
            IF(RDIS.LE.CEFF)ENEW = ENEW + BATT(RDIS)
215      CONTINUE
      ENDIF
200   CONTINUE
      EOLD = 0.0D0
      DO J = 3, N
         DO K = 1, J - 2
            I1 = IMOL + J
            I2 = IMOL + K
            DX = X(I1)-X(I2)
            DY = Y(I1)-Y(I2)
            DZ = Z(I1)-Z(I2)
            DX = DX - DNINT(DX)
            DY = DY - DNINT(DY)
            DZ = DZ - DNINT(DZ)
            RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
            EOLD = EOLD + BATT(RDIS)
         ENDDO
      ENDDO
      DO J = 1, N
         IBD = IMOL + J
         CALL OLDEN1(I,J)
         XATT = X(IBD)
         YATT = Y(IBD)
         ZATT = Z(IBD)
         DSX = XATT - DNINT(XATT)
         DSY = YATT - DNINT(YATT)
         DSZ = ZATT - DNINT(ZATT)
         RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
         EOLD = EOLD + ENER + BBATT(RCP)
      ENDDO
      IF(EOLD.LE.ENEW)THEN
         BOLTZ = DEXP(EOLD-ENEW)
         IF(BOLTZ.LE.RAN(NSEED))RETURN
      ENDIF
!     MOVE ACCEPTED
      ISUC = 0
      RETURN
      END
!
      SUBROUTINE REPT(I,ISUC)
!     PERFORMS A REPTATION MOVE ON LINEAR CHAIN I
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER (NBMAX=10000,NMAX=40)
      DIMENSION X(NBMAX),Y(NBMAX),Z(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      COMMON /POS/ X,Y,Z,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /BOX/ AL,BDIA
      COMMON /INTVAR/ NMOL,N,NVL
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /BIGBEAD/ BBEPS,RSP,RSPL2
      COMMON /SIGS/ DEFF,CEFF

      ISUC = 1
      IMOL = (I-1)*N
      ENEW = 0.0D0

!     CUT OFF END AND ATTACH TO BEAD 1
      CALL RUV(DELX,DELY,DELZ)
      XITR(1) = X(IMOL+1) + DELX /AL
      YITR(1) = Y(IMOL+1) + DELY /AL
      ZITR(1) = Z(IMOL+1) + DELZ /AL
!     CHECK FOR OVERLAP WITH CENTRAL SPHERE
      XATT = XITR(1)
      YATT = YITR(1)
      ZATT = ZITR(1)
      DSX = XATT - DNINT(XATT)
      DSY = YATT - DNINT(YATT)
      DSZ = ZATT - DNINT(ZATT)
      RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
      IF (RCP.LT.RSPL2)RETURN
!     CHECK FOR INTER CHAIN OVERLAP
      CALL NEWEN1(I,1)
      IF(NVL.EQ.1)RETURN
      ENEW = ENEW + ENER + BBATT(RCP)
!     CHECK FOR INTRA CHAIN OVERLAP
      DO 320 II = 2, N-1
         I1 = IMOL + II
         DX = XITR(1)-X(I1)
         DY = YITR(1)-Y(I1)
         DZ = ZITR(1)-Z(I1)
         DX = DX - DNINT(DX)
         DY = DY - DNINT(DY)
         DZ = DZ - DNINT(DZ)
         RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
         IF(RDIS.LT.DEFF)RETURN
         IF(RDIS.LE.CEFF)ENEW = ENEW + BATT(RDIS)
320   CONTINUE
!     CALCULATE OLD ENERGY
      EOLD = 0.0D0
      I2 = IMOL+N
      DO II = 1, N-2
         I1 = IMOL + II
         DX = X(I1)-X(I2)
         DY = Y(I1)-Y(I2)
         DZ = Z(I1)-Z(I2)
         DX = DX - DNINT(DX)
         DY = DY - DNINT(DY)
         DZ = DZ - DNINT(DZ)
         RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
         EOLD = EOLD + BATT(RDIS)
      ENDDO
      CALL OLDEN1(I,N)
      XATT = X(I2)
      YATT = Y(I2)
      ZATT = Z(I2)
      DSX = XATT - DNINT(XATT)
      DSY = YATT - DNINT(YATT)
      DSZ = ZATT - DNINT(ZATT)
      RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
      EOLD = EOLD + ENER + BBATT(RCP)
      IF(EOLD.LE.ENEW)THEN
         BOLTZ = DEXP(EOLD-ENEW)
         IF(BOLTZ.LE.RAN(NSEED))RETURN
      ENDIF
!     MOVE ACCEPTED
!     RESET NEW CHAIN BEADS
      DO J = 2, N
         I1 = IMOL + J-1
         XITR(J) = X(I1)
         YITR(J) = Y(I1)
         ZITR(J) = Z(I1)
      ENDDO
      ISUC = 0
      RETURN
      END
!
      SUBROUTINE CCB(I,ISUC)
!     PERFORMS A CCB MOVE ON CHAIN I
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RAN
      PARAMETER (PI=3.141592653589793D0)
      PARAMETER ( NSAMP=15 )
      PARAMETER (NBMAX=10000,NMAX=40)
      DIMENSION X(NBMAX),Y(NBMAX),Z(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      DIMENSION XT(NSAMP),YT(NSAMP),ZT(NSAMP),ET(NSAMP),
     C  WS(NMAX),STX(NMAX),STY(NMAX),STZ(NMAX),CO(NMAX),
     C  XN(NMAX),YN(NMAX),ZN(NMAX)
      COMMON /POS/ X,Y,Z,XITR,YITR,ZITR
      COMMON /SEED/ NSEED
      COMMON /BOX/ AL,BDIA
      COMMON /INTVAR/ NMOL,N,NVL
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /BIGBEAD/ BBEPS,RSP,RSPL2
      COMMON /SIGS/ DEFF,CEFF
      ISUC = 1
      WN = 1.0D0
      WO = 1.0D0

      IMOL = (I-1)*N
      DO J = 1, N
         I1 = IMOL + J
         XN(J) = X(I1)
         YN(J) = Y(I1)
         ZN(J) = Z(I1)
      ENDDO
      ICUT = INT((N-1)*RAN(NSEED)) + 2

      DO 50 J = ICUT, N
         SUM = 0.0D0
         DO 20 K = 1, NSAMP
            ENI = 0.0D0
            ET(K) = 1.0D0
            CALL RUV(DELX,DELY,DELZ)
            XT(K) = XN(J-1) + DELX/AL
            YT(K) = YN(J-1) + DELY/AL
            ZT(K) = ZN(J-1) + DELZ/AL
            XITR(J) = XT(K)
            YITR(J) = YT(K)
            ZITR(J) = ZT(K)
!     CHECK FOR OVERLAP WITH CENTRAL SPHERE
            XATT = XITR(J)
            YATT = YITR(J)
            ZATT = ZITR(J)
            DSX = XATT - DNINT(XATT)
            DSY = YATT - DNINT(YATT)
            DSZ = ZATT - DNINT(ZATT)
            RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
            IF (RCP.LT.RSPL2)RETURN
!     CHECK FOR INTER CHAIN OVERLAP
            CALL NEWEN1(I,J)
            IF(NVL.EQ.1)THEN
               ET(K) = 0.0D0
               GOTO 19
            ENDIF
            ENI = ENI + ENER + BBATT(RCP)
!     CHECK FOR INTRA CHAIN OVERLAP
            IF(J.GE.3)THEN
               DO 16 II = 1, J - 2
                  DX = XT(K)-XN(II)
                  DY = YT(K)-YN(II)
                  DZ = ZT(K)-ZN(II)
                  DX = DX - DNINT(DX)
                  DY = DY - DNINT(DY)
                  DZ = DZ - DNINT(DZ)
                  RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
                  IF(RDIS.LT.DEFF)THEN
                     ET(K) = 0.0D0
                     GOTO 19
                  ENDIF
                  IF(RDIS.LE.CEFF)ENI = ENI + BATT(RDIS)
16             CONTINUE
            ENDIF
            ET(K) = DEXP(-ENI)
19          CONTINUE
            SUM = SUM + ET(K)
20       CONTINUE
         IF(SUM.LT.1.0D-10)RETURN
         SUM = 1.0D0/SUM
         DO 25 K = 1, NSAMP
            ET(K) = ET(K)*SUM
25       CONTINUE
         XRAN = RAN(NSEED)
         S = 0.0D0
         DO 30 K = 1, NSAMP
            S = S + ET(K)
            IF(XRAN.LT.S)GOTO 35
30       CONTINUE
35       CONTINUE
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
!  COMPUTE WEIGHT OF EXISTING CHAIN
      DO 65 J = 1, N
         I1 = IMOL + J
         XN(J) = X(I1)
         YN(J) = Y(I1)
         ZN(J) = Z(I1)
65    CONTINUE

      DO 100 J = ICUT, N
         SUM = 0.0D0
         DO 80 K = 1, NSAMP
            ENO = 0.0D0
            IF(K.EQ.1)THEN
               CALL OLDEN1(I,J)
               I1 = IMOL + J
               XATT = X(I1)
               YATT = Y(I1)
               ZATT = Z(I1)
               DSX = XATT - DNINT(XATT)
               DSY = YATT - DNINT(YATT)
               DSZ = ZATT - DNINT(ZATT)
               RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
               ENO = ENO + ENER + BBATT(RCP)
               IF(J.GE.3)THEN
                  DO II = 1, J - 2
                     I1 = IMOL + J
                     I2 = IMOL + II
                     DX = X(I1)-X(I2)
                     DY = Y(I1)-Y(I2)
                     DZ = Z(I1)-Z(I2)
                     DX = DX - DNINT(DX)
                     DY = DY - DNINT(DY)
                     DZ = DZ - DNINT(DZ)
                     RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
                     IF(RDIS.LT.DEFF)WRITE(6,*)'ERROR1 CHECK OLD CONFIG'
                     IF(RDIS.LE.CEFF)ENO = ENO + BATT(RDIS)
                  ENDDO
               ENDIF
               ET(K) = DEXP(-ENO)
               GOTO 79
            ENDIF
            CALL RUV(DELX,DELY,DELZ)
            XT(K) = XN(J-1) + DELX /AL
            YT(K) = YN(J-1) + DELY /AL
            ZT(K) = ZN(J-1) + DELZ /AL
            XITR(J) = XT(K)
            YITR(J) = YT(K)
            ZITR(J) = ZT(K)
!     CHECK FOR OVERLAP WITH CENTRAL SPHERE
            XATT = XITR(J)
            YATT = YITR(J)
            ZATT = ZITR(J)
            DSX = XATT - DNINT(XATT)
            DSY = YATT - DNINT(YATT)
            DSZ = ZATT - DNINT(ZATT)
            RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
            IF (RCP.LT.RSPL2)RETURN
!      CHECK FOR INTER CHAIN OVERLAP
            CALL NEWEN1(I,J)
            IF(NVL.EQ.1)THEN
               ET(K) = 0.0D0
               GOTO 79
            ENDIF
            ENO = ENO + ENER + BBATT(RCP)
!     CHECK FOR INTRA CHAIN OVERLAP
            IF(J.GE.3)THEN
               DO 75 II = 1, J - 2
                  DX = XT(K)-XN(II)
                  DY = YT(K)-YN(II)
                  DZ = ZT(K)-ZN(II)
                  DX = DX - DNINT(DX)
                  DY = DY - DNINT(DY)
                  DZ = DZ - DNINT(DZ)
                  RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
                  IF( RDIS .LT. DEFF )THEN
                     ET(K) = 0.0D0
                     GOTO 79
                  ENDIF
                  IF(RDIS.LE.CEFF)ENO = ENO + BATT(RDIS)
75             CONTINUE
            ENDIF
            ET(K) = DEXP(-ENO)
79          CONTINUE
            SUM = SUM + ET(K)
80       CONTINUE
         ET(1) = ET(1)/SUM
         WO = WO*ET(1)
100   CONTINUE
!     CALCULATE AND VERIFY ENERGIES
      DO J = 1, N
         XITR(J) = STX(J)
         YITR(J) = STY(J)
         ZITR(J) = STZ(J)
      ENDDO
      EOLD = 0.0D0
      ENEW = 0.0D0
      DO J = 1, N
         CALL NEWEN1(I,J)
         IF(NVL.EQ.1)WRITE(6,*)'Screwed in Linear CCB'
         XATT = XITR(J)
         YATT = YITR(J)
         ZATT = ZITR(J)
         DSX = XATT - DNINT(XATT)
         DSY = YATT - DNINT(YATT)
         DSZ = ZATT - DNINT(ZATT)
         RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
         ENEW = ENEW + ENER + BBATT(RCP)
         CALL OLDEN1(I,J)
         IMOLJ = IMOL + J
         XATT = X(IMOLJ)
         YATT = Y(IMOLJ)
         ZATT = Z(IMOLJ)
         DSX = XATT - DNINT(XATT)
         DSY = YATT - DNINT(YATT)
         DSZ = ZATT - DNINT(ZATT)
         RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8

         EOLD = EOLD + ENER + BBATT(RCP)
         IF(J.GE.3)THEN
            DO K = 1, J - 2
               DX = XITR(J)-XITR(K)
               DY = YITR(J)-YITR(K)
               DZ = ZITR(J)-ZITR(K)
               DX = DX - DNINT(DX)
               DY = DY - DNINT(DY)
               DZ = DZ - DNINT(DZ)
               RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
               IF(RDIS.LE.CEFF)ENEW = ENEW + BATT(RDIS)
               IJ = IMOL + J
               IK = IMOL + K
               DX = X(IJ)-X(IK)
               DY = Y(IJ)-Y(IK)
               DZ = Z(IJ)-Z(IK)
               DX = DX - DNINT(DX)
               DY = DY - DNINT(DY)
               DZ = DZ - DNINT(DZ)
               RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
               IF(RDIS.LE.CEFF)EOLD = EOLD + BATT(RDIS)
            ENDDO
         ENDIF
      ENDDO
      BOLTZ = WO/WN*DEXP(EOLD-ENEW)
      IF(BOLTZ.LT.RAN(NSEED))RETURN
      ISUC = 0
      RETURN
      END
!
      SUBROUTINE RUV(X,Y,Z)
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
         X = 2.D0 * B1 * BH
         Y = 2.D0 * B2 * BH
         Z = 1.D0 - 2.D0*BSQ
      ENDIF
      RETURN
      END
!
      DOUBLE PRECISION FUNCTION BBATT(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /BIGBEAD/ BBEPS,RSP,RSPL2
      COMMON /BOX/ AL,BDIA
      S = X*AL**2
      RIJ = DSQRT(S)
      RSPSQ = DSQRT(RSPL2)
      IF (RIJ.LT.RSPSQ)BBATT = 0
      RETURN
      BBATT = BBEPS*(DEXP(-2.5D0*(RIJ-RSPSQ))/RIJ)
      RETURN
      END
!
      DOUBLE PRECISION FUNCTION BATT(X)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /BOX/ AL,BDIA
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /SIGS/ DEFF,CEFF
      S = X*AL**2
      RIJ = DSQRT(S)
      BATT = BEPS*(DEXP(-2.5D0*(RIJ-1.0D0))/RIJ)
      RETURN
      END
!
      SUBROUTINE ENERCAL(ENERGY)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NBMAX=10000,NMAX=40)
      DIMENSION X(NBMAX),Y(NBMAX),Z(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      COMMON /POS/ X,Y,Z,XITR,YITR,ZITR
      COMMON /BOX/ AL,BDIA
      COMMON /INTVAR/ NMOL,N,NVL
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /BIGBEAD/ BBEPS,RSP,RSPL2
      COMMON /SIGS/ DEFF,CEFF
!     INTERMOLECULAR ENERGY
      ENERGY = 0.0D0
      DO I = 1, NMOL
         DO J = 1, N
            CALL OLDEN1(I,J)
            ENERGY = ENERGY + ENER
         ENDDO
      ENDDO
      ENERGY = ENERGY*0.50D0
      DO IMOL = 1, NMOL
         DO J = 1, N
            IMOLJ = (IMOL-1)*N + J
            XATT = X(IMOLJ)
            YATT = Y(IMOLJ)
            ZATT = Z(IMOLJ)
            DSX = XATT - DNINT(XATT)
            DSY = YATT - DNINT(YATT)
            DSZ = ZATT - DNINT(ZATT)
            RCP=DSX**2 + DSY**2 + DSZ**2 + 1.D-8
            ENERGY = ENERGY + BBATT(RCP)
         ENDDO
      ENDDO
!     INTRAMOLECULAR ENERGY
      DO I = 1, NMOL
         IMOL = (I-1)*N
         DO J = 3, N
            DO K = 1, J - 2
               IJ = IMOL + J
               IK = IMOL + K
               DX = X(IK)-X(IJ)
               DY = Y(IJ)-Y(IK)
               DZ = Z(IJ)-Z(IK)
               DX = DX - DNINT(DX)
               DY = DY - DNINT(DY)
               DZ = DZ - DNINT(DZ)
               RDIS = DX*DX+DY*DY+DZ*DZ + 1.D-8
               ENERGY = ENERGY + BATT(RDIS)
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END
!
      SUBROUTINE OLDEN1(ITRY,JTRY)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NBMAX=10000,NMAX=40,MCMAX=25)
      DIMENSION X(NBMAX),Y(NBMAX),Z(NBMAX),JNEAR(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      COMMON /POS/ X,Y,Z,XITR,YITR,ZITR
      COMMON /INTVAR/ NMOL,N,NVL
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /SIGS/ DEFF,CEFF
      COMMON /BOX/ AL,BDIA

      NVL = 1
      ENER = 0.0D0
      IMN = (ITRY-1)*N + JTRY
      X1 = X(IMN)
      Y1 = Y(IMN)
      Z1 = Z(IMN)

      DO K = 1, NMOL
         DO J = 1, N
            KJ = (K-1)*N + J
            IF(K.NE.ITRY) THEN
               DX = X(KJ) - X1
               DY = Y(KJ) - Y1
               DZ = Z(KJ) - Z1
               DX = DX - DNINT(DX)
               DY = DY - DNINT(DY)
               DZ = DZ - DNINT(DZ)
               RDIS = DX*DX + DY*DY + DZ*DZ + 1.D-8
               ENER = ENER + BATT(RDIS)
           ENDIF
         ENDDO
      ENDDO
      NVL = 0
      RETURN
      END
!
      SUBROUTINE NEWEN1(ITRY,JTRY)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (NBMAX=10000,NMAX=40,MCMAX=25)
      DIMENSION X(NBMAX),Y(NBMAX),Z(NBMAX),JNEAR(NBMAX)
      DIMENSION XITR(NMAX),YITR(NMAX),ZITR(NMAX)
      COMMON /POS/ X,Y,Z,XITR,YITR,ZITR
      COMMON /INTVAR/ NMOL,N,NVL
      COMMON /ENS/ BEPS,ENER,EOLD,ENEW
      COMMON /SIGS/ DEFF,CEFF
      COMMON /BOX/ AL,BDIA

      NVL = 1
      ENER = 0.0D0
      X1 = XITR(JTRY)
      Y1 = YITR(JTRY)
      Z1 = ZITR(JTRY)

      DO K = 1, NMOL
         DO J = 1, N
            KJ = (K-1)*N + J
            IF(K.NE.ITRY) THEN
               DX = X(KJ) - X1
               DY = Y(KJ) - Y1
               DZ = Z(KJ) - Z1
               DX = DX - DNINT(DX)
               DY = DY - DNINT(DY)
               DZ = DZ - DNINT(DZ)
               RDIS = DX*DX + DY*DY + DZ*DZ + 1.D-8
               IF(RDIS.LT.DEFF)RETURN
               ENER = ENER + BATT(RDIS)
            ENDIF
         ENDDO
      ENDDO
      NVL = 0
      RETURN
      END
