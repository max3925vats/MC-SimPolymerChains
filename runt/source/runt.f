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
