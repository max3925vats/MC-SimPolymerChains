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
