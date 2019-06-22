!     Polyfj Move subroutines: Dickman, Reptation and CCB

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
!     DISPLACE FIRST SEGMENT
      XITR(1) = X1(IBD) + ZS1*DLR
      YITR(1) = Y1(IBD) + ZS2*DLR
      ZITR(1) = Z1(IBD) + ZS3*DLR
!     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
      XT = XITR(1) - DNINT(XITR(1))
      YT = YITR(1) - DNINT(YITR(1))
      ZT = ZITR(1) - DNINT(ZITR(1))
      RDIST = XT**2+YT**2+ZT**2
      IF(RDIST.LT.RSPL2)RETURN
!     CHECK FOR INTERMOLECULAR OVERLAP
      CALL OVER(I,1)
      IF(NVL.EQ.1)RETURN
!     DISPLACE SUBSEQUENT SEGMENTS
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
211      CONTINUE
         XITR(II) = XITR(II-1) + D1
         YITR(II) = YITR(II-1) + D2
         ZITR(II) = ZITR(II-1) + D3
!     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
         XT = XITR(II) - DNINT(XITR(II))
         YT = YITR(II) - DNINT(YITR(II))
         ZT = ZITR(II) - DNINT(ZITR(II))
         RDIST = XT**2+YT**2+ZT**2
         IF(RDIST.LT.RSPL2)RETURN
!      CHECK FOR INTERMOLECULAR OVERLAP
         CALL OVER(I,II)
         IF(NVL.EQ.1)RETURN
!      CHECK FOR INTRA OVERLAP
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
215         CONTINUE
         ENDIF
200   CONTINUE
!     MOVE ACCEPTED
      ISUC = 0
      RETURN
      END

      SUBROUTINE REPT(I,ISUC)
!     PERFORMS A REPTATION MOVE ON LINEAR CHAIN I
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
!
      ISUC = 1
      IF(RAN(NSEED).GT.0.5)CALL CVERT(I)
      IMOL = (I-1)*N
!
!     CUT OFF END AND ATTACH TO BEAD 1
      CALL RUV(DELX,DELY,DELZ)
      XITR(1) = X1(IMOL+1) + DELX /AL
      YITR(1) = Y1(IMOL+1) + DELY /AL
      ZITR(1) = Z1(IMOL+1) + DELZ /AL
!     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
      XT = XITR(1) - DNINT(XITR(1))
      YT = YITR(1) - DNINT(YITR(1))
      ZT = ZITR(1) - DNINT(ZITR(1))
      RDIST = XT**2+YT**2+ZT**2
      IF(RDIST.LT.RSPL2)RETURN
!     CHECK FOR INTER CHAIN OVERLAP
      CALL OVER(I,1)
      IF(NVL.EQ.1)RETURN
!     CHECK FOR INTRA CHAIN OVERLAP
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
!     MOVE ACCEPTED
!     RESET NEW CHAIN BEADS
      DO J = 2, N
         I1 = IMOL + J-1
         XITR(J) = X1(I1)
         YITR(J) = Y1(I1)
         ZITR(J) = Z1(I1)
      ENDDO
      ISUC = 0
      RETURN
      END

      SUBROUTINE CCB(I,ISUC)
!
!     PERFORMS A CCB MOVE ON CHAIN I
!
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
!
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
!      write(*,*)ICUT
      DO 50 J = ICUT, N
         SUM = 0.0D0
         DO 20 K = 1, NSAMP
            ET(K) = 1.0D0
            CALL RUV(DELX,DELY,DELZ)
            XT(K) = XN(J-1) + DELX/AL
            YT(K) = YN(J-1) + DELY/AL
            ZT(K) = ZN(J-1) + DELZ/AL
!     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
            XT1 = XT(K) - DNINT(XT(K))
            YT1 = YT(K) - DNINT(YT(K))
            ZT1 = ZT(K) - DNINT(ZT(K))
            RDIST = XT1**2+YT1**2+ZT1**2
            IF(RDIST.LT.RSPL2)THEN
               ET(K) = 0.0D0
               GOTO 19
            ENDIF
!      CHECK FOR INTER CHAIN OVERLAP
            XITR(J) = XT(K)
            YITR(J) = YT(K)
            ZITR(J) = ZT(K)
            CALL OVER(I,J)
            IF(NVL.EQ.1)THEN
               ET(K) = 0.0D0
               GOTO 19
            ENDIF
!     CHECK FOR INTRA CHAIN OVERLAP
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
16             CONTINUE
            ENDIF
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
!     CHECK FIRST FOR OVERLAP WITH LARGE SPHERE
            XT1 = XT(K) - DNINT(XT(K))
            YT1 = YT(K) - DNINT(YT(K))
            ZT1 = ZT(K) - DNINT(ZT(K))
            RDIST = XT1**2+YT1**2+ZT1**2
            IF(RDIST.LT.RSPL2)THEN
               ET(K) = 0.0D0
               GOTO 79
            ENDIF
!      CHECK FOR INTER CHAIN OVERLAP
            XITR(J) = XT(K)
            YITR(J) = YT(K)
            ZITR(J) = ZT(K)
            CALL OVER(I,J)
            IF(NVL.EQ.1)THEN
               ET(K) = 0.0D0
               GOTO 79
            ENDIF
!     CHECK FOR INTRA CHAIN OVERLAP
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
75             CONTINUE
            ENDIF
79       CONTINUE
         SUM = SUM + ET(K)
80       CONTINUE
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
      ENDIF  
      RETURN
      END

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
!
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
      ENDDO
      RETURN
      END

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
10          CONTINUE
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
!
      NVL = 0
      RETURN
      END
