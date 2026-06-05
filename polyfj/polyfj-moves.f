! =============================================================================
! FILE: polyfj-moves.f
! CONTENTS: MC move and overlap routines for an ATHERMAL hard-chain /
!           hard-sphere system.
!
! ACCEPTANCE MODEL: a move is accepted iff it creates NO hard overlap.
!   There is no energy, no Metropolis criterion -- purely steric exclusion.
!
! UNIT CONVENTION: all coordinates are stored in REDUCED units (coord / AL,
!   where AL is the box length). The bead diameter sigma = 1 in real units,
!   so the hard-core threshold in reduced units is DEFF = (1/AL)^2.
!   Minimum-image distances are computed with DNINT (round-to-nearest).
!   The central sphere contact condition uses RSPL2 = (R_sphere + 0.5)^2 /
!   AL^2 (sphere radius plus one bead radius, in reduced units).
!
! OVERLAP TEST IN USE: OVER (brute-force O(N*NMOL) scan) -- this is the
!   routine actually called by all move routines.
! DEAD CODE: OVER1 (linked-cell version) is NEVER called anywhere. See the
!   NOTE block above OVER1 for details.
!
! ROUTINE INDEX:
!   DICK  (I, ISUC) -- Dickman move (translate bead 1, regrow rest)
!   REPT  (I, ISUC) -- reptation (remove tail, grow new head)
!   CCB   (I, ISUC) -- configurational-bias regrowth from random cut
!   OVER  (ITRY, JTRY) -- brute-force inter-chain overlap test [LIVE]
!   CVERT (I)       -- reverse chain in place, update linked-cell list
!   OVER1 (ITRY, JTRY) -- cell-list overlap test [DEAD CODE -- never called]
! =============================================================================
!     Polyfj Move subroutines: Dickman, Reptation and CCB

! -----------------------------------------------------------------------------
! SUBROUTINE DICK(I, ISUC)
! PURPOSE: Dickman move on chain I.  Displaces bead 1 by a random step, then
!   regrowing beads 2..N using the existing bond directions with a 25% chance
!   to perturb each bond before renormalizing to unit length.
! ARGUMENTS:
!   I    (in)  -- chain index (1-based)
!   ISUC (out) -- success flag: 0 = accepted, 1 = rejected (overlap found)
! ALGORITHM:
!   1. Translate bead 1 by (ZS1,ZS2,ZS3)*DLR; reject on central-sphere or
!      inter-chain overlap.
!   2. For each subsequent bead: 25% chance to perturb the bond vector by
!      DINT*random_unit, then renormalize to unit length in reduced coords.
!   3. Place bead; reject on central-sphere, inter-chain (OVER), or
!      intra-chain overlap (all non-bonded pairs separated by >=2 beads).
!   4. If all N beads placed without overlap, set ISUC=0 (accepted).
! -----------------------------------------------------------------------------
      SUBROUTINE DICK(I,ISUC)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RANF
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
      ZS1 = RANF(NSEED) - .5
      ZS2 = RANF(NSEED) - .5
      ZS3 = RANF(NSEED) - .5
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
         IF(RANF(NSEED).GT.0.25)GOTO 211
!        25% chance: perturb bond direction by a random unit vector scaled by
!        DINT, then renormalize the result to unit length in reduced coords
         CALL RUV(ZS1,ZS2,ZS3)
         D1 = D1 + DINT*ZS1
         D2 = D2 + DINT*ZS2
         D3 = D3 + DINT*ZS3
         CNM = D1*D1 + D2*D2 + D3*D3
!        CN = 1/(|D|*AL): normalizes the bond to unit real-space length (=1/AL
!        in reduced units)
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

! -----------------------------------------------------------------------------
! SUBROUTINE REPT(I, ISUC)
! PURPOSE: Reptation move on chain I.  Removes the tail bead and grows a new
!   head bead at the other end using a uniformly random bond direction.
! ARGUMENTS:
!   I    (in)  -- chain index (1-based)
!   ISUC (out) -- success flag: 0 = accepted, 1 = rejected (overlap found)
! ALGORITHM:
!   1. With 50% probability, reverse the chain (CVERT) so either end can be
!      the growing head.
!   2. Grow a new bead at a random unit-vector step from the current head.
!   3. Reject on central-sphere, inter-chain, or intra-chain overlap.
!   4. On acceptance, shift all bead positions: old bead J becomes new bead
!      J+1, and the new head occupies position 1.
! -----------------------------------------------------------------------------
      SUBROUTINE REPT(I,ISUC)
!     PERFORMS A REPTATION MOVE ON LINEAR CHAIN I
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RANF
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
      IF(RANF(NSEED).GT.0.5)CALL CVERT(I)
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

! -----------------------------------------------------------------------------
! SUBROUTINE CCB(I, ISUC)
! PURPOSE: Configurational-bias (Rosenbluth) regrowth of the tail of chain I
!   starting from a random cut point ICUT.
! ARGUMENTS:
!   I    (in)  -- chain index (1-based)
!   ISUC (out) -- success flag: 0 = accepted, 1 = rejected
! ALGORITHM:
!   1. Optionally reverse the chain (50% chance) so either end can be the tail.
!   2. Choose a random cut point ICUT in [2, N]; beads 1..ICUT-1 are kept.
!   3. NEW CHAIN: for each bead J from ICUT to N, generate NSAMP=15 random
!      trial directions; compute Rosenbluth weight ET(K)=1 if trial K has no
!      overlap, 0 otherwise.  WN accumulates the product of normalised weights.
!   4. Select one trial direction with probability proportional to ET(K);
!      store that bead position.
!   5. OLD CHAIN: for each bead J from ICUT to N, treat the existing
!      position as trial K=1 and generate NSAMP-1 additional random trials;
!      compute WO from the old-chain Rosenbluth weights (ET(1)/SUM at each J).
!   6. Accept if WO/WN > uniform random (athermal Rosenbluth criterion).
! NOTE ON ROSENBLUTH WEIGHTS: ET(K) is 1 for a valid (non-overlapping)
!   direction and 0 for an overlapping one.  SUM = number of valid trials.
!   The selected bead's normalised weight WS(J) = 1/SUM if chosen; WN is the
!   product over all regrown beads, WO likewise for the old chain fragment.
! -----------------------------------------------------------------------------
      SUBROUTINE CCB(I,ISUC)
!
!     PERFORMS A CCB MOVE ON CHAIN I
!
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RANF
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
      IF(RANF(NSEED).GT.0.5)CALL CVERT(I)
      IMOL = (I-1)*N
      DO J = 1, N
         I1 = IMOL + J
         XN(J) = X1(I1)
         YN(J) = Y1(I1)
         ZN(J) = Z1(I1)
      ENDDO
      ICUT = INT((N-1)*RANF(NSEED)) + 2
!      write(*,*)ICUT
!     --- GROW NEW CHAIN from bead ICUT to N ---
      DO 50 J = ICUT, N
         SUM = 0.0D0
!        Generate NSAMP trial bond directions; ET(K)=1 if valid, 0 if overlap
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
!        If all NSAMP directions overlap, the chain cannot be grown; reject
         IF(SUM.LT.1.0D-10)RETURN
!        Normalise ET(K) so they sum to 1 (selection probabilities)
         SUM = 1.0D0/SUM
         DO 25 K = 1, NSAMP
            ET(K) = ET(K)*SUM
25       CONTINUE
!        Select one trial direction by cumulative-probability roulette wheel
         XRAN = RANF(NSEED)
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
!        Accumulate new-chain Rosenbluth weight (product of 1/SUM at each bead)
         WN = WN*WS(J)
50    CONTINUE
!     Store the newly grown chain fragment in STX/STY/STZ
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

!     --- COMPUTE OLD-CHAIN ROSENBLUTH WEIGHT (WO) ---
!     For each bead J, the existing position is K=1; NSAMP-1 additional random
!     trials are generated.  ET(1)/SUM gives the old bead's Rosenbluth factor.
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
!        ET(1) (the old bead's weight) is normalised by the sum over all trials
         ET(1) = ET(1)/SUM
         WO = WO*ET(1)
100   CONTINUE
!     Athermal Rosenbluth acceptance: accept if WO/WN > uniform random
      BOLTZ = WO/WN
      IF(BOLTZ.GT.RANF(NSEED))THEN
         DO J = 1, N
            XITR(J) = STX(J)
            YITR(J) = STY(J)
            ZITR(J) = STZ(J)
         ENDDO
         ISUC = 0
      ENDIF  
      RETURN
      END

! -----------------------------------------------------------------------------
! SUBROUTINE OVER(ITRY, JTRY)
! PURPOSE: Brute-force inter-chain overlap test.  This is the LIVE overlap
!   routine used by all move subroutines (DICK, REPT, CCB).
! ARGUMENTS:
!   ITRY (in)  -- index of the chain being moved (excluded from the check)
!   JTRY (in)  -- index of the trial bead within the moving chain (in XITR/
!                 YITR/ZITR) to be tested
! RETURNS (via COMMON /INTVAR/):
!   NVL = 1 if trial bead JTRY overlaps any bead of any chain other than ITRY
!   NVL = 0 if no overlap found
! ALGORITHM:
!   O(N * NMOL) scan: loop over all molecules and all beads; skip molecule
!   ITRY; compute minimum-image distance squared and compare to DEFF.
! -----------------------------------------------------------------------------
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

! -----------------------------------------------------------------------------
! SUBROUTINE CVERT(I)
! PURPOSE: Reverse chain I in place (bead 1 becomes bead N, etc.) AND update
!   the linked-cell list (HEAD/LIST) to reflect the new bead positions.
! ARGUMENTS:
!   I (in) -- chain index (1-based)
! ALGORITHM:
!   1. Copy chain I in reversed order into the trial arrays XITR/YITR/ZITR.
!   2. For each bead slot IJ (position in the global bead arrays):
!      a. Find the cell ICO of the old position X1(IJ) and remove IJ from
!         that cell's linked list.
!      b. Find the cell ICN of the new position XITR(J) and append IJ to
!         that cell's linked list.
!   3. Write the reversed coordinates back into X1/Y1/Z1.
! NOTE: Cell indices are computed by shifting coordinates to [0,1) via
!   minimum-image then multiplying by CELLI (cells per box length).
! -----------------------------------------------------------------------------
      SUBROUTINE CVERT(I)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*4 RANF
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
!     Copy chain I reversed into trial arrays: XITR(J) = X1(IMOL + N-J+1)
      DO J = 1, N
         K = N - J + 1
         IK = IMOL + K
         XITR(J) = X1(IK)
         YITR(J) = Y1(IK)
         ZITR(J) = Z1(IK)
      ENDDO
!     Update linked-cell list: for each bead slot IJ, remove the old position
!     and insert the reversed position into the appropriate cell
      DO J = 1, N
         IJ = IMOL + J
!        Map old position to cell index ICO (shift to [0,1) then scale)
         XO = X1(IJ) - DNINT(X1(IJ)) + 0.50D0
         YO = Y1(IJ) - DNINT(Y1(IJ)) + 0.50D0
         ZO = Z1(IJ) - DNINT(Z1(IJ)) + 0.50D0
!        Map new (reversed) position to cell index ICN
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

! =============================================================================
! NOTE: OVER1 IS DEAD CODE -- NEVER CALLED
! ----------------------------------------
! OVER1 is a linked-cell (HEAD/MAP/LIST) implementation of the inter-chain
! overlap test, intended as a faster O(1) alternative to the brute-force OVER.
! However, OVER1 is NEVER called anywhere in this program.  All move routines
! (DICK, REPT, CCB) call OVER instead.
!
! Additionally, OVER1 reads IDI(IM) to identify which chain a bead belongs to,
! but IDI is never populated anywhere in the codebase.  As written, IDI will
! always be zero (default COMMON initialisation), causing the chain-exclusion
! filter (IDI(IM).EQ.ITRY) to never match, which would incorrectly reject
! every move.
!
! This routine should be treated as vestigial/prototype code.  Do not call it
! without first populating IDI and verifying the MAP initialisation.
! =============================================================================
! -----------------------------------------------------------------------------
! SUBROUTINE OVER1(ITRY, JTRY)
! PURPOSE: Linked-cell inter-chain overlap test (VESTIGIAL -- see NOTE above).
!   Intended to replace OVER with an O(1) cell-list scan over the 27
!   neighbouring cells of trial bead JTRY.
! ARGUMENTS:
!   ITRY (in)  -- chain index to exclude from the overlap check
!   JTRY (in)  -- trial bead index (in XITR/YITR/ZITR)
! RETURNS (via COMMON /INTVAR/):
!   NVL = 1 if overlap found, NVL = 0 if clear
! STATUS: DEAD CODE.  Never called.  IDI (chain-membership array) is never
!   populated, so chain exclusion logic is broken.  Use OVER instead.
! -----------------------------------------------------------------------------
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
!
      REAL*4 FUNCTION RANF(ISEED)
!     PORTABLE UNIFORM (0,1) GENERATOR -- PARK-MILLER "MINIMAL STANDARD"
!     LCG WITH SCHRAGE'S METHOD, UPDATING ISEED IN PLACE. REPLACES THE
!     LEGACY RAN(NSEED) INTRINSIC, WHICH UNDER gfortran RETURNS A CONSTANT
!     AND NEVER ADVANCES THE SEED (BREAKING ALL MOVE SELECTION/SAMPLING).
      INTEGER ISEED, IA, IM, IQ, IR, K
      PARAMETER (IA=16807, IM=2147483647, IQ=127773, IR=2836)
      IF (ISEED.LE.0) ISEED = 2937
      K = ISEED/IQ
      ISEED = IA*(ISEED - K*IQ) - IR*K
      IF (ISEED.LT.0) ISEED = ISEED + IM
      RANF = REAL(ISEED) * (1.0/REAL(IM))
      RETURN
      END
