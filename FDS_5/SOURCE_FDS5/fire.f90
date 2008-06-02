MODULE FIRE
 
! Compute combustion 
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: SECOND
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: fireid='$Id$'
CHARACTER(255), PARAMETER :: firerev='$Revision$'
CHARACTER(255), PARAMETER :: firedate='$Date$'

TYPE(REACTION_TYPE), POINTER :: RN

PUBLIC COMBUSTION, GET_REV_fire
 
CONTAINS
 

SUBROUTINE COMBUSTION(NM)

INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()

CALL POINT_TO_MESH(NM)

! Choose between mixture fraction formulation or finite-rate chemistry

IF (MIXTURE_FRACTION) THEN
   CALL COMBUSTION_MF
ELSE
   CALL COMBUSTION_FR
ENDIF

! Mirror Q across external boundaries -- for Smokeview visualization only

CALL COMBUSTION_BC

TUSED(10,NM)=TUSED(10,NM)+SECOND()-TNOW

CONTAINS


SUBROUTINE COMBUSTION_MF

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_MOLECULAR_WEIGHT
REAL(EB) :: Y_FU_0,A,ETRM,Y_O2_0,Y_CO_0,DYF,DX_FDT,HFAC_F,DTT,DELTA2,& 
            Y_O2_MAX,TMP_MIN,Y_O2_CORR,Q_NEW,Q_OLD,F_TO_CO,DELTAH_CO,DYCO,HFAC_CO,RHOX, &
            X_FU,X_O2,X_FU_0,X_O2_0,X_FU_S,X_O2_S,X_FU_N,X_O2_N,CO_TO_O2, &
            Y_FU_MAX,TMP_F_MIN,Y_F_CORR,Z_2_MIN,WGT,OMWGT
REAL(EB), PARAMETER :: Y_FU_MIN=1.E-10_EB,Y_O2_MIN=1.E-10_EB,X_O2_MIN=1.E-16_EB,X_FU_MIN=1.E-16_EB,Y_CO_MIN=1.E-10_EB
INTEGER :: NODETS,N,I,J,K,II,JJ,KK,IOR,IC,IW,IWA(-3:3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: Y_O2,R_SUM_DILUENTS,MIX_TIME

! Weighting factor for time-averaging of local HRR

WGT   = MIN(1._EB,DT/Q_AVG_TIME)
OMWGT = 1._EB - WGT

! Misc initializations

Y_O2     => WORK1
Y_O2     =  0._EB
MIX_TIME => WORK2
MIX_TIME =  DT
Q        =  0._EB

! Compute and save O2 in all cells

IF (.NOT.CO_PRODUCTION) THEN
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            CALL GET_MASS_FRACTION(YY(I,J,K,I_FUEL),0._EB,YY(I,J,K,I_PROG_F),O2_INDEX,Y_SUM(I,J,K),Y_O2(I,J,K))   
         ENDDO
      ENDDO
   ENDDO
ELSE
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            CALL GET_MASS_FRACTION(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),O2_INDEX,Y_SUM(I,J,K),Y_O2(I,J,K))   
         ENDDO
      ENDDO
   ENDDO
ENDIF

! Compute the (fast) reaction of fuel to either CO or CO2

RN => REACTION(1)
HFAC_F  = RN%HEAT_OF_COMBUSTION/DT

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE
         Y_FU_0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))*RN%Y_F_INLET
         IF (Y_FU_0<=Y_FU_MIN) CYCLE
         Y_O2_0  = Y_O2(I,J,K)
         IF (Y_O2_0<=Y_O2_MIN) CYCLE

         IF_SUPPRESSION: IF (SUPPRESSION) THEN  ! Get maximum O2 in the current and neighboring cells to see if flame viable

            Y_O2_MAX  = 0._EB
            Y_FU_MAX  = 0._EB
            TMP_MIN   = TMPA
            TMP_F_MIN = TMPA

            SEARCH_LOOP: DO IOR=-3,3

               II = I
               JJ = J
               KK = K
               SELECT CASE(IOR)
                  CASE(-1)
                     II = I-1
                  CASE( 1)
                     II = I+1
                  CASE(-2)
                     JJ = J-1
                  CASE( 2)
                     JJ = J+1
                  CASE(-3)
                     KK = K-1
                  CASE( 3)
                     KK = K+1
               END SELECT

               IW = WALL_INDEX(IC,IOR)
               IF (IW==0 .OR. BOUNDARY_TYPE(IW)==OPEN_BOUNDARY .OR. BOUNDARY_TYPE(IW)==INTERPOLATED_BOUNDARY) THEN
                  IF (Y_O2(II,JJ,KK)>Y_O2_MAX) THEN
                     Y_O2_MAX = Y_O2(II,JJ,KK)
                     TMP_MIN  = TMP(II,JJ,KK)
                  ENDIF
                  IF (YY(II,JJ,KK,I_FUEL)>Y_FU_MAX) THEN
                     Y_FU_MAX   = YY(II,JJ,KK,I_FUEL)
                     TMP_F_MIN = TMP(II,JJ,KK)
                  ENDIF
               ENDIF

            ENDDO SEARCH_LOOP

            ! Evaluate empirical extinction criteria

            Y_O2_CORR = RN%Y_O2_LL*(RN%CRIT_FLAME_TMP-TMP_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
            Y_F_CORR  = RN%Y_F_LFL*(RN%CRIT_FLAME_TMP-TMP_F_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
            IF (Y_O2_MAX < Y_O2_CORR .OR. Y_FU_MAX*RN%Y_F_INLET < Y_F_CORR) CYCLE

         ENDIF IF_SUPPRESSION

         DYF = MIN(Y_FU_0,Y_O2_0/RN%O2_F_RATIO)
         IF (LES .AND. EDDY_BREAKUP) THEN
            IF (.NOT.TWO_D) THEN
               DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
            ELSE
               DELTA2 = DX(I)*DZ(K)
            ENDIF
            MIX_TIME(I,J,K) = MIX_TIME_FACTOR*SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
         ENDIF
         Q_NEW = MIN((Q_AVG_MAX-OMWGT*Q_AVG(I,J,K))/WGT,Q_UPPER,DYF*RHO(I,J,K)*HFAC_F*MIN(1._EB,DT/MIX_TIME(I,J,K)))
         Q_AVG(I,J,K) = OMWGT*Q_AVG(I,J,K) + WGT*Q_NEW
         DYF = Q_NEW /(RHO(I,J,K)*HFAC_F*RN%Y_F_INLET)
         Q(I,J,K)   = Q_NEW
         YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
         IF (CO_PRODUCTION) THEN
            YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) + DYF
         ELSE
            YY(I,J,K,I_PROG_F)  = YY(I,J,K,I_PROG_F)  + DYF
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Optional second (slow) reaction to convert CO to CO_2
   
CONVERT_CO: IF (CO_PRODUCTION) THEN 

   RN => REACTION(2)
   DELTAH_CO = (REACTION(2)%HEAT_OF_COMBUSTION - REACTION(1)%HEAT_OF_COMBUSTION) * &
                REACTION(1)%MW_FUEL/((REACTION(1)%NU_CO-REACTION(2)%NU_CO)*MW_CO)
   F_TO_CO   = REACTION(1)%MW_FUEL/(REACTION(1)%NU_CO*MW_CO)  
   HFAC_CO   = DELTAH_CO/DT
   CO_TO_O2  = MW_O2/MW_O2*0.5_EB
   A  = RN%BOF 
   NODETS = 20
   DTT    = DT/REAL(NODETS,EB)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

            Y_O2_0  = Y_O2(I,J,K)
            Y_CO_0  = MAX(0._EB,YY(I,J,K,I_PROG_CO))*RN%Y_F_INLET / F_TO_CO 
            IF (Y_CO_0<=Y_CO_MIN .OR. Y_O2_0<=Y_O2_MIN) CYCLE

            ! Get max conversion allowed

            Z_2_MIN = REACTION(2)%NU_CO / (REACTION(2)%NU_CO+REACTION(2)%NU_CO2) * (YY(I,J,K,I_PROG_CO)+YY(I,J,K,I_PROG_F))
            IF (YY(I,J,K,I_PROG_CO) < Z_2_MIN) CYCLE
            IF((TMP(I,J,K) < 500._EB .AND. Q(I,J,K)==0._EB) .OR. Q(I,J,K)>=Q_UPPER) CYCLE

            ! Compute slow reaction

            IF (Q(I,J,K)/=0._EB) THEN
               DYCO = MIN(Y_CO_0,Y_O2_0/CO_TO_O2)*MIN(1._EB,DT/MIX_TIME(I,J,K))
            ELSE
               RHOX = 1000._EB*RHO(I,J,K)
               X_FU_0 = Y_CO_0*RHOX/MW_CO*1.E-6_EB
               X_O2_0 = Y_O2_0*RHOX/MW_O2*1.E-6_EB
               X_FU   = X_FU_0
               X_O2   = X_O2_0
               ETRM = EXP(-RN%E/(R0*TMP(I,J,K)))   
               ODE_LOOP2: DO II=1,NODETS
                  IF (X_FU<=X_FU_MIN .OR. X_O2<=X_O2_MIN) EXIT ODE_LOOP2
                  DX_FDT= -A*ETRM*X_FU*X_O2
                  X_FU_S = X_FU + DTT*DX_FDT
                  X_O2_S = X_O2 + DTT*DX_FDT*0.5_EB
                  IF (X_O2_S<=X_O2_MIN) THEN
                     X_FU = MAX(0.0_EB,X_FU-X_O2*2._EB)
                     EXIT ODE_LOOP2
                  ENDIF
                  IF (X_FU_S<=X_FU_MIN) THEN
                     X_FU = X_FU_MIN
                     EXIT ODE_LOOP2
                  ENDIF
                  DX_FDT = -A*ETRM*X_FU_S*X_O2_S
                  X_FU_N = .5_EB*(X_FU+X_FU_S+DTT*DX_FDT)
                  X_O2_N = .5_EB*(X_O2+X_O2_S+DTT*DX_FDT*0.5_EB)
                  IF (X_O2_N<=X_O2_MIN) THEN
                     X_FU = MAX(0.0_EB,0.5_EB*((X_FU+X_FU_S)-(X_O2+X_O2_S)*2._EB))
                     EXIT ODE_LOOP2
                  ENDIF
                  IF (X_FU_N<=X_FU_MIN) THEN
                     X_FU = X_FU_MIN
                     EXIT ODE_LOOP2
                  ENDIF
                  X_FU  = X_FU_N
                  X_O2  = X_O2_N
               ENDDO ODE_LOOP2
               DYCO = MIN(X_FU_0,X_FU_0 - X_FU)*MW_CO/RHOX*1.E6
            ENDIF
            Q_OLD = Q(I,J,K)
            Q_NEW = MIN(Q_UPPER-Q_OLD,DYCO*RHO(I,J,K)*HFAC_CO)
            DYCO = Q_NEW/(RHO(I,J,K)*HFAC_CO*RN%Y_F_INLET) * F_TO_CO
            IF (YY(I,J,K,I_PROG_CO) - DYCO < Z_2_MIN) THEN
               Q_NEW = Q_NEW*(YY(I,J,K,I_PROG_CO)-Z_2_MIN)/DYCO
               DYCO  = YY(I,J,K,I_PROG_CO) - Z_2_MIN
            ENDIF
            Q(I,J,K)   = Q_OLD + Q_NEW
            YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) - DYCO
            YY(I,J,K,I_PROG_F)  = YY(I,J,K,I_PROG_F)  + DYCO
         ENDDO
      ENDDO
   ENDDO

ENDIF CONVERT_CO

! Sum up diluent mass fractions

IF (N_SPEC_DILUENTS > 0) THEN
   R_SUM_DILUENTS => WORK4
   R_SUM_DILUENTS =  0._EB
   DO N=1,N_SPECIES
      IF (SPECIES(N)%MODE==GAS_SPECIES) R_SUM_DILUENTS(:,:,:) = R_SUM_DILUENTS(:,:,:) + SPECIES(N)%RCON*YY(:,:,:,N)
   ENDDO
ENDIF

! Compute new mixture molecular weight

IF (CO_PRODUCTION) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            CALL GET_MOLECULAR_WEIGHT(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
            RSUM(I,J,K) = R0/RSUM(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ELSE
   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            CALL GET_MOLECULAR_WEIGHT(YY(I,J,K,I_FUEL),0._EB,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
            RSUM(I,J,K) = R0/RSUM(I,J,K)
         ENDDO
      ENDDO
   ENDDO
ENDIF


IF (N_SPEC_DILUENTS > 0) RSUM = RSUM*(1._EB-Y_SUM) + R_SUM_DILUENTS

END SUBROUTINE COMBUSTION_MF



SUBROUTINE COMBUSTION_FR

! Finite-Rate Combustion

REAL(EB) :: TMPD,ETRM,&
            DX_FDT,DTT,YYNEW,WFAC,QFAC,&
            X_O2_MIN,X_FU_MIN,X_FU_S,X_O2_S
REAL(EB), DIMENSION(1:N_REACTIONS,0:N_SPECIES) :: MPUE
REAL(EB), DIMENSION(1:N_REACTIONS) :: HFAC_F, A, ETRM1, X_FU,X_O2, &
           DYF,Q_NR,X_FU_N,X_O2_N
INTEGER :: NODETS,N,I,J,K,II,NR
LOGICAL :: NO_REACTION

Q =  0._EB

DO NR=1,N_REACTIONS
   RN => REACTION(NR)
   HFAC_F(NR) = RN%HEAT_OF_COMBUSTION/DT
   DO N=1,N_SPECIES
      MPUE(NR,N) = SPECIES(N)%MW*RN%NU(N)/ABS(RN%EPUMO2*SPECIES(RN%I_OXIDIZER)%MW*RN%NU(RN%I_OXIDIZER))
   ENDDO
   A(NR) = RN%BOF
ENDDO 
NODETS = 20
DTT    = DT/REAL(NODETS,EB)
X_O2_MIN = 1.E-14_EB
X_FU_MIN = 1.E-14_EB

DO K=1,KBAR
   DO J=1,JBAR
      ILOOP3: DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOP3
         DYF = 0._EB
         !Convert mass fractions to concentrations
         TMPD = TMP(I,J,K)
         DO NR = 1, N_REACTIONS
            ETRM1(NR) = A(NR)*EXP(-REACTION(NR)%E/(R0*TMPD))
         ENDDO
         NO_REACTION = .FALSE.
         ODE_LOOP: DO II=1,NODETS
            X_O2_N  = 0._EB
            X_FU_N = 0._EB
            REACTION_LOOP: DO NR=1,N_REACTIONS
               Q_NR(NR) = 0._EB
               RN => REACTION(NR)
               X_O2(NR) = YY(I,J,K,RN%I_OXIDIZER)*1000._EB*RHO(I,J,K)/SPECIES(RN%I_OXIDIZER)%MW*1.E-6_EB
               X_FU(NR) = YY(I,J,K,RN%I_FUEL)*1000._EB    *RHO(I,J,K)/RN%MW_FUEL*1.E-6_EB
               IF (X_FU(NR)<=X_FU_MIN .OR. X_O2(NR)<=X_O2_MIN) THEN
                  IF (NR == 1) THEN
                     NO_REACTION = .TRUE.
                  ELSE
                     NO_REACTION = NO_REACTION .AND. .TRUE.
                  ENDIF
                  CYCLE REACTION_LOOP
               ELSE
                  NO_REACTION = .FALSE.
               ENDIF
               ETRM = ETRM1(NR)

               ! Local flame extinction criteria

               SPEC_LOOP: DO N=1,N_SPECIES
                  IF (N==RN%I_FUEL .OR. N==RN%I_OXIDIZER .OR. RN%N(N)==-999._EB) CYCLE SPEC_LOOP
                  IF (YY(I,J,K,N) < 1.E-20_EB) CYCLE SPEC_LOOP
                  ETRM = ETRM * (YY(I,J,K,N)*1000._EB*RHO(I,J,K)/SPECIES(N)%MW*1.E-6_EB)**RN%N(N)
               ENDDO SPEC_LOOP

               ! Solve the simple ODE to deplete fuel and oxidizer due to reaction

               DX_FDT= -ETRM*X_FU(NR)**RN%N_F*X_O2(NR)**RN%N_O
               X_FU_S = X_FU(NR) + DTT*DX_FDT
               X_O2_S = X_O2(NR) + DTT*DX_FDT*RN%NU(RN%I_OXIDIZER)/RN%NU(RN%I_FUEL)
               IF (X_O2_S<X_O2_MIN) THEN
                  X_O2_S = X_O2_MIN 
                  X_FU_S = MAX(0.0_EB,X_FU(NR)-(X_O2(NR)-X_O2_S)*RN%NU(RN%I_FUEL)/RN%NU(RN%I_OXIDIZER))
               ENDIF  
               IF (X_FU_S<X_FU_MIN) THEN
                  X_FU_S = X_FU_MIN
                  X_O2_S = MAX(0.0_EB,X_O2(NR)-(X_FU(NR)-X_FU_S)*RN%NU(RN%I_OXIDIZER)/RN%NU(RN%I_FUEL))
               ENDIF
               DX_FDT= -ETRM*X_FU_S**RN%N_F*X_O2_S**RN%N_O
               X_FU_N(NR) = .5_EB*(X_FU(NR)+X_FU_S+DTT*DX_FDT)
               X_O2_N(NR) = .5_EB*(X_O2(NR)+X_O2_S+DTT*DX_FDT*RN%NU(RN%I_OXIDIZER)/RN%NU(RN%I_FUEL))
               IF (X_O2_N(NR)<X_O2_MIN) &
                  X_FU_N(NR) = MAX(0.0_EB,0.5_EB*((X_FU(NR)+X_FU_S)-(X_O2(NR)+X_O2_S)*RN%NU(RN%I_FUEL)/RN%NU(RN%I_OXIDIZER)))
               IF (X_FU_N(NR)<X_FU_MIN) X_FU_N(NR) = X_FU_MIN
               DYF(NR) = MAX(-X_FU(NR),X_FU_N(NR)-X_FU(NR))*RN%MW_FUEL*0.001_EB*1.E6_EB
               Q_NR(NR) = -DYF(NR) * HFAC_F(NR)
               IF (Q(I,J,K) + Q_NR(NR) > Q_UPPER) THEN
                  QFAC = (1._EB - (Q(I,J,K) + Q_NR(NR) - Q_UPPER) / Q_NR(NR))
                  DYF(NR) = DYF(NR) * QFAC
                  Q_NR(NR) = Q_NR(NR) * QFAC
                  Q(I,J,K) = Q_UPPER
               ENDIF   
               DO N=1,N_SPECIES
                  YYNEW = YY(I,J,K,N) + DT*MPUE(NR,N)*Q_NR(NR)/RHO(I,J,K)
                  YY(I,J,K,N) = MAX(YYMIN(N),YYNEW)
               ENDDO
               IF (Q(I,J,K)==Q_UPPER) EXIT ODE_LOOP
               Q(I,J,K) = Q(I,J,K) + Q_NR(NR)
            ENDDO REACTION_LOOP
            IF (NO_REACTION) EXIT ODE_LOOP
            NO_REACTION = .FALSE.

         ENDDO ODE_LOOP
      ENDDO ILOOP3
   ENDDO
ENDDO

! Adjust the average molecular weight term, R*Sum(Yi/Mi)

RSUM = SPECIES(0)%RCON
SLOOP: DO N=1,N_SPECIES
   IF (SPECIES(N)%MODE/=GAS_SPECIES) CYCLE SLOOP         
   WFAC = SPECIES(N)%RCON - SPECIES(0)%RCON
   RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YY(:,:,:,N)
ENDDO SLOOP

END SUBROUTINE COMBUSTION_FR


SUBROUTINE COMBUSTION_BC

INTEGER :: II,JJ,KK,IIG,JJG,KKG,IW

DO IW=1,NEWC
   IF (BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY .AND. BOUNDARY_TYPE(IW)/=OPEN_BOUNDARY) CYCLE
   II  = IJKW(1,IW)
   JJ  = IJKW(2,IW)
   KK  = IJKW(3,IW)
   IIG = IJKW(6,IW)
   JJG = IJKW(7,IW)
   KKG = IJKW(8,IW)
   Q(II,JJ,KK) = Q(IIG,JJG,KKG)
ENDDO

END SUBROUTINE COMBUSTION_BC

END SUBROUTINE COMBUSTION



SUBROUTINE GET_REV_fire(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') firerev(INDEX(firerev,':')+1:LEN_TRIM(firerev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') firedate

END SUBROUTINE GET_REV_fire
 
END MODULE FIRE

