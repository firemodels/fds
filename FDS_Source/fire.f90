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
REAL(EB) :: Q_UPPER

PUBLIC COMBUSTION, GET_REV_fire
 
CONTAINS
 

SUBROUTINE COMBUSTION(NM)

INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()

CALL POINT_TO_MESH(NM)

! Upper bounds on local HRR per unit volume

Q_UPPER = HRRPUA_SHEET/CELL_SIZE

! Choose between mixture fraction formulation or finite-rate chemistry

IF (MIXTURE_FRACTION) THEN
   IF (COMBUSTION2) THEN
      CALL COMBUSTION_MF2
   ELSE
      CALL COMBUSTION_MF   
   ENDIF
ELSE
   IF (COMBUSTION2) THEN
      CALL COMBUSTION_FR2
   ELSE
      CALL COMBUSTION_FR
   ENDIF
ENDIF

! Mirror Q across external boundaries -- for Smokeview visualization only

CALL COMBUSTION_BC

TUSED(10,NM)=TUSED(10,NM)+SECOND()-TNOW

CONTAINS

SUBROUTINE COMBUSTION_MF

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_SPECIFIC_GAS_CONSTANT
REAL(EB) :: Y_FU_0,A,ETRM,Y_O2_0,Y_CO_0,DYF,DX_FDT,HFAC_F,DTT,DELTA2,& 
            Y_O2_MAX,TMP_MIN,Y_O2_CORR,Q_NEW,Q_OLD,F_TO_CO,DELTAH_CO,DYCO,HFAC_CO,RHOX, &
            X_FU,X_O2,X_FU_0,X_O2_0,X_FU_S,X_O2_S,X_FU_N,X_O2_N,CO_TO_O2, &
            Y_FU_MAX,TMP_F_MIN,Y_F_CORR,Z_2_MIN,Z_2_MIN_FAC,WGT,OMWGT,Q_BOUND_1,Q_BOUND_2,Q_BOUND_3,YY_GET(1:N_SPECIES)
REAL(EB), PARAMETER :: Y_FU_MIN=1.E-10_EB,Y_O2_MIN=1.E-10_EB,X_O2_MIN=1.E-16_EB,X_FU_MIN=1.E-16_EB,Y_CO_MIN=1.E-10_EB
INTEGER :: NODETS,I,J,K,II,JJ,KK,IOR,IC,IW,IWA(-3:3)
REAL(EB), POINTER, DIMENSION(:,:,:) :: Y_O2,Y_O2_NEW,MIX_TIME

! Weighting factor for time-averaging of local HRR

WGT   = MIN(1._EB,DT/HRRPUV_AVERAGE_TIME)
OMWGT = 1._EB - WGT

! Misc initializations

Y_O2     => WORK1
Y_O2     =  0._EB
MIX_TIME => WORK2
MIX_TIME =  DT
Q        =  0._EB

! Compute and save O2 in all cells


   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            YY_GET(:) = YY(I,J,K,:)
            CALL GET_MASS_FRACTION(YY_GET,O2_INDEX,Y_O2(I,J,K))   
         ENDDO
      ENDDO
   ENDDO
IF (CO_PRODUCTION) THEN
   Y_O2_NEW => WORK3
   Y_O2_NEW =  Y_O2
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

            ! Get O2 and temperature in nearest neighbor cells

            IF (SUPPRESSION_SEARCH) THEN

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

            ENDIF

            ! Evaluate empirical extinction criteria
            
            Y_O2_CORR = RN%Y_O2_LL*(RN%CRIT_FLAME_TMP-TMP_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
            Y_F_CORR  = RN%Y_F_LFL*(RN%CRIT_FLAME_TMP-TMP_F_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
            IF (Y_O2_MAX < Y_O2_CORR .OR. Y_FU_MAX*RN%Y_F_INLET < Y_F_CORR) CYCLE

         ENDIF IF_SUPPRESSION
         DYF = MIN(Y_FU_0,Y_O2_0/RN%O2_F_RATIO)
         IF (LES .AND. EDDY_DISSIPATION) THEN
            IF (.NOT.TWO_D) THEN
               DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
            ELSE
               DELTA2 = DX(I)*DZ(K)
            ENDIF
            MIX_TIME(I,J,K) = C_EDC*SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
         ENDIF
         Q_BOUND_1 = DYF*RHO(I,J,K)*HFAC_F*MIN(1._EB,DT/MIX_TIME(I,J,K))
         Q_BOUND_2 = Q_UPPER
         IF (LES) THEN
            Q_BOUND_3 = (HRRPUV_AVERAGE-OMWGT*Q_AVG(I,J,K))/WGT
            Q_NEW = MIN(Q_BOUND_1,Q_BOUND_2,Q_BOUND_3)
         ELSE
            Q_NEW = MIN(Q_BOUND_1,Q_BOUND_2)
         ENDIF
         DYF = Q_NEW /(RHO(I,J,K)*HFAC_F*RN%Y_F_INLET)         
         Q(I,J,K)   = Q_NEW
         YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
         IF (CO_PRODUCTION) THEN
            IF (SOOT_DEPOSITION) THEN
               YY(I,J,K,I_PROG_CO)   = YY(I,J,K,I_PROG_CO)   + DYF * (1._EB - RN%SOOT_YIELD)
               YY(I,J,K,I_PROG_SOOT) = YY(I,J,K,I_PROG_SOOT) + DYF * RN%SOOT_YIELD
            ELSE
               YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) + DYF            
            ENDIF
            Y_O2_NEW(I,J,K) = Y_O2_NEW(I,J,K) - DYF * RN%O2_F_RATIO
         ELSE
            IF (SOOT_DEPOSITION) THEN
               YY(I,J,K,I_PROG_F)    = YY(I,J,K,I_PROG_F)    + DYF * (1._EB - RN%SOOT_YIELD)
               YY(I,J,K,I_PROG_SOOT) = YY(I,J,K,I_PROG_SOOT) + DYF * RN%SOOT_YIELD
            ELSE
               YY(I,J,K,I_PROG_F)  = YY(I,J,K,I_PROG_F)  + DYF            
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Optional second (slow) reaction to convert CO to CO_2
   CONVERT_CO: IF (CO_PRODUCTION) THEN 
   RN => REACTION(2)
   DELTAH_CO = (REACTION(2)%HEAT_OF_COMBUSTION - REACTION(1)%HEAT_OF_COMBUSTION) * &
                REACTION(1)%MW_FUEL/((REACTION(1)%NU(CO_INDEX)-REACTION(2)%NU(CO_INDEX))*MW_CO)
   F_TO_CO   = REACTION(1)%MW_FUEL/(REACTION(1)%NU(CO_INDEX)*MW_CO)  
   HFAC_CO   = DELTAH_CO/DT
   CO_TO_O2  = MW_O2/(MW_CO*2._EB)
   Z_2_MIN_FAC = F_TO_CO * REACTION(2)%CO_YIELD / REACTION(1)%Y_F_INLET
   A  = RN%BOF 
   NODETS = 20
   DTT    = DT/REAL(NODETS,EB)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

            Y_O2_0  = Y_O2_NEW(I,J,K)
            Y_CO_0  = MAX(0._EB,YY(I,J,K,I_PROG_CO))*RN%Y_F_INLET / F_TO_CO 
            IF (Y_CO_0<=Y_CO_MIN .OR. Y_O2_0<=Y_O2_MIN) CYCLE

            ! Get max conversion allowed

            Z_2_MIN = Z_2_MIN_FAC * (YY(I,J,K,I_PROG_CO)+YY(I,J,K,I_PROG_F))
            IF (YY(I,J,K,I_PROG_CO) < Z_2_MIN) CYCLE
            IF (LES) THEN
               Q_BOUND_2 = MIN((HRRPUV_AVERAGE-OMWGT*Q_AVG(I,J,K))/WGT,Q_UPPER)
            ELSE
               Q_BOUND_2 = Q_UPPER
            ENDIF
            IF((TMP(I,J,K) < 500._EB .AND. Q(I,J,K)==0._EB) .OR. Q(I,J,K)>=Q_BOUND_2) CYCLE

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
            Q_NEW = MIN(Q_BOUND_2-Q_OLD,DYCO*RHO(I,J,K)*HFAC_CO)
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

IF (LES) Q_AVG = OMWGT*Q_AVG + WGT*Q

! Compute new mixture molecular weight

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         YY_GET(:) = YY(I,J,K,:)
         CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE COMBUSTION_MF


SUBROUTINE COMBUSTION_MF2

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_SPECIFIC_GAS_CONSTANT,GET_AVERAGE_SPECIFIC_HEAT
REAL(EB) :: Y_FU_0,A,ETRM,Y_O2_0,Y_CO_0,DYF,DX_FDT,HFAC_F,DTT,DELTA2,& 
            Y_O2_MAX,TMP_MIN,Y_O2_CORR,Q_NEW,Q_OLD,F_TO_CO,DELTAH_CO,DELTAH_F,DYCO,HFAC_CO,RHOX,H_F_0,H_G_0,H_G_N, &
            X_FU,X_O2,X_FU_0,X_O2_0,X_FU_S,X_O2_S,X_FU_N,X_O2_N,CO_TO_O2, &
            Y_FU_MAX,TMP_F_MIN,Y_F_CORR,Z_2_MIN,Z_2_MIN_FAC,WGT,OMWGT,Q_BOUND_1,Q_BOUND_2,Q_BOUND_3,YY_GET(1:N_SPECIES)
REAL(EB), PARAMETER :: Y_FU_MIN=1.E-10_EB,Y_O2_MIN=1.E-10_EB,X_O2_MIN=1.E-16_EB,X_FU_MIN=1.E-16_EB,Y_CO_MIN=1.E-10_EB
INTEGER :: NODETS,I,J,K,II,JJ,KK,IOR,IC,IW,IWA(-3:3),ITMP,ICFT
REAL(EB), POINTER, DIMENSION(:,:,:) :: Y_O2,Y_O2_NEW,MIX_TIME

! Weighting factor for time-averaging of local HRR

WGT   = MIN(1._EB,DT/HRRPUV_AVERAGE_TIME)
OMWGT = 1._EB - WGT

! Misc initializations

Y_O2     => WORK1
Y_O2     =  0._EB
MIX_TIME => WORK2
MIX_TIME =  DT
Q        =  0._EB

! Compute and save O2 in all cells


   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            YY_GET(:) = YY(I,J,K,:)
            CALL GET_MASS_FRACTION(YY_GET,O2_INDEX,Y_O2(I,J,K))   
         ENDDO
      ENDDO
   ENDDO
IF (CO_PRODUCTION) THEN
   Y_O2_NEW => WORK3
   Y_O2_NEW =  Y_O2
ENDIF

! Compute the (fast) reaction of fuel to either CO or CO2

RN => REACTION(1)
HFAC_F  = RN%HEAT_OF_COMBUSTION/DT
IF (CO_PRODUCTION) THEN
   DELTAH_F = REACTION(2)%HEAT_OF_COMBUSTION
ELSE
   DELTAH_F = REACTION(1)%HEAT_OF_COMBUSTION
ENDIF
ICFT = NINT(MIN(5000._EB,RN%CRIT_FLAME_TMP))
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE
         Y_FU_0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))*RN%Y_F_INLET
         IF (Y_FU_0<=Y_FU_MIN) CYCLE
         Y_O2_0  = Y_O2(I,J,K)
         IF (Y_O2_0<=Y_O2_MIN) CYCLE
         DYF = MIN(Y_FU_0,Y_O2_0/RN%O2_F_RATIO)         
         IF_SUPPRESSION: IF (SUPPRESSION) THEN  ! Determine if combustion can yield critical flame temperature
            ITMP = NINT(MIN(5000._EB,TMP(I,J,K)))
            YY_GET = 0._EB
            YY_GET(I_FUEL) = 1._EB
            CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_F_0,ITMP)
            YY_GET = YY(I,J,K,:)
            YY_GET(I_FUEL) = 0._EB
            YY_GET = YY_GET / (1 - Y_FU_0)
            CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_0,ITMP)
            YY_GET = YY(I,J,K,:)
            YY_GET(I_FUEL) = YY_GET(I_FUEL) - DYF
            YY_GET(I_PROG_F) = YY_GET(I_PROG_F) + DYF
            CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_N,ICFT)
            IF ( (DYF*H_F_0 + (1._EB-DYF)*MIN(1._EB,DYF * RN%O2_F_RATIO/Y_O2_0)*H_G_0)*TMP(I,J,K) + DYF*DELTAH_F < &
             (DYF+(1._EB-DYF)*MIN(1._EB,DYF * RN%O2_F_RATIO/Y_O2_0))*H_G_N*RN%CRIT_FLAME_TMP) CYCLE
         ENDIF IF_SUPPRESSION
         IF (LES .AND. EDDY_DISSIPATION) THEN
            IF (.NOT.TWO_D) THEN
               DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
            ELSE
               DELTA2 = DX(I)*DZ(K)
            ENDIF
            MIX_TIME(I,J,K) = C_EDC*SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
         ENDIF
         Q_BOUND_1 = DYF*RHO(I,J,K)*HFAC_F*MIN(1._EB,DT/MIX_TIME(I,J,K))
         Q_BOUND_2 = Q_UPPER
         IF (LES) THEN
            Q_BOUND_3 = (HRRPUV_AVERAGE-OMWGT*Q_AVG(I,J,K))/WGT
            Q_NEW = MIN(Q_BOUND_1,Q_BOUND_2,Q_BOUND_3)
         ELSE
            Q_NEW = MIN(Q_BOUND_1,Q_BOUND_2)
         ENDIF
         DYF = Q_NEW /(RHO(I,J,K)*HFAC_F*RN%Y_F_INLET)         
         Q(I,J,K)   = Q_NEW
         YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
         IF (CO_PRODUCTION) THEN
            IF (SOOT_DEPOSITION) THEN
               YY(I,J,K,I_PROG_CO)   = YY(I,J,K,I_PROG_CO)   + DYF * (1._EB - RN%SOOT_YIELD)
               YY(I,J,K,I_PROG_SOOT) = YY(I,J,K,I_PROG_SOOT) + DYF * RN%SOOT_YIELD
            ELSE
               YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) + DYF            
            ENDIF
            Y_O2_NEW(I,J,K) = Y_O2_NEW(I,J,K) - DYF * RN%O2_F_RATIO
         ELSE
            IF (SOOT_DEPOSITION) THEN
               YY(I,J,K,I_PROG_F)    = YY(I,J,K,I_PROG_F)    + DYF * (1._EB - RN%SOOT_YIELD)
               YY(I,J,K,I_PROG_SOOT) = YY(I,J,K,I_PROG_SOOT) + DYF * RN%SOOT_YIELD
            ELSE
               YY(I,J,K,I_PROG_F)  = YY(I,J,K,I_PROG_F)  + DYF            
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Optional second (slow) reaction to convert CO to CO_2
   CONVERT_CO: IF (CO_PRODUCTION) THEN 
   RN => REACTION(2)
   DELTAH_CO = (REACTION(2)%HEAT_OF_COMBUSTION - REACTION(1)%HEAT_OF_COMBUSTION) * &
                REACTION(1)%MW_FUEL/((REACTION(1)%NU(CO_INDEX)-REACTION(2)%NU(CO_INDEX))*MW_CO)
   F_TO_CO   = REACTION(1)%MW_FUEL/(REACTION(1)%NU(CO_INDEX)*MW_CO)  
   HFAC_CO   = DELTAH_CO/DT
   CO_TO_O2  = MW_O2/(MW_CO*2._EB)
   Z_2_MIN_FAC = F_TO_CO * REACTION(2)%CO_YIELD / REACTION(1)%Y_F_INLET
   A  = RN%BOF 
   NODETS = 20
   DTT    = DT/REAL(NODETS,EB)

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE

            Y_O2_0  = Y_O2_NEW(I,J,K)
            Y_CO_0  = MAX(0._EB,YY(I,J,K,I_PROG_CO))*RN%Y_F_INLET / F_TO_CO 
            IF (Y_CO_0<=Y_CO_MIN .OR. Y_O2_0<=Y_O2_MIN) CYCLE

            ! Get max conversion allowed

            Z_2_MIN = Z_2_MIN_FAC * (YY(I,J,K,I_PROG_CO)+YY(I,J,K,I_PROG_F))
            IF (YY(I,J,K,I_PROG_CO) < Z_2_MIN) CYCLE
            IF (LES) THEN
               Q_BOUND_2 = MIN((HRRPUV_AVERAGE-OMWGT*Q_AVG(I,J,K))/WGT,Q_UPPER)
            ELSE
               Q_BOUND_2 = Q_UPPER
            ENDIF
            IF((TMP(I,J,K) < 500._EB .AND. Q(I,J,K)==0._EB) .OR. Q(I,J,K)>=Q_BOUND_2) CYCLE

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
            Q_NEW = MIN(Q_BOUND_2-Q_OLD,DYCO*RHO(I,J,K)*HFAC_CO)
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

IF (LES) Q_AVG = OMWGT*Q_AVG + WGT*Q

! Compute new mixture molecular weight

DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         YY_GET(:) = YY(I,J,K,:)
         CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE COMBUSTION_MF2


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

SUBROUTINE COMBUSTION_FR2
USE PHYSICAL_FUNCTIONS, ONLY:GET_SPECIFIC_GAS_CONSTANT
! Finite-Rate Combustion

REAL(EB) :: TMPD,ETRM1(N_REACTIONS),Y_MIN,Y_MIN_MIN,X_REAC(N_REAC_SPECIES),YY_GET(N_SPECIES),RHOP,Q_TMP
INTEGER :: NODETS,N,I,J,K,II,NR,NS
LOGICAL :: NO_REACTION

Q =  0._EB
Y_MIN     = 1.E-14_EB

DO K=1,KBAR
   DO J=1,JBAR
      ILOOP3: DO I=1,IBAR
         IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOP3
         TMPD = TMP(I,J,K)
         CALL RATE_CONSTANT(TMPD,ETRM1)
         NO_REACTION = .TRUE.
! See if all rate constants or all reactant species are near zero         
         REACLOOP: DO NR = 1, N_REACTIONS
            IF (ETRM1(NR) < 1.E-30_EB) CYCLE REACLOOP
            Y_MIN_MIN = Y_MIN
            RN => REACTION(NR)
            SPECLOOP: DO NS = 1, RN%N_SPECIES
               IF (RN%NU(NS) < 0._EB) EXIT SPECLOOP
               Y_MIN_MIN = MIN(Y_MIN_MIN,YY(I,J,K,RN%SPECIES(NS)))
            ENDDO SPECLOOP
            IF (Y_MIN_MIN >=Y_MIN) NO_REACTION=.FALSE.
         ENDDO REACLOOP
         IF (NO_REACTION) CYCLE ILOOP3
         RHOP = RHO(I,J,K)
         YY_GET = YY(I,J,K,:)
         DO NS = 1,N_REAC_SPECIES
            X_REAC(NS) = YY_GET(REAC_SPECIES(NS)%SPEC_INDEX)/SPECIES(REAC_SPECIES(NS)%SPEC_INDEX)%MW*RHOP
         ENDDO
         CALL RK5AS(X_REAC,ETRM1,Q_TMP,DT)
         DO NS = 1,N_REAC_SPECIES
            YY_GET(REAC_SPECIES(NS)%SPEC_INDEX) = X_REAC(NS) * SPECIES(REAC_SPECIES(NS)%SPEC_INDEX)%MW / RHOP
         ENDDO
         YY(I,J,K,:) = YY_GET
         Q(I,J,K) = Q_TMP
      ENDDO ILOOP3
   ENDDO
ENDDO

! Adjust the average molecular weight term, R*Sum(Yi/Mi)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         YY_GET(:) = YY(I,J,K,:)
         CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
      ENDDO
   ENDDO
ENDDO

END SUBROUTINE COMBUSTION_FR2

SUBROUTINE RK5AS(X,KR,Q,DT)
!Adaptive time step, 5th order RK, chemical kinetics scheme
REAL(EB), INTENT(INOUT) :: X(N_SPECIES)
REAL(EB), INTENT(IN) :: KR(N_REACTIONS)
REAL(EB), PARAMETER :: B21=0.2_EB,&
            B31=0.075_EB,B32=0.225_EB,&
            B41=0.3_EB,B42=-0.9_EB,B43=1.2_EB,&
            B51=-11._EB/54._EB,B52=2.5_EB,B53=-70._EB/27._EB,B54=35._EB/27._EB,&
            B61=1631._EB/55296._EB,B62=175._EB/512._EB,B63=575._EB/13824._EB,B64=44275._EB/110592._EB,B65=253._EB/4096._EB
REAL(EB),PARAMETER :: C1=37._EB/378._EB,C2=0._EB,C3=250._EB/621._EB,C4=125._EB/594._EB,C5=0._EB,C6=512._EB/1771._EB
REAL(EB) :: DC1=C1-2825._EB/27648._EB,DC2=0._EB,DC3=C3-18575._EB/48384._EB,DC4=C4-13525._EB/55296._EB,DC5=C5-277._EB/14336._EB,&
            DC6=C6-0.25_EB
REAL(EB) :: EPS = 1.E-4_EB,POWERUP=-0.2_EB,POWERDOWN=-0.25_EB,RELAX=0.9_EB        
REAL(EB) :: XERR(N_REAC_SPECIES),XSCALE(N_REAC_SPECIES),XO(N_REAC_SPECIES),XN(N_REAC_SPECIES),K(6,N_REAC_SPECIES),DQ(6)
REAL(EB) :: DTSTEP,TSTEP,ERRMAX,MINSTEP,DT,Q,DQDT
INTEGER :: NS

DTSTEP = DT/20._EB
MINSTEP = 0.001_EB * DT
TSTEP = 0._EB

DO NS=1,N_REAC_SPECIES
   XO(NS) = X(REAC_SPECIES(NS)%SPEC_INDEX)
ENDDO
Q = 0._EB
DQ = 0._EB

TIMELOOP: DO 
  
   XN = XO
   CALL DERIV(XN,KR,DQDT)
   K(1,:) = DTSTEP*XN
   DQ(1) = DTSTEP*DQDT

   XN = XO+B21*K(1,:)
   CALL DERIV(XN,KR,DQDT)
   K(2,:) = DTSTEP*XN
   DQ(2) = DTSTEP*DQDT

   XN = XO+B31*K(1,:)+B32*K(2,:)
   CALL DERIV(XN,KR,DQDT)
   K(3,:) = DTSTEP*XN
   DQ(3) = DTSTEP*DQDT

   XN = XO+B41*K(1,:)+B42*K(2,:)+B43*K(3,:)
   CALL DERIV(XN,KR,DQDT)
   K(4,:) = DTSTEP*XN
   DQ(4) = DTSTEP*DQDT

   XN = XO+B51*K(1,:)+B52*K(2,:)+B53*K(3,:)+B54*K(4,:)
   CALL DERIV(XN,KR,DQDT)
   K(5,:) = DTSTEP*XN
   DQ(5) = DTSTEP*DQDT

   XN = XO+B61*K(1,:)+B62*K(2,:)+B63*K(3,:)+B64*K(4,:)
   CALL DERIV(XN,KR,DQDT)
   K(6,:) = DTSTEP*XN
   DQ(6) = DTSTEP*DQDT
   
   XN  = XO + C1 * K(1,:) +  C2 * K(2,:) +  C3 * K(3,:) +  C4 * K(4,:) +  C5  * K(5,:) +  C6  * K(6,:)
   XERR =   DC1 * K(1,:) + DC2 * K(2,:) + DC3 * K(3,:) + DC4 * K(4,:) + DC5  * K(5,:) + DC6  * K(6,:)

   XSCALE = ABS(X) + ZERO_P! + ABS(K(1,:)) + ZERO_P
   ERRMAX = MAXVAL(ABS(XERR) / XSCALE) / EPS

   IF (ERRMAX > 1._EB .AND. DTSTEP > MINSTEP) THEN      
      DTSTEP = MAX(MINSTEP,0.1_EB*DTSTEP,RELAX*DTSTEP*ERRMAX**POWERDOWN)
   ELSE
      XO = XN
      Q = Q + C1 * DQ(1) + C2 * DQ(2) + C3 * DQ(3) + C4 * DQ(4) + C5 * DQ(5) + C6 * DQ(6)
      TSTEP = TSTEP + DTSTEP
      IF (TSTEP >= DT) EXIT TIMELOOP
      IF (ERRMAX < 1._EB) DTSTEP = DTSTEP*MIN(5._EB,RELAX*ERRMAX**POWERUP)
      DTSTEP = MIN(DTSTEP,DT-TSTEP)
   ENDIF

ENDDO TIMELOOP

DO NS=1,N_REAC_SPECIES
   X(REAC_SPECIES(NS)%SPEC_INDEX) = XN(NS)
ENDDO

END SUBROUTINE RK5AS

SUBROUTINE DERIV(DXDT,KR,DQDT)
REAL(EB),INTENT(INOUT) :: DXDT(N_REAC_SPECIES),DQDT
REAL(EB),INTENT(IN) :: KR(N_REACTIONS)
REAL(EB) :: COLLISION_SUM,RATE(N_REACTIONS)
INTEGER :: NR,NS
TYPE (REAC_SPECIES_TYPE), POINTER :: RS
TYPE (REACTION_TYPE), POINTER :: RN

RATE = 0._EB
DQDT = 0._EB
DO NR = 1, N_REACTIONS
   RN => REACTION(NR)
   IF (RN%COLLISION) THEN
      COLLISION_SUM = 0._EB
      DO NS = 1, RN%N_COLLISION_SPECIES
         COLLISION_SUM = COLLISION_SUM + DXDT(RN%COLLISION_SPECIES(NS))
      ENDDO
      COLLISION_SUM = COLLISION_SUM**RN%COLLISION_EXPONENT
   ELSE
      COLLISION_SUM = 1._EB
   ENDIF
   RATE(NR) = KR(NR) * COLLISION_SUM
   DO NS = 1,RN%N_SPECIES
      RATE(NR) = RATE(NR) * X(RN%SPECIES(NS))**RN%N(NS)
   ENDDO
   DQDT = DQDT + RATE(NR)*RN%HEAT_OF_COMBUSTION
ENDDO

DO NS = 1,N_REAC_SPECIES
   DXDT(NS) = 0._EB
   RS => REAC_SPECIES(NS)
   DO NR = 1, RS%N_REACTIONS
      DXDT(NS) = DXDT(NS) + RATE(RS%REACTION_POINTER(NR,1)) * REACTION(RS%REACTION_POINTER(NR,1))%NU(RS%REACTION_POINTER(NR,2))
   ENDDO 
ENDDO

END SUBROUTINE DERIV

SUBROUTINE RATE_CONSTANT(TMP,KR)
REAL(EB), INTENT(IN) :: TMP
REAL(EB), INTENT(INOUT) :: KR(N_REACTIONS)
INTEGER :: NR
TYPE (REACTION_TYPE), POINTER :: RN

DO NR = 1, N_REACTIONS
   RN => REACTION(NR)
   KR(NR) = RN%BOF*EXP(-RN%E/(R0*TMP))
END DO

END SUBROUTINE RATE_CONSTANT

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

