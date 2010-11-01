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

TYPE(REACTION_TYPE), POINTER :: RN=>NULL()
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

Q_UPPER = HRRPUA_SHEET/CELL_SIZE**HRRPUA_SHEET_EXPONENT + HRRPUV_AVERAGE

! Choose between mixture fraction formulation or finite-rate chemistry

IF (MIXTURE_FRACTION) THEN
   CALL COMBUSTION_MF   
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

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_SPECIFIC_GAS_CONSTANT,GET_AVERAGE_SPECIFIC_HEAT,GET_CONDUCTIVITY
REAL(EB) :: Y_FU_0,A,ETRM,Y_O2_0,Y_CO_0,DYF,DX_FDT,HFAC_F,DTT,DELTA,DELTA2,ACCEL, & 
            Y_O2_MAX,TMP_MIN,Y_O2_CORR,Q_NEW,Q_OLD,F_TO_CO,DELTAH_CO,DYCO,HFAC_CO,RHOX, &
            X_FU,X_O2,X_FU_0,X_O2_0,X_FU_S,X_O2_S,X_FU_N,X_O2_N,CO_TO_O2,CRIT_FLAME_TMPA, &
            Y_FU_MAX,TMP_F_MIN,Y_F_CORR,Z_2_MIN,Z_2_MIN_FAC,WGT,OMWGT,Q_BOUND_1,Q_BOUND_2,Q_BOUND_3,YY_GET(1:N_SPECIES), &
            ZETA,CS2,H_F_0,H_F_N,H_G_0,H_G_N,DYAIR,DELTAH_F,TAU_D,TAU_U,TAU_G,EPSK,KSGS,KP,CP,S_L,TAU_CHEM, &
            DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,SS2,S12,S13,S23
REAL(EB), PARAMETER :: Y_FU_MIN=1.E-10_EB,Y_O2_MIN=1.E-10_EB,X_O2_MIN=1.E-16_EB,X_FU_MIN=1.E-16_EB,Y_CO_MIN=1.E-10_EB, &
                       M_MIN=0.1_EB,M_MAX=0.3_EB
INTEGER :: NODETS,I,J,K,II,JJ,KK,IOR,IC,IW,IWA(-3:3),ITMP,ICFT
REAL(EB), POINTER, DIMENSION(:,:,:) :: Y_O2=>NULL(),Y_O2_NEW=>NULL(),MIX_TIME=>NULL(), &
                                       UU=>NULL(),VV=>NULL(),WW=>NULL()

! Misc initializations

Y_O2     => WORK1
MIX_TIME => WORK2 ! need an array because MIX_TIME is reused in CO_PRODUCTION
!$OMP PARALLEL SHARED(Y_O2,RSUM)
!$OMP WORKSHARE
Y_O2     =  0._EB
MIX_TIME =  DT
Q        =  0._EB
!$OMP END WORKSHARE

! Compute and save O2 in all cells

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,YY_GET)
   DO K=0,KBP1
      DO J=0,JBP1
         DO I=0,IBP1
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            YY_GET(:) = YY(I,J,K,:)
            CALL GET_MASS_FRACTION(YY_GET,O2_INDEX,Y_O2(I,J,K))   
         ENDDO
      ENDDO
   ENDDO
!$OMP END DO
IF (CO_PRODUCTION) THEN
   !$OMP SINGLE
   Y_O2_NEW => WORK3
   !$OMP END SINGLE
   !$OMP WORKSHARE
   Y_O2_NEW =  Y_O2
   !$OMP END WORKSHARE
ENDIF
!$OMP SINGLE

! Compute the (fast) reaction of fuel to either CO or CO2
IF (CO_PRODUCTION) THEN
   DELTAH_F = REACTION(2)%HEAT_OF_COMBUSTION
ELSE
   DELTAH_F = REACTION(1)%HEAT_OF_COMBUSTION
ENDIF


RN => REACTION(1)
HFAC_F  = RN%HEAT_OF_COMBUSTION/DT

IF (NEW_MIX_TIME) THEN
   UU => US
   VV => VS
   WW => WS
ENDIF

!$OMP END SINGLE

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,IC,Y_FU_0,Y_O2_0,Y_O2_MAX,Y_FU_MAX,TMP_MIN,TMP_F_MIN,IOR,II,JJ,KK,IW) &
!$OMP PRIVATE(Y_O2_CORR,Y_F_CORR,DYF,DELTA,Q_BOUND_1,Q_BOUND_2,Q_NEW,CRIT_FLAME_TMPA,YY_GET,DYAIR) &
!$OMP PRIVATE(TAU_D,KP,CP,S_L,TAU_CHEM,DUDX,DUDY,DUDZ,DVDX,DVDY,DVDZ,DWDX,DWDY,DWDZ,SS2,S12,S13,S23,TAU_U,TAU_G,EPSK,KSGS)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IC = CELL_INDEX(I,J,K)
         IF (SOLID(IC)) CYCLE
         IF (TMP(I,J,K) < RN%AIT) CYCLE
         Y_FU_0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))*RN%Y_F_INLET
         IF (Y_FU_0<=Y_FU_MIN) CYCLE
         Y_O2_0  = Y_O2(I,J,K)
         IF (Y_O2_0<=Y_O2_MIN) CYCLE
         
         IF_SUPPRESSION: IF (SUPPRESSION) THEN  ! Get maximum O2 in the current and neighboring cells to see if flame viable

            Y_O2_MAX  = 0._EB
            Y_FU_MAX  = 0._EB
            TMP_MIN   = TMPA
            TMP_F_MIN = TMPA
            CRIT_FLAME_TMPA = 20._EB + TMPM

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
            IF (EXTINCTION2) THEN
               DYF = MIN(Y_FU_0,Y_O2_0/RN%O2_F_RATIO) 
               YY_GET = 0._EB
               YY_GET(I_FUEL) = 1._EB
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_F_0,TMP(I,J,K))
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_F_N,RN%CRIT_FLAME_TMP)
               YY_GET = YY(I,J,K,:)
               YY_GET(I_FUEL) = 0._EB
               YY_GET = YY_GET / (1._EB - Y_FU_0)
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_0,TMP(I,J,K))            
               CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,H_G_N,RN%CRIT_FLAME_TMP)
               DYAIR = DYF * (1._EB - Y_FU_0) / Y_O2_0 * RN%O2_F_RATIO
               IF ( (DYF*H_F_0 + DYAIR*H_G_0)*TMP(I,J,K) + DYF*DELTAH_F < (DYF*H_F_N + DYAIR*H_G_N)*RN%CRIT_FLAME_TMP) CYCLE
            ELSE
               Y_O2_CORR = RN%Y_O2_LL*(RN%CRIT_FLAME_TMP-TMP_MIN)/(RN%CRIT_FLAME_TMP-CRIT_FLAME_TMPA)
               Y_F_CORR  = RN%Y_F_LFL*(RN%CRIT_FLAME_TMP-TMP_F_MIN)/(RN%CRIT_FLAME_TMP-CRIT_FLAME_TMPA)
               IF (Y_O2_MAX < Y_O2_CORR .OR. Y_FU_MAX*RN%Y_F_INLET < Y_F_CORR) CYCLE
            ENDIF

         ENDIF IF_SUPPRESSION
         
         LES_IF: IF (LES) THEN
            
            IF (USE_MAX_FILTER_WIDTH) THEN
               DELTA=MAX(DX(I),DY(J),DZ(K))
            ELSE
               IF (.NOT.TWO_D) THEN
                  DELTA = (DX(I)*DY(J)*DZ(K))**ONTH
               ELSE
                  DELTA = SQRT(DX(I)*DZ(K))
               ENDIF
            ENDIF
 
            EXPERIMENTAL_IF: IF (NEW_MIX_TIME) THEN
               ! experimental
               TAU_D = SC*RHO(I,J,K)*DELTA**2/MU(I,J,K)   ! diffusive time scale
               
               ! compute local filtered strain
               DUDX = RDX(I)*(UU(I,J,K)-UU(I-1,J,K))
               DVDY = RDY(J)*(VV(I,J,K)-VV(I,J-1,K))
               DWDZ = RDZ(K)*(WW(I,J,K)-WW(I,J,K-1))
               DUDY = 0.25_EB*RDY(J)*(UU(I,J+1,K)-UU(I,J-1,K)+UU(I-1,J+1,K)-UU(I-1,J-1,K))
               DUDZ = 0.25_EB*RDZ(K)*(UU(I,J,K+1)-UU(I,J,K-1)+UU(I-1,J,K+1)-UU(I-1,J,K-1)) 
               DVDX = 0.25_EB*RDX(I)*(VV(I+1,J,K)-VV(I-1,J,K)+VV(I+1,J-1,K)-VV(I-1,J-1,K))
               DVDZ = 0.25_EB*RDZ(K)*(VV(I,J,K+1)-VV(I,J,K-1)+VV(I,J-1,K+1)-VV(I,J-1,K-1))
               DWDX = 0.25_EB*RDX(I)*(WW(I+1,J,K)-WW(I-1,J,K)+WW(I+1,J,K-1)-WW(I-1,J,K-1))
               DWDY = 0.25_EB*RDY(J)*(WW(I,J+1,K)-WW(I,J-1,K)+WW(I,J+1,K-1)-WW(I,J-1,K-1))
               S12 = 0.5_EB*(DUDY+DVDX)
               S13 = 0.5_EB*(DUDZ+DWDX)
               S23 = 0.5_EB*(DVDZ+DWDY)
               SS2 = 2._EB*(DUDX**2 + DVDY**2 + DWDZ**2 + 2._EB*(S12**2 + S13**2 + S23**2))
               
               EPSK = MU(I,J,K)/RHO(I,J,K)*SS2       ! ke dissipation rate, assumes production=dissipation
               KSGS = 2.25_EB*(EPSK*DELTA/PI)**TWTH  ! estimate of subgrid ke, from Kolmogorov spectrum

               TAU_U = DELTA/SQRT(2._EB*KSGS+1.E-10_EB)   ! advective time scale
               TAU_G = SQRT(2._EB*DELTA/(GRAV+1.E-10_EB)) ! acceleration time scale
               MIX_TIME(I,J,K)=MIN(TAU_D,TAU_U,TAU_G)
            ELSE EXPERIMENTAL_IF
               ! FDS 5 default
               MIX_TIME(I,J,K)=C_EDC*SC*RHO(I,J,K)*DELTA**2/MU(I,J,K)
            ENDIF EXPERIMENTAL_IF
            
         ENDIF LES_IF
         
         ! chemical time scale
         CHEM_IF: IF (CHECK_CHEMICAL_TIME_SCALE) THEN
            TAU_CHEM = 0._EB ! infinitely fast chemistry
            YY_GET(:) = YY(I,J,K,:)
            CALL GET_CONDUCTIVITY(YY_GET,KP,TMP(I,J,K))
            CALL GET_AVERAGE_SPECIFIC_HEAT(YY_GET,CP,TMP(I,J,K))  
            S_L = LAMINAR_FLAME_SPEED(TMPA,RN%O2_F_RATIO/(Y_O2_0/Y_FU_0))
            IF (S_L>0._EB) TAU_CHEM = KP/(RHO(I,J,K)*CP*S_L*S_L)
            MIX_TIME(I,J,K)=MAX(TAU_CHEM,MIX_TIME(I,J,K))           
         ENDIF CHEM_IF
         
         IF (FIXED_MIX_TIME>0._EB) MIX_TIME(I,J,K)=FIXED_MIX_TIME
         
         NEW_DYF_IF: IF (NEW_MIX_TIME) THEN
            IF (Y_FU_0 < Y_O2_0/RN%O2_F_RATIO) THEN
               DYF = Y_FU_0 * (1._EB -EXP(-DT/MIX_TIME(I,J,K)))
            ELSE
               DYF = Y_O2_0/RN%O2_F_RATIO * (1._EB -EXP(-DT/MIX_TIME(I,J,K)))
            ENDIF
            Q_BOUND_1 = DYF*RHO(I,J,K)*HFAC_F
         ELSE NEW_DYF_IF
            ! FDS 5 default
            DYF = MIN(Y_FU_0,Y_O2_0/RN%O2_F_RATIO)
            Q_BOUND_1 = DYF*RHO(I,J,K)*HFAC_F*MIN(1._EB,DT/MIX_TIME(I,J,K))
         ENDIF NEW_DYF_IF
         
         Q_BOUND_2 = Q_UPPER
         Q_NEW = MIN(Q_BOUND_1,Q_BOUND_2)
         DYF = Q_NEW /(RHO(I,J,K)*HFAC_F*RN%Y_F_INLET)
         
         Q(I,J,K)  = Q_NEW
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
!$OMP END DO

! Optional second (slow) reaction to convert CO to CO_2
CONVERT_CO: IF (CO_PRODUCTION) THEN 
   !$OMP SINGLE
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
   !$OMP END SINGLE

   !$OMP DO COLLAPSE(3) PRIVATE(K,J,I,Y_O2_0,Y_CO_0,Z_2_MIN,DYCO,RHOX,X_FU_0,X_O2_0,X_FU,X_O2,ETRM) &
   !$OMP PRIVATE(II,DX_FDT,X_FU_S,X_O2_S,X_FU_N,X_O2_N,Q_OLD,Q_NEW)
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
            IF ((TMP(I,J,K) < 500._EB .AND. Q(I,J,K)==0._EB) .OR. Q(I,J,K)>=Q_UPPER) CYCLE

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
   !$OMP END DO

ENDIF CONVERT_CO

! Compute new mixture molecular weight

!$OMP DO COLLAPSE(3) PRIVATE(YY_GET)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         YY_GET(:) = YY(I,J,K,:)
         CALL GET_SPECIFIC_GAS_CONSTANT(YY_GET,RSUM(I,J,K))
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

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

!$OMP PARALLEL
!$OMP WORKSHARE
Q =  0._EB
!$OMP END WORKSHARE

!$OMP SINGLE
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
!$OMP END SINGLE

!$OMP DO COLLAPSE(3) PRIVATE(K,J,I,DYF,TMPD,NR,ETRM1,NO_REACTION,II,X_O2_N,X_FU_N,Q_NR,RN,X_O2,X_FU,ETRM,N) &
!$OMP PRIVATE(DX_FDT,X_FU_S,X_O2_S,QFAC,YYNEW)
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
!$OMP END DO

! Adjust the average molecular weight term, R*Sum(Yi/Mi)

!$OMP WORKSHARE
RSUM = SPECIES(0)%RCON
!$OMP END WORKSHARE
!$OMP END PARALLEL
SLOOP: DO N=1,N_SPECIES
   IF (SPECIES(N)%MODE/=GAS_SPECIES) CYCLE SLOOP         
   WFAC = SPECIES(N)%RCON - SPECIES(0)%RCON
   !$OMP PARALLEL WORKSHARE
   RSUM(:,:,:) = RSUM(:,:,:) + WFAC*YY(:,:,:,N)
   !$OMP END PARALLEL WORKSHARE
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
TYPE (REAC_SPECIES_TYPE), POINTER :: RS=>NULL()
TYPE (REACTION_TYPE), POINTER :: RN=>NULL()

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
TYPE (REACTION_TYPE), POINTER :: RN=>NULL()

DO NR = 1, N_REACTIONS
   RN => REACTION(NR)
   KR(NR) = RN%BOF*EXP(-RN%E/(R0*TMP))
END DO

END SUBROUTINE RATE_CONSTANT

SUBROUTINE COMBUSTION_BC

INTEGER :: II,JJ,KK,IIG,JJG,KKG,IW

!$OMP PARALLEL DO PRIVATE(II,JJ,KK,IIG,JJG,KKG)
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
!$OMP END PARALLEL DO

END SUBROUTINE COMBUSTION_BC


REAL(EB) FUNCTION LAMINAR_FLAME_SPEED(TMP_U,EQ_RAT)
IMPLICIT NONE

REAL(EB), INTENT(IN) :: TMP_U     ! temperature of unburned gases
REAL(EB), INTENT(IN) :: EQ_RAT    ! equivalence ratio (F/A)/(F/A)_stoic
REAL(EB), PARAMETER :: TMP_REF=298._EB
REAL(EB) :: S_L_REF,GAMMA,PHI

! Reference:
! Stephen Turns, 'An Introduction to Combustion: Concepts and Applications', Chapter 8
! see Eq. (8.33)
!
! Metghalachi and Keck, Burning Velocities at High P and T, Combustion and Flame, 48:191-210 (1982).
!
! Comments:
! For now, we base S_L on the data for Propane. Pressue is assumed to be 1 atm.

! data for Propane from Turns
! ---------------------------------------------
REAL(EB), PARAMETER :: PHI_M = 1.08_EB
REAL(EB), PARAMETER :: B_M   = 0.3422_EB  ! m/s
REAL(EB), PARAMETER :: B2    = -1.3865_EB ! m/s
REAL(EB), PARAMETER :: PHI_MIN = 0.8_EB ! lower limit of M&K data
REAL(EB), PARAMETER :: PHI_MAX = 1.4_EB ! upper limit of M&K data
! ---------------------------------------------

PHI = MIN(MAX(PHI_MIN,EQ_RAT),PHI_MAX)
S_L_REF = B_M + B2*(PHI-PHI_M)**2
GAMMA = 2.18_EB - 0.8_EB*(PHI-1._EB)

LAMINAR_FLAME_SPEED = S_L_REF*(MAX(TMP_U,TMP_REF)/TMP_REF)**GAMMA

END FUNCTION LAMINAR_FLAME_SPEED


END SUBROUTINE COMBUSTION



SUBROUTINE GET_REV_fire(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') firerev(INDEX(firerev,':')+1:LEN_TRIM(firerev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') firedate

END SUBROUTINE GET_REV_fire
 
END MODULE FIRE

