MODULE FIRE
 
! Compute combustion 
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: SECOND
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: fireid='$Id: fire.f90 3764 2009-04-13 19:04:27Z mcgratta $'
CHARACTER(255), PARAMETER :: firerev='$Revision: 3764 $'
CHARACTER(255), PARAMETER :: firedate='$Date: 2009-04-13 15:04:27 -0400 (Mon, 13 Apr 2009) $'

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

Q_UPPER = HRRPUA_SHEET/CELL_SIZE

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

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_SPECIFIC_GAS_CONSTANT,GET_SPECIFIC_ENTHALPY
REAL(EB) :: A,ETRM,DYF,DYFO2,DX_FDT,DTT,DELTA2,& 
            Y_O2_MAX,TMP_MIN,Y_O2_CORR,Q_NEW,Q_OLD,DYCO,RHOX, &
            DELTAH_F,DELTAH_CO,HFAC_F,HFAC_CO, &
            F_TO_CO,F_TO_CO2,CO_TO_O2,CO_CO2,OM_CO_CO2, &
            Y_FU_MAX,TMP_F_MIN,Y_F_CORR,Z_2_MIN,Z_2_MIN_FAC,WGT,OMWGT,Q_BOUND_1,Q_BOUND_2,YY_GET(1:N_SPECIES),&
            H_GAS,H_LHS,H_CFL,RHO_G, &
            Y_FU_0,Y_O2_0,Y_CO_0,Y_CO2_0,Y_CO_NC_0, &
            X_FU_0,X_O2_0,X_CO_0,X_CO2_0,X_CO_NC_0,X_FU,X_O2,X_CO,X_CO2,X_CO_NC,X_FU_1,X_O2_1,X_CO_1,X_CO2_1,X_CO_NC_1, &
            R11,R12,R13,R21,R22,R23,QTEMP,HDTT
REAL(EB), PARAMETER :: Y_FU_MIN=1.E-10_EB,Y_O2_MIN=1.E-10_EB,Y_CO_MIN=1.E-10_EB, &
                       X_FU_MIN=1.E-16_EB,X_O2_MIN=1.E-16_EB,X_CO_MIN=1.E-16_EB
INTEGER :: NODETS,I,J,K,II,JJ,KK,IC,ITMP,N
REAL(EB), POINTER, DIMENSION(:,:,:) :: Y_O2,MIX_TIME
TYPE(REACTION_TYPE), POINTER,DIMENSION(:) :: RN
! Weighting factor for time-averaging of local HRR

WGT   = MIN(1._EB,DT/HRRPUV_AVERAGE_TIME)
OMWGT = 1._EB - WGT

! Misc initializations

Y_O2     => WORK1
Y_O2     =  0._EB
MIX_TIME => WORK2
MIX_TIME =  DT / 24._EB
Q        =  0._EB
RN       => REACTION
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

! Compute the (fast) reaction of fuel to either CO or CO2

HFAC_F  = RN(1)%HEAT_OF_COMBUSTION/DT

makeco: IF (.NOT. CO_PRODUCTION) THEN
   DO K=1,KBAR
      DO J=1,JBAR
         singlestep: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE singlestep
            Y_FU_0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))*RN(1)%Y_F_INLET
            IF (Y_FU_0<=Y_FU_MIN) CYCLE singlestep
            Y_O2_0  = Y_O2(I,J,K)
            IF (Y_O2_0<=Y_O2_MIN) CYCLE singlestep
            DYFO2 = Y_O2_0/RN(1)%O2_F_RATIO
            ! Evaluate empirical extinction criteria (Jukka's idea)
            IF (SUPPRESSION) THEN
               RHO_G = RHO(I,J,K)   
               ITMP = NINT(TMP(I,J,K))
               YY_GET = YY(I,J,K,:)
               CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_GAS,ITMP)
               ITMP = NINT(RN(1)%CRIT_FLAME_TMP)
               IF (DYFO2 > Y_FU_0) THEN
                  YY_GET(I_FUEL) = YY_GET(I_FUEL) + YY(I,J,K,I_FUEL)- Y_FU_0
                  YY_GET(I_PROG_F) = YY_GET(I_PROG_F) + Y_FU_0
                  CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_CFL,ITMP)
                  H_LHS = (1._EB+Y_FU_0/DYFO2*(1._EB-YY(I,J,K,I_FUEL)))*(H_CFL-H_GAS)
               ELSE
                  YY_GET(I_FUEL) = YY_GET(I_FUEL) + YY(I,J,K,I_FUEL) - DYFO2
                  YY_GET(I_PROG_F) = YY_GET(I_PROG_F) + DYFO2
                  CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_CFL,ITMP)
                  H_LHS = (1._EB+1._EB/DYFO2*(1._EB-YY(I,J,K,I_FUEL)))*(H_CFL-H_GAS)
               ENDIF
               IF (H_LHS > RN(1)%HEAT_OF_COMBUSTION) THEN
                  CYCLE singlestep
               ENDIF
            ENDIF 
            IF (LES .AND. EDDY_DISSIPATION) THEN
               IF (.NOT.TWO_D) THEN
                  DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
               ELSE
                  DELTA2 = DX(I)*DZ(K)
               ENDIF
               MIX_TIME(I,J,K) = C_EDC*SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
            ENDIF
            IF (Y_FU_0 < DYFO2) THEN
               DYF = Y_FU_0 * (1._EB - EXP(-DT/MIX_TIME(I,J,K)))
            ELSE
               DYF = DYFO2  * (1._EB - EXP(-DT/MIX_TIME(I,J,K)))
            ENDIF
            Q_BOUND_1 = DYF*RHO(I,J,K)*HFAC_F
            Q_BOUND_2 = Q_UPPER
            IF (LES) Q_BOUND_2 = MIN(Q_BOUND_2,(HRRPUV_AVERAGE-OMWGT*Q_AVG(I,J,K))/WGT)
            Q_NEW = MIN(Q_BOUND_1,Q_BOUND_2)

            DYF = Q_BOUND_1/(Q_NEW*RN(1)%Y_F_INLET)
            Q(I,J,K)   = Q_NEW
            YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
            IF (SOOT_DEPOSITION) THEN
               YY(I,J,K,I_PROG_F)    = YY(I,J,K,I_PROG_F)    + DYF * (1._EB - RN(1)%SOOT_YIELD)
               YY(I,J,K,I_PROG_SOOT) = YY(I,J,K,I_PROG_SOOT) + DYF * RN(1)%SOOT_YIELD
            ELSE
               YY(I,J,K,I_PROG_F)  = YY(I,J,K,I_PROG_F)  + DYF            
            ENDIF
         ENDDO singlestep
      ENDDO
   ENDDO 

ELSE makeco
   
   DELTAH_F  = RN(1)%HEAT_OF_COMBUSTION * RN(1)%MW_FUEL / 1000._EB !J/mol
   DELTAH_CO = (CO2_HEAT_OF_FORMATION - CO_HEAT_OF_FORMATION) !J/mol
!   WRITE(*,*) DELTAH_F,DELTAH_CO,RN(1)%HEAT_OF_COMBUSTION,RN(2)%HEAT_OF_COMBUSTION
   F_TO_CO   = RN(1)%MW_FUEL/(RN(1)%NU(CO_INDEX)*MW_CO)  
   F_TO_CO2  = RN(1)%MW_FUEL/(RN(1)%NU(CO_INDEX)*MW_CO2) 
   HFAC_CO   = DELTAH_CO/DT
   CO_TO_O2  = MW_O2/(MW_CO*2._EB)   
   Z_2_MIN_FAC = F_TO_CO * RN(2)%CO_YIELD / RN(1)%Y_F_INLET
   A  = RN(2)%BOF 
   NODETS = 24
   DTT    = DT/REAL(NODETS,EB)
!   WRITE(*,*) DT, ' DT'
   HDTT   = 0.5_EB * DTT
   Q_BOUND_2 = Q_UPPER
   CO_CO2 = RN(2)%NU(CO2_INDEX)/RN(1)%NU(CO_INDEX)
   OM_CO_CO2 = 1._EB - CO_CO2   
   DO K=1,KBAR
      DO J=1,JBAR
         twostep: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE twostep
            Y_O2_0  = Y_O2(I,J,K)
            IF (Y_O2_0<=Y_O2_MIN) CYCLE twostep
            DYFO2 = Y_O2_0/RN(2)%O2_F_RATIO
            Y_FU_0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))*RN(1)%Y_F_INLET
            Y_CO_0  = MAX(0._EB,YY(I,J,K,I_PROG_CO))*RN(1)%Y_F_INLET / F_TO_CO 
            IF (Y_FU_0<=Y_FU_MIN .AND. Y_CO_0 <= Y_CO_MIN) CYCLE twostep
            ! Evaluate empirical extinction criteria (Jukka's idea)
            RHO_G = RHO(I,J,K)   
            ITMP = NINT(TMP(I,J,K))
            YY_GET = YY(I,J,K,:)
            CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_GAS,ITMP)
            ITMP = NINT(RN(1)%CRIT_FLAME_TMP)
            IF (DYFO2 > Y_FU_0) THEN
               YY_GET(I_FUEL) = YY_GET(I_FUEL) + YY(I,J,K,I_FUEL)- Y_FU_0
               YY_GET(I_PROG_F) = YY_GET(I_PROG_F) + Y_FU_0
               CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_CFL,ITMP)
               H_LHS = (1._EB+Y_FU_0/DYFO2*(1._EB-YY(I,J,K,I_FUEL)))*(H_CFL-H_GAS)
            ELSE
               YY_GET(I_FUEL) = YY_GET(I_FUEL) + YY(I,J,K,I_FUEL) - DYFO2
               YY_GET(I_PROG_F) = YY_GET(I_PROG_F) + DYFO2
               CALL GET_SPECIFIC_ENTHALPY(YY_GET,H_CFL,ITMP)
               H_LHS = (1._EB+1._EB/DYFO2*(1._EB-YY(I,J,K,I_FUEL)))*(H_CFL-H_GAS)
            ENDIF
     
            IF (H_LHS > RN(2)%HEAT_OF_COMBUSTION) THEN 
!IF(J==20 .AND. K==1)              WRITE(*,*) I,J,K,NINT(TMP(I,J,K)),ITMP,Y_FU_0,Y_O2_0,H_LHS
!IF(J==20 .AND. K==1)              WRITE(*,*) YY(I,J,K,:)
               IF (TMP(I,J,K) < 500._EB .OR. Y_CO_0<= Z_2_MIN_FAC * (Y_CO_0+YY(I,J,K,I_PROG_F))) CYCLE twostep
               MIX_TIME(I,J,K) = -1._EB
               ETRM = EXP(-RN(2)%E/(R0*TMP(I,J,K))) 
               Y_FU_0=0._EB
            ELSE
               IF (LES) THEN
                  IF (.NOT.TW    O_D) THEN
                     DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
                  ELSE
                     DELTA2 = DX(I)*DZ(K)
                  ENDIF
                  MIX_TIME(I,J,K) = C_EDC*SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
               ENDIF
            ENDIF
            Y_CO2_0  = MAX(0._EB,YY(I,J,K,I_PROG_F))*RN(1)%Y_F_INLET/F_TO_CO2
            Y_CO_NC_0 = Z_2_MIN_FAC*(YY(I,J,K,I_PROG_CO)+YY(I,J,K,I_PROG_F))*RN(1)%Y_F_INLET/F_TO_CO
            Y_CO_0 = Y_CO_0 - Y_CO_NC_0
!IF(J==20 .AND. K==1)            WRITE(*,*) YY(I,J,K,:),I,J,K,'A'
            IF (LES) Q_BOUND_2 = MIN((HRRPUV_AVERAGE-OMWGT*Q_AVG(I,J,K))/WGT,Q_BOUND_2)
            RHOX = 1000._EB*RHO(I,J,K)                  
            X_FU_0 = Y_FU_0*RHOX/RN(1)%MW_FUEL*1.E-6_EB
            X_O2_0 = Y_O2_0*RHOX/MW_O2*1.E-6_EB
            X_CO_0 = Y_CO_0*RHOX/MW_CO*1.E-6_EB
            X_CO_NC_0 = Y_CO_NC_0*RHOX/MW_CO*1.E-6_EB
            X_CO2_0 = Y_CO2_0*RHOX/MW_CO2*1.E-6_EB
            X_FU=X_FU_0
            X_O2=X_O2_0
            X_CO=X_CO_0
            X_CO_NC=X_CO_NC_0
            X_CO2=X_CO2_0
!IF(J==20 .AND. K==1) WRITE(*,'(5(E14.8,1X))') X_FU,X_O2,X_CO,X_CO_NC,X_CO2            
!IF(J==20 .AND. K==1) WRITE(*,'(1X,4(E17.10,1X))') RHOX,MIX_TIME(I,J,K)
            Q_NEW = 0._EB
            DYF = 0._EB
            DYCO = 0._EB
            co_rk: DO N=1,NODETS
!RK Step 1
              IF (MIX_TIME(I,J,K)<0._EB) THEN
                 R11 = 0._EB
                 R21 = -A*ETRM*(X_CO+X_CO_NC_1)*X_O2
              ELSE
                 R11 = RN(1)%NU(CO_INDEX) * MIN(X_FU,X_O2/RN(1)%NU(O2_INDEX)) / MIX_TIME(I,J,K)
                 R21 = -1._EB*MIN(X_CO+X_CO_NC,2._EB*X_O2) / MIX_TIME(I,J,K) * 0.2_EB
              ENDIF
              X_FU_1 = X_FU - 1._EB/RN(1)%NU(CO_INDEX) * R11* HDTT
              X_O2_1 = X_O2 - (RN(1)%NU(O2_INDEX)/RN(1)%NU(CO_INDEX) * R11 - 0.5_EB*CO_CO2*R21)* HDTT 
              X_CO_1 = X_CO + (R11+R21) * HDTT
              X_CO_NC_1 = X_CO_NC - OM_CO_CO2*R21 * HDTT
!IF(J==20 .AND. K==1) WRITE(*,'(4(E14.8,1X))') X_FU_1,X_O2_1,X_CO_1,X_CO_NC_1
!RK Step 2
              IF (MIX_TIME(I,J,K)<0._EB) THEN
                 R12 = 0._EB
                 R22 = -A*ETRM*(X_CO_1+X_CO_NC_1)*X_O2_1
              ELSE
                 R12 = RN(1)%NU(CO_INDEX) * MIN(X_FU_1,X_O2_1/RN(1)%NU(O2_INDEX)) / MIX_TIME(I,J,K)
                 R22 = -1._EB*MIN(X_CO_1+X_CO_NC_1,2._EB*X_O2_1) / MIX_TIME(I,J,K) * 0.2_EB
              ENDIF
              X_FU_1 = X_FU_1 - 1._EB/RN(1)%NU(CO_INDEX) * (2._EB*R12-R11)* DTT
              X_O2_1 = X_O2_1 - (RN(1)%NU(O2_INDEX)/RN(1)%NU(CO_INDEX) * (2._EB*R12-R11) - 0.5_EB*CO_CO2*(2._EB*R22-R21))* DTT 
              X_CO_1 = X_CO_1 + (2._EB*R12-R11+(2._EB*R22-R21)) * DTT
              X_CO_NC_1 = X_CO_NC_1 - OM_CO_CO2*(2*R22-R21) * DTT
!IF(J==20 .AND. K==1) WRITE(*,'(4(E17.11,1X))') X_FU_1,X_O2_1,X_CO_1,X_CO_NC_1
 !RK Step 3
              IF (MIX_TIME(I,J,K)<0._EB) THEN
                 R13 = 0._EB
                 R23 = -A*ETRM*(X_CO_1+X_CO_NC_1)*X_O2_1
              ELSE
                 R13 = RN(1)%NU(CO_INDEX) * MIN(X_FU_1,X_O2_1/RN(1)%NU(O2_INDEX)) / MIX_TIME(I,J,K)
                 R23 = -1._EB*MIN(X_CO_1+X_CO_NC_1,2._EB*X_O2_1) / MIX_TIME(I,J,K) * 0.2_EB
              ENDIF
!IF(J==20 .AND. K==1) WRITE(*,*) R11,R12,R13,' 1'
!IF(J==20 .AND. K==1) WRITE(*,*) R21,R22,R23,' 2'
              R11 = ONSI * (R11 + 4._EB*R12 + R13) * DTT
              R21 = ONSI * (R21 + 4._EB*R22 + R23) * DTT
!IF(J==20 .AND. K==1) WRITE(*,*) R11/DTT,R21/DTT,' 3'
              QTEMP = (R11/RN(1)%NU(CO_INDEX) * DELTAH_F - CO_CO2 * R21 * DELTAH_CO)*1.E6_EB
              IF ((Q_NEW + QTEMP)/DT > Q_BOUND_2) THEN
                 R11 = DT*(Q_BOUND_2-Q_NEW/DT)/QTEMP * R11
                 R21 = DT*(Q_BOUND_2-Q_NEW/DT)/QTEMP * R21
!  WRITE(*,*) Q_BOUND_2,Q_NEW,QTEMP
!  WRITE(*,*) R11/DTT,R21/DTT
              ENDIF
              IF (X_CO+R11 + R21 < 0._EB) &
                 R21 = -(X_CO+R11)
              QTEMP = (R11/RN(1)%NU(CO_INDEX) * DELTAH_F - CO_CO2 * R21 * DELTAH_CO)*1.E6_EB
              X_FU = X_FU  - 1._EB/RN(1)%NU(CO_INDEX) *  R11 
              X_O2 = X_O2  - (RN(1)%NU(O2_INDEX)/RN(1)%NU(CO_INDEX) * R11 - 0.5_EB*CO_CO2*R21) 
              X_CO = X_CO + R11 + R21
              X_CO_NC = X_CO_NC - OM_CO_CO2*R21
              X_CO2 = X_CO2 - CO_CO2*R21
!IF(J==20 .AND. K==1) WRITE(*,'(5(E14.8,1X))') X_FU,X_O2,X_CO,X_CO_NC,X_CO2                          
!IF(J==20 .AND. K==1) WRITE(*,*) 'BBB',QTEMP,Q_NEW,Q_NEW/DT,Q_BOUND_2
              Q_NEW = Q_NEW + QTEMP
              IF (Q_NEW >= Q_BOUND_2) EXIT co_rk
              IF (X_FU < X_FU_MIN .AND. X_CO < X_CO_MIN) EXIT co_rk
            ENDDO co_rk
!IF(J==20 .AND. K==1)           WRITE(*,*) Q_NEW
            Q(I,J,K)   = Q_NEW / DT 
!IF(J==20 .AND. K==1) WRITE(*,*) I,J,K,YY(I,J,K,:)
            YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL)-(X_FU_0-X_FU)*RN(1)%MW_FUEL/RHOX*1.E6/RN(1)%Y_F_INLET                
            YY(I,J,K,I_PROG_CO) = (X_CO+X_CO_NC)*MW_CO/RHOX*1.E6/RN(1)%Y_F_INLET*F_TO_CO
            YY(I,J,K,I_PROG_F) = X_CO2*MW_CO2/RHOX*1.E6/RN(1)%Y_F_INLET*F_TO_CO2
!IF(J==20 .AND. K==1)           WRITE(*,*) YY(I,J,K,:),'B'
            !IF (SOOT_DEPOSITION) THEN
            !   YY(I,J,K,I_PROG_CO)   = YY(I,J,K,I_PROG_CO)   + DYF * (1._EB - RN(1)%SOOT_YIELD)
            !   YY(I,J,K,I_PROG_SOOT) = YY(I,J,K,I_PROG_SOOT) + DYF * RN(1)%SOOT_YIELD
            !ELSE
            !   YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) + DYF            
            !ENDIF
         ENDDO twostep
      ENDDO
   ENDDO
ENDIF makeco

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

RETURN

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

