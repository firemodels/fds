MODULE FIRE
 
! Compute combustion 
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY: SECOND
 
IMPLICIT NONE
PRIVATE
CHARACTER(255), PARAMETER :: fireid='$Id: fire.f90,v 1.11 2007/04/09 12:15:51 dr_jfloyd Exp $'
TYPE(REACTION_TYPE), POINTER :: RN

PUBLIC COMBUSTION
 
CONTAINS
 
SUBROUTINE COMBUSTION(NM)

INTEGER, INTENT(IN) :: NM
REAL(EB) :: TNOW

IF( EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()

CALL POINT_TO_MESH(NM)

IF (MIXTURE_FRACTION) THEN
   CALL COMBUSTION_MF(NM)
ELSE
   CALL COMBUSTION_FR(NM)
ENDIF

TUSED(10,NM)=TUSED(10,NM)+SECOND()-TNOW

CONTAINS

SUBROUTINE COMBUSTION_MF(NM)
USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_MOLECULAR_WEIGHT
REAL(EB) :: TMPD,YFU0,A,ETRM,YCOMIN,YO2MIN,YFUMIN,YO20,YCO0,&
            DYF,DX_FDT,HFAC_F,DTT,& 
            YYNEW,Z_F,Y_O2_MAX,TMP_MIN,Y_O2_CORR, WFAC,Q_NEW,Q_OLD,&
            F_TO_CO,DELTAH_CO,DYCO,HFAC_CO,Z_2,RHOX, &
            X_FU,X_O,X_FU_0,X_O_0,X_FU_S,X_O_S,X_FU_N,X_O_N,X_O_MIN,X_FU_MIN,COTOO2
REAL(EB), DIMENSION(0:N_SPECIES,1:N_REACTIONS) :: MPUE
INTEGER :: NODETS,N,I,J,K,II,JJ,KK,NR,SCAN_DISTANCE,REAC_COUNT,IC
LOGICAL :: BURN
INTEGER, INTENT(IN) :: NM
REAL(EB), POINTER, DIMENSION(:,:,:) :: YO2Z,QT,R_SUM_DILUENTS
LOGICAL, POINTER, DIMENSION(:,:,:) :: IGNITE
PRODUCE_CO: IF (.NOT. CO_PRODUCTION) THEN !Combustion without CO formation and destruction
   YO2Z => WORK1
   Q =  0._EB
   Z_2 = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         ILOOPA: DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOPA
            CALL GET_MASS_FRACTION (YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),O2_INDEX,Y_SUM(I,J,K),Z_F,YO2Z(I,J,K))   
         ENDDO ILOOPA
      ENDDO
   ENDDO
   YO2MIN = 1.E-10_EB
   YFUMIN = 1.E-10_EB
   !Loop and do F -> CO
   RN => REACTION(1)
   HFAC_F   = RN%HEAT_OF_COMBUSTION/DT
   SUPPRESSIONIF: IF (SUPPRESSION) THEN !TMP and Y_O2 dependent combustion
      SCAN_DISTANCE=1  
      DO K=1,KBAR
         DO J=1,JBAR
            ILOOPB: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOPB
               YO20  = YO2Z(I,J,K)
               YFU0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))
               IF (YFU0<=YFUMIN .OR. YO20<=YO2MIN) CYCLE ILOOPB
               !Get min O2
               Y_O2_MAX = 0._EB   
               BURN = .FALSE.
               BURN_LOOP: DO KK=MAX(1,K-SCAN_DISTANCE),MIN(KBAR,K+SCAN_DISTANCE)
                  DO JJ=MAX(1,J-SCAN_DISTANCE),MIN(JBAR,J+SCAN_DISTANCE)
                     DO II=MAX(1,I-SCAN_DISTANCE),MIN(IBAR,I+SCAN_DISTANCE)
                        IF (YO2Z(II,JJ,KK)>Y_O2_MAX .AND. .NOT.SOLID(CELL_INDEX(II,JJ,KK)) ) THEN
                           Y_O2_MAX = YO2Z(II,JJ,KK)
                           TMP_MIN = TMP(II,JJ,KK)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO BURN_LOOP
               Y_O2_CORR = RN%Y_O2_LL*(RN%CRIT_FLAME_TMP-TMP_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
               IF (Y_O2_MAX<Y_O2_CORR) CYCLE ILOOPB
               DYF = MIN(YFU0,YO20/RN%O2_F_RATIO)
               Q_NEW = MIN(Q_UPPER,DYF*RHO(I,J,K)*HFAC_F)
               DYF = Q_NEW /(RHO(I,J,K)*HFAC_F)
               Q(I,J,K) = Q_NEW
               YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
               YY(I,J,K,I_PROG_F) = YY(I,J,K,I_PROG_F) + DYF
            ENDDO ILOOPB
         ENDDO
      ENDDO
      
   ELSE SUPPRESSIONIF !No suppression
      DO K=1,KBAR
         DO J=1,JBAR
            ILOOPC: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOPC
               YO20  = YO2Z(I,J,K)
               YFU0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))
               IF (YFU0<=YFUMIN .OR. YO20<=YO2MIN) CYCLE ILOOPC
               DYF = MIN(YFU0,YO20/RN%O2_F_RATIO)
               Q_NEW = MIN(Q_UPPER,DYF*RHO(I,J,K)*HFAC_F)
               DYF = Q_NEW /(RHO(I,J,K)*HFAC_F)
               Q(I,J,K) = Q_NEW
               YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
               YY(I,J,K,I_PROG_F) = YY(I,J,K,I_PROG_F) + DYF
            ENDDO ILOOPC
         ENDDO
      ENDDO
   ENDIF SUPPRESSIONIF
       
ELSE PRODUCE_CO !Combustion with suppression and CO production
   YO2Z   => WORK1
!   IGNITE => LOGICAL_WORK
!   IGNITE = .FALSE.
   QT     => WORK3
   DO K=1,KBAR
      DO J=1,JBAR
         ILOOPX: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE ILOOPX
!            IF (Q(I,J,K) > 0._EB) IGNITE(I,J,K) = .TRUE.
!            IF (MAXVAL(WALL_INDEX(IC,:)) > 0) THEN
!               DO N = -3,3
!                  IF (WALL_INDEX(IC,N) > 0) THEN
!                     IF (SURFACE(IJKW(5,WALL_INDEX(IC,N)))%HRRPUA>0._EB)  IGNITE(I,J,K) = .TRUE.
!                  ENDIF
!               END DO
!            ENDIF
            CALL GET_MASS_FRACTION(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),O2_INDEX,Y_SUM(I,J,K),Z_F,YO2Z(I,J,K))
         ENDDO ILOOPX
      ENDDO
   ENDDO
   QT = 0._EB
   Q =  0._EB
   RN => REACTION(1)
   HFAC_F   = RN%HEAT_OF_COMBUSTION/DT
   F_TO_CO = RN%MW_FUEL/(RN%NU_CO*MW_CO)  
   YO2MIN = 1.E-10_EB
   YFUMIN = 1.E-10_EB
   !Loop and do F -> CO
   SCAN_DISTANCE=1
   DO K=1,KBAR
      DO J=1,JBAR
         ILOOPY: DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOPY
            YO20  = YO2Z(I,J,K)
            YFU0  = MIN(1._EB,MAX(0._EB,YY(I,J,K,I_FUEL)))
            IF (YFU0<=YFUMIN .OR. YO20<=YO2MIN) CYCLE ILOOPY
            !Get min O2
               Y_O2_MAX = 0._EB
               DO KK=MAX(1,K-SCAN_DISTANCE),MIN(KBAR,K+SCAN_DISTANCE)
                  DO JJ=MAX(1,J-SCAN_DISTANCE),MIN(JBAR,J+SCAN_DISTANCE)
                     DO II=MAX(1,I-SCAN_DISTANCE),MIN(IBAR,I+SCAN_DISTANCE)
                        IF (.NOT. SOLID(CELL_INDEX(II,JJ,KK))) THEN
!                           IF (IGNITE(II,JJ,KK)) BURN = .TRUE.
                           IF (YO2Z(II,JJ,KK)>Y_O2_MAX) THEN
                              Y_O2_MAX = YO2Z(II,JJ,KK)
                              TMP_MIN = TMP(II,JJ,KK)
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
               Y_O2_CORR = RN%Y_O2_LL*(RN%CRIT_FLAME_TMP-TMP_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
!               IF (Y_O2_MAX<Y_O2_CORR .OR. (TMP(I,J,K) < 700._EB .AND. .NOT. BURN)) CYCLE ILOOPY
               IF (Y_O2_MAX<Y_O2_CORR) CYCLE ILOOPY
               DYF = MIN(YFU0,YO20/RN%O2_F_RATIO)
            Q_NEW = MIN(Q_UPPER,DYF*RHO(I,J,K)*HFAC_F)
            DYF = Q_NEW /(RHO(I,J,K)*HFAC_F)
            Q(I,J,K) = Q_NEW
            YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
            YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) + DYF
            YO2Z(I,J,K) = YO2Z(I,J,K) - DYF * RN%O2_F_RATIO
         ENDDO ILOOPY
      ENDDO
   ENDDO
   !Loop and do CO -> CO2  
   RN => REACTION(2)
   DELTAH_CO = (REACTION(2)%HEAT_OF_COMBUSTION - REACTION(1)%HEAT_OF_COMBUSTION) * &
               REACTION(1)%MW_FUEL/((REACTION(1)%NU_CO-REACTION(2)%NU_CO)*MW_CO)
   HFAC_CO   = DELTAH_CO/DT
   COTOO2 = MW_O2/MW_O2*0.5_EB
   A  = RN%BOF 
   NODETS = 20
   DTT    = DT/REAL(NODETS,EB)
   X_O_MIN = 1.E-16_EB
   X_FU_MIN = 1.E-16_EB
   YCOMIN = 1.E-10_EB
   DO K=1,KBAR
      DO J=1,JBAR
         ILOOPY1: DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOPY1
            YO20  = YO2Z(I,J,K)
            YCO0  = MAX(0._EB,YY(I,J,K,I_PROG_CO)) / F_TO_CO 
            IF (YCO0<=YCOMIN .OR. YO20<=YO2MIN) CYCLE ILOOPY1
            !Get max CO2
            IF((TMP(I,J,K) < 500._EB .AND. Q(I,J,K)==0._EB) .OR. Q(I,J,K)>=Q_UPPER) CYCLE ILOOPY1
            IF (Q(I,J,K)/=0._EB) THEN
               DYCO = MIN(YCO0,YO20/COTOO2)    
            ELSE
               RHOX = 1000._EB*RHO(I,J,K)
               X_FU_0 = YCO0*RHOX/MW_CO*1.E-6_EB
               X_O_0  = YO20*RHOX/MW_O2*1.E-6_EB
               X_FU  = X_FU_0
               X_O   = X_O_0
               ETRM = EXP(-RN%E/(R0*TMP(I,J,K)))   
               ODE_LOOP2: DO II=1,NODETS
                  IF (X_FU<=X_FU_MIN .OR. X_O<=X_O_MIN) EXIT ODE_LOOP2
                  DX_FDT= -A*ETRM*X_FU*X_O
                  X_FU_S = X_FU + DTT*DX_FDT
                  X_O_S = X_O + DTT*DX_FDT*0.5_EB
                  IF (X_O_S<=X_O_MIN) THEN
                     X_FU = MAX(0.0_EB,X_FU-X_O*2._EB)
                     EXIT ODE_LOOP2
                  ENDIF
                  IF (X_FU_S<=X_FU_MIN) THEN
                     X_FU = X_FU_MIN
                     EXIT ODE_LOOP2
                  ENDIF
                  DX_FDT= -A*ETRM*X_FU_S*X_O_S
                  X_FU_N = .5_EB*(X_FU+X_FU_S+DTT*DX_FDT)
                  X_O_N  = .5_EB*(X_O+X_O_S+DTT*DX_FDT*0.5_EB)
                  IF (X_O_N<=X_O_MIN) THEN
                     X_FU = MAX(0.0_EB,0.5_EB*((X_FU+X_FU_S)-(X_O+X_O_S)*2._EB))
                     EXIT ODE_LOOP2
                  ENDIF
                  IF (X_FU_N<=X_FU_MIN) THEN
                     X_FU = X_FU_MIN
                     EXIT ODE_LOOP2
                  ENDIF
                  X_FU  = X_FU_N
                  X_O   = X_O_N
               ENDDO ODE_LOOP2
               DYCO = MIN(X_FU_0,X_FU_0 - X_FU)*MW_CO/RHOX*1.E6
            ENDIF
            Q_OLD = Q(I,J,K)
            Q_NEW = MIN(Q_UPPER-Q_OLD,DYCO*RHO(I,J,K)*HFAC_CO)
            DYCO = Q_NEW/(RHO(I,J,K)*HFAC_CO)
            Q(I,J,K) = Q_OLD + Q_NEW
            YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) - DYCO * F_TO_CO
            YY(I,J,K,I_PROG_F) = YY(I,J,K,I_PROG_F) + DYCO * F_TO_CO
         ENDDO ILOOPY1
      ENDDO
   ENDDO
ENDIF PRODUCE_CO
IF (N_SPEC_DILUENTS > 0) THEN
   R_SUM_DILUENTS => WORK4
   R_SUM_DILUENTS =  0._EB
   DO N=1,N_SPECIES
      IF (SPECIES(N)%MODE==GAS_SPECIES) R_SUM_DILUENTS(:,:,:) = R_SUM_DILUENTS(:,:,:) + SPECIES(N)%RCON*YYS(:,:,:,N)
   ENDDO
ENDIF
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CO_PRODUCTION) THEN
            Z_2 = YYS(I,J,K,I_PROG_CO)
         ELSE
            Z_2 = 0._EB
         ENDIF
         CALL GET_MOLECULAR_WEIGHT(YYS(I,J,K,I_FUEL),Z_2,YYS(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
         RSUM(I,J,K) = R0/RSUM(I,J,K)
      ENDDO
   ENDDO
ENDDO

IF (N_SPEC_DILUENTS > 0) RSUM = RSUM*(1._EB-Y_SUM) + R_SUM_DILUENTS

!FORALL (I=0:IBP1,J=0:JBP1,K=0:KBP1) TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
!TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

END SUBROUTINE COMBUSTION_MF

SUBROUTINE COMBUSTION_FR(NM)
REAL(EB) :: TMPD,YFU0,ETRM,&
            DX_FDT,DTT,YYNEW,WFAC,Q_NEW,QFAC,RHOX, RATE_SUM,&
            X_O_MIN,X_FU_MIN,EXPON,X_FU_S,X_O_S
REAL(EB), DIMENSION(1:N_REACTIONS,0:N_SPECIES) :: MPUE
REAL(EB), DIMENSION(1:N_REACTIONS) :: HFAC_F, A, ETRM1, X_FU,X_O,X_FU_0,X_O_0,&
           DYF,Q_NR,X_FU_N,X_O_N
INTEGER :: NODETS,N,I,J,K,II,NR,REAC_COUNT
INTEGER, INTENT(IN) :: NM
LOGICAL :: NO_REACTION
Q =  0._EB
DO NR=1,N_REACTIONS
   RN => REACTION(NR)
   HFAC_F(NR) = RN%HEAT_OF_COMBUSTION/DT
   EXPON = 0._EB
   DO N=0,N_SPECIES
      MPUE(NR,N) = SPECIES(N)%MW*RN%NU(N)/ABS(RN%EPUMO2*SPECIES(RN%I_OXIDIZER)%MW*RN%NU(RN%I_OXIDIZER))
   ENDDO
   A(NR) = RN%BOF
ENDDO 
NODETS = 20
DTT    = DT/REAL(NODETS,EB)
X_O_MIN = 1.E-14_EB
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
            X_O_N  = 0._EB
            X_FU_N = 0._EB
            REACTION_LOOP: DO NR=1,N_REACTIONS
               Q_NR(NR) = 0._EB
               RN => REACTION(NR)
               X_O(NR)  = YY(I,J,K,RN%I_OXIDIZER)*1000._EB*RHO(I,J,K)/SPECIES(RN%I_OXIDIZER)%MW*1.E-6_EB
               X_FU(NR) = YY(I,J,K,RN%I_FUEL)*1000._EB*RHO(I,J,K)/RN%MW_FUEL*1.E-6_EB
!               IF(I>=8 .AND. J==1 .AND. K==7 .AND. TMPD>1200.) &
!               WRITE(*,'(1X,I3,4(1X,E14.5))') NR,X_O(NR),X_FU(NR),YY(I,J,K,RN%I_OXIDIZER),YY(I,J,K,RN%I_FUEL)
               IF (X_FU(NR)<=X_FU_MIN .OR. X_O(NR)<=X_O_MIN) THEN
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
!             IF(I>=8 .AND. J==1 .AND. K==7.AND. TMPD>1200.) WRITE(*,*) NR,ETRM,TMPD,'B'
! Solve the simple ODE to deplete fuel and oxidizer due to reaction
               DX_FDT= -ETRM*X_FU(NR)**RN%N_F*X_O(NR)**RN%N_O
               X_FU_S = X_FU(NR) + DTT*DX_FDT
               X_O_S = X_O(NR) + DTT*DX_FDT*RN%N_O/RN%N_F
!               IF(I>=8 .AND. J==1 .AND. K==7) WRITE(*,*) II,NR,NO_REACTION,X_FU_S,'C1'
               IF (X_O_S<X_O_MIN) THEN
                  X_O_S = X_O_MIN 
                  X_FU_S = MAX(0.0_EB,X_FU(NR)-(X_O(NR)-X_O_S)*RN%NU(RN%I_FUEL)/RN%NU(RN%I_OXIDIZER))
               ENDIF  
               IF (X_FU_S<X_FU_MIN) THEN
                  X_FU_S = X_FU_MIN
                  X_O_S = MAX(0.0_EB,X_O(NR)-(X_FU(NR)-X_FU_S)*RN%NU(RN%I_OXIDIZER)/RN%NU(RN%I_FUEL))
               ENDIF
               DX_FDT= -ETRM*X_FU_S**RN%N_F*X_O_S**RN%N_O
               X_FU_N(NR) = .5_EB*(X_FU(NR)+X_FU_S+DTT*DX_FDT)
               X_O_N(NR)  = .5_EB*(X_O(NR)+X_O_S+DTT*DX_FDT*RN%NU(RN%I_OXIDIZER)/RN%NU(RN%I_FUEL))
!               IF(I>=8 .AND. J==1 .AND. K==7) WRITE(*,*) II,NR,NO_REACTION,X_FU_N(NR),'C2'               
               IF (X_O_N(NR)<X_O_MIN) &
                  X_FU_N(NR) = MAX(0.0_EB,0.5_EB*((X_FU(NR)+X_FU_S)-(X_O(NR)+X_O_S)*RN%NU(RN%I_FUEL)/RN%NU(RN%I_OXIDIZER)))
               IF (X_FU_N(NR)<X_FU_MIN) X_FU_N(NR) = X_FU_MIN
!               IF(I>=8 .AND. J==1 .AND. K==7.AND. TMPD>1200.) WRITE(*,*) II,NR,NO_REACTION,X_FU_N(NR),X_FU(NR),'C'
               DYF(NR) = MAX(-X_FU(NR),X_FU_N(NR)-X_FU(NR))*RN%MW_FUEL*0.001_EB*1.E6_EB
               Q_NR(NR) = -DYF(NR) * HFAC_F(NR)
!               IF(I>=8 .AND. J==1 .AND. K==7.AND. TMPD>1200.) WRITE(*,*) I,Q_NR(NR),DYF(NR),'D'             
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
! Update HRR and species mass fractions
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

!FORALL (I=0:IBP1,J=0:JBP1,K=0:KBP1) TMP(I,J,K) = PBAR(K,PRESSURE_ZONE(I,J,K))/(RSUM(I,J,K)*RHO(I,J,K))
!TMP = MAX(TMPMIN,MIN(TMPMAX,TMP))

RETURN

END SUBROUTINE COMBUSTION_FR

END SUBROUTINE COMBUSTION
 
END MODULE FIRE


