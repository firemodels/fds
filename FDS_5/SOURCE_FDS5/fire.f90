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

IF( EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()

CALL POINT_TO_MESH(NM)

IF (MIXTURE_FRACTION) THEN
   CALL COMBUSTION_MF
ELSE
   CALL COMBUSTION_FR
ENDIF

TUSED(10,NM)=TUSED(10,NM)+SECOND()-TNOW

CONTAINS

SUBROUTINE COMBUSTION_MF
USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_MOLECULAR_WEIGHT
REAL(EB) :: YFU0,A,ETRM,YCOMIN,YO2MIN,YFUMIN,YO2W,YO20,YCO0,Y_SUM_W,&
            DYF,DX_FDT,HFAC_F,DTT,DELTA2,& 
            Y_O2_MAX,TMP_MIN,Y_O2_CORR,Q_NEW,Q_OLD,&
            F_TO_CO,DELTAH_CO,DYCO,HFAC_CO,Z_2,RHOX, &
            X_FU,X_O,X_FU_0,X_O_0,X_FU_S,X_O_S,X_FU_N,X_O_N,X_O_MIN,X_FU_MIN,COTOO2, &
            Y_F_MAX,TMP_F_MIN,Y_F_CORR,Z2_MIN
INTEGER :: NODETS,N,I,J,K,II,IC,IW,IWA(-3:3)
!LOGICAL :: BURN
REAL(EB), POINTER, DIMENSION(:,:,:) :: YO2Z,QT,R_SUM_DILUENTS,MIX_TIME
!LOGICAL, POINTER, DIMENSION(:,:,:) :: IGNITE

PRODUCE_CO: IF (.NOT. CO_PRODUCTION) THEN !Combustion without CO formation and destruction
   YO2Z => WORK1
   YO2Z = 0._EB
   MIX_TIME => WORK2
   MIX_TIME = DT
   Q =  0._EB
   Z_2 = 0._EB
   DO K=1,KBAR
      DO J=1,JBAR
         ILOOPA: DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOPA
            CALL GET_MASS_FRACTION(YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),O2_INDEX,Y_SUM(I,J,K),YO2Z(I,J,K))   
         ENDDO ILOOPA
      ENDDO
   ENDDO
   YO2MIN = 1.E-10_EB
   YFUMIN = 1.E-10_EB
   !Loop and do F -> CO
   RN => REACTION(1)
   HFAC_F   = RN%HEAT_OF_COMBUSTION/DT
   SUPPRESSIONIF: IF (SUPPRESSION) THEN !TMP and Y_O2 dependent combustion
      DO K=1,KBAR
         DO J=1,JBAR
            ILOOPB: DO I=1,IBAR
               IC = CELL_INDEX(I,J,K)
               IF (SOLID(IC)) CYCLE ILOOPB
               IWA = WALL_INDEX(IC,:)
               YO20  = YO2Z(I,J,K)
               YFU0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))*RN%Y_F_INLET
               IF (YFU0<=YFUMIN .OR. YO20<=YO2MIN) CYCLE ILOOPB

               !Get min O2
               Y_O2_MAX = YO20
               Y_F_MAX = YFU0/RN%Y_F_INLET
               TMP_MIN = TMP(I,J,K)
               TMP_F_MIN = TMP(I,J,K)
               !Check neighboring cells for fuel and oxygen
               !X direction
               IF (IWA(-1)==0) THEN
                  IF (YO2Z(I-1,J,K)>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2Z(I-1,J,K)
                     TMP_MIN = TMP(I-1,J,K)
                  ENDIF
                  IF (YY(I-1,J,K,I_FUEL)>Y_F_MAX) THEN
                     Y_F_MAX = YY(I-1,J,K,I_FUEL)
                     TMP_F_MIN = TMP(I-1,J,K)
                  ENDIF
               ELSE
                  IW = IWA(-1)
                  IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN
                     Y_SUM_W = 0._EB           
                     DO N=1,N_SPECIES
                        IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                           Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                        ENDIF
                     ENDDO
                     CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),Z_2,YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                     IF (YO2W>Y_O2_MAX) THEN
                        Y_O2_MAX = YO2W
                        TMP_MIN = TMP_F(IW)
                     ENDIF
                     IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                        Y_F_MAX = YY_W(IW,I_FUEL)
                        TMP_F_MIN = TMP_F(IW)
                     ENDIF
                  ENDIF
               ENDIF
               IF (IWA(1)==0) THEN
                  IF (YO2Z(I+1,J,K)>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2Z(I+1,J,K)
                     TMP_MIN = TMP(I+1,J,K)
                  ENDIF
                  IF (YY(I+1,J,K,I_FUEL)>Y_F_MAX) THEN
                     Y_F_MAX = YY(I+1,J,K,I_FUEL)
                     TMP_F_MIN = TMP(I+1,J,K)
                  ENDIF
               ELSE
                  IW = IWA(1)
                  IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN    
                     Y_SUM_W = 0._EB           
                     DO N=1,N_SPECIES
                        IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                           Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                        ENDIF
                     ENDDO
                     CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),Z_2,YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                     IF (YO2W>Y_O2_MAX) THEN
                        Y_O2_MAX = YO2W
                        TMP_MIN = TMP_F(IW)
                     ENDIF
                     IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                        Y_F_MAX = YY_W(IW,I_FUEL)
                        TMP_F_MIN = TMP_F(IW)
                     ENDIF
                  ENDIF
               ENDIF
               !Y direction            
               IF (IWA(-2)==0) THEN
                  IF (YO2Z(I,J-1,K)>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2Z(I,J-1,K)
                     TMP_MIN = TMP(I,J-1,K)
                  ENDIF
                  IF (YY(I,J-1,K,I_FUEL)>Y_F_MAX) THEN
                     Y_F_MAX = YY(I,J-1,K,I_FUEL)
                     TMP_F_MIN = TMP(I,J-1,K)
                  ENDIF
               ELSE
                  IW = IWA(-2)
                  IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN               
                     Y_SUM_W = 0._EB           
                     DO N=1,N_SPECIES
                        IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                           Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                        ENDIF
                     ENDDO
                     CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),Z_2,YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                     IF (YO2W>Y_O2_MAX) THEN
                        Y_O2_MAX = YO2W
                        TMP_MIN = TMP_F(IW)
                     ENDIF
                     IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                        Y_F_MAX = YY_W(IW,I_FUEL)
                        TMP_F_MIN = TMP_F(IW)
                     ENDIF
                  ENDIF
               ENDIF
               IF (IWA(2)==0) THEN
                  IF (YO2Z(I,J+1,K)>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2Z(I,J+1,K)
                     TMP_MIN = TMP(I,J+1,K)
                  ENDIF
                  IF (YY(I,J+1,K,I_FUEL)>Y_F_MAX) THEN
                     Y_F_MAX = YY(I,J+1,K,I_FUEL)
                     TMP_F_MIN = TMP(I,J+1,K)
                  ENDIF
               ELSE
                  IW = IWA(2)
                  IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN               
                     Y_SUM_W = 0._EB           
                     DO N=1,N_SPECIES
                        IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                           Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                        ENDIF
                     ENDDO
                     CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),Z_2,YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                     IF (YO2W>Y_O2_MAX) THEN
                        Y_O2_MAX = YO2W
                        TMP_MIN = TMP_F(IW)
                     ENDIF
                     IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                        Y_F_MAX = YY_W(IW,I_FUEL)
                        TMP_F_MIN = TMP_F(IW)
                     ENDIF
                  ENDIF
               ENDIF                        
               !Z direction
               IF (IWA(-3)==0) THEN
                  IF (YO2Z(I,J,K-1)>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2Z(I,J,K-1)
                     TMP_MIN = TMP(I,J,K-1)
                  ENDIF
                  IF (YY(I,J,K-1,I_FUEL)>Y_F_MAX) THEN
                     Y_F_MAX = YY(I,J,K-1,I_FUEL)
                     TMP_F_MIN = TMP(I,J,K-1)
                  ENDIF
               ELSE
                  IW = IWA(-3)
                  IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN               
                     Y_SUM_W = 0._EB           
                     DO N=1,N_SPECIES
                        IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                           Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                        ENDIF
                     ENDDO
                     CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),Z_2,YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                     IF (YO2W>Y_O2_MAX) THEN
                        Y_O2_MAX = YO2W
                        TMP_MIN = TMP_F(IW)
                     ENDIF
                     IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                        Y_F_MAX = YY_W(IW,I_FUEL)
                        TMP_F_MIN = TMP_F(IW)
                     ENDIF
                  ENDIF
               ENDIF
               IF (IWA(3)==0) THEN
                  IF (YO2Z(I,J,K+1)>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2Z(I,J,K+1)
                     TMP_MIN = TMP(I,J,K+1)
                  ENDIF
                  IF (YY(I,J,K+1,I_FUEL)>Y_F_MAX) THEN
                     Y_F_MAX = YY(I,J,K+1,I_FUEL)
                     TMP_F_MIN = TMP(I,J,K+1)
                  ENDIF
               ELSE
                  IW = IWA(3)
                  IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN
                     Y_SUM_W = 0._EB           
                     DO N=1,N_SPECIES
                        IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                           Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                        ENDIF
                     ENDDO
                     CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),Z_2,YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                     IF (YO2W>Y_O2_MAX) THEN
                        Y_O2_MAX = YO2W
                        TMP_MIN = TMP_F(IW)
                     ENDIF
                     IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                        Y_F_MAX = YY_W(IW,I_FUEL)
                        TMP_F_MIN = TMP_F(IW)
                     ENDIF
                  ENDIF
               ENDIF
               Y_O2_CORR = RN%Y_O2_LL*(RN%CRIT_FLAME_TMP-TMP_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
               Y_F_CORR  = RN%Y_F_LFL*(RN%CRIT_FLAME_TMP-TMP_F_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
               IF (Y_O2_MAX < Y_O2_CORR .OR. Y_F_MAX*RN%Y_F_INLET < Y_F_CORR) CYCLE ILOOPB
               DYF = MIN(YFU0,YO20/RN%O2_F_RATIO)
               IF (LES .AND. EDDY_BREAKUP) THEN
                  IF (.NOT.TWO_D) THEN
                     DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
                  ELSE
                     DELTA2 = DX(I)*DZ(K)
                  ENDIF
                  MIX_TIME(I,J,K) = SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
               ENDIF
               Q_NEW = MIN(Q_UPPER,DYF*RHO(I,J,K)*HFAC_F*MIN(1._EB,DT/MIX_TIME(I,J,K)))
               DYF = Q_NEW /(RHO(I,J,K)*HFAC_F*RN%Y_F_INLET)
               Q(I,J,K) = Q_NEW
               YY(I,J,K,I_FUEL)   = YY(I,J,K,I_FUEL)   - DYF
               YY(I,J,K,I_PROG_F) = YY(I,J,K,I_PROG_F) + DYF
            ENDDO ILOOPB
         ENDDO
      ENDDO
      
   ELSE SUPPRESSIONIF  ! No suppression

      DO K=1,KBAR
         DO J=1,JBAR
            ILOOPC: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ILOOPC
               YO20  = YO2Z(I,J,K)
               YFU0  = MAX(0._EB,MIN(1._EB,YY(I,J,K,I_FUEL)))*RN%Y_F_INLET
               IF (YFU0<=YFUMIN .OR. YO20<=YO2MIN) CYCLE ILOOPC
               DYF = MIN(YFU0,YO20/RN%O2_F_RATIO)
               IF (LES .AND. EDDY_BREAKUP) THEN
                  IF (.NOT.TWO_D) THEN
                     DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
                  ELSE
                     DELTA2 = DX(I)*DZ(K)
                  ENDIF
                  MIX_TIME(I,J,K) = SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
               ENDIF
               Q_NEW = MIN(Q_UPPER,DYF*RHO(I,J,K)*HFAC_F*MIN(1._EB,DT/MIX_TIME(I,J,K)))
               DYF = Q_NEW /(RHO(I,J,K)*HFAC_F*RN%Y_F_INLET)
               Q(I,J,K) = Q_NEW
               YY(I,J,K,I_FUEL)   = YY(I,J,K,I_FUEL)   - DYF
               YY(I,J,K,I_PROG_F) = YY(I,J,K,I_PROG_F) + DYF
            ENDDO ILOOPC
         ENDDO
      ENDDO
   ENDIF SUPPRESSIONIF
       
ELSE PRODUCE_CO  ! Combustion with suppression and CO production

   YO2Z   => WORK1
   YO2Z   = 0._EB
   MIX_TIME => WORK2
   MIX_TIME = 0._EB
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
            CALL GET_MASS_FRACTION(YY(I,J,K,I_FUEL),YY(I,J,K,I_PROG_CO),YY(I,J,K,I_PROG_F),O2_INDEX,Y_SUM(I,J,K),YO2Z(I,J,K))
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
   DO K=1,KBAR
      DO J=1,JBAR
         ILOOPY: DO I=1,IBAR
            IC = CELL_INDEX(I,J,K)
            IF (SOLID(IC)) CYCLE ILOOPY
            IWA = WALL_INDEX(IC,:)
            YO20  = YO2Z(I,J,K)
            YFU0  = MIN(1._EB,MAX(0._EB,YY(I,J,K,I_FUEL)))*RN%Y_F_INLET
            IF (YFU0<=YFUMIN .OR. YO20<=YO2MIN) CYCLE ILOOPY
            !Get min O2
            Y_O2_MAX = YO20
            Y_F_MAX = YFU0/RN%Y_F_INLET
            TMP_MIN = TMP(I,J,K)
            TMP_F_MIN = TMP(I,J,K)
            !Check neighboring cells for fuel and oxygen
            !X direction
            IF (IWA(-1)==0) THEN
               IF (YO2Z(I-1,J,K)>Y_O2_MAX) THEN
                  Y_O2_MAX = YO2Z(I-1,J,K)
                  TMP_MIN = TMP(I-1,J,K)
               ENDIF
               IF (YY(I-1,J,K,I_FUEL)>Y_F_MAX) THEN
                  Y_F_MAX = YY(I-1,J,K,I_FUEL)
                  TMP_F_MIN = TMP(I-1,J,K)
               ENDIF
            ELSE
               IW = IWA(-1)
               IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN
                  Y_SUM_W = 0._EB           
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                        Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                     ENDIF
                  ENDDO
                  CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),YY_W(IW,I_PROG_CO),YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                  IF (YO2W>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2W
                     TMP_MIN = TMP_F(IW)
                  ENDIF
                  IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                     Y_F_MAX = YY_W(IW,I_FUEL)
                     TMP_F_MIN = TMP_F(IW)
                  ENDIF
               ENDIF
            ENDIF
            IF (IWA(1)==0) THEN
               IF (YO2Z(I+1,J,K)>Y_O2_MAX) THEN
                  Y_O2_MAX = YO2Z(I+1,J,K)
                  TMP_MIN = TMP(I+1,J,K)
               ENDIF
               IF (YY(I+1,J,K,I_FUEL)>Y_F_MAX) THEN
                  Y_F_MAX = YY(I+1,J,K,I_FUEL)
                  TMP_F_MIN = TMP(I+1,J,K)
               ENDIF
            ELSE
               IW = IWA(1)
               IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN               
                  Y_SUM_W = 0._EB           
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                        Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                     ENDIF
                  ENDDO
                  CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),YY_W(IW,I_PROG_CO),YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                  IF (YO2W>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2W
                     TMP_MIN = TMP_F(IW)
                  ENDIF
                  IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                     Y_F_MAX = YY_W(IW,I_FUEL)
                     TMP_F_MIN = TMP_F(IW)
                  ENDIF
               ENDIF
            ENDIF
            !Y direction            
            IF (IWA(-2)==0) THEN
               IF (YO2Z(I,J-1,K)>Y_O2_MAX) THEN
                  Y_O2_MAX = YO2Z(I,J-1,K)
                  TMP_MIN = TMP(I,J-1,K)
               ENDIF
               IF (YY(I,J-1,K,I_FUEL)>Y_F_MAX) THEN
                  Y_F_MAX = YY(I,J-1,K,I_FUEL)
                  TMP_F_MIN = TMP(I,J-1,K)
               ENDIF
            ELSE
               IW = IWA(-2)
               IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN               
                  Y_SUM_W = 0._EB           
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                        Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                     ENDIF
                  ENDDO
                  CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),YY_W(IW,I_PROG_CO),YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                  IF (YO2W>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2W
                     TMP_MIN = TMP_F(IW)
                  ENDIF
                  IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                     Y_F_MAX = YY_W(IW,I_FUEL)
                     TMP_F_MIN = TMP_F(IW)
                  ENDIF
               ENDIF
            ENDIF
            IF (IWA(2)==0) THEN
               IF (YO2Z(I,J+1,K)>Y_O2_MAX) THEN
                  Y_O2_MAX = YO2Z(I,J+1,K)
                  TMP_MIN = TMP(I,J+1,K)
               ENDIF
               IF (YY(I,J+1,K,I_FUEL)>Y_F_MAX) THEN
                  Y_F_MAX = YY(I,J+1,K,I_FUEL)
                  TMP_F_MIN = TMP(I,J+1,K)
               ENDIF
            ELSE
               IW = IWA(2)
               IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN               
                  Y_SUM_W = 0._EB           
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                        Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                     ENDIF
                  ENDDO
                  CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),YY_W(IW,I_PROG_CO),YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                  IF (YO2W>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2W
                     TMP_MIN = TMP_F(IW)
                  ENDIF
                  IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                     Y_F_MAX = YY_W(IW,I_FUEL)
                     TMP_F_MIN = TMP_F(IW)
                  ENDIF
               ENDIF
            ENDIF       
            !Z direction
            IF (IWA(-3)==0) THEN
               IF (YO2Z(I,J,K-1)>Y_O2_MAX) THEN
                  Y_O2_MAX = YO2Z(I,J,K-1)
                  TMP_MIN = TMP(I,J,K-1)
               ENDIF
               IF (YY(I,J,K-1,I_FUEL)>Y_F_MAX) THEN
                  Y_F_MAX = YY(I,J,K-1,I_FUEL)
                  TMP_F_MIN = TMP(I,J,K-1)
               ENDIF
            ELSE
               IW = IWA(-3)
               IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN               
                  Y_SUM_W = 0._EB           
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                        Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                     ENDIF
                  ENDDO
                  CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),YY_W(IW,I_PROG_CO),YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                  IF (YO2W>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2W
                     TMP_MIN = TMP_F(IW)
                  ENDIF
                  IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                     Y_F_MAX = YY_W(IW,I_FUEL)
                     TMP_F_MIN = TMP_F(IW)
                  ENDIF
               ENDIF
            ENDIF
            IF (IWA(3)==0) THEN
               IF (YO2Z(I,J,K+1)>Y_O2_MAX) THEN
                  Y_O2_MAX = YO2Z(I,J,K+1)
                  TMP_MIN = TMP(I,J,K+1)
               ENDIF
               IF (YY(I,J,K+1,I_FUEL)>Y_F_MAX) THEN
                  Y_F_MAX = YY(I,J,K+1,I_FUEL)
                  TMP_F_MIN = TMP(I,J,K+1)
               ENDIF
            ELSE
               IW = IWA(3)
               IF (SURFACE(IJKW(5,IW))%SPECIES_BC_INDEX/=NO_MASS_FLUX) THEN
                  Y_SUM_W = 0._EB           
                  DO N=1,N_SPECIES
                     IF (SPECIES(N)%MODE==GAS_SPECIES) THEN
                        Y_SUM_W = Y_SUM_W + YY_W(IW,N)
                     ENDIF
                  ENDDO
                  CALL GET_MASS_FRACTION(YY_W(IW,I_FUEL),YY_W(IW,I_PROG_CO),YY_W(IW,I_PROG_F),O2_INDEX,Y_SUM_W,YO2W)   
                  IF (YO2W>Y_O2_MAX) THEN
                     Y_O2_MAX = YO2W
                     TMP_MIN = TMP_F(IW)
                  ENDIF
                  IF (MAX(0._EB,MIN(1._EB,YY_W(IW,I_FUEL)))>Y_F_MAX) THEN
                     Y_F_MAX = YY_W(IW,I_FUEL)
                     TMP_F_MIN = TMP_F(IW)
                  ENDIF
               ENDIF
            ENDIF
            Y_O2_CORR = RN%Y_O2_LL*(RN%CRIT_FLAME_TMP-TMP_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
            Y_F_CORR  = RN%Y_F_LFL*(RN%CRIT_FLAME_TMP-TMP_F_MIN)/(RN%CRIT_FLAME_TMP-TMPA)
            IF (Y_O2_MAX < Y_O2_CORR .OR. Y_F_MAX*RN%Y_F_INLET < Y_F_CORR) CYCLE ILOOPY
            DYF = MIN(YFU0,YO20/RN%O2_F_RATIO)
            IF (LES .AND. EDDY_BREAKUP) THEN
               IF (.NOT.TWO_D) THEN
                  DELTA2 = (DX(I)*DY(J)*DZ(K))**TWTH
               ELSE
                  DELTA2 = DX(I)*DZ(K)
               ENDIF
               MIX_TIME(I,J,K) = SC*RHO(I,J,K)*DELTA2/MU(I,J,K)
            ENDIF
            Q_NEW = MIN(Q_UPPER,DYF*RHO(I,J,K)*HFAC_F*MIN(1._EB,DT/MIX_TIME(I,J,K)))
            DYF = Q_NEW /(RHO(I,J,K)*HFAC_F*RN%Y_F_INLET)
            Q(I,J,K) = Q_NEW
            YY(I,J,K,I_FUEL) = YY(I,J,K,I_FUEL) - DYF
            YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) + DYF
            YO2Z(I,J,K) = YO2Z(I,J,K) - DYF * RN%Y_F_INLET * RN%O2_F_RATIO
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
            YCO0  = MAX(0._EB,YY(I,J,K,I_PROG_CO))*RN%Y_F_INLET / F_TO_CO 
            IF (YCO0<=YCOMIN .OR. YO20<=YO2MIN) CYCLE ILOOPY1
            !Get max conversion allowed
            Z2_MIN = REACTION(2)%NU_CO / (REACTION(2)%NU_CO+REACTION(2)%NU_CO2) * (YY(I,J,K,I_PROG_CO)+YY(I,J,K,I_PROG_F))
            IF (YY(I,J,K,I_PROG_CO) < Z2_MIN) CYCLE ILOOPY1
            !Get max CO2
            IF((TMP(I,J,K) < 500._EB .AND. Q(I,J,K)==0._EB) .OR. Q(I,J,K)>=Q_UPPER) CYCLE ILOOPY1
            IF (Q(I,J,K)/=0._EB) THEN
               DYCO = MIN(YCO0,YO20/COTOO2)*MIN(1._EB,DT/MIX_TIME(I,J,K))
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
            DYCO = Q_NEW/(RHO(I,J,K)*HFAC_CO*RN%Y_F_INLET) * F_TO_CO
            IF (YY(I,J,K,I_PROG_CO) - DYCO < Z2_MIN) THEN
               Q_NEW = (YY(I,J,K,I_PROG_CO) - Z2_MIN) / DYCO * Q_NEW
               DYCO = YY(I,J,K,I_PROG_CO) - Z2_MIN
            ENDIF
            Q(I,J,K) = Q_OLD + Q_NEW
            YY(I,J,K,I_PROG_CO) = YY(I,J,K,I_PROG_CO) - DYCO
            YY(I,J,K,I_PROG_F) = YY(I,J,K,I_PROG_F) + DYCO
         ENDDO ILOOPY1
      ENDDO
   ENDDO
ENDIF PRODUCE_CO
IF (N_SPEC_DILUENTS > 0) THEN
   R_SUM_DILUENTS => WORK4
   R_SUM_DILUENTS =  0._EB
   DO N=1,N_SPECIES
      IF (SPECIES(N)%MODE==GAS_SPECIES) R_SUM_DILUENTS(:,:,:) = R_SUM_DILUENTS(:,:,:) + SPECIES(N)%RCON*YY(:,:,:,N)
   ENDDO
ENDIF
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF (CO_PRODUCTION) THEN
            Z_2 = YY(I,J,K,I_PROG_CO)
         ELSE
            Z_2 = 0._EB
         ENDIF
         CALL GET_MOLECULAR_WEIGHT(YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),RSUM(I,J,K))
         RSUM(I,J,K) = R0/RSUM(I,J,K)
      ENDDO
   ENDDO
ENDDO

IF (N_SPEC_DILUENTS > 0) RSUM = RSUM*(1._EB-Y_SUM) + R_SUM_DILUENTS

END SUBROUTINE COMBUSTION_MF

SUBROUTINE COMBUSTION_FR
REAL(EB) :: TMPD,ETRM,&
            DX_FDT,DTT,YYNEW,WFAC,QFAC,&
            X_O_MIN,X_FU_MIN,X_FU_S,X_O_S
REAL(EB), DIMENSION(1:N_REACTIONS,0:N_SPECIES) :: MPUE
REAL(EB), DIMENSION(1:N_REACTIONS) :: HFAC_F, A, ETRM1, X_FU,X_O,&
           DYF,Q_NR,X_FU_N,X_O_N
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

! Solve the simple ODE to deplete fuel and oxidizer due to reaction
               DX_FDT= -ETRM*X_FU(NR)**RN%N_F*X_O(NR)**RN%N_O
               X_FU_S = X_FU(NR) + DTT*DX_FDT
               X_O_S = X_O(NR) + DTT*DX_FDT*RN%NU(RN%I_OXIDIZER)/RN%NU(RN%I_FUEL)
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
               IF (X_O_N(NR)<X_O_MIN) &
                  X_FU_N(NR) = MAX(0.0_EB,0.5_EB*((X_FU(NR)+X_FU_S)-(X_O(NR)+X_O_S)*RN%NU(RN%I_FUEL)/RN%NU(RN%I_OXIDIZER)))
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

RETURN

END SUBROUTINE COMBUSTION_FR

END SUBROUTINE COMBUSTION

SUBROUTINE GET_REV_fire(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') firerev(INDEX(firerev,':')+1:LEN_TRIM(firerev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') firedate

END SUBROUTINE GET_REV_fire
 
END MODULE FIRE


