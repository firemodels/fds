MODULE RAD
 
! Radiation heat transfer
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_VARIABLES
USE RADCONS

IMPLICIT NONE
PRIVATE

CHARACTER(255), PARAMETER :: radiid='$Id$'
CHARACTER(255), PARAMETER :: radirev='$Revision$'
CHARACTER(255), PARAMETER :: radidate='$Date$'

PUBLIC INIT_RADIATION,COMPUTE_RADIATION,NSB,NRA,UIIDIM,NRT,RSA,NRP,RTMPMAX,RTMPMIN,GET_REV_radi
 
CONTAINS
 
 
SUBROUTINE INIT_RADIATION
USE MEMORY_FUNCTIONS, ONLY : CHKMEMERR
USE MIEV
USE RADCALV
USE DEVICE_VARIABLES, ONLY : DEVICE, GAS_CELL_RAD_FLUX, GAS_CELL_RAD_DEVC_INDEX, N_GAS_CELL_RAD_DEVC
REAL(EB) :: THETAUP,THETALOW,PHIUP,PHILOW,F_THETA,MW_RADCAL,PLANCK_C2,KSI,LT, &
            RCRHO,YY,BBF,AP0,AMEAN,MTOT,XLENG,YLENG,ZLENG
INTEGER  :: N,I,J,K,IZERO,NN,NI,II,JJ,IIM,JJM,IBND,NS,I_RADCAL
TYPE(SPECIES_TYPE), POINTER :: SS
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC
 
! Determine the number of polar angles (theta)
 
NRA = NUMBER_RADIATION_ANGLES
IF (CYLINDRICAL) THEN
   NRT = NINT(SQRT(REAL(NRA)))
ELSEIF (TWO_D) THEN
   NRT = 1
ELSE
   NRT = 2*NINT(0.5_EB*1.17*REAL(NRA)**(1._EB/2.26))
ENDIF      
 
ALLOCATE(NRP(1:NRT),STAT=IZERO)
CALL ChkMemErr('INIT','NRP',IZERO)
 
! Determine number of azimuthal angles (phi)
 
N = 0
DO I=1,NRT
   IF (CYLINDRICAL) THEN
      NRP(I) = NINT(REAL(NRA)/(REAL(NRT)))
      ELSEIF (TWO_D) THEN
      NRP(I) = 4*NINT(0.25_EB*REAL(NRA))
      ELSE
      THETALOW = PI*REAL(I-1)/REAL(NRT)
      THETAUP  = PI*REAL(I)/REAL(NRT)
      NRP(I) = NINT(0.5_EB*REAL(NRA)*(COS(THETALOW)-COS(THETAUP)))
      NRP(I) = MAX(4,NRP(I))
      NRP(I) = 4*NINT(0.25_EB*REAL(NRP(I)))
      ENDIF
   N = N + NRP(I)
ENDDO
NRA = N
NUMBER_RADIATION_ANGLES = NRA
 
! Set the opening angle of the cylindrical geometry equal to the azimuthal angle
 
IF (CYLINDRICAL) DPHI0 = PI/REAL(NRP(1))
 
NSB = NUMBER_SPECTRAL_BANDS
 
ALLOCATE(RSA(1:NRA),STAT=IZERO)
CALL ChkMemErr('INIT','RSA',IZERO)
ALLOCATE(DLX(1:NRA),STAT=IZERO)
CALL ChkMemErr('INIT','DLX',IZERO)
ALLOCATE(DLY(1:NRA),STAT=IZERO)
CALL ChkMemErr('INIT','DLY',IZERO)
ALLOCATE(DLZ(1:NRA),STAT=IZERO)
CALL ChkMemErr('INIT','DLZ',IZERO)
IF (CYLINDRICAL) THEN
   ALLOCATE(DLB(1:NRA),STAT=IZERO)
   CALL ChkMemErr('INIT','DLB',IZERO)
ENDIF
ALLOCATE(DLN(-3:3,1:NRA),STAT=IZERO)
CALL ChkMemErr('INIT','DLN',IZERO)
ALLOCATE(DLM(1:NRA,3),STAT=IZERO)
CALL ChkMemErr('INIT','DLM',IZERO)

! Setup gas phase rad flux devices 

IF (GAS_CELL_RAD_FLUX) THEN
   DO N=1, N_GAS_CELL_RAD_DEVC
      ALLOCATE(DEVICE(GAS_CELL_RAD_DEVC_INDEX(N))%ILW(NRA,NSB),STAT=IZERO)
      CALL ChkMemErr('INIT','DV%ILW',IZERO)
      DEVICE(GAS_CELL_RAD_DEVC_INDEX(N))%ILW = 0._EB
   ENDDO
ENDIF
 
! Determine mean direction normals and sweeping orders
 
N = 0
DO I=1,NRT
   DO J=1,NRP(I)
      N = N + 1
      THETALOW  = PI*REAL(I-1)/REAL(NRT)
      THETAUP   = PI*REAL(I)/REAL(NRT)
      F_THETA   = 0.5_EB*(THETAUP-THETALOW  - COS(THETAUP)*SIN(THETAUP) + COS(THETALOW)*SIN(THETALOW))
      IF (CYLINDRICAL) THEN
         PHILOW = PI*REAL(J-1)/REAL(NRP(I))
         PHIUP  = PI*REAL(J)/REAL(NRP(I))
      ELSEIF (TWO_D) THEN
         PHILOW = TWOPI*REAL(J-1)/REAL(NRP(I)) + PIO2
         PHIUP  = TWOPI*REAL(J)/REAL(NRP(I))   + PIO2
      ELSE
         PHILOW = TWOPI*REAL(J-1)/REAL(NRP(I))
         PHIUP  = TWOPI*REAL(J)/REAL(NRP(I))
      ENDIF
      RSA(N) = (PHIUP-PHILOW)*(COS(THETALOW)-COS(THETAUP))
      IF (CYLINDRICAL) THEN
         DLX(N) = 2._EB*SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW)) *F_THETA
         DLY(N) =  (-SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW))  +COS(DPHI0/2.)*(COS(PHILOW)-COS(PHIUP)))*F_THETA
         DLB(N) =  (-SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW))  -COS(DPHI0/2.)*(COS(PHILOW)-COS(PHIUP)))*F_THETA
         DLZ(N)    = 0.5_EB*(PHIUP-PHILOW)   * ((SIN(THETAUP))**2-(SIN(THETALOW))**2)
      ELSEIF (TWO_D) THEN
         DLX(N) = (SIN(PHIUP)-SIN(PHILOW))*F_THETA
         DLY(N) = 0._EB
         DLZ(N) = (COS(PHILOW)-COS(PHIUP))*F_THETA
      ELSE
         DLX(N) = (SIN(PHIUP)-SIN(PHILOW))*F_THETA
         DLY(N) = (COS(PHILOW)-COS(PHIUP))*F_THETA
         DLZ(N)    = 0.5_EB*(PHIUP-PHILOW)      * ((SIN(THETAUP))**2-(SIN(THETALOW))**2)
      ENDIF
   ENDDO
ENDDO
 
! Set (wall normal)*(angle vector) value
 
DO N = 1,NRA
   DLN(-1,N) = -DLX(N)
   DLN( 1,N) =  DLX(N)
   DLN(-2,N) = -DLY(N)
   DLN( 2,N) =  DLY(N)
   DLN(-3,N) = -DLZ(N)
   DLN( 3,N) =  DLZ(N)
ENDDO
 
! Calculate mirroring matrix
 
N = 0
DO I=1,NRT
   DO J=1,NRP(I)
      N = N + 1
      DO K=1,3
         IF (TWO_D .AND. .NOT.CYLINDRICAL) THEN
            SELECT CASE(K)
               CASE(1)             ! X-surfaces
                  IIM = 1
                  JJM = NRP(I) - J + 1
               CASE(2)             ! Y-surfaces
                  IIM = 1
                  JJM = J
               CASE(3)             ! Z-surfaces
                  IIM = 1
                  JJM = NRP(I)/2 - J + 1
            END SELECT
            JJM = MODULO(JJM,NRP(I))
            IF (JJM==0) JJM = NRP(I)
         ELSE
            SELECT CASE(K)
               CASE(1)             ! X-surfaces
                  IIM = I
                  JJM = NRP(I)/2 - J + 1
               CASE(2)             ! Y-surfaces
                  IIM = I
                  JJM = NRP(I) - J + 1
               CASE(3)             ! Z-surfaces
                  IIM = NRT - I + 1
                  JJM = J
            END SELECT
            IIM = MODULO(IIM,NRT)
            JJM = MODULO(JJM,NRP(I))
            IF (IIM==0) IIM = NRT
            IF (JJM==0) JJM = NRP(I)
         ENDIF
 
         NN = 0
         DO II = 1,IIM
            DO JJ = 1,NRP(II)
               NN = NN + 1
               IF ((II==IIM).AND.(JJ==JJM)) NI = NN
            ENDDO
         ENDDO
         DLM(N,K) = NI
      ENDDO
   ENDDO
ENDDO
 
!-----------------------------------------------------
!
!            Radiative properties computation
!
!-----------------------------------------------------

!
! General parameters
!
RTMPMAX = 2400._EB       ! Maximum temperature for property tables
RTMPMIN = 300._EB      ! Minimum temperature for property tables
! 
! Setup spectral information
!
INIT_WIDE_BAND: IF (WIDE_BAND_MODEL) THEN
 
! Fraction of blackbody emission in a wavelength interval
 
   PLANCK_C2 = 14387.69_EB            ! micron.K
   NLAMBDAT  = 4000
   LTSTEP    = 25.0_EB           ! maximum LAMBDA*T = NLANBDAT*LTSTEP
   ALLOCATE(BBFRAC(0:NLAMBDAT),STAT=IZERO)
   CALL ChkMemErr('INIT','BBFRAC',IZERO)
   BBFRAC = 0._EB
   LT = 0._EB
   DO I = 1,NLAMBDAT
      LT =  LT + LTSTEP
      KSI = PLANCK_C2/LT
      DO J = 1,50
         BBFRAC(I) = BBFRAC(I) + EXP(-KSI*REAL(J))/REAL(J) * (KSI**3 + 3.*KSI**2/REAL(J) + 6.*KSI/REAL(J)**2 + 6./REAL(J)**3)
      ENDDO
   ENDDO
   BBFRAC =  BBFRAC * 15._EB/PI**4
 
! Define band limit wave lengths in micrometers
 
   ALLOCATE(WL_LOW(1:NSB),STAT=IZERO)
   CALL ChkMemErr('INIT','WL_LOW',IZERO)
   ALLOCATE(WL_HIGH(1:NSB),STAT=IZERO)
   CALL ChkMemErr('INIT','WL_HIGH',IZERO)
   IF (CH4_BANDS) THEN
      WL_LOW(1:NSB) =(/1.00_EB, 2.63_EB, 2.94_EB, 3.57_EB, 4.17_EB, 4.6_EB, 7.00_EB, 8.62_EB, 10.0_EB /)
      WL_HIGH(1:NSB)=(/2.63_EB, 2.94_EB, 3.57_EB, 4.17_EB, 4.60_EB, 7.0_EB, 8.62_EB, 10.0_EB, 200._EB /) 
   ELSE
      WL_LOW(1:NSB) =(/1.00_EB, 2.63_EB, 2.94_EB, 4.17_EB, 4.6_EB, 10.0_EB /)
      WL_HIGH(1:NSB)=(/2.63_EB, 2.94_EB, 4.17_EB, 4.6_EB, 10.0_EB, 200.0_EB /)
   ENDIF
 
ENDIF INIT_WIDE_BAND
 
!----------------------------------------------------------------------------
!
!     Tables for gas phase absorption coefficient
!
!     CONTROLLING PROGRAM FOR SUBROUTINE "RADCAL", A NARROW-BAND
!     MODEL FOR CALCULATING SPECTRAL INTENSITY (W/M-2/SR/MICRON) AND
!     SPECTRAL TRANSMITTANCE VERSUS WAVELENGTH (MICRONS) IN A NONISO-
!     THERMAL, VARIABLE COMPOSITION  MIXTURE OF CO2, H2O, CO, N2, O2,
!     CH4, AND SOOT. FOR A HOMOGENEOUS PATH, THE PROGRAM ALSO COMPUTES
!     THE PLANCK-MEAN ABSORPTION COEF, AP0, THE INCIDENT-MEAN ABSORPTION
!     COEFFICIENT, AIWALL, AND THE EFFECTIVE-MEAN ABSORPTION COEFFICIENT,
!     AMEAN, ALL IN UNITS OF INVERSE METERS.
!
!     INPUT PARAMETERS:
!          NPT=NUMBER OF HOMOGENEOUS ELEMENTS
!          DD(J)=THICKNESS OF J TH ELEMENT, M
!          RCT(J)=TEMPERATURE OF J TH ELEMENT, K.
!          P(I,J)=PARTIAL PRESSURE OF GASEOUS COMPONENTS, kPa:
!                  I   GASEOUS SPECIES
!                  1        CO2
!                  2        H2O
!                  3        CH4
!                  4        CO
!                  5        O2
!                  6        N2
!          SVF(J)=SOOT VOLUME FRACTION OF J TH ELEMENT
!          OMMIN=MINIMUM WAVE NUMBER IN SPECTRUM, CM-1.
!          OMMAX=MAXIMUM WAVE NUMBER IN SPECTRUM, CM-1.
!
!-------------------------------------------------------------------------
 
CALL RCALLOC  ! Allocate arrays for RadCal
 
! 20% of mean beam length, Eq 8-51, Holman, 7th Ed. Heat Transfer. Length = 3.6*Volume/Area
 
XLENG = MESHES(1)%XF-MESHES(1)%XS
YLENG = MESHES(1)%YF-MESHES(1)%YS
ZLENG = MESHES(1)%ZF-MESHES(1)%ZS
IF (PATH_LENGTH < 0.0_EB) THEN  ! default was -1.0
   IF (TWO_D) THEN ! calculate based on the geometry
      PATH_LENGTH = MIN( 10._EB , 0.1_EB*3.6_EB*XLENG*ZLENG/(XLENG+ZLENG) )
   ELSE
      PATH_LENGTH = MIN( 10._EB , 0.1_EB*3.6_EB*XLENG*YLENG*ZLENG/(XLENG*YLENG+XLENG*ZLENG+YLENG*ZLENG) )
   ENDIF
ENDIF
DD(1) = MAX(PATH_LENGTH,1.0E-4_EB)
 
! Using RadCal, create look-up tables for the absorption coefficients for all gas species, mixture fraction or aerosols

SPECIES_LOOP: DO NS=1,N_SPECIES

   SS => SPECIES(NS)
   IF (.NOT. SS%ABSORBING) CYCLE SPECIES_LOOP
 
   GAS_TYPE: SELECT CASE (SS%MODE) 

      CASE (GAS_SPECIES) GAS_TYPE

         SS%NKAP_TEMP = 21
         SS%NKAP_MASS = 200
         SS%MAXMASS=1._EB
         ALLOCATE(SS%KAPPA(0:SS%NKAP_MASS,0:SS%NKAP_TEMP,1:NSB),STAT=IZERO)
         CALL ChkMemErr('RADI','KAPPA',IZERO)
         IF (NS==I_CO2) THEN
            I_RADCAL = 1
            MW_RADCAL = MW_CO2
         ELSEIF (NS==I_CO) THEN
            I_RADCAL = 4
            MW_RADCAL = MW_CO
         ELSEIF (NS==I_WATER) THEN
            I_RADCAL = 2
            MW_RADCAL = MW_H2O
         ELSE
            IF (SS%ABSORBING) THEN
               I_RADCAL = 3
               MW_RADCAL = 16
            ELSE
               I_RADCAL = 5
               MW_RADCAL = 32
            ENDIF
         END IF
         BAND_LOOP_FR1: DO IBND = 1,NSB
            IF (NSB>1) THEN
               OMMIN = REAL(NINT(1.E4_EB/WL_HIGH(IBND)),EB)
               OMMAX = REAL(NINT(1.E4_EB/WL_LOW(IBND)),EB)
            ELSE
               OMMIN = 50._EB
               OMMAX = 10000._EB
            ENDIF
            CALL INIT_RADCAL
            T_LOOP_MF2: DO K = 0,SS%NKAP_TEMP
               RCT(1) = RTMPMIN + K*(RTMPMAX-RTMPMIN)/SS%NKAP_TEMP
               Z_LOOP_MF2: DO J=0,SS%NKAP_MASS
                  YY = SS%MAXMASS*(REAL(J)/REAL(SS%NKAP_MASS))**4
                  MTOT = YY/MW_RADCAL+(1._EB-YY)/MW_N2
                  SPECIE = 0._EB
                  SPECIE(I_RADCAL) = YY
                  P(:,1) = 0._EB
                  SVF = 0._EB
                  P(I_RADCAL,1) = YY/MW_RADCAL/MTOT
                  P(6,1) = (1-YY)/MW_N2/MTOT
                  CALL RADCAL(AMEAN,AP0)
                  IF (NSB==1 .AND. PATH_LENGTH > 0._EB) THEN
                     SS%KAPPA(J,K,IBND) = MIN(AMEAN,AP0)
                  ELSE
                     IF (NSB==1) THEN
                        BBF = 1._EB
                     ELSE
                        BBF = BLACKBODY_FRACTION(WL_LOW(1),WL_HIGH(1),RCT(1))
                     ENDIF
                     SS%KAPPA(J,K,IBND) = AP0/BBF
                  ENDIF
              ENDDO Z_LOOP_MF2
            ENDDO T_LOOP_MF2
         ENDDO BAND_LOOP_FR1

      CASE (AEROSOL_SPECIES) GAS_TYPE

         SS%NKAP_TEMP = 21
         SS%NKAP_MASS = 200  
         SS%MAXMASS=0.2_EB   
         ALLOCATE(SS%KAPPA(0:SS%NKAP_MASS,0:SS%NKAP_TEMP,1:NSB),STAT=IZERO)
         CALL ChkMemErr('RADI','KAPPA',IZERO)
         BAND_LOOP_FR2: DO IBND = 1,NSB
            IF (NSB>1) THEN
               OMMIN = REAL(NINT(1.E4_EB/WL_HIGH(IBND)),EB)
               OMMAX = REAL(NINT(1.E4_EB/WL_LOW(IBND)),EB)
            ELSE
               OMMIN = 50._EB
               OMMAX = 10000._EB
            ENDIF
            CALL INIT_RADCAL
            I_RADCAL=5
            T_LOOP_MF3: DO K = 0,SS%NKAP_TEMP
               RCT(1) = RTMPMIN + K*(RTMPMAX-RTMPMIN)/SS%NKAP_TEMP
               Z_LOOP_MF3: DO J=0,SS%NKAP_MASS
                  YY = SS%MAXMASS*(REAL(J)/REAL(SS%NKAP_MASS))**4
                  RCRHO = 29._EB * P_INF/(R0*RCT(1)) !Good enough
                  SPECIE = 0._EB               
                  P(:,1) = 0._EB               
                  P(6,1) = 1._EB
                  SPECIE(I_RADCAL) = YY*RCRHO/RHO_SOOT                  
                  SVF(1) = YY*RCRHO/RHO_SOOT                                    
                  CALL RADCAL(AMEAN,AP0)
                  IF (NSB==1 .AND. PATH_LENGTH > 0._EB) THEN
                     SS%KAPPA(J,K,IBND) = MIN(AMEAN,AP0)
                  ELSE
                     IF (NSB==1) THEN
                        BBF = 1._EB
                     ELSE
                        BBF = BLACKBODY_FRACTION(WL_LOW(1),WL_HIGH(1),RCT(1))
                     ENDIF
                     SS%KAPPA(J,K,IBND) = AP0/BBF
                  ENDIF
               ENDDO Z_LOOP_MF3
            ENDDO T_LOOP_MF3
         ENDDO BAND_LOOP_FR2

   END SELECT GAS_TYPE

ENDDO SPECIES_LOOP

IF (MIXTURE_FRACTION) THEN
   Z2KAPPA_T = 42
   Z2KAPPA_M = 50
   ALLOCATE (Z2KAPPA_M4(5),STAT=IZERO)
   CALL ChkMemErr('RADI','Z2KAPPA_M4',IZERO)
   Z2KAPPA_M4(1) = REAL(Z2KAPPA_M,EB)**4/REACTION(1)%MW_FUEL
   Z2KAPPA_M4(2) = REAL(Z2KAPPA_M,EB)**4/MW_CO2
   Z2KAPPA_M4(3) = REAL(Z2KAPPA_M,EB)**4/MW_CO
   Z2KAPPA_M4(4) = REAL(Z2KAPPA_M,EB)**4/MW_H2O
   Z2KAPPA_M4(5) = REAL(Z2KAPPA_M,EB)**4*5._EB
   ALLOCATE (Z2KAPPA(5,0:Z2KAPPA_M,0:Z2KAPPA_T,NSB),STAT=IZERO)
   CALL ChkMemErr('RADI','Z2KAPPA',IZERO)
   Z2KAPPA = 0._EB
   BBF = 1._EB
   OMMIN = 50._EB
   OMMAX = 10000._EB
   BAND_LOOP_Z: DO IBND = 1,NSB
      IF (NSB>1) THEN
       OMMIN = REAL(NINT(1.E4_EB/WL_HIGH(IBND)),EB)
       OMMAX = REAL(NINT(1.E4_EB/WL_LOW(IBND)),EB)
      ENDIF
      CALL INIT_RADCAL 
      T_LOOP_Z: DO K = 0,Z2KAPPA_T
         RCT(1) = RTMPMIN + K*(RTMPMAX-RTMPMIN)/Z2KAPPA_T         
         IF (NSB>1) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),RCT(1))
         Y_LOOP_Z: DO J=0,Z2KAPPA_M
            YY = (REAL(J,EB)/REAL(Z2KAPPA_M,EB))**4
            SVF(1) = 0._EB
            !FUEL
            SPECIE = 0._EB
            P = 0._EB
            SPECIE(3) = YY
            P(3,1) = YY
            P(6,1) = (1._EB-YY)
            CALL RADCAL(AMEAN,AP0)
            IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
               Z2KAPPA(1,J,K,IBND) = MIN(AMEAN,AP0)
            ELSE
               Z2KAPPA(1,J,K,IBND) = AP0/BBF
            ENDIF
            !CO2
            SPECIE = 0._EB
            P = 0._EB
            SPECIE(1) = YY
            P(1,1) = YY
            P(6,1) = (1._EB-YY)
            CALL RADCAL(AMEAN,AP0)
            IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
               Z2KAPPA(2,J,K,IBND) = MIN(AMEAN,AP0)
            ELSE
               Z2KAPPA(2,J,K,IBND) = AP0/BBF
            ENDIF
            !CO
            SPECIE = 0._EB
            P = 0._EB
            SPECIE(4) = YY
            P(4,1) = YY
            P(6,1) = (1._EB-YY)
            CALL RADCAL(AMEAN,AP0)
            IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
               Z2KAPPA(3,J,K,IBND) = MIN(AMEAN,AP0)
            ELSE
               Z2KAPPA(3,J,K,IBND) = AP0/BBF
            ENDIF
            !H2O
            SPECIE = 0._EB
            P = 0._EB
            SPECIE(2) = YY
            P(2,1) = YY
            P(6,1) = (1._EB-YY)
            CALL RADCAL(AMEAN,AP0)
            IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
               Z2KAPPA(4,J,K,IBND) = MIN(AMEAN,AP0)
            ELSE
               Z2KAPPA(4,J,K,IBND) = AP0/BBF
            ENDIF
            !Soot
            RCRHO = MW_AIR*P_INF/(R0*RCT(1))
            YY = YY * 0.2_EB
            SPECIE = 0._EB
            P = 0._EB
            SPECIE(5) = YY*RCRHO/RHO_SOOT
            P(6,1) = 1._EB
            SVF(1) = YY*RCRHO/RHO_SOOT
            CALL RADCAL(AMEAN,AP0)
            IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
               Z2KAPPA(5,J,K,IBND) = MIN(AMEAN,AP0)
            ELSE
               Z2KAPPA(5,J,K,IBND) = AP0/BBF
            ENDIF
         ENDDO Y_LOOP_Z
      ENDDO T_LOOP_Z
   ENDDO BAND_LOOP_Z
ENDIF

CALL RCDEALLOC  ! Deallocate RadCal arrays
 
 
!-----------------------------------------------------
!
!     Tables for droplet absorption coefficients
! 
!-----------------------------------------------------

DROPLETS: IF (N_EVAP_INDICIES>0) THEN
   GET_PC_RADI: DO J=1,N_PART
      PC => PARTICLE_CLASS(J)
      IF ((.NOT. PC%FUEL) .AND. (.NOT. PC%WATER)) CYCLE GET_PC_RADI
      IF (PC%FUEL) THEN  ! Tuntomo, Tien and Park, Comb. Sci. and Tech., 84, 133-140, 1992
         PC%WQABS(:,1) = (/0,10,16,28,52,98,191,386,792,1629,3272,6163, &
            10389,15588,20807,23011,22123,22342,22200,22241,21856, &
            22795,23633,24427,25285,26207,27006,27728,28364,28866, &
            29260/)!,29552,29748,30000/)
         PC%KWR(0) = 0._EB
         DO I=1,NDG
            PC%KWR(I)=EXP(I/2.5_EB-4._EB)
         ENDDO
         PC%WQABS=PC%WQABS/10000._EB
         PC%KWR=PC%KWR/1000000._EB
      ENDIF
      IF (PC%WATER) CALL MEAN_CROSS_SECTIONS(J)
   ENDDO GET_PC_RADI
ENDIF DROPLETS
 
! A few miscellaneous constants

FOUR_SIGMA = 4._EB*SIGMA
RPI_SIGMA  = RPI*SIGMA
 
! In axially symmetric case, each angle represents two symmetric angles. So weight the intensities by two.

W_AXI = 1._EB
IF (CYLINDRICAL) W_AXI = 2._EB
 
END SUBROUTINE INIT_RADIATION


 
SUBROUTINE COMPUTE_RADIATION(NM)

! Call radiation routine or simply specify the radiative loss

USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY : SECOND  
REAL(EB) :: TNOW
INTEGER, INTENT(IN) :: NM 

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()
 
CALL POINT_TO_MESH(NM)

IF (RADIATION) THEN
   CALL RADIATION_FVM(NM)
ELSE
   IF (N_REACTIONS>0) QR = -RADIATIVE_FRACTION*Q
ENDIF

TUSED(9,NM)=TUSED(9,NM)+SECOND()-TNOW

CONTAINS 
 
SUBROUTINE RADIATION_FVM(NM)
USE MIEV
USE MATH_FUNCTIONS, ONLY : INTERPOLATE1D
USE DEVICE_VARIABLES, ONLY : DEVICE_TYPE,DEVICE, GAS_CELL_RAD_FLUX, GAS_CELL_RAD_DEVC_INDEX, N_GAS_CELL_RAD_DEVC
REAL(EB) :: ZZ, RAP, AX, AXU, AXD, AY, AYU, AYD, AZ, VC, RU, RD, RP, &
            ILXU, ILYU, ILZU, QVAL, BBF, BBFA, NCSDROP, RSA_RAT, WAXIDLN, KAPPA_1, Z_2, COSINE, &
            Q_SUM,K_SUM,U_SUM,KAPPA_CORRECTOR,VOL
INTEGER  :: N, NN,IIG,JJG,KKG,I,J,K,IW,II,JJ,KK,IOR,IC,IWUP,IWDOWN, &
            ISTART, IEND, ISTEP, JSTART, JEND, JSTEP, &
            KSTART, KEND, KSTEP, NSTART, NEND, NSTEP, &
            I_UIID, N_UPDATES, IBND, TYY, NOM, NS, IBC,EVAP_INDEX
LOGICAL :: UPDATE_INTENSITY
REAL(EB), POINTER, DIMENSION(:,:,:) :: KFST4, IL, UIIOLD, KAPPAW, KFST4W, EXTCOE, SCAEFF
REAL(EB), POINTER, DIMENSION(:)     :: OUTRAD_W, INRAD_W
INTEGER, INTENT(IN) :: NM
TYPE (OMESH_TYPE), POINTER :: M2
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(DEVICE_TYPE), POINTER :: DV
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC

KFST4    => WORK1
IL       => WORK2
UIIOLD   => WORK3
EXTCOE   => WORK4
KAPPAW   => WORK5
SCAEFF   => WORK6
KFST4W   => WORK7
OUTRAD_W => WALL_WORK1
INRAD_W  => WALL_WORK2
 

! Ratio of solid angle, used in scattering
 
RSA_RAT = 1._EB/(1._EB-1._EB/NRA)
 
! Check if it time to update radiation intensity field
 
RAD_CALL_COUNTER  = RAD_CALL_COUNTER+1
IF ( MOD(RAD_CALL_COUNTER,TIME_STEP_INCREMENT)==1 .OR. TIME_STEP_INCREMENT==1) THEN
   UPDATE_INTENSITY = .TRUE.
ELSE
   UPDATE_INTENSITY = .FALSE.
ENDIF

IF (WIDE_BAND_MODEL) QR = 0._EB
IF (UPDATE_INTENSITY) QRADIN = 0._EB
 
! Loop over spectral bands
 
BAND_LOOP: DO IBND = 1,NUMBER_SPECTRAL_BANDS
 
   KAPPAW = 0._EB
   KFST4W = 0._EB
   SCAEFF = 0._EB
 
   ! Calculate fraction on ambient black body radiation
 
   IF (NUMBER_SPECTRAL_BANDS==1) THEN
      BBFA = 1._EB
   ELSE
      BBFA = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMPA)
   ENDIF      

   ! Generate water absorption and scattering coefficients
 
   IF_DROPLETS_INCLUDED: IF (NLP>0 .AND. N_EVAP_INDICIES>0) THEN
      IF (NUMBER_SPECTRAL_BANDS==1) THEN
         BBF = 1._EB
      ELSE
         BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),RADTMP)
      ENDIF      
      DO K=1,KBAR
         DO J=1,JBAR
            ZLOOPM: DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ZLOOPM
                  PC_LOOP: DO N = 1,N_PART
                     PC => PARTICLE_CLASS(N)
                     EVAP_INDEX = PC%EVAP_INDEX
                     IF (EVAP_INDEX==0) CYCLE PC_LOOP
                     IF (AVG_DROP_DEN(I,J,K,EVAP_INDEX)==0._EB) CYCLE PC_LOOP
                     NCSDROP = THFO*AVG_DROP_DEN(I,J,K,EVAP_INDEX)/ (PC%DENSITY*AVG_DROP_RAD(I,J,K,EVAP_INDEX))
                     ! Absorption and scattering efficiency
                     CALL INTERPOLATE1D(PC%KWR,PC%WQABS(:,IBND), 0.95_EB*AVG_DROP_RAD(I,J,K,EVAP_INDEX),QVAL) 
                     KAPPAW(I,J,K) = KAPPAW(I,J,K) + NCSDROP*QVAL
                     IF (PC%WATER) THEN
                        CALL INTERPOLATE1D(PC%KWR,PC%WQSCA(:,IBND),0.95_EB*AVG_DROP_RAD(I,J,K,EVAP_INDEX),QVAL)
                        SCAEFF(I,J,K) = SCAEFF(I,J,K) + NCSDROP*QVAL
                     ENDIF
                     KFST4W(I,J,K) = KFST4W(I,J,K)+ BBF*KAPPAW(I,J,K)*FOUR_SIGMA*AVG_DROP_TMP(I,J,K,EVAP_INDEX)**4
                  END DO PC_LOOP
            ENDDO ZLOOPM
         ENDDO
      ENDDO
      QR_W = 0._EB
   ENDIF IF_DROPLETS_INCLUDED
 
   ! Compute absorption coefficient KAPPA and source term KAPPA*4*SIGMA*TMP**4
 
   BBF   = 1._EB
   KAPPA = KAPPA0
   U_SUM = 0._EB
   K_SUM = 0._EB
   Q_SUM = 0._EB

   DO K=1,KBAR
      DO J=1,JBAR
         ZLOOP: DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE ZLOOP
            SUM_SPECIES: DO NS=1,N_SPECIES
               IF (SPECIES(NS)%MODE == MIXTURE_FRACTION_SPECIES) CYCLE SUM_SPECIES
               IF (.NOT. SPECIES(NS)%ABSORBING)                  CYCLE SUM_SPECIES
               TYY = NINT(SPECIES(NS)%NKAP_TEMP*(TMP(I,J,K)-RTMPMIN)/(RTMPMAX-RTMPMIN))
               TYY = MAX(0,MIN(SPECIES(NS)%NKAP_TEMP,TYY))
               ZZ = MIN(1._EB,MAX(0._EB,YY(I,J,K,NS)))
               KAPPA(I,J,K) = KAPPA(I,J,K) + YY2KAPPA(ZZ,TYY,IBND,NS)
            ENDDO SUM_SPECIES
            IF (MIXTURE_FRACTION) THEN
               TYY = NINT(Z2KAPPA_T*(TMP(I,J,K)-RTMPMIN)/(RTMPMAX-RTMPMIN))
               TYY = MAX(0,MIN(Z2KAPPA_T,TYY))
               IF (CO_PRODUCTION) THEN
                  Z_2 = YY(I,J,K,I_PROG_CO)
               ELSE
                  Z_2 = 0._EB
               ENDIF
               CALL GET_KAPPA(YY(I,J,K,I_FUEL),Z_2,YY(I,J,K,I_PROG_F),Y_SUM(I,J,K),KAPPA_1,TYY,IBND)      
               KAPPA(I,J,K) = KAPPA(I,J,K) + KAPPA_1
            ENDIF
            IF (WIDE_BAND_MODEL)  BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMP(I,J,K))
            KFST4(I,J,K) = BBF*KAPPA(I,J,K)*FOUR_SIGMA*TMP(I,J,K)**4
            IF (RADIATIVE_FRACTION*Q(I,J,K)>0._EB) THEN
               KFST4(I,J,K) = MAX(KFST4(I,J,K),BBF*RADIATIVE_FRACTION*Q(I,J,K))
               VOL = R(I)*DX(I)*DY(J)*DZ(K)
               IF (WIDE_BAND_MODEL) THEN
                  U_SUM = U_SUM + KAPPA(I,J,K)*UIID(I,J,K,IBND)*VOL
               ELSE
                  U_SUM = U_SUM + KAPPA(I,J,K)*UII(I,J,K)*VOL
               ENDIF
               K_SUM = K_SUM + KFST4(I,J,K)*VOL
               Q_SUM = Q_SUM + BBF*RADIATIVE_FRACTION*Q(I,J,K)*VOL
            ENDIF
         ENDDO ZLOOP
      ENDDO
   ENDDO

   ! Add a corrective factor to the radiative source term to achieve desired RADIATIVE_FRACTION

   IF (RADIATIVE_FRACTION>0._EB .AND. K_SUM>0._EB) THEN
      KAPPA_CORRECTOR = (Q_SUM+U_SUM)/K_SUM
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (Q(I,J,K)>0._EB) KFST4(I,J,K) = KFST4(I,J,K)*KAPPA_CORRECTOR
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Calculate extinction coefficient
 
   EXTCOE = KAPPA + KAPPAW + SCAEFF*RSA_RAT

   ! Update intensity field
 
   INTENSITY_UPDATE: IF (UPDATE_INTENSITY) THEN
 
      IF (WIDE_BAND_MODEL) THEN
         UIIOLD = UIID(:,:,:,IBND)
      ELSE
         UIIOLD = UII
      ENDIF
      UII = 0._EB

      ! Compute boundary condition intensity emissivity*sigma*Tw**4/pi or emissivity*QRADOUT/pi for wall with internal radiation
 
      BBF = 1.0_EB
      DO IW = 1,NWC
         IF (WIDE_BAND_MODEL) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMP_F(IW))
         IBC = IJKW(5,IW)
         SF  => SURFACE(IBC)
         IF (.NOT. SF%INTERNAL_RADIATION) THEN
            QRADOUT(IW) = E_WALL(IW)*SIGMA*TMP_F(IW)**4
         ENDIF
         OUTRAD_W(IW) = BBF*RPI*QRADOUT(IW)
      ENDDO

      ! Compute boundary condition term incoming radiation integral
 
      DO IW = 1,NWC
         IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) CYCLE
         IOR = IJKW(4,IW)
         INRAD_W(IW) = SUM(-W_AXI*DLN(IOR,:)* WALL(IW)%ILW(:,IBND),1, DLN(IOR,:)<0._EB)
      ENDDO
 
      ! If updating intensities first time, sweep ALL angles
 
      N_UPDATES = 1
      IF (RAD_CALL_COUNTER==1) N_UPDATES = ANGLE_INCREMENT

      UPDATE_LOOP: DO I_UIID = 1,N_UPDATES
 
      ! Update counters inside the radiation routine
 
         ANGLE_INC_COUNTER = MOD(ANGLE_INC_COUNTER,ANGLE_INCREMENT) + 1
         IF (WIDE_BAND_MODEL) THEN
            UIID(:,:,:,IBND) = 0._EB
         ELSE
            UIID(:,:,:,ANGLE_INC_COUNTER) = 0._EB
         ENDIF
 
         ! Set the bounds and increment for the angleloop. Step downdard because in cylindrical case the Nth angle 
         ! boundary condition comes from (N+1)th angle.
 
         NSTART    = NRA - ANGLE_INC_COUNTER + 1
         NEND      = 1
         NSTEP     = -ANGLE_INCREMENT
         IL(:,:,:) = BBFA*RPI_SIGMA*TMPA4

         ANGLE_LOOP: DO N = NSTART,NEND,NSTEP  ! Sweep through control angles
 
            ! Boundary conditions: Intensities leaving the boundaries.
 
            WALL_LOOP1: DO IW=1,NWC
               IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY .OR. BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP1
               IOR = IJKW(4,IW)
               IF (DLN(IOR,N) < 0._EB) CYCLE WALL_LOOP1
               II  = IJKW(1,IW)
               JJ  = IJKW(2,IW)
               KK  = IJKW(3,IW)
               IF (.NOT.TWO_D .OR. ABS(IOR)/=2) THEN
                  SELECT CASE (BOUNDARY_TYPE(IW))
                     CASE (OPEN_BOUNDARY) 
                           IL(II,JJ,KK) = BBFA*RPI_SIGMA*TMPA4
                     CASE (MIRROR_BOUNDARY) 
                        WALL(IW)%ILW(N,IBND) = WALL(IW)%ILW(DLM(N,ABS(IOR)),IBND)
                        IL(II,JJ,KK) = WALL(IW)%ILW(N,IBND)
                     CASE (INTERPOLATED_BOUNDARY) 
                        IL(II,JJ,KK) = WALL(IW)%ILW(N,IBND)
                     CASE DEFAULT ! solid wall
                        WALL(IW)%ILW(N,IBND) = OUTRAD_W(IW) + RPI*(1._EB-E_WALL(IW))* INRAD_W(IW)
                  END SELECT
               ELSEIF (CYLINDRICAL) THEN
                  IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) CYCLE WALL_LOOP1
                  IL(II,JJ,KK) = WALL(IW)%ILW(N,IBND)
               ENDIF
            ENDDO WALL_LOOP1
 
            ! Determine sweep direction in physical space
 
            ISTART = 1
            JSTART = 1
            KSTART = 1
            IEND   = IBAR
            JEND   = JBAR
            KEND   = KBAR
            ISTEP  = 1
            JSTEP  = 1
            KSTEP  = 1
            IF (DLX(N) < 0._EB) THEN
               ISTART = IBAR
               IEND   = 1
               ISTEP  = -1
            ENDIF
            IF (DLY(N) < 0._EB) THEN
               JSTART = JBAR
               JEND   = 1
               JSTEP  = -1
            ENDIF
            IF (DLZ(N) < 0._EB) THEN
               KSTART = KBAR
               KEND   = 1
               KSTEP  = -1
            ENDIF
 
            GEOMETRY: IF (CYLINDRICAL) THEN  ! Sweep in axisymmetric geometry

               J = 1
               CKLOOP: DO K=KSTART,KEND,KSTEP
                  CILOOP: DO I=ISTART,IEND,ISTEP
                     IC = CELL_INDEX(I,J,K)
                     IF (SOLID(IC)) CYCLE CILOOP
                     ILXU = IL(I-ISTEP,J,K)
                     ILYU = IL(I,J-JSTEP,K)
                     ILZU = IL(I,J,K-KSTEP)
                     IF (DLX(N)>=0._EB) THEN
                        RU  = R(I-1)
                        RD  = R(I)
                     ELSE
                        RU  = R(I)
                        RD  = R(I-1)
                     ENDIF
                     RP  = SQRT(0.5_EB*(RU**2+RD**2))
                     VC  = DX(I)  * RP*DPHI0 * DZ(K)
                     AXU =          RU       * DZ(K) * ABS(DLX(N))
                     AXD =          RD       * DZ(K) * ABS(DLX(N))
                     AYU = DX(I)             * DZ(K) * ABS(DLB(N))
                     AYD = DX(I)             * DZ(K) * ABS(DLY(N))
                     AZ  = DX(I)  * RP*DPHI0         * ABS(DLZ(N))
                     IF (MODULO(N,NRP(1))==1) AYD = 0._EB  ! Zero out the terms involving symmetric overhang
                     IF (MODULO(N,NRP(1))==0) AYU = 0._EB
                     IF (IC/=0) THEN
                        IW = WALL_INDEX(IC,-ISTEP)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILXU = WALL(IW)%ILW(N,IBND)
                           AYU = 0.5*AYU
                           AZ = 0.5*AZ
                        ENDIF
                        IW = WALL_INDEX(IC,-JSTEP*2)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILYU = WALL(IW)%ILW(N,IBND)
                           AXU = 0.5*AXU
                           AZ = 0.5*AZ
                        ENDIF
                        IW = WALL_INDEX(IC,-KSTEP*3)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILZU = WALL(IW)%ILW(N,IBND)
                           AXU = 0.5*AXU
                           AYU = 0.5*AYU
                        ENDIF
                     ENDIF
                     RAP = 1._EB/(AXD+AYD+AZ+EXTCOE(I,J,K)*VC*RSA(N))
                     IL(I,J,K) = MAX(0._EB, RAP * (AXU*ILXU + AYU*ILYU + AZ*ILZU +  &
                                 VC*RSA(N)*RFPI*( KFST4(I,J,K)+KFST4W(I,J,K) +RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
                  ENDDO CILOOP
               ENDDO CKLOOP

            ELSEIF (TWO_D) THEN GEOMETRY  ! Sweep in 2D cartesian geometry

               J = 1
               K2LOOP: DO K=KSTART,KEND,KSTEP
                  I2LOOP: DO I=ISTART,IEND,ISTEP
                     IC = CELL_INDEX(I,J,K)
                     IF (SOLID(IC)) CYCLE I2LOOP
                     ILXU  = IL(I-ISTEP,J,K)
                     ILZU  = IL(I,J,K-KSTEP)
                     VC  = DX(I) * DZ(K)
                     AX  =         DZ(K) * ABS(DLX(N))
                     AZ  = DX(I)         * ABS(DLZ(N))
                     IF (IC/=0) THEN
                        IW = WALL_INDEX(IC,-ISTEP)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILXU = WALL(IW)%ILW(N,IBND)
                           AZ = 0.5*AZ
                        ENDIF
                        IW = WALL_INDEX(IC,-KSTEP*3)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILZU = WALL(IW)%ILW(N,IBND)
                           AX = 0.5*AX
                        ENDIF
                     ENDIF
                     RAP = 1._EB/(AX+AZ+EXTCOE(I,J,K)*VC*RSA(N))
                     IL(I,J,K) = MAX(0._EB, RAP * (AX*ILXU + AZ*ILZU + &
                                 VC*RSA(N)*RFPI*(KFST4(I,J,K)+KFST4W(I,J,K) +  RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) ) 
                  ENDDO I2LOOP
               ENDDO K2LOOP

            ELSE GEOMETRY  ! Sweep in 3D cartesian geometry

               KLOOP: DO K=KSTART,KEND,KSTEP
                  JLOOP: DO J=JSTART,JEND,JSTEP
                     ILOOP: DO I=ISTART,IEND,ISTEP
                        IC = CELL_INDEX(I,J,K)
                        IF (SOLID(IC)) CYCLE ILOOP
                        ILXU  = IL(I-ISTEP,J,K)
                        ILYU  = IL(I,J-JSTEP,K)
                        ILZU  = IL(I,J,K-KSTEP)
                        VC  = DX(I) * DY(J) * DZ(K)
                        AX  =         DY(J) * DZ(K) * ABS(DLX(N))
                        AY  = DX(I)         * DZ(K) * ABS(DLY(N))
                        AZ  = DX(I) * DY(J)         * ABS(DLZ(N))
                        IF (IC/=0) THEN
                           IW = WALL_INDEX(IC,-ISTEP)
                           IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                              AY = 0.5*AY
                              AZ = 0.5*AZ
                              ILXU = WALL(IW)%ILW(N,IBND)
                           ENDIF
                           IW = WALL_INDEX(IC,-JSTEP*2)
                           IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                              AX = 0.5*AX
                              AZ = 0.5*AZ
                              ILYU = WALL(IW)%ILW(N,IBND)
                           ENDIF
                           IW = WALL_INDEX(IC,-KSTEP*3)
                           IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                              AX = 0.5*AX
                              AY = 0.5*AY
                              ILZU = WALL(IW)%ILW(N,IBND)
                           ENDIF
                        ENDIF
                        RAP = 1._EB/(AX+AY+AZ+EXTCOE(I,J,K)*VC*RSA(N))
                        IL(I,J,K) = MAX(0._EB, RAP * ( AX*ILXU + AY*ILYU + AZ*ILZU + &
                                    VC*RSA(N)*RFPI*( KFST4(I,J,K)+KFST4W(I,J,K) + RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
                     ENDDO ILOOP
                  ENDDO JLOOP
               ENDDO KLOOP
 
            ENDIF GEOMETRY
 
            ! Boundary values: Incoming radiation
 
            WALL_LOOP2: DO IW=1,NWC
               IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY)   CYCLE WALL_LOOP2     
               IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY)   CYCLE WALL_LOOP2  
               IF (BOUNDARY_TYPE(IW)==POROUS_BOUNDARY) CYCLE WALL_LOOP2  
               IOR = IJKW(4,IW)
               IF (TWO_D .AND. .NOT.CYLINDRICAL  .AND. ABS(IOR)==2) CYCLE WALL_LOOP2  ! 2-D non cylindrical
               IF (DLN(IOR,N)>=0._EB) CYCLE WALL_LOOP2     ! outgoing
               WAXIDLN = - W_AXI*DLN(IOR,N)
               IIG = IJKW(6,IW)
               JJG = IJKW(7,IW)
               KKG = IJKW(8,IW)
               INRAD_W(IW) = INRAD_W(IW) - WAXIDLN * WALL(IW)%ILW(N,IBND) ! update incoming radiation,step 1
               WALL(IW)%ILW(N,IBND) = IL(IIG,JJG,KKG)
               INRAD_W(IW) = INRAD_W(IW) + WAXIDLN * WALL(IW)%ILW(N,IBND) ! update incoming radiation,step 2
            ENDDO WALL_LOOP2
 
            ! Copy the Y-downwind intensities to Y-upwind in cylindrical case
 
            IF (CYLINDRICAL) THEN
               J=1
               DO K=1,KBAR
                  DO I=1,IBAR
                     IWUP   = WALL_INDEX(CELL_INDEX(I,J,K),-2)
                     IWDOWN = WALL_INDEX(CELL_INDEX(I,J,K), 2)
                     WALL(IWUP)%ILW(MAX(1,N-1),IBND) = WALL(IWDOWN)%ILW(N,IBND)
                  ENDDO
               ENDDO
            ENDIF
 
            ! Calculate integrated intensity UIID
 
            IF (WIDE_BAND_MODEL) THEN
               UIID(:,:,:,IBND) = UIID(:,:,:,IBND) + W_AXI*RSA(N)*IL
            ELSE
               UIID(:,:,:,ANGLE_INC_COUNTER) = UIID(:,:,:,ANGLE_INC_COUNTER) + W_AXI*RSA(N)*IL
            ENDIF
 
            ! Interpolate boundary intensities onto other meshes
 
            INTERPOLATION_LOOP: DO NOM=1,NMESHES
               IF (NM==NOM) CYCLE INTERPOLATION_LOOP
               IF (EVACUATION_ONLY(NOM)) CYCLE INTERPOLATION_LOOP
               IF (NIC(NOM,NM)==0) CYCLE INTERPOLATION_LOOP
               M2=>OMESH(NOM)
               OTHER_WALL_LOOP: DO IW=1,MESHES(NOM)%NEWC
                  IF (M2%IJKW(9,IW)/=NM .OR. M2%BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE OTHER_WALL_LOOP
                  IOR = M2%IJKW(4,IW)
                  IF (DLN(IOR,N)<=0._EB) CYCLE OTHER_WALL_LOOP
                  M2%WALL(IW)%ILW(N,IBND)=IL(M2%IJKW(10,IW),M2%IJKW(11,IW),M2%IJKW(12,IW))
               ENDDO OTHER_WALL_LOOP
            ENDDO INTERPOLATION_LOOP
            IF (GAS_CELL_RAD_FLUX) THEN
               DEVC_LOOP: DO NN = 1, N_GAS_CELL_RAD_DEVC
                  DV => DEVICE(GAS_CELL_RAD_DEVC_INDEX(NN))
                  IF (DV%MESH /= NM) CYCLE DEVC_LOOP
                  COSINE = DV%ORIENTATION(1)*DLX(N)+DV%ORIENTATION(2)*DLY(N)+DV%ORIENTATION(3)*DLZ(N)
                  IF (COSINE <0._EB) DV%ILW(N,IBND) = -COSINE * IL(DV%I,DV%J,DV%K)
               ENDDO DEVC_LOOP
            ENDIF
         ENDDO ANGLE_LOOP
      ENDDO UPDATE_LOOP
 
      ! Compute incoming flux on walls 

      DO IW=1,NWC
         IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) CYCLE 
         IBC = IJKW(5,IW)
         QRADIN(IW)  = QRADIN(IW) + E_WALL(IW)*(INRAD_W(IW)+SURFACE(IBC)%EXTERNAL_FLUX/NUMBER_SPECTRAL_BANDS)
      ENDDO 
   ENDIF INTENSITY_UPDATE
 
   ! Save source term for the energy equation (QR = -DIV Q)

   IF (WIDE_BAND_MODEL) THEN
      QR = QR + KAPPA*UIID(:,:,:,IBND)-KFST4
      IF (NLP>0 .AND. N_EVAP_INDICIES>0) QR_W = QR_W + KAPPAW*UIID(:,:,:,IBND) - KFST4W
   ENDIF

ENDDO BAND_LOOP

! Sum up the parts of the intensity

IF (UPDATE_INTENSITY) UII = SUM(UIID, DIM = 4)

! Save source term for the energy equation (QR = -DIV Q). Done only in one-band (gray gas) case.

IF (.NOT. WIDE_BAND_MODEL) THEN
   QR  = KAPPA*UII - KFST4
   IF (NLP>0 .AND. N_EVAP_INDICIES>0) QR_W = QR_W + KAPPAW*UII - KFST4W
ENDIF

! Check for VIRTUAL devices

IF (GAS_CELL_RAD_FLUX) THEN
   DEVC_LOOP2: DO NN = 1, N_GAS_CELL_RAD_DEVC
      DV => DEVICE(GAS_CELL_RAD_DEVC_INDEX(NN))
      IF (DV%MESH /= NM) CYCLE DEVC_LOOP2
      IW = DV%VIRTUAL_WALL_INDEX
      IF (IW>0) THEN
         QRADIN(IW) = E_WALL(IW)*SUM(DV%ILW)
      ENDIF
   ENDDO DEVC_LOOP2
ENDIF

END SUBROUTINE RADIATION_FVM

END SUBROUTINE COMPUTE_RADIATION



REAL(EB) FUNCTION YY2KAPPA(YY_IN,TYY,IBND,NS)

! Calculate KAPPA as a function of YY

REAL(EB), INTENT(IN) :: YY_IN
REAL(EB) :: MAX_MASS,YY_K,INT_FAC
INTEGER, INTENT(IN) :: TYY,IBND,NS
INTEGER  :: LBND,UBND,NIYY
REAL(EB), POINTER, DIMENSION(:,:,:) :: KAPPA

KAPPA => SPECIES(NS)%KAPPA
NIYY = SPECIES(NS)%NKAP_MASS
MAX_MASS = SPECIES(NS)%MAXMASS

YY_K = MAX(0._EB,MIN(YY_IN,MAX_MASS))
INT_FAC = REAL(NIYY,EB)*(YY_K/MAX_MASS)**0.25_EB

LBND = INT(INT_FAC)
INT_FAC = INT_FAC - LBND   
LBND = MIN(LBND,NIYY)
UBND = MIN(LBND+1,NIYY)

YY2KAPPA = KAPPA(LBND,TYY,IBND)
YY2KAPPA = YY2KAPPA + INT_FAC*(KAPPA(UBND,TYY,IBND)-YY2KAPPA)

END FUNCTION YY2KAPPA

 
REAL(EB) FUNCTION BLACKBODY_FRACTION(L1,L2,TEMP)
 
! Calculates the fraction of black body radiation between wavelengths L1 and L2 (micron) in Temperature TEMP
 
REAL(EB) :: L1,L2,TEMP,LT1,LT2,BBFLOW,BBFHIGH
INTEGER  :: IYY

LT1    =   L1 * TEMP/LTSTEP
LT2    =   L2 * TEMP/LTSTEP
IYY = MIN(NLAMBDAT,MAX(0,NINT(LT1)))
BBFLOW = BBFRAC(IYY)
IYY = MIN(NLAMBDAT,MAX(0,NINT(LT2)))
BBFHIGH = BBFRAC(IYY)
BLACKBODY_FRACTION = BBFHIGH - BBFLOW

END FUNCTION BLACKBODY_FRACTION

SUBROUTINE GET_KAPPA(Z1,Z2,Z3,YY_SUM,KAPPA,TYY,IBND)
USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION_ALL,GET_MOLECULAR_WEIGHT
! GET_KAPPA returns the radiative absorption
REAL(EB) :: YY_SUM,KAPPA,Z1,Z2,Z3,Y_MF(9),MWA
REAL(EB) :: K_FUEL,K_CO2,K_CO,K_H2O,K_SOOT,INT_FAC
INTEGER :: IBND,TYY,LBND,UBND

IF (YY_SUM >=1._EB) THEN
   KAPPA = 0._EB
   RETURN
ELSE
   CALL GET_MASS_FRACTION_ALL(Z1,Z2,Z3,YY_SUM,Y_MF)
   CALL GET_MOLECULAR_WEIGHT(Z1,Z2,Z3,YY_SUM,MWA)
   INT_FAC = MAX(0._EB,(Y_MF(FUEL_INDEX)*MWA*Z2KAPPA_M4(1)))**0.25_EB
   LBND = INT(INT_FAC)
   INT_FAC = INT_FAC - LBND   
   LBND = MIN(LBND,Z2KAPPA_M)
   UBND = MIN(LBND+1,Z2KAPPA_M)
   K_FUEL = Z2KAPPA (1,LBND,TYY,IBND)
   K_FUEL = K_FUEL+INT_FAC*(Z2KAPPA (1,UBND,TYY,IBND)-K_FUEL)

   INT_FAC = MAX(0._EB,(Y_MF(CO2_INDEX)*MWA*Z2KAPPA_M4(2)))**0.25_EB
   LBND = INT(INT_FAC)
   INT_FAC = INT_FAC - LBND  
   LBND = MIN(LBND,Z2KAPPA_M)
   UBND = MIN(LBND+1,Z2KAPPA_M)
   K_CO2 = Z2KAPPA (2,LBND,TYY,IBND)
   K_CO2 = K_CO2 + INT_FAC*(Z2KAPPA (2,UBND,TYY,IBND)-K_CO2)

   INT_FAC = MAX(0._EB,(Y_MF(CO_INDEX)*MWA*Z2KAPPA_M4(3)))**0.25_EB
   LBND = INT(INT_FAC)
   INT_FAC = INT_FAC - LBND   
   LBND = MIN(LBND,Z2KAPPA_M)
   UBND = MIN(LBND+1,Z2KAPPA_M)
   K_CO = Z2KAPPA (3,LBND,TYY,IBND)
   K_CO = K_CO + INT_FAC*(Z2KAPPA (3,UBND,TYY,IBND)-K_CO)
   
   INT_FAC = MAX(0._EB,(Y_MF(H2O_INDEX)*MWA*Z2KAPPA_M4(4)))**0.25_EB
   LBND = INT(INT_FAC)
   INT_FAC = INT_FAC - LBND  
   LBND = MIN(LBND,Z2KAPPA_M)
   UBND = MIN(LBND+1,Z2KAPPA_M)
   K_H2O = Z2KAPPA (4,LBND,TYY,IBND)
   K_H2O = K_H2O + INT_FAC*(Z2KAPPA (4,UBND,TYY,IBND)-K_H2O)

   INT_FAC = MAX(0._EB,(Y_MF(SOOT_INDEX)*Z2KAPPA_M4(5)))**0.25_EB
   LBND = INT(INT_FAC)
   INT_FAC = INT_FAC - LBND   
   LBND = MIN(LBND,Z2KAPPA_M)
   UBND = MIN(LBND+1,Z2KAPPA_M)
   K_SOOT = Z2KAPPA (5,LBND,TYY,IBND)
   K_SOOT = K_SOOT + INT_FAC*(Z2KAPPA (5,UBND,TYY,IBND)-K_SOOT)
   KAPPA = K_FUEL + K_CO2 + K_CO + K_H2O + K_SOOT
ENDIF

END SUBROUTINE GET_KAPPA 

SUBROUTINE GET_REV_radi(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') radirev(INDEX(radirev,':')+1:LEN_TRIM(radirev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') radidate

END SUBROUTINE GET_REV_radi
 
END MODULE RAD

