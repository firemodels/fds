MODULE RAD

! Radiation heat transfer

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE MESH_VARIABLES
USE RADCONS

IMPLICIT NONE
PRIVATE

PUBLIC INIT_RADIATION,COMPUTE_RADIATION,BLACKBODY_FRACTION

CONTAINS


SUBROUTINE INIT_RADIATION

! Meanings of some variables defined here:
!
! THETAUP       Upper limit of spherical polar angle THETA on a solid angle interval
! THETALOW      Lower limit of spherical polar angle THETA on a solid angle interval
! PHIUP         Upper limit of spherical polar angle PHI on a solid angle interval
! PHILOW        Lower limit of spherical polar angle PHI on a solid angle interval
! PLANCK_C2     Second Planck radiation constant, 14387.69 micron.K (or 1.438769 cm.K)
! N, I, J, K    Integers used as indices of arrays in DO loops
! NRA           An integer that specifies the number of discrete radiation angles in the FVM
! NSB           Number of spectral bands in the piecewise constant absorption-emission spectrum
! IPC           An integer for looping over the Lagrangian particle classes


USE MEMORY_FUNCTIONS, ONLY : CHKMEMERR
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
USE MIEV
USE RADCAL_CALC
REAL(EB) :: THETAUP,THETALOW,PHIUP,PHILOW,F_THETA,PLANCK_C2,KSI,LT,RCRHO,YY,YY2,BBF,AP0,AMEAN,RADIANCE,TRANSMISSIVITY,X_N2,NRM
INTEGER  :: N,I,J,K,IPC,IZERO,NN,NI,II,JJ,IIM,JJM,IBND,NS,NS2,NRA,NSB,RADCAL_TEMP(16)=0,RCT_SKIP=-1,OR_IN,I1,I2,NM,NSTEPS
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE (RAD_FILE_TYPE), POINTER :: RF
REAL(EB), ALLOCATABLE, DIMENSION(:) :: COSINE_ARRAY,COSINE_ARRAY_2

! A few miscellaneous constants

FOUR_SIGMA = 4._EB*SIGMA            ! Four times the Stefan-Boltzmann constant
RPI_SIGMA  = RPI*SIGMA              ! Stefan-Boltzmann constant divided by PI (RPI = reciprocal pi)

NRA = NUMBER_RADIATION_ANGLES
NSB = NUMBER_SPECTRAL_BANDS

! Set the opening angle of the cylindrical geometry equal to the azimuthal angle

IF (CYLINDRICAL) DPHI0 = PI/REAL(NRP(1))

ALLOCATE(RSA(1:NRA),STAT=IZERO)
CALL ChkMemErr('RADI','RSA',IZERO)
ALLOCATE(DLX(1:NRA),STAT=IZERO)
CALL ChkMemErr('RADI','DLX',IZERO)
ALLOCATE(DLY(1:NRA),STAT=IZERO)
CALL ChkMemErr('RADI','DLY',IZERO)
ALLOCATE(DLZ(1:NRA),STAT=IZERO)
CALL ChkMemErr('RADI','DLZ',IZERO)
IF (CYLINDRICAL) THEN
   ALLOCATE(DLB(1:NRA),STAT=IZERO)
   CALL ChkMemErr('RADI','DLB',IZERO)
ENDIF
ALLOCATE(DLN(-3:3,1:NRA),STAT=IZERO)
CALL ChkMemErr('RADI','DLN',IZERO)
ALLOCATE(DLM(1:NRA,3),STAT=IZERO)
CALL ChkMemErr('RADI','DLM',IZERO)

! Determine mean direction normals and sweeping orders
! as described in the FDS Tech. Ref. Guide Vol. 1 Sec. 6.2.2.

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
         DLX(N) =  (SIN(PHIUP)-SIN(PHILOW)) *F_THETA
         DLY(N) =  (-SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW))  +COS(DPHI0/2.)*(COS(PHILOW)-COS(PHIUP)))*F_THETA
         DLB(N) =  (-SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW))  -COS(DPHI0/2.)*(COS(PHILOW)-COS(PHIUP)))*F_THETA
         DLZ(N)    = 0.5_EB*(PHIUP-PHILOW)   * ((SIN(THETAUP))**2-(SIN(THETALOW))**2)
         IF (N==1000000) WRITE(LU_ERR,'(A)') 'This line should never get executed. It is here only to prevent optimization.'
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
   DLN( 0,N) = 0._EB !prevent undefined variable errors
   DLN(-1,N) = -DLX(N)
   DLN( 1,N) =  DLX(N)
   DLN(-2,N) = -DLY(N)
   DLN( 2,N) =  DLY(N)
   DLN(-3,N) = -DLZ(N)
   DLN( 3,N) =  DLZ(N)
ENDDO

! In axially symmetric case, each angle represents two symmetric angles. So weight the intensities by two.

WEIGH_CYL = 1._EB
IF (CYLINDRICAL) THEN
   WEIGH_CYL = 2._EB
   ! Wall direction cosines are only used for flux integrations, so they can by multiplied in advance.
   DLN = WEIGH_CYL * DLN
ENDIF

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


! General parameters

RTMPMAX = 2470._EB     ! Maximum temperature for property tables
RTMPMIN = 270._EB      ! Minimum temperature for property tables

! Setup spectral information

INIT_WIDE_BAND: IF (WIDE_BAND_MODEL) THEN

   ! Fraction of blackbody emission in a wavelength interval

   PLANCK_C2 = 14387.69_EB       ! Value of the 2nd Planck radiation constant in micron.K
   NLAMBDAT  = 4000
   LTSTEP    = 25.0_EB           ! maximum LAMBDA*T = NLANBDAT*LTSTEP
   ALLOCATE(BBFRAC(0:NLAMBDAT),STAT=IZERO)
   CALL ChkMemErr('INIT','BBFRAC',IZERO)

   BBFRAC = 0._EB
   LT     = 0._EB

   DO I = 1,NLAMBDAT
      LT  = LT + LTSTEP
      KSI = PLANCK_C2/LT

      DO J = 1,50
         BBFRAC(I) = BBFRAC(I) + EXP(-KSI*REAL(J))/REAL(J) * (KSI**3 + 3.*KSI**2/REAL(J) + 6.*KSI/REAL(J)**2 + 6./REAL(J)**3)
      ENDDO
   ENDDO

   BBFRAC =  BBFRAC * 15._EB/PI**4

   ! Define band limit wave lengths in micrometers

   IF (.NOT.ALLOCATED(WL_LOW).OR. .NOT.ALLOCATED(WL_HIGH)) THEN

      ALLOCATE(WL_LOW(1:NSB),STAT=IZERO)
      CALL ChkMemErr('INIT','WL_LOW',IZERO)
      ALLOCATE(WL_HIGH(1:NSB),STAT=IZERO)
      CALL ChkMemErr('INIT','WL_HIGH',IZERO)

      ! Define the band limits as function of the fuel used
      ! Use (SPECIES(FUEL_INDEX)%RADCAL_ID (FUEL_INDEX global variable defined in module GLOBAL_CONSTANTS

      IF (FUEL_INDEX>0) THEN

         SELECT CASE(SPECIES(FUEL_INDEX)%RADCAL_ID)
            CASE('METHANE')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.630_EB, 2.940_EB, 4.170_EB, 4.600_EB,  10.000_EB /)
               WL_HIGH(1:NSB) = (/ 2.630_EB, 2.940_EB, 4.170_EB, 4.600_EB, 10.00_EB, 200.000_EB /)
            CASE('ETHANE')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.632_EB, 2.985_EB, 4.000_EB, 6.061_EB,   9.174_EB /)
               WL_HIGH(1:NSB) = (/ 2.632_EB, 2.985_EB, 4.000_EB, 6.061_EB, 9.174_EB, 200.000_EB /)
            CASE('ETHYLENE')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.632_EB, 2.963_EB, 3.571_EB, 6.061_EB,  12.821_EB /)
               WL_HIGH(1:NSB) = (/ 2.632_EB, 2.963_EB, 3.571_EB, 6.061_EB,12.821_EB, 200.000_EB /)
            CASE('PROPANE')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.632_EB, 2.985_EB, 3.922_EB, 6.061_EB,   8.511_EB /)
               WL_HIGH(1:NSB) = (/ 2.632_EB, 2.985_EB, 3.922_EB, 6.061_EB, 8.511_EB, 200.000_EB /)
            CASE('PROPYLENE')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.632_EB, 3.077_EB, 3.846_EB, 5.128_EB,   8.511_EB /)
               WL_HIGH(1:NSB) = (/ 2.632_EB, 3.077_EB, 3.846_EB, 5.128_EB, 8.511_EB, 200.000_EB /)
            CASE('N-HEPTANE')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.632_EB, 3.077_EB, 3.922_EB, 5.634_EB,   9.091_EB /)
               WL_HIGH(1:NSB) = (/ 2.632_EB, 3.077_EB, 3.922_EB, 5.634_EB, 9.091_EB, 200.000_EB /)
            CASE('TOLUENE')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.632_EB, 3.125_EB, 3.922_EB, 4.878_EB,   8.333_EB /)
               WL_HIGH(1:NSB) = (/ 2.632_EB, 3.125_EB, 3.922_EB, 4.878_EB, 8.333_EB, 200.000_EB /)
            CASE('METHANOL')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.614_EB, 3.125_EB, 3.846_EB, 5.970_EB,   8.889_EB /)
               WL_HIGH(1:NSB) = (/ 2.614_EB, 3.125_EB, 3.846_EB, 5.970_EB, 8.889_EB, 200.000_EB /)
            CASE('MMA')
               WL_LOW(1:NSB)  = (/ 1.000_EB, 2.632_EB, 3.077_EB, 3.774_EB, 4.878_EB,   8.000_EB /)
               WL_HIGH(1:NSB) = (/ 2.632_EB, 3.077_EB, 3.774_EB, 4.878_EB, 8.000_EB, 200.000_EB /)
         END SELECT

      ELSE

         ! Use the methane bands if there is no fuel species

         WL_LOW(1:NSB)  = (/ 1.000_EB, 2.630_EB, 2.940_EB, 4.170_EB, 4.600_EB,  10.000_EB /)
         WL_HIGH(1:NSB) = (/ 2.630_EB, 2.940_EB, 4.170_EB, 4.600_EB, 10.00_EB, 200.000_EB /)

      ENDIF

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
!          SEGMENT_LENGTH_M(J)=THICKNESS OF J TH ELEMENT, M
!          TEMP_GAS(J)=TEMPERATURE OF J TH ELEMENT, K.
!          PARTIAL_PRESSURES_ATM(I,J)=PARTIAL PRESSURE OF GASEOUS COMPONENTS I, kPa:
!          OMMIN=MINIMUM WAVE NUMBER IN SPECTRUM, CM-1.
!          OMMAX=MAXIMUM WAVE NUMBER IN SPECTRUM, CM-1.
!
!-------------------------------------------------------------------------

MAKE_KAPPA_ARRAYS: IF (.NOT.SOLID_PHASE_ONLY .AND. ANY(SPECIES%RADCAL_ID/='null')) THEN

   ! Check for valid RADCAL species and setup arrays from ZZ to RADCAL_YY

   N_RADCAL_ARRAY_SIZE = 0
   GET_RADCAL_SPECIES: DO NS=1,N_SPECIES
      SELECT CASE (SPECIES(NS)%RADCAL_ID)
         CASE('CARBON DIOXIDE')
            IF (RADCAL_TEMP(1)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=1
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='CARBON DIOXIDE'
               RADCAL_TEMP(1)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(1)
         CASE('WATER VAPOR')
            IF (RADCAL_TEMP(2)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=2
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='WATER VAPOR'
               SPECIES(NS)%RADCAL_INDEX=N_RADCAL_ARRAY_SIZE
               RADCAL_TEMP(2)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(2)
         CASE('CARBON MONOXIDE')
            IF (RADCAL_TEMP(3)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=3
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='CARBON MONOXIDE'
               RADCAL_TEMP(3)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(3)
         CASE('METHANE')
            IF (RADCAL_TEMP(4)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=4
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='METHANE'
               RADCAL_TEMP(4)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(4)
         CASE('ETHYLENE')
            IF (RADCAL_TEMP(5)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=5
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='ETHYLENE'
               RADCAL_TEMP(5)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(5)
         CASE('ETHANE')
            IF (RADCAL_TEMP(6)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=6
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='ETHANE'
               RADCAL_TEMP(6)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(6)
         CASE('PROPYLENE')
            IF (RADCAL_TEMP(7)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=7
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='PROPYLENE'
               RADCAL_TEMP(7)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(7)
         CASE('PROPANE')
            IF (RADCAL_TEMP(8)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=8
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='PROPANE'
               RADCAL_TEMP(8)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(8)
         CASE('TOLUENE')
            IF (RADCAL_TEMP(9)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=9
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='TOLUENE'
               RADCAL_TEMP(9)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(9)
         CASE('N-HEPTANE')
            IF (RADCAL_TEMP(10)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=10
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='N-HEPTANE'
               RADCAL_TEMP(10)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(10)
         CASE('METHANOL')
            IF (RADCAL_TEMP(11)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=11
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='METHANOL'
               RADCAL_TEMP(11)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(11)
         CASE('MMA')
            IF (RADCAL_TEMP(12)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=12
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='MMA'
               RADCAL_TEMP(12)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(12)
         CASE('SOOT')
            IF (RADCAL_TEMP(16)==0) THEN
               N_RADCAL_ARRAY_SIZE=N_RADCAL_ARRAY_SIZE+1
               RADCAL_SPECIES_INDEX(N_RADCAL_ARRAY_SIZE)=16
               RADCAL_SPECIES_ID(N_RADCAL_ARRAY_SIZE)='SOOT'
               RADCAL_TEMP(16)=N_RADCAL_ARRAY_SIZE
            ENDIF
            SPECIES(NS)%RADCAL_INDEX=RADCAL_TEMP(16)
      END SELECT
   END DO GET_RADCAL_SPECIES

   BUILD_KAPPA_ARRAY: IF (N_RADCAL_ARRAY_SIZE>0) THEN

      KAPPA_ARRAY=.TRUE.

      ALLOCATE(Z2RADCAL_SPECIES(N_RADCAL_ARRAY_SIZE,1:N_TRACKED_SPECIES),STAT=IZERO)
      CALL ChkMemErr('RADI','ZZ2RADCAL_SPECIES',IZERO)
      Z2RADCAL_SPECIES = 0._EB

      DO NS=1,N_TRACKED_SPECIES
         DO NS2=1,N_SPECIES
            IF (SPECIES(NS2)%RADCAL_INDEX > 0) THEN
               IF (SPECIES(NS2)%RADCAL_ID/='SOOT') THEN
                  Z2RADCAL_SPECIES(SPECIES(NS2)%RADCAL_INDEX,NS) = Z2RADCAL_SPECIES(SPECIES(NS2)%RADCAL_INDEX,NS) + &
                                                                   REAL(N_KAPPA_Y,EB)**4 / SPECIES(NS2)%MW * Z2Y(NS2,NS)
               ELSE
                  Z2RADCAL_SPECIES(SPECIES(NS2)%RADCAL_INDEX,NS) = Z2RADCAL_SPECIES(SPECIES(NS2)%RADCAL_INDEX,NS) + &
                                                                   REAL(N_KAPPA_Y,EB)**4 * 5._EB * Z2Y(NS2,NS) * &
                                                                   SPECIES(SOOT_INDEX)%DENSITY_SOLID/SPECIES(NS2)%DENSITY_SOLID
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      ! Allocate arrays for RadCal

      CALL RCALLOC

      ! Set the Mean Beam Length to 10 cm unless the user desires otherwise

      IF (PATH_LENGTH < 0._EB) PATH_LENGTH = 0.1_EB
      ALLOCATE(SEGMENT_LENGTH_M(1))
      ALLOCATE(TOTAL_PRESSURE_ATM(1))
      ALLOCATE(TEMP_GAS(1))
      ALLOCATE(PARTIAL_PRESSURES_ATM(16,1))
      SEGMENT_LENGTH_M(1) = MAX(PATH_LENGTH,1.0E-4_EB)
      TOTAL_PRESSURE_ATM(1)=P_INF/P_STP
      NPT = 1
      TWALL = RADTMP
      LAMBDAMIN = -1.1E+4_EB
      LAMBDAMAX = -1.0E+4_EB

      ! Using RadCal, create look-up tables for the absorption coefficients for all gas species, mixture fraction or aerosols

      ALLOCATE (RADCAL_SPECIES2KAPPA(N_RADCAL_ARRAY_SIZE,0:N_KAPPA_Y,0:N_KAPPA_T,NSB),STAT=IZERO)
      CALL ChkMemErr('RADI','RADCAL_SPECIES2KAPPA',IZERO)
      RADCAL_SPECIES2KAPPA = 0._EB
      BBF = 1._EB
      OMMIN = 50._EB
      OMMAX = 10000._EB
      BAND_LOOP_Z: DO IBND = 1,NSB
         IF (NSB>1) THEN
            OMMIN = REAL(NINT(1.E4_EB/WL_HIGH(IBND)),EB)
            OMMAX = REAL(NINT(1.E4_EB/WL_LOW(IBND)),EB)
         ENDIF
         CALL INIT_RADCAL
         T_LOOP_Z: DO K = 0,N_KAPPA_T
            TEMP_GAS(1) = RTMPMIN + K*(RTMPMAX-RTMPMIN)/N_KAPPA_T
            ! AMEAN will not be calculated close to RADTMP, where it cannot be solved
            IF (ABS(TEMP_GAS(1)-RADTMP)<=0.4_EB*(RTMPMAX-RTMPMIN)/N_KAPPA_T .AND. PATH_LENGTH > 0.0_EB) THEN
               RCT_SKIP = K
               CYCLE T_LOOP_Z
            ENDIF
            RCRHO = MW_AIR*P_INF/(R0*TEMP_GAS(1))
            IF (NSB>1) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TEMP_GAS(1))
            Y_LOOP_Z: DO J=1,N_KAPPA_Y
               YY = (REAL(J,EB)/REAL(N_KAPPA_Y,EB))**4
               YY2 = 1._EB-YY
               X_N2 = YY2/28._EB
               N = 0
               RADCAL_SPECIES_LOOP: DO NS = 1, N_RADCAL_ARRAY_SIZE
                  PARTIAL_PRESSURES_ATM = 0._EB
                  SELECT CASE(RADCAL_SPECIES_INDEX(NS))
                     CASE(1) ! CARBON DIOXIDE
                        PARTIAL_PRESSURES_ATM(1,1) = YY/(YY+44._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(1,1)
                     CASE(2) ! WATER VAPOR
                        PARTIAL_PRESSURES_ATM(2,1) = YY/(YY+18._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(2,1)
                     CASE(3) ! CARBON MONOXIDE
                        PARTIAL_PRESSURES_ATM(3,1) = YY/(YY+28._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(3,1)
                     CASE(4) ! METHANE
                        PARTIAL_PRESSURES_ATM(4,1) = YY/(YY+16._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(4,1)
                     CASE(5) ! EHTYLENE
                        PARTIAL_PRESSURES_ATM(5,1) = YY/(YY+28._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(5,1)
                     CASE(6) ! ETHANE
                        PARTIAL_PRESSURES_ATM(6,1) = YY/(YY+30._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(6,1)
                     CASE(7) ! PROPYLENE
                        PARTIAL_PRESSURES_ATM(7,1) = YY/(YY+42._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(7,1)
                     CASE(8) ! PROPANE
                        PARTIAL_PRESSURES_ATM(8,1) = YY/(YY+44._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(8,1)
                     CASE(9) ! TOLUENE
                        PARTIAL_PRESSURES_ATM(9,1) = YY/(YY+92._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(9,1)
                     CASE(10) ! N-HEPTANE
                        PARTIAL_PRESSURES_ATM(10,1) = YY/(YY+100._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(10,1)
                     CASE(11) ! METHANOL
                        PARTIAL_PRESSURES_ATM(11,1) = YY/(YY+32._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(11,1)
                     CASE(12) ! MMA
                        PARTIAL_PRESSURES_ATM(12,1) = YY/(YY+100._EB*X_N2)
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB - PARTIAL_PRESSURES_ATM(12,1)
                     CASE(16) ! SOOT
                        YY2 = 0.2_EB*YY
                        PARTIAL_PRESSURES_ATM(16,1) = YY2*RCRHO/SPECIES(SOOT_INDEX)%DENSITY_SOLID
                        PARTIAL_PRESSURES_ATM(14,1) = 1._EB
                  END SELECT
                  CALL SUB_RADCAL(AMEAN,AP0,RADIANCE,TRANSMISSIVITY)
                  IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
                     RADCAL_SPECIES2KAPPA(NS,J,K,1) = MIN(AMEAN,AP0)
                  ELSE
                     RADCAL_SPECIES2KAPPA(NS,J,K,IBND) = AP0/BBF
                  ENDIF
               END DO RADCAL_SPECIES_LOOP
            ENDDO Y_LOOP_Z
         ENDDO T_LOOP_Z
         ! Interpolate KAPPA at RADTMP
         IF (RCT_SKIP == 0) THEN
            RADCAL_SPECIES2KAPPA(:,:,0,IBND) = RADCAL_SPECIES2KAPPA(:,:,1,IBND)
         ELSEIF (RCT_SKIP == N_KAPPA_T) THEN
            RADCAL_SPECIES2KAPPA(:,:,N_KAPPA_T,IBND) = RADCAL_SPECIES2KAPPA(:,:,N_KAPPA_T-1,IBND)
         ELSEIF (RCT_SKIP > 0) THEN
            RADCAL_SPECIES2KAPPA(:,:,RCT_SKIP,IBND) = 0.5_EB*(RADCAL_SPECIES2KAPPA(:,:,RCT_SKIP-1,IBND)+ &
                                                              RADCAL_SPECIES2KAPPA(:,:,RCT_SKIP+1,IBND))
         ENDIF
         CALL RCDEALLOC2  ! Deallocate RadCal wavelength dependent arrays
      ENDDO BAND_LOOP_Z
      !Adjust values from /cm to /m
      RADCAL_SPECIES2KAPPA =  RADCAL_SPECIES2KAPPA * 100._EB
      CALL RCDEALLOC  ! Deallocate RadCal arrays

   ENDIF BUILD_KAPPA_ARRAY

   ! Trap any errors

   IF (ANY(RADCAL_SPECIES2KAPPA<0._EB)) CALL SHUTDOWN('ERROR: KAPPA < 0 in RADCAL')

ENDIF MAKE_KAPPA_ARRAYS

! Tables for PARTICLE absorption coefficients

DO IPC=1,N_LAGRANGIAN_CLASSES
   LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
   IF (LPC%LIQUID_DROPLET) CALL MEAN_CROSS_SECTIONS(IPC)
ENDDO

! Determine angle factors for Lagrangian particles with ORIENTATION

IF (SOLID_PARTICLES) THEN
   ALLOCATE(ORIENTATION_FACTOR(NRA,N_ORIENTATION_VECTOR))
   ORIENTATION_FACTOR = 0._EB
   PARTICLE_CLASS_LOOP: DO IPC=1,N_LAGRANGIAN_CLASSES
      LPC => LAGRANGIAN_PARTICLE_CLASS(IPC)
      IF (LPC%N_ORIENTATION==0) CYCLE PARTICLE_CLASS_LOOP
      I1 = LPC%ORIENTATION_INDEX
      I2 = LPC%ORIENTATION_INDEX+LPC%N_ORIENTATION-1
      ALLOCATE(COSINE_ARRAY(I1:I2))
      ALLOCATE(COSINE_ARRAY_2(NRA))
      LPC%SOLID_ANGLE = 0._EB
      ANGLE_LOOP: DO N=1,NRA
         ORIENTATION_LOOP: DO OR_IN=I1,I2
            COSINE_ARRAY(OR_IN) = ORIENTATION_VECTOR(1,OR_IN)*DLX(N) + &
                                  ORIENTATION_VECTOR(2,OR_IN)*DLY(N) + &
                                  ORIENTATION_VECTOR(3,OR_IN)*DLZ(N)
         ENDDO ORIENTATION_LOOP
         COSINE_ARRAY_2(N) = COSINE_ARRAY(I1)
         ORIENTATION_FACTOR(N,I1+MINLOC(COSINE_ARRAY)-1) = 1._EB
         LPC%SOLID_ANGLE(MINLOC(COSINE_ARRAY)) = LPC%SOLID_ANGLE(MINLOC(COSINE_ARRAY)) + RSA(N)
      ENDDO ANGLE_LOOP
      DO N=1,NRA
         ORIENTATION_FACTOR(N,I1:I2) = ORIENTATION_FACTOR(N,I1:I2)*PI/LPC%SOLID_ANGLE(1:I2-I1+1)
      ENDDO
      LPC%NEAREST_RAD_ANGLE_INDEX = MAXLOC(COSINE_ARRAY_2,DIM=1)
      DEALLOCATE(COSINE_ARRAY)
      DEALLOCATE(COSINE_ARRAY_2)
   ENDDO PARTICLE_CLASS_LOOP
ENDIF

! Initialize radiation file (RADF)

DO NM=1,NMESHES

   IF (MYID/=PROCESS(NM)) CYCLE

   CALL POINT_TO_MESH(NM)

   DO N=1,N_RADF
      RF => RAD_FILE(N)
      ALLOCATE(RF%IL_SAVE(RF%I1:RF%I2,RF%J1:RF%J2,RF%K1:RF%K2,NUMBER_RADIATION_ANGLES))
      IF (.NOT.APPEND) THEN
         OPEN(LU_RADF(N,NM),FILE=FN_RADF(N,NM),FORM='FORMATTED',STATUS='REPLACE')
         WRITE(LU_RADF(N,NM),'(A)') 'NSTEPS'
         NSTEPS = INT((T_RADF_END-T_RADF_BEGIN)/DT_RADF) + 1
         WRITE(LU_RADF(N,NM),'(I4)') NSTEPS
         WRITE(LU_RADF(N,NM),'(/A)') 'TIMES'
         DO NN=1,NSTEPS
            WRITE(LU_RADF(N,NM),'(F8.2)') T_RADF_BEGIN + (NN-1)*DT_RADF
         ENDDO
         WRITE(LU_RADF(N,NM),'(/A)') 'NP'
         WRITE(LU_RADF(N,NM),'(I8)') RF%N_POINTS
         WRITE(LU_RADF(N,NM),'(/A)') 'XYZ_INTENSITIES'
         DO K=RF%K1,RF%K2,RF%K_STEP
            DO J=RF%J1,RF%J2,RF%J_STEP
               DO I=RF%I1,RF%I2,RF%I_STEP
                  WRITE(LU_RADF(N,NM),'(3F8.3)') XC(I),YC(J),ZC(K)
               ENDDO
            ENDDO
         ENDDO
         WRITE(LU_RADF(N,NM),'(/A)') 'NI'
         WRITE(LU_RADF(N,NM),'(I4)') NUMBER_RADIATION_ANGLES
         WRITE(LU_RADF(N,NM),'(/A)') 'XYZ_DIRECTIONS'
         DO NN=1,NRA
            NRM = NORM2([DLX(NN),DLY(NN),DLZ(NN)])
            WRITE(LU_RADF(N,NM),'(3F7.3)') DLX(NN)/NRM,DLY(NN)/NRM,DLZ(NN)/NRM
         ENDDO
      ELSE
         OPEN(LU_RADF(N,NM),FILE=FN_RADF(N,NM),FORM='FORMATTED',STATUS='OLD')
      ENDIF
   ENDDO

ENDDO

END SUBROUTINE INIT_RADIATION



SUBROUTINE COMPUTE_RADIATION(T,NM,RAD_ITER)

! Call radiation routine or simply specify the radiative loss

USE COMP_FUNCTIONS, ONLY : CURRENT_TIME
USE COMPLEX_GEOMETRY
REAL(EB) :: TNOW,T
INTEGER, INTENT(IN) :: NM,RAD_ITER

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=CURRENT_TIME()

CALL POINT_TO_MESH(NM)

IF (RADIATION) THEN
   RADIATION_COMPLETED = .FALSE.
   CALL RADIATION_FVM(T,NM,RAD_ITER)
ELSE
   RADIATION_COMPLETED = .TRUE.
   IF (N_REACTIONS>0) QR = -CHI_R*Q
ENDIF

T_USED(9)=T_USED(9)+CURRENT_TIME()-TNOW

CONTAINS


SUBROUTINE RADIATION_FVM(T,NM,RAD_ITER)
USE MIEV
USE MATH_FUNCTIONS, ONLY : INTERPOLATE1D, EVALUATE_RAMP
USE TRAN, ONLY : GET_IJK
USE COMPLEX_GEOMETRY, ONLY : IBM_IDRA,IBM_CGSC,IBM_SOLID
USE MPI
REAL(EB) :: T, RAP, AX, AXU, AXD, AY, AYU, AYD, AZ, VC, RU, RD, RP, &
            ILXU, ILYU, ILZU, QVAL, BBF, BBFA, NCSDROP, RSA_RAT,EFLUX,TYY_FAC, &
            AIU_SUM,A_SUM,VOL,VC1,AY1,AZ1,COSINE,AFX,AFY,AFZ,ILXU_AUX,ILYU_AUX,ILZU_AUX,AFX_AUX,AFY_AUX,AFZ_AUX
INTEGER  :: N,NN,IIG,JJG,KKG,I,J,K,IW,ICF,II,JJ,KK,IOR,IC,IWUP,IWDOWN, &
            ISTART, IEND, ISTEP, JSTART, JEND, JSTEP, &
            KSTART, KEND, KSTEP, NSTART, NEND, NSTEP, &
            I_UIID, N_UPDATES, IBND, TYY, NOM, ARRAY_INDEX,NRA, &
            IMIN, JMIN, KMIN, IMAX, JMAX, KMAX, N_SLICE, M_IJK, IJK, LL
INTEGER  :: IADD,ICR,IFA
INTEGER, ALLOCATABLE :: IJK_SLICE(:,:)
REAL(EB) :: XID,YJD,ZKD,AREA_VOLUME_RATIO,DLF,DLA(3),TSI,TMP_EXTERIOR,DUMMY
REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_GET
INTEGER :: IID,JJD,KKD,IP
LOGICAL :: UPDATE_INTENSITY, UPDATE_QRW2
REAL(EB), POINTER, DIMENSION(:,:,:) :: IL,UIIOLD,KAPPA_PART,KFST4_GAS,KFST4_PART,EXTCOE,SCAEFF,IL_UP
REAL(EB), POINTER, DIMENSION(:)     :: OUTRAD_W,INRAD_W,OUTRAD_F,INRAD_F
INTEGER, INTENT(IN) :: NM,RAD_ITER
TYPE (OMESH_TYPE), POINTER :: M2
TYPE(SURFACE_TYPE), POINTER :: SF
TYPE(VENTS_TYPE), POINTER :: VT
TYPE(RAD_FILE_TYPE), POINTER :: RF
TYPE(LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC
TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER :: LP
CHARACTER(20) :: FORMT
ALLOCATE( IJK_SLICE(3, IBAR*KBAR) )

KFST4_GAS => WORK1
IL       => WORK2
UIIOLD   => WORK3
EXTCOE   => WORK4
KAPPA_PART => WORK5
SCAEFF   => WORK6
KFST4_PART   => WORK7
IL_UP    => WORK8
OUTRAD_W => WALL_WORK1
INRAD_W  => WALL_WORK2
IF (CC_IBM) THEN
   OUTRAD_F => FACE_WORK1
   INRAD_F  => FACE_WORK2
   INRAD_F  = 0._EB
ENDIF

! Ratio of solid angle, used in scattering

NRA     = NUMBER_RADIATION_ANGLES
RSA_RAT = 1._EB/(1._EB-1._EB/NRA)

! Check if it time to update radiation intensity field

IF ( MOD(RAD_CALL_COUNTER,TIME_STEP_INCREMENT)==0 .OR. INITIALIZATION_PHASE .OR. ICYC==1) THEN
   UPDATE_INTENSITY   = .TRUE.
   EXCHANGE_RADIATION = .TRUE.
ELSE
   UPDATE_INTENSITY   = .FALSE.
   EXCHANGE_RADIATION = .FALSE.
ENDIF

IF (RAD_ITER==RADIATION_ITERATIONS) RAD_CALL_COUNTER  = RAD_CALL_COUNTER + 1

IF (WIDE_BAND_MODEL) THEN
   QR = 0._EB
ENDIF

! Zero out radiation flux to wall, particles, facets if the intensity is to be updated

IF (UPDATE_INTENSITY) THEN
   DO IW=1,N_INTERNAL_WALL_CELLS+N_EXTERNAL_WALL_CELLS
      WALL(IW)%ONE_D%Q_RAD_IN = 0._EB
   ENDDO
   DO IP=1,NLP
      LAGRANGIAN_PARTICLE(IP)%ONE_D%Q_RAD_IN = 0._EB
   ENDDO
   DO ICF=1,N_CFACE_CELLS
      CFACE(ICF)%ONE_D%Q_RAD_IN = 0._EB
   ENDDO
ENDIF

UPDATE_QRW2 = .FALSE.

! Loop over spectral bands

BAND_LOOP: DO IBND = 1,NUMBER_SPECTRAL_BANDS

   KAPPA_PART = 0._EB
   KFST4_GAS  = 0._EB
   KFST4_PART = 0._EB
   SCAEFF = 0._EB

   ! Calculate fraction on ambient black body radiation

   IF (NUMBER_SPECTRAL_BANDS==1) THEN
      BBFA = 1._EB
   ELSE
      BBFA = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMPA)
   ENDIF

   ! Generate water absorption and scattering coefficients

   IF_PARTICLES_INCLUDED: IF (NLP>0 .AND. N_LP_ARRAY_INDICES>0) THEN

      IF (NUMBER_SPECTRAL_BANDS==1) THEN
         BBF = 1._EB
      ELSE
         BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),RADTMP)
      ENDIF

      PC_LOOP: DO N=1,N_LAGRANGIAN_CLASSES
         LPC => LAGRANGIAN_PARTICLE_CLASS(N)
         IF (.NOT.LPC%LIQUID_DROPLET) CYCLE PC_LOOP
         ARRAY_INDEX = LPC%ARRAY_INDEX
         IF (ARRAY_INDEX==0) CYCLE PC_LOOP
         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  IF (ABS(AVG_DROP_AREA(I,J,K,ARRAY_INDEX))<TWO_EPSILON_EB) CYCLE
                  NCSDROP = AVG_DROP_AREA(I,J,K,ARRAY_INDEX)
                  CALL INTERPOLATE1D(LPC%R50,LPC%WQABS(:,IBND),AVG_DROP_RAD(I,J,K,ARRAY_INDEX),QVAL)
                  KAPPA_PART(I,J,K) = KAPPA_PART(I,J,K) + NCSDROP*QVAL
                  KFST4_PART(I,J,K) = KFST4_PART(I,J,K)+ BBF*NCSDROP*QVAL*FOUR_SIGMA*AVG_DROP_TMP(I,J,K,ARRAY_INDEX)**4
                  CALL INTERPOLATE1D(LPC%R50,LPC%WQSCA(:,IBND),AVG_DROP_RAD(I,J,K,ARRAY_INDEX),QVAL)
                  SCAEFF(I,J,K) = SCAEFF(I,J,K) + NCSDROP*QVAL
               ENDDO
            ENDDO
         ENDDO
      ENDDO PC_LOOP

      QR_W = 0._EB

   ENDIF IF_PARTICLES_INCLUDED

   ! Compute the absorption coefficient, KAPPA_PART, for a collection of solid particles

   IF (NLP>0 .AND. SOLID_PARTICLES) THEN
      DO IP = 1,NLP
         LP => LAGRANGIAN_PARTICLE(IP)
         LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
         IF (.NOT.LPC%SOLID_PARTICLE) CYCLE
         CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XID,YJD,ZKD,IID,JJD,KKD)
         AREA_VOLUME_RATIO = LPC%SHAPE_FACTOR*LP%PWT*LP%ONE_D%AREA*RDX(IID)*RRN(IID)*RDY(JJD)*RDZ(KKD)
         KAPPA_PART(IID,JJD,KKD) = KAPPA_PART(IID,JJD,KKD) + AREA_VOLUME_RATIO
         KFST4_PART(IID,JJD,KKD) = KFST4_PART(IID,JJD,KKD) + FOUR_SIGMA*AREA_VOLUME_RATIO*LP%ONE_D%TMP_F**4
      ENDDO
   ENDIF

   ! Compute absorption coefficient of the gases, KAPPA_GAS

   IF (KAPPA0>=0._EB) THEN

      KAPPA_GAS = KAPPA0

   ELSEIF (KAPPA_ARRAY) THEN

      TYY_FAC = N_KAPPA_T / (RTMPMAX-RTMPMIN)
      !$OMP PARALLEL PRIVATE(ZZ_GET, TYY)
      ALLOCATE(ZZ_GET(1:N_TRACKED_SPECIES))
      !$OMP DO SCHEDULE(STATIC)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               TYY = MAX(0 , MIN(N_KAPPA_T,INT((TMP(I,J,K) - RTMPMIN) * TYY_FAC)))
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
               KAPPA_GAS(I,J,K) = GET_KAPPA(ZZ_GET,TYY,IBND)
            ENDDO
         ENDDO
      ENDDO
      !$OMP END DO
      DEALLOCATE(ZZ_GET)
      !$OMP END PARALLEL
   ENDIF

   ! Compute source term KAPPA*4*SIGMA*TMP**4

   WIDE_BAND_MODEL_IF: IF (WIDE_BAND_MODEL) THEN

      ! Wide band model

      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMP(I,J,K))
               KFST4_GAS(I,J,K) = BBF*KAPPA_GAS(I,J,K)*FOUR_SIGMA*TMP(I,J,K)**4
            ENDDO
         ENDDO
      ENDDO

   ELSE WIDE_BAND_MODEL_IF

      ! Gray gas model

      RTE_SOURCE_CORRECTION_IF: IF (RTE_SOURCE_CORRECTION) THEN ! default RTE_SOURCE_CORRECTION=.TRUE.

         ! Only apply the correction to KFST4_GAS for gray gas model

         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  KFST4_GAS(I,J,K) = KAPPA_GAS(I,J,K)*FOUR_SIGMA*TMP(I,J,K)**4
                  IF (CHI_R(I,J,K)*Q(I,J,K)>QR_CLIP) THEN
                     VOL = R(I)*DX(I)*DY(J)*DZ(K)
                     RAD_Q_SUM = RAD_Q_SUM + (CHI_R(I,J,K)*Q(I,J,K)+KAPPA_GAS(I,J,K)*UII(I,J,K))*VOL
                     KFST4_SUM = KFST4_SUM + KFST4_GAS(I,J,K)*VOL
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         ! Correct the source term in the RTE based on user-specified RADIATIVE_FRACTION on REAC

         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  IF (CHI_R(I,J,K)*Q(I,J,K)>QR_CLIP) KFST4_GAS(I,J,K) = KFST4_GAS(I,J,K)*RTE_SOURCE_CORRECTION_FACTOR
               ENDDO
            ENDDO
         ENDDO

      ELSE RTE_SOURCE_CORRECTION_IF  ! OPTICALLY_THIN

         ! Use specified radiative fraction

         DO K=1,KBAR
            DO J=1,JBAR
               DO I=1,IBAR
                  IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
                  KFST4_GAS(I,J,K) = CHI_R(I,J,K)*Q(I,J,K)+KAPPA_GAS(I,J,K)*UII(I,J,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF RTE_SOURCE_CORRECTION_IF

   ENDIF WIDE_BAND_MODEL_IF

   ! Calculate extinction coefficient

   EXTCOE = KAPPA_GAS + KAPPA_PART + SCAEFF*RSA_RAT

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
      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         IF (WALL(IW)%BOUNDARY_TYPE == OPEN_BOUNDARY) THEN
            BBF = BBFA
         ELSE
            IF (WIDE_BAND_MODEL) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),WALL(IW)%ONE_D%TMP_F)
            SF  => SURFACE(WALL(IW)%SURF_INDEX)
            IF (.NOT. SF%INTERNAL_RADIATION) WALL(IW)%ONE_D%Q_RAD_OUT = WALL(IW)%ONE_D%EMISSIVITY*SIGMA*WALL(IW)%ONE_D%TMP_F**4
         ENDIF
         OUTRAD_W(IW) = BBF*RPI*WALL(IW)%ONE_D%Q_RAD_OUT
      ENDDO

      BBF = 1.0_EB
      DO ICF = 1,N_CFACE_CELLS
         IF (WIDE_BAND_MODEL) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),CFACE(ICF)%ONE_D%TMP_F)
         SF => SURFACE(CFACE(ICF)%SURF_INDEX)
         IF (.NOT. SF%INTERNAL_RADIATION) CFACE(ICF)%ONE_D%Q_RAD_OUT = CFACE(ICF)%ONE_D%EMISSIVITY*SIGMA*CFACE(ICF)%ONE_D%TMP_F**4
         OUTRAD_F(ICF) = BBF*RPI*CFACE(ICF)%ONE_D%Q_RAD_OUT
      ENDDO

      ! Compute boundary condition term incoming radiation integral

      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
         INRAD_W(IW) = SUM(-DLN(WALL(IW)%ONE_D%IOR,:)*WALL(IW)%ONE_D%BAND(IBND)%ILW(:), 1, DLN(WALL(IW)%ONE_D%IOR,:)<0._EB)
      ENDDO

      DO ICF = 1,N_CFACE_CELLS
         DO N=1,NRA
            DLA = (/DLX(N),DLY(N),DLZ(N)/)
            DLF = DOT_PRODUCT(CFACE(ICF)%NVEC,DLA) ! face normal * radiation angle
            IF (DLF<0._EB) INRAD_F(ICF) = INRAD_F(ICF) - DLF*CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N)
         ENDDO
      ENDDO

      ! If updating intensities first time, sweep ALL angles

      N_UPDATES = 1
      IF (INITIALIZATION_PHASE .OR. ICYC==1) N_UPDATES = ANGLE_INCREMENT

      UPDATE_LOOP: DO I_UIID = 1,N_UPDATES

         ! Update counters inside the radiation routine

         ANGLE_INC_COUNTER = MOD(ANGLE_INC_COUNTER,ANGLE_INCREMENT) + 1

         ! If this is the last set of angles to update, indicate that the radiation routine has finished a full update

         IF (ANGLE_INC_COUNTER==ANGLE_INCREMENT) RADIATION_COMPLETED = .TRUE.

         ! Zero out UIID, the integrated intensity

         IF (WIDE_BAND_MODEL) THEN
            UIID(:,:,:,IBND) = 0._EB
         ELSE
            UIID(:,:,:,ANGLE_INC_COUNTER) = 0._EB
         ENDIF

         DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
            IF (WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) WALL(IW)%ONE_D%BAND(IBND)%ILW(ANGLE_INC_COUNTER) = 0._EB
         ENDDO

         ! Set the bounds and increment for the angleloop. Step downdard because in cylindrical case the Nth angle
         ! boundary condition comes from (N+1)th angle.

         NSTART    = NRA - ANGLE_INC_COUNTER + 1
         NEND      = 1
         NSTEP     = -ANGLE_INCREMENT

         IL(:,:,:) = BBFA*RPI_SIGMA*TMPA4

         ANGLE_LOOP: DO N = NSTART,NEND,NSTEP  ! Sweep through control angles

            ! Boundary conditions: Intensities leaving the boundaries.

            !$OMP PARALLEL DO PRIVATE(IOR, II, JJ, KK, LL, NOM) SCHEDULE(GUIDED)
            WALL_LOOP1: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
               IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE WALL_LOOP1
               IOR = WALL(IW)%ONE_D%IOR
               IF (DLN(IOR,N) < 0._EB) CYCLE WALL_LOOP1
               II  = WALL(IW)%ONE_D%II
               JJ  = WALL(IW)%ONE_D%JJ
               KK  = WALL(IW)%ONE_D%KK
               IF (.NOT.TWO_D .OR. ABS(IOR)/=2) THEN
                  SELECT CASE (WALL(IW)%BOUNDARY_TYPE)
                     CASE (OPEN_BOUNDARY)
                        VT => VENTS(WALL(IW)%VENT_INDEX)
                        IF (VT%TMP_EXTERIOR>0._EB) THEN
                           TSI = T - T_BEGIN
                           TMP_EXTERIOR = TMP_0(KK)+EVALUATE_RAMP(TSI,DUMMY,VT%TMP_EXTERIOR_RAMP_INDEX)*(VT%TMP_EXTERIOR-TMP_0(KK))
                           IL(II,JJ,KK) = BBFA*RPI_SIGMA*TMP_EXTERIOR**4
                        ELSE
                           IL(II,JJ,KK) = BBFA*RPI_SIGMA*TMPA4
                        ENDIF
                     CASE (MIRROR_BOUNDARY)
                        WALL(IW)%ONE_D%BAND(IBND)%ILW(N) = WALL(IW)%ONE_D%BAND(IBND)%ILW(DLM(N,ABS(IOR)))
                        IL(II,JJ,KK) = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                     CASE (INTERPOLATED_BOUNDARY)
                        ! IL_R holds the intensities from mesh NOM in the ghost cells of mesh NM.
                        ! IL(II,JJ,KK) is the average of the intensities from the other mesh.
                        NOM = EXTERNAL_WALL(IW)%NOM
                        IL(II,JJ,KK) = 0._EB
                        DO LL=EXTERNAL_WALL(IW)%NIC_MIN,EXTERNAL_WALL(IW)%NIC_MAX
                           IL(II,JJ,KK) = IL(II,JJ,KK) + OMESH(NOM)%IL_R(LL,N,IBND)
                        ENDDO
                        IL(II,JJ,KK) = IL(II,JJ,KK)/REAL(EXTERNAL_WALL(IW)%NIC_MAX-EXTERNAL_WALL(IW)%NIC_MIN+1,EB)
                     CASE DEFAULT ! solid wall
                        WALL(IW)%ONE_D%BAND(IBND)%ILW(N) = OUTRAD_W(IW) + RPI*(1._EB-WALL(IW)%ONE_D%EMISSIVITY)*INRAD_W(IW)
                  END SELECT
               ELSEIF (CYLINDRICAL) THEN
                  IF (WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY) CYCLE WALL_LOOP1
                  IL(II,JJ,KK) = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
               ENDIF
            ENDDO WALL_LOOP1
            !$OMP END PARALLEL DO

            CFACE_LOOP1: DO ICF=1,N_CFACE_CELLS
               IF (CFACE(ICF)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE CFACE_LOOP1
               CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N) = OUTRAD_F(ICF) + RPI*(1._EB-CFACE(ICF)%ONE_D%EMISSIVITY)*INRAD_F(ICF)
            ENDDO CFACE_LOOP1

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
            IMIN = ISTART
            JMIN = JSTART
            KMIN = KSTART
            IMAX = IEND
            JMAX = JEND
            KMAX = KEND
            IF (DLX(N) < 0._EB) THEN
               ISTART = IBAR
               IEND   = 1
               ISTEP  = -1
               IMIN = IEND
               IMAX = ISTART
            ENDIF
            IF (DLY(N) < 0._EB) THEN
               JSTART = JBAR
               JEND   = 1
               JSTEP  = -1
               JMIN = JEND
               JMAX = JSTART
            ENDIF
            IF (DLZ(N) < 0._EB) THEN
               KSTART = KBAR
               KEND   = 1
               KSTEP  = -1
               KMIN = KEND
               KMAX = KSTART
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
                     AXU = 2._EB*SIN(DPHI0/2.)*RU       * DZ(K) * ABS(DLX(N))
                     AXD = 2._EB*SIN(DPHI0/2.)*RD       * DZ(K) * ABS(DLX(N))
                     AYU = DX(I)             * DZ(K) * ABS(DLB(N))
                     AYD = DX(I)             * DZ(K) * ABS(DLY(N))
                     AZ  = DX(I)  * RP*DPHI0         * ABS(DLZ(N))
                     IF (MODULO(N,NRP(1))==1) AYD = 0._EB  ! Zero out the terms involving symmetric overhang
                     IF (MODULO(N,NRP(1))==0) AYU = 0._EB
                     IF (IC/=0) THEN
                        IW = WALL_INDEX(IC,-ISTEP)
                        IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILXU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                        IW = WALL_INDEX(IC,-JSTEP*2)
                        IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILYU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                        IW = WALL_INDEX(IC,-KSTEP*3)
                        IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILZU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                     ENDIF
                     AIU_SUM = AXU*ILXU + AYU*ILYU + AZ*ILZU
                     A_SUM = AXD + AYD + AZ
                     RAP = 1._EB/(A_SUM + EXTCOE(I,J,K)*VC*RSA(N))
                     IL(I,J,K) = MAX(0._EB, RAP * (AIU_SUM + VC*RSA(N)*RFPI* &
                                 ( KFST4_GAS(I,J,K) + KFST4_PART(I,J,K) + RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
                     IF (SOLID_PARTICLES) IL_UP(I,J,K) = MAX(0._EB,AIU_SUM/A_SUM)
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
                        IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILXU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                        IW = WALL_INDEX(IC,-KSTEP*3)
                        IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILZU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                     ENDIF
                     AIU_SUM = AX*ILXU + AZ*ILZU
                     A_SUM = AX + AZ
                     RAP = 1._EB/(A_SUM + EXTCOE(I,J,K)*VC*RSA(N))
                     IL(I,J,K) = MAX(0._EB, RAP * (AIU_SUM + VC*RSA(N)*RFPI* &
                                    (KFST4_GAS(I,J,K) + KFST4_PART(I,J,K) + RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
                     IF (SOLID_PARTICLES) IL_UP(I,J,K) = MAX(0._EB,AIU_SUM/A_SUM)
                  ENDDO I2LOOP
               ENDDO K2LOOP

            ELSE GEOMETRY  ! Sweep in 3D cartesian geometry

              DO N_SLICE = ISTEP*ISTART + JSTEP*JSTART + KSTEP*KSTART, &
                          ISTEP*IEND + JSTEP*JEND + KSTEP*KEND
                M_IJK = 0
                DO K = KMIN, KMAX
                  IF (ISTEP*JSTEP > 0) THEN ! I STARTS HIGH
                    JSTART = MAX(JMIN, JSTEP*(N_SLICE - KSTEP*K - ISTEP*IMAX))
                    JEND   = MIN(JMAX, JSTEP*(N_SLICE - KSTEP*K - ISTEP*IMIN))
                  ELSE IF (ISTEP*JSTEP < 0) THEN ! I STARTS LOW
                    JSTART = MAX(JMIN, JSTEP*(N_SLICE - KSTEP*K - ISTEP*IMIN))
                    JEND   = MIN(JMAX, JSTEP*(N_SLICE - KSTEP*K - ISTEP*IMAX))
                  ENDIF
                  IF (JSTART > JEND) THEN
                    CYCLE
                  ENDIF
                  DO J = JSTART, JEND
                    I = ISTEP * (N_SLICE - J*JSTEP - K*KSTEP)
                    M_IJK = M_IJK+1
                    IJK_SLICE(:,M_IJK) = (/I,J,K/)
                  ENDDO
                ENDDO

                 !$OMP PARALLEL DO SCHEDULE(GUIDED) &
                 !$OMP& PRIVATE(I, J, K, AY1, AX, VC1, AZ1, IC, ILXU, ILYU, &
                 !$OMP& ILZU, VC, AY, AZ, IW, A_SUM, AIU_SUM, RAP, AFX, AFY, AFZ, &
                 !$OMP& AFX_AUX, AFY_AUX, AFZ_AUX, ILXU_AUX, ILYU_AUX, ILZU_AUX, &
                 !$OMP& ICF, IADD, ICR, IFA )
                 SLICELOOP: DO IJK = 1, M_IJK
                   I = IJK_SLICE(1,IJK)
                   J = IJK_SLICE(2,IJK)
                   K = IJK_SLICE(3,IJK)

                   AY1 = DZ(K) * ABS(DLY(N))
                   AX  = DY(J) * DZ(K) * ABS(DLX(N))
                   VC1 = DY(J) * DZ(K)
                   AZ1 = DY(J) * ABS(DLZ(N))
                   IC = CELL_INDEX(I,J,K)
                   IF (SOLID(IC)) CYCLE SLICELOOP
                   ILXU  = IL(I-ISTEP,J,K)
                   ILYU  = IL(I,J-JSTEP,K)
                   ILZU  = IL(I,J,K-KSTEP)
                   VC  = DX(I) * VC1
                   AY  = DX(I) * AY1
                   AZ  = DX(I) * AZ1
                   IF (IC/=0) THEN
                       IW = WALL_INDEX(IC,-ISTEP)
                       IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILXU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                       IW = WALL_INDEX(IC,-JSTEP*2)
                       IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILYU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                       IW = WALL_INDEX(IC,-KSTEP*3)
                       IF (WALL(IW)%BOUNDARY_TYPE==SOLID_BOUNDARY) ILZU = WALL(IW)%ONE_D%BAND(IBND)%ILW(N)
                   ENDIF
                   IF (CC_IBM) THEN
                      IF (CCVAR(I,J,K,IBM_CGSC) == IBM_SOLID) CYCLE SLICELOOP
                      AFX_AUX  = 0._EB; AFY_AUX  = 0._EB; AFZ_AUX  = 0._EB
                      ILXU_AUX = 0._EB; ILYU_AUX = 0._EB; ILZU_AUX = 0._EB
                      ! X axis
                      IADD= -(1+ISTEP)/2
                      ICR = FCVAR(I+IADD,J,K,IBM_IDRA,IAXIS) ! List of CFACES assigned to upwind X face.
                      DO IFA=1,RAD_CFACE(ICR)%N_ASSIGNED_CFACES_RADI
                         ICF=RAD_CFACE(ICR)%ASSIGNED_CFACES_RADI(IFA)
                         IF (REAL(ISTEP,EB)*CFACE(ICF)%NVEC(IAXIS)>0._EB) THEN
                            AFX      = ABS(CFACE(ICF)%NVEC(IAXIS))*CFACE(ICF)%AREA/(DY(J)*DZ(K))
                            AFX_AUX  = AFX_AUX  + AFX
                            ILXU_AUX = ILXU_AUX + CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N)*AFX
                         ENDIF
                      ENDDO
                      ! Y axis
                      IADD= -(1+JSTEP)/2
                      ICR = FCVAR(I,J+IADD,K,IBM_IDRA,JAXIS) ! List of CFACES assigned to upwind Y face.
                      DO IFA=1,RAD_CFACE(ICR)%N_ASSIGNED_CFACES_RADI
                         ICF=RAD_CFACE(ICR)%ASSIGNED_CFACES_RADI(IFA)
                         IF (REAL(JSTEP,EB)*CFACE(ICF)%NVEC(JAXIS)>0._EB) THEN
                            AFY      = ABS(CFACE(ICF)%NVEC(JAXIS))*CFACE(ICF)%AREA/(DX(I)*DZ(K))
                            AFY_AUX  = AFY_AUX  + AFY
                            ILYU_AUX = ILYU_AUX + CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N)*AFY
                         ENDIF
                      ENDDO
                      ! Z axis
                      IADD= -(1+KSTEP)/2
                      ICR = FCVAR(I,J,K+IADD,IBM_IDRA,KAXIS) ! List of CFACES assigned to upwind Z face.
                      DO IFA=1,RAD_CFACE(ICR)%N_ASSIGNED_CFACES_RADI
                         ICF=RAD_CFACE(ICR)%ASSIGNED_CFACES_RADI(IFA)
                         IF (REAL(KSTEP,EB)*CFACE(ICF)%NVEC(KAXIS)>0._EB) THEN
                            AFZ      = ABS(CFACE(ICF)%NVEC(KAXIS))*CFACE(ICF)%AREA/(DX(I)*DY(J))
                            AFZ_AUX  = AFZ_AUX  + AFZ
                            ILZU_AUX = ILZU_AUX + CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N)*AFZ
                         ENDIF
                      ENDDO
                      ILXU = ILXU*(1._EB-AFX_AUX) + ILXU_AUX
                      ILYU = ILYU*(1._EB-AFY_AUX) + ILYU_AUX
                      ILZU = ILZU*(1._EB-AFZ_AUX) + ILZU_AUX
                   ENDIF
                   A_SUM = AX + AY + AZ
                   AIU_SUM = AX*ILXU + AY*ILYU + AZ*ILZU
                   IF (SOLID_PARTICLES) IL_UP(I,J,K) = MAX(0._EB,AIU_SUM/A_SUM)
                   RAP = 1._EB/(A_SUM + EXTCOE(I,J,K)*VC*RSA(N))
                   IL(I,J,K) = MAX(0._EB, RAP * (AIU_SUM + VC*RSA(N)*RFPI* &
                                   ( KFST4_GAS(I,J,K) + KFST4_PART(I,J,K) + RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
                 ENDDO SLICELOOP
                 !$OMP END PARALLEL DO

               ENDDO ! IPROP

            ENDIF GEOMETRY

            ! Copy the Y-downwind intensities to Y-upwind in cylindrical case

            IF (CYLINDRICAL) THEN
               J=1
               CKLOOP2: DO K=1,KBAR
               CILOOP2: DO I=1,IBAR
                  IC = CELL_INDEX(I,J,K)
                  IF (SOLID(IC)) CYCLE CILOOP2
                  IWUP   = WALL_INDEX(CELL_INDEX(I,J,K),-2)
                  IWDOWN = WALL_INDEX(CELL_INDEX(I,J,K), 2)
                  IF (IWUP /=0 .AND. IWDOWN /= 0) THEN
                     IF (MODULO(N,NRP(1))==1) THEN
                        WALL(IWUP)%ONE_D%BAND(IBND)%ILW(N)   = WALL(IWDOWN)%ONE_D%BAND(IBND)%ILW(N)
                     ELSE
                        WALL(IWUP)%ONE_D%BAND(IBND)%ILW(N-1) = WALL(IWDOWN)%ONE_D%BAND(IBND)%ILW(N)
                     ENDIF
                  ENDIF
               ENDDO CILOOP2
               ENDDO CKLOOP2
            ENDIF

            ! Boundary values: Incoming radiation

            !$OMP PARALLEL PRIVATE(IOR, IIG, JJG, KKG)
            !$OMP DO SCHEDULE(GUIDED)
            WALL_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
               IF (WALL(IW)%BOUNDARY_TYPE==NULL_BOUNDARY)   CYCLE WALL_LOOP2
               IF (WALL(IW)%BOUNDARY_TYPE==OPEN_BOUNDARY)   CYCLE WALL_LOOP2
               IOR = WALL(IW)%ONE_D%IOR
               IF (TWO_D .AND. .NOT.CYLINDRICAL  .AND. ABS(IOR)==2) CYCLE WALL_LOOP2  ! 2-D non cylindrical
               IF (DLN(IOR,N)>=0._EB) CYCLE WALL_LOOP2     ! outgoing
               IIG = WALL(IW)%ONE_D%IIG
               JJG = WALL(IW)%ONE_D%JJG
               KKG = WALL(IW)%ONE_D%KKG
               INRAD_W(IW) = INRAD_W(IW) + DLN(IOR,N) * WALL(IW)%ONE_D%BAND(IBND)%ILW(N) ! update incoming radiation,step 1
               WALL(IW)%ONE_D%BAND(IBND)%ILW(N) = IL(IIG,JJG,KKG)
               INRAD_W(IW) = INRAD_W(IW) - DLN(IOR,N) * WALL(IW)%ONE_D%BAND(IBND)%ILW(N) ! update incoming radiation,step 2
            ENDDO WALL_LOOP2
            !$OMP END DO

            !$OMP DO SCHEDULE(GUIDED)
            WALL_LOOP3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
               IF (WALL(IW)%BOUNDARY_TYPE/=OPEN_BOUNDARY)   CYCLE WALL_LOOP3
               IOR = WALL(IW)%ONE_D%IOR
               IF (DLN(IOR,N)>=0._EB) CYCLE WALL_LOOP3     ! outgoing
               IIG = WALL(IW)%ONE_D%IIG
               JJG = WALL(IW)%ONE_D%JJG
               KKG = WALL(IW)%ONE_D%KKG
               WALL(IW)%ONE_D%BAND(IBND)%ILW(ANGLE_INC_COUNTER) = WALL(IW)%ONE_D%BAND(IBND)%ILW(ANGLE_INC_COUNTER) - &
                                                                  DLN(IOR,N)*IL(IIG,JJG,KKG)
            ENDDO WALL_LOOP3
            !$OMP END DO
            !$OMP END PARALLEL

            CFACE_LOOP2: DO ICF=1,N_CFACE_CELLS
               IF (CFACE(ICF)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE CFACE_LOOP2
               DLA = (/DLX(N),DLY(N),DLZ(N)/)
               DLF = DOT_PRODUCT(CFACE(ICF)%NVEC,DLA) ! face normal * radiation angle
               IF (DLF>=0._EB) CYCLE CFACE_LOOP2      ! outgoing
               II = CFACE(ICF)%ONE_D%II
               JJ = CFACE(ICF)%ONE_D%JJ
               KK = CFACE(ICF)%ONE_D%KK
               INRAD_F(ICF) = INRAD_F(ICF) + DLF * CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N) ! update incoming radiation,step 1
               CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N) = IL(II,JJ,KK)
               INRAD_F(ICF) = INRAD_F(ICF) - DLF * CFACE(ICF)%ONE_D%BAND(IBND)%ILW(N) ! update incoming radiation,step 2
            ENDDO CFACE_LOOP2


            ! Calculate integrated intensity UIID

            IF (WIDE_BAND_MODEL) THEN
               UIID(:,:,:,IBND) = UIID(:,:,:,IBND) + WEIGH_CYL*RSA(N)*IL
            ELSE
               UIID(:,:,:,ANGLE_INC_COUNTER) = UIID(:,:,:,ANGLE_INC_COUNTER) + WEIGH_CYL*RSA(N)*IL
            ENDIF

            ! Interpolate boundary intensities onto other meshes.
            ! IL_S is an array holding the intensities IL of cells just outside of mesh NOM.

            INTERPOLATION_LOOP: DO NOM=1,NMESHES
               IF (NM==NOM) CYCLE INTERPOLATION_LOOP
               IF (EVACUATION_ONLY(NOM)) CYCLE INTERPOLATION_LOOP
               M2=>OMESH(NOM)
               IF (M2%NIC_S==0) CYCLE INTERPOLATION_LOOP
               OTHER_WALL_LOOP: DO LL=1,M2%NIC_S
                  M2%IL_S(LL,N,IBND) = IL(M2%IIO_S(LL),M2%JJO_S(LL),M2%KKO_S(LL))
               ENDDO OTHER_WALL_LOOP
            ENDDO INTERPOLATION_LOOP

            ! Compute projected intensity on particles

            IF (SOLID_PARTICLES) THEN
               PARTICLE_RADIATION_LOOP: DO IP=1,NLP
                  LP => LAGRANGIAN_PARTICLE(IP)
                  LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
                  SELECT CASE(LPC%N_ORIENTATION)
                     CASE(0)
                        CYCLE PARTICLE_RADIATION_LOOP
                     CASE(1)
                        COSINE = ORIENTATION_VECTOR(1,LP%ORIENTATION_INDEX)*DLX(N) + &
                                 ORIENTATION_VECTOR(2,LP%ORIENTATION_INDEX)*DLY(N) + &
                                 ORIENTATION_VECTOR(3,LP%ORIENTATION_INDEX)*DLZ(N)
                        IF (COSINE<0._EB) THEN
                           IF (LPC%MASSLESS_TARGET) THEN
                              LP%ONE_D%BAND(IBND)%ILW(N) = -COSINE * IL(LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)
                              IF (N==LPC%NEAREST_RAD_ANGLE_INDEX) &
                                 LP%ONE_D%IL(IBND) = IL(LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)
                           ELSE
                              ! IL_UP does not account for the absorption of radiation within the cell occupied by the particle
                              LP%ONE_D%BAND(IBND)%ILW(N) = -COSINE * IL_UP(LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)
                           ENDIF
                        ENDIF
                     CASE(2:)
                        LP%ONE_D%BAND(IBND)%ILW(N) = ORIENTATION_FACTOR(N,LP%ORIENTATION_INDEX)*RSA(N)* &
                                                     IL(LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG)
                  END SELECT
               ENDDO PARTICLE_RADIATION_LOOP
            ENDIF

            ! Save the intensities for the radiation file (RADF)

            DO NN=1,N_RADF
               RF => RAD_FILE(NN)
               DO K=RF%K1,RF%K2,RF%K_STEP
                  DO J=RF%J1,RF%J2,RF%J_STEP
                     DO I=RF%I1,RF%I2,RF%I_STEP
                        RF%IL_SAVE(I,J,K,N) = IL(I,J,K)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO

         ENDDO ANGLE_LOOP

      ENDDO UPDATE_LOOP

      ! Compute incoming flux on walls and particles

      DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         IF (WALL(IW)%BOUNDARY_TYPE/=SOLID_BOUNDARY) CYCLE
         SF  => SURFACE(WALL(IW)%SURF_INDEX)
         EFLUX = EVALUATE_RAMP(T,SF%TAU(TIME_EFLUX),SF%RAMP_INDEX(TIME_EFLUX))*SF%EXTERNAL_FLUX
         WALL(IW)%ONE_D%Q_RAD_IN  = WALL(IW)%ONE_D%Q_RAD_IN + WALL(IW)%ONE_D%EMISSIVITY*(INRAD_W(IW)+BBFA*EFLUX)
      ENDDO

      DO ICF=1,N_CFACE_CELLS
         IF (CFACE(ICF)%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE
         SF  => SURFACE(CFACE(ICF)%SURF_INDEX)
         EFLUX = EVALUATE_RAMP(T,SF%TAU(TIME_EFLUX),SF%RAMP_INDEX(TIME_EFLUX))*SF%EXTERNAL_FLUX
         CFACE(ICF)%ONE_D%Q_RAD_IN  = CFACE(ICF)%ONE_D%Q_RAD_IN + CFACE(ICF)%ONE_D%EMISSIVITY*(INRAD_F(ICF)+BBFA*EFLUX)
      ENDDO

   ENDIF INTENSITY_UPDATE

   ! Save source term for the energy equation (QR = -DIV Q)

   IF (WIDE_BAND_MODEL) THEN
      QR = QR + KAPPA_GAS*UIID(:,:,:,IBND)-KFST4_GAS
      IF (NLP>0 .AND. N_LP_ARRAY_INDICES>0) THEN
         QR_W = QR_W + KAPPA_PART*UIID(:,:,:,IBND) - KFST4_PART
      ENDIF
   ENDIF

ENDDO BAND_LOOP

! Sum up intensities and compute incoming flux at open boundaries

IF (UPDATE_INTENSITY) THEN

   UII = SUM(UIID, DIM = 4)

   DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      IF (WALL(IW)%BOUNDARY_TYPE/=OPEN_BOUNDARY) CYCLE
      WALL(IW)%ONE_D%Q_RAD_IN = 0._EB
      DO IBND=1,NUMBER_SPECTRAL_BANDS
         WALL(IW)%ONE_D%Q_RAD_IN  = WALL(IW)%ONE_D%Q_RAD_IN + SUM(WALL(IW)%ONE_D%BAND(IBND)%ILW(1:NUMBER_RADIATION_ANGLES))
      ENDDO
   ENDDO

ENDIF

! Save source term for the energy equation (QR = -DIV Q). Done only in one-band (gray gas) case.

IF (.NOT. WIDE_BAND_MODEL) THEN
   QR  = KAPPA_GAS*UII - KFST4_GAS
   IF (NLP>0 .AND. N_LP_ARRAY_INDICES>0) QR_W = QR_W + KAPPA_PART*UII - KFST4_PART
ENDIF

! Calculate the incoming radiative flux onto the solid particles

IF (SOLID_PARTICLES .AND. UPDATE_INTENSITY) THEN
   PARTICLE_LOOP: DO IP=1,NLP
      LP => LAGRANGIAN_PARTICLE(IP)
      LPC => LAGRANGIAN_PARTICLE_CLASS(LP%CLASS_INDEX)
      SF => SURFACE(LPC%SURF_INDEX)
      IF (LPC%SOLID_PARTICLE .OR. LPC%MASSLESS_TARGET) THEN
         EFLUX = EVALUATE_RAMP(T,SF%TAU(TIME_EFLUX),SF%RAMP_INDEX(TIME_EFLUX))*SF%EXTERNAL_FLUX
         IF (LP%ORIENTATION_INDEX>0) THEN
            LP%ONE_D%Q_RAD_IN = 0._EB
            DO IBND=1,NUMBER_SPECTRAL_BANDS
               LP%ONE_D%Q_RAD_IN = LP%ONE_D%Q_RAD_IN + LP%ONE_D%EMISSIVITY * &
                                 (WEIGH_CYL*SUM(LP%ONE_D%BAND(IBND)%ILW(1:NUMBER_RADIATION_ANGLES)) + EFLUX)
            ENDDO
         ELSE
            LP%ONE_D%Q_RAD_IN = LP%ONE_D%EMISSIVITY*(0.25_EB*UII(LP%ONE_D%IIG,LP%ONE_D%JJG,LP%ONE_D%KKG) + EFLUX)
         ENDIF
         IF (LPC%SOLID_PARTICLE) LP%ONE_D%Q_RAD_OUT = LP%ONE_D%EMISSIVITY*SIGMA*LP%ONE_D%TMP_F**4
      ENDIF
   ENDDO PARTICLE_LOOP
ENDIF

DEALLOCATE(IJK_SLICE)

! Write out intensities to the radiation file (RADF)

IF (N_RADF>0 .AND. T>=RADF_CLOCK(NM)) THEN

   DO N=1,N_RADF
      RF => RAD_FILE(N)
      WRITE(LU_RADF(N,NM),'(/A)') 'TIME'
      WRITE(LU_RADF(N,NM),'(F8.2)') T
      WRITE(LU_RADF(N,NM),'(/A)') 'INTENSITIES'
      WRITE(FORMT,'(A,I4,A)') '(',NUMBER_RADIATION_ANGLES+2,'F12.2)'
      DO K=RF%K1,RF%K2,RF%K_STEP
         DO J=RF%J1,RF%J2,RF%J_STEP
            DO I=RF%I1,RF%I2,RF%I_STEP
               WRITE(LU_RADF(N,NM),FORMT) TMP(I,J,K),0._EB,(RF%IL_SAVE(I,J,K,NN),NN=1,NUMBER_RADIATION_ANGLES)
            ENDDO
         ENDDO
      ENDDO
   ENDDO

   RADF_CLOCK(NM) = RADF_CLOCK(NM) + DT_RADF

ENDIF

END SUBROUTINE RADIATION_FVM


END SUBROUTINE COMPUTE_RADIATION


REAL(EB) FUNCTION BLACKBODY_FRACTION(L1,L2,TEMP)

! Calculates the fraction of black body radiation between wavelengths L1 and L2 (micron) in Temperature TEMP

USE MATH_FUNCTIONS, ONLY: INTERPOLATE1D_UNIFORM
REAL(EB),INTENT(IN) :: L1,L2,TEMP
REAL(EB) :: LT1,LT2,BBFLOW,BBFHIGH

LT1    =   L1 * TEMP/LTSTEP
CALL INTERPOLATE1D_UNIFORM(LBOUND(BBFRAC,1),BBFRAC,LT1,BBFLOW)

LT2    =   L2 * TEMP/LTSTEP
CALL INTERPOLATE1D_UNIFORM(LBOUND(BBFRAC,1),BBFRAC,LT2,BBFHIGH)

BLACKBODY_FRACTION = BBFHIGH - BBFLOW

END FUNCTION BLACKBODY_FRACTION


FUNCTION GET_KAPPA(Z_IN,TYY,IBND)

! Returns the radiative absorption

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION_ALL,GET_MOLECULAR_WEIGHT
REAL(EB), INTENT(INOUT) :: Z_IN(1:N_TRACKED_SPECIES)
REAL(EB) :: KAPPA_TEMP,INT_FAC,GET_KAPPA,SCALED_Y_RADCAL_SPECIES,MWA
INTEGER, INTENT(IN) :: IBND,TYY
INTEGER :: LBND,UBND,N

GET_KAPPA = 0._EB

CALL GET_MOLECULAR_WEIGHT(Z_IN,MWA)

DO N = 1, N_RADCAL_ARRAY_SIZE
   SCALED_Y_RADCAL_SPECIES = DOT_PRODUCT(Z2RADCAL_SPECIES(N,:),Z_IN)
   IF (SCALED_Y_RADCAL_SPECIES<TWO_EPSILON_EB) CYCLE
   IF (RADCAL_SPECIES_INDEX(N)==16) THEN
      INT_FAC = MAX(0._EB,SCALED_Y_RADCAL_SPECIES)**0.25_EB
   ELSE
      INT_FAC = MAX(0._EB,SCALED_Y_RADCAL_SPECIES*MWA)**0.25_EB
   ENDIF
   LBND = INT(INT_FAC)
   INT_FAC = INT_FAC - LBND
   LBND = MIN(LBND,N_KAPPA_Y)
   UBND = MIN(LBND+1,N_KAPPA_Y)
   KAPPA_TEMP = RADCAL_SPECIES2KAPPA(N,LBND,TYY,IBND)
   GET_KAPPA = GET_KAPPA + KAPPA_TEMP + INT_FAC*(RADCAL_SPECIES2KAPPA(N,UBND,TYY,IBND)-KAPPA_TEMP)
ENDDO

END FUNCTION GET_KAPPA

END MODULE RAD
