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

PUBLIC INIT_RADIATION,COMPUTE_RADIATION,GET_KAPPA,GET_REV_radi
 
CONTAINS
 
 
SUBROUTINE INIT_RADIATION

USE MEMORY_FUNCTIONS, ONLY : CHKMEMERR
USE MIEV
USE RADCALV
USE DEVICE_VARIABLES, ONLY : DEVICE, GAS_CELL_RAD_FLUX, GAS_CELL_RAD_DEVC_INDEX, N_GAS_CELL_RAD_DEVC
REAL(EB) :: THETAUP,THETALOW,PHIUP,PHILOW,F_THETA,PLANCK_C2,KSI,LT,RCRHO,YY,BBF,AP0,AMEAN
INTEGER  :: N,I,J,K,IZERO,NN,NI,II,JJ,IIM,JJM,IBND,NS,NRA,NSB
 
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
         DLX(N) =  (SIN(PHIUP)-SIN(PHILOW)) *F_THETA
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
 
RTMPMAX = 2400._EB     ! Maximum temperature for property tables
RTMPMIN = 300._EB      ! Minimum temperature for property tables
  
! Setup spectral information
 
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
!          SVF=SOOT VOLUME FRACTION OF J TH ELEMENT
!          OMMIN=MINIMUM WAVE NUMBER IN SPECTRUM, CM-1.
!          OMMAX=MAXIMUM WAVE NUMBER IN SPECTRUM, CM-1.
!
!-------------------------------------------------------------------------

MAKE_KAPPA_ARRAYS: IF (SOOT_INDEX /= 0 .OR. CO_INDEX /= 0 .OR. FUEL_INDEX /= 0 .OR. CO2_INDEX /= 0 .OR. H2O_INDEX /=0) THEN 

   KAPPA_ARRAY = .TRUE.

   ! Allocate arrays for RadCal

   CALL RCALLOC
 
   ! Set the Mean Beam Length to 5 times the smallest cell dimension unless the user desires otherwise
 
   IF (PATH_LENGTH < 0._EB) PATH_LENGTH = MIN( 10._EB , 5._EB*CHARACTERISTIC_CELL_SIZE )
   DD = MAX(PATH_LENGTH,1.0E-4_EB)
 
   ! Using RadCal, create look-up tables for the absorption coefficients for all gas species, mixture fraction or aerosols

   N_KAPPA_ARRAY=0
      IF (FUEL_INDEX /= 0) N_KAPPA_ARRAY = N_KAPPA_ARRAY + 1
      IF (CO2_INDEX /= 0)  N_KAPPA_ARRAY = N_KAPPA_ARRAY + 1
      IF (CO_INDEX /= 0)   N_KAPPA_ARRAY = N_KAPPA_ARRAY + 1
      IF (H2O_INDEX /= 0)  N_KAPPA_ARRAY = N_KAPPA_ARRAY + 1
      IF (SOOT_INDEX /= 0) N_KAPPA_ARRAY = N_KAPPA_ARRAY + 1
   Z2KAPPA_T = 42
   Z2KAPPA_M = 50
   NS = 0
   ALLOCATE (Z2KAPPA_M4(N_KAPPA_ARRAY),STAT=IZERO)
   CALL ChkMemErr('RADI','Z2KAPPA_M4',IZERO)
   ALLOCATE (KAPPA_INDEX(N_KAPPA_ARRAY),STAT=IZERO)
   CALL ChkMemErr('RADI','KAPPA_INDEX',IZERO)
   IF (FUEL_INDEX > 0) THEN
      NS = NS + 1
      Z2KAPPA_M4(NS) = REAL(Z2KAPPA_M,EB)**4/SPECIES(FUEL_INDEX)%MW
      KAPPA_INDEX(NS) = FUEL_INDEX
   ENDIF
   IF (CO2_INDEX > 0) THEN
      NS = NS + 1
      Z2KAPPA_M4(NS) = REAL(Z2KAPPA_M,EB)**4/SPECIES(CO2_INDEX)%MW
      KAPPA_INDEX(NS) = CO2_INDEX
   ENDIF
   IF (CO_INDEX > 0) THEN
      NS = NS + 1
      Z2KAPPA_M4(NS) = REAL(Z2KAPPA_M,EB)**4/SPECIES(CO_INDEX)%MW
      KAPPA_INDEX(NS) = CO_INDEX
   ENDIF      
   IF (H2O_INDEX > 0) THEN
      NS = NS + 1
      Z2KAPPA_M4(NS) = REAL(Z2KAPPA_M,EB)**4/SPECIES(H2O_INDEX)%MW
      KAPPA_INDEX(NS) = H2O_INDEX
   ENDIF      
   IF (SOOT_INDEX > 0) THEN
      NS = NS + 1
      Z2KAPPA_M4(NS) = REAL(Z2KAPPA_M,EB)**4*5._EB
      KAPPA_INDEX(NS) = SOOT_INDEX
   ENDIF
   ALLOCATE (Z2KAPPA(N_KAPPA_ARRAY,0:Z2KAPPA_M,0:Z2KAPPA_T,NSB),STAT=IZERO)
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
         RCT = RTMPMIN + K*(RTMPMAX-RTMPMIN)/Z2KAPPA_T         
         IF (NSB>1) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),RCT)
         Y_LOOP_Z: DO J=0,Z2KAPPA_M
            YY = (REAL(J,EB)/REAL(Z2KAPPA_M,EB))**4
            N = 0
            KAPPA_SPECIES: DO NS = 1, 5
               SELECT CASE(NS)
                  CASE(1) ! FUEL
                     IF(FUEL_INDEX > 0) THEN
                        N = N + 1
                        SVF    = 0._EB
                        SPECIE = 0._EB
                        P      = 0._EB
                        SPECIE(3) = YY
                        P(3)      = YY
                        P(6)      = (1._EB-YY)
                        CALL RADCAL(AMEAN,AP0)
                        IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
                           Z2KAPPA(N,J,K,IBND) = MIN(AMEAN,AP0)
                        ELSE
                           Z2KAPPA(N,J,K,IBND) = AP0/BBF
                        ENDIF
                     ENDIF
                  CASE(2) ! CO2
                     IF(CO2_INDEX > 0) THEN
                        N = N + 1
                        SVF    = 0._EB
                        SPECIE = 0._EB
                        P      = 0._EB
                        SPECIE(1) = YY
                        P(1)      = YY
                        P(6)      = (1._EB-YY)
                        CALL RADCAL(AMEAN,AP0)
                        IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
                           Z2KAPPA(N,J,K,IBND) = MIN(AMEAN,AP0)
                        ELSE
                           Z2KAPPA(N,J,K,IBND) = AP0/BBF
                        ENDIF
                     ENDIF
                  CASE(3) ! CO
                     IF(CO_INDEX > 0) THEN
                        N = N + 1
                        SVF    = 0._EB
                        SPECIE = 0._EB
                        P      = 0._EB
                        SPECIE(4) = YY
                        P(4)      = YY
                        P(6)      = (1._EB-YY)
                        CALL RADCAL(AMEAN,AP0)
                        IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
                           Z2KAPPA(N,J,K,IBND) = MIN(AMEAN,AP0)
                        ELSE
                           Z2KAPPA(N,J,K,IBND) = AP0/BBF
                        ENDIF
                     ENDIF
                  CASE(4) ! H2O
                     IF(H2O_INDEX > 0) THEN
                        N = N + 1
                        SVF    = 0._EB
                        SPECIE = 0._EB
                        P      = 0._EB
                        SPECIE(2) = YY
                        P(2)      = YY
                        P(6)      = (1._EB-YY)
                        CALL RADCAL(AMEAN,AP0)
                        IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
                           Z2KAPPA(N,J,K,IBND) = MIN(AMEAN,AP0)
                        ELSE
                           Z2KAPPA(N,J,K,IBND) = AP0/BBF
                        ENDIF
                     ENDIF
                  CASE(5) !Soot
                     IF(SOOT_INDEX > 0) THEN
                        N = N + 1
                        RCRHO = MW_AIR*P_INF/(R0*RCT)
                        YY = YY * 0.2_EB
                        SPECIE = 0._EB
                        P      = 0._EB
                        SPECIE(5) = YY*RCRHO/RHO_SOOT
                        P(6) = 1._EB
                        SVF = YY*RCRHO/RHO_SOOT
                        CALL RADCAL(AMEAN,AP0)                        
                        IF (NSB==1 .AND. PATH_LENGTH > 0.0_EB) THEN
                           Z2KAPPA(N,J,K,IBND) = MIN(AMEAN,AP0)
                        ELSE
                           Z2KAPPA(N,J,K,IBND) = AP0/BBF
                        ENDIF
                     ENDIF
               END SELECT
            END DO KAPPA_SPECIES
         ENDDO Y_LOOP_z
      ENDDO T_LOOP_Z
   ENDDO BAND_LOOP_Z

   CALL RCDEALLOC  ! Deallocate RadCal arrays

ENDIF MAKE_KAPPA_ARRAYS
 
!-----------------------------------------------------
!
!     Tables for droplet absorption coefficients
! 
!-----------------------------------------------------

DROPLETS: IF (N_EVAP_INDICES>0) THEN
   GET_PC_RADI: DO J=1,N_PART
      IF (.NOT.PARTICLE_CLASS(J)%MASSLESS) CALL MEAN_CROSS_SECTIONS(J)
   ENDDO GET_PC_RADI
ENDIF DROPLETS
 
! A few miscellaneous constants

FOUR_SIGMA = 4._EB*SIGMA
RPI_SIGMA  = RPI*SIGMA
 

END SUBROUTINE INIT_RADIATION


 
SUBROUTINE COMPUTE_RADIATION(T,NM)

! Call radiation routine or simply specify the radiative loss

USE MESH_POINTERS
USE COMP_FUNCTIONS, ONLY : SECOND  
REAL(EB) :: TNOW,T
INTEGER, INTENT(IN) :: NM 

IF (EVACUATION_ONLY(NM)) RETURN

TNOW=SECOND()
 
CALL POINT_TO_MESH(NM)

IF (RADIATION) THEN
   CALL RADIATION_FVM(T,NM)
ELSE
   IF (N_REACTIONS>0) QR = -RADIATIVE_FRACTION*Q
ENDIF

TUSED(9,NM)=TUSED(9,NM)+SECOND()-TNOW

CONTAINS 
 
SUBROUTINE RADIATION_FVM(T,NM)
USE MIEV
USE MATH_FUNCTIONS, ONLY : INTERPOLATE1D, EVALUATE_RAMP 
USE DEVICE_VARIABLES, ONLY : DEVICE_TYPE,DEVICE, GAS_CELL_RAD_FLUX, GAS_CELL_RAD_DEVC_INDEX, N_GAS_CELL_RAD_DEVC
REAL(EB) :: T, RAP, AX, AXU, AXD, AY, AYU, AYD, AZ, VC, RU, RD, RP, &
            ILXU, ILYU, ILZU, QVAL, BBF, BBFA, NCSDROP, RSA_RAT,COSINE,EFLUX,TYY_FAC, &
            AIU_SUM, A_SUM,RSUM,HRR,RTE_SOURCE_CORRECTION_FACTOR,VOL
INTEGER  :: N, NN,IIG,JJG,KKG,I,J,K,IW,II,JJ,KK,IOR,IC,IWUP,IWDOWN, &
            ISTART, IEND, ISTEP, JSTART, JEND, JSTEP, &
            KSTART, KEND, KSTEP, NSTART, NEND, NSTEP, &
            I_UIID, N_UPDATES, IBND, TYY, NOM, IBC,EVAP_INDEX,NRA, I_DROP
REAL(EB) :: XID,YJD,ZKD,ZZ_GET(0:N_TRACKED_SPECIES),KAPPA_PART,SURFACE_AREA
INTEGER :: IPC,IID,JJD,KKD,ID
LOGICAL :: UPDATE_INTENSITY
REAL(EB), POINTER, DIMENSION(:,:,:) :: KFST4=>NULL(), IL=>NULL(), UIIOLD=>NULL(), KAPPAW=>NULL(), &
                                       KFST4W=>NULL(), EXTCOE=>NULL(), SCAEFF=>NULL(),IL_UP=>NULL(),UIIP=>NULL()
REAL(EB), POINTER, DIMENSION(:)     :: OUTRAD_W=>NULL(), INRAD_W=>NULL()
INTEGER, INTENT(IN) :: NM
TYPE (OMESH_TYPE), POINTER :: M2=>NULL()
TYPE(SURFACE_TYPE), POINTER :: SF=>NULL()
TYPE(DEVICE_TYPE), POINTER :: DV=>NULL()
TYPE(PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()
TYPE(DROPLET_TYPE), POINTER :: DR=>NULL()

KFST4    => WORK1
IL       => WORK2
UIIOLD   => WORK3
EXTCOE   => WORK4
KAPPAW   => WORK5
SCAEFF   => WORK6
KFST4W   => WORK7
IL_UP    => WORK8
OUTRAD_W => WALL_WORK1
INRAD_W  => WALL_WORK2
 
! Ratio of solid angle, used in scattering
 
NRA     = NUMBER_RADIATION_ANGLES
RSA_RAT = 1._EB/(1._EB-1._EB/NRA)
 
! Check if it time to update radiation intensity field
 
RAD_CALL_COUNTER  = RAD_CALL_COUNTER+1
IF ( MOD(RAD_CALL_COUNTER,TIME_STEP_INCREMENT)==1 .OR. TIME_STEP_INCREMENT==1) THEN
   UPDATE_INTENSITY = .TRUE.
ELSE
   UPDATE_INTENSITY = .FALSE.
ENDIF

IF (WIDE_BAND_MODEL) THEN
   QR = 0._EB
ENDIF
IF (UPDATE_INTENSITY) THEN
   QRADIN = 0._EB
ENDIF
 
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
 
   IF_DROPLETS_INCLUDED: IF (NLP>0 .AND. N_EVAP_INDICES>0) THEN
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
                     IF (PC%TREE .OR. PC%SURF_INDEX>0) CYCLE PC_LOOP
                     EVAP_INDEX = PC%EVAP_INDEX
                     IF (EVAP_INDEX==0) CYCLE PC_LOOP
                     IF (ABS(AVG_DROP_AREA(I,J,K,EVAP_INDEX))<ZERO_P) CYCLE PC_LOOP
                     NCSDROP = AVG_DROP_AREA(I,J,K,EVAP_INDEX)
                     ! Absorption and scattering efficiency
                     CALL INTERPOLATE1D(PC%R50,PC%WQABS(:,IBND),AVG_DROP_RAD(I,J,K,EVAP_INDEX),QVAL) 
                     KAPPAW(I,J,K) = KAPPAW(I,J,K) + NCSDROP*QVAL
                     CALL INTERPOLATE1D(PC%R50,PC%WQSCA(:,IBND),AVG_DROP_RAD(I,J,K,EVAP_INDEX),QVAL)
                     SCAEFF(I,J,K) = SCAEFF(I,J,K) + NCSDROP*QVAL
                     KFST4W(I,J,K) = KFST4W(I,J,K)+ BBF*KAPPAW(I,J,K)*FOUR_SIGMA*AVG_DROP_TMP(I,J,K,EVAP_INDEX)**4
                  ENDDO PC_LOOP
            ENDDO ZLOOPM
         ENDDO
      ENDDO

      QR_W = 0._EB

   ENDIF IF_DROPLETS_INCLUDED

   ! Compute raised vegetation absorption and emission source term 

   IF_VEG_INCLUDED: IF (NLP>0 .AND. TREE_MESH(NM)) THEN
      VEG_ELEM_LOOP: DO ID = 1,NLP
         IPC = DROPLET(ID)%CLASS
         IF (.NOT. PARTICLE_CLASS(IPC)%TREE) CYCLE VEG_ELEM_LOOP
         XID = CELLSI(FLOOR((DROPLET(ID)%X-XS)*RDXINT))
         YJD = CELLSJ(FLOOR((DROPLET(ID)%Y-YS)*RDYINT))
         ZKD = CELLSK(FLOOR((DROPLET(ID)%Z-ZS)*RDZINT))
         IID = FLOOR(XID+1._EB)
         JJD = FLOOR(YJD+1._EB)
         KKD = FLOOR(ZKD+1._EB)
         KAPPAW(IID,JJD,KKD) = KAPPAW(IID,JJD,KKD) + DROPLET(ID)%VEG_KAPPA
         KFST4W(IID,JJD,KKD) = KFST4W(IID,JJD,KKD) + DROPLET(ID)%VEG_EMISS
      ENDDO VEG_ELEM_LOOP
      QR_W = 0._EB
   ENDIF IF_VEG_INCLUDED

   ! Virtual particles
   
   IF (NLP>0 .AND. VIRTUAL_PARTICLES) THEN
      DO ID = 1,NLP
         DR => DROPLET(ID)
         IPC = DR%CLASS
         PC => PARTICLE_CLASS(IPC)
         IF (PC%SURF_INDEX<1) CYCLE
         SF => SURFACE(PC%SURF_INDEX)
         XID = CELLSI(FLOOR((DR%X-XS)*RDXINT))
         YJD = CELLSJ(FLOOR((DR%Y-YS)*RDYINT))
         ZKD = CELLSK(FLOOR((DR%Z-ZS)*RDZINT))
         IID = FLOOR(XID+1._EB)
         JJD = FLOOR(YJD+1._EB)
         KKD = FLOOR(ZKD+1._EB)
         SELECT CASE(SF%GEOMETRY)
            CASE(SURF_CARTESIAN)
               SURFACE_AREA = SF%LENGTH*SF%WIDTH  
            CASE(SURF_CYLINDRICAL)
               SURFACE_AREA = SF%LENGTH*TWOPI*DR%R
            CASE(SURF_SPHERICAL)
               SURFACE_AREA = 4._EB*PI*DR%R**2
         END SELECT
         KAPPA_PART = 0.25_EB*DR%PWT*SURFACE_AREA*RDX(IID)*RRN(IID)*RDY(JJD)*RDZ(KKD)
         KAPPAW(IID,JJD,KKD) = KAPPAW(IID,JJD,KKD) + KAPPA_PART
         KFST4W(IID,JJD,KKD) = KFST4W(IID,JJD,KKD) + FOUR_SIGMA*KAPPA_PART*DR%TMP**4
      ENDDO
      QR_W = 0._EB
   ENDIF
 
   ! Compute absorption coefficient KAPPA
 
   KAPPA = KAPPA0

   IF (KAPPA_ARRAY) THEN
      TYY_FAC = Z2KAPPA_T / (RTMPMAX-RTMPMIN)
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               TYY = MAX(0 , MIN(Z2KAPPA_T,INT((TMP(I,J,K) - RTMPMIN) * TYY_FAC)))
               ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(I,J,K,1:N_TRACKED_SPECIES)
               KAPPA(I,J,K) = KAPPA(I,J,K) + GET_KAPPA(ZZ_GET,TYY,IBND)
            ENDDO
         ENDDO
      ENDDO
   ENDIF

   ! Compute source term KAPPA*4*SIGMA*TMP**4

   HRR  = 0._EB
   RSUM = 0._EB
   BBF  = 1._EB

   DO K=1,KBAR
      DO J=1,JBAR
         DO I=1,IBAR
            IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
            IF (WIDE_BAND_MODEL) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMP(I,J,K))
            KFST4(I,J,K) = BBF*KAPPA(I,J,K)*FOUR_SIGMA*TMP(I,J,K)**4
            IF (RADIATIVE_FRACTION>0._EB .AND. Q(I,J,K)>0._EB) THEN
               VOL  = R(I)*DX(I)*DY(J)*DZ(K)
               HRR  = HRR  + Q(I,J,K)*VOL
               RSUM = RSUM + KFST4(I,J,K)*VOL
            ENDIF
         ENDDO
      ENDDO
   ENDDO

   ! Correct the source term in the RTE based on user-specified RADIATIVE_FRACTION

   IF (RSUM>0._EB) THEN
      RTE_SOURCE_CORRECTION_FACTOR = RADIATIVE_FRACTION*HRR/RSUM
      IF (WIDE_BAND_MODEL) THEN
         UIIP => UIID(:,:,:,IBND)
      ELSE
         UIIP => UII
      ENDIF
      DO K=1,KBAR
         DO J=1,JBAR
            DO I=1,IBAR
               IF (SOLID(CELL_INDEX(I,J,K))) CYCLE
               IF (RADIATIVE_FRACTION>0._EB .AND. Q(I,J,K)>0._EB) THEN
                  KFST4(I,J,K) = KFST4(I,J,K)*RTE_SOURCE_CORRECTION_FACTOR + KAPPA(I,J,K)*UIIP(I,J,K)
               ENDIF
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
      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS+N_VIRTUAL_WALL_CELLS
         IF (BOUNDARY_TYPE(IW) == OPEN_BOUNDARY) THEN
            BBF = BBFA
         ELSE
            IF (WIDE_BAND_MODEL) BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMP_F(IW))
            IBC = IJKW(5,IW)
            SF  => SURFACE(IBC)
            IF (.NOT. SF%INTERNAL_RADIATION) QRADOUT(IW) = E_WALL(IW)*SIGMA*TMP_F(IW)**4
         ENDIF
         OUTRAD_W(IW) = BBF*RPI*QRADOUT(IW)
      ENDDO

      ! Compute boundary condition term incoming radiation integral
 
      DO IW = 1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) CYCLE
         IOR = IJKW(4,IW)
         INRAD_W(IW) = SUM(-DLN(IOR,:)* WALL(IW)%ILW(:,IBND),1, DLN(IOR,:)<0._EB)
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

         DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
            IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY) WALL(IW)%ILW(ANGLE_INC_COUNTER,IBND) = 0._EB
         ENDDO
 
         ! Set the bounds and increment for the angleloop. Step downdard because in cylindrical case the Nth angle 
         ! boundary condition comes from (N+1)th angle.
          
         NSTART    = NRA - ANGLE_INC_COUNTER + 1
         NEND      = 1
         NSTEP     = -ANGLE_INCREMENT

         IL(:,:,:) = BBFA*RPI_SIGMA*TMPA4

         ANGLE_LOOP: DO N = NSTART,NEND,NSTEP  ! Sweep through control angles
 
            ! Boundary conditions: Intensities leaving the boundaries.
            
            WALL_LOOP1: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
               IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY) CYCLE WALL_LOOP1
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
                     AXU = 2._EB*SIN(DPHI0/2.)*RU       * DZ(K) * ABS(DLX(N))
                     AXD = 2._EB*SIN(DPHI0/2.)*RD       * DZ(K) * ABS(DLX(N))
                     AYU = DX(I)             * DZ(K) * ABS(DLB(N))
                     AYD = DX(I)             * DZ(K) * ABS(DLY(N))
                     AZ  = DX(I)  * RP*DPHI0         * ABS(DLZ(N))
                     IF (MODULO(N,NRP(1))==1) AYD = 0._EB  ! Zero out the terms involving symmetric overhang
                     IF (MODULO(N,NRP(1))==0) AYU = 0._EB
                     IF (IC/=0) THEN
                        IW = WALL_INDEX(IC,-ISTEP)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILXU = WALL(IW)%ILW(N,IBND)
                        ENDIF
                        IW = WALL_INDEX(IC,-JSTEP*2)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILYU = WALL(IW)%ILW(N,IBND)
                        ENDIF
                        IW = WALL_INDEX(IC,-KSTEP*3)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILZU = WALL(IW)%ILW(N,IBND)
                        ENDIF
                     ENDIF
                     AIU_SUM = AXU*ILXU + AYU*ILYU + AZ*ILZU
                     A_SUM = AXD + AYD + AZ
                     RAP = 1._EB/(A_SUM + EXTCOE(I,J,K)*VC*RSA(N))
                     IL(I,J,K) = MAX(0._EB, RAP * (AIU_SUM + VC*RSA(N)*RFPI* &
                                 ( KFST4(I,J,K)+KFST4W(I,J,K) +RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
                     IF (VIRTUAL_PARTICLES) IL_UP(I,J,K) = MAX(0._EB,AIU_SUM/A_SUM)
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
                        ENDIF
                        IW = WALL_INDEX(IC,-KSTEP*3)
                        IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                           ILZU = WALL(IW)%ILW(N,IBND)
                        ENDIF
                     ENDIF
                     AIU_SUM = AX*ILXU + AZ*ILZU 
                     A_SUM = AX + AZ
                     RAP = 1._EB/(A_SUM + EXTCOE(I,J,K)*VC*RSA(N))
                     IL(I,J,K) = MAX(0._EB, RAP * (AIU_SUM + VC*RSA(N)*RFPI* &
                                    (KFST4(I,J,K)+KFST4W(I,J,K) +  RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) ) 
                     IF (VIRTUAL_PARTICLES) IL_UP(I,J,K) = MAX(0._EB,AIU_SUM/A_SUM)
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
                              ILXU = WALL(IW)%ILW(N,IBND)
                           ENDIF
                           IW = WALL_INDEX(IC,-JSTEP*2)
                           IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                              ILYU = WALL(IW)%ILW(N,IBND)
                           ENDIF
                           IW = WALL_INDEX(IC,-KSTEP*3)
                           IF (BOUNDARY_TYPE(IW)==SOLID_BOUNDARY) THEN
                              ILZU = WALL(IW)%ILW(N,IBND)
                           ENDIF
                        ENDIF
                        A_SUM = AX + AY + AZ
                        AIU_SUM = AX*ILXU + AY*ILYU + AZ*ILZU 
                        RAP = 1._EB/(A_SUM + EXTCOE(I,J,K)*VC*RSA(N))
                        IL(I,J,K) = MAX(0._EB, RAP * (AIU_SUM + VC*RSA(N)*RFPI* &
                                       ( KFST4(I,J,K)+KFST4W(I,J,K) + RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
                        IF (VIRTUAL_PARTICLES) IL_UP(I,J,K) = MAX(0._EB,AIU_SUM/A_SUM)
                     ENDDO ILOOP
                  ENDDO JLOOP
               ENDDO KLOOP
 
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
                  IF (MODULO(N,NRP(1))==1) THEN
                     WALL(IWUP)%ILW(N,IBND) = WALL(IWDOWN)%ILW(N,IBND)
                  ELSE
                     WALL(IWUP)%ILW(N-1,IBND) = WALL(IWDOWN)%ILW(N,IBND)
                  ENDIF
               ENDDO CILOOP2
               ENDDO CKLOOP2
            ENDIF
 
            ! Boundary values: Incoming radiation
            
            WALL_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
               IF (BOUNDARY_TYPE(IW)==NULL_BOUNDARY)   CYCLE WALL_LOOP2     
               IF (BOUNDARY_TYPE(IW)==OPEN_BOUNDARY)   CYCLE WALL_LOOP2  
               IOR = IJKW(4,IW)
               IF (TWO_D .AND. .NOT.CYLINDRICAL  .AND. ABS(IOR)==2) CYCLE WALL_LOOP2  ! 2-D non cylindrical
               IF (DLN(IOR,N)>=0._EB) CYCLE WALL_LOOP2     ! outgoing
               IIG = IJKW(6,IW)
               JJG = IJKW(7,IW)
               KKG = IJKW(8,IW)
               INRAD_W(IW) = INRAD_W(IW) + DLN(IOR,N) * WALL(IW)%ILW(N,IBND) ! update incoming radiation,step 1
               WALL(IW)%ILW(N,IBND) = IL(IIG,JJG,KKG)
               INRAD_W(IW) = INRAD_W(IW) - DLN(IOR,N) * WALL(IW)%ILW(N,IBND) ! update incoming radiation,step 2
            ENDDO WALL_LOOP2
 
            WALL_LOOP3: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
               IF (BOUNDARY_TYPE(IW)/=OPEN_BOUNDARY)   CYCLE WALL_LOOP3 
               IOR = IJKW(4,IW)
               IF (DLN(IOR,N)>=0._EB) CYCLE WALL_LOOP3     ! outgoing
               IIG = IJKW(6,IW)
               JJG = IJKW(7,IW)
               KKG = IJKW(8,IW)
               WALL(IW)%ILW(ANGLE_INC_COUNTER,IBND) = WALL(IW)%ILW(ANGLE_INC_COUNTER,IBND)-DLN(IOR,N)*IL(IIG,JJG,KKG)
            ENDDO WALL_LOOP3

            ! Calculate integrated intensity UIID
 
            IF (WIDE_BAND_MODEL) THEN
               UIID(:,:,:,IBND) = UIID(:,:,:,IBND) + WEIGH_CYL*RSA(N)*IL
            ELSE
               UIID(:,:,:,ANGLE_INC_COUNTER) = UIID(:,:,:,ANGLE_INC_COUNTER) + WEIGH_CYL*RSA(N)*IL
            ENDIF
 
            ! Interpolate boundary intensities onto other meshes
 
            INTERPOLATION_LOOP: DO NOM=1,NMESHES
               IF (NM==NOM) CYCLE INTERPOLATION_LOOP
               IF (EVACUATION_ONLY(NOM)) CYCLE INTERPOLATION_LOOP
               M2=>OMESH(NOM)
               IF (M2%NIC_R==0) CYCLE INTERPOLATION_LOOP
               OTHER_WALL_LOOP: DO IW=1,MESHES(NOM)%N_EXTERNAL_WALL_CELLS
                  IF (M2%IJKW(9,IW)/=NM .OR. M2%BOUNDARY_TYPE(IW)/=INTERPOLATED_BOUNDARY) CYCLE OTHER_WALL_LOOP
                  IOR = M2%IJKW(4,IW)
                  IF (DLN(IOR,N)<=0._EB) CYCLE OTHER_WALL_LOOP
                  M2%WALL(IW)%ILW(N,IBND)=IL(M2%IJKW(10,IW),M2%IJKW(11,IW),M2%IJKW(12,IW))
               ENDDO OTHER_WALL_LOOP
            ENDDO INTERPOLATION_LOOP

            ! Compute projected intensity on gas phase devices measuring radiative flux

            IF (GAS_CELL_RAD_FLUX) THEN
               DEVC_LOOP: DO NN = 1, N_GAS_CELL_RAD_DEVC
                  DV => DEVICE(GAS_CELL_RAD_DEVC_INDEX(NN))
                  IF (DV%MESH /= NM) CYCLE DEVC_LOOP
                  COSINE = DV%ORIENTATION(1)*DLX(N)+DV%ORIENTATION(2)*DLY(N)+DV%ORIENTATION(3)*DLZ(N)
                  IF (COSINE <0._EB) DV%ILW(N,IBND) = -COSINE * IL(DV%I,DV%J,DV%K)
               ENDDO DEVC_LOOP
            ENDIF

            IF (VIRTUAL_PARTICLES) THEN
               DROPLET_RADIATION_LOOP: DO I_DROP=1,NLP
                  DR => DROPLET(I_DROP)
                  PC => PARTICLE_CLASS(DR%CLASS)
                  IF (PC%N_SPLIT < 2) CYCLE DROPLET_RADIATION_LOOP
                  IBC = PC%SURF_INDEX
                  IF (IBC<1) CYCLE DROPLET_RADIATION_LOOP
                  IW = DR%WALL_INDEX
                  COSINE = PC%ORIENTATION(DR%SPLIT_IOR,1)*DLX(N) + &
                           PC%ORIENTATION(DR%SPLIT_IOR,2)*DLY(N) + &
                           PC%ORIENTATION(DR%SPLIT_IOR,3)*DLZ(N)
                  IF (COSINE <0._EB) WALL(IW)%ILW(N,IBND) = -COSINE * IL_UP(IJKW(1,IW),IJKW(2,IW),IJKW(3,IW))
               ENDDO DROPLET_RADIATION_LOOP
            ENDIF

         ENDDO ANGLE_LOOP

      ENDDO UPDATE_LOOP
 
      ! Compute incoming flux on walls 

      DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
         IF (BOUNDARY_TYPE(IW)/=SOLID_BOUNDARY) CYCLE 
         IBC = IJKW(5,IW)
         SF  => SURFACE(IBC)      
         EFLUX = EVALUATE_RAMP(T,SF%TAU(TIME_EFLUX),SF%RAMP_INDEX(TIME_EFLUX))*SF%EXTERNAL_FLUX
         QRADIN(IW)  = QRADIN(IW) + E_WALL(IW)*(INRAD_W(IW)+BBFA*EFLUX)
      ENDDO 

   ENDIF INTENSITY_UPDATE
 
   ! Save source term for the energy equation (QR = -DIV Q)

   IF (WIDE_BAND_MODEL) THEN
      QR = QR + KAPPA*UIID(:,:,:,IBND)-KFST4
      IF (NLP>0 .AND. N_EVAP_INDICES>0) THEN
         QR_W = QR_W + KAPPAW*UIID(:,:,:,IBND) - KFST4W
      ENDIF
   ENDIF

ENDDO BAND_LOOP

! Sum up intensities and compute incoming flux at open boundaries

IF (UPDATE_INTENSITY) THEN

   UII = SUM(UIID, DIM = 4)

   DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
      IF (BOUNDARY_TYPE(IW)/=OPEN_BOUNDARY) CYCLE 
      QRADIN(IW)  = SUM(WALL(IW)%ILW)
   ENDDO 

ENDIF

! Save source term for the energy equation (QR = -DIV Q). Done only in one-band (gray gas) case.

IF (.NOT. WIDE_BAND_MODEL) THEN
   QR  = KAPPA*UII - KFST4
   IF (NLP>0 .AND. N_EVAP_INDICES>0) QR_W = QR_W + KAPPAW*UII - KFST4W

   IF_VEG_INCLUDED2: IF (NLP>0 .AND. TREE_MESH(NM)) THEN
      VEG_ELEM_LOOP2: DO ID = 1,NLP
         IPC = DROPLET(ID)%CLASS
         IF (.NOT. PARTICLE_CLASS(IPC)%TREE) CYCLE VEG_ELEM_LOOP2
         DROPLET(ID)%VEG_DIVQR = 0.0_EB
         XID = CELLSI(FLOOR((DROPLET(ID)%X-XS)*RDXINT))
         YJD = CELLSJ(FLOOR((DROPLET(ID)%Y-YS)*RDYINT))
         ZKD = CELLSK(FLOOR((DROPLET(ID)%Z-ZS)*RDZINT))
         IID = FLOOR(XID+1._EB)
         JJD = FLOOR(YJD+1._EB)
         KKD = FLOOR(ZKD+1._EB)
         DROPLET(ID)%VEG_DIVQR = DROPLET(ID)%VEG_KAPPA*UII(IID,JJD,KKD) - DROPLET(ID)%VEG_EMISS
      ENDDO VEG_ELEM_LOOP2
   ENDIF IF_VEG_INCLUDED2

ENDIF

! Check for VIRTUAL devices and particles

IF (GAS_CELL_RAD_FLUX) THEN
   DEVC_LOOP2: DO NN = 1, N_GAS_CELL_RAD_DEVC
      DV => DEVICE(GAS_CELL_RAD_DEVC_INDEX(NN))
      IF (DV%MESH /= NM) CYCLE DEVC_LOOP2
      IW = DV%VIRTUAL_WALL_INDEX
      IF (IW>0) THEN
         QRADIN(IW) = E_WALL(IW)*WEIGH_CYL*SUM(DV%ILW)
      ENDIF
   ENDDO DEVC_LOOP2
ENDIF

IF (VIRTUAL_PARTICLES) THEN
   DROPLET_LOOP: DO I=1,NLP
      DR => DROPLET(I)
      PC => PARTICLE_CLASS(DR%CLASS)
      IF (PC%SURF_INDEX>0) THEN
         IW = DR%WALL_INDEX
         SF => SURFACE(PC%SURF_INDEX)
         EFLUX = EVALUATE_RAMP(T,SF%TAU(TIME_EFLUX),SF%RAMP_INDEX(TIME_EFLUX))*SF%EXTERNAL_FLUX
         IF(PC%N_SPLIT>1) THEN
            QRADIN(IW) = E_WALL(IW)*(WEIGH_CYL*SUM(WALL(IW)%ILW) + EFLUX)
         ELSE
            QRADIN(IW) = E_WALL(IW)*(0.25_EB*UII(IJKW(1,IW),IJKW(2,IW),IJKW(3,IW)) + EFLUX)
         ENDIF   
      ENDIF
   ENDDO DROPLET_LOOP
ENDIF

END SUBROUTINE RADIATION_FVM

END SUBROUTINE COMPUTE_RADIATION
 
REAL(EB) FUNCTION BLACKBODY_FRACTION(L1,L2,TEMP)
 
! Calculates the fraction of black body radiation between wavelengths L1 and L2 (micron) in Temperature TEMP
 
REAL(EB) :: L1,L2,TEMP,LT1,LT2,BBFLOW,BBFHIGH,F
INTEGER  :: IYY

LT1    =   L1 * TEMP/LTSTEP
IF(LT1 > NLAMBDAT) THEN
   BBFLOW  = BBFRAC(NLAMBDAT)
ELSE
   IYY = MIN(NLAMBDAT-1,INT(LT1))
   F = LT1-REAL(IYY,EB)
   BBFLOW  = (1._EB-F)*BBFRAC(IYY) + F*BBFRAC(IYY+1)
ENDIF

LT2    =   L2 * TEMP/LTSTEP
IF(LT2 > NLAMBDAT) THEN
   BBFHIGH  = BBFRAC(NLAMBDAT)
ELSE
   IYY = MIN(NLAMBDAT-1,INT(LT2))
   F = LT2-REAL(IYY,EB)
   BBFHIGH  = (1._EB-F)*BBFRAC(IYY) + F*BBFRAC(IYY+1)
ENDIF
      
BLACKBODY_FRACTION = BBFHIGH - BBFLOW      

END FUNCTION BLACKBODY_FRACTION



FUNCTION GET_KAPPA(Z_IN,TYY,IBND)

! Returns the radiative absorption

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION_ALL,GET_MOLECULAR_WEIGHT
REAL(EB), INTENT(INOUT) :: Z_IN(0:N_TRACKED_SPECIES)
REAL(EB) :: KAPPA_TEMP,INT_FAC,GET_KAPPA,Y_MF(1:N_SPECIES),MWA
INTEGER, INTENT(IN) :: IBND,TYY
INTEGER :: LBND,UBND,N

CALL GET_MASS_FRACTION_ALL(Z_IN,Y_MF)
CALL GET_MOLECULAR_WEIGHT(Z_IN,MWA)
GET_KAPPA = 0._EB
DO N = 1, N_KAPPA_ARRAY
   IF (KAPPA_INDEX(N)==SOOT_INDEX) THEN
      INT_FAC = MAX(0._EB,(Y_MF(KAPPA_INDEX(N))*Z2KAPPA_M4(N)))**0.25_EB
   ELSE
      INT_FAC = MAX(0._EB,(Y_MF(KAPPA_INDEX(N))*MWA*Z2KAPPA_M4(N)))**0.25_EB
   ENDIF
   LBND = INT(INT_FAC)
   INT_FAC = INT_FAC - LBND   
   LBND = MIN(LBND,Z2KAPPA_M)
   UBND = MIN(LBND+1,Z2KAPPA_M)
   KAPPA_TEMP = Z2KAPPA(N,LBND,TYY,IBND)
   GET_KAPPA = GET_KAPPA + KAPPA_TEMP + INT_FAC*(Z2KAPPA(N,UBND,TYY,IBND)-KAPPA_TEMP)
ENDDO

END FUNCTION GET_KAPPA 

SUBROUTINE GET_REV_radi(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') radirev(INDEX(radirev,':')+1:LEN_TRIM(radirev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') radidate

END SUBROUTINE GET_REV_radi
 
END MODULE RAD

