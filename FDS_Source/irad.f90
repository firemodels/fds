!------------------------------------------------------------------------------
! This file contains 5 modules
! 1 - RADCONS
! 2 - RADCAL_VAR
! 3 - RADACL   
! 4 - SPECDATA
! 5 - MIEV
!------------------------------------------------------------------------------


MODULE RADCONS
!------------------------------------------------------------------------------
! Radiation Parameters

USE PRECISION_PARAMETERS
IMPLICIT NONE

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: DLN,ORIENTATION_FACTOR
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: BBFRAC, WL_LOW, WL_HIGH
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: DLX, DLY, DLZ, DLB, RSA

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: DLM
INTEGER, ALLOCATABLE, DIMENSION(:)   :: NRP

REAL(EB) :: RADTMP, PATH_LENGTH, RADIATIVE_FRACTION
REAL(EB) :: DGROUP_A, DGROUP_B, WEIGH_CYL
REAL(EB) :: DPHI0, FOUR_SIGMA, RPI_SIGMA, LTSTEP, RTMPMAX, RTMPMIN
REAL(EB) :: MIE_MINIMUM_DIAMETER,MIE_MAXIMUM_DIAMETER

INTEGER :: TIME_STEP_INCREMENT,NMIEANG
INTEGER :: NRDMIE, NLMBDMIE, MIE_NDG
INTEGER :: NRT,NCO,UIIDIM,NLAMBDAT,NKAPPAT,NKAPPAZ

LOGICAL :: WIDE_BAND_MODEL

INTEGER :: N_RADCAL_ARRAY_SIZE,RADCAL_SPECIES_INDEX(16),N_KAPPA_T=44,N_KAPPA_Y=50
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: Z2RADCAL_SPECIES
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: RADCAL_SPECIES2KAPPA
CHARACTER(LABEL_LENGTH) :: RADCAL_SPECIES_ID(16)='null'


!------------------------------------------------------------------------------
!
!     BBFRAC    Fraction of blackbody radiation
!     DLX       Mean X-component of the control angle ray vector
!     DLY       Mean Y-component of the control angle ray vector
!     DLZ       Mean Z-component of the control angle ray vector
!     DLB       Mean Bottom-component of rayn vector (cylindrical case)
!     DLM       Mirroring indexes
!     DLN       Wall normal matrix
!     DPHI0     Opening angle of the cylindrical domain
!     E_WALL    Wall emissivity
!     ILW       Radiation intensities on solid mirrors and mesh interfaces.
!               Intensity integrals (band specific or angle increment) for solid and open walls
!     INRAD_W   Incident radiative heat flux on a cell (QRADIN = E_WALL*INRAD_W)
!     R50       Array of PARTICLE radii corresponding to the median diameters 
!               of the distributions used in the generation of WQABS and WQSCA arrays.
!     MIE_NDG   Number of PARTICLE radii in WQABS and WQSCA arrays
!     NLMBDMIE  Number of wave lengths in Mie calculations
!     NMIEANG   Number of angle bins in forward scattering integration
!     NUMBER_RADIATION_ANGLES
!     NRA       Total number of radiation control angles
!     NRDMIE    Number of PARTICLE radii in Mie calculations
!     NRT       Number of radiation theta angles
!     NRP       Number of radiation phi angles on each theta band
!     NUMBER_SPECTRA_BANDS
!     NSB       Number of spectral bands (1=gray, 6=wide band, 9=wide band w. CH4)
!     OUTRAD_W  Emitted intensity from a wall (OUTRAD_W = QRADOUT/PI)
!     ORIENTATION_FACTOR
!               0.5 if the particle orientation most closely aligns with the solid angle; 0 otherwise
!     PHIUP     Upper limit of solid angle component PHI
!     PHILOW    Lower limit of solid angle component PHI
!     RADTMP    Radiation temperature for absorption properties (Mie)
!     QRADIN    Absorbed radiative heat flux into a surface cell (solid wall or open)
!     QRADOUT   Emitted radiative heat flux from a surface (solid wall or open) 
!     RSA       Array of solid angles
!     RTMPMAX   Maximum temperature for tabulation of radiative properties
!     RTMPMIN   Minimum temperature for tabulation of radiative properties
!     THETAUP   Upper limit of solid angle component THETA
!     THETALOW  Lower limit of solid angle component THETA
!     UII       Integrated intensity
!     UIID      Parts of UII  if WIDE_BAND_MODEL = TRUE, UIID contains the band specific intensity
!                             if WIDE_BAND_MODEL/= TRUE, UIID contains the ANGLE_INCREMENTs of intensity
!     WEIGH_CYL In cylindrical coordinates, all intensities represent two actual control angles
!     WL_LOW    Lower wavelength limit of the spectral bands
!     WL_HIGH   Upper wavelength limit of the spectral bands
!     WQABS     Absorption efficiency factor array 
!     WQSCA     Scattering efficiency factor array
!
!     PATH_LENGTH             Mean path length for the gray gas abs. coef.
!     ANGLE_INCREMENT         How many angles are skipped on each update
!     TIME_STEP_INCREMENT     How often is the radiation solver called
!
!------------------------------------------------------------------------------

END MODULE RADCONS


MODULE RADCAL_VAR

!------------------------------------------------------------------------------
! VARIABLES
! OMMIN : (REAL) MINIMUM WAVE NUMBER IN SPECTRUM, CM-1
! OMMAX : (REAL) MAXIMUM WAVE NUMBER IN SPECTRUM, CM-1
! NOM   : (INTEGER) NUMBER OF WAVELENGTH INTERVALS
! GC    : (REAL) COLLISION BROADENED HALF-WIDTH
! LAMBDA: (REAL) WAVELENGTH, IN MICROMETER 
!------------------------------------------------------------------------------
! Module wrapper for RadCal subroutine

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS, ONLY: AL2O3
IMPLICIT NONE

CHARACTER(255) :: MESSAGE
INTEGER, PARAMETER :: n_radcal_species = 16      ! number of radcal species (including soot) plus nitrogen and oxygen
INTEGER, PARAMETER :: n_radcal_species_gas=n_radcal_species-1 ! number of radcal species (including soot) plus nitrogen and oxygen

TYPE elements
   CHARACTER(32)           :: radcal_id
   CHARACTER(32)           :: id
   CHARACTER(2048)         :: comments
   REAL(EB), POINTER, DIMENSION(:) :: bands
   REAL(EB), POINTER, DIMENSION(:) :: temp_EXP
END TYPE elements

TYPE(elements), DIMENSION(n_radcal_species) :: radcal_species

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: partial_pressures_atm
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: temp_gas
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: segment_length_m
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: total_pressure_atm

CHARACTER(255) :: CHID_RADCAL  ! radcal CASE id
CHARACTER(255) :: TITLE_RADCAL ! radcal CASE title

REAL(EB) :: ommin, ommax, lambdamin, lambdamax
REAL(EB) :: twall, tau, xtot, fv

INTEGER:: npt

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: gamma, sd15, sd, sd7, sd3

INTEGER, PARAMETER :: n_temp_c3h6=7,  n_band_c3h6=3, &
                      n_temp_c3h8=7,  n_band_c3h8=2, &
                      n_temp_c7h16=7, n_band_c7h16=2, &
                      n_temp_c7h8=7,  n_band_c7h8=5, &
                      n_temp_ch4=23,  n_band_ch4=3, &
                      n_temp_ch3oh=7, n_band_ch3oh=4, &
                      n_temp_c5h8o2=7,n_band_c5h8o2=6, &
                      n_temp_c2h6=7,  n_band_c2h6=3, &
                      n_temp_c2h4=7,  n_band_c2h4=4

! 1: goody, 2: malkmus, 3: elsasser 
INTEGER, PARAMETER :: i_model_c3h6=2,i_model_c3h8=2,i_model_c7h16=2,i_model_c7h8=2,i_model_ch3oh=2,i_model_c5h8o2=2, &
                      i_model_c2h6=2,i_model_c2h4=2,i_model_co2=1,i_model_h2o=1,i_model_co=1,i_model_ch4=3

REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_c3h6_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_c3h8_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_c7h16_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_c7h8_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_ch3oh_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_c5h8o2_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_c2h6_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_ch4_temp
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: sd_c2h4_temp

REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_c3h6
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_c3h8
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_c7h16
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_c2h6
REAL(EB), ALLOCATABLE, DIMENSION(:)    :: be_c2h4

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_ch4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_c3h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_c3h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_c7h16
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_c2h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: om_bnd_c2h4

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_ch4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_ch4

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_c3h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_c3h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd3_c3h6

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_c3h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_c3h8

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_c7h16
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_c7h16

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd3_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd4_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd5_c7h8

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd3_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd4_ch3oh

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd3_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd4_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd5_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd6_c5h8o2

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_c2h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_c2h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd3_c2h6

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd1_c2h4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd2_c2h4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd3_c2h4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: sd4_c2h4

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_c3h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_c3h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad3_c3h6

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_c3h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_c3h8

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_c7h16
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_c7h16

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad3_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad4_c7h8
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad5_c7h8

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad3_ch3oh
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad4_ch3oh

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad3_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad4_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad5_c5h8o2
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad6_c5h8o2

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_c2h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_c2h6
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad3_c2h6

REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad1_c2h4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad2_c2h4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad3_c2h4
REAL(EB), ALLOCATABLE, DIMENSION(:,:)  :: gammad4_c2h4

REAL(EB), ALLOCATABLE, DIMENSION(:)  :: ab, wave_number, lambda
REAL(EB), ALLOCATABLE, DIMENSION(:)  :: incident_radiance

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: ttau

INTEGER :: nom

REAL(EB), PARAMETER :: m_to_cm = 100.0_EB  ! convertion factor meters -> centimeters
REAL(EB), PARAMETER :: cm_to_m = 0.01_EB   ! convertion factor centimeters -> meters

CHARACTER(30) :: radcal_id
REAL(EB) :: pfuel

INTEGER, PARAMETER :: i_co=3, i_co2=1, i_h2o=2, i_n2=14, i_o2=15, i_fv=16, i_c2h4=5, i_c2h6=6, i_c3h6=7, i_c3h8=8, &
                      i_c7h8=9, i_c7h16=10, i_ch3oh=11, i_ch4=4, i_ch4_old=13, i_mma=12

CHARACTER(255), PARAMETER :: iradid='$Id$'
CHARACTER(255), PARAMETER :: iradrev='$Revision$'
CHARACTER(255), PARAMETER :: iraddate='$Date$'
   
CONTAINS
 
!==============================================================================
SUBROUTINE GET_REV_irad(MODULE_REV,MODULE_DATE)
!==============================================================================
! Variables passed in

INTEGER, INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE
!------------------------------------------------------------------------------

WRITE(MODULE_DATE,'(A)') iradrev(INDEX(iradrev,':')+2:LEN_TRIM(iradrev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') iraddate

!------------------------------------------------------------------------------
END SUBROUTINE GET_REV_irad   

END MODULE RADCAL_VAR

   
MODULE RADCAL_CALC
   
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS, ONLY: SIGMA
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
USE RADCONS, ONLY: RPI_SIGMA
USE RADCAL_VAR

PUBLIC INIT_RADCAL,RCALLOC,RCDEALLOC,PLANCK,SUB_RADCAL

CONTAINS

!==============================================================================
SUBROUTINE init_radcal
!==============================================================================
IMPLICIT NONE
!------------------------------------------------------------------------------
! local variables

REAL(EB) :: DOM, omega

INTEGER :: i_wavenumb

IF((ommax<1100._EB).and.(ommin<ommax)) THEN 
   nom = INT((ommax-ommin)/5._EB)
ELSEIF(ommin>5000._EB) THEN
   nom = INT((ommax-ommin)/50._EB)
ELSEIF(ommin<1100._EB.and.ommax>5000._EB) THEN
   nom = INT((1100._EB-ommin)/5._EB)+INT((5000._EB-1100._EB)/25._EB) +INT((ommax-5000._EB)/50._EB)
ELSEIF(ommin<1100._EB) THEN
   nom = INT((1100._EB-ommin)/5._EB)+INT((ommax-1100._EB)/25._EB)
ELSEIF(ommax>5000._EB) THEN
   nom = INT((5000._EB-ommin)/25._EB)+INT((ommax-5000._EB)/50._EB)
ELSE
   nom = INT((ommax-ommin)/25._EB)
ENDIF

ALLOCATE(incident_radiance(nom)); incident_radiance   = 0.0_EB
ALLOCATE(ttau(0:npt,nom))       ; ttau                = 1.0_EB

ALLOCATE(lambda(nom))      ; lambda      = 0.0_EB
ALLOCATE(ab(nom))          ; ab          = 0.0_EB
ALLOCATE(wave_number(nom)) ; wave_number = 0.0_EB

! populate vector wave_number and lambda 

DOM    = 5.0_EB
omega  = ommin-DOM

DO i_wavenumb = 1, nom
   omega=omega+DOM
   IF(omega > 1100._EB) omega = omega + 20._EB
   IF(omega > 5000._EB) omega = omega + 25._EB

   wave_number(i_wavenumb) = omega            ! wavenumber in cm^-1
   lambda(i_wavenumb)      = 10000.0_EB/omega ! wavelength in micron
ENDDO

!------------------------------------------------------------------------------
END SUBROUTINE init_radcal


!==============================================================================
SUBROUTINE sub_radcal(effective_absorption,planck_mean_absorption,radiance,    &
                     total_transmissivity)
!==============================================================================
!
! input
!------
! io   : (INTEGER) ouput file unit

! output
!-------
! effective_absorption   : (REAL) effective absorption coefficient,   units: 1/cm
! planck_mean_absorption : (REAL) planck mean absorption coefficient, units: 1/cm
! radiance               : (REAL) total incident power per unit of area per 
!                                 unit of solid angle, units: w/m2/str
! total_transmissivity   : (REAL) total transmissivity, computed ONLY when twall 
!                                 (blackbody source) is strictly greater than 0
!                                 DIMENSIONless number
!
! local
!------
! a_collision : (REAL) collision broadened fine structure PARAMETER
!               (line_half_width/line_spacing)
! a_doppler   : (REAL) DOppler broadened fine structure PARAMETER
!               (line_half_width/line_spacing)
! azotemp     : (REAL) ratio reference temperature (273k) over local temperature
! omega       : (REAL) wavenumber
! optical_thickness: (REAL) optical thickness of the ith species,  units: cmstp
! ptot        : (REAL) total pressure in atm
! temp4       : (REAL) local temperature raised to the power 4
! 
! x_collision : (REAL) ! optical depth for a pure collision curve of growth
! x_doppler   : (REAL) ! optical depth for a pure DOppler curve of growth
! x_particle  : (REAL) ! optical depth for particle
!------------------------------------------------------------------------------

IMPLICIT NONE

! variables passed output
REAL(EB), INTENT(OUT) :: effective_absorption, planck_mean_absorption,       &
                        radiance, total_transmissivity

! local variables

REAL(EB), ALLOCATABLE, DIMENSION(:)   :: path_length_cm, azotemp

REAL(EB), ALLOCATABLE, DIMENSION(:)   :: taus, taul, x_particle, ab_planck,  &
                                       ab_tau, bb_spect, bb_gas, bb_wall

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: curtis_xstar, curtis_acollision,    &
                                       curtis_aDOppler, gc,                &
                                       optical_thickness, optical_depth
 
REAL(EB) :: rsl, rss, omega, wave_length,  dambda, lterm, temp4, & 
            avg_temp, total_length_cm, int_bb_gas, int_bb_wall

INTEGER :: i_wavenumb, kmax, kmin, i_path

! initialization

int_bb_gas            = 0.0_EB
int_bb_wall           = 0.0_EB
incident_radiance     = 0.0_EB
radiance              = 0.0_EB
total_transmissivity  = 1.0_EB  ! initialization values to arbitrary 
                                 ! large number so that user knows when it 
                                 ! has not been calculated
!------------------------------------------------------------------------------
! ALLOCATE arrays needed for calculation of radiation properties at a given 
! wavenumber: gc, species_pressure_atm, optical_thickness

IF(ALLOCATED(curtis_xstar))         DEALLOCATE(curtis_xstar)
ALLOCATE(curtis_xstar(1:n_radcal_species,0:npt));      curtis_xstar      = 0.0_EB

IF(ALLOCATED(curtis_acollision))    DEALLOCATE(curtis_acollision)
ALLOCATE(curtis_acollision(1:n_radcal_species,0:npt)); curtis_acollision = 0.0_EB

IF(ALLOCATED(curtis_aDOppler))      DEALLOCATE(curtis_aDOppler)
ALLOCATE(curtis_aDOppler(1:n_radcal_species,0:npt));   curtis_aDOppler   = 0.0_EB

IF(ALLOCATED(optical_thickness))    DEALLOCATE(optical_thickness)
ALLOCATE(optical_thickness(1:n_radcal_species,1:npt)); optical_thickness = 0.0_EB

IF(ALLOCATED(optical_depth))         DEALLOCATE(optical_depth)
ALLOCATE(optical_depth(1:n_radcal_species,0:npt));      optical_depth    = 0.0_EB
  
IF(ALLOCATED(gc))                   DEALLOCATE(gc)
ALLOCATE(gc(1:n_radcal_species,1:npt));                 gc               = 0.0_EB

IF(ALLOCATED(path_length_cm))       DEALLOCATE(path_length_cm)
ALLOCATE(path_length_cm(1:npt))

! To compute the effective absorption coefficient
IF (ALLOCATED(bb_wall)) DEALLOCATE(bb_wall)
ALLOCATE(bb_wall(nom)); bb_wall = 0.0_EB

IF (ALLOCATED(bb_gas)) DEALLOCATE(bb_gas)
ALLOCATE(bb_gas(nom));  bb_gas  = 0.0_EB

! convert path length from meters to cm
path_length_cm = m_to_cm*segment_length_m

! compute ratio standard temperature over considered temperature
IF(ALLOCATED(azotemp))              DEALLOCATE(azotemp)
ALLOCATE(azotemp(1:npt))
azotemp = 273._EB/temp_gas

!------------------------------------------------------------------------------
! compute the segments properties that are indepENDent of the frequency
! store them in arrays with DIMENSIONs (npt,n_radcal_species)

init_path: DO i_path = 1, npt  ! loop over elements of the path length. 

!------------------------------------------------------------------------------
! test. IF path_length_cm set to zero (bad segment), CYCLE
! ELSE, proceed with calculation

   IF (path_length_cm(i_path)<=TWO_EPSILON_EB) CYCLE

!------------------------------------------------------------------------------
! compute species optical thickness and collision broadening

   optical_thickness(:,i_path) = azotemp(i_path)*path_length_cm(i_path)*               &
                                                partial_pressures_atm(:,i_path)

   gc(:,i_path)                = collision_broadening( partial_pressures_atm(:,i_path), &
                                                      total_pressure_atm(i_path),      &
                                                      azotemp(i_path)                )

ENDDO init_path

! compute total length of the line of sight
total_length_cm = SUM(path_length_cm)

! compute length weighted average temp4 and avg_temp along line of sight

temp4    = SUM(temp_gas**4*path_length_cm)/total_length_cm
avg_temp = SUM(temp_gas*path_length_cm)/total_length_cm

!------------------------------------------------------------------------------
! loop spectrum computes each spectral contribution,
! loop over each wave number intervals
! this loop can be parallized with OPENmp

!!                   below are OPENmp instructions
!$OMP PARALLEL PRIVATE (i_wavenumb, i_path) &
!$ firstPRIVATE(curtis_xstar,curtis_acollision,curtis_aDOppler,optical_depth,xtot)

!$OMP DO
lspectrum: DO i_wavenumb = 1, nom
   curtis_xstar      = 0.0_EB
   curtis_aDOppler   = 0.0_EB
   optical_depth     = 0.0_EB
   curtis_acollision = 0.0_EB
!------------------------------------------------------------------------------
! loop over the dIFferent segments of the line of sight

   path: DO i_path = 1,npt  ! loop over elements of the path length. 
! evaluate the combined spectral transmissivity and radiance
! arguments out: x (vector)
      CALL species_optical_depth(path_length_cm(i_path), partial_pressures_atm(:,i_path),   & 
                                 total_pressure_atm(i_path), gc(:,i_path), temp_gas(i_path),&
                                 wave_number(i_wavenumb), optical_thickness(:,i_path),      &
                                 curtis_xstar(:,i_path-1:i_path),                           &
                                 curtis_acollision(:,i_path-1:i_path),                      &
                                 curtis_aDOppler(:,i_path-1:i_path),                        &
                                 optical_depth(:,i_path-1:i_path))
      xtot = SUM(optical_depth(:,i_path))   
      ttau(i_path,i_wavenumb) = EXP(-xtot)
   ENDDO path

   ! ab : weakline limits optical depth per unit length, units in cm-1
   IF (total_length_cm.gt.0.0_EB) THEN
      ab(i_wavenumb) = SUM(curtis_xstar(:,npt))/total_length_cm
   ELSE
      ab(i_wavenumb) = 0.0_EB
   ENDIF

! integrate the incident radiance over the whole line of sight
   incident_radiance(i_wavenumb) = 0.0_EB

   DO i_path = 1, npt
      ! reCALL: ttau(0,i_wavenumb) = 1.0_EB

      incident_radiance(i_wavenumb) = incident_radiance(i_wavenumb) +                       &
                                    (ttau(i_path-1,i_wavenumb)-ttau(i_path,i_wavenumb))*  &
                                    planck_wn(temp_gas(i_path),wave_number(i_wavenumb))
   ENDDO

   incident_radiance(i_wavenumb) = incident_radiance(i_wavenumb) +  &
                                 ttau(npt,i_wavenumb)*planck_wn(twall,wave_number(i_wavenumb))
 
   bb_gas(i_wavenumb)  = planck_wn(avg_temp,wave_number(i_wavenumb))
   bb_wall(i_wavenumb) = planck_wn(twall,   wave_number(i_wavenumb))
ENDDO lspectrum
!$OMP END DO
!$OMP END PARALLEL

!------------------------------------------------------------------------------   
! integrate the radiance over the spectrum

radiance    = integration(wave_number,incident_radiance)
int_bb_gas  = integration(wave_number,bb_gas)
int_bb_wall = integration(wave_number,bb_wall)

!------------------------------------------------------------------------------
! soot contribution
! determine soot radiance for short (between ommax and 25000) and 
! long wavelengths (between 5 and ommin)
! this corresponds to the ranges between 400 nm (visible) to ommax (user depENDent)
! and ommin (user depENDent) and 2000 microns
     
rsl  = 0._EB
rss  = 0._EB

kmax = INT(ommin)
kmin = INT(ommax)

!------------------------------------------------------------------------------
long: DO i_wavenumb = 5, kmax, 5 ! loop over long wavelengths

   omega       = REAL(i_wavenumb,EB)
   wave_length = 10000._EB/omega
   dambda      = 10000._EB/(omega-2.5_EB)-10000._EB/(omega+2.5_EB)

! initialize taul: transmissivity for long wavelengths

   IF (ALLOCATED(taul)) DEALLOCATE(taul)
   ALLOCATE(taul(0:npt));     taul       = 1.0_EB
 
   IF (ALLOCATED(x_particle)) DEALLOCATE(x_particle)
   ALLOCATE(x_particle(npt)); x_particle = 0.0_EB

   ! loop over the line of sight elements
   DO i_path = 1, npt
    
      CALL pod(omega,partial_pressures_atm(i_fv,i_path),path_length_cm(i_path), x_particle(i_path))  

      taul(i_path) = EXP(-SUM(x_particle(1:i_path)))      
      rsl = rsl + (taul(i_path-1)-taul(i_path))*planck(temp_gas(i_path),wave_length)*dambda

   ENDDO   

! add contribution of the incident radiance from the wall
   rsl = rsl + taul(npt)*planck(twall,wave_length)*dambda

! calculate contribution to int_bb_wall & int_bb_gas
   int_bb_wall = int_bb_wall + planck(twall,    wave_length)*dambda
   int_bb_gas  = int_bb_gas  + planck(avg_temp, wave_length)*dambda

ENDDO long

!------------------------------------------------------------------------------
short: DO i_wavenumb = kmin, 25000, 100 ! loop over short wavelengths

   omega       = REAL(i_wavenumb,EB)
   wave_length = 10000._EB/omega
   dambda      = 10000._EB/(omega-50._EB)-10000._EB/(omega+50._EB)

! initialize taus: transmissivity for short wavelengths

   IF (ALLOCATED(taus)) DEALLOCATE(taus)
   ALLOCATE(taus(0:npt));     taus       = 1.0_EB
 
   IF (ALLOCATED(x_particle)) DEALLOCATE(x_particle)
   ALLOCATE(x_particle(npt)); x_particle = 0.0_EB

   ! loop over the line of sight elements

   DO i_path = 1, npt
    
      CALL pod(omega,partial_pressures_atm(i_fv,i_path),path_length_cm(i_path), x_particle(i_path))  

      taus(i_path) = EXP(-SUM(x_particle(1:i_path)))      
      rss = rss + (taus(i_path-1)-taus(i_path))*planck(temp_gas(i_path),wave_length)*dambda

   ENDDO   

! add contribution of the incident radiance from the wall
   rss = rss + taus(npt)*planck(twall,wave_length)*dambda

! calculate contribution to int_bb_wall & int_bb_gas
   int_bb_wall = int_bb_wall + planck(twall,    wave_length)*dambda
   int_bb_gas  = int_bb_gas  + planck(avg_temp, wave_length)*dambda

ENDDO short

! END loop over species
!------------------------------------------------------------------------------

! add the contribution of soot over larger wave length DOMain
radiance = radiance + rss + rsl

!------------------------------------------------------------------------------
! DEALLOCATE some variables (gc,x,optical_thickness)

IF(ALLOCATED(bb_gas))            DEALLOCATE(bb_gas)
IF(ALLOCATED(bb_wall))           DEALLOCATE(bb_wall)
IF(ALLOCATED(gc))                DEALLOCATE(gc)
IF(ALLOCATED(optical_depth))     DEALLOCATE(optical_depth)
IF(ALLOCATED(optical_thickness)) DEALLOCATE(optical_thickness)
IF(ALLOCATED(curtis_xstar))      DEALLOCATE(curtis_xstar)
IF(ALLOCATED(curtis_acollision)) DEALLOCATE(curtis_acollision)
IF(ALLOCATED(curtis_aDOppler))   DEALLOCATE(curtis_aDOppler)
IF(ALLOCATED(x_particle))        DEALLOCATE(x_particle)
IF(ALLOCATED(taus))              DEALLOCATE(taus)
IF(ALLOCATED(taul))              DEALLOCATE(taul)
!------------------------------------------------------------------------------

IF (total_length_cm > 0.0_EB) THEN

   IF(ABS(twall-avg_temp)<=spacing(twall)) THEN 
      lterm = 1.0_EB
   ELSE
      lterm = MIN(1.0_EB,MAX(1e-30_EB,(radiance-int_bb_gas)/(int_bb_wall-int_bb_gas)))
   ENDIF
   !IF (lterm.gt.1.0_EB) THEN
   !   WRITE(*,*) "lterm: ", lterm
   !   WRITE(*,*) "radiance/rpi_sigma ", radiance/rpi_sigma
   !   WRITE(*,*) "twall**4 ", twall**4
   !   WRITE(*,*) "temp**4  ", temp4
   !ENDIF
   
   effective_absorption = -1._EB/total_length_cm*LOG(lterm)
ELSE
   effective_absorption = 0.0_EB
ENDIF 

!------------------------------------------------------------------------------
! computes the mean absorption coefficients over the whole pathlength
! integrate planck_mean_absorption*planck over lambda

IF (ALLOCATED(ab_planck)) DEALLOCATE(ab_planck)
ALLOCATE(ab_planck(nom))
 
DO i_wavenumb = 1, nom
   ab_planck(i_wavenumb) = ab(i_wavenumb)*planck(avg_temp,lambda(i_wavenumb))
ENDDO

planck_mean_absorption = integration(lambda,ab_planck)

! divide by sigma*pi*temp**4 to get final coefficient planck_mean_absorption 

IF (temp4 > TWO_EPSILON_EB) THEN
   planck_mean_absorption = planck_mean_absorption/(rpi_sigma*temp4)
ELSE
   planck_mean_absorption = 0.0_EB
ENDIF

!------------------------------------------------------------------------------
! computes the total transmissivity over the whole pathlength. it is calculated 
! ONLY IF twall is strictly greater than 0
! integrate tau*planck over lambda

IF (twall**4 > TWO_EPSILON_EB) THEN

   IF (ALLOCATED(ab_tau)) DEALLOCATE(ab_tau)
   ALLOCATE(ab_tau(nom))
 
   IF (ALLOCATED(bb_spect)) DEALLOCATE(bb_spect)
   ALLOCATE(bb_spect(nom))


   DO i_wavenumb = 1, nom
      bb_spect(i_wavenumb) = planck_wn(twall,wave_number(i_wavenumb))
      ab_tau(i_wavenumb)   = ttau(npt,i_wavenumb)*bb_spect(i_wavenumb)
   ENDDO

   total_transmissivity = integration(wave_number,ab_tau)/integration(wave_number,bb_spect)

ENDIF

!------------------------------------------------------------------------------
! DEALLOCATE ab_planck & ab_tau
IF (ALLOCATED(ab_planck)) DEALLOCATE(ab_planck)
IF (ALLOCATED(ab_tau))    DEALLOCATE(ab_tau)
IF (ALLOCATED(bb_spect))  DEALLOCATE(bb_spect)

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE sub_radcal


!==============================================================================
FUNCTION collision_broadening(species_pressure_atm,ptot,azotemp) RESULT(gc)
!==============================================================================
!! this FUNCTION computes the collision broadening half-width at half-height. 
!! (calculation proceeds in accordance with the slg model, table 5-18, in nasa sp-3080._EB)
!! based on eq 5-34 from nasa report sp-3080
!! RETURNs vector gc that CONTAINS the collision broadening half-width at 
!! half-height including foreign and self-broadening collision coefficients.
!! gc units in cm-1
! checked on 02/11/13
!------------------------------------------------------------------------------

! arguments in
REAL(EB), DIMENSION(:), INTENT(IN) :: species_pressure_atm ! species partial pressure, units: atm
REAL(EB), INTENT(IN) :: ptot         ! total pressure, units: atm  
REAL(EB), INTENT(IN) :: azotemp      ! reCALL azotemp = 273/t, units: NONE

! argument out
REAL(EB), ALLOCATABLE, DIMENSION(:) :: gc ! units in cm-1

! local variables
REAL(EB) :: p_fuel

INTEGER, DIMENSION(6) :: i_indices
INTEGER :: i, ii, i_fuel

i_fuel = i_c2h4

! populate i_indices used for computations of gc
! [note: the total intensity calculated is that which leaves interval j=1.
! p(i,j) is partial pressure, atm, of species i in  interval j.
! i=1,2,3,4,5, or 6 implies species is co2, h2o, fuel, co, o2, or n2, resp.]

i_indices(1) = i_co2
i_indices(2) = i_h2o 
i_indices(3) = i_fuel     ! fuel
i_indices(4) = i_co
i_indices(5) = i_o2
i_indices(6) = i_n2

! [note: the total intensity calculated is that which leaves interval j=1.
! p(i,j) is partial pressure, atm, of species i in  interval j.
! i=1,2,3,4,5, or 6 implies species is co2, h2o, ch4, co, o2, or n2, resp.]

IF(ALLOCATED(gc)) DEALLOCATE(gc)
ALLOCATE(gc(n_radcal_species))
gc = 0.0_EB

! compute fuel pressure. remove contribution of non fuel species: 
! co2, co, h2o, o2, n2

p_fuel = ptot - (species_pressure_atm(i_co2)+species_pressure_atm(i_co)+      &
               species_pressure_atm(i_h2o)+species_pressure_atm(i_n2)+      &
               species_pressure_atm(i_o2))

! compute collisional broadening half-width at half height, eq 5-34 nasa report sp-3080
! initialization

! co2(ii=1), h2o(ii=2), fuel(ii=3), co(ii=4), o2 (ii=5), n2 (ii=6) 
DO i = 1, 4

! IF the species is not present, i.e. species_pressure_atm(i_indices(i))==0.0_EB
! THEN set gc to zero and CYCLE. ELSE proceed with calculation

   IF (species_pressure_atm(i_indices(i))==0.0_EB.and.i/=3) CYCLE

! include non-resonant foreign and self-broadening collisions
   DO ii = 1, 6
      IF (i_indices(ii) == i_fuel) THEN
         gc(i_indices(i)) = gc(i_indices(i))+gamma(i,ii)*p_fuel*SQRT(azotemp)
      ELSE
         gc(i_indices(i)) = gc(i_indices(i))+gamma(i,ii)*species_pressure_atm(i_indices(ii))*SQRT(azotemp)
      ENDIF
   ENDDO
! include resonant self-broadening collisions
   IF (i_indices(i) /= i_fuel) THEN
      gc(i_indices(i))=gc(i_indices(i))+gamma(i,7)*species_pressure_atm(i_indices(i))*azotemp     
   ELSE
      gc(i_indices(i))=gc(i_indices(i))+gamma(i,7)*p_fuel*azotemp
   ENDIF
ENDDO

! compute the collision broadening coefficient for the fuel

DO i = 1, n_radcal_species
   IF (i/=i_co2.and.i/=i_o2.and.i/=i_co.and.i/=i_h2o.and.i/=i_n2.and.i/=i_fv) &
      gc(i) = gc(i_fuel)
ENDDO

!------------------------------------------------------------------------------
END FUNCTION collision_broadening

!==============================================================================
SUBROUTINE species_optical_depth(path_length_cm,p_atm,ptot,gc,temp,omega,          &
                              optical_thickness, curtis_xstar,                  &
                              curtis_acollision, curtis_aDOppler,optical_depth)
!==============================================================================
!! this SUBROUTINE computes the species optical_depth for all the species 
!! given a wavenumber omega and species optical_thickness.
!! RETURNs the optical depth for each species using the curtis godson approximation 
! 
! 
! arguments in:
! ------------
! path_length_cm   : scalar, physical length, in cm
! p_atm            : vector, (nradcal_species) species partial pressure, in atm
! ptot             : scalar, total pressure, in atm
! gc               : vector, collision-broadened half-width values, in cm-1
! temp             : scalar, local temperature, in k
! omega            : scalar, wavenumber, in cm-1
! optical_thickness: vector (nradcal_species) elements, units in stp.cm.atm
!
! arguments inout:
!-----------------
! varibles used to compute the curtis-godson approximation PARAMETERs
! the first values (i.e.) (nradcal_species,1) was either computed previously
! or is equal to 0 (for the first point)
! 
! curtis_xstar      : array of values of xstart.     DIMENSION: (nradcal_species,2).
! curtis_acollision : array of values of acollision. DIMENSION: (nradcal_species,2).
! curtis_aDOppler   : array of values of aDOppler.   DIMENSION: (nradcal_species,2).
! optical_depth     : array of values of elements optical depth, DIMENSIONless
!                     DIMENSION: (nradcal_species,2).
!
!------------------------------------------------------------------------------

IMPLICIT NONE

! variables passed in
REAL(EB), DIMENSION(:), INTENT(IN) :: p_atm, optical_thickness, gc
REAL(EB), INTENT(IN) :: ptot, temp, omega, path_length_cm

! variables passed inout
REAL(EB), DIMENSION(:,:), INTENT(INOUT) :: curtis_xstar, curtis_acollision,   &
                                          curtis_aDOppler, optical_depth

! local
REAL(EB) :: x_doppler, x_collision, xstar
REAL(EB) :: sdweak,    a_collision, a_doppler

INTEGER  :: i_model ! model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser
INTEGER  :: i

! initialize x_doppler, x_collision
x_doppler   = 0.0_EB
x_collision = 0.0_EB

! initialize curtis_xstar(:,2), curtis_acollision(:,2), 
!            curtis_aDOppler(:,2), optical_depth(:,2)
! these values will be updated IF sdweak of a given species is greater than TINY_EB

curtis_xstar(:,2)      = curtis_xstar(:,1)
curtis_acollision(:,2) = curtis_acollision(:,1)  
curtis_aDOppler(:,2)   = curtis_aDOppler(:,1)
optical_depth(:,2)     = optical_depth(:,1)

! loop species computes the contribution of each species to tau

lspecies: DO i = 1, n_radcal_species
   IF (p_atm(i) <= TWO_EPSILON_EB) CYCLE lspecies

   i_model = 1 ! enforce goody model when not specIFied  

   sdweak      = TINY_EB  ! average spectral absorption  coefficient. units: (atm^-1.cm^-1)
   a_collision = 1.0_EB  ! fine structure PARAMETER. DIMENSIONless
   a_doppler   = 1.0_EB

   ! IF the species is present in this segment, get sdweak, a_collision, a_doppler.
   ! otherwise USE 0.0_EB, 1.0_EB, 1.0_EB
 
   SELECT CASE (i)
      CASE (i_co2)
         CALL co2(    omega,temp,gc(i_co2), sdweak,a_collision,a_doppler,i_model)
      CASE (i_h2o)
         CALL h2o(    omega,temp,gc(i_h2o), sdweak,a_collision,a_doppler,i_model)
      CASE (i_co)
         CALL co(     omega,temp,gc(i_co),  sdweak,a_collision,a_doppler,i_model)
      CASE (i_ch4_old)
         CALL ch4_old(omega,temp,p_atm(i_ch4_old),  ptot,gc(i_ch4_old),sdweak,a_collision,a_doppler,i_model)
      CASE (i_ch4)
         CALL ch4(    omega,temp,p_atm(i_ch4),  ptot,gc(i_ch4),sdweak,a_collision,a_doppler,i_model)
      CASE (i_c3h8)
         CALL c3h8(   omega,temp,p_atm(i_c3h8), ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_c7h16)
         CALL c7h16(  omega,temp,p_atm(i_c7h16),ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_ch3oh)
         CALL ch3oh(  omega,temp,p_atm(i_ch3oh),ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_c7h8)
         CALL c7h8(   omega,temp,p_atm(i_c7h8), ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_c3h6)
         CALL c3h6(   omega,temp,p_atm(i_c3h6), ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_mma)
         CALL c5h8o2( omega,temp,p_atm(i_mma),  ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_c2h6)
         CALL c2h6(   omega,temp,p_atm(i_c2h6), ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_c2h4)
         CALL c2h4(   omega,temp,p_atm(i_c2h4), ptot,sdweak,a_collision,a_doppler,i_model)
      CASE (i_fv)
      CASE DEFAULT
         CYCLE lspecies
      END SELECT

   ! IF (sdweak <= TINY_EB) THEN the species DOes not participate on this segments. skip the calculation of 
   ! optical_depth(:,2). IF not, proceed with calculation of x_collision, x_doppler, etc...

   IF (i/=i_fv) THEN

      IF(sdweak <= TINY_EB) CYCLE lspecies

      xstar             = sdweak*optical_thickness(i) ! optical depth for the weak line limit.  sp-3080 eq. 5-27
      curtis_xstar(i,2) = curtis_xstar(i,1) + xstar
      IF (curtis_xstar(i,2)>=TINY_EB) THEN
      ! collision broadening
         curtis_acollision(i,2) = 1.0_EB/curtis_xstar(i,2)*(curtis_acollision(i,1)*curtis_xstar(i,1)+ &
         a_collision*xstar)
      ! DOppler broadening
         curtis_aDOppler(i,2) = 1.0_EB/curtis_xstar(i,2)*(curtis_aDOppler(i,1)*curtis_xstar(i,1)+ &
         a_doppler*xstar)
      ELSE
         curtis_acollision(i,2) = 1.0_EB
         curtis_aDOppler(i,2)   = 1.0_EB
      ENDIF

      IF (curtis_xstar(i,2) >= 1.e-20_EB .AND. curtis_acollision(i,2)>0._EB) THEN
!------------------------------------------------------------------------------
! compute the optical depth for the ith species using ranDOM band model
! composed of lines of DOppler-lorentz shape with:
! ranDOM distribution of either equal strength lines for DOppler or decaying
! USE appropriate model for lorentz line (elsasser,goody,malkmus)
! compute optical depth for pure DOppler curve of growth 
         x_doppler = growth_doppler(curtis_xstar(i,2),curtis_aDOppler(i,2),i_model)
 
! lorentz optical depth 
! compute optical depth for pure collision curve of growth 
         SELECT CASE (i_model)
            CASE(1)  ! goody model
               x_collision = goody(curtis_xstar(i,2),curtis_acollision(i,2))
            CASE(2)  ! malkmus model
               x_collision = malkmus(curtis_xstar(i,2),curtis_acollision(i,2))
            CASE(3)  ! elsasser model
               x_collision = elsasser(curtis_xstar(i,2),curtis_acollision(i,2))
            CASE DEFAULT
               x_collision = goody(curtis_xstar(i,2),curtis_acollision(i,2))
         END SELECT

      ! compute optical depth
         optical_depth(i,2) = combined_lines(x_collision, x_doppler, curtis_xstar(i,2))

      ELSE ! xstar < 1.e-20_EB very weak line
         optical_depth(i,2) = curtis_xstar(i,2)
      ENDIF ! condition (xstar >= 1.e-20_EB)

   ! determine optical depth of soot
   ELSE ! treatment for soot volume fraction
      CALL pod(omega, p_atm(i),path_length_cm,xstar)
      curtis_xstar(i,2)  = curtis_xstar(i,1) + xstar
      optical_depth(i,2) = curtis_xstar(i,2)

   ENDIF

END DO lspecies

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE species_optical_depth

!==============================================================================
ELEMENTAL FUNCTION growth_doppler(xstar,a_doppler,i_model) RESULT(x_doppler)
!==============================================================================
!! this FUNCTION computes the DOppler curve of growth.
!! two CASEs
!! 1) equally intense, ranDOMly spaced DOppler line (i_model =1 or 3)
!! 2) geometriCALLy decaying strength, ranDOMly spaced DOppler line (i_model=2)
!
! variables in:
! xstar: optical pathlength * mean absorption coefficient
! a_doppler: DOppler fine structure PARAMETER
! i_model: model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser
!
! variable out
! x_doppler

! variables in:
REAL(EB), INTENT(IN) :: xstar
REAL(EB), INTENT(IN) :: a_doppler
INTEGER,  INTENT(IN) :: i_model

! variable out
REAL(EB) :: x_doppler

! source: ! sp-3080, chap 3.2.5 e and f. see page 105.

SELECT CASE(i_model)
CASE(1) ! goody model: equally intense, ranDOMly spaced DOppler line 
      x_doppler = 1.7_EB*a_doppler*SQRT(LOG(1._EB+(0.58824_EB*xstar/a_doppler)**2)) 
CASE(2) ! malkmus: geometriCALLy decaying strength, ranDOMly spaced DOppler line
      x_doppler = 0.937_EB*a_doppler*(LOG(1.0_EB+(1.07_EB*xstar/a_doppler)**(0.6666_EB)))* &
                                 SQRT(LOG(1.0_EB+(1.07_EB*xstar/a_doppler)**(0.6666_EB)))
CASE(3)
      x_doppler = 1.7_EB*a_doppler*SQRT(LOG(1._EB+(0.58824_EB*xstar/a_doppler)**2)) 
CASE DEFAULT
      x_doppler = 1.7_EB*a_doppler*SQRT(LOG(1._EB+(0.58824_EB*xstar/a_doppler)**2)) 
END SELECT
!------------------------------------------------------------------------------
END FUNCTION growth_doppler

!==============================================================================
ELEMENTAL FUNCTION combined_lines(x_collision, x_doppler, xstar) RESULT(optical_depth)
!==============================================================================
!! compute the combined collision and DOppler optical depths and 
! the optical depth based on EXPression given in sp-3080 eq 5-26 & eq 5-25 (p220)
!
! variables in:
! x_collision : lorenz curve of growth
! x_doppler: DOppler curve of growth
! xstar: optical pathlength * mean absorption coefficient
!
! variable out
! optical_depth

! variables passed in
REAL(EB), INTENT(IN) :: x_collision
REAL(EB), INTENT(IN) :: x_doppler
REAL(EB), INTENT(IN) :: xstar

! variable passed output
REAL(EB) :: optical_depth

! local variables
REAL(EB) :: y_doppler, y_collision, y_combined

! compute the combined collision and DOppler optical depths y_combined

y_doppler   = MAX(1._EB-(x_doppler/xstar)**2,TWO_EPSILON_EB)
y_collision = MAX(1._EB-(x_collision/xstar)**2,TWO_EPSILON_EB)
! sp-3080 eq. 5-26
y_combined = MAX(1._EB/y_collision**2+1._EB/y_doppler**2-1._EB,1._EB)
! finally compute optical depth
! sp-3080 eq. 5-25
optical_depth = xstar*SQRT(1._EB-(y_combined**(-0.5_EB))) 

!------------------------------------------------------------------------------
END FUNCTION combined_lines


!==============================================================================
ELEMENTAL FUNCTION goody(xstar,a_collision) RESULT(x_goody)
!==============================================================================
!! FUNCTION computes the equivalent line width over the average line spacing
!! using the goody statistical model
! 
! variables in:
! xstar       : optical pathlength * mean absorption coefficient
! a_collision : fine structure PARAMETER
!
! RETURNs x_goody
  
IMPLICIT NONE

REAL(EB), INTENT(IN) :: xstar
REAL(EB), INTENT(IN) :: a_collision

REAL(EB)  :: x_goody

IF (xstar<1.0e-15_EB) THEN ! weak line limit
x_goody = xstar
ELSEIF(xstar>1.0e+5_EB) THEN ! strong line limit
x_goody = SQRT(4.0_EB*xstar*a_collision)
ELSE                        ! formulation to avoid division by zero
x_goody = SQRT(a_collision)*xstar/SQRT(a_collision+0.25_EB*xstar)
ENDIF

RETURN
!------------------------------------------------------------------------------
END FUNCTION goody

!==============================================================================
ELEMENTAL FUNCTION malkmus(xstar,a_collision) RESULT(x_malkmus)
!==============================================================================
!
!! FUNCTION computes the equivalent line width over the average line spacing
!! using the malkmus statistical model
! 
! variables in:
! xstar       : optical pathlength * mean absorption coefficient
! a_collision : fine structure PARAMETER
!
! RETURNs x_malkmus (optical depth for a pure collision curve)
  
IMPLICIT NONE

! variables passed in
REAL(EB), INTENT(IN) :: xstar
REAL(EB), INTENT(IN) :: a_collision

REAL(EB)  :: x_malkmus

IF    (xstar < 1.0e-15_EB) THEN ! weak line limit
x_malkmus = xstar
ELSEIF(xstar > 1.0e+5_EB ) THEN ! strong line limit
x_malkmus = SQRT(4.0_EB*xstar*a_collision)
ELSE             ! formulation to avoid division by a_collision IF equal to zero
x_malkmus = SQRT(4.0_EB*a_collision**2+4.0_EB*a_collision*xstar)-2.0_EB*a_collision
ENDIF

RETURN
!------------------------------------------------------------------------------
END FUNCTION malkmus

!==============================================================================
FUNCTION elsasser(xstar,a_collision) RESULT(x_elsasser)
USE MATH_FUNCTIONS, ONLY : ERF
!==============================================================================
!! FUNCTION computes the equivalent line width over the average line spacing
!! using the elsasser statistical model
! 
! variables in:
! xstar       : optical pathlength * mean absorption coefficient
! a_collision : fine structure PARAMETER
!
! RETURNs x_elsasser
!------------------------------------------------------------------------------

IMPLICIT NONE

! variables passed in
REAL(EB), INTENT(IN) :: xstar
REAL(EB), INTENT(IN) :: a_collision

! variables passed out
REAL(EB)  :: x_elsasser

! local variables

REAL(EB) :: x_collision

! the following computes the optical depth, x_collision, for any species represented 
! by the elsasser model using the godson equation and an approximation to the ladenberg-reiche
! FUNCTION as recommENDed by brosmer and tien (jqsrt 33,p 521), eq. 2, page 523 

x_collision = SQRT(4.0_EB*a_collision)*xstar/SQRT(16.0_EB/pi*a_collision+xstar)

x_elsasser  = -LOG(1.0_EB-ERF(x_collision)+TWO_EPSILON_EB)

!------------------------------------------------------------------------------
RETURN
END FUNCTION elsasser

!==============================================================================
pure FUNCTION partition_FUNCTION(transition_wn,degeneracy,temp)
!==============================================================================
!! This FUNCTION calculate the partition_FUNCTION of an harmonic 
!! quantum oscillator knowing the dIFferent frequencies and degeneracies
!------------------------------------------------------------------------------

IMPLICIT NONE

! Variable passed in 
REAL(EB), DIMENSION(:), INTENT(IN) :: transition_wn ! Wavenumber associated with
                                                   ! a vibration mode

INTEGER, DIMENSION(:), INTENT(IN) :: degeneracy  ! Denegeracy of energy level
                                                ! associated with a vibration
                                                ! mode

REAL(EB), INTENT(IN) :: temp                     ! Temperature in Kelvin

! Variable passed out
REAL(EB) :: partition_FUNCTION

! Local variables
REAL(EB), PARAMETER :: q2 = 1.4388_EB   ! speed_of_light*planck_cns/boltzmann
REAL(EB) :: q2_over_T
INTEGER  :: n_mode
INTEGER  :: i_mode

q2_over_T = q2/temp 
n_mode    = SIZE(transition_wn)

partition_FUNCTION = 1.0_EB

DO i_mode = 1, n_mode
partition_FUNCTION = partition_FUNCTION* &
         (1.0_EB-EXP(-q2_over_T*transition_wn(i_mode)))**(degeneracy(i_mode))
ENDDO

!------------------------------------------------------------------------------
RETURN
END FUNCTION partition_FUNCTION


!==============================================================================
pure SUBROUTINE approx_vib_rot(transition_matrix,line_strength,om,temp,omega,  &
                              bcnt,be,gc1, gd, gdinv,gddinv,sdweak,line_spacing)
!==============================================================================
!! just overlapping spectral lines
! 
! 
!------------------------------------------------------------------------------

IMPLICIT NONE

! Variable passed in
INTEGER, DIMENSION(:,:), INTENT(IN) :: transition_matrix ! Square matrix
REAL(EB), DIMENSION(:),   INTENT(IN) :: bcnt ! Band center in wavenumbers (cm-1)
REAL(EB), DIMENSION(:),   INTENT(IN) :: om
REAL(EB), DIMENSION(:),   INTENT(IN) :: line_strength

REAL(EB), INTENT(IN) :: omega ! Wavenumber of sought properties
REAL(EB), INTENT(IN) :: be    ! Rotational constant
REAL(EB), INTENT(IN) :: temp  ! Temperature in Kelvin
REAL(EB), INTENT(IN) :: gc1
REAL(EB), INTENT(IN) :: gd
REAL(EB), INTENT(IN), optional :: line_spacing ! Optional

! Variables passed out
REAL(EB), INTENT(OUT) :: sdweak
REAL(EB), INTENT(OUT) :: gdinv
REAL(EB), INTENT(OUT) :: gddinv

! local Variables
REAL(EB), DIMENSION(:), ALLOCATABLE :: transition_wn, atot
REAL(EB), PARAMETER :: q2 = 1.4388_EB   ! speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: t0 = 300._EB     ! reference temperature, in kelvin

REAL(EB) :: t0ot
REAL(EB) :: q2ot
REAL(EB) :: dinv

INTEGER  :: n_transitions
INTEGER  :: n_vibrations

INTEGER :: i_transition

sdweak  =  0.0_EB
t0ot    =  t0/temp
q2ot    = -q2/temp

n_transitions = SIZE(transition_matrix,2)
n_vibrations  = SIZE(transition_matrix,1)

IF (.not.present(line_spacing)) THEN 
dinv = 1.0_EB/(4.0_EB*be)
ELSE
dinv = line_spacing
ENDIF

! compute the associated transition wavenumbers
! ALLOCATE variable transition_wn

IF (ALLOCATED(transition_wn)) DEALLOCATE(transition_wn)
ALLOCATE(transition_wn(n_vibrations))

! ALLOCATE variable atot

IF (ALLOCATED(atot)) DEALLOCATE(atot)
ALLOCATE(atot(n_vibrations))

! Calculate the vector of tranistion wavenumbers
transition_wn = MATMUL(transition_matrix,om)

! Compute the dIFferent line intensities at temperature temp

DO i_transition = 1, n_transitions
atot(i_transition) = line_strength(i_transition)*t0ot*(1.0_EB-EXP(q2ot*transition_wn(i_transition)))/ & 
partition_FUNCTION(om,transition_matrix(i_transition,:),temp)
ENDDO

! compute the mean absorption coefficient using eq 11-44 from penner in quantitative molecular spectrometry, 1955
! sdweak = SUM(atot*propability) 
! just overlapping spectral line

DO i_transition = 1, n_transitions
sdweak = sdweak + atot(i_transition)*(-q2ot/(4._EB*be))*ABS(omega-bcnt(i_transition))* & 
                  EXP(q2ot/(4._EB*be)*(omega-bcnt(i_transition))**2)
ENDDO

gdinv  = gc1*dinv
gddinv = gd*dinv
!EXPress s/d at stp, as is in nasa sp-3080
sdweak = sdweak*temp/273._EB

!------------------------------------------------------------------------------
RETURN 
END SUBROUTINE Approx_Vib_Rot

!==============================================================================
ELEMENTAL SUBROUTINE co2(omega,temp,gc1,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser
IMPLICIT NONE

! variables passed in
REAL(EB), INTENT(IN) :: omega ! wavenumber, units in cm-1
REAL(EB), INTENT(IN) :: temp  ! temperature, units in kelvin
REAL(EB), INTENT(IN) :: gc1   ! half width at half height for lorentz line, units in cm-1

! variables passed out
REAL(EB), INTENT(OUT):: sdweak  ! narrow band mean ABS(ortion coefficient, units in cm-1
REAL(EB), INTENT(OUT):: gdinv   ! line width to line spacing ration for lorentz lines
REAL(EB), INTENT(OUT):: gddinv  ! line width to line spacing ration for DOppler
INTEGER,  INTENT(OUT):: i_model ! model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! local variables
INTEGER  :: i,j,k,l

REAL(EB), PARAMETER :: wm_co2 = 44._EB      ! molar mass (in gram) of co2
REAL(EB), PARAMETER :: be     = 0.391635_EB ! rotational constant co2
REAL(EB), PARAMETER :: q2     = 1.4388_EB   ! speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: t0     = 300._EB     ! reference temperature, in kelvin

REAL(EB), PARAMETER :: om1 = 1354.91_EB  ! fundamental frequency for transition (000)->(100), units cm-1
REAL(EB), PARAMETER :: om2 =  673.0_EB   ! fundamental frequency for transition (000)->(010), units cm-1
REAL(EB), PARAMETER :: om3 = 2396.49_EB  ! fundamental frequency for transition (000)->(001), units cm-1

REAL(EB), DIMENSION(3) :: atot, bcnt
REAL(EB), DIMENSION(:), ALLOCATABLE  :: om, line_strength
INTEGER, DIMENSION(:,:), ALLOCATABLE :: transition_matrix

REAL(EB) :: aa,bb,cc,qq,ee,ff,gg, &
         sminus,splus,sdstrg,gd, &
         x13,x23,x33,xbar,om12,alpha,omprim,v3,gam,  &
         omvv3,delta,v,omvbar,f1,f2,unflo1,unflo2,unflo3,test, &
         vbar1,oma,omb,ttemp,tt,t1,tw,ww,temp1,temp2,temp3,    &
         dinv,a,b,d,g,w1,dinv1,dinv2,dinv3,q2ot,t0ot,q2ot0

i_model = i_model_co2 

! initialize sdweak, gdinv, gddinv

sdweak = 0._EB
gdinv  = 1._EB
gddinv = 1._EB

IF (omega > 5725.0_EB) THEN
sdweak = 0.0_EB
gdinv  = 1.0_EB
gddinv = 1.0_EB
ELSE
  
! compute DOppler broadening half-widths sp-3080 eq: 5-35
gd = 5.94e-6_EB*omega*SQRT(temp/(273.0_EB*wm_co2))

IF(omega>4550.0_EB) THEN
! contribution to 2.0 micron band from:
! (000)-(041)
! (000)-(121),
! (000)-(201) transitions.

! band center 

   bcnt(1) = 4860.5_EB ! (000)-(041)
   bcnt(2) = 4983.5_EB ! (000)-(121)
   bcnt(3) = 5109.0_EB ! (000)-(201)
      
   IF(ALLOCATED(transition_matrix)) DEALLOCATE(transition_matrix)
   ALLOCATE(transition_matrix(3,3))

   IF(ALLOCATED(line_strength)) DEALLOCATE(line_strength)
   ALLOCATE(line_strength(3))

   IF(ALLOCATED(om)) DEALLOCATE(om)
   ALLOCATE(om(3))

   ! ! compute the associated transition wavenumbers
     
   transition_matrix(1,1:3) = [0,4,1] ! (000)-(041)
   transition_matrix(2,1:3) = [1,2,1] ! (000)-(121)
   transition_matrix(3,1:3) = [2,0,1] ! (000)-(201)

   om(1) = om1
   om(2) = om2
   om(3) = om3

!! compute the integrated intensity of the band
!! 0.272: integrated intensity (in cm-2-atm-1) of band bcnt(1) = 4860.5_EB ((041)-(000)) at t=300k
!! 1.01 : integrated intensity (in cm-2-atm-1) of band bcnt(2) = 4983.5_EB ((121)-(000)) at t=300k
!! 0.426: integrated intensity (in cm-2-atm-1) of band bcnt(3) = 5109  _EB ((201)-(000)) at t=300k
!! values from penner and varanasi jqsrt vol 4 pp 799-806, 1964
      
   line_strength = [0.272_EB,1.01_EB,0.426_EB]

   CALL approx_vib_rot(transition_matrix,line_strength,om,temp,omega,bcnt,be,gc1, &
                        gd,gdinv,gddinv,sdweak)
 
ELSEIF((omega<=4550._EB).and.(omega>3800._EB)) THEN
   sdweak = 0.0_EB
   gdinv  = 1.0_EB
   gddinv = 1.0_EB

ELSEIF((omega<=3800._EB).and.(omega>3050._EB)) THEN
   b   = 0.391635_EB
   a   = 0.0030875_EB

   x13 = -19.37_EB
   x23 = -12.53_EB
   x33 = -12.63_EB

   t0ot = t0/temp

   q2ot  = -q2/temp
   q2ot0 = -q2/t0

   xbar  = 0.5_EB*(.5_EB*x13+x23)
   om12  = 0.5_EB*(.5_EB*om1+om2)

   sdweak = 0._EB
   sdstrg = 0._EB

!calculate absorption coef. and line spacing PARAMETER for 2.7 micron band
   l = 1
!contribution to 2.7 micron band from (000)-(021) and (010)-(031) trans.
   alpha  = 28.5_EB
   omprim = 2.0_EB*om2+om3

   l120: DO
      aa = alpha*b*q2/(a*(1._EB-EXP(om3*q2ot0))*(1._EB-EXP(om12*q2ot0))**3*(1._EB+EXP(om12*q2ot0))*(1._EB-EXP(omprim*q2ot0)))
      bb = (1._EB-EXP(q2ot*omega))*(1._EB-EXP(q2ot*om3))* (1._EB-EXP(om12*q2ot))**3*(1._EB+EXP(om12*q2ot)) &
            *(1._EB-EXP(q2ot*omprim))
      cc = aa*bb*omega*t0/temp**2

      l102: DO j = 1, 20
         v = REAL(j-1,EB)
         IF ((j/2*2)==j) g = (v+1.0_EB)*(v+3.0_EB)/4.0_EB ! replace with modulo 2
         IF ((j/2*2)/=j) g = (v+2.0_EB)*(v+2.0_EB)/4.0_EB ! replace with modulo 2

         vbar1 = -1.0_EB+(v+3.0_EB)*(v+4.0_EB)/(v+2.0_EB)/6.0_EB

         IF(j/2*2==j) vbar1 = -1.0_EB+(v+5.0_EB)/6.0_EB ! replace with modulo 2

         l101: DO k = 1, 10
            v3  = REAL(k-1,EB)
            qq  = (v3+1)*g*EXP((v3*om3+v*om12)*q2ot)*(vbar1+1.0_EB)
            gam = b-a*(v3+1._EB)

            IF(l==2) THEN
               omvv3 = 3728.0_EB-5.0_EB*v-47._EB*v3
               IF (v<=TWO_EPSILON_EB) omvv3 = 3715._EB-47._EB*v3
            ELSE
               omvv3 = 3598.0_EB-18.0_EB*v-47.0_EB*v3
               IF (v<=TWO_EPSILON_EB) omvv3 = 3613.0_EB-47.0_EB*v3
            ENDIF
 
            delta = a*(omega-omvv3)

            IF ((gam**2)<=delta) CYCLE l102

            d      = 2.0_EB*SQRT(gam**2-delta)
            omvbar = omvv3*(1.0_EB-EXP(omvv3*q2ot))
            f1     = gam - 0.5_EB*d
            f2     = gam + 0.5_EB*d
            ee     = q2*gam/(a**2*temp)
            unflo1 = ee*delta*(1.0_EB + 0.5_EB*a/gam)

            IF(unflo1<=-78.0_EB) CYCLE l102

            unflo2 = ee*2.0_EB*gam*f1

            IF(unflo2>=78.0_EB) CYCLE l102

            ff     = EXP(ee*delta*(1.0_EB+0.5_EB*a/gam))
            sminus = cc*qq/omvbar*ABS(f1)*ff*EXP(-ee*2.0_EB*gam*f1)
            unflo3 = ee*2.0_EB*gam*f2

            IF(unflo3>=78.0_EB) THEN
               splus = 0.0_EB
            ELSE
               splus = cc*qq/omvbar*ABS(f2)*ff*EXP(-ee*2.0_EB*gam*f2)
            ENDIF

            gg     = sdweak
            sdweak = (sminus+splus)/d+sdweak
            test   = (sdweak-gg)/sdweak
            IF(test<.0001_EB) CYCLE l102
            sdstrg = SQRT(0.5_EB*g)*(SQRT(sminus)+SQRT(splus))/d+sdstrg
         ENDDO l101
      ENDDO l102

      IF(l==2) exit l120
!contribution to 2.7 micron band from (000)-(101) and (010)-(111) trans.
      alpha  = 42.3_EB
      omprim = om1+om3
      l      = 2
   ENDDO l120
!calculate absorption coef and line spacing PARAMETER for 4.3 micron band
   IF(sdweak<=TINY_EB) THEN
      sdweak = 0.0_EB
      gdinv  = 1.0_EB
      gddinv = 1.0_EB
   ELSE
      dinv   = sdstrg*sdstrg/sdweak
      gdinv  = gc1*dinv
      gddinv = gd*dinv
!***EXPress s/d at stp, as is k in nasa sp-3080
      sdweak = sdweak*temp/273.0_EB
!         ENDIF
   ENDIF

ELSEIF((omega<=3050.0_EB).and.(omega>2474.0_EB)) THEN
   sdweak = 0.0_EB
   gdinv  = 1.0_EB
   gddinv = 1.0_EB

ELSEIF((omega<=2474.0_EB).and.(omega>1975.0_EB)) THEN

b = .391635_EB
a = .0030875_EB

x13 = -19.37_EB
x23 = -12.53_EB
x33 = -12.63_EB

xbar = 0.5_EB*(0.5_EB*x13+x23)
om12 = 0.5_EB*(0.5_EB*om1+om2)

sdweak=0.0_EB
sdstrg=0.0_EB

! dipole approximation

IF(omega<=2395.0_EB) THEN
      alpha  = 2700.0_EB
      omprim = om3
      aa     = alpha*b*q2/(a*(1.0_EB-EXP(-om3*q2/t0))*(1.0_EB-EXP(-om12*q2 /t0))**3* & 
                              (1.0_EB+EXP(-om12*q2/t0))*(1.0_EB-EXP(-omprim*q2/t0)))
      bb     = (1._EB-EXP(-q2*omega/temp))*(1.0_EB-EXP(-q2*om3/temp))*(1.0_EB-EXP(-om12*q2/temp))**3*&
               (1._EB+EXP(-om12*q2/temp))*(1.0_EB-EXP(-q2*omprim/temp))
      cc     = aa*bb*omega*t0/temp**2

      l202a: DO j = 1, 20
         v  = REAL(j-1,EB)

         IF(j/2*2==j) g = 0.25_EB*(v+1.0_EB)*(v+3.0_EB) ! replace with modulo
         IF(j/2*2/=j) g = 0.25_EB*(v+2.0_EB)*(v+2.0_EB) ! replace with modulo

         l201a: DO k=1,10

            v3    = REAL(k-1,EB)
            qq    = (v3+1.0_EB)*g*EXP(-(v3*om3+v*om12)*q2/temp)
            gam   = b-a*(v3+1.0_EB)
            omvv3 = om3 + 0.5_EB*x13 + x23 + 2.0_EB*x33 + xbar*v + 2.0_EB*x33*v3
            delta = a*(omega-omvv3)

            IF(gam*gam<=delta) CYCLE l202a

            d      = 2.0_EB*SQRT(gam**2-delta)
            omvbar = omvv3*(1._EB-EXP(-omvv3*q2/temp))

            f1     = gam-0.5_EB*d
            f2     = gam+0.5_EB*d
            ee     = q2*gam/(a*a*temp)
            unflo1 = ee*delta*(1.0_EB + 0.5_EB*a/gam)

            IF(unflo1<=-78.0_EB) CYCLE l202a
            unflo2 = ee*2.0_EB*gam*f1

            IF(unflo2>=78.0_EB) CYCLE l202a

            ff     = EXP(ee*delta*(1.0_EB + 0.5_EB*a/gam))
            sminus = cc*qq/omvbar*ABS(f1)*ff*EXP(-ee*2.0_EB*gam*f1)
            unflo3 = ee*2.0_EB*gam*f2

            IF (unflo3 >= 78.0_EB) THEN
               splus=0.0_EB
            ELSE
               splus=cc*qq/omvbar*ABS(f2)*ff*EXP(-ee*2.0_EB*gam*f2)
            ENDIF

            gg     = sdweak
            sdweak = (sminus+splus)/d+sdweak
            test   = (sdweak-gg)/sdweak

            IF(test<1.0e-4_EB) CYCLE l202a

            sdstrg=SQRT(0.5_EB*g)*(SQRT(sminus)+SQRT(splus))/d+sdstrg

         ENDDO l201a
      ENDDO l202a

      IF(sdweak<=TINY_EB) THEN
         sdweak = 0.0_EB
         gdinv  = 1.0_EB
         gddinv = 1.0_EB
      ELSE
         dinv   = sdstrg**2/sdweak
         gdinv  = gc1*dinv
         gddinv = gd*dinv
!***  EXPress s/d at stp, as is k in nasa sp-3080
         sdweak = sdweak*temp/273.0_EB
      ENDIF
   ELSE
!calculate absorption coef. and line spacing PARAMETER for 2.7 micron band
      l=1
!contribution to 2.7 micron band from (000)-(021) and (010)-(031) trans.
      alpha  = 28.5_EB
      omprim = 2.0_EB*om2+om3

      l120a: DO
         aa = alpha*b*q2/(a*(1.0_EB-EXP(-om3*q2/t0))*(1.0_EB-EXP(-om12*q2 /t0))**3* & 
                              (1.0_EB+EXP(-om12*q2/t0))*(1.0_EB-EXP(-omprim*q2/t0)))
         bb = (1._EB-EXP(-q2*omega/temp))*(1.0_EB-EXP(-q2*om3/temp))*    & 
               (1._EB-EXP(-om12*q2/temp))**3*(1.0_EB+EXP(-om12*q2/temp))* &
               (1._EB-EXP(-q2*omprim/temp))

         cc = aa*bb*omega/temp*t0/temp

         l102a: DO j = 1, 20
            v=REAL(j-1,EB)

            IF(j/2*2==j) g = 0.25_EB*(v+1.0_EB)*(v+3.0_EB) ! USE modulo
            IF(j/2*2/=j) g = 0.25_EB*(v+2.0_EB)*(v+2.0_EB) ! USE modulo

            vbar1 = -1.0_EB + (v+3.0_EB)*(v+4.0_EB)/(6.0_EB*(v+2.0_EB))
            IF(j/2*2==j) vbar1 = -1.0_EB+(v+5.0_EB)/6.0_EB ! USE modulo

            l101a: DO k = 1, 10

               v3  = REAL(k-1,EB)
               qq  = (v3+1)*g*EXP( -(v3*om3 + v*om12)*q2/temp )*(vbar1+1.0_EB)
               gam = b-a*(v3+1.0_EB)

               IF(l==2) THEN
                  omvv3 = 3728.0_EB - 5.0_EB*v - 47.0_EB*v3
                  IF(v<=TWO_EPSILON_EB) omvv3 = 3715.0_EB - 47.0_EB*v3
               ELSE
                  omvv3 = 3598.0_EB - 18.0_EB*v - 47.0_EB*v3
                  IF(v<=TWO_EPSILON_EB) omvv3 = 3613.0_EB - 47.0_EB*v3
               ENDIF

               delta = a*(omega-omvv3)
               IF((gam**2)<=delta) CYCLE l102a

               d      = 2.0_EB*SQRT(gam**2-delta)
               omvbar = omvv3*(1.0_EB-EXP(-omvv3*q2/temp))
               f1     = gam-0.5_EB*d
               f2     = gam+0.5_EB*d

               ee     = q2*gam/(a**2*temp)
               unflo1 = ee*delta*(1.0_EB+0.5_EB*a/gam)

               IF(unflo1<= -78.0_EB) CYCLE l102a

               unflo2 = ee*2.0_EB*gam*f1

               IF(unflo2>= 78._EB) CYCLE l102a

               ff     = EXP(ee*delta*(1.0_EB+0.5_EB*a/gam))
               sminus = cc*qq/omvbar*ABS(f1)*ff*EXP(-ee*2.0_EB*gam*f1)
               unflo3 = ee*2.0_EB*gam*f2

               IF(unflo3>=78._EB) THEN
                  splus = 0._EB
               ELSE
                  splus = cc*qq/omvbar*ABS(f2)*ff*EXP(-ee*2.0_EB*gam*f2)
               ENDIF

               gg     = sdweak
               sdweak = (sminus+splus)/d+sdweak
               test   = (sdweak-gg)/sdweak

               IF(test<1.0e-4_EB) CYCLE l102a

               sdstrg = SQRT(0.5_EB*g)*(SQRT(sminus)+SQRT(splus))/d+sdstrg

            ENDDO l101a
         ENDDO l102a

         IF(l==2) exit l120a
!contribution to 2.7 micron band from (000)-(101) and (010)-(111) trans.
         alpha  = 42.3_EB
         omprim = om1+om3
         l      = 2
      ENDDO l120a
!calculate absorption coef and line spacing PARAMETER for 4.3 micron band
      IF(sdweak<=TINY_EB) THEN
         sdweak = 0.0_EB
         gdinv  = 1.0_EB
         gddinv = 1.0_EB
      ELSE
         dinv   = sdstrg**2/sdweak
         gdinv  = gc1*dinv
         gddinv = gd*dinv
!***EXPress s/d at stp, as is k in nasa sp-3080
         sdweak = sdweak*temp/273.0_EB
      ENDIF
   ENDIF

ELSEIF((omega<=1975.0_EB).and.(omega>1100.0_EB)) THEN
   sdweak = 0.0_EB
   gdinv  = 1.0_EB
   gddinv = 1.0_EB

! band 10 microns_
ELSEIF((omega<=1100.0_EB).and.(omega>880.0_EB)) THEN
!contribution to 10.0 micron band from (100)-(001) and (020)-(001) trans.
   bcnt(1) = 960.8_EB
   bcnt(2) = 1063.6_EB

   oma = om3
   omb = (om1+2.0_EB*om2)/2.0_EB

   atot(1)= 0.0219_EB
   atot(2)= 0.0532_EB

   DO k=1,2
      atot(k)=t0/temp*atot(k)*EXP(q2*omb*(1.0_EB/t0-1.0_EB/temp))*(1.0_EB-EXP(-q2*(oma-omb)/temp)) &
                  /( (1.0_EB-EXP(-q2*oma/temp))*(1.0_EB-EXP(-omb*q2/temp)) )
   ENDDO
   sdweak = 0.0_EB

   DO i = 1, 2
      sdweak = sdweak + atot(i)*q2/(4.0_EB*be*temp)*ABS(omega-bcnt(i))*EXP(-q2/(4.0_EB*be*temp)*(omega-bcnt(i))**2)
   ENDDO

   dinv   = 1.0_EB/(4.0_EB*be)
   gdinv  = gc1*dinv
   gddinv = gd*dinv
!***EXPress s/d at stp, as is in nasa sp-3080
   sdweak = sdweak*temp/273.0_EB

! tablulated band: 15 micron
   ELSEIF((omega<=880.0_EB).and.(omega>500.0_EB))  THEN
!contribution to 15.0 micron band from (000)-(010) trans.
   ttemp = temp
   j  = (floor(omega)-495)/5
   w1 = 495.0_EB + 5.0_EB*REAL(j,EB)
   ww = (omega-w1)/5

   IF(temp>=2400._EB) ttemp = 2399.99_EB
   IF(temp < 300._EB) ttemp =  300.00_EB

   i = ttemp/300._EB
 
   SELECT CASE(i)
      CASE(3)
         i=2
         tt=(ttemp-600._EB)/600._EB
      CASE(6:7) 
         i=5
         tt=(ttemp-1800._EB)/600._EB
      CASE DEFAULT
         t1=REAL(i,EB)*300._EB
         tt=(ttemp-t1)/300._EB
         IF (i>=4) i=i-1     
   END SELECT

   tw=tt*ww
   sdweak=sd15(i,j)*(1._EB-tt-ww+tw)+sd15(i+1,j)*(tt-tw)  +sd15(i,j+1)*(ww-tw)+sd15(i+1,j+1)*tw

   IF(sdweak<=TINY_EB) THEN
      sdweak = 0.0_EB
      gdinv  = 1.0_EB
      gddinv = 1.0_EB
   ELSE
!calculate line spacing PARAMETER for 15.0 micron band
      dinv1 =  1.2_EB
      dinv2 =  8.0_EB
      dinv3 = 30.0_EB

      temp1 = 300.0_EB
      temp2 = 550.0_EB
      temp3 = 830.0_EB

      dinv   = dinv1*(ttemp-temp2)*(ttemp-temp3)/(temp1-temp2)  /(temp1-temp3)+dinv2*(ttemp-temp1)*(ttemp-temp3)&
               /(temp2-temp1)/(temp2-temp3)+dinv3*(ttemp-temp1)  *(ttemp-temp2)/(temp3-temp1)/(temp3-temp2)
      gdinv  = gc1*dinv
      gddinv = gd*dinv
   ENDIF

ELSE
   sdweak = 0.0_EB
   gdinv  = 1.0_EB
   gddinv = 1.0_EB
ENDIF
ENDIF

!------------------------------------------------------------------------------
RETURN
END SUBROUTINE co2

!==============================================================================
ELEMENTAL SUBROUTINE h2o(omega,temp,gc2,sdweak,gdinv,gddinv,i_model)
!==============================================================================
!! this SUBROUTINE RETURNs the mean spectral absorption coefficient (sdweak),
!! the collision broadening fine structure PARAMETER (gdinv),
!! the DOppler broadening fine structure PARAMETER (gddinv), and the 
!! narrow band model to be used (i_model)
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser
!------------------------------------------------------------------------------

IMPLICIT NONE

! variables passed in

REAL(EB), INTENT(IN) :: omega ! wavenumber,  in cm-1
REAL(EB), INTENT(IN) :: temp  ! temperature, in kelvin
REAL(EB), INTENT(IN) :: gc2   ! half width at half height, in cm-1

! variables passed out

REAL(EB), INTENT(OUT) :: sdweak  ! spectral absorption coefficient, in cm-1
REAL(EB), INTENT(OUT) :: gddinv  ! line width to line spacing ration for DOppler
REAL(EB), INTENT(OUT) :: gdinv   ! line width to line spacing ratio for lorentz boradening
INTEGER,  INTENT(OUT) :: i_model ! model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! local variables

INTEGER  :: i,j
REAL(EB) :: w1, ww, t1, tt, tw, d, b, dinv, ttemp, gd
REAL(EB), PARAMETER :: wm_h2o = 18._EB

i_model = i_model_h2o

sdweak = 0._EB
gdinv  = 1._EB
gddinv = 1._EB

IF (omega>=50._EB.and.omega<9300) THEN
! compute DOpller half-width. eq 5-35
gd   = 5.94e-6_EB*omega*SQRT(temp/(273._EB*wm_h2o))

j    = (omega-25._EB)/25._EB
ttemp= temp

IF(temp>=2500._EB) ttemp = 2499.99_EB
IF(temp<300._EB)   ttemp = 300._EB

i = floor(ttemp/500._EB) + 1

IF (i==2.and.ttemp<600._EB) i = 1

w1 = 25._EB+25._EB*REAL(j,EB)
ww = (omega-w1)/25._EB

IF(i>2) THEN
   t1=REAL(i-1,EB)*500._EB
   tt=(ttemp-t1)/500._EB
ELSE
   IF(i==1) tt = (ttemp-300._EB)/300._EB
   IF(i==2) tt = (ttemp-600._EB)/400._EB
ENDIF

tw     = tt*ww
! perform interpolation
sdweak = sd(i,j)*(1._EB-tt-ww+tw)+sd(i+1,j)*(tt-tw)+sd(i,j+1) *(ww-tw)+sd(i+1,j+1)*tw

d      = -2.294_EB+.3004e-02_EB*ttemp-.366e-06_EB*temp**2
b      = dsin(.0036_EB*omega-8.043_EB)
dinv   = EXP(.7941_EB*b+d) ! from nasa ir handbook)

!     dinv=EXP(0.00106*temp-1.21)

! update argument out variables
gdinv  = gc2*dinv
gddinv = gd*dinv

ENDIF

!------------------------------------------------------------------------------
RETURN
END SUBROUTINE h2o


!==============================================================================
SUBROUTINE co(omega,temp,gc4,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser
INTEGER j
REAL(EB) omega,temp,gc4,sdweak,gdinv,gddinv,aa,bb,cc,qq,ee,ff,gg, &
      sminus,splus,sdstrg,b,alpha,a,ome,wx,wy,omprim,t0, &
      q2,wm,gd,v,gam,omv,delta,d,omvbar,f1,f2,test,dinv,q2ot,toaz

INTEGER, INTENT(OUT) :: i_model

i_model = i_model_co

IF(omega<1600._EB .or. omega>2400._EB) THEN
sdweak = 0._EB
gdinv  = 1._EB
gddinv = 1._EB
ELSE
b     = 1.93139_EB
alpha = 260._EB
a     = .017485_EB
ome   = 2170.21_EB
wx    = 13.461_EB
wy    = .0308_EB
omprim= ome-2._EB*wx+3.25_EB*wy

t0 = 300._EB
q2 = 1.4388_EB

toaz = temp/273._EB
q2ot = q2/temp

wm     = 28._EB
! DOppler broadening half-width
gd     = 5.94e-6_EB*omega*SQRT(toaz/wm)
sdweak = 1.e-99_EB
sdstrg = 1.e-99_EB

aa = alpha*b*q2/(a*(1._EB-EXP(-omprim*q2/t0))**2)
bb = (1._EB-EXP(-omega*q2ot))*(1._EB-EXP(-omprim*q2ot))**2
cc = aa*bb*omega*t0/temp**2

l101: DO j=1,20
   v    = REAL(j-1,EB)
   qq   = (v+1._EB)*EXP(-v*ome*q2ot)
   gam  = b-a*(v+1._EB)
   omv  = ome-2._EB*(v+1._EB)*wx+(3._EB*(v+1._EB)*(v+1._EB)+.25_EB)*wy
   delta= a*(omega-omv)

   IF(gam**2<=delta) exit l101

   d      = 2._EB*SQRT(gam*gam-delta)
   omvbar = omv*(1._EB-EXP(-omv*q2ot))
   f1     = gam-0.5_EB*d
   f2     = gam+0.5_EB*d
   ee     = q2*gam/(a*a*temp)
      
   ff     = EXP(ee*delta*(1._EB+.5_EB*a/gam))
   sminus = cc*qq/omvbar*ABS(f1)*ff*EXP(-ee*2._EB*gam*f1)
   splus  = cc*qq/omvbar*ABS(f2)*ff*EXP(-ee*2._EB*gam*f2)
   gg     = sdweak
   sdweak = (sminus+splus)/d+sdweak
   test   = (sdweak-gg)/sdweak
      
   IF(test<.0001_EB) exit l101

   sdstrg=(SQRT(sminus)+SQRT(splus))/d+sdstrg
ENDDO l101

dinv   = sdstrg*sdstrg/sdweak
gdinv  = gc4*dinv
gddinv = gd*dinv
!***EXPress s/d at stp, as is k in nasa sp-3080
sdweak = sdweak*toaz
ENDIF
!------------------------------------------------------------------------------
END SUBROUTINE co


!==============================================================================
!ELEMENTAL SUBROUTINE pod(omega,soot_volume_frac,path_length_cm,x_particle)
SUBROUTINE pod(omega,soot_volume_frac,path_length_cm,x_particle)
!==============================================================================
!! pod calculates particle optical depth, x_particle, of the volume fraction of 
!! soot particles in gas cloud.  rin and rik are
!! the REAL and imaginary parts of the index of refraction. the particles are 
!! assumed to be in the rayleigh limit.
! 
! argument in:
! omega: wavenumber, units: cm-1
! soot_volume_frac: soot volume fraction. no units
! path_length_cm: physical path length. units in cm
! temp: temperature in kelvin
!
! argument out:
! x_particle: opticle depth due to the presence of soot particle. no units.
!*-----------------------------------------------------------------------------

IMPLICIT NONE

! variables passed in:

REAL(EB), INTENT(IN) :: omega, path_length_cm, soot_volume_frac

! varibles passed out: x_particle: particle optical depth

REAL(EB), INTENT(OUT) :: x_particle

! locals

REAL(EB) :: abco,ff,lambda,rin,rik

lambda = 10000._EB/omega
 
rin = 1.6_EB
rik = 0.5_EB

!ff = 36.0_EB*pi*rin*rik/lambda/((rin**2-rik**2+2.0_EB)**2+(2.0_EB*rin*rik)**2)

! absorption coef. is based upon measurements of dalzell and 
! sarofim

ff         = 7.0_EB/lambda
abco       = ff*soot_volume_frac*1.e6_EB
x_particle = abco*path_length_cm*cm_to_m

!------------------------------------------------------------------------------
RETURN
END SUBROUTINE pod


!==============================================================================
SUBROUTINE ch4_old(omega,temp,pch4,ptot,gc3,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser
INTEGER i,j
REAL(EB) omega,temp,pch4,ptot,gc3,sdweak,gdinv,gddinv,be,q2, &
   wm,gd,om1,om2,om3,om4,com1,com2,com3,com4,dinv,pe,w1,sdb, &
   sda,sdc,q2ot,azot,toaz

INTEGER, INTENT(OUT) :: i_model

REAL(EB), DIMENSION(4) :: atot, bcnt

i_model = i_model_ch4

IF(omega>5000._EB .or. omega<1125._EB) THEN

sdweak=0.0_EB
gdinv=1._EB
gddinv=1._EB

ELSE

be=5.2412_EB
q2=1.4388_EB
wm=16._EB
q2ot = -q2/temp
azot = 273._EB/temp
toaz = temp/273._EB

gd=5.94e-6_EB*omega*SQRT(toaz)/4._EB

IF(omega>=3400._EB) THEN

! contribution to 2.4 micron band from (0000)-(0110), (0000)-(0011),
! (0000)-(1001), and (0000)-(0102) trans.  the integrated band intensities
! of vincent-geisse (annales de physique ser.12, v. 10, 1955) have
! been multiplied by a factor of 4 and the line spacing is that
! of v4 from gray and penner (jqsrt v. 5, 1965).

   om1=2914.2_EB
   om2=1526.0_EB
   om3=3020.3_EB
   om4=1306.2_EB

   bcnt(1) = 4123.0_EB
   bcnt(2) = 4216.3_EB
   bcnt(3) = 4313.2_EB
   bcnt(4) = 4546.0_EB

   com1=om2+2._EB*om4
   com2=om1+om4
   com3=om3+om4
   com4=om2+om3

   atot(1)=0.64_EB*azot*(1._EB-EXP(q2ot*com1))/((1._EB-EXP(q2ot*om2))*(1._EB-EXP(q2ot*om4))**2)
   atot(2)=17.6_EB*azot*(1._EB-EXP(q2ot*com2))/((1._EB-EXP(q2ot*om1))*(1._EB-EXP(q2ot*om4)))
   atot(3)=14.8_EB*azot*(1._EB-EXP(q2ot*com3))/((1._EB-EXP(q2ot*om3))*(1._EB-EXP(q2ot*om4)))
   atot(4)=5.04_EB*azot*(1._EB-EXP(q2ot*com4))/((1._EB-EXP(q2ot*om2))*(1._EB-EXP(q2ot*om3)))

   dinv=1._EB/5.74_EB
   gdinv=gc3*dinv
   gddinv=gd*dinv
   sdweak=0.0_EB

   DO i=1,4
      sdweak=sdweak+2._EB*(omega-bcnt(i))**2*(-q2ot*be)**1.5_EB*atot(i)/sqrtpi*dinv**3*EXP(q2ot*be*dinv**2 &
      *(omega-bcnt(i))**2)
   ENDDO
   sdweak=sdweak*toaz
ELSE
   pe=ptot+.3_EB*pch4

   IF(omega>=2625._EB) THEN
! contribution to 3.3 micron band from (0000)-(0010) trans.
! refer to brosmer and tien, jqsrt v. 33, p. 521

      gdinv  = .00734_EB*pe*SQRT(azot)*EXP(1.02_EB*(toaz-1._EB))
      gddinv = gd/9.4_EB

      j  = (omega-2600._EB)/25._EB
      w1 = 2600._EB+25._EB*REAL(j,EB)
      sdb= sd3(2,j)+(omega-w1)/25._EB*(sd3(2,j+1)-sd3(2,j))

      IF(temp>600._EB) THEN
         sdc=sd3(3,j)+(omega-w1)/25._EB*(sd3(3,j+1)-sd3(3,j))
         sdweak=sdb+(temp-600._EB)/250._EB*(sdc-sdb)
         IF(sdweak<0._EB)sdweak=0._EB
      ELSE
         sda=sd3(1,j)+(omega-w1)/25._EB*(sd3(1,j+1)-sd3(1,j))
         sdweak=sda+(temp-290._EB)/310._EB*(sdb-sda)
         IF(sdweak<0._EB)sdweak=0._EB
      ENDIF

   ELSEIF(omega>1450._EB) THEN

      sdweak=0.0_EB
      gdinv=1._EB
      gddinv=1._EB

   ELSE
! contribution to 7.7 micron band from (0000)-(0001) trans.
! refer to brosmer and tien, jqsrt v. 33, p. 521.
      gdinv  = .0243_EB*pe*(toaz)**.8_EB
      gddinv = gd/5.1_EB

      j   = (omega-1100._EB)/25._EB
      w1  = 1100._EB+25._EB*REAL(j,EB)
      sdb = sd7(2,j)+(omega-w1)/25._EB*(sd7(2,j+1)-sd7(2,j))

      IF(temp>600._EB) THEN
         sdc = sd7(3,j)+(omega-w1)/25._EB*(sd7(3,j+1)-sd7(3,j))
         sdweak = sdb+(temp-600._EB)/250._EB*(sdc-sdb)
         IF(sdweak<0._EB)sdweak=0._EB
      ELSE
         sda = sd7(1,j)+(omega-w1)/25._EB*(sd7(1,j+1)-sd7(1,j))
         sdweak = sda+(temp-290._EB)/310._EB*(sdb-sda)
         IF(sdweak<0._EB)sdweak=0._EB
      ENDIF
   ENDIF
ENDIF
ENDIF
!------------------------------------------------------------------------------
END SUBROUTINE ch4_old


!==============================================================================
SUBROUTINE ch4(omega,temp,pch4,ptot,gc3,sdweak,gdinv,gddinv,i_model)
!==============================================================================
!! compute methane optical properties using single line group
!!
!! band #1: 1150 cm-1 - 1600 cm-1 
!! band #2: 2700 cm-1 - 3250 cm-1 
!! band #3: 3400 cm-1 - 5000 cm-1
! 
! variables passed in
! omega: (REAL) wavenumber in cm-1
! temp : (REAL) temperature in k
! pch4 : (REAL) partial pressure of methane in atm
! ptot : (REAL) total pressure in atm
! gc3  : (REAL) collisional broadening half-width at half heigh
! 
! vraibles passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser
! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! be_ch4             : (REAL) rotational constant for ch4 [cm^-1]
! dinv_ch4           : inverse line spacing [cm]

REAL(EB), INTENT(IN)  :: omega, temp, pch4, ptot, gc3
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT)  :: i_model

REAL(EB) :: gd, pressure_effective, q2ot, azot, toaz, fact1

REAL(EB), PARAMETER :: q2     = 1.4388_EB ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_ch4 = 16.0425_EB

!line spacing d2_ch4=5.74 cm^-1 for the v4 fundamental of methane gray & penner,612
!footnote at bottom
REAL(EB), PARAMETER :: be_ch4 = 5.248_EB   ! rotational constants [cm^-1]
REAL(EB), DIMENSION(3), PARAMETER :: dinv_ch4 = (/1._EB/5.1_EB,1._EB/9.4_EB,1._EB/5.74_EB/)

REAL(EB), DIMENSION(4), PARAMETER :: &
               om_ch4  = (/2914.2_EB, 1526.0_EB, 3020.3_EB, 1306.2_EB/),                                               &
               com_ch4 = (/1526.0_EB+2._EB*1306.2_EB, 2914.2_EB+1306.2_EB, 3020.3_EB+1306.2_EB, 1526.0_EB+3020.3_EB/), &
               s2_ch4  = (/0.64_EB,17.6_EB,14.8_EB,5.04_EB/)

REAL(EB), DIMENSION(4) :: atot

INTEGER i

i_model = i_model_ch4

! initialize output

sdweak = 0.0_EB  !spectral absorption coefficient. 
gdinv  = 1._EB   !line width to line spacing ratio: fine structure PARAMETER
gddinv = 1._EB   !line width to line spacing ratio for DOppler boradening fine structure PARAMETER

! initialize some coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOpller half-width. eq 5-35
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_ch4) !DOppler half width [cm^-1]. nasa,222

! compute effective pressure. brosmer & tien,524 eq. 7 (specIFic to methane)
pressure_effective = ptot+.3_EB*pch4

!------------------------------------------------------------------------------
! computed properties 2.4 micron band, band #3: 3400 cm-1 - 5000 cm-1
IF((om_bnd_ch4(n_band_ch4,1)<=omega).and.(omega<=om_bnd_ch4(n_band_ch4,2))) THEN

   ! contribution of the v1+v4 combinaison band (2.4 micron)
   ! contribution to 2.4 micron band from
   ! (0000)-(0102)
   ! (0000)-(1001)
   ! (0000)-(0011)
   ! (0000)-(0110) transitions
   ! the integrated band intensities (s2_ch4)
   ! of vincent-geisse (annales de physique ser.12, v. 10, 1955) have
   ! been multiplied by a factor of 4 
   ! the line spacing (1/dinv_ch4) is that
   ! of v4 from gray and penner (jqsrt v. 5, 1965).                 

   atot(1) = s2_ch4(1)*azot*(1._EB-EXP(q2ot*com_ch4(1)))/((1._EB-EXP(q2ot*om_ch4(2)))*(1._EB-EXP(q2ot*om_ch4(4)))**2)
   atot(2) = s2_ch4(2)*azot*(1._EB-EXP(q2ot*com_ch4(2)))/((1._EB-EXP(q2ot*om_ch4(1)))*(1._EB-EXP(q2ot*om_ch4(4))))
   atot(3) = s2_ch4(3)*azot*(1._EB-EXP(q2ot*com_ch4(3)))/((1._EB-EXP(q2ot*om_ch4(3)))*(1._EB-EXP(q2ot*om_ch4(4))))
   atot(4) = s2_ch4(4)*azot*(1._EB-EXP(q2ot*com_ch4(4)))/((1._EB-EXP(q2ot*om_ch4(2)))*(1._EB-EXP(q2ot*om_ch4(3))))

   
   fact1  = q2ot*be_ch4*dinv_ch4(3)**2

   DO i=1,4
      sdweak = sdweak+(omega-com_ch4(i))**2*atot(i)*EXP(fact1*(omega-com_ch4(i))**2)
   ENDDO

   ! eq 6 in gray and penner jqsrt, vol 5, page 611-620, 1965
   sdweak=sdweak*2._EB*(-q2ot*be_ch4)**1.5_EB/sqrtpi*dinv_ch4(3)**3

   !***EXPress s/d at standard temperature and pressure, as is in nasa sp-3080
   sdweak = sdweak*toaz

   gdinv  = gc3*dinv_ch4(3)
   gddinv = gd*dinv_ch4(3) 

!------------------------------------------------------------------------------
! tabulated properties
ELSE IF((om_bnd_ch4(2,1)<=omega).and.(omega<om_bnd_ch4(2,2))) THEN
      
   !------------------------------------------------------------------------------
   ! contribution to 3.3 micron band, band #2: 2700 cm-1 - 3250 cm-1
   ! (0000)-(0010) transition (v3 fundamental)
   ! 9.4 is the average line position for the 3.3 micron band (cm^-1)
   ! source: brosmer & tien, jqsrt v. 33, p. 525 (specIFic to methane)

   sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch4_temp,om_bnd_ch4(2,:),sd2_ch4)
   ! line shape PARAMETER brosmer & tien, jqsrt v. 33, p. 525 eq. 10 (specIFic to methane)
   gdinv  = .00734_EB*pressure_effective*SQRT(azot)*EXP(1.02_EB*(toaz-1._EB)) 
   gddinv = gd*dinv_ch4(2)
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
   
ELSE IF((om_bnd_ch4(1,1)<=omega).and.(omega<om_bnd_ch4(1,2))) THEN
   !------------------------------------------------------------------------------
   ! contribution to 7.7 micron band, band #1: 1150 cm-1 - 1600 cm-1 
   ! (0000)-(0001) transition (v4 fundamental)
   ! 5.1 (cm-1) is the average line position for the 7.7 micron band 
   ! source: brosmer & tien, jqsrt v. 33, p. 525 (specIFic to methane)

   sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch4_temp,om_bnd_ch4(1,:),sd1_ch4)
   ! line shape PARAMETER brosmer & tien, jqsrt v. 33, 525 eq. 11 (specIFic to methane)
   gdinv  = .0243_EB*pressure_effective*(toaz)**.8_EB
   gddinv = gd*dinv_ch4(1)
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF 

!------------------------------------------------------------------------------
END SUBROUTINE ch4

!==============================================================================
SUBROUTINE c3h6(omega,temp,pc3h6,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
!! compute propylene optical properties using single line group
!! 
!! absorption bands: 
!! there are 3 bands for propylene
!!
!! band #1: 1250 cm-1 - 1950 cm-1 
!! band #2: 2700 cm-1 - 3200 cm-1 
!! band #3: 700 cm-1 - 1150 cm-1 
!
! variables passed in
! omega  : (REAL) wavenumber in cm-1
! temp   : (REAL) temperature in k
! pc3h6 : (REAL) partial pressure of propylene
! ptot   : (REAL) total pressure
! gc3    : (REAL) collisional broadening half-width at half heigh
! 
! vraibles passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_c3h6         : inverse line spacing [cm] for propylene

REAL(EB), INTENT(IN)  :: omega, temp, pc3h6, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_c3h6 

REAL(EB), PARAMETER :: q2       = 1.4388_EB  ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_c3h6  = 42.0797_EB ! nist wEBbook data

INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_c3h6
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_c3h6) !DOppler half width [cm^-1]. nasa,222

! compute average line spacing
! model propylene as a prolate symmetric top. the line spacing is approximated to
! 2*rotational_constant b. reCALL that rotational_constant b = rotational_constant c
! rotational_constant_b obtained from cccbdv wEBsite, using EXPerimental values  
! propylene belongs to group cs
dinv_c3h6 = 1._EB/(2._EB*0.312_EB)

! IF omega is within absorption band, THEN compute sdweak from tabulated data
! propylene has three bands
! band #1: 775 cm-1 - 1150 cm-1 
! band #2: 1225 cm-1 - 1975 cm-1 
! band #3: 2650 cm-1 - 3275 cm-1  

! get the index of the band
band: DO i = 1, n_band_c3h6
   IF ((om_bnd_c3h6(i,1)<=omega).and.(omega<om_bnd_c3h6(i,2))) THEN
      in_band = .true.
      i_band  = i      
      exit band
   ENDIF
ENDDO band

IF (in_band) THEN
   
   pressure_effective = ptot+(be_c3h6(i_band)-1.0_EB)*pc3h6 
   
   SELECT CASE (i_band)
      CASE (1) ! band #1: 775 cm-1 - 1150 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h6_temp,om_bnd_c3h6(i_band,:) ,sd1_c3h6)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h6_temp,om_bnd_c3h6(i_band,:) ,gammad1_c3h6)
      CASE (2) ! band #2: 1225 cm-1 - 1975 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h6_temp,om_bnd_c3h6(i_band,:) ,sd2_c3h6)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h6_temp,om_bnd_c3h6(i_band,:) ,gammad2_c3h6)
      CASE (3) ! band #3: 2650 cm-1 - 3275 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h6_temp,om_bnd_c3h6(i_band,:) ,sd3_c3h6)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h6_temp,om_bnd_c3h6(i_band,:) ,gammad3_c3h6)
   END SELECT

   ! compute properties
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_c3h6
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE c3h6

!==============================================================================
SUBROUTINE c3h8(omega,temp,pc3h8,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! compute propane optical properties using single line group
! 
! absorption bands: 
! there are 2 bands for propane

! band #1: 1175 cm-1 - 1675 cm-1 
! band #2: 2550 cm-1 - 3375 cm-1 
!
! variables passed in
! omega  : (REAL) wavenumber in cm-1
! temp   : (REAL) temperature in k
! pc3h8 : (REAL) partial pressure of propane
! ptot   : (REAL) total pressure
! gc3    : (REAL) collisional broadening half-width at half heigh
! 
! vraibles passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_c3h8         : inverse line spacing [cm] for propane

REAL(EB), INTENT(IN)  :: omega, temp, pc3h8, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_c3h8 

REAL(EB), PARAMETER :: q2       = 1.4388_EB  ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_c3h8  = 44.0956_EB ! nist wEBbook data

INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_c3h8
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_c3h8) !DOppler half width [cm^-1]. nasa,222

! compute average line spacing
! propane is a prolate symmetric top. the line spacing is approximated to
! 2*rotational_constant b. reCALL that rotational_constant b = rotational_constant c
! rotational_constant_b obtained from cccbdv wEBsite, using EXPerimental values  
! propane belongs to group c2v
dinv_c3h8 = 1._EB/(2._EB*0.28173_EB)

! IF omega is within absorption band, THEN compute sdweak from tabulated data
! propane has three bands
! band #1: 1175 cm-1 - 1675 cm-1 
! band #2: 2550 cm-1 - 3375 cm-1 

! get the index of the band
band: DO i = 1, n_band_c3h8
   IF ((om_bnd_c3h8(i,1)<=omega).and.(omega<om_bnd_c3h8(i,2))) THEN
      in_band = .true.
      i_band  = i
      exit band
   ENDIF
ENDDO band

IF (in_band) THEN

   pressure_effective = ptot+(be_c3h8(i_band)-1.0_EB)*pc3h8 

   SELECT CASE (i_band)
      CASE (1) ! band #1: 1175 cm-1 - 1675 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h8_temp,om_bnd_c3h8(i_band,:) ,sd1_c3h8)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h8_temp,om_bnd_c3h8(i_band,:) ,gammad1_c3h8)
      CASE (2) ! band #2: 2550 cm-1 - 3375 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h8_temp,om_bnd_c3h8(i_band,:) ,sd2_c3h8)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c3h8_temp,om_bnd_c3h8(i_band,:) ,gammad2_c3h8)
   END SELECT

   ! compute properties  
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_c3h8
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE c3h8

!==============================================================================
SUBROUTINE c7h16(omega,temp,pc7h16,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! compute heptane optical properties using single line group
! 
! absorption bands: 
! there are 2 bands for heptane

! band #1: 1100 cm-1 - 1800 cm-1 
! band #2: 2550 cm-1 - 3275 cm-1 
!
! variables passed in
! omega  : (REAL) wavenumber in cm-1
! temp   : (REAL) temperature in k
! pc7h16 : (REAL) partial pressure of heptane
! ptot   : (REAL) total pressure
! gc3    : (REAL) collisional broadening half-width at half heigh
! 
! vraibles passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_c7h16         : inverse line spacing [cm] for heptane

REAL(EB), INTENT(IN)  :: omega, temp, pc7h16, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_c7h16 

REAL(EB), PARAMETER :: q2       = 1.4388_EB  ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_c7h16 = 100.2019_EB ! nist wEBbook data

INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_c7h16
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_c7h16) !DOppler half width [cm^-1]. nasa,222

! compute average line spacing
! n_heptane is a prolate symmetric top. the line spacing is approximated to
! 2*rotational_constant b. reCALL that rotational_constant b = rotational_constant c
! rotational_constant_b obtained from cccbdv wEBsite, using EXPerimental values  
! n-heptane belongs to group c2v
dinv_c7h16=1._EB/(2.0_EB*0.024_EB)

! IF omega is within absorption band, THEN compute sdweak from tabulated data
! heptane has three bands
! band #1: 1100 cm-1 - 1800 cm-1 
! band #2: 2550 cm-1 - 3275 cm-1 

! get the index of the band
band: DO i = 1, n_band_c7h16
   IF ((om_bnd_c7h16(i,1)<=omega).and.(omega<om_bnd_c7h16(i,2))) THEN
      in_band = .true.
      i_band  = i
      exit band
   ENDIF
ENDDO band

IF (in_band) THEN

   pressure_effective = ptot+(be_c7h16(i_band)-1.0_EB)*pc7h16 

   SELECT CASE (i_band)
      CASE (1) ! band #1: 1175 cm-1 - 1675 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h16_temp,om_bnd_c7h16(i_band,:),sd1_c7h16)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h16_temp,om_bnd_c7h16(i_band,:),gammad1_c7h16)
      CASE (2) ! band #2: 2550 cm-1 - 3375 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h16_temp,om_bnd_c7h16(i_band,:),sd2_c7h16)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h16_temp,om_bnd_c7h16(i_band,:),gammad2_c7h16)
   END SELECT

   ! compute properties  
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_c7h16
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE c7h16

!==============================================================================
SUBROUTINE c7h8(omega,temp,pc7h8,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! compute toluene optical properties using single line group
! 
! absorption bands: 
! there are 5 bands for toluene

! band #1: 700 cm-1 - 825 cm-1 
! band #2: 975 cm-1 - 1175 cm-1 
! band #3: 1275 cm-1 - 1675 cm-1 
! band #4: 1650 cm-1 - 2075 cm-1 
! band #5: 2675 cm-1 - 3225 cm-1 
!
! variables passed in
! omega  : (REAL) wavenumber in cm-1
! temp   : (REAL) temperature in k
! pc7h8 : (REAL) partial pressure of toluene
! ptot   : (REAL) total pressure
! gc3    : (REAL) collisional broadening half-width at half heigh
! 
! vraibles passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_c7h8         : inverse line spacing [cm] for toluene

REAL(EB), INTENT(IN)  :: omega, temp, pc7h8, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_c7h8 

REAL(EB), PARAMETER :: q2       = 1.4388_EB  ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_c7h8  = 92.1384_EB ! nist wEBbook data

INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_c7h8
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_c7h8) !DOppler half width [cm^-1]. nasa,222

! toluene belongs to group cs 
! toluene is model as a prolate symmetric top. the line spacing is approximated to
! 2*rotational_constant b. reCALL that rotational_constant b = rotational_constant c
! rotational_constant_b obtained from cccbdv wEBsite.  using EXPerimental values  

dinv_c7h8=1._EB/(2.0_EB*0.084_EB)


! IF omega is within absorption band, THEN compute sdweak from tabulated data
! toluene has 5 bands
! band #1: 700 cm-1 - 825 cm-1 
! band #2: 975 cm-1 - 1175 cm-1 
! band #3: 1275 cm-1 - 1675 cm-1 
! band #4: 1650 cm-1 - 2075 cm-1 
! band #5: 2675 cm-1 - 3225 cm-1 

! get the index of the band
band: DO i = 1, n_band_c7h8
   IF ((om_bnd_c7h8(i,1)<=omega).and.(omega<om_bnd_c7h8(i,2))) THEN
      in_band = .true.
      i_band  = i
      exit band
   ENDIF
ENDDO band


IF (in_band) THEN

   pressure_effective = ptot+(be_c7h8(i_band)-1.0_EB)*pc7h8 

   SELECT CASE (i_band)
      CASE (1) ! band #1: 700 cm-1 - 825 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),sd1_c7h8)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),gammad1_c7h8)
      CASE (2) ! band #2: 975 cm-1 - 1175 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),sd2_c7h8)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),gammad2_c7h8)
      CASE (3) ! band #3: 1275 cm-1 - 1675 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),sd3_c7h8)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),gammad3_c7h8)
      CASE (4) ! band #4: 1650 cm-1 - 2075 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),sd4_c7h8)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),gammad4_c7h8)
      CASE (5) ! band #5: 2675 cm-1 - 3225 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),sd5_c7h8)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c7h8_temp,om_bnd_c7h8(i_band,:),gammad5_c7h8)
   END SELECT

   ! compute properties  
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_c7h8
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE c7h8


!==============================================================================
SUBROUTINE ch3oh(omega,temp,pch3oh,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! compute methanol optical properties using single line group
! 
! absorption bands: 
! there are 4 bands for methanol

! band #1: 825 cm-1 - 1150 cm-1 
! band #2: 1125 cm-1 - 1700 cm-1 
! band #3: 2600 cm-1 - 3225 cm-1 
! band #4: 3525 cm-1 - 3850 cm-1
!
! variables passed in
! omega  : (REAL) wavenumber in cm-1
! temp   : (REAL) temperature in k
! pch3oh : (REAL) partial pressure of methanol
! ptot   : (REAL) total pressure
! 
! variables passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_ch3oh         : inverse line spacing [cm] for methanol

REAL(EB), INTENT(IN)  :: omega, temp, pch3oh, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_ch3oh 

REAL(EB), PARAMETER :: q2       = 1.4388_EB  ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_ch3oh = 32.0419_EB ! nist wEBbook data
INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_ch3oh
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_ch3oh) !DOppler half width [cm^-1]. nasa,222

! compute average line spacing
! methanol belongs to point group cs. for now the line spacing is approximated to
! 2*rotational_constant b since rotational_constant b ~ rotational_constant c
! rotational_constant_b obtained from cccbdv wEBsite, using EXPerimental values  
dinv_ch3oh=1._EB/(2.0_EB*0.82338_EB)

! IF omega is within absorption band, THEN compute sdweak from tabulated data
! methanol has 4 bands
! band #1: 825 cm-1 - 1150 cm-1 
! band #2: 1125 cm-1 - 1700 cm-1 
! band #3: 2600 cm-1 - 3225 cm-1 
! band #4: 3525 cm-1 - 3850 cm-1 

! get the index of the band
band: DO i = 1, n_band_ch3oh
   IF ((om_bnd_ch3oh(i,1)<=omega).and.(omega<om_bnd_ch3oh(i,2))) THEN
      in_band = .true.
      i_band  = i
      exit band
   ENDIF
ENDDO band


IF (in_band) THEN

   pressure_effective = ptot+(be_ch3oh(i_band)-1.0_EB)*pch3oh 

   SELECT CASE (i_band)
      CASE (1) ! band #1: 825 cm-1 - 1150 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,sd1_ch3oh)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,gammad1_ch3oh)
      CASE (2) ! band #2: 1125 cm-1 - 1700 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,sd2_ch3oh)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,gammad2_ch3oh)
      CASE (3) ! band #3: 2600 cm-1 - 3225 cm-1  
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,sd3_ch3oh)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,gammad3_ch3oh)
      CASE (4) ! band #4: 3525 cm-1 - 3850 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,sd4_ch3oh)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_ch3oh_temp,om_bnd_ch3oh(i_band,:) ,gammad4_ch3oh)
   END SELECT

   ! compute properties  
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_ch3oh
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE ch3oh

!==============================================================================
SUBROUTINE c5h8o2(omega,temp,pc5h8o2,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! compute mma optical properties using single line group
! 
! absorption bands: 
! there are 6 bands for mma
! band #1: 750 cm-1 - 900 cm-1 
! band #2: 875 cm-1 - 1075 cm-1 
! band #3: 1050 cm-1 - 1275 cm-1 
! band #4: 1250 cm-1 - 1575 cm-1 
! band #5: 1550 cm-1 - 1975 cm-1 
! band #6: 2650 cm-1 - 3275 cm-1
!
! variables passed in
! omega   : (REAL) wavenumber in cm-1
! temp    : (REAL) temperature in k
! pc5h8o2 : (REAL) partial pressure of mma
! ptot    : (REAL) total pressure
! 
! variables passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_c5h8o2         : inverse line spacing [cm] for mma

REAL(EB), INTENT(IN)  :: omega, temp, pc5h8o2, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_c5h8o2 

REAL(EB), PARAMETER :: q2        = 1.4388_EB   ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_c5h8o2 = 100.1158_EB ! nist wEBbook data
INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_c5h8o2
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_c5h8o2) !DOppler half width [cm^-1]. nasa,222

! compute average line spacing
! mma belongs to point group ? for now using data from acetylacetone (c5h8o2)
dinv_c5h8o2=1._EB/(2.0_EB*0.067_EB)

! IF omega is within absorption band, THEN compute sdweak from tabulated data
! there are 6 bands for mma
! band #1: 750 cm-1 - 900 cm-1 
! band #2: 875 cm-1 - 1075 cm-1 
! band #3: 1050 cm-1 - 1275 cm-1 
! band #4: 1250 cm-1 - 1575 cm-1 
! band #5: 1550 cm-1 - 1975 cm-1 
! band #6: 2650 cm-1 - 3275 cm-1 

! get the index of the band
band: DO i = 1, n_band_c5h8o2
   IF ((om_bnd_c5h8o2(i,1)<=omega).and.(omega<om_bnd_c5h8o2(i,2))) THEN
      in_band = .true.
      i_band  = i
      exit band
   ENDIF
ENDDO band

IF (in_band) THEN

   pressure_effective = ptot+(be_c5h8o2(i_band)-1.0_EB)*pc5h8o2 

   SELECT CASE (i_band)
      CASE (1) ! band #1: 750 cm-1 - 900 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),sd1_c5h8o2)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),gammad1_c5h8o2)
      CASE (2) ! band #2: 875 cm-1 - 1075 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),sd2_c5h8o2)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),gammad2_c5h8o2)
      CASE (3) ! band #3: 1050 cm-1 - 1275 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),sd3_c5h8o2)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),gammad3_c5h8o2)
      CASE (4) ! band #4: 1250 cm-1 - 1575 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),sd4_c5h8o2)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),gammad4_c5h8o2)
      CASE (5) ! band #5: 1550 cm-1 - 1975 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),sd5_c5h8o2)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),gammad5_c5h8o2)
      CASE (6) ! band #6: 2650 cm-1 - 3275 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),sd6_c5h8o2)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c5h8o2_temp,om_bnd_c5h8o2(i_band,:),gammad6_c5h8o2)
   END SELECT

   ! compute properties  
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_c5h8o2
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE c5h8o2


!==============================================================================
SUBROUTINE c2h6(omega,temp,pc2h6,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! compute ethane optical properties using single line group
! 
! absorption bands: 
! there are 3 bands for ethane
! band #1: 750 cm-1 - 1100 cm-1 
! band #2: 1250 cm-1 - 1700 cm-1 
! band #3: 2550 cm-1 - 3375 cm-1 

!
! variables passed in
! omega   : (REAL) wavenumber in cm-1
! temp    : (REAL) temperature in k
! pc2h6 : (REAL) partial pressure of ethane
! ptot    : (REAL) total pressure
! 
! variables passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_c2h6         : inverse line spacing [cm] for ethane

REAL(EB), INTENT(IN)  :: omega, temp, pc2h6, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_c2h6 

REAL(EB), PARAMETER :: q2      = 1.4388_EB  ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_c2h6 = 30.0690_EB ! nist wEBbook data
INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_c2h6
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_c2h6) !DOppler half width [cm^-1]. nasa,222

! compute average line spacing
! ethane belongs to point group d3d with rotational constant a = 2.51967 cm-1
! b = 0.68341 cm-1, and c = 0.68341 cm-1.
dinv_c2h6=1._EB/(2.0_EB*0.68341_EB)

! IF omega is within absorption band, THEN compute sdweak from tabulated data
! there are 3 bands for ethane
! band #1: 750 cm-1 - 1100 cm-1 
! band #2: 1250 cm-1 - 1700 cm-1 
! band #3: 2550 cm-1 - 3375 cm-1 

! get the index of the band
band: DO i = 1, n_band_c2h6
   IF ((om_bnd_c2h6(i,1)<=omega).and.(omega<om_bnd_c2h6(i,2))) THEN
      in_band = .true.
      i_band  = i
      exit band
   ENDIF
ENDDO band

IF (in_band) THEN

   pressure_effective = ptot+(be_c2h6(i_band)-1.0_EB)*pc2h6 

   SELECT CASE (i_band)
      CASE (1)  ! band #1:  750 cm-1 - 1100 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h6_temp,om_bnd_c2h6(i_band,:),sd1_c2h6)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h6_temp,om_bnd_c2h6(i_band,:),gammad1_c2h6)
      CASE (2) ! band #2: 1250 cm-1 - 1700 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h6_temp,om_bnd_c2h6(i_band,:),sd2_c2h6)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h6_temp,om_bnd_c2h6(i_band,:),gammad2_c2h6)
      CASE (3) ! band #3: 2550 cm-1 - 3375 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h6_temp,om_bnd_c2h6(i_band,:),sd3_c2h6)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h6_temp,om_bnd_c2h6(i_band,:),gammad3_c2h6)
   END SELECT

   ! compute properties  
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_c2h6
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE c2h6


!==============================================================================
SUBROUTINE c2h4(omega,temp,pc2h4,ptot,sdweak,gdinv,gddinv,i_model)
!==============================================================================
! compute ethylene optical properties using single line group
! 
! absorption bands: 
! there are 4 bands for ethylene
! band #1: 800 cm-1  - 1250 cm-1 
! band #2: 1325 cm-1 - 1575 cm-1 
! band #3: 1775 cm-1 - 2050 cm-1 
! band #4: 2825 cm-1 - 3350 cm-1 
!
! variables passed in
! omega   : (REAL) wavenumber in cm-1
! temp    : (REAL) temperature in k
! pc2h4 : (REAL) partial pressure of ethylene
! ptot    : (REAL) total pressure
! 
! variables passed out
! sdweak : (REAL) spectral absorption coefficient
! gdinv  : (REAL) line width to line spacing ratio, fine structure PARAMETER
! gddinv : (REAL) line width to line spacing ratio for DOppler boradening,
!                 fine structure PARAMETER
! i_model: (INTEGER) model used for snb fitting: 1: goody, 2 malkmus, 3 elsasser

! locals
! pressure_effective : (REAL) effective pressure, (formely pe)
! dinv_c2h4          : inverse line spacing [cm] for ethylene

REAL(EB), INTENT(IN)  :: omega, temp, pc2h4, ptot
REAL(EB), INTENT(OUT) :: sdweak, gdinv, gddinv

INTEGER, INTENT(OUT) :: i_model

REAL(EB) :: q2ot, azot, toaz
REAL(EB) :: gd, pressure_effective
REAL(EB) :: dinv_c2h4 

REAL(EB), PARAMETER :: q2      = 1.4388_EB  ! q2 = speed_of_light*planck_cns/boltzmann
REAL(EB), PARAMETER :: wm_c2h4 = 28.0532_EB ! nist wEBbook data
INTEGER :: i_band, i

LOGICAL :: in_band ! true IF omega is within some tabulated band, false otherwise

! get model for snb fitting

i_model = i_model_c2h4
i_band  = 0
in_band = .false.

! set initial values to sdweak, gdinv, gddinv - values are RETURNed IF omega is not
! within absorption bands
sdweak = 0.0_EB
gdinv  = 1._EB 
gddinv = 1._EB 

! compute thermal coefficients
q2ot = -q2/temp       
azot = 273._EB/temp
toaz = temp/273._EB

! compute DOppler broadening half width gd
gd   = 5.94e-6_EB*omega*SQRT(toaz/wm_c2h4) !DOppler half width [cm^-1]. nasa,222

! compute average line spacing
! ethylene belongs to point group d2h with rotational constant a = 4.828 cm-1
! b = 1.0012 cm-1, and c = 0.8282 cm-1. we make the assumption b ~ c
dinv_c2h4=1._EB/(2.0_EB*1.0012_EB)

! IF omega is within absorption band, THEN compute sdweak from tabulated data
! there are 4 bands for ethylene
! band #1: 800 cm-1  - 1250 cm-1 
! band #2: 1325 cm-1 - 1575 cm-1 
! band #3: 1775 cm-1 - 2050 cm-1 
! band #4: 2825 cm-1 - 3350 cm-1 

! get the index of the band
band: DO i = 1, n_band_c2h4
   IF ((om_bnd_c2h4(i,1)<=omega).and.(omega<om_bnd_c2h4(i,2))) THEN
      in_band = .true.
      i_band  = i
      exit band
   ENDIF
ENDDO band

IF (in_band) THEN

   pressure_effective = ptot+(be_c2h4(i_band)-1.0_EB)*pc2h4 

   SELECT CASE (i_band)
      CASE (1)  ! band #1: 800 cm-1  - 1250 cm-1 
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),sd1_c2h4)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),gammad1_c2h4)
      CASE (2) ! band #2: 1325 cm-1 - 1575 cm-1  
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),sd2_c2h4)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),gammad2_c2h4)
      CASE (3) ! band #3: 1775 cm-1 - 2050 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),sd3_c2h4)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),gammad3_c2h4)
      CASE (4) ! band #4: 2825 cm-1 - 3350 cm-1
         sdweak = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),sd4_c2h4)
         gdinv  = GET_SPECTRAL_ABSORPTION(omega,temp,sd_c2h4_temp,om_bnd_c2h4(i_band,:),gammad4_c2h4)
   END SELECT

   ! compute properties  
   gdinv  = gdinv*pressure_effective
   gddinv = gd*dinv_c2h4
   !***EXPress s/d at stp, as is k in nasa sp-3080
   sdweak = sdweak*toaz
ENDIF

RETURN
!------------------------------------------------------------------------------
END SUBROUTINE c2h4

!==============================================================================
pure REAL(EB) FUNCTION GET_SPECTRAL_ABSORPTION(omega,temp,sd_temp,bounds,absorption_data)
!==============================================================================
! RETURNs interpolated values of spectral absorption
! note assumed shape
!
! input:
! omega  : wavenumber, unit = 1/cm
! temp   : gas temperature, unit = kelvin
! sd_temp: vector of temperature measurements, DIMENSION(*) unit = kelvin
! bounds : vector of wavenumber bounds and spacing for a given band, DIMENSIONs(3)
! absorption_data: array of species spectral measurements, DIMENSIONs(*,*)
!                 note = first DIMENSION equals to DIMENSION of sd_temp
!
! action performed: bilinear interpolation of absorption_data
!------------------------------------------------------------------------------
! variables passed in

REAL(EB), INTENT(IN) :: omega
REAL(EB), INTENT(IN) :: temp
REAL(EB), DIMENSION(:), INTENT(IN)   :: sd_temp, bounds
REAL(EB), DIMENSION(:,:), INTENT(IN) :: absorption_data 

REAL(EB) :: delta_t
REAL(EB) :: delta_w
REAL(EB) :: w1
REAL(EB) :: omega_min, omega_max, delta_omega
REAL(EB) :: ttemp
REAL(EB) :: interpolated_value

INTEGER :: n_temp, n_bounds, n_abs_data1, n_abs_data2
INTEGER :: i_omega, i_temp, i_omega1, i_omega2

LOGICAL :: cross

cross = .false.

! cvl need to link MODULE that defines errmsg
! part of the code commented for now
!  CHARACTER(255) messag
!  LOGICAL fatal

n_temp      = ubound(sd_temp,1)  ! get number of temperature measurements
n_bounds    = ubound(bounds,1)   ! get number of elements in bounds. should be 3
n_abs_data1 = ubound(absorption_data,1) ! get number of temperature data
n_abs_data2 = ubound(absorption_data,2) ! get number of wavenumber data

! ! cvl need to link MODULE that defines errmsg
! ! part of the code commented for now
! !-----------------------------test DIMENSIONs of bounds and absorption_data----
! ! n_bounds should be 3
! 
!  IF (n_bounds/=3) THEN 
!     messag='error in radcal. n_bounds should be equal to 3'
!     fatal = .true.
!     CALL errmsg( TRIM(messag), fatal )
!  ENDIF
! 
! ! n_abs_data1 should be equal to n_temp
! 
!  IF (n_abs_data1/=n_temp) THEN 
!     messag='error in radcal. first DIMENSION of absorption_data should be equal to n_temp'
!     fatal = .true.
!     CALL errmsg( TRIM(messag), fatal )
!  ENDIF

!------------------------------------------------------------------------------
omega_min   = bounds(1)
omega_max   = bounds(2)
delta_omega = bounds(3)
 
! treat special CASE when band crosses the 1100 cm-1 mark.
! changes of stride
! below 1100 cm-1, delta_omega = 5 cm-1
! above 1100 cm-1, delta_omega = 25 cm-1

IF ((omega_min<=1100.0_EB).and.(omega_max>1100.0_EB)) THEN
cross = .true.
IF (omega>1100._EB) delta_omega = 25.0_EB
ENDIF

! bounds ttemp so that MIN(sd_temp) <= ttemp <= MAX(sd_temp)
! IF ttemp < MIN(sd_temp) THEN ttemp := MIN(sd_temp)
! IF MAX(sd_temp) < ttemp THEN ttemp := MAX(sd_temp)
ttemp       = MIN(maxval(sd_temp),MAX(temp,minval(sd_temp)))  

IF (ttemp <= minval(sd_temp)) THEN
   i_temp  = 1
   delta_t = 0._EB
ELSE IF (ttemp >= maxval(sd_temp)) THEN
   i_temp  = n_temp-1
   delta_t = 1.0_EB
ELSE
! find index i_temp such that sd_temp(i_temp)<ttemp<=sd_temp(i_temp+1)
   DO i_temp = 1, n_temp-1
      IF ((sd_temp(i_temp)<ttemp).and.(ttemp<=sd_temp(i_temp+1))) exit
   ENDDO

   delta_t = (ttemp-sd_temp(i_temp))/(sd_temp(i_temp+1)-sd_temp(i_temp))
ENDIF

IF ((omega_min<=omega).and.(omega<=omega_max)) THEN 

   IF(.not.(cross)) THEN ! the band is on either side of the 1100 cm-1 mark

      IF(omega<omega_max) THEN
! find i_omega such that omega(i_omega)<=omega<omega(i_omega+1)
         i_omega = floor((omega-omega_min)/delta_omega)+1
         i_omega = MIN(i_omega,n_abs_data2-1)
         w1      = omega_min+delta_omega*(REAL(i_omega-1,EB))
         delta_w = (omega-w1)/delta_omega
      ELSE
! CASE i_omega = omega_max => assume THEN that i_omega = n_abs_data2 -1 
         i_omega = n_abs_data2 - 1
         w1      = omega_min+delta_omega*(REAL(i_omega-1,EB))
         delta_w = (omega-w1)/delta_omega
      ENDIF

   ELSE ! the band is on both side of 1100 cm-1
      ! test IF omega <= 1100
      IF(omega<=1100) THEN
! find i_omega such that omega(i_omega)<=omega<omega(i_omega+1)
         i_omega = floor((omega-omega_min)/delta_omega)+1
         i_omega = MIN(i_omega,n_abs_data2-1)
         w1      = omega_min+delta_omega*(REAL(i_omega-1,EB))
         delta_w = (omega-w1)/delta_omega
      ELSE
! CASE omega > 1100 cm-1
         i_omega1 = floor((1100-omega_min)/5.0_EB)+1
         i_omega2 = floor((omega-1100)/delta_omega)
         i_omega  = MIN(i_omega1+i_omega2-1,n_abs_data2-1)
         w1      = omega_min+5.0_EB*(REAL(i_omega1,EB))+delta_omega*(REAL(i_omega2-1,EB))
         delta_w = (omega-w1)/delta_omega
      ENDIF

   ENDIF

! perform bilinear interpolation to obtain values of sdweak
   interpolated_value = (1._EB-delta_w)*((1._EB-delta_t)*absorption_data(i_temp,i_omega)+ & 
                        absorption_data(i_temp+1,i_omega)*delta_t) +                      &
                        delta_w*(absorption_data(i_temp,i_omega+1)*(1._EB-delta_t)+       &
                        absorption_data(i_temp+1,i_omega+1)*delta_t)

ELSE
   interpolated_value = 0._EB
ENDIF

! ensure positivity
GET_SPECTRAL_ABSORPTION = MAX(0._EB,interpolated_value)
!------------------------------------------------------------------------------
END FUNCTION GET_SPECTRAL_ABSORPTION


!==============================================================================
ELEMENTAL REAL(EB) FUNCTION planck(temp,lambda)
!==============================================================================
! computes blackbody FUNCTION in units of w/m^2/micron/sr
! input:
! temp: temperature, unit = kelvin
! lambda: wavelength, unit = micrometers
!------------------------------------------------------------------------------
! variables passed in
REAL (EB), INTENT(IN):: temp, lambda

! local variables
REAL (EB), PARAMETER :: q1 = 1.19088e8_EB ! q1 = 2*speed_of_light^2*planck_cnst
REAL (EB), PARAMETER :: q2 = 14388._EB    ! q2 = speed_of_light*planck_cns/boltzmann
REAL (EB) :: c   ! c = lambda * temp

IF(ABS(temp)<TWO_EPSILON_EB) THEN
   planck = 0._EB
ELSE
   c = temp * lambda
   IF (q2/c > 38._EB) THEN
      planck = 0._EB
   ELSE
      planck=q1*(lambda**(-5))/(EXP(q2/c)-1._EB)
END IF
ENDIF
!------------------------------------------------------------------------------
END FUNCTION planck

!==============================================================================
ELEMENTAL REAL(EB) FUNCTION planck_wn(temp,omega)
!==============================================================================
! computes blackbody FUNCTION in units of w/m^2/cm^-1/sr
! input:
! temp : temperature, unit = kelvin
! omega: wavenumber,  unit = cm-1
! checked by vlecous1 on april 15th 2013
!------------------------------------------------------------------------------
! variables passed in
REAL (EB), INTENT(IN):: temp, omega

! local variables
REAL (EB), PARAMETER :: q1 = 1.19088e-8_EB ! q1 = 2*speed_of_light^2*planck_cnst 
REAL (EB), PARAMETER :: q2 = 1.4388_EB     
! q2 = speed_of_light*planck_cns/boltzmann
REAL (EB) :: c   ! c = temp/omega

IF(ABS(temp)<TWO_EPSILON_EB) THEN
   planck_wn = 0._EB
ELSE
   c = temp/omega
   IF (q2/c > 38._EB) THEN
      planck_wn = 0._EB
   ELSE
      planck_wn=q1*(omega**3)/(EXP(q2/c)-1._EB)
END IF
ENDIF
!------------------------------------------------------------------------------
END FUNCTION planck_wn


!==============================================================================
REAL(EB) FUNCTION integration(x,y)
!==============================================================================
! numerical integration of the vector y over the range given by vector x
! DOes not assume any regularity of x
! integration based on Simpson rule over non regular ABS(cissa
! based on lagrangian interpolation using 3 point (quadratic interpolation)
! note: io is unit of output file. needed for printing error message
!------------------------------------------------------------------------------
! variables passed in 
REAL(EB), INTENT(IN), DIMENSION(:) :: x,y

! local variable
REAL(EB) :: alpha, beta, gamma, b, a, segment_int
 
INTEGER  :: i_point, n_points, i, j ,k

! get number of elements CONTAINS in x and y
 
n_points = SIZE(y)

! test whether SIZE(x) = SIZE(y). STOP IF condition not met
IF (SIZE(x)/=n_points) THEN
   WRITE(MESSAGE, '(A)') 'Error 2: (internal) vector x and y for integration should have same SIZE(.'
   CALL SHUTDOWN(MESSAGE)
ENDIF

! test whether n_points is greater or equal to 3. IF not, print error msg and STOP radcal 
IF (n_points.lt.2) THEN
   WRITE(MESSAGE, '(A)') 'Error 3: (internal) not enough points for integration.' 
   CALL SHUTDOWN(MESSAGE)
ENDIF
! starts integration
integration = 0.0_EB

! treat first and last segment in the loop 
DO i_point = 1,n_points
 
   IF (i_point == 1) THEN
      i = i_point
      j = i_point+1
      k = i_point+2

      a = x(i)
      b = 0.5_EB*(x(i)+x(j))
   ELSEIF(i_point == n_points) THEN
      i = i_point-2
      j = i_point-1
      k = i_point

      a = 0.5_EB*(x(j)+x(k))
      b = x(k)
   ELSE
      i = i_point-1
      j = i_point
      k = i_point+1

      a = 0.5_EB*(x(i)+x(j))
      b = 0.5_EB*(x(j)+x(k))
   ENDIF

   alpha = y(i)/((x(i)-x(j))*(x(i)-x(k)))
   beta  = y(j)/((x(j)-x(i))*(x(j)-x(k)))
   gamma = y(k)/((x(k)-x(i))*(x(k)-x(j)))

   segment_int = 1.0_EB/3.0_EB*(b**3-a**3)*(alpha+beta+gamma)
   segment_int = segment_int - 0.5_EB*(b**2-a**2)*((x(j)+x(k))*alpha+beta*(x(i)+x(k))+gamma*(x(i)+x(j)))
   segment_int = segment_int + (b-a)*(alpha*x(j)*x(k)+beta*x(k)*x(i)+gamma*x(j)*x(i))

   integration = integration + segment_int

END DO 

! Finally integration = sign(b-a)*integration. This accounts for 
! descENDing order vector x.
 
integration = sign(1.0_EB,b-a)*integration

RETURN

!------------------------------------------------------------------------------
END FUNCTION integration


!==============================================================================
SUBROUTINE RCALLOC
!==============================================================================

! Old Radcal variables
ALLOCATE(gamma(4,7))
ALLOCATE(sd15(6,80))
ALLOCATE(sd(6,376))
ALLOCATE(sd7(3,16))
ALLOCATE(sd3(3,32))

!-------------------------methane data-------------------


! there are 3 bands for methane

! band #1: 1150 cm-1 - 1600 cm-1 
! band #2: 2700 cm-1 - 3250 cm-1 
! band #3: 3400 cm-1 - 5000 cm-1

ALLOCATE(sd_ch4_temp(n_temp_ch4)) 

! initialize bands wavenumber bounds for methane ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., THEN the band is calculated and not tabulated.

ALLOCATE(om_bnd_ch4(n_band_ch4,3)) 

om_bnd_ch4 = RESHAPE((/ &
   1150._EB, 2700._EB, 3400._EB,&
   1600._EB, 3250._EB, 5000._EB,&
   25._EB,   25._EB,    0._EB/),(/n_band_ch4,3/)) 

sd_ch4_temp = (/ &
   300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
   600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
   900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
   1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(sd1_ch4(n_temp_ch4,19)) 

! band #1: 1150 cm-1 - 1600 cm-1 

sd1_ch4(1:23,1:8) = RESHAPE((/ &  ! 1150-1325 cm-1
   2.674955e-03_EB, 4.122616e-03_EB, 6.116504e-03_EB, 8.394889e-03_EB, 1.066588e-02_EB, 1.270929e-02_EB, &
   1.431368e-02_EB, 1.543500e-02_EB, 1.607798e-02_EB, 1.623394e-02_EB, 1.594103e-02_EB, 1.537005e-02_EB, &
   1.456109e-02_EB, 1.356591e-02_EB, 1.253062e-02_EB, 1.145338e-02_EB, 1.037156e-02_EB, 9.319748e-03_EB, &
   8.296535e-03_EB, 7.359273e-03_EB, 6.512254e-03_EB, 5.764075e-03_EB, 5.034367e-03_EB,  &
   1.056047e-02_EB, 1.689522e-02_EB, 2.382293e-02_EB, 3.013931e-02_EB, 3.518013e-02_EB, 3.873081e-02_EB, &
   4.056505e-02_EB, 4.099451e-02_EB, 4.031689e-02_EB, 3.867720e-02_EB, 3.630281e-02_EB, 3.360423e-02_EB, &
   3.072697e-02_EB, 2.769650e-02_EB, 2.480529e-02_EB, 2.203561e-02_EB, 1.945492e-02_EB, 1.710571e-02_EB, &
   1.490918e-02_EB, 1.297063e-02_EB, 1.127359e-02_EB, 9.806206e-03_EB, 8.433200e-03_EB,  &
   6.535368e-02_EB, 8.634735e-02_EB, 1.035285e-01_EB, 1.149338e-01_EB, 1.211259e-01_EB, 1.231160e-01_EB, &
   1.212403e-01_EB, 1.167139e-01_EB, 1.104665e-01_EB, 1.027695e-01_EB, 9.409726e-02_EB, 8.537322e-02_EB, &
   7.679556e-02_EB, 6.836183e-02_EB, 6.060003e-02_EB, 5.339354e-02_EB, 4.679438e-02_EB, 4.089038e-02_EB, &
   3.545789e-02_EB, 3.072468e-02_EB, 2.660506e-02_EB, 2.306236e-02_EB, 1.979016e-02_EB,  &
   2.797310e-01_EB, 2.878774e-01_EB, 2.844637e-01_EB, 2.718418e-01_EB, 2.544158e-01_EB, 2.351542e-01_EB, &
   2.144282e-01_EB, 1.940107e-01_EB, 1.744123e-01_EB, 1.554283e-01_EB, 1.374171e-01_EB, 1.210226e-01_EB, &
   1.060935e-01_EB, 9.238361e-02_EB, 8.033186e-02_EB, 6.960028e-02_EB, 6.015950e-02_EB, 5.186871e-02_EB, &
   4.448682e-02_EB, 3.815557e-02_EB, 3.272683e-02_EB, 2.816055e-02_EB, 2.399102e-02_EB,  &
   6.055294e-01_EB, 5.085203e-01_EB, 4.314546e-01_EB, 3.661128e-01_EB, 3.120083e-01_EB, 2.675996e-01_EB, &
   2.295433e-01_EB, 1.976798e-01_EB, 1.706085e-01_EB, 1.470613e-01_EB, 1.263463e-01_EB, 1.086762e-01_EB, &
   9.339644e-02_EB, 7.991483e-02_EB, 6.850499e-02_EB, 5.859531e-02_EB, 5.007218e-02_EB, 4.275497e-02_EB, &
   3.636644e-02_EB, 3.094626e-02_EB, 2.637845e-02_EB, 2.255529e-02_EB, 1.911152e-02_EB,  &
   7.291248e-01_EB, 6.247580e-01_EB, 5.567927e-01_EB, 5.047634e-01_EB, 4.622435e-01_EB, 4.253591e-01_EB, &
   3.901550e-01_EB, 3.566718e-01_EB, 3.248608e-01_EB, 2.936112e-01_EB, 2.628808e-01_EB, 2.344160e-01_EB, &
   2.078641e-01_EB, 1.829068e-01_EB, 1.606407e-01_EB, 1.403181e-01_EB, 1.221762e-01_EB, 1.060864e-01_EB, &
   9.154534e-02_EB, 7.892764e-02_EB, 6.808643e-02_EB, 5.882322e-02_EB, 5.031768e-02_EB,  &
   1.503577e+00_EB, 1.216438e+00_EB, 1.009787e+00_EB, 8.445781e-01_EB, 7.110515e-01_EB, 6.021463e-01_EB, &
   5.099874e-01_EB, 4.327438e-01_EB, 3.680691e-01_EB, 3.125303e-01_EB, 2.644163e-01_EB, 2.242054e-01_EB, &
   1.897708e-01_EB, 1.602130e-01_EB, 1.354708e-01_EB, 1.143770e-01_EB, 9.654761e-02_EB, 8.152145e-02_EB, &
   6.854917e-02_EB, 5.774320e-02_EB, 4.875464e-02_EB, 4.129146e-02_EB, 3.467888e-02_EB,  &
   1.113787e+00_EB, 8.666106e-01_EB, 6.937600e-01_EB, 5.636603e-01_EB, 4.651973e-01_EB, 3.899245e-01_EB, &
   3.293625e-01_EB, 2.804438e-01_EB, 2.403481e-01_EB, 2.062563e-01_EB, 1.767199e-01_EB, 1.518103e-01_EB, &
   1.304026e-01_EB, 1.116263e-01_EB, 9.571162e-02_EB, 8.192030e-02_EB, 7.003835e-02_EB, 5.986022e-02_EB, &
   5.092725e-02_EB, 4.337782e-02_EB, 3.700459e-02_EB, 3.164349e-02_EB, 2.682395e-02_EB/),(/23,8/))

sd1_ch4(1:23,9:16) = RESHAPE((/ &  ! 1350-1525 cm-1
   5.859637e-01_EB, 5.741640e-01_EB, 5.446284e-01_EB, 5.016955e-01_EB, 4.533548e-01_EB, 4.052375e-01_EB, &
   3.578427e-01_EB, 3.138042e-01_EB, 2.741405e-01_EB, 2.379703e-01_EB, 2.050667e-01_EB, 1.765003e-01_EB, &
   1.515652e-01_EB, 1.294551e-01_EB, 1.106049e-01_EB, 9.430003e-02_EB, 8.028460e-02_EB, 6.831086e-02_EB, &
   5.788843e-02_EB, 4.905327e-02_EB, 4.167658e-02_EB, 3.549325e-02_EB, 2.997796e-02_EB,  &
   5.143640e-02_EB, 7.394302e-02_EB, 9.480931e-02_EB, 1.108135e-01_EB, 1.211319e-01_EB, 1.262273e-01_EB, &
   1.262909e-01_EB, 1.227527e-01_EB, 1.167145e-01_EB, 1.087709e-01_EB, 9.943737e-02_EB, 9.002538e-02_EB, &
   8.060527e-02_EB, 7.137067e-02_EB, 6.299760e-02_EB, 5.512164e-02_EB, 4.801455e-02_EB, 4.169764e-02_EB, &
   3.593587e-02_EB, 3.095498e-02_EB, 2.664278e-02_EB, 2.295572e-02_EB, 1.959682e-02_EB,  &
   7.616134e-03_EB, 9.957002e-03_EB, 1.281279e-02_EB, 1.589717e-02_EB, 1.894916e-02_EB, 2.174268e-02_EB, &
   2.395920e-02_EB, 2.555557e-02_EB, 2.648904e-02_EB, 2.671274e-02_EB, 2.629477e-02_EB, 2.540283e-02_EB, &
   2.418621e-02_EB, 2.265546e-02_EB, 2.101542e-02_EB, 1.928221e-02_EB, 1.752831e-02_EB, 1.582992e-02_EB, &
   1.416926e-02_EB, 1.261250e-02_EB, 1.121041e-02_EB, 9.951538e-03_EB, 8.728757e-03_EB,  &
   1.277402e-02_EB, 1.327177e-02_EB, 1.322108e-02_EB, 1.267471e-02_EB, 1.186190e-02_EB, 1.093562e-02_EB, &
   9.920972e-03_EB, 8.920419e-03_EB, 7.971129e-03_EB, 7.057028e-03_EB, 6.197372e-03_EB, 5.429023e-03_EB, &
   4.732635e-03_EB, 4.092017e-03_EB, 3.548360e-03_EB, 3.062048e-03_EB, 2.633809e-03_EB, 2.263255e-03_EB, &
   1.935254e-03_EB, 1.656795e-03_EB, 1.416023e-03_EB, 1.215289e-03_EB, 1.031774e-03_EB,  &
   8.014464e-03_EB, 7.194752e-03_EB, 6.553713e-03_EB, 6.000398e-03_EB, 5.531684e-03_EB, 5.121393e-03_EB, &
   4.722456e-03_EB, 4.343916e-03_EB, 3.983544e-03_EB, 3.621247e-03_EB, 3.268422e-03_EB, 2.936068e-03_EB, &
   2.617151e-03_EB, 2.322363e-03_EB, 2.051315e-03_EB, 1.804840e-03_EB, 1.578471e-03_EB, 1.378963e-03_EB, &
   1.197596e-03_EB, 1.036958e-03_EB, 9.007021e-04_EB, 7.820103e-04_EB, 6.693842e-04_EB,  &
   2.792665e-03_EB, 2.575454e-03_EB, 2.537125e-03_EB, 2.578233e-03_EB, 2.648059e-03_EB, 2.713453e-03_EB, &
   2.741224e-03_EB, 2.722787e-03_EB, 2.662395e-03_EB, 2.559402e-03_EB, 2.415338e-03_EB, 2.254521e-03_EB, &
   2.077537e-03_EB, 1.891657e-03_EB, 1.708933e-03_EB, 1.533948e-03_EB, 1.367874e-03_EB, 1.212162e-03_EB, &
   1.064882e-03_EB, 9.339250e-04_EB, 8.169687e-04_EB, 7.180516e-04_EB, 6.220503e-04_EB,  &
   8.239408e-04_EB, 8.207289e-04_EB, 8.769232e-04_EB, 9.471812e-04_EB, 1.013068e-03_EB, 1.072041e-03_EB, &
   1.109555e-03_EB, 1.125480e-03_EB, 1.120373e-03_EB, 1.094883e-03_EB, 1.043855e-03_EB, 9.887558e-04_EB, &
   9.221581e-04_EB, 8.517148e-04_EB, 7.786200e-04_EB, 7.086512e-04_EB, 6.363394e-04_EB, 5.709486e-04_EB, &
   5.037963e-04_EB, 4.478632e-04_EB, 3.933607e-04_EB, 3.461225e-04_EB, 3.036065e-04_EB,  &
   1.787674e-02_EB, 1.651301e-02_EB, 1.538618e-02_EB, 1.427021e-02_EB, 1.321135e-02_EB, 1.220556e-02_EB, &
   1.120017e-02_EB, 1.022271e-02_EB, 9.293408e-03_EB, 8.381487e-03_EB, 7.491250e-03_EB, 6.673691e-03_EB, &
   5.914957e-03_EB, 5.201979e-03_EB, 4.566356e-03_EB, 3.988961e-03_EB, 3.472458e-03_EB, 3.017625e-03_EB, &
   2.605493e-03_EB, 2.248445e-03_EB, 1.941774e-03_EB, 1.676885e-03_EB, 1.438512e-03_EB/),(/23,8/))

sd1_ch4(1:23,17:19) = RESHAPE((/ &  ! 1550-1600 cm-1
   3.842192e-03_EB, 4.909077e-03_EB, 5.876585e-03_EB, 6.637230e-03_EB, 7.159498e-03_EB, 7.449619e-03_EB, &
   7.502474e-03_EB, 7.374704e-03_EB, 7.105819e-03_EB, 6.720456e-03_EB, 6.241450e-03_EB, 5.729337e-03_EB, &
   5.213056e-03_EB, 4.685851e-03_EB, 4.183214e-03_EB, 3.714659e-03_EB, 3.276286e-03_EB, 2.879051e-03_EB, &
   2.508044e-03_EB, 2.182116e-03_EB, 1.899894e-03_EB, 1.653926e-03_EB, 1.423668e-03_EB,  &
   1.587677e-03_EB, 1.446538e-03_EB, 1.552901e-03_EB, 1.804439e-03_EB, 2.128742e-03_EB, 2.458565e-03_EB, &
   2.743898e-03_EB, 2.955974e-03_EB, 3.088835e-03_EB, 3.134690e-03_EB, 3.095698e-03_EB, 3.002805e-03_EB, &
   2.862955e-03_EB, 2.683299e-03_EB, 2.486550e-03_EB, 2.282038e-03_EB, 2.073847e-03_EB, 1.873143e-03_EB, &
   1.674361e-03_EB, 1.489429e-03_EB, 1.321804e-03_EB, 1.173554e-03_EB, 1.026702e-03_EB,  &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB/),(/23,3/))

ALLOCATE(sd2_ch4(n_temp_ch4,23)) 

! band #2: 2700 cm-1 - 3250 cm-1 

sd2_ch4(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
   1.217530e-02_EB, 1.117112e-02_EB, 1.018335e-02_EB, 9.113159e-03_EB, 8.069237e-03_EB, 7.089143e-03_EB, &
   6.164984e-03_EB, 5.339743e-03_EB, 4.605161e-03_EB, 3.952813e-03_EB, 3.369126e-03_EB, 2.873535e-03_EB, &
   2.442004e-03_EB, 2.069306e-03_EB, 1.752468e-03_EB, 1.483929e-03_EB, 1.255931e-03_EB, 1.056565e-03_EB, &
   8.920522e-04_EB, 7.508758e-04_EB, 6.317838e-04_EB, 5.352971e-04_EB, 4.530464e-04_EB,  &
   2.513916e-02_EB, 2.114376e-02_EB, 1.780513e-02_EB, 1.491293e-02_EB, 1.247243e-02_EB, 1.043356e-02_EB, &
   8.698933e-03_EB, 7.256269e-03_EB, 6.066735e-03_EB, 5.057876e-03_EB, 4.209670e-03_EB, 3.508303e-03_EB, &
   2.922322e-03_EB, 2.428483e-03_EB, 2.028794e-03_EB, 1.689558e-03_EB, 1.409261e-03_EB, 1.176715e-03_EB, &
   9.788265e-04_EB, 8.167162e-04_EB, 6.823326e-04_EB, 5.730157e-04_EB, 4.786727e-04_EB,  &
   2.557782e-02_EB, 2.133345e-02_EB, 1.832454e-02_EB, 1.586992e-02_EB, 1.378920e-02_EB, 1.202910e-02_EB, &
   1.045539e-02_EB, 9.073223e-03_EB, 7.878108e-03_EB, 6.813401e-03_EB, 5.858711e-03_EB, 5.036798e-03_EB, &
   4.322566e-03_EB, 3.689893e-03_EB, 3.157259e-03_EB, 2.690401e-03_EB, 2.292743e-03_EB, 1.951361e-03_EB, &
   1.658405e-03_EB, 1.405787e-03_EB, 1.195120e-03_EB, 1.018165e-03_EB, 8.624014e-04_EB,  &
   3.267234e-02_EB, 2.794974e-02_EB, 2.414894e-02_EB, 2.075791e-02_EB, 1.780983e-02_EB, 1.531517e-02_EB, &
   1.312860e-02_EB, 1.127375e-02_EB, 9.704346e-03_EB, 8.352011e-03_EB, 7.159191e-03_EB, 6.165400e-03_EB, &
   5.304681e-03_EB, 4.552437e-03_EB, 3.918329e-03_EB, 3.371374e-03_EB, 2.897155e-03_EB, 2.490584e-03_EB, &
   2.139129e-03_EB, 1.834560e-03_EB, 1.579448e-03_EB, 1.362473e-03_EB, 1.164728e-03_EB,  &
   4.729257e-02_EB, 4.020146e-02_EB, 3.655117e-02_EB, 3.486375e-02_EB, 3.451183e-02_EB, 3.486976e-02_EB, &
   3.526911e-02_EB, 3.549943e-02_EB, 3.539419e-02_EB, 3.477334e-02_EB, 3.359574e-02_EB, 3.209631e-02_EB, &
   3.033091e-02_EB, 2.826035e-02_EB, 2.617425e-02_EB, 2.398793e-02_EB, 2.182750e-02_EB, 1.975692e-02_EB, &
   1.770980e-02_EB, 1.581216e-02_EB, 1.409330e-02_EB, 1.256203e-02_EB, 1.106706e-02_EB,  &
   3.967221e-02_EB, 3.764359e-02_EB, 3.949747e-02_EB, 4.277841e-02_EB, 4.623172e-02_EB, 4.903094e-02_EB, &
   5.058120e-02_EB, 5.097742e-02_EB, 5.029388e-02_EB, 4.863988e-02_EB, 4.605751e-02_EB, 4.308689e-02_EB, &
   3.983596e-02_EB, 3.638363e-02_EB, 3.298206e-02_EB, 2.967115e-02_EB, 2.650179e-02_EB, 2.356866e-02_EB, &
   2.077236e-02_EB, 1.827032e-02_EB, 1.605444e-02_EB, 1.411180e-02_EB, 1.226114e-02_EB,  &
   7.223868e-02_EB, 8.544183e-02_EB, 9.962661e-02_EB, 1.104863e-01_EB, 1.171420e-01_EB, 1.196787e-01_EB, &
   1.182118e-01_EB, 1.139658e-01_EB, 1.078438e-01_EB, 1.002471e-01_EB, 9.161471e-02_EB, 8.298437e-02_EB, &
   7.452297e-02_EB, 6.622584e-02_EB, 5.857678e-02_EB, 5.147962e-02_EB, 4.505849e-02_EB, 3.932179e-02_EB, &
   3.408714e-02_EB, 2.948172e-02_EB, 2.551286e-02_EB, 2.211846e-02_EB, 1.897647e-02_EB,  &
   1.759866e-01_EB, 2.015363e-01_EB, 2.180490e-01_EB, 2.232613e-01_EB, 2.203020e-01_EB, 2.118597e-01_EB, &
   1.990075e-01_EB, 1.841205e-01_EB, 1.683887e-01_EB, 1.521963e-01_EB, 1.360646e-01_EB, 1.209307e-01_EB, &
   1.069128e-01_EB, 9.382435e-02_EB, 8.206146e-02_EB, 7.155641e-02_EB, 6.219310e-02_EB, 5.392479e-02_EB, &
   4.648837e-02_EB, 4.008493e-02_EB, 3.458414e-02_EB, 2.988156e-02_EB, 2.558176e-02_EB/),(/23,8/))

sd2_ch4(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
   3.311709e-01_EB, 3.344762e-01_EB, 3.267663e-01_EB, 3.099531e-01_EB, 2.886161e-01_EB, 2.659881e-01_EB, &
   2.423456e-01_EB, 2.191913e-01_EB, 1.973904e-01_EB, 1.764196e-01_EB, 1.564157e-01_EB, 1.382986e-01_EB, &
   1.218590e-01_EB, 1.067421e-01_EB, 9.331013e-02_EB, 8.132822e-02_EB, 7.068471e-02_EB, 6.139407e-02_EB, &
   5.297295e-02_EB, 4.571219e-02_EB, 3.948284e-02_EB, 3.416733e-02_EB, 2.931231e-02_EB,  &
   8.235613e-01_EB, 7.238968e-01_EB, 6.324968e-01_EB, 5.466365e-01_EB, 4.712117e-01_EB, 4.058971e-01_EB, &
   3.486088e-01_EB, 2.993153e-01_EB, 2.575278e-01_EB, 2.209771e-01_EB, 1.890168e-01_EB, 1.618217e-01_EB, &
   1.384064e-01_EB, 1.180022e-01_EB, 1.008139e-01_EB, 8.601549e-02_EB, 7.332027e-02_EB, 6.253544e-02_EB, &
   5.308446e-02_EB, 4.515135e-02_EB, 3.847721e-02_EB, 3.288247e-02_EB, 2.786602e-02_EB,  &
   6.821008e-01_EB, 5.401363e-01_EB, 4.353503e-01_EB, 3.533680e-01_EB, 2.892836e-01_EB, 2.389427e-01_EB, &
   1.981961e-01_EB, 1.653239e-01_EB, 1.387030e-01_EB, 1.163771e-01_EB, 9.777404e-02_EB, 8.236369e-02_EB, &
   6.949908e-02_EB, 5.853014e-02_EB, 4.942064e-02_EB, 4.176497e-02_EB, 3.531768e-02_EB, 2.985698e-02_EB, &
   2.516886e-02_EB, 2.128060e-02_EB, 1.802439e-02_EB, 1.532867e-02_EB, 1.292454e-02_EB,  &
   5.823907e-01_EB, 4.331161e-01_EB, 3.353248e-01_EB, 2.655224e-01_EB, 2.144294e-01_EB, 1.761169e-01_EB, &
   1.472744e-01_EB, 1.241615e-01_EB, 1.057407e-01_EB, 9.045740e-02_EB, 7.745806e-02_EB, 6.667193e-02_EB, &
   5.744892e-02_EB, 4.942743e-02_EB, 4.263289e-02_EB, 3.676096e-02_EB, 3.166692e-02_EB, 2.731256e-02_EB, &
   2.343114e-02_EB, 2.014747e-02_EB, 1.734533e-02_EB, 1.497211e-02_EB, 1.281431e-02_EB,  &
   3.210821e+00_EB, 2.793045e+00_EB, 2.475891e+00_EB, 2.198508e+00_EB, 1.954455e+00_EB, 1.740419e+00_EB, &
   1.543274e+00_EB, 1.366589e+00_EB, 1.208668e+00_EB, 1.063917e+00_EB, 9.309610e-01_EB, 8.135837e-01_EB, &
   7.091172e-01_EB, 6.150707e-01_EB, 5.334009e-01_EB, 4.611724e-01_EB, 3.980142e-01_EB, 3.431784e-01_EB, &
   2.944666e-01_EB, 2.526486e-01_EB, 2.171496e-01_EB, 1.870189e-01_EB, 1.596229e-01_EB,  &
   3.517151e-01_EB, 2.541996e-01_EB, 1.936159e-01_EB, 1.517283e-01_EB, 1.221750e-01_EB, 1.004308e-01_EB, &
   8.351777e-02_EB, 7.016371e-02_EB, 5.949379e-02_EB, 5.053966e-02_EB, 4.298934e-02_EB, 3.667062e-02_EB, &
   3.130112e-02_EB, 2.669085e-02_EB, 2.279490e-02_EB, 1.946235e-02_EB, 1.660600e-02_EB, 1.417075e-02_EB, &
   1.204697e-02_EB, 1.025668e-02_EB, 8.756988e-03_EB, 7.490243e-03_EB, 6.351877e-03_EB,  &
   7.921358e-01_EB, 6.055366e-01_EB, 4.798618e-01_EB, 3.867414e-01_EB, 3.174199e-01_EB, 2.642038e-01_EB, &
   2.216893e-01_EB, 1.874171e-01_EB, 1.594862e-01_EB, 1.359448e-01_EB, 1.159382e-01_EB, 9.908354e-02_EB, &
   8.473452e-02_EB, 7.233655e-02_EB, 6.187191e-02_EB, 5.286391e-02_EB, 4.518233e-02_EB, 3.863453e-02_EB, &
   3.288079e-02_EB, 2.803528e-02_EB, 2.394608e-02_EB, 2.051761e-02_EB, 1.743587e-02_EB,  &
   1.334409e+00_EB, 1.095178e+00_EB, 9.120448e-01_EB, 7.631562e-01_EB, 6.425576e-01_EB, 5.444728e-01_EB, &
   4.621443e-01_EB, 3.938260e-01_EB, 3.368579e-01_EB, 2.880892e-01_EB, 2.459378e-01_EB, 2.104005e-01_EB, &
   1.799327e-01_EB, 1.535802e-01_EB, 1.314063e-01_EB, 1.121907e-01_EB, 9.581696e-02_EB, 8.185813e-02_EB, &
   6.964125e-02_EB, 5.934825e-02_EB, 5.067442e-02_EB, 4.339843e-02_EB, 3.686411e-02_EB/),(/23,8/))

sd2_ch4(1:23,17:23) = RESHAPE((/ &  ! 3100-3250 cm-1
   9.674073e-01_EB, 9.013091e-01_EB, 8.256467e-01_EB, 7.423397e-01_EB, 6.619262e-01_EB, 5.879801e-01_EB, &
   5.191373e-01_EB, 4.573467e-01_EB, 4.022801e-01_EB, 3.523929e-01_EB, 3.072101e-01_EB, 2.674592e-01_EB, &
   2.326179e-01_EB, 2.012004e-01_EB, 1.742107e-01_EB, 1.505433e-01_EB, 1.298551e-01_EB, 1.119451e-01_EB, &
   9.604141e-02_EB, 8.242372e-02_EB, 7.089244e-02_EB, 6.110245e-02_EB, 5.221546e-02_EB,  &
   3.658653e-01_EB, 4.184364e-01_EB, 4.441884e-01_EB, 4.457809e-01_EB, 4.319013e-01_EB, 4.088715e-01_EB, &
   3.793819e-01_EB, 3.474065e-01_EB, 3.151904e-01_EB, 2.830904e-01_EB, 2.516607e-01_EB, 2.228145e-01_EB, &
   1.964195e-01_EB, 1.718178e-01_EB, 1.502324e-01_EB, 1.307846e-01_EB, 1.135071e-01_EB, 9.842138e-02_EB, &
   8.481592e-02_EB, 7.309981e-02_EB, 6.307130e-02_EB, 5.451710e-02_EB, 4.669207e-02_EB,  &
   8.373191e-02_EB, 1.229238e-01_EB, 1.573970e-01_EB, 1.825905e-01_EB, 1.975753e-01_EB, 2.040052e-01_EB, &
   2.024706e-01_EB, 1.957505e-01_EB, 1.855755e-01_EB, 1.726368e-01_EB, 1.580202e-01_EB, 1.433211e-01_EB, &
   1.287520e-01_EB, 1.144879e-01_EB, 1.013973e-01_EB, 8.925648e-02_EB, 7.818868e-02_EB, 6.826146e-02_EB, &
   5.920229e-02_EB, 5.129713e-02_EB, 4.441206e-02_EB, 3.851923e-02_EB, 3.306979e-02_EB,  &
   1.457007e-02_EB, 2.736439e-02_EB, 4.340085e-02_EB, 6.002981e-02_EB, 7.506875e-02_EB, 8.727601e-02_EB, &
   9.581284e-02_EB, 1.007875e-01_EB, 1.027285e-01_EB, 1.017913e-01_EB, 9.839048e-02_EB, 9.365514e-02_EB, &
   8.789001e-02_EB, 8.123320e-02_EB, 7.449732e-02_EB, 6.766345e-02_EB, 6.100922e-02_EB, 5.470686e-02_EB, &
   4.865358e-02_EB, 4.309088e-02_EB, 3.810993e-02_EB, 3.369467e-02_EB, 2.945916e-02_EB,  &
   1.191095e-03_EB, 1.552045e-03_EB, 2.059789e-03_EB, 2.649547e-03_EB, 3.236595e-03_EB, 3.762561e-03_EB, &
   4.165025e-03_EB, 4.448429e-03_EB, 4.600555e-03_EB, 4.628224e-03_EB, 4.534948e-03_EB, 4.374710e-03_EB, &
   4.151710e-03_EB, 3.880527e-03_EB, 3.599880e-03_EB, 3.300161e-03_EB, 2.998249e-03_EB, 2.715486e-03_EB, &
   2.425276e-03_EB, 2.164309e-03_EB, 1.929373e-03_EB, 1.713087e-03_EB, 1.506010e-03_EB,  &
   4.538291e-04_EB, 4.599393e-04_EB, 4.587246e-04_EB, 4.526653e-04_EB, 4.375098e-04_EB, 4.167315e-04_EB, &
   3.918634e-04_EB, 3.644121e-04_EB, 3.335314e-04_EB, 3.072189e-04_EB, 2.753771e-04_EB, 2.453914e-04_EB, &
   2.183408e-04_EB, 1.968316e-04_EB, 1.704817e-04_EB, 1.507858e-04_EB, 1.319531e-04_EB, 1.134243e-04_EB, &
   9.801919e-05_EB, 8.900765e-05_EB, 7.481087e-05_EB, 6.511635e-05_EB, 5.390619e-05_EB,  &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, &
   0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB, 0.000000e+00_EB/),(/23,7/))

!-------------------------ethane data-------------------


! there are 3 bands for ethane

! band #1: 730 cm-1 - 1095 cm-1 
! band #2: 1250 cm-1 - 1700 cm-1 
! band #3: 2550 cm-1 - 3375 cm-1 

ALLOCATE(sd_c2h6_temp(n_temp_c2h6)) 

! initialize bands wavenumber bounds for ethane ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_c2h6(n_band_c2h6,3)) 
ALLOCATE(be_c2h6(n_band_c2h6)) 

om_bnd_c2h6 = RESHAPE((/ &
   730._EB, 1250._EB, 2550._EB, &
   1095._EB, 1700._EB, 3375._EB, &
   5._EB, 25._EB, 25._EB/),(/n_band_c2h6,3/)) 

sd_c2h6_temp = (/ &
   296._EB, 400._EB, 450._EB, 500._EB, 600._EB, 800._EB,&
   1000._EB/)

be_c2h6 = (/ &
   1.000_EB, 1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_c2h6(n_temp_c2h6,74)) 

! band #1: 730 cm-1 - 1095 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.69096 % 

sd1_c2h6(1:7,1:8) = RESHAPE((/ &  ! 730-765 cm-1
   2.794554e-04_EB, 3.019036e-04_EB, 3.080672e-04_EB, 2.762484e-04_EB, 2.579677e-04_EB, 3.235741e-04_EB, &
   3.596358e-04_EB, &
   3.869441e-04_EB, 3.054332e-04_EB, 4.065149e-04_EB, 3.649137e-04_EB, 3.706177e-04_EB, 3.794478e-04_EB, &
   3.418980e-04_EB, &
   3.961398e-04_EB, 3.859808e-04_EB, 4.060006e-04_EB, 3.536294e-03_EB, 3.591171e-04_EB, 3.754312e-04_EB, &
   3.669844e-04_EB, &
   7.822607e-04_EB, 3.453771e-04_EB, 2.262803e-02_EB, 1.021724e-02_EB, 1.653952e-02_EB, 6.505243e-03_EB, &
   3.164663e-04_EB, &
   1.410615e-02_EB, 3.210967e-04_EB, 6.360977e-02_EB, 1.166115e-02_EB, 1.932412e+00_EB, 1.575805e-02_EB, &
   3.405423e-04_EB, &
   1.497481e-02_EB, 3.592094e-04_EB, 1.309596e-02_EB, 7.475906e-03_EB, 1.892167e+00_EB, 1.483042e-02_EB, &
   2.983060e-04_EB, &
   4.020853e-02_EB, 4.386772e-04_EB, 7.714318e-01_EB, 2.671185e-02_EB, 1.498279e+00_EB, 5.703491e-03_EB, &
   4.203063e-04_EB, &
   5.124663e-02_EB, 1.551881e-01_EB, 3.195236e-01_EB, 2.525727e-02_EB, 1.027149e+00_EB, 1.160139e-02_EB, &
   3.688548e-04_EB/),(/7,8/))

sd1_c2h6(1:7,9:16) = RESHAPE((/ &  ! 770-805 cm-1
   6.559937e-02_EB, 3.631953e-01_EB, 1.309011e-01_EB, 5.307552e-02_EB, 1.606841e-01_EB, 1.891214e-02_EB, &
   3.202870e-04_EB, &
   1.060830e-01_EB, 2.002917e-01_EB, 3.031810e-01_EB, 6.540746e-02_EB, 2.598325e-01_EB, 3.968316e-01_EB, &
   3.244087e-04_EB, &
   1.068733e-01_EB, 2.931690e-01_EB, 1.405611e-01_EB, 5.788228e-02_EB, 8.277802e-01_EB, 9.351959e-01_EB, &
   3.515174e-04_EB, &
   1.515345e-01_EB, 2.276208e-01_EB, 1.256607e-01_EB, 4.591616e-01_EB, 2.628483e+00_EB, 8.399321e-01_EB, &
   3.078236e-04_EB, &
   1.710208e-01_EB, 7.250193e-01_EB, 9.620347e-02_EB, 1.004229e-01_EB, 1.608952e+00_EB, 9.263801e-01_EB, &
   4.007032e-04_EB, &
   2.129896e-01_EB, 3.313998e-01_EB, 1.381292e-01_EB, 2.290180e-01_EB, 1.075200e+00_EB, 7.619251e-01_EB, &
   4.125753e-04_EB, &
   2.536966e-01_EB, 2.974149e-01_EB, 1.376301e-01_EB, 2.140445e-01_EB, 1.112248e+00_EB, 1.347530e+00_EB, &
   3.626265e-04_EB, &
   2.932385e-01_EB, 2.623213e-01_EB, 1.391961e-01_EB, 1.394845e-01_EB, 4.998654e-01_EB, 1.363158e+00_EB, &
   3.577957e-04_EB/),(/7,8/))

sd1_c2h6(1:7,17:24) = RESHAPE((/ &  ! 810-845 cm-1
   3.183925e-01_EB, 2.940497e-01_EB, 1.942798e-01_EB, 1.955972e-01_EB, 2.139000e-01_EB, 1.687942e+00_EB, &
   4.048250e-04_EB, &
   3.562908e-01_EB, 2.514123e-01_EB, 1.948319e-01_EB, 2.290250e-01_EB, 6.104920e-01_EB, 1.456760e+00_EB, &
   3.298797e-04_EB, &
   3.457713e-01_EB, 2.350119e-01_EB, 2.001751e-01_EB, 3.161435e-01_EB, 5.531987e-01_EB, 1.962497e-01_EB, &
   3.488040e-04_EB, &
   3.720387e-01_EB, 2.431526e-01_EB, 2.192578e-01_EB, 1.640348e-01_EB, 9.663389e-01_EB, 6.727370e-02_EB, &
   3.578073e-04_EB, &
   3.527769e-01_EB, 2.743749e-01_EB, 2.170241e-01_EB, 1.992689e-01_EB, 6.608995e-01_EB, 6.062421e-02_EB, &
   3.421209e-04_EB, &
   3.345994e-01_EB, 2.376781e-01_EB, 2.379620e-01_EB, 2.238005e-01_EB, 1.419365e-01_EB, 6.753641e-02_EB, &
   2.941644e-04_EB, &
   3.110441e-01_EB, 2.342208e-01_EB, 1.884871e-01_EB, 2.033274e-01_EB, 3.626503e-01_EB, 7.846996e-02_EB, &
   3.468572e-03_EB, &
   2.854259e-01_EB, 2.338848e-01_EB, 2.014852e-01_EB, 1.903667e-01_EB, 4.180199e-01_EB, 1.112186e-01_EB, &
   1.298106e-02_EB/),(/7,8/))

sd1_c2h6(1:7,25:32) = RESHAPE((/ &  ! 850-885 cm-1
   2.751807e-01_EB, 2.333134e-01_EB, 1.748430e-01_EB, 1.437778e-01_EB, 2.133839e-01_EB, 1.309889e-01_EB, &
   2.178703e-02_EB, &
   2.593350e-01_EB, 2.501009e-01_EB, 1.770492e-01_EB, 1.695213e-01_EB, 1.244732e-01_EB, 1.113414e-01_EB, &
   4.564412e-02_EB, &
   2.216701e-01_EB, 1.761779e-01_EB, 1.650476e-01_EB, 1.312825e-01_EB, 1.553354e-01_EB, 1.429574e-01_EB, &
   3.560235e-02_EB, &
   2.022044e-01_EB, 1.510443e-01_EB, 2.176619e-01_EB, 1.572246e-01_EB, 1.135620e-01_EB, 7.703252e-02_EB, &
   2.142633e-02_EB, &
   1.812600e-01_EB, 1.353210e-01_EB, 1.613117e-01_EB, 2.734768e-01_EB, 3.022567e-01_EB, 6.876571e-02_EB, &
   2.406288e-02_EB, &
   1.343889e-01_EB, 1.259396e-01_EB, 1.578672e-01_EB, 1.373682e-01_EB, 2.232547e-01_EB, 6.066181e-02_EB, &
   5.509858e-02_EB, &
   1.095603e-01_EB, 8.840111e-02_EB, 1.614927e-01_EB, 1.024326e-01_EB, 4.054561e-01_EB, 5.622874e-02_EB, &
   8.081915e-02_EB, &
   8.280573e-02_EB, 1.072322e-01_EB, 1.083568e-01_EB, 8.976527e-02_EB, 7.891590e-02_EB, 8.667161e-02_EB, &
   1.032621e-01_EB/),(/7,8/))

sd1_c2h6(1:7,33:40) = RESHAPE((/ &  ! 890-925 cm-1
   6.315700e-02_EB, 1.839109e-01_EB, 9.563206e-02_EB, 1.558355e-01_EB, 1.009801e-01_EB, 4.350811e-01_EB, &
   2.056632e-01_EB, &
   4.243703e-02_EB, 1.021804e-01_EB, 5.212402e-02_EB, 1.014279e-01_EB, 1.625632e-01_EB, 6.953883e-01_EB, &
   7.419429e-02_EB, &
   3.637108e-02_EB, 3.530674e-01_EB, 3.834881e-02_EB, 4.227429e-01_EB, 1.077391e-01_EB, 6.239971e-01_EB, &
   2.133274e-01_EB, &
   2.862217e-02_EB, 6.251917e-02_EB, 2.879530e-02_EB, 3.363107e-02_EB, 3.988006e-02_EB, 8.289923e-01_EB, &
   2.326386e-01_EB, &
   5.105750e-02_EB, 2.945599e-03_EB, 2.656328e-02_EB, 2.709003e-01_EB, 7.963168e-01_EB, 9.731125e-01_EB, &
   2.237348e-01_EB, &
   2.094406e-02_EB, 2.532325e-04_EB, 8.603313e-02_EB, 1.258408e-01_EB, 6.987088e-02_EB, 9.704437e-01_EB, &
   7.049559e-01_EB, &
   1.382001e-02_EB, 2.164374e-04_EB, 5.300388e-03_EB, 4.040185e-01_EB, 5.065259e-01_EB, 1.779278e+00_EB, &
   9.229669e-02_EB, &
   1.013101e-02_EB, 3.460394e-04_EB, 2.566832e-03_EB, 6.388567e-02_EB, 2.147372e-01_EB, 1.237580e+00_EB, &
   1.062539e-01_EB/),(/7,8/))

sd1_c2h6(1:7,41:48) = RESHAPE((/ &  ! 930-965 cm-1
   2.805218e-03_EB, 3.631322e-04_EB, 2.361065e-03_EB, 2.853649e-03_EB, 1.391981e-02_EB, 1.135302e+00_EB, &
   1.345327e-01_EB, &
   2.616753e-03_EB, 3.519962e-04_EB, 2.459903e-04_EB, 2.436123e-01_EB, 3.966290e-03_EB, 1.048639e+00_EB, &
   5.178554e-02_EB, &
   5.158410e-04_EB, 3.832391e-04_EB, 2.145139e-04_EB, 4.559609e-02_EB, 7.809468e-04_EB, 1.207831e+00_EB, &
   8.442327e-02_EB, &
   2.986755e-04_EB, 3.753807e-04_EB, 2.668659e-04_EB, 2.070056e-02_EB, 2.165559e-04_EB, 1.341963e+00_EB, &
   1.118404e-01_EB, &
   4.009616e-04_EB, 3.667449e-04_EB, 3.857376e-04_EB, 2.884928e-04_EB, 2.499297e-04_EB, 1.164307e+00_EB, &
   1.179168e-01_EB, &
   3.855251e-04_EB, 4.011809e-04_EB, 4.248678e-04_EB, 4.035946e-04_EB, 9.176744e-04_EB, 1.155933e+00_EB, &
   1.250155e-01_EB, &
   2.767518e-04_EB, 3.807194e-04_EB, 3.312430e-04_EB, 3.468054e-04_EB, 3.815901e-04_EB, 5.697397e-01_EB, &
   1.342966e-01_EB, &
   3.522300e-04_EB, 3.241922e-04_EB, 2.938979e-04_EB, 3.128230e-04_EB, 3.702311e-04_EB, 1.288695e-01_EB, &
   1.274616e-01_EB/),(/7,8/))

sd1_c2h6(1:7,49:56) = RESHAPE((/ &  ! 970-1005 cm-1
   3.734602e-04_EB, 3.927916e-04_EB, 3.867028e-04_EB, 3.488727e-04_EB, 3.129810e-04_EB, 5.653561e-03_EB, &
   1.114115e-01_EB, &
   3.090640e-04_EB, 3.182551e-04_EB, 3.234946e-04_EB, 3.977928e-04_EB, 3.092056e-04_EB, 3.223704e-04_EB, &
   2.110842e-01_EB, &
   3.694516e-04_EB, 3.386979e-04_EB, 3.942745e-04_EB, 4.268219e-04_EB, 4.823033e-04_EB, 3.062293e-04_EB, &
   1.371742e-01_EB, &
   3.449803e-04_EB, 3.878879e-04_EB, 4.122692e-04_EB, 3.550996e-04_EB, 3.586092e-04_EB, 3.305017e-04_EB, &
   4.857086e-01_EB, &
   3.341938e-04_EB, 4.611333e-04_EB, 2.824504e-04_EB, 3.440987e-04_EB, 3.632731e-04_EB, 3.190974e-04_EB, &
   1.036874e-01_EB, &
   3.732784e-04_EB, 3.561815e-04_EB, 3.078580e-04_EB, 4.006334e-04_EB, 3.049231e-04_EB, 3.978448e-04_EB, &
   4.974197e-02_EB, &
   2.979406e-04_EB, 3.116971e-04_EB, 3.889874e-04_EB, 3.786843e-04_EB, 3.072701e-04_EB, 2.867540e-04_EB, &
   4.862621e-02_EB, &
   3.325163e-04_EB, 3.138932e-04_EB, 4.250213e-04_EB, 3.119730e-04_EB, 3.732476e-04_EB, 3.717123e-04_EB, &
   4.454824e-02_EB/),(/7,8/))

sd1_c2h6(1:7,57:64) = RESHAPE((/ &  ! 1010-1045 cm-1
   3.350681e-04_EB, 3.530821e-04_EB, 4.528042e-04_EB, 3.516044e-04_EB, 3.923830e-04_EB, 3.532259e-04_EB, &
   6.542840e-02_EB, &
   3.204957e-04_EB, 3.836292e-04_EB, 3.138970e-04_EB, 3.100875e-04_EB, 3.990551e-04_EB, 3.191793e-04_EB, &
   7.406912e-02_EB, &
   3.472691e-04_EB, 3.920082e-04_EB, 3.234488e-04_EB, 4.026391e-04_EB, 3.384464e-04_EB, 2.968447e-04_EB, &
   2.322588e-01_EB, &
   3.347127e-04_EB, 3.546956e-04_EB, 2.920812e-04_EB, 4.081782e-04_EB, 3.386853e-04_EB, 3.307871e-04_EB, &
   4.463887e-01_EB, &
   3.248741e-04_EB, 4.609033e-04_EB, 3.231879e-04_EB, 3.456361e-04_EB, 3.298612e-04_EB, 3.329381e-04_EB, &
   1.035782e+00_EB, &
   4.690522e-04_EB, 3.215237e-04_EB, 4.686739e-04_EB, 3.787632e-04_EB, 3.650086e-04_EB, 3.038461e-04_EB, &
   6.376563e-01_EB, &
   3.712414e-04_EB, 3.428482e-04_EB, 3.459720e-04_EB, 3.001890e-04_EB, 3.653620e-04_EB, 3.093595e-04_EB, &
   1.783476e-01_EB, &
   3.413817e-04_EB, 3.624958e-04_EB, 3.637267e-04_EB, 3.228750e-04_EB, 2.920456e-04_EB, 3.181545e-04_EB, &
   2.620076e-01_EB/),(/7,8/))

sd1_c2h6(1:7,65:72) = RESHAPE((/ &  ! 1050-1085 cm-1
   3.325905e-04_EB, 3.330314e-04_EB, 3.617122e-04_EB, 3.505622e-04_EB, 3.869601e-04_EB, 3.838024e-04_EB, &
   5.425202e-01_EB, &
   3.314643e-04_EB, 3.409742e-04_EB, 2.922515e-04_EB, 2.908445e-04_EB, 3.423906e-04_EB, 2.970385e-04_EB, &
   6.238502e-01_EB, &
   3.492482e-04_EB, 3.606574e-04_EB, 3.274995e-04_EB, 3.213303e-04_EB, 3.179339e-04_EB, 3.480095e-04_EB, &
   4.163639e-02_EB, &
   3.199232e-04_EB, 3.955618e-04_EB, 2.798966e-04_EB, 3.155772e-04_EB, 2.916205e-04_EB, 3.756559e-04_EB, &
   1.227547e-02_EB, &
   3.705657e-04_EB, 3.404134e-04_EB, 3.590301e-04_EB, 3.440815e-04_EB, 3.347934e-04_EB, 3.248160e-04_EB, &
   1.130962e-02_EB, &
   3.884066e-04_EB, 3.017404e-04_EB, 3.954654e-04_EB, 3.017435e-04_EB, 2.609871e-04_EB, 2.759450e-04_EB, &
   2.756706e-03_EB, &
   3.676700e-04_EB, 3.019733e-04_EB, 3.533386e-04_EB, 3.287691e-04_EB, 4.033519e-04_EB, 3.370772e-04_EB, &
   3.213983e-04_EB, &
   3.262600e-04_EB, 3.615855e-04_EB, 3.142415e-04_EB, 3.309319e-04_EB, 4.228572e-04_EB, 3.478703e-04_EB, &
   3.234945e-04_EB/),(/7,8/))

sd1_c2h6(1:7,73:74) = RESHAPE((/ &  ! 1090-1095 cm-1
   2.711945e-04_EB, 2.496329e-04_EB, 3.150593e-04_EB, 2.603387e-04_EB, 2.992890e-04_EB, 3.589836e-04_EB, &
   2.656451e-04_EB, &
   2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, &
   2.145139e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_c2h6(n_temp_c2h6,74)) 

! band #1: 730 cm-1 - 1095 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.69096 % 

! print fine structure array gamma_d 

gammad1_c2h6(1:7,1:8) = RESHAPE((/ &  ! 730-765 cm-1
   2.006565e-04_EB, 2.337401e-04_EB, 2.177856e-04_EB, 2.048546e-04_EB, 1.869983e-04_EB, 2.713666e-04_EB, &
   2.905961e-04_EB, &
   2.306868e-04_EB, 2.167323e-04_EB, 3.350447e-04_EB, 2.526474e-04_EB, 2.261799e-04_EB, 2.480626e-04_EB, &
   2.521419e-04_EB, &
   2.854020e-04_EB, 3.329316e-04_EB, 1.598485e-04_EB, 3.267121e-02_EB, 2.234431e-04_EB, 2.957501e-04_EB, &
   2.173206e-04_EB, &
   2.137117e-04_EB, 2.266418e-04_EB, 4.728511e-05_EB, 2.245423e-01_EB, 3.064851e-05_EB, 8.222848e-03_EB, &
   1.859738e-04_EB, &
   3.462857e-03_EB, 2.519057e-04_EB, 3.480152e-05_EB, 3.925807e-02_EB, 1.572683e-05_EB, 3.101866e-01_EB, &
   1.995997e-04_EB, &
   9.741100e-03_EB, 2.710087e-04_EB, 2.101240e-03_EB, 7.469797e-02_EB, 2.264929e-05_EB, 1.275440e-01_EB, &
   1.730011e-04_EB, &
   7.831446e-03_EB, 2.356564e-04_EB, 2.537222e-05_EB, 5.290921e-02_EB, 3.544110e-05_EB, 1.081723e-01_EB, &
   2.635023e-04_EB, &
   1.830629e-02_EB, 5.427184e-06_EB, 1.638831e-04_EB, 1.590424e-02_EB, 4.986288e-05_EB, 2.169141e-01_EB, &
   2.505642e-04_EB/),(/7,8/))

gammad1_c2h6(1:7,9:16) = RESHAPE((/ &  ! 770-805 cm-1
   3.746173e-02_EB, 1.736644e-05_EB, 6.085991e-04_EB, 4.704301e-01_EB, 6.606423e-04_EB, 7.502293e-02_EB, &
   2.108317e-04_EB, &
   2.746203e-02_EB, 1.840834e-04_EB, 5.921780e-04_EB, 2.043773e-02_EB, 6.107989e-04_EB, 1.388557e-04_EB, &
   1.893995e-04_EB, &
   2.420680e-01_EB, 4.191667e-04_EB, 3.132099e-03_EB, 1.886145e+00_EB, 3.070676e-04_EB, 1.079068e-04_EB, &
   2.691653e-04_EB, &
   8.925252e-02_EB, 1.124913e-03_EB, 8.738915e-03_EB, 1.559750e-03_EB, 2.395539e-04_EB, 1.072791e-04_EB, &
   1.752395e-04_EB, &
   1.067081e-01_EB, 6.895830e-04_EB, 4.690182e-02_EB, 4.174376e-02_EB, 4.527850e-04_EB, 2.137964e-04_EB, &
   2.588101e-04_EB, &
   8.746003e-02_EB, 2.440063e-03_EB, 2.130236e-02_EB, 6.354666e-03_EB, 6.671299e-04_EB, 3.378605e-04_EB, &
   2.343097e-04_EB, &
   9.852102e-02_EB, 4.561646e-03_EB, 3.949280e-02_EB, 1.123806e-02_EB, 8.233218e-04_EB, 2.522466e-04_EB, &
   2.226954e-04_EB, &
   1.207618e-01_EB, 8.941124e-03_EB, 1.085484e-01_EB, 3.349804e-02_EB, 1.725479e-03_EB, 3.076476e-04_EB, &
   1.938371e-04_EB/),(/7,8/))

gammad1_c2h6(1:7,17:24) = RESHAPE((/ &  ! 810-845 cm-1
   1.241365e-01_EB, 9.920833e-03_EB, 2.540899e-02_EB, 1.887039e-02_EB, 6.826400e-03_EB, 2.126909e-04_EB, &
   2.773335e-04_EB, &
   1.417207e-01_EB, 2.020016e-02_EB, 3.589721e-02_EB, 1.758195e-02_EB, 1.968415e-03_EB, 2.289210e-04_EB, &
   2.228792e-04_EB, &
   1.684374e-01_EB, 2.315876e-02_EB, 3.432510e-02_EB, 1.023915e-02_EB, 2.040480e-03_EB, 2.528208e-03_EB, &
   2.461503e-04_EB, &
   2.993636e-01_EB, 4.679863e-02_EB, 4.993611e-02_EB, 5.799012e-01_EB, 1.375382e-03_EB, 5.540907e+00_EB, &
   2.729510e-04_EB, &
   3.763363e-01_EB, 3.147273e-02_EB, 5.066712e-02_EB, 5.649368e-02_EB, 2.734156e-03_EB, 8.088974e+00_EB, &
   2.352079e-04_EB, &
   3.923696e-01_EB, 4.535395e-02_EB, 3.848111e-02_EB, 2.837120e-02_EB, 3.903926e-02_EB, 4.193583e+00_EB, &
   1.708561e-04_EB, &
   6.116601e-01_EB, 3.911633e-02_EB, 9.404679e-02_EB, 3.598464e-02_EB, 4.985811e-03_EB, 2.536374e+00_EB, &
   1.671554e-02_EB, &
   6.488105e-01_EB, 3.434365e-02_EB, 5.620492e-02_EB, 4.579674e-02_EB, 3.831981e-03_EB, 2.040147e-02_EB, &
   1.535947e-01_EB/),(/7,8/))

gammad1_c2h6(1:7,25:32) = RESHAPE((/ &  ! 850-885 cm-1
   2.185282e-01_EB, 2.624089e-02_EB, 8.281180e-02_EB, 1.337021e-01_EB, 1.164389e-02_EB, 1.441667e-02_EB, &
   1.326192e+00_EB, &
   1.054436e-01_EB, 1.564415e-02_EB, 4.232910e-02_EB, 2.542335e-02_EB, 9.182689e-02_EB, 2.979608e-02_EB, &
   3.412059e+00_EB, &
   1.037085e-01_EB, 2.307616e-02_EB, 2.791608e-02_EB, 7.390912e-02_EB, 1.792121e-02_EB, 8.709794e-03_EB, &
   9.123564e-01_EB, &
   6.557347e-02_EB, 2.400623e-02_EB, 1.184967e-02_EB, 1.759675e-02_EB, 4.025549e-02_EB, 1.650755e-01_EB, &
   4.126380e-01_EB, &
   3.854672e-02_EB, 1.545381e-02_EB, 1.176956e-02_EB, 5.353178e-03_EB, 3.790587e-03_EB, 2.857260e+00_EB, &
   4.896543e-02_EB, &
   3.823510e-02_EB, 1.163165e-02_EB, 7.527145e-03_EB, 1.044536e-02_EB, 4.277846e-03_EB, 2.782661e+00_EB, &
   9.862962e-03_EB, &
   2.698248e-02_EB, 1.682716e-02_EB, 4.532103e-03_EB, 2.286498e-02_EB, 1.772522e-03_EB, 2.075720e+00_EB, &
   2.124687e-02_EB, &
   2.066380e-02_EB, 4.723103e-03_EB, 6.343220e-03_EB, 2.397806e-02_EB, 2.607826e-02_EB, 7.702274e-03_EB, &
   1.637205e-02_EB/),(/7,8/))

gammad1_c2h6(1:7,33:40) = RESHAPE((/ &  ! 890-925 cm-1
   1.658402e-02_EB, 9.717357e-04_EB, 4.999287e-03_EB, 3.360172e-03_EB, 7.640958e-03_EB, 5.866049e-04_EB, &
   4.203039e-03_EB, &
   9.306719e-03_EB, 1.059841e-03_EB, 1.554113e-02_EB, 4.017085e-03_EB, 1.562439e-03_EB, 2.629554e-04_EB, &
   2.723941e-02_EB, &
   6.324845e-03_EB, 1.323745e-04_EB, 2.015294e-02_EB, 4.021824e-04_EB, 2.200845e-03_EB, 1.909660e-04_EB, &
   2.384092e-03_EB, &
   6.442116e-03_EB, 1.976663e-04_EB, 1.301740e-02_EB, 7.108478e-03_EB, 9.667620e-03_EB, 1.093856e-04_EB, &
   1.651473e-03_EB, &
   6.177417e-04_EB, 5.252625e-04_EB, 3.292456e-03_EB, 9.638064e-05_EB, 1.157061e-04_EB, 9.272832e-05_EB, &
   1.772000e-03_EB, &
   3.119150e-03_EB, 1.688032e-04_EB, 1.157769e-04_EB, 4.344200e-05_EB, 6.951343e-04_EB, 9.547409e-05_EB, &
   5.716592e-04_EB, &
   3.033519e-03_EB, 1.700977e-04_EB, 1.923551e-03_EB, 1.441479e-05_EB, 9.010921e-05_EB, 4.649518e-05_EB, &
   3.993180e-03_EB, &
   4.505246e-03_EB, 2.206088e-04_EB, 4.709851e-04_EB, 4.044407e-05_EB, 9.427052e-05_EB, 4.271051e-05_EB, &
   4.728030e-03_EB/),(/7,8/))

gammad1_c2h6(1:7,41:48) = RESHAPE((/ &  ! 930-965 cm-1
   8.515529e-04_EB, 2.523157e-04_EB, 3.215017e-03_EB, 8.896712e-05_EB, 2.653323e-03_EB, 2.395564e-05_EB, &
   5.128282e-03_EB, &
   1.598460e-03_EB, 2.391773e-04_EB, 2.024249e-04_EB, 4.570502e-06_EB, 1.197295e-03_EB, 2.420553e-05_EB, &
   7.176382e-01_EB, &
   5.405162e-04_EB, 2.754095e-04_EB, 1.688326e-04_EB, 1.901197e-05_EB, 1.882451e-04_EB, 2.561929e-05_EB, &
   1.225580e+00_EB, &
   2.092858e-04_EB, 2.662552e-04_EB, 2.026153e-04_EB, 1.070822e-05_EB, 1.684702e-04_EB, 2.442314e-05_EB, &
   2.444828e+01_EB, &
   2.529925e-04_EB, 2.765314e-04_EB, 2.642207e-04_EB, 2.200703e-04_EB, 1.995939e-04_EB, 2.380865e-05_EB, &
   2.685955e+01_EB, &
   2.173973e-04_EB, 2.813395e-04_EB, 2.874743e-04_EB, 3.035869e-04_EB, 6.183918e-04_EB, 2.202463e-05_EB, &
   2.519809e+00_EB, &
   1.980501e-04_EB, 2.678571e-04_EB, 2.137226e-04_EB, 2.556760e-04_EB, 2.222342e-04_EB, 1.634015e-05_EB, &
   3.727180e-02_EB, &
   2.678728e-04_EB, 2.424340e-04_EB, 1.979944e-04_EB, 2.135470e-04_EB, 2.146239e-04_EB, 1.781477e-05_EB, &
   7.049181e-03_EB/),(/7,8/))

gammad1_c2h6(1:7,49:56) = RESHAPE((/ &  ! 970-1005 cm-1
   2.364062e-04_EB, 2.712640e-04_EB, 2.842134e-04_EB, 2.567382e-04_EB, 2.256683e-04_EB, 2.498984e-05_EB, &
   4.004321e-03_EB, &
   2.157659e-04_EB, 2.064345e-04_EB, 2.288569e-04_EB, 2.738793e-04_EB, 2.050141e-04_EB, 2.197032e-04_EB, &
   1.028210e-03_EB, &
   2.197468e-04_EB, 2.679068e-04_EB, 2.304862e-04_EB, 2.905492e-04_EB, 2.725375e-04_EB, 2.159399e-04_EB, &
   1.419923e-03_EB, &
   2.335319e-04_EB, 2.681889e-04_EB, 3.239557e-04_EB, 2.569741e-04_EB, 2.238526e-04_EB, 2.007748e-04_EB, &
   2.537471e-04_EB, &
   2.263569e-04_EB, 3.026600e-04_EB, 2.089141e-04_EB, 2.243273e-04_EB, 2.228419e-04_EB, 2.279698e-04_EB, &
   7.317457e-04_EB, &
   2.753690e-04_EB, 2.292202e-04_EB, 2.318102e-04_EB, 2.597276e-04_EB, 2.082393e-04_EB, 2.415048e-04_EB, &
   2.144539e-03_EB, &
   2.130814e-04_EB, 2.160667e-04_EB, 3.005135e-04_EB, 2.457929e-04_EB, 2.262402e-04_EB, 2.271924e-04_EB, &
   1.982460e-03_EB, &
   2.260545e-04_EB, 2.078552e-04_EB, 2.252198e-04_EB, 2.292910e-04_EB, 2.617146e-04_EB, 2.719938e-04_EB, &
   9.539658e-03_EB/),(/7,8/))

gammad1_c2h6(1:7,57:64) = RESHAPE((/ &  ! 1010-1045 cm-1
   2.185131e-04_EB, 2.784174e-04_EB, 2.724240e-04_EB, 2.166217e-04_EB, 2.454826e-04_EB, 2.343101e-04_EB, &
   4.681303e-03_EB, &
   2.140851e-04_EB, 2.908903e-04_EB, 2.238573e-04_EB, 1.936346e-04_EB, 2.929472e-04_EB, 1.997998e-04_EB, &
   5.199463e-03_EB, &
   2.704701e-04_EB, 2.762782e-04_EB, 2.520038e-04_EB, 2.774566e-04_EB, 2.294524e-04_EB, 1.897369e-04_EB, &
   1.101146e-03_EB, &
   2.261254e-04_EB, 2.439262e-04_EB, 2.109405e-04_EB, 2.752945e-04_EB, 2.380712e-04_EB, 2.219084e-04_EB, &
   4.729845e-04_EB, &
   2.442952e-04_EB, 3.469938e-04_EB, 2.322461e-04_EB, 2.311834e-04_EB, 2.364839e-04_EB, 2.042171e-04_EB, &
   2.219015e-04_EB, &
   2.634592e-04_EB, 2.110525e-04_EB, 3.161689e-04_EB, 2.414846e-04_EB, 2.423626e-04_EB, 2.110280e-04_EB, &
   3.145896e-04_EB, &
   2.659955e-04_EB, 2.470635e-04_EB, 2.700914e-04_EB, 2.258166e-04_EB, 2.553061e-04_EB, 2.172449e-04_EB, &
   8.521657e-04_EB, &
   2.210763e-04_EB, 2.410264e-04_EB, 2.690630e-04_EB, 2.296574e-04_EB, 1.980775e-04_EB, 2.236759e-04_EB, &
   5.582951e-04_EB/),(/7,8/))

gammad1_c2h6(1:7,65:72) = RESHAPE((/ &  ! 1050-1085 cm-1
   2.066800e-04_EB, 2.209692e-04_EB, 2.738395e-04_EB, 2.085763e-04_EB, 2.168987e-04_EB, 2.377505e-04_EB, &
   1.629293e-04_EB, &
   2.355869e-04_EB, 2.212229e-04_EB, 2.006984e-04_EB, 2.157291e-04_EB, 2.203006e-04_EB, 2.060948e-04_EB, &
   9.639746e-05_EB, &
   2.268691e-04_EB, 2.800409e-04_EB, 1.984453e-04_EB, 2.045951e-04_EB, 1.977593e-04_EB, 2.051878e-04_EB, &
   9.621759e-04_EB, &
   2.273679e-04_EB, 2.692054e-04_EB, 2.117256e-04_EB, 2.039894e-04_EB, 1.777981e-04_EB, 2.609597e-04_EB, &
   1.607451e-03_EB, &
   2.293047e-04_EB, 2.421597e-04_EB, 2.833909e-04_EB, 2.164108e-04_EB, 2.199777e-04_EB, 2.398280e-04_EB, &
   1.282936e-03_EB, &
   2.497990e-04_EB, 2.311811e-04_EB, 2.716080e-04_EB, 2.074448e-04_EB, 1.857482e-04_EB, 1.792751e-04_EB, &
   9.092886e-04_EB, &
   2.440601e-04_EB, 2.141159e-04_EB, 2.414717e-04_EB, 2.388990e-04_EB, 2.516826e-04_EB, 2.158533e-04_EB, &
   2.049667e-04_EB, &
   2.220025e-04_EB, 2.592544e-04_EB, 1.988151e-04_EB, 1.978750e-04_EB, 2.337158e-04_EB, 2.300715e-04_EB, &
   1.991380e-04_EB/),(/7,8/))

gammad1_c2h6(1:7,73:74) = RESHAPE((/ &  ! 1090-1095 cm-1
   1.847719e-04_EB, 1.895118e-04_EB, 2.362824e-04_EB, 1.933361e-04_EB, 2.035927e-04_EB, 2.735696e-04_EB, &
   2.160458e-04_EB, &
   1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, &
   1.688326e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_c2h6(n_temp_c2h6,19)) 

! band #2: 1250 cm-1 - 1700 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 1.536 % 

sd2_c2h6(1:7,1:8) = RESHAPE((/ &  ! 1250-1425 cm-1
   2.890736e-04_EB, 2.831224e-04_EB, 2.500860e-04_EB, 2.757140e-04_EB, 3.408382e-04_EB, 8.798483e-04_EB, &
   9.996940e-04_EB, &
   3.173085e-04_EB, 3.511457e-04_EB, 3.615941e-04_EB, 3.532323e-04_EB, 3.063903e-03_EB, 5.665886e-01_EB, &
   6.009682e-01_EB, &
   7.762194e-03_EB, 3.886970e-04_EB, 3.445713e-04_EB, 3.346007e-04_EB, 9.791971e-03_EB, 3.101724e-01_EB, &
   3.854181e-01_EB, &
   5.755406e-03_EB, 2.965005e-03_EB, 3.003864e-03_EB, 2.921568e-02_EB, 2.167033e-02_EB, 2.481666e-01_EB, &
   4.975186e-01_EB, &
   2.729491e-02_EB, 4.272374e-02_EB, 4.078658e-02_EB, 5.957540e-02_EB, 5.312857e-02_EB, 7.121161e-02_EB, &
   4.478193e-02_EB, &
   1.116554e-01_EB, 1.211306e-01_EB, 1.109669e-01_EB, 1.148103e-01_EB, 9.457225e-02_EB, 9.397124e-02_EB, &
   6.487302e-02_EB, &
   2.433578e-01_EB, 2.251434e-01_EB, 1.980541e-01_EB, 1.663919e-01_EB, 1.429235e-01_EB, 1.679094e-01_EB, &
   1.043870e-01_EB, &
   3.085264e-01_EB, 2.668615e-01_EB, 2.316373e-01_EB, 1.816814e-01_EB, 1.675634e-01_EB, 1.397869e-01_EB, &
   9.888312e-02_EB/),(/7,8/))

sd2_c2h6(1:7,9:16) = RESHAPE((/ &  ! 1450-1625 cm-1
   4.739873e-01_EB, 3.630012e-01_EB, 2.939668e-01_EB, 2.366840e-01_EB, 2.000756e-01_EB, 1.538908e-01_EB, &
   1.749130e-01_EB, &
   5.226241e-01_EB, 3.793992e-01_EB, 3.109342e-01_EB, 2.476162e-01_EB, 1.999044e-01_EB, 1.438092e-01_EB, &
   3.486559e-01_EB, &
   4.230308e-01_EB, 3.178107e-01_EB, 2.777402e-01_EB, 1.894468e-01_EB, 1.805602e-01_EB, 1.554106e-01_EB, &
   1.422521e-01_EB, &
   2.761096e-01_EB, 2.597293e-01_EB, 2.257682e-01_EB, 1.514748e-01_EB, 1.401495e-01_EB, 1.253687e-01_EB, &
   7.720757e-02_EB, &
   1.111228e-01_EB, 1.672047e-01_EB, 1.760444e-01_EB, 8.528705e-02_EB, 1.082694e-01_EB, 9.903350e-02_EB, &
   6.563257e-02_EB, &
   3.538374e-02_EB, 1.038124e-01_EB, 5.979240e-02_EB, 7.230511e-02_EB, 5.579668e-02_EB, 6.549601e-02_EB, &
   1.002546e-01_EB, &
   4.798824e-03_EB, 1.584507e-02_EB, 1.391213e-02_EB, 7.436389e-01_EB, 2.920638e-02_EB, 1.145160e-01_EB, &
   1.008028e-01_EB, &
   3.514365e-04_EB, 3.078667e-04_EB, 6.289715e-04_EB, 7.024049e-03_EB, 1.079463e-02_EB, 2.595491e-01_EB, &
   1.608240e-02_EB/),(/7,8/))

sd2_c2h6(1:7,17:19) = RESHAPE((/ &  ! 1650-1700 cm-1
   3.438068e-04_EB, 3.506371e-04_EB, 3.282171e-04_EB, 3.346269e-04_EB, 1.558845e-03_EB, 6.506371e-02_EB, &
   2.338883e-01_EB, &
   2.659394e-04_EB, 3.048186e-04_EB, 2.844540e-04_EB, 2.662972e-04_EB, 3.138724e-04_EB, 3.077023e-03_EB, &
   7.859120e-03_EB, &
   2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, &
   2.145139e-04_EB/),(/7,3/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_c2h6(n_temp_c2h6,19)) 

! band #2: 1250 cm-1 - 1700 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 1.536 % 

! print fine structure array gamma_d 

gammad2_c2h6(1:7,1:8) = RESHAPE((/ &  ! 1250-1425 cm-1
   2.119975e-04_EB, 2.146327e-04_EB, 1.903252e-04_EB, 2.022541e-04_EB, 2.602782e-04_EB, 4.275791e-04_EB, &
   2.592334e-04_EB, &
   2.197390e-04_EB, 2.683644e-04_EB, 2.472885e-04_EB, 2.304333e-04_EB, 8.267459e-03_EB, 1.628893e-05_EB, &
   3.112078e-05_EB, &
   1.421119e-02_EB, 2.763312e-04_EB, 2.321580e-04_EB, 2.291520e-04_EB, 1.958198e-02_EB, 6.419678e-05_EB, &
   1.177209e-04_EB, &
   5.134518e-02_EB, 1.180015e-03_EB, 5.665049e-03_EB, 9.959465e-05_EB, 7.884425e-02_EB, 2.789112e-04_EB, &
   2.062441e-04_EB, &
   6.493055e-01_EB, 8.533910e-03_EB, 1.523050e-02_EB, 5.127453e-03_EB, 6.635840e-01_EB, 1.697909e-02_EB, &
   1.007998e+00_EB, &
   6.644469e+00_EB, 5.703315e-02_EB, 7.526340e-02_EB, 3.982739e-02_EB, 2.373472e+00_EB, 2.681325e-02_EB, &
   6.481877e-01_EB, &
   4.838580e+01_EB, 1.082873e-01_EB, 1.011916e-01_EB, 2.605553e-01_EB, 2.237840e+00_EB, 1.611963e-02_EB, &
   3.893042e-01_EB, &
   7.997133e+01_EB, 1.230787e-01_EB, 1.231629e-01_EB, 3.381738e+00_EB, 2.828329e-01_EB, 5.186004e-02_EB, &
   4.989001e-02_EB/),(/7,8/))

gammad2_c2h6(1:7,9:16) = RESHAPE((/ &  ! 1450-1625 cm-1
   7.999999e+01_EB, 1.520030e-01_EB, 1.960215e-01_EB, 3.752618e-01_EB, 2.469986e-01_EB, 4.882190e-02_EB, &
   1.110526e-02_EB, &
   7.998354e+01_EB, 1.855013e-01_EB, 1.836005e-01_EB, 2.792553e-01_EB, 2.412743e-01_EB, 4.914911e-02_EB, &
   2.585712e-03_EB, &
   8.000000e+01_EB, 1.653039e-01_EB, 1.276875e-01_EB, 4.450008e+01_EB, 1.282813e-01_EB, 1.923521e-02_EB, &
   9.527892e-03_EB, &
   8.000000e+01_EB, 7.676701e-02_EB, 7.528859e-02_EB, 3.303007e+01_EB, 3.079957e-01_EB, 2.422473e-02_EB, &
   3.925627e-02_EB, &
   7.999946e+01_EB, 2.618451e-02_EB, 2.459894e-02_EB, 1.811547e+01_EB, 7.028942e-02_EB, 2.413638e-02_EB, &
   3.482649e-02_EB, &
   3.034575e+00_EB, 3.670718e-03_EB, 3.712350e-02_EB, 1.721291e-02_EB, 2.149261e+00_EB, 1.626810e-02_EB, &
   2.799184e-03_EB, &
   2.512595e-02_EB, 8.026413e-04_EB, 5.409205e-02_EB, 1.139540e-04_EB, 7.373019e-01_EB, 1.105238e-03_EB, &
   1.244382e-03_EB, &
   2.370622e-04_EB, 2.109883e-04_EB, 4.709529e-04_EB, 1.239090e-03_EB, 2.096932e-01_EB, 1.101878e-04_EB, &
   6.601257e-03_EB/),(/7,8/))

gammad2_c2h6(1:7,17:19) = RESHAPE((/ &  ! 1650-1700 cm-1
   2.490089e-04_EB, 2.329594e-04_EB, 2.331607e-04_EB, 2.360833e-04_EB, 1.239584e-03_EB, 1.730171e-04_EB, &
   6.965657e-05_EB, &
   1.969514e-04_EB, 2.180366e-04_EB, 2.121351e-04_EB, 1.905845e-04_EB, 2.199848e-04_EB, 9.700224e-04_EB, &
   1.623935e-03_EB, &
   1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, &
   1.688326e-04_EB/),(/7,3/))

!---------------------------------------------------------------------------
ALLOCATE(sd3_c2h6(n_temp_c2h6,34)) 

! band #3: 2550 cm-1 - 3375 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.3765 % 

sd3_c2h6(1:7,1:8) = RESHAPE((/ &  ! 2550-2725 cm-1
   1.273212e-03_EB, 2.818357e-04_EB, 2.919712e-04_EB, 2.725676e-04_EB, 3.224456e-04_EB, 1.576641e-03_EB, &
   7.098264e-02_EB, &
   1.760425e-02_EB, 2.470202e-03_EB, 6.517460e-04_EB, 1.488758e-03_EB, 2.133431e-03_EB, 7.058130e-03_EB, &
   5.189907e-01_EB, &
   1.544358e-02_EB, 7.033132e-03_EB, 6.497614e-03_EB, 7.444432e-03_EB, 7.481133e-03_EB, 9.709862e-03_EB, &
   1.866477e-01_EB, &
   2.776353e-02_EB, 1.788161e-02_EB, 1.440452e-02_EB, 1.497498e-02_EB, 1.098338e-02_EB, 1.398847e-02_EB, &
   1.830505e-01_EB, &
   4.053879e-02_EB, 2.592830e-02_EB, 1.797176e-02_EB, 1.619226e-02_EB, 1.327884e-02_EB, 1.766187e-02_EB, &
   1.018847e-02_EB, &
   4.388545e-02_EB, 3.290936e-02_EB, 2.841386e-02_EB, 1.850778e-02_EB, 1.703843e-02_EB, 2.108276e-02_EB, &
   2.254688e-02_EB, &
   3.693893e-02_EB, 3.278040e-02_EB, 3.020784e-02_EB, 2.524477e-02_EB, 2.353619e-02_EB, 2.930089e-02_EB, &
   3.312656e-02_EB, &
   5.710175e-02_EB, 5.583766e-02_EB, 4.445106e-02_EB, 4.789108e-02_EB, 3.783170e-02_EB, 4.950327e-02_EB, &
   5.137114e-02_EB/),(/7,8/))

sd3_c2h6(1:7,9:16) = RESHAPE((/ &  ! 2750-2925 cm-1
   9.437750e-02_EB, 7.826799e-02_EB, 6.481721e-02_EB, 7.020402e-02_EB, 6.017769e-02_EB, 7.963757e-02_EB, &
   9.068663e-02_EB, &
   1.449953e-01_EB, 1.110552e-01_EB, 9.612305e-02_EB, 1.060636e-01_EB, 8.633211e-02_EB, 1.043293e-01_EB, &
   1.280871e-01_EB, &
   1.001382e-01_EB, 8.904423e-02_EB, 8.899341e-02_EB, 1.080298e-01_EB, 9.748937e-02_EB, 1.507712e-01_EB, &
   2.152193e-01_EB, &
   1.789013e-01_EB, 1.576201e-01_EB, 1.688876e-01_EB, 1.893222e-01_EB, 2.046234e-01_EB, 2.778477e-01_EB, &
   3.326868e-01_EB, &
   4.519875e-01_EB, 4.722980e-01_EB, 5.002552e-01_EB, 5.176398e-01_EB, 5.259894e-01_EB, 5.581873e-01_EB, &
   5.217227e-01_EB, &
   1.942344e+00_EB, 1.479985e+00_EB, 1.319915e+00_EB, 1.189402e+00_EB, 1.020515e+00_EB, 8.684235e-01_EB, &
   7.596269e-01_EB, &
   2.588837e+00_EB, 1.924667e+00_EB, 1.726463e+00_EB, 1.589788e+00_EB, 1.374889e+00_EB, 1.152320e+00_EB, &
   9.331464e-01_EB, &
   4.349487e+00_EB, 3.315472e+00_EB, 2.890903e+00_EB, 2.568487e+00_EB, 2.077789e+00_EB, 1.524372e+00_EB, &
   1.129152e+00_EB/),(/7,8/))

sd3_c2h6(1:7,17:24) = RESHAPE((/ &  ! 2950-3125 cm-1
   4.686779e+00_EB, 3.462217e+00_EB, 2.991646e+00_EB, 2.656886e+00_EB, 2.157561e+00_EB, 1.584726e+00_EB, &
   1.167425e+00_EB, &
   5.661792e+00_EB, 4.018616e+00_EB, 3.411886e+00_EB, 2.975545e+00_EB, 2.380752e+00_EB, 1.690561e+00_EB, &
   1.242062e+00_EB, &
   4.918064e+00_EB, 3.441837e+00_EB, 2.903757e+00_EB, 2.525786e+00_EB, 2.015014e+00_EB, 1.429458e+00_EB, &
   1.053831e+00_EB, &
   3.342233e+00_EB, 2.430853e+00_EB, 2.070409e+00_EB, 1.791861e+00_EB, 1.453571e+00_EB, 1.041375e+00_EB, &
   8.022528e-01_EB, &
   1.298165e+00_EB, 1.255308e+00_EB, 1.160176e+00_EB, 1.048498e+00_EB, 9.329619e-01_EB, 7.444407e-01_EB, &
   5.802852e-01_EB, &
   2.926016e-01_EB, 4.494914e-01_EB, 4.813448e-01_EB, 4.684675e-01_EB, 5.012802e-01_EB, 4.679462e-01_EB, &
   4.003903e-01_EB, &
   6.131343e-02_EB, 1.274991e-01_EB, 1.601441e-01_EB, 1.663262e-01_EB, 2.229073e-01_EB, 2.672821e-01_EB, &
   2.381398e-01_EB, &
   4.322732e-02_EB, 6.015505e-02_EB, 6.743728e-02_EB, 6.674353e-02_EB, 1.110929e-01_EB, 1.464795e-01_EB, &
   1.569146e-01_EB/),(/7,8/))

sd3_c2h6(1:7,25:32) = RESHAPE((/ &  ! 3150-3325 cm-1
   3.604160e-02_EB, 1.194845e-01_EB, 4.058106e-02_EB, 4.315323e-02_EB, 9.519510e-02_EB, 8.917098e-02_EB, &
   1.405300e-01_EB, &
   4.315334e-02_EB, 4.712174e-01_EB, 4.312845e-02_EB, 3.219296e-02_EB, 3.525944e-01_EB, 1.060124e-01_EB, &
   5.046058e-01_EB, &
   1.257731e-01_EB, 8.432473e-01_EB, 9.486353e-02_EB, 2.860309e-01_EB, 4.209545e-01_EB, 1.897258e-01_EB, &
   1.586487e+00_EB, &
   5.478996e-02_EB, 5.751983e-01_EB, 1.238009e-01_EB, 1.250768e-01_EB, 9.468278e-01_EB, 3.490872e-01_EB, &
   1.256817e+00_EB, &
   7.050923e-02_EB, 6.595339e-01_EB, 1.143623e-01_EB, 3.234956e-01_EB, 7.050330e-01_EB, 3.260248e-01_EB, &
   5.549485e-01_EB, &
   4.360067e-02_EB, 5.450163e-01_EB, 1.152844e-01_EB, 2.285287e-01_EB, 6.637506e-01_EB, 2.501009e-01_EB, &
   3.535923e-01_EB, &
   1.872369e-01_EB, 3.564558e-01_EB, 4.090535e-02_EB, 1.306462e-01_EB, 4.246963e-01_EB, 7.148910e-02_EB, &
   1.438647e-01_EB, &
   5.015709e-03_EB, 1.912379e-01_EB, 9.285764e-03_EB, 1.392927e-03_EB, 2.040005e-01_EB, 9.061367e-02_EB, &
   9.395617e-04_EB/),(/7,8/))

sd3_c2h6(1:7,33:34) = RESHAPE((/ &  ! 3350-3375 cm-1
   2.644024e-04_EB, 1.032051e-02_EB, 3.153028e-03_EB, 3.073098e-04_EB, 1.028061e-02_EB, 1.242207e-03_EB, &
   2.965961e-04_EB, &
   2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, 2.145139e-04_EB, &
   2.145139e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad3_c2h6(n_temp_c2h6,34)) 

! band #3: 2550 cm-1 - 3375 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.3765 % 

! print fine structure array gamma_d 

gammad3_c2h6(1:7,1:8) = RESHAPE((/ &  ! 2550-2725 cm-1
   5.906363e-04_EB, 2.102426e-04_EB, 2.132346e-04_EB, 2.065179e-04_EB, 2.655299e-04_EB, 8.306361e-04_EB, &
   3.355503e-05_EB, &
   2.516050e-04_EB, 3.355121e-03_EB, 4.172370e-04_EB, 1.011463e-03_EB, 2.596114e-02_EB, 1.182928e-03_EB, &
   2.503340e-05_EB, &
   4.163269e-03_EB, 1.125552e-01_EB, 4.980295e-03_EB, 2.696283e-03_EB, 3.105873e-02_EB, 4.712694e-03_EB, &
   6.033084e-05_EB, &
   2.623691e-02_EB, 1.923864e-02_EB, 3.655389e-02_EB, 4.141023e-03_EB, 3.609314e-02_EB, 3.824795e-03_EB, &
   1.219453e-04_EB, &
   5.265718e-02_EB, 8.436915e-03_EB, 1.960254e-02_EB, 1.875421e-03_EB, 2.111866e-02_EB, 3.776296e-03_EB, &
   3.441962e-02_EB, &
   3.937315e-02_EB, 8.289911e-03_EB, 4.840661e-03_EB, 4.421232e-03_EB, 7.814993e-02_EB, 1.662021e-02_EB, &
   1.455835e-01_EB, &
   3.720731e-02_EB, 8.723819e-03_EB, 7.153107e-03_EB, 7.013544e-03_EB, 7.319968e-02_EB, 1.767570e-01_EB, &
   8.122963e-01_EB, &
   4.868277e-02_EB, 1.269629e-02_EB, 3.095084e-02_EB, 1.285202e-02_EB, 3.311119e-01_EB, 4.745323e-02_EB, &
   5.448600e-01_EB/),(/7,8/))

gammad3_c2h6(1:7,9:16) = RESHAPE((/ &  ! 2750-2925 cm-1
   7.027332e-02_EB, 3.651191e-02_EB, 1.128884e-01_EB, 2.630950e-02_EB, 2.289043e-01_EB, 3.168499e-02_EB, &
   2.263308e-02_EB, &
   1.007340e-01_EB, 9.210111e-02_EB, 1.687337e-01_EB, 2.959829e-02_EB, 1.810977e-01_EB, 9.553982e-02_EB, &
   3.918733e-02_EB, &
   4.509352e-02_EB, 1.086112e-01_EB, 7.647679e-02_EB, 2.814778e-02_EB, 8.981709e-01_EB, 1.273344e-01_EB, &
   4.408137e-02_EB, &
   9.590449e-02_EB, 2.770484e-01_EB, 1.477733e-01_EB, 1.049910e-01_EB, 3.105849e-01_EB, 2.918868e-01_EB, &
   1.291524e-01_EB, &
   3.138946e-01_EB, 5.413508e-01_EB, 3.949461e-01_EB, 3.999787e-01_EB, 6.266594e-01_EB, 4.691785e-01_EB, &
   3.754321e-01_EB, &
   1.106384e+00_EB, 1.162295e+00_EB, 1.094260e+00_EB, 1.051491e+00_EB, 1.215748e+00_EB, 7.154832e-01_EB, &
   4.545049e-01_EB, &
   1.517339e+00_EB, 1.415673e+00_EB, 1.368748e+00_EB, 1.312400e+00_EB, 1.453900e+00_EB, 9.649236e-01_EB, &
   7.807490e-01_EB, &
   4.303118e+00_EB, 2.952375e+00_EB, 2.831306e+00_EB, 2.384206e+00_EB, 2.387748e+00_EB, 1.208513e+00_EB, &
   8.645917e-01_EB/),(/7,8/))

gammad3_c2h6(1:7,17:24) = RESHAPE((/ &  ! 2950-3125 cm-1
   3.315145e+00_EB, 2.599752e+00_EB, 2.743646e+00_EB, 2.629391e+00_EB, 2.506127e+00_EB, 1.333128e+00_EB, &
   9.626906e-01_EB, &
   1.456169e+00_EB, 1.207547e+00_EB, 1.363886e+00_EB, 1.509996e+00_EB, 1.699983e+00_EB, 1.288544e+00_EB, &
   8.734672e-01_EB, &
   1.081840e+00_EB, 1.007445e+00_EB, 1.186416e+00_EB, 1.295045e+00_EB, 1.467006e+00_EB, 1.089437e+00_EB, &
   7.294066e-01_EB, &
   1.769766e+00_EB, 1.424872e+00_EB, 1.588295e+00_EB, 1.745279e+00_EB, 1.453194e+00_EB, 9.836635e-01_EB, &
   4.884795e-01_EB, &
   8.210794e-01_EB, 7.811090e-01_EB, 9.322482e-01_EB, 1.373898e+00_EB, 1.084482e+00_EB, 6.956616e-01_EB, &
   5.118669e-01_EB, &
   1.964022e-01_EB, 3.039042e-01_EB, 3.956607e-01_EB, 1.328063e+00_EB, 5.899273e-01_EB, 6.132437e-01_EB, &
   5.047836e-01_EB, &
   2.041406e-02_EB, 8.312137e-02_EB, 1.318700e-01_EB, 7.573599e+00_EB, 2.512530e-01_EB, 2.496958e-01_EB, &
   1.051791e+00_EB, &
   5.023107e-03_EB, 9.491129e-03_EB, 3.053563e-02_EB, 9.379598e-02_EB, 3.673982e-02_EB, 7.455316e-02_EB, &
   2.454365e-01_EB/),(/7,8/))

gammad3_c2h6(1:7,25:32) = RESHAPE((/ &  ! 3150-3325 cm-1
   6.352743e-03_EB, 7.943354e-04_EB, 1.083748e-02_EB, 9.588981e-03_EB, 4.904702e-03_EB, 1.679226e-02_EB, &
   1.304495e-02_EB, &
   6.766775e-03_EB, 1.596208e-04_EB, 6.130026e-03_EB, 7.152650e-03_EB, 3.279571e-04_EB, 1.794526e-03_EB, &
   4.823723e-04_EB, &
   1.510899e-03_EB, 1.196785e-04_EB, 1.094313e-03_EB, 1.868976e-04_EB, 1.950517e-04_EB, 3.851491e-04_EB, &
   8.135632e-05_EB, &
   1.592258e-02_EB, 2.062254e-04_EB, 8.766044e-04_EB, 4.062596e-04_EB, 1.031744e-04_EB, 1.363841e-04_EB, &
   4.489251e-05_EB, &
   5.336883e-03_EB, 1.621021e-04_EB, 8.122952e-04_EB, 9.938396e-05_EB, 1.057126e-04_EB, 8.283648e-05_EB, &
   2.421683e-05_EB, &
   7.894809e-03_EB, 1.025986e-04_EB, 4.657427e-04_EB, 6.874881e-05_EB, 6.535834e-05_EB, 3.906150e-05_EB, &
   4.299378e-05_EB, &
   1.168996e-04_EB, 6.996908e-05_EB, 5.712580e-04_EB, 3.973347e-05_EB, 4.446145e-05_EB, 8.103887e-05_EB, &
   3.469019e-05_EB, &
   1.835174e-04_EB, 4.182929e-05_EB, 1.449024e-03_EB, 6.981809e-04_EB, 3.965878e-05_EB, 3.732167e-05_EB, &
   1.777260e-02_EB/),(/7,8/))

gammad3_c2h6(1:7,33:34) = RESHAPE((/ &  ! 3350-3375 cm-1
   1.966771e-04_EB, 1.267047e-04_EB, 5.949690e-04_EB, 2.287935e-04_EB, 1.200025e-04_EB, 2.518698e-04_EB, &
   2.116261e-04_EB, &
   1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, 1.688326e-04_EB, &
   1.688326e-04_EB/),(/7,2/))

!-------------------------ethylene data-------------------


! there are 4 bands for ethylene

! band #1: 750 cm-1 - 1250 cm-1 
! band #2: 1300 cm-1 - 1600 cm-1 
! band #3: 1750 cm-1 - 2075 cm-1 
! band #4: 2800 cm-1 - 3400 cm-1 

ALLOCATE(sd_c2h4_temp(n_temp_c2h4)) 

! initialize bands wavenumber bounds for ethylene ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_c2h4(n_band_c2h4,3)) 
ALLOCATE(be_c2h4(n_band_c2h4)) 

om_bnd_c2h4 = RESHAPE((/ &
   750._EB, 1300._EB, 1750._EB, 2800._EB, &
   1250._EB, 1600._EB, 2075._EB, 3400._EB, &
   5._EB, 25._EB, 25._EB, 25._EB/),(/n_band_c2h4,3/)) 

sd_c2h4_temp = (/ &
   296._EB, 400._EB, 450._EB, 500._EB, 601._EB, 801._EB,&
   1000._EB/)

be_c2h4 = (/ &
   1.000_EB, 1.000_EB, 1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_c2h4(n_temp_c2h4,77)) 

! band #1: 750 cm-1 - 1250 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.6588 % 

sd1_c2h4(1:7,1:8) = RESHAPE((/ &  ! 750-785 cm-1
   3.061616e-04_EB, 2.767365e-04_EB, 2.573770e-04_EB, 3.123088e-04_EB, 2.802708e-04_EB, 3.214089e-04_EB, &
   1.645179e-04_EB, &
   3.863930e-04_EB, 3.214740e-04_EB, 3.599651e-04_EB, 2.875226e-04_EB, 3.220027e-04_EB, 3.442021e-04_EB, &
   3.429771e-04_EB, &
   3.479571e-04_EB, 3.988382e-04_EB, 5.410044e-04_EB, 3.711756e-04_EB, 4.049933e-04_EB, 3.780624e-04_EB, &
   1.581304e-04_EB, &
   3.363011e-04_EB, 4.180798e-04_EB, 3.197069e-04_EB, 3.653986e-04_EB, 4.187833e-04_EB, 7.101448e-03_EB, &
   3.809923e-04_EB, &
   3.533603e-04_EB, 3.413308e-04_EB, 3.254621e-04_EB, 4.091800e-04_EB, 3.812469e-04_EB, 2.897062e-01_EB, &
   2.857766e-04_EB, &
   3.471021e-04_EB, 3.981441e-04_EB, 3.391170e-04_EB, 3.788735e-04_EB, 3.733369e-04_EB, 7.154619e-01_EB, &
   2.973615e-04_EB, &
   3.515529e-04_EB, 3.867219e-04_EB, 3.372267e-04_EB, 4.293329e-04_EB, 3.549482e-03_EB, 2.273679e+00_EB, &
   2.916009e-04_EB, &
   4.086122e-04_EB, 3.034723e-04_EB, 3.481328e-04_EB, 4.240499e-04_EB, 2.168522e-02_EB, 8.658210e-01_EB, &
   2.764510e-04_EB/),(/7,8/))

sd1_c2h4(1:7,9:16) = RESHAPE((/ &  ! 790-825 cm-1
   3.370619e-04_EB, 3.449516e-04_EB, 3.068753e-03_EB, 2.166061e-03_EB, 3.067872e-01_EB, 2.126871e+00_EB, &
   3.160891e-04_EB, &
   2.931244e-04_EB, 2.735279e-03_EB, 6.531861e-01_EB, 6.019424e-03_EB, 5.639027e-01_EB, 1.299202e+00_EB, &
   3.396657e-04_EB, &
   4.877063e-04_EB, 9.280827e-03_EB, 1.579979e-02_EB, 1.496897e-02_EB, 4.992397e-01_EB, 1.406058e+00_EB, &
   3.043444e-04_EB, &
   4.001483e-04_EB, 1.891644e-02_EB, 9.491147e-03_EB, 2.267553e-02_EB, 3.421534e-01_EB, 6.382093e-02_EB, &
   1.852766e-04_EB, &
   1.290823e-03_EB, 2.799153e-02_EB, 1.446431e-02_EB, 2.843022e-02_EB, 8.398144e-02_EB, 4.512817e-02_EB, &
   1.294866e-02_EB, &
   2.070200e-01_EB, 3.632895e-02_EB, 3.296335e-02_EB, 4.820045e-02_EB, 1.668313e-01_EB, 1.323383e-01_EB, &
   3.901722e+00_EB, &
   2.001955e-01_EB, 4.074161e-02_EB, 5.504791e-02_EB, 2.033738e-01_EB, 1.370229e-01_EB, 1.403502e+00_EB, &
   9.703462e-01_EB, &
   2.718536e-02_EB, 7.068372e-02_EB, 8.035231e-02_EB, 1.066961e-01_EB, 1.330522e-01_EB, 2.467964e-01_EB, &
   4.709421e+00_EB/),(/7,8/))

sd1_c2h4(1:7,17:24) = RESHAPE((/ &  ! 830-865 cm-1
   8.016587e-02_EB, 1.232301e-01_EB, 1.080896e-01_EB, 1.876693e-01_EB, 1.436722e-01_EB, 1.806775e-01_EB, &
   1.955761e-01_EB, &
   8.134800e-02_EB, 1.769076e-01_EB, 1.680321e-01_EB, 1.930321e-01_EB, 1.757502e-01_EB, 2.398550e-01_EB, &
   2.039868e-01_EB, &
   1.441756e-01_EB, 2.209548e-01_EB, 2.092624e-01_EB, 2.537575e-01_EB, 2.144884e-01_EB, 2.788351e-01_EB, &
   7.848435e-01_EB, &
   1.963952e-01_EB, 2.703514e-01_EB, 2.622330e-01_EB, 3.903483e-01_EB, 3.265592e-01_EB, 3.009724e-01_EB, &
   6.260775e+00_EB, &
   2.991975e-01_EB, 3.352004e-01_EB, 3.564035e-01_EB, 3.796975e-01_EB, 4.082076e-01_EB, 3.402704e-01_EB, &
   2.111753e+00_EB, &
   3.860882e-01_EB, 4.216740e-01_EB, 4.053700e-01_EB, 4.209808e-01_EB, 5.214527e-01_EB, 5.131750e-01_EB, &
   3.472334e-01_EB, &
   4.648316e-01_EB, 4.632385e-01_EB, 4.893760e-01_EB, 5.410936e-01_EB, 5.227714e-01_EB, 3.928599e-01_EB, &
   3.841886e-01_EB, &
   6.022275e-01_EB, 5.990898e-01_EB, 5.570012e-01_EB, 6.192645e-01_EB, 6.317220e-01_EB, 4.457061e-01_EB, &
   4.504311e-01_EB/),(/7,8/))

sd1_c2h4(1:7,25:32) = RESHAPE((/ &  ! 870-905 cm-1
   7.041966e-01_EB, 6.621798e-01_EB, 6.254956e-01_EB, 6.363171e-01_EB, 6.107430e-01_EB, 5.497340e-01_EB, &
   5.235314e-01_EB, &
   8.729129e-01_EB, 7.949380e-01_EB, 7.594193e-01_EB, 7.573864e-01_EB, 6.279590e-01_EB, 5.238851e-01_EB, &
   6.195839e-01_EB, &
   1.076407e+00_EB, 9.480760e-01_EB, 8.612630e-01_EB, 8.823757e-01_EB, 7.150564e-01_EB, 6.412975e-01_EB, &
   5.713233e-01_EB, &
   1.149326e+00_EB, 9.164652e-01_EB, 8.407015e-01_EB, 7.805755e-01_EB, 7.361358e-01_EB, 4.885371e-01_EB, &
   6.246788e-01_EB, &
   1.285145e+00_EB, 1.001830e+00_EB, 9.098842e-01_EB, 8.665704e-01_EB, 7.069087e-01_EB, 5.784402e-01_EB, &
   5.420318e-01_EB, &
   1.656594e+00_EB, 1.243694e+00_EB, 1.132558e+00_EB, 1.004866e+00_EB, 8.573835e-01_EB, 6.721572e-01_EB, &
   5.651920e-01_EB, &
   1.512695e+00_EB, 1.063450e+00_EB, 9.534704e-01_EB, 8.623387e-01_EB, 7.340704e-01_EB, 5.647638e-01_EB, &
   5.138506e-01_EB, &
   1.503391e+00_EB, 1.046223e+00_EB, 9.264262e-01_EB, 8.216305e-01_EB, 7.214293e-01_EB, 4.921282e-01_EB, &
   9.505467e-01_EB/),(/7,8/))

sd1_c2h4(1:7,33:40) = RESHAPE((/ &  ! 910-945 cm-1
   2.136589e+00_EB, 1.441371e+00_EB, 1.265404e+00_EB, 1.131971e+00_EB, 8.629746e-01_EB, 5.974312e-01_EB, &
   6.931755e-01_EB, &
   1.937708e+00_EB, 1.292014e+00_EB, 1.117250e+00_EB, 9.960421e-01_EB, 7.945183e-01_EB, 6.856903e-01_EB, &
   6.607081e-01_EB, &
   1.854239e+00_EB, 1.153175e+00_EB, 9.838961e-01_EB, 9.179450e-01_EB, 6.997797e-01_EB, 5.900011e-01_EB, &
   1.136208e+00_EB, &
   1.823799e+00_EB, 1.165735e+00_EB, 1.007681e+00_EB, 9.292605e-01_EB, 7.896203e-01_EB, 5.208819e-01_EB, &
   1.763502e+00_EB, &
   1.757932e+00_EB, 1.135853e+00_EB, 9.853861e-01_EB, 9.026759e-01_EB, 7.075480e-01_EB, 6.179452e-01_EB, &
   6.802571e-01_EB, &
   1.474799e+00_EB, 1.023497e+00_EB, 9.336669e-01_EB, 8.778123e-01_EB, 8.131435e-01_EB, 6.616974e-01_EB, &
   7.567051e-01_EB, &
   2.289039e+00_EB, 1.607837e+00_EB, 1.440647e+00_EB, 1.347166e+00_EB, 1.239509e+00_EB, 1.038545e+00_EB, &
   1.100592e+00_EB, &
   3.685332e+00_EB, 2.419819e+00_EB, 2.144192e+00_EB, 1.965314e+00_EB, 1.747995e+00_EB, 1.497634e+00_EB, &
   1.667794e+00_EB/),(/7,8/))

sd1_c2h4(1:7,41:48) = RESHAPE((/ &  ! 950-985 cm-1
   9.810246e+00_EB, 7.439967e+00_EB, 6.716823e+00_EB, 6.081792e+00_EB, 4.841636e+00_EB, 3.239145e+00_EB, &
   2.178622e+00_EB, &
   2.747963e+00_EB, 2.185259e+00_EB, 2.110959e+00_EB, 2.104012e+00_EB, 1.967093e+00_EB, 1.887167e+00_EB, &
   1.592410e+00_EB, &
   1.878302e+00_EB, 1.339111e+00_EB, 1.209909e+00_EB, 1.180489e+00_EB, 1.063613e+00_EB, 1.043418e+00_EB, &
   9.802449e-01_EB, &
   1.846689e+00_EB, 1.233103e+00_EB, 1.078269e+00_EB, 1.010819e+00_EB, 7.843025e-01_EB, 7.220401e-01_EB, &
   6.485561e-01_EB, &
   2.135746e+00_EB, 1.274099e+00_EB, 1.072632e+00_EB, 9.574284e-01_EB, 7.244360e-01_EB, 6.000363e-01_EB, &
   5.922618e-01_EB, &
   1.598509e+00_EB, 1.039908e+00_EB, 9.058024e-01_EB, 8.461204e-01_EB, 6.619076e-01_EB, 5.572688e-01_EB, &
   6.188948e-01_EB, &
   2.103202e+00_EB, 1.300976e+00_EB, 1.111634e+00_EB, 9.847661e-01_EB, 7.744278e-01_EB, 5.718997e-01_EB, &
   6.346134e-01_EB, &
   2.261675e+00_EB, 1.454471e+00_EB, 1.243161e+00_EB, 1.093485e+00_EB, 8.306949e-01_EB, 6.753502e-01_EB, &
   5.723775e-01_EB/),(/7,8/))

sd1_c2h4(1:7,49:56) = RESHAPE((/ &  ! 990-1025 cm-1
   2.079467e+00_EB, 1.263897e+00_EB, 1.076642e+00_EB, 9.600844e-01_EB, 7.479525e-01_EB, 5.904677e-01_EB, &
   5.384319e-01_EB, &
   1.966589e+00_EB, 1.324982e+00_EB, 1.148642e+00_EB, 1.030295e+00_EB, 8.285130e-01_EB, 6.028164e-01_EB, &
   5.615524e-01_EB, &
   1.704129e+00_EB, 1.136312e+00_EB, 1.028789e+00_EB, 9.378779e-01_EB, 7.696240e-01_EB, 6.407307e-01_EB, &
   5.740712e-01_EB, &
   2.099061e+00_EB, 1.480171e+00_EB, 1.279590e+00_EB, 1.150617e+00_EB, 8.965259e-01_EB, 6.108346e-01_EB, &
   5.651648e-01_EB, &
   1.532195e+00_EB, 1.072161e+00_EB, 9.586967e-01_EB, 8.739703e-01_EB, 7.409633e-01_EB, 5.463573e-01_EB, &
   5.543106e-01_EB, &
   1.449718e+00_EB, 1.130569e+00_EB, 1.010466e+00_EB, 9.314315e-01_EB, 8.044840e-01_EB, 5.741548e-01_EB, &
   5.103097e-01_EB, &
   1.361629e+00_EB, 1.048216e+00_EB, 9.484550e-01_EB, 8.875882e-01_EB, 7.583217e-01_EB, 6.224803e-01_EB, &
   5.104751e-01_EB, &
   1.215443e+00_EB, 9.557761e-01_EB, 8.676531e-01_EB, 8.218519e-01_EB, 6.763815e-01_EB, 5.340442e-01_EB, &
   4.917932e-01_EB/),(/7,8/))

sd1_c2h4(1:7,57:64) = RESHAPE((/ &  ! 1030-1065 cm-1
   1.134122e+00_EB, 9.472716e-01_EB, 8.758677e-01_EB, 8.193195e-01_EB, 7.086094e-01_EB, 5.766945e-01_EB, &
   6.808968e-01_EB, &
   8.950378e-01_EB, 7.770052e-01_EB, 7.253665e-01_EB, 7.245858e-01_EB, 6.204068e-01_EB, 5.655605e-01_EB, &
   6.332963e-01_EB, &
   6.844087e-01_EB, 6.493525e-01_EB, 6.230393e-01_EB, 6.224226e-01_EB, 5.363550e-01_EB, 4.974590e-01_EB, &
   7.583943e-01_EB, &
   7.556324e-01_EB, 7.049278e-01_EB, 6.726842e-01_EB, 6.598185e-01_EB, 5.886987e-01_EB, 5.424892e-01_EB, &
   5.990280e-01_EB, &
   6.289721e-01_EB, 6.069297e-01_EB, 6.011876e-01_EB, 5.779180e-01_EB, 5.326335e-01_EB, 4.750991e-01_EB, &
   4.588251e-01_EB, &
   4.745725e-01_EB, 4.627839e-01_EB, 4.803008e-01_EB, 4.911406e-01_EB, 4.593009e-01_EB, 3.894394e-01_EB, &
   4.927539e-01_EB, &
   4.243657e-01_EB, 4.444157e-01_EB, 4.520903e-01_EB, 4.522320e-01_EB, 4.526266e-01_EB, 3.758722e-01_EB, &
   4.237227e-01_EB, &
   2.905835e-01_EB, 3.554843e-01_EB, 3.584158e-01_EB, 3.801254e-01_EB, 4.065901e-01_EB, 4.115057e-01_EB, &
   5.026136e-01_EB/),(/7,8/))

sd1_c2h4(1:7,65:72) = RESHAPE((/ &  ! 1070-1125 cm-1
   3.642080e-01_EB, 4.087106e-01_EB, 4.100121e-01_EB, 4.077213e-01_EB, 4.343523e-01_EB, 3.642286e-01_EB, &
   4.076571e-01_EB, &
   2.285227e-01_EB, 2.715118e-01_EB, 2.835998e-01_EB, 3.085575e-01_EB, 3.428083e-01_EB, 2.858531e-01_EB, &
   3.538602e-01_EB, &
   2.072180e-01_EB, 2.508628e-01_EB, 2.687615e-01_EB, 3.061642e-01_EB, 3.126805e-01_EB, 3.288883e-01_EB, &
   3.309101e-01_EB, &
   1.485077e-01_EB, 1.994016e-01_EB, 2.138919e-01_EB, 2.352105e-01_EB, 2.665822e-01_EB, 2.680282e-01_EB, &
   4.173662e-01_EB, &
   1.233539e-01_EB, 1.739268e-01_EB, 2.005570e-01_EB, 2.264850e-01_EB, 2.387994e-01_EB, 2.855644e-01_EB, &
   3.039316e-01_EB, &
   1.208930e-01_EB, 1.604253e-01_EB, 1.886652e-01_EB, 1.954718e-01_EB, 2.106737e-01_EB, 2.480406e-01_EB, &
   2.674441e-01_EB, &
   9.300476e-02_EB, 1.298575e-01_EB, 1.529873e-01_EB, 1.727849e-01_EB, 2.062593e-01_EB, 2.191541e-01_EB, &
   2.754971e-01_EB, &
   2.843120e-02_EB, 4.750139e-02_EB, 6.454680e-02_EB, 7.833683e-02_EB, 1.007079e-01_EB, 1.334723e-01_EB, &
   1.952035e-01_EB/),(/7,8/))

sd1_c2h4(1:7,73:77) = RESHAPE((/ &  ! 1150-1250 cm-1
   1.746127e-02_EB, 7.841172e-03_EB, 2.585562e-02_EB, 3.153723e-02_EB, 3.990568e-02_EB, 6.617133e-02_EB, &
   2.077741e-01_EB, &
   3.428479e-04_EB, 3.716895e-04_EB, 5.481270e-03_EB, 9.380871e-03_EB, 1.423446e-02_EB, 3.592881e-02_EB, &
   5.844943e-01_EB, &
   3.844693e-04_EB, 4.065786e-04_EB, 5.088321e-04_EB, 3.670823e-04_EB, 1.012681e-03_EB, 1.177143e-02_EB, &
   1.441741e+00_EB, &
   2.863454e-04_EB, 2.697429e-04_EB, 2.663683e-04_EB, 2.847322e-04_EB, 2.814940e-04_EB, 4.490652e-04_EB, &
   2.765929e-03_EB, &
   2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, &
   1.035699e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_c2h4(n_temp_c2h4,77)) 

! band #1: 750 cm-1 - 1250 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.6588 % 

! print fine structure array gamma_d 

gammad1_c2h4(1:7,1:8) = RESHAPE((/ &  ! 750-785 cm-1
   2.416099e-04_EB, 1.972502e-04_EB, 2.014288e-04_EB, 2.567319e-04_EB, 2.171052e-04_EB, 2.409977e-04_EB, &
   1.211331e-04_EB, &
   2.468780e-04_EB, 2.281487e-04_EB, 2.453616e-04_EB, 1.924299e-04_EB, 2.048382e-04_EB, 2.401308e-04_EB, &
   2.496562e-04_EB, &
   2.424362e-04_EB, 2.577800e-04_EB, 3.990244e-04_EB, 2.948000e-04_EB, 2.703709e-04_EB, 2.386217e-04_EB, &
   1.097937e-04_EB, &
   2.310838e-04_EB, 3.105661e-04_EB, 2.545394e-04_EB, 2.425587e-04_EB, 2.350566e-04_EB, 5.850300e-04_EB, &
   2.717545e-04_EB, &
   2.188209e-04_EB, 2.694507e-04_EB, 2.395597e-04_EB, 2.852713e-04_EB, 3.015723e-04_EB, 1.785307e-05_EB, &
   1.970987e-04_EB, &
   2.681455e-04_EB, 3.065002e-04_EB, 2.277458e-04_EB, 2.640085e-04_EB, 2.242100e-04_EB, 5.823579e-05_EB, &
   2.107112e-04_EB, &
   2.864384e-04_EB, 2.901989e-04_EB, 2.177055e-04_EB, 2.318611e-04_EB, 8.020695e-03_EB, 6.053782e-05_EB, &
   1.973837e-04_EB, &
   2.864942e-04_EB, 1.888603e-04_EB, 2.876981e-04_EB, 2.918797e-04_EB, 4.025819e-04_EB, 1.500744e-04_EB, &
   1.850902e-04_EB/),(/7,8/))

gammad1_c2h4(1:7,9:16) = RESHAPE((/ &  ! 790-825 cm-1
   2.028540e-04_EB, 2.332145e-04_EB, 8.835388e-04_EB, 1.557804e-03_EB, 6.921471e-05_EB, 1.608572e-04_EB, &
   2.284767e-04_EB, &
   1.964413e-04_EB, 1.354319e-03_EB, 3.099296e-05_EB, 1.230936e-01_EB, 9.177521e-05_EB, 3.735015e-04_EB, &
   2.567094e-04_EB, &
   3.523993e-04_EB, 4.495628e-03_EB, 1.170244e-01_EB, 1.834956e-01_EB, 2.119215e-04_EB, 4.987881e-04_EB, &
   2.324660e-04_EB, &
   2.481552e-04_EB, 3.297622e-03_EB, 3.408731e-01_EB, 1.941870e-02_EB, 5.655524e-04_EB, 6.632397e+00_EB, &
   1.455497e-04_EB, &
   6.200957e-02_EB, 5.552681e-03_EB, 1.304796e+00_EB, 9.167035e-01_EB, 1.036969e-02_EB, 4.543832e+00_EB, &
   4.663762e-02_EB, &
   4.857487e-05_EB, 7.312172e-03_EB, 6.397405e-01_EB, 1.805633e+00_EB, 7.876697e-03_EB, 1.052298e+01_EB, &
   2.200552e-04_EB, &
   1.266088e-04_EB, 3.671828e-01_EB, 1.052101e+00_EB, 4.874822e-03_EB, 1.186820e-02_EB, 1.436321e-03_EB, &
   2.304645e-03_EB, &
   7.234540e-02_EB, 1.423099e-01_EB, 4.250396e+00_EB, 4.471836e-01_EB, 1.394765e-01_EB, 2.304974e-02_EB, &
   9.320014e-04_EB/),(/7,8/))

gammad1_c2h4(1:7,17:24) = RESHAPE((/ &  ! 830-865 cm-1
   5.320200e-03_EB, 3.524362e-02_EB, 1.338375e+01_EB, 2.690068e-02_EB, 2.237749e+01_EB, 3.550085e+01_EB, &
   4.352470e+01_EB, &
   1.362377e+00_EB, 4.799018e-02_EB, 3.645279e-01_EB, 1.401350e-01_EB, 2.192962e+01_EB, 1.052708e-01_EB, &
   4.352470e+01_EB, &
   4.653039e-02_EB, 8.421067e-02_EB, 1.759985e+00_EB, 1.658517e-01_EB, 5.613535e+01_EB, 4.860610e+01_EB, &
   9.812101e-03_EB, &
   2.321729e-01_EB, 2.414939e-01_EB, 1.567072e+01_EB, 7.788555e-02_EB, 1.585202e-01_EB, 2.057238e+01_EB, &
   2.123765e-03_EB, &
   1.555629e-01_EB, 8.669673e-01_EB, 5.965418e-01_EB, 3.491496e-01_EB, 1.276794e-01_EB, 4.796946e-01_EB, &
   7.827500e-03_EB, &
   2.604872e-01_EB, 5.734264e-01_EB, 3.245414e+01_EB, 2.473726e+00_EB, 1.188363e-01_EB, 8.314367e-02_EB, &
   4.352145e+01_EB, &
   3.170365e-01_EB, 1.356303e+00_EB, 6.470392e-01_EB, 2.532365e-01_EB, 1.744398e-01_EB, 7.021255e-01_EB, &
   4.352322e+01_EB, &
   5.909101e-01_EB, 1.099398e+00_EB, 4.090840e+00_EB, 4.223559e-01_EB, 2.001671e-01_EB, 4.828266e+01_EB, &
   4.352310e+01_EB/),(/7,8/))

gammad1_c2h4(1:7,25:32) = RESHAPE((/ &  ! 870-905 cm-1
   6.696240e-01_EB, 9.061473e-01_EB, 1.096877e+00_EB, 6.942716e-01_EB, 5.833627e-01_EB, 7.942786e-01_EB, &
   3.236155e-01_EB, &
   7.130009e-01_EB, 9.362301e-01_EB, 7.015500e-01_EB, 5.223698e-01_EB, 6.248457e+00_EB, 4.863162e+01_EB, &
   1.544108e-01_EB, &
   1.043070e+00_EB, 7.601537e-01_EB, 8.543522e-01_EB, 5.009053e-01_EB, 8.115324e-01_EB, 2.773377e-01_EB, &
   1.734668e+00_EB, &
   9.830346e-01_EB, 1.107606e+00_EB, 7.646673e-01_EB, 7.719984e-01_EB, 3.261238e-01_EB, 2.851084e+00_EB, &
   2.517361e-01_EB, &
   1.198560e+00_EB, 1.081208e+00_EB, 7.434343e-01_EB, 6.010115e-01_EB, 9.163865e-01_EB, 3.554812e-01_EB, &
   2.317131e-01_EB, &
   9.814423e-01_EB, 1.013604e+00_EB, 6.693517e-01_EB, 9.522875e-01_EB, 8.114383e-01_EB, 3.910559e-01_EB, &
   3.337523e+01_EB, &
   7.620617e-01_EB, 7.803621e-01_EB, 5.747705e-01_EB, 6.345028e-01_EB, 5.888869e-01_EB, 4.300627e-01_EB, &
   4.310793e+01_EB, &
   8.991018e-01_EB, 8.909633e-01_EB, 6.980425e-01_EB, 1.113145e+00_EB, 7.124349e-01_EB, 4.862871e+01_EB, &
   6.096945e-02_EB/),(/7,8/))

gammad1_c2h4(1:7,33:40) = RESHAPE((/ &  ! 910-945 cm-1
   5.715478e-01_EB, 5.814464e-01_EB, 5.080659e-01_EB, 4.983276e-01_EB, 8.896438e-01_EB, 3.707252e+00_EB, &
   2.418005e-01_EB, &
   4.564794e-01_EB, 5.408890e-01_EB, 5.281411e-01_EB, 6.161207e-01_EB, 9.304299e-01_EB, 1.822560e-01_EB, &
   1.433739e-01_EB, &
   4.474194e-01_EB, 5.933999e-01_EB, 5.926916e-01_EB, 4.188992e-01_EB, 8.916726e-01_EB, 2.166682e-01_EB, &
   3.172885e-02_EB, &
   6.679241e-01_EB, 7.624342e-01_EB, 6.139160e-01_EB, 4.989378e-01_EB, 3.780157e-01_EB, 4.858759e+01_EB, &
   2.070389e-02_EB, &
   1.572323e+00_EB, 1.507380e+00_EB, 1.025688e+00_EB, 7.931157e-01_EB, 2.355979e+00_EB, 4.083523e-01_EB, &
   4.042678e-01_EB, &
   1.840928e+00_EB, 1.516332e+00_EB, 9.674401e-01_EB, 7.614660e-01_EB, 5.264076e-01_EB, 2.038597e+00_EB, &
   6.604494e-01_EB, &
   2.359084e+00_EB, 1.464881e+00_EB, 1.174436e+00_EB, 1.055262e+00_EB, 8.221130e-01_EB, 1.578212e+00_EB, &
   5.803317e-01_EB, &
   1.353327e+00_EB, 1.007304e+00_EB, 8.550921e-01_EB, 9.242096e-01_EB, 8.460605e-01_EB, 1.145876e+00_EB, &
   3.951015e-01_EB/),(/7,8/))

gammad1_c2h4(1:7,41:48) = RESHAPE((/ &  ! 950-985 cm-1
   1.026068e+00_EB, 7.509719e-01_EB, 6.828323e-01_EB, 7.181245e-01_EB, 8.830267e-01_EB, 1.535879e+00_EB, &
   3.739161e+00_EB, &
   1.923367e+00_EB, 1.440812e+00_EB, 1.103009e+00_EB, 1.030902e+00_EB, 1.337752e+00_EB, 1.153093e+00_EB, &
   2.040631e+00_EB, &
   1.610198e+00_EB, 1.196060e+00_EB, 1.104793e+00_EB, 8.780307e-01_EB, 1.156102e+00_EB, 6.886495e-01_EB, &
   2.238055e+00_EB, &
   1.383929e+00_EB, 1.280808e+00_EB, 1.097115e+00_EB, 7.326487e-01_EB, 5.498608e+01_EB, 4.452870e-01_EB, &
   4.352469e+01_EB, &
   1.188786e+00_EB, 1.159934e+00_EB, 9.631713e-01_EB, 7.953825e-01_EB, 1.129330e+01_EB, 6.701708e-01_EB, &
   5.733374e-01_EB, &
   1.242802e+00_EB, 1.350551e+00_EB, 1.172729e+00_EB, 7.385986e-01_EB, 3.096362e+01_EB, 8.122837e-01_EB, &
   2.708279e-01_EB, &
   6.907316e-01_EB, 7.891104e-01_EB, 7.404165e-01_EB, 7.153440e-01_EB, 1.513590e+00_EB, 8.412093e-01_EB, &
   2.566696e-01_EB, &
   7.900363e-01_EB, 8.478862e-01_EB, 7.713692e-01_EB, 7.316301e-01_EB, 1.686388e+00_EB, 3.127574e-01_EB, &
   9.749592e-01_EB/),(/7,8/))

gammad1_c2h4(1:7,49:56) = RESHAPE((/ &  ! 990-1025 cm-1
   5.092459e-01_EB, 6.454904e-01_EB, 6.133995e-01_EB, 5.823764e-01_EB, 1.459319e+00_EB, 3.862803e-01_EB, &
   4.351923e+01_EB, &
   7.621436e-01_EB, 8.099637e-01_EB, 7.550770e-01_EB, 6.563547e-01_EB, 8.599084e-01_EB, 1.265807e+00_EB, &
   4.352461e+01_EB, &
   8.135765e-01_EB, 9.399794e-01_EB, 6.891443e-01_EB, 6.724564e-01_EB, 1.443782e+00_EB, 4.675502e-01_EB, &
   6.832951e+00_EB, &
   1.261934e+00_EB, 1.187896e+00_EB, 1.050548e+00_EB, 8.976709e-01_EB, 1.264652e+00_EB, 4.862525e+01_EB, &
   8.286417e+00_EB, &
   8.043161e-01_EB, 8.668835e-01_EB, 7.890285e-01_EB, 7.711386e-01_EB, 1.050733e+00_EB, 4.856478e+01_EB, &
   4.350007e+01_EB, &
   1.357377e+00_EB, 1.080455e+00_EB, 1.033605e+00_EB, 1.033196e+00_EB, 8.256082e-01_EB, 8.647109e+00_EB, &
   4.350847e+01_EB, &
   8.412833e-01_EB, 8.456811e-01_EB, 8.447705e-01_EB, 7.603788e-01_EB, 9.443837e-01_EB, 2.903532e-01_EB, &
   4.352240e+01_EB, &
   9.318939e-01_EB, 8.628581e-01_EB, 7.702967e-01_EB, 5.812089e-01_EB, 1.478510e+00_EB, 7.524972e-01_EB, &
   4.352454e+01_EB/),(/7,8/))

gammad1_c2h4(1:7,57:64) = RESHAPE((/ &  ! 1030-1065 cm-1
   8.629044e-01_EB, 8.975496e-01_EB, 7.798708e-01_EB, 6.608358e-01_EB, 6.759715e-01_EB, 2.711107e-01_EB, &
   1.220647e-01_EB, &
   7.370489e-01_EB, 8.019656e-01_EB, 1.000196e+00_EB, 5.336084e-01_EB, 1.804552e+00_EB, 2.787416e-01_EB, &
   1.281215e-01_EB, &
   7.787285e-01_EB, 7.334174e-01_EB, 8.721264e-01_EB, 5.887939e-01_EB, 3.552536e+01_EB, 8.474746e-01_EB, &
   8.392367e-02_EB, &
   5.450501e-01_EB, 6.247731e-01_EB, 7.400770e-01_EB, 5.561437e-01_EB, 9.526879e-01_EB, 3.143742e-01_EB, &
   1.138369e-01_EB, &
   5.144688e-01_EB, 6.585084e-01_EB, 5.903669e-01_EB, 7.810567e-01_EB, 8.401970e-01_EB, 6.417739e-01_EB, &
   4.420604e-01_EB, &
   3.665531e-01_EB, 8.719562e-01_EB, 5.223315e-01_EB, 4.375375e-01_EB, 5.743448e-01_EB, 2.121049e+00_EB, &
   2.496608e-01_EB, &
   4.667181e-01_EB, 7.081388e-01_EB, 6.911857e-01_EB, 6.081257e-01_EB, 3.737837e-01_EB, 2.130304e+00_EB, &
   6.904458e-01_EB, &
   3.306129e-01_EB, 3.244247e-01_EB, 1.436316e+00_EB, 6.770932e-01_EB, 3.549036e-01_EB, 2.672621e-01_EB, &
   1.358031e-01_EB/),(/7,8/))

gammad1_c2h4(1:7,65:72) = RESHAPE((/ &  ! 1070-1125 cm-1
   2.919422e-01_EB, 3.828297e-01_EB, 5.395788e-01_EB, 5.826860e-01_EB, 2.566111e-01_EB, 2.715471e+01_EB, &
   2.745597e-01_EB, &
   3.398729e-01_EB, 4.374885e-01_EB, 1.141266e+00_EB, 3.939052e-01_EB, 1.637264e-01_EB, 4.863163e+01_EB, &
   4.352422e+01_EB, &
   1.420571e-01_EB, 4.454635e-01_EB, 1.365866e+00_EB, 2.363442e-01_EB, 2.876851e-01_EB, 4.033493e-01_EB, &
   4.352427e+01_EB, &
   2.987867e-01_EB, 2.946281e-01_EB, 3.722689e+00_EB, 5.611341e-01_EB, 3.323989e-01_EB, 4.184566e+01_EB, &
   8.215267e-02_EB, &
   6.464375e-02_EB, 2.415560e-01_EB, 5.270835e-01_EB, 2.591978e-01_EB, 1.623247e+00_EB, 2.083192e-01_EB, &
   2.703531e+00_EB, &
   1.652449e-01_EB, 7.624040e-01_EB, 3.783911e-01_EB, 6.227587e-01_EB, 4.813467e+00_EB, 4.088673e-01_EB, &
   4.352470e+01_EB, &
   4.426509e-02_EB, 2.377002e-01_EB, 5.194564e+00_EB, 4.416648e-01_EB, 1.768354e-01_EB, 6.614723e-01_EB, &
   8.659489e-01_EB, &
   1.448604e-02_EB, 6.642067e-02_EB, 6.463703e-01_EB, 1.052108e+00_EB, 1.675080e-01_EB, 3.912344e-01_EB, &
   3.760628e-01_EB/),(/7,8/))

gammad1_c2h4(1:7,73:77) = RESHAPE((/ &  ! 1150-1250 cm-1
   2.129075e-04_EB, 2.822799e-02_EB, 2.461798e-02_EB, 8.986262e-02_EB, 1.247997e+00_EB, 1.467934e-01_EB, &
   1.669303e-02_EB, &
   2.363189e-04_EB, 2.307838e-04_EB, 6.030981e-03_EB, 4.177046e-03_EB, 8.698422e-02_EB, 1.673678e-02_EB, &
   1.814471e-03_EB, &
   2.536687e-04_EB, 2.591352e-04_EB, 3.402664e-04_EB, 2.425216e-04_EB, 8.814904e-04_EB, 1.309772e-03_EB, &
   1.809337e-04_EB, &
   2.020888e-04_EB, 2.044892e-04_EB, 1.950677e-04_EB, 2.120114e-04_EB, 2.042227e-04_EB, 3.416414e-04_EB, &
   5.953816e-02_EB, &
   1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, &
   8.125358e-05_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_c2h4(n_temp_c2h4,13)) 

! band #2: 1300 cm-1 - 1600 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.84783 % 

sd2_c2h4(1:7,1:8) = RESHAPE((/ &  ! 1300-1475 cm-1
   2.882468e-04_EB, 2.872390e-04_EB, 2.844740e-04_EB, 2.830853e-04_EB, 2.869806e-04_EB, 2.911547e-04_EB, &
   1.572859e-04_EB, &
   3.652126e-04_EB, 3.386829e-04_EB, 3.628235e-04_EB, 3.693377e-04_EB, 3.974359e-04_EB, 5.740813e-04_EB, &
   1.406431e+00_EB, &
   3.864016e-04_EB, 3.824772e-04_EB, 4.379786e-03_EB, 2.882465e-03_EB, 2.712238e-02_EB, 1.918484e-02_EB, &
   7.860096e-01_EB, &
   1.282101e-02_EB, 2.679804e-02_EB, 4.946063e-02_EB, 5.785843e-02_EB, 7.524531e-02_EB, 1.070631e-01_EB, &
   1.846526e-01_EB, &
   2.572462e-01_EB, 2.610537e-01_EB, 2.762334e-01_EB, 2.568286e-01_EB, 2.524152e-01_EB, 2.204080e-01_EB, &
   1.917598e-01_EB, &
   5.165021e-01_EB, 3.183373e-01_EB, 2.741011e-01_EB, 2.245903e-01_EB, 1.734576e-01_EB, 1.131539e-01_EB, &
   1.085824e-01_EB, &
   4.620901e-01_EB, 3.014428e-01_EB, 2.695816e-01_EB, 2.383870e-01_EB, 1.998293e-01_EB, 1.716093e-01_EB, &
   1.585133e-01_EB, &
   3.857202e-01_EB, 2.649185e-01_EB, 2.402621e-01_EB, 2.174361e-01_EB, 1.792406e-01_EB, 1.252422e-01_EB, &
   1.649973e-01_EB/),(/7,8/))

sd2_c2h4(1:7,9:13) = RESHAPE((/ &  ! 1500-1600 cm-1
   5.926638e-02_EB, 7.277371e-02_EB, 8.111314e-02_EB, 7.867862e-02_EB, 9.699567e-02_EB, 8.865234e-02_EB, &
   7.947370e-02_EB, &
   1.150376e-03_EB, 2.245129e-03_EB, 7.274074e-03_EB, 7.462653e-03_EB, 2.221317e-02_EB, 4.638888e-02_EB, &
   2.712345e-02_EB, &
   3.976051e-04_EB, 3.647587e-04_EB, 3.852124e-04_EB, 3.887963e-04_EB, 3.797778e-04_EB, 1.070583e-03_EB, &
   3.011592e-04_EB, &
   3.339180e-04_EB, 3.021121e-04_EB, 2.677955e-04_EB, 2.991550e-04_EB, 2.786140e-04_EB, 3.017892e-04_EB, &
   1.540746e-04_EB, &
   2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, &
   1.035699e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_c2h4(n_temp_c2h4,13)) 

! band #2: 1300 cm-1 - 1600 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.84783 % 

! print fine structure array gamma_d 

gammad2_c2h4(1:7,1:8) = RESHAPE((/ &  ! 1300-1475 cm-1
   2.158499e-04_EB, 2.198834e-04_EB, 2.077167e-04_EB, 2.022446e-04_EB, 2.019156e-04_EB, 2.108599e-04_EB, &
   1.180395e-04_EB, &
   2.584797e-04_EB, 2.404707e-04_EB, 2.365446e-04_EB, 2.427590e-04_EB, 2.779223e-04_EB, 3.857109e-04_EB, &
   4.413016e-05_EB, &
   2.651063e-04_EB, 2.657195e-04_EB, 1.443449e-03_EB, 2.228190e-03_EB, 2.843764e-04_EB, 1.574315e-02_EB, &
   3.976791e-04_EB, &
   3.442081e-02_EB, 6.159804e-02_EB, 2.276418e-02_EB, 3.598936e-02_EB, 1.165571e-01_EB, 3.704352e-01_EB, &
   3.636410e-02_EB, &
   4.705384e-01_EB, 6.099726e-01_EB, 2.288867e-01_EB, 5.066479e-01_EB, 3.305657e-01_EB, 3.501155e-01_EB, &
   1.491325e-01_EB, &
   9.295562e-01_EB, 6.218514e-01_EB, 2.607721e-01_EB, 5.342728e-01_EB, 6.235127e-01_EB, 1.582834e+00_EB, &
   1.082596e-01_EB, &
   4.753247e-01_EB, 4.602367e-01_EB, 2.660304e-01_EB, 3.252024e-01_EB, 4.963318e-01_EB, 1.168146e-01_EB, &
   1.189208e-01_EB, &
   5.580701e-01_EB, 5.558061e-01_EB, 2.447941e-01_EB, 1.720602e-01_EB, 1.909313e-01_EB, 2.351203e+00_EB, &
   2.093321e-02_EB/),(/7,8/))

gammad2_c2h4(1:7,9:13) = RESHAPE((/ &  ! 1500-1600 cm-1
   4.272966e-01_EB, 1.013776e+00_EB, 9.918181e-02_EB, 1.052649e+00_EB, 4.393797e-02_EB, 3.996882e-02_EB, &
   9.062687e-02_EB, &
   2.701330e-03_EB, 2.355961e-02_EB, 2.364569e-03_EB, 1.244040e-02_EB, 4.075644e-03_EB, 1.508431e-03_EB, &
   3.184016e-01_EB, &
   2.534164e-04_EB, 2.480212e-04_EB, 2.545401e-04_EB, 2.448742e-04_EB, 2.630816e-04_EB, 2.165252e-04_EB, &
   2.247869e-04_EB, &
   2.330947e-04_EB, 2.194445e-04_EB, 2.014090e-04_EB, 2.218115e-04_EB, 2.076602e-04_EB, 2.253786e-04_EB, &
   1.141333e-04_EB, &
   1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, &
   8.125358e-05_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(sd3_c2h4(n_temp_c2h4,14)) 

! band #3: 1750 cm-1 - 2075 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.45198 % 

sd3_c2h4(1:7,1:8) = RESHAPE((/ &  ! 1750-1925 cm-1
   2.815672e-04_EB, 3.085650e-04_EB, 3.011997e-04_EB, 2.507052e-04_EB, 3.301713e-04_EB, 2.824203e-04_EB, &
   1.401316e-04_EB, &
   3.584064e-04_EB, 3.532685e-04_EB, 3.752903e-04_EB, 3.459819e-04_EB, 6.155777e-04_EB, 3.497840e-03_EB, &
   1.410505e-01_EB, &
   3.713646e-04_EB, 3.521998e-04_EB, 5.743608e-04_EB, 8.300439e-04_EB, 2.211093e-02_EB, 1.662667e-02_EB, &
   7.102650e-01_EB, &
   1.552528e-02_EB, 1.444417e-02_EB, 2.574225e-02_EB, 2.639928e-02_EB, 4.136403e-02_EB, 4.346707e-02_EB, &
   8.056475e-01_EB, &
   1.377934e-01_EB, 1.025802e-01_EB, 1.048520e-01_EB, 9.315740e-02_EB, 9.685325e-02_EB, 7.199117e-02_EB, &
   6.245123e-01_EB, &
   1.901714e-01_EB, 1.029944e-01_EB, 9.811801e-02_EB, 8.417027e-02_EB, 7.035037e-02_EB, 5.536099e-02_EB, &
   4.099564e-01_EB, &
   2.542475e-01_EB, 1.411985e-01_EB, 1.285533e-01_EB, 1.090056e-01_EB, 9.009588e-02_EB, 7.442026e-02_EB, &
   1.611310e-01_EB, &
   1.937759e-01_EB, 1.535567e-01_EB, 1.503998e-01_EB, 1.340973e-01_EB, 1.208789e-01_EB, 9.651057e-02_EB, &
   1.654443e-01_EB/),(/7,8/))

sd3_c2h4(1:7,9:14) = RESHAPE((/ &  ! 1950-2075 cm-1
   1.548981e-02_EB, 2.175487e-02_EB, 3.117281e-02_EB, 3.251084e-02_EB, 4.215815e-02_EB, 5.301204e-02_EB, &
   1.143451e-01_EB, &
   4.016633e-04_EB, 3.491488e-04_EB, 4.867031e-04_EB, 4.500772e-04_EB, 2.320822e-03_EB, 7.486336e-03_EB, &
   2.340736e-01_EB, &
   3.822381e-04_EB, 3.883093e-04_EB, 3.278614e-04_EB, 3.723909e-04_EB, 3.584776e-04_EB, 3.541946e-04_EB, &
   2.753666e-01_EB, &
   3.861390e-04_EB, 3.578653e-04_EB, 3.623971e-04_EB, 3.869958e-04_EB, 3.611574e-04_EB, 3.348041e-04_EB, &
   5.612807e-03_EB, &
   2.666887e-04_EB, 2.879471e-04_EB, 2.979300e-04_EB, 2.891663e-04_EB, 3.019136e-04_EB, 2.928543e-04_EB, &
   1.433099e-04_EB, &
   2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, &
   1.035699e-04_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(gammad3_c2h4(n_temp_c2h4,14)) 

! band #3: 1750 cm-1 - 2075 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.45198 % 

! print fine structure array gamma_d 

gammad3_c2h4(1:7,1:8) = RESHAPE((/ &  ! 1750-1925 cm-1
   2.062935e-04_EB, 2.213095e-04_EB, 2.125965e-04_EB, 1.924384e-04_EB, 2.245648e-04_EB, 2.001139e-04_EB, &
   1.052265e-04_EB, &
   2.525656e-04_EB, 2.298492e-04_EB, 2.561819e-04_EB, 2.386624e-04_EB, 3.763751e-04_EB, 6.701805e-04_EB, &
   3.650678e-05_EB, &
   2.518994e-04_EB, 2.486015e-04_EB, 3.999398e-04_EB, 6.636791e-04_EB, 2.149851e-04_EB, 1.018456e-02_EB, &
   1.523459e-04_EB, &
   9.583444e-03_EB, 5.252649e-02_EB, 2.359414e-02_EB, 1.373517e-01_EB, 1.580336e-02_EB, 3.863765e-01_EB, &
   4.400337e-04_EB, &
   1.817495e-01_EB, 1.535393e+00_EB, 1.403613e-01_EB, 5.133053e-01_EB, 5.338970e-02_EB, 9.422243e-01_EB, &
   8.608663e-04_EB, &
   1.604957e-01_EB, 1.152431e+00_EB, 9.900652e-02_EB, 1.034304e-01_EB, 1.658726e-01_EB, 9.190427e-01_EB, &
   1.104201e-03_EB, &
   1.875157e-01_EB, 6.534721e-01_EB, 1.215243e-01_EB, 1.499339e-01_EB, 1.669372e-01_EB, 1.540611e-01_EB, &
   5.452210e-03_EB, &
   2.586428e-01_EB, 3.173253e-01_EB, 1.281338e-01_EB, 2.222569e-01_EB, 1.791652e-01_EB, 9.692458e-01_EB, &
   1.136010e-02_EB/),(/7,8/))

gammad3_c2h4(1:7,9:14) = RESHAPE((/ &  ! 1950-2075 cm-1
   1.624541e-02_EB, 2.116364e-02_EB, 2.358393e-02_EB, 5.479449e-02_EB, 3.978462e-02_EB, 5.533594e-02_EB, &
   7.054981e-03_EB, &
   2.704645e-04_EB, 2.392112e-04_EB, 3.409077e-04_EB, 3.457806e-04_EB, 1.714364e-03_EB, 1.531894e-03_EB, &
   1.624673e-04_EB, &
   2.502376e-04_EB, 2.651107e-04_EB, 2.268500e-04_EB, 2.246834e-04_EB, 2.477431e-04_EB, 2.372743e-04_EB, &
   1.590970e-05_EB, &
   2.596193e-04_EB, 2.512472e-04_EB, 2.549796e-04_EB, 2.535211e-04_EB, 2.549967e-04_EB, 2.385259e-04_EB, &
   1.055504e-04_EB, &
   1.974440e-04_EB, 2.070414e-04_EB, 2.178990e-04_EB, 2.070114e-04_EB, 2.210429e-04_EB, 2.175134e-04_EB, &
   1.061814e-04_EB, &
   1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, &
   8.125358e-05_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(sd4_c2h4(n_temp_c2h4,25)) 

! band #4: 2800 cm-1 - 3400 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.39903 % 

sd4_c2h4(1:7,1:8) = RESHAPE((/ &  ! 2800-2975 cm-1
   2.787454e-04_EB, 2.798900e-04_EB, 3.048047e-04_EB, 2.953254e-04_EB, 3.008118e-04_EB, 2.854943e-04_EB, &
   1.499675e-04_EB, &
   1.974546e-03_EB, 3.576021e-04_EB, 3.665839e-04_EB, 3.719708e-04_EB, 4.008233e-04_EB, 3.722047e-04_EB, &
   3.066454e-04_EB, &
   4.675051e-03_EB, 3.770033e-04_EB, 1.005255e-03_EB, 5.231568e-04_EB, 3.365147e-03_EB, 4.031952e-03_EB, &
   1.001045e-02_EB, &
   1.579559e-02_EB, 4.036503e-04_EB, 9.793374e-03_EB, 4.364856e-03_EB, 1.125997e-02_EB, 1.687883e-02_EB, &
   4.109982e-02_EB, &
   1.901642e-02_EB, 7.948388e-03_EB, 2.110029e-02_EB, 1.680041e-02_EB, 3.681206e-02_EB, 5.911015e-02_EB, &
   9.045692e-02_EB, &
   7.057398e-02_EB, 8.615787e-02_EB, 1.147527e-01_EB, 1.077993e-01_EB, 1.296121e-01_EB, 1.341954e-01_EB, &
   1.504629e-01_EB, &
   4.303610e-01_EB, 3.433046e-01_EB, 3.232359e-01_EB, 2.803463e-01_EB, 2.416656e-01_EB, 1.843691e-01_EB, &
   1.701492e-01_EB, &
   6.023500e-01_EB, 3.794196e-01_EB, 3.422745e-01_EB, 2.887488e-01_EB, 2.637113e-01_EB, 2.107903e-01_EB, &
   2.126134e-01_EB/),(/7,8/))

sd4_c2h4(1:7,9:16) = RESHAPE((/ &  ! 3000-3175 cm-1
   7.798028e-01_EB, 4.757862e-01_EB, 3.974023e-01_EB, 3.382752e-01_EB, 2.973165e-01_EB, 2.380506e-01_EB, &
   2.451699e-01_EB, &
   7.172871e-01_EB, 5.831037e-01_EB, 5.314972e-01_EB, 4.800679e-01_EB, 4.358886e-01_EB, 3.600701e-01_EB, &
   3.346981e-01_EB, &
   5.311509e-01_EB, 4.669769e-01_EB, 4.432862e-01_EB, 4.085659e-01_EB, 3.930117e-01_EB, 3.358462e-01_EB, &
   2.929137e-01_EB, &
   8.649072e-01_EB, 5.873164e-01_EB, 4.935644e-01_EB, 4.177742e-01_EB, 3.508014e-01_EB, 2.332914e-01_EB, &
   2.128007e-01_EB, &
   7.478903e-01_EB, 4.640090e-01_EB, 3.775139e-01_EB, 3.190243e-01_EB, 2.690863e-01_EB, 1.975084e-01_EB, &
   1.921879e-01_EB, &
   1.077261e+00_EB, 6.909648e-01_EB, 5.870020e-01_EB, 4.958355e-01_EB, 4.101225e-01_EB, 3.093396e-01_EB, &
   2.599685e-01_EB, &
   7.637944e-01_EB, 5.720239e-01_EB, 5.128777e-01_EB, 4.357737e-01_EB, 3.822630e-01_EB, 3.059218e-01_EB, &
   2.471299e-01_EB, &
   3.843373e-01_EB, 3.060968e-01_EB, 2.831408e-01_EB, 2.396863e-01_EB, 2.332271e-01_EB, 1.957615e-01_EB, &
   1.637279e-01_EB/),(/7,8/))

sd4_c2h4(1:7,17:24) = RESHAPE((/ &  ! 3200-3375 cm-1
   2.099610e-01_EB, 1.791031e-01_EB, 1.774091e-01_EB, 1.486224e-01_EB, 1.585195e-01_EB, 1.392503e-01_EB, &
   1.207946e-01_EB, &
   1.199548e-01_EB, 1.020644e-01_EB, 9.502921e-02_EB, 8.145885e-02_EB, 1.125359e-01_EB, 1.424681e-01_EB, &
   8.669709e-02_EB, &
   4.616169e-01_EB, 7.858449e-02_EB, 7.272250e-02_EB, 2.867532e-02_EB, 1.168091e-01_EB, 1.330542e-01_EB, &
   5.570436e-02_EB, &
   8.691472e-02_EB, 9.220433e-03_EB, 1.896226e-02_EB, 3.974691e-03_EB, 3.169841e-01_EB, 2.218115e-02_EB, &
   5.082209e-02_EB, &
   2.867450e-02_EB, 3.903775e-04_EB, 3.229916e-04_EB, 3.302777e-04_EB, 1.149402e-01_EB, 6.959691e-02_EB, &
   1.770145e-02_EB, &
   3.690812e-03_EB, 3.360817e-04_EB, 3.487733e-04_EB, 4.020878e-04_EB, 3.692208e-04_EB, 3.501982e-04_EB, &
   8.401837e-03_EB, &
   3.318610e-04_EB, 3.541529e-04_EB, 3.588894e-04_EB, 3.595732e-04_EB, 3.703241e-04_EB, 3.649087e-04_EB, &
   2.903059e-04_EB, &
   3.061877e-04_EB, 2.641235e-04_EB, 2.931529e-04_EB, 2.985364e-04_EB, 2.827758e-04_EB, 2.741559e-04_EB, &
   1.355576e-04_EB/),(/7,8/))

sd4_c2h4(1:7,25:25) = RESHAPE((/ &  ! 3400-3400 cm-1
   2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, 2.172015e-04_EB, &
   1.035699e-04_EB/),(/7,1/))

!---------------------------------------------------------------------------
ALLOCATE(gammad4_c2h4(n_temp_c2h4,25)) 

! band #4: 2800 cm-1 - 3400 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.39903 % 

! print fine structure array gamma_d 

gammad4_c2h4(1:7,1:8) = RESHAPE((/ &  ! 2800-2975 cm-1
   2.053802e-04_EB, 2.055037e-04_EB, 2.250544e-04_EB, 2.225113e-04_EB, 2.136803e-04_EB, 2.087893e-04_EB, &
   1.105053e-04_EB, &
   4.253229e-04_EB, 2.343264e-04_EB, 2.233710e-04_EB, 2.452992e-04_EB, 2.751557e-04_EB, 2.508501e-04_EB, &
   2.222309e-04_EB, &
   1.276922e-03_EB, 2.646545e-04_EB, 7.364197e-04_EB, 3.348640e-04_EB, 3.561604e-03_EB, 2.031177e-03_EB, &
   2.949697e-01_EB, &
   1.385628e-02_EB, 2.885267e-04_EB, 9.042876e-01_EB, 2.692734e-03_EB, 9.713740e-03_EB, 4.653156e-03_EB, &
   2.653808e+00_EB, &
   2.397896e-02_EB, 8.918801e-03_EB, 4.080122e-01_EB, 2.625635e-01_EB, 3.017722e-02_EB, 6.429571e-02_EB, &
   6.032175e+00_EB, &
   9.755870e-02_EB, 2.422956e+00_EB, 1.631270e+01_EB, 6.661998e+00_EB, 3.034249e-01_EB, 5.074180e+00_EB, &
   6.298484e+00_EB, &
   4.898871e-01_EB, 2.215583e+00_EB, 6.480333e+01_EB, 1.549944e+01_EB, 1.573448e+00_EB, 4.766076e-01_EB, &
   4.349295e+01_EB, &
   4.836954e-01_EB, 1.451915e+00_EB, 6.485194e+01_EB, 3.798161e+01_EB, 6.188424e-01_EB, 1.109668e+00_EB, &
   4.352421e+01_EB/),(/7,8/))

gammad4_c2h4(1:7,9:16) = RESHAPE((/ &  ! 3000-3175 cm-1
   3.345606e-01_EB, 5.278129e-01_EB, 6.469998e+01_EB, 1.288707e+01_EB, 8.242915e-01_EB, 1.035867e+00_EB, &
   4.352463e+01_EB, &
   1.563175e+00_EB, 2.625564e+00_EB, 6.454977e+01_EB, 6.093823e+01_EB, 2.035803e+00_EB, 1.210794e+00_EB, &
   4.352445e+01_EB, &
   8.952736e-01_EB, 1.070318e+00_EB, 1.932614e+00_EB, 1.478721e+01_EB, 1.380362e+00_EB, 4.608167e-01_EB, &
   4.352008e+01_EB, &
   1.218284e+00_EB, 7.511695e-01_EB, 9.571422e-01_EB, 3.171139e+00_EB, 1.952744e+00_EB, 1.419092e+00_EB, &
   4.351690e+01_EB, &
   1.188336e+00_EB, 5.437254e-01_EB, 7.645704e-01_EB, 1.220082e+00_EB, 1.364237e+00_EB, 2.953192e-01_EB, &
   3.944037e+01_EB, &
   1.069749e+00_EB, 7.136400e-01_EB, 5.480781e-01_EB, 7.511426e-01_EB, 8.072618e-01_EB, 2.442109e-01_EB, &
   1.336082e+01_EB, &
   5.123621e-01_EB, 4.456214e-01_EB, 3.621148e-01_EB, 6.473524e-01_EB, 5.364253e-01_EB, 1.595842e-01_EB, &
   2.523820e+01_EB, &
   2.185443e-01_EB, 1.605263e-01_EB, 1.238794e-01_EB, 2.209229e-01_EB, 1.353909e-01_EB, 5.702117e-02_EB, &
   3.092779e+01_EB/),(/7,8/))

gammad4_c2h4(1:7,17:24) = RESHAPE((/ &  ! 3200-3375 cm-1
   5.244654e-02_EB, 7.872155e-02_EB, 5.156657e-02_EB, 8.595139e-02_EB, 5.856846e-02_EB, 2.220866e-02_EB, &
   1.517203e-01_EB, &
   3.762586e-03_EB, 1.216532e-02_EB, 1.917467e-02_EB, 2.846263e-02_EB, 2.089697e-02_EB, 4.969412e-03_EB, &
   1.387417e-01_EB, &
   4.291422e-05_EB, 7.807504e-04_EB, 1.341258e-03_EB, 6.152627e-03_EB, 2.309479e-03_EB, 1.073258e-03_EB, &
   3.311750e-01_EB, &
   3.584912e-05_EB, 2.671120e-04_EB, 1.957471e-04_EB, 1.342098e-03_EB, 1.278685e-04_EB, 3.324302e-03_EB, &
   3.684340e-03_EB, &
   8.216113e-05_EB, 2.682532e-04_EB, 2.247056e-04_EB, 2.195120e-04_EB, 7.774207e-05_EB, 3.815559e-05_EB, &
   2.872761e-03_EB, &
   1.294855e-04_EB, 2.265159e-04_EB, 2.506069e-04_EB, 2.536762e-04_EB, 2.424375e-04_EB, 2.435020e-04_EB, &
   1.033795e-04_EB, &
   2.206754e-04_EB, 2.418204e-04_EB, 2.208317e-04_EB, 2.403401e-04_EB, 2.291933e-04_EB, 2.353414e-04_EB, &
   2.133339e-04_EB, &
   2.205017e-04_EB, 1.924693e-04_EB, 2.146102e-04_EB, 2.118201e-04_EB, 2.108670e-04_EB, 2.100997e-04_EB, &
   9.661135e-05_EB/),(/7,8/))

gammad4_c2h4(1:7,25:25) = RESHAPE((/ &  ! 3400-3400 cm-1
   1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, 1.709157e-04_EB, &
   8.125358e-05_EB/),(/7,1/))

!-------------------------heptane data-------------------


! there are 2 bands for heptane

! band #1: 1100 cm-1 - 1800 cm-1 
! band #2: 2550 cm-1 - 3275 cm-1 

ALLOCATE(sd_c7h16_temp(n_temp_c7h16)) 

! initialize bands wavenumber bounds for heptane ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_c7h16(n_band_c7h16,3)) 
ALLOCATE(be_c7h16(n_band_c7h16)) 

om_bnd_c7h16 = RESHAPE((/ &
   1100._EB, 2550._EB, &
   1800._EB, 3275._EB, &
   25._EB, 25._EB/),(/n_band_c7h16,3/)) 

sd_c7h16_temp = (/ &
   293._EB, 400._EB, 450._EB, 490._EB, 593._EB, 794._EB,&
   1000._EB/)

be_c7h16 = (/ &
   1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_c7h16(n_temp_c7h16,29)) 

! band #1: 1100 cm-1 - 1800 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.72543 % 

sd1_c7h16(1:7,1:8) = RESHAPE((/ &  ! 1100-1275 cm-1
   1.426744e-03_EB, 1.201770e-03_EB, 1.274320e-03_EB, 1.347928e-03_EB, 1.656416e-03_EB, 1.325255e-03_EB, &
   2.057861e-03_EB, &
   1.811130e-03_EB, 2.521922e-03_EB, 1.872758e-03_EB, 1.707256e-03_EB, 1.927722e-03_EB, 1.955475e-03_EB, &
   2.516760e-01_EB, &
   1.747928e-03_EB, 2.120166e-03_EB, 2.046218e-03_EB, 1.649154e-03_EB, 1.666396e-03_EB, 1.992630e-03_EB, &
   2.005874e+00_EB, &
   1.624614e-03_EB, 2.021985e-03_EB, 4.956340e-03_EB, 1.877145e-03_EB, 1.905600e-03_EB, 1.746241e-03_EB, &
   3.502611e-01_EB, &
   1.930396e-03_EB, 1.894347e-03_EB, 1.299287e-02_EB, 1.476436e-02_EB, 4.597890e-03_EB, 2.359195e-02_EB, &
   4.013873e+00_EB, &
   6.004547e-02_EB, 4.017701e-02_EB, 3.934126e-02_EB, 3.083190e-02_EB, 2.179721e-01_EB, 4.014522e-02_EB, &
   2.600484e+00_EB, &
   7.613398e-02_EB, 5.549015e-02_EB, 3.003541e-01_EB, 1.344358e-01_EB, 4.411584e-01_EB, 7.993576e-02_EB, &
   4.665685e+00_EB, &
   1.373972e-01_EB, 1.249248e-01_EB, 1.230316e-01_EB, 1.155939e-01_EB, 3.749324e-01_EB, 1.429788e-01_EB, &
   2.876741e+00_EB/),(/7,8/))

sd1_c7h16(1:7,9:16) = RESHAPE((/ &  ! 1300-1475 cm-1
   1.960350e-01_EB, 7.251248e-01_EB, 4.660543e-01_EB, 2.468984e-01_EB, 3.746538e-01_EB, 1.932187e-01_EB, &
   1.906827e+00_EB, &
   2.182462e-01_EB, 1.037076e+00_EB, 4.532828e-01_EB, 3.341447e-01_EB, 4.124866e-01_EB, 3.832806e-01_EB, &
   1.538802e+00_EB, &
   5.511687e-01_EB, 9.486430e-01_EB, 6.841490e-01_EB, 5.497692e-01_EB, 5.893967e-01_EB, 6.049837e-01_EB, &
   6.095137e+00_EB, &
   2.169296e+00_EB, 1.762238e+00_EB, 1.496288e+00_EB, 1.153363e+00_EB, 9.332138e-01_EB, 6.553914e-01_EB, &
   1.559757e+00_EB, &
   1.152586e+00_EB, 1.398012e+00_EB, 1.008290e+00_EB, 8.691367e-01_EB, 7.857883e-01_EB, 5.751618e-01_EB, &
   3.554067e+00_EB, &
   2.809542e-01_EB, 2.417925e+00_EB, 4.381657e-01_EB, 5.312573e-01_EB, 5.715073e-01_EB, 5.154692e-01_EB, &
   1.149283e+00_EB, &
   2.520765e+00_EB, 2.328593e+00_EB, 1.895157e+00_EB, 1.577062e+00_EB, 1.328041e+00_EB, 1.069522e+00_EB, &
   2.814896e+00_EB, &
   4.096478e+00_EB, 3.000085e+00_EB, 2.158821e+00_EB, 1.809690e+00_EB, 1.347468e+00_EB, 8.991423e-01_EB, &
   2.342350e+00_EB/),(/7,8/))

sd1_c7h16(1:7,17:24) = RESHAPE((/ &  ! 1500-1675 cm-1
   4.388670e-01_EB, 4.751103e-01_EB, 5.711364e-01_EB, 4.125599e-01_EB, 4.895056e-01_EB, 3.320026e-01_EB, &
   2.177032e+00_EB, &
   1.319705e-01_EB, 1.769183e-01_EB, 1.924812e-01_EB, 1.563493e-01_EB, 2.730796e-01_EB, 1.797482e-01_EB, &
   1.084963e+00_EB, &
   1.057772e-01_EB, 1.102722e-01_EB, 1.166527e-01_EB, 7.750779e-02_EB, 1.371680e-01_EB, 7.530555e-02_EB, &
   4.673835e+00_EB, &
   1.595659e-03_EB, 8.558479e-01_EB, 3.816360e-02_EB, 1.698589e-01_EB, 3.558832e-02_EB, 1.126854e+00_EB, &
   6.274235e+00_EB, &
   1.754518e-03_EB, 3.490672e-01_EB, 1.460628e-01_EB, 3.092481e-01_EB, 1.777143e-03_EB, 5.457342e+00_EB, &
   4.901413e+00_EB, &
   1.574552e-03_EB, 2.022684e-03_EB, 1.771614e-03_EB, 8.573879e-02_EB, 1.949107e-03_EB, 8.452931e-01_EB, &
   6.243195e+00_EB, &
   1.827339e-03_EB, 2.129922e-03_EB, 2.108249e-03_EB, 4.442893e-03_EB, 1.724490e-03_EB, 3.038490e-03_EB, &
   3.195029e+00_EB, &
   1.956178e-03_EB, 2.075598e-03_EB, 1.738792e-03_EB, 1.976314e-03_EB, 1.668323e-03_EB, 1.791602e-03_EB, &
   4.975893e+00_EB/),(/7,8/))

sd1_c7h16(1:7,25:29) = RESHAPE((/ &  ! 1700-1800 cm-1
   1.635422e-03_EB, 2.360211e-03_EB, 1.966190e-03_EB, 1.838431e-03_EB, 1.806395e-03_EB, 2.025722e-03_EB, &
   3.595214e+00_EB, &
   1.685412e-03_EB, 1.845260e-03_EB, 1.904451e-03_EB, 1.947180e-03_EB, 1.880024e-03_EB, 1.570401e-03_EB, &
   4.352547e+00_EB, &
   1.693693e-03_EB, 2.153411e-03_EB, 2.016532e-03_EB, 1.811016e-03_EB, 1.826653e-03_EB, 1.799395e-03_EB, &
   1.496407e+00_EB, &
   1.300712e-03_EB, 1.598790e-03_EB, 1.518162e-03_EB, 1.529646e-03_EB, 1.291412e-03_EB, 1.245281e-03_EB, &
   1.384262e-01_EB, &
   9.565465e-04_EB, 9.572222e-04_EB, 9.564413e-04_EB, 9.564850e-04_EB, 9.561277e-04_EB, 9.585852e-04_EB, &
   1.922450e-03_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_c7h16(n_temp_c7h16,29)) 

! band #1: 1100 cm-1 - 1800 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.72543 % 

! print fine structure array gamma_d 

gammad1_c7h16(1:7,1:8) = RESHAPE((/ &  ! 1100-1275 cm-1
   1.079377e-03_EB, 8.321757e-04_EB, 1.106526e-03_EB, 8.589890e-04_EB, 1.218045e-03_EB, 9.003631e-04_EB, &
   1.652303e-03_EB, &
   1.219534e-03_EB, 1.678161e-03_EB, 1.027749e-03_EB, 1.087240e-03_EB, 1.300295e-03_EB, 1.381191e-03_EB, &
   2.766162e-01_EB, &
   1.117750e-03_EB, 1.217712e-03_EB, 1.293043e-03_EB, 1.061960e-03_EB, 1.158654e-03_EB, 1.485883e-03_EB, &
   6.838761e-04_EB, &
   1.138939e-03_EB, 1.252798e-03_EB, 2.807650e-03_EB, 1.167864e-03_EB, 1.225067e-03_EB, 1.314096e-03_EB, &
   2.282961e-02_EB, &
   1.420120e-03_EB, 1.020475e-03_EB, 6.385377e-04_EB, 2.240175e-03_EB, 1.390330e-03_EB, 3.503565e-01_EB, &
   2.532075e-04_EB, &
   1.303079e-02_EB, 1.684651e-02_EB, 1.656972e-02_EB, 6.785015e-03_EB, 1.293076e-04_EB, 1.692435e-02_EB, &
   9.530208e-04_EB, &
   1.370385e-02_EB, 1.586469e-02_EB, 1.596226e-04_EB, 6.447976e-04_EB, 2.699714e-04_EB, 4.064881e-01_EB, &
   7.587172e-04_EB, &
   1.038621e-02_EB, 6.544712e-03_EB, 6.378851e-03_EB, 9.284895e-03_EB, 4.414444e-04_EB, 4.010176e-01_EB, &
   1.276912e-03_EB/),(/7,8/))

gammad1_c7h16(1:7,9:16) = RESHAPE((/ &  ! 1300-1475 cm-1
   9.663644e-02_EB, 3.809049e-04_EB, 7.060483e-04_EB, 2.346751e-03_EB, 1.601612e-03_EB, 2.284259e-01_EB, &
   3.174484e-03_EB, &
   2.648317e-02_EB, 4.915970e-04_EB, 1.642486e-03_EB, 3.803913e-03_EB, 3.050190e-03_EB, 1.055777e-02_EB, &
   4.320127e-03_EB, &
   3.165256e-02_EB, 3.335451e-03_EB, 6.233122e-03_EB, 1.162381e-02_EB, 1.333659e-02_EB, 1.163372e-02_EB, &
   8.774580e-04_EB, &
   7.316298e-02_EB, 3.067612e-02_EB, 2.893574e-02_EB, 4.067615e-02_EB, 5.248547e-02_EB, 4.649271e-01_EB, &
   1.044691e-02_EB, &
   4.533824e-02_EB, 1.035218e-02_EB, 1.958798e-02_EB, 2.401292e-02_EB, 2.064596e-02_EB, 3.122001e-02_EB, &
   3.284470e-03_EB, &
   2.887489e-02_EB, 4.618211e-04_EB, 1.291555e-02_EB, 9.570390e-03_EB, 1.139081e-02_EB, 4.131392e-02_EB, &
   2.376013e-02_EB, &
   9.435371e-02_EB, 3.559624e-02_EB, 4.477145e-02_EB, 5.321919e-02_EB, 5.410966e-02_EB, 4.965735e-02_EB, &
   3.060972e-03_EB, &
   1.434765e-01_EB, 5.756528e-02_EB, 7.898010e-02_EB, 6.856662e-02_EB, 6.369824e-02_EB, 6.660056e-02_EB, &
   3.552779e-03_EB/),(/7,8/))

gammad1_c7h16(1:7,17:24) = RESHAPE((/ &  ! 1500-1675 cm-1
   6.164439e-02_EB, 4.097206e-02_EB, 2.017864e-02_EB, 5.078487e-02_EB, 2.150953e-02_EB, 6.143261e-01_EB, &
   2.140967e-03_EB, &
   4.999099e-02_EB, 6.175675e-01_EB, 5.926015e-01_EB, 5.735879e-01_EB, 3.127071e-03_EB, 8.449651e-02_EB, &
   4.359548e-03_EB, &
   2.308160e-01_EB, 2.876922e+00_EB, 2.678573e+00_EB, 3.088436e-01_EB, 1.022662e-03_EB, 4.348389e-01_EB, &
   5.138371e-04_EB, &
   1.136554e-03_EB, 2.312863e-05_EB, 1.191809e+00_EB, 5.400494e-04_EB, 8.452069e-04_EB, 3.414922e-05_EB, &
   2.732241e-04_EB, &
   1.236428e-03_EB, 7.575429e-06_EB, 1.804662e-05_EB, 1.168575e-04_EB, 1.132712e-03_EB, 1.348325e-05_EB, &
   2.605356e-04_EB, &
   1.005482e-03_EB, 1.099106e-03_EB, 1.146793e-03_EB, 1.793408e-04_EB, 1.319446e-03_EB, 2.932722e-05_EB, &
   2.217645e-04_EB, &
   1.351400e-03_EB, 1.302077e-03_EB, 1.237548e-03_EB, 2.270156e-03_EB, 1.159496e-03_EB, 2.121847e-03_EB, &
   2.405564e-04_EB, &
   1.421491e-03_EB, 1.320464e-03_EB, 1.069063e-03_EB, 1.283897e-03_EB, 9.514195e-04_EB, 1.174814e-03_EB, &
   1.274543e-04_EB/),(/7,8/))

gammad1_c7h16(1:7,25:29) = RESHAPE((/ &  ! 1700-1800 cm-1
   1.029018e-03_EB, 1.522013e-03_EB, 1.190609e-03_EB, 1.137355e-03_EB, 1.227873e-03_EB, 1.521396e-03_EB, &
   4.984208e-05_EB, &
   1.062872e-03_EB, 1.240560e-03_EB, 1.115838e-03_EB, 1.071608e-03_EB, 1.073014e-03_EB, 1.091421e-03_EB, &
   1.043998e-04_EB, &
   1.264467e-03_EB, 1.452995e-03_EB, 1.276660e-03_EB, 1.062100e-03_EB, 1.296189e-03_EB, 1.314713e-03_EB, &
   1.161390e-04_EB, &
   9.490158e-04_EB, 1.088570e-03_EB, 1.017566e-03_EB, 1.036075e-03_EB, 9.453163e-04_EB, 9.054199e-04_EB, &
   1.849978e-04_EB, &
   7.632406e-04_EB, 7.638486e-04_EB, 7.631088e-04_EB, 7.630990e-04_EB, 7.628205e-04_EB, 7.652755e-04_EB, &
   1.536272e-03_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_c7h16(n_temp_c7h16,30)) 

! band #2: 2550 cm-1 - 3275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 5.5204 % 

sd2_c7h16(1:7,1:8) = RESHAPE((/ &  ! 2550-2725 cm-1
   1.230745e-03_EB, 1.059249e-02_EB, 5.258556e-03_EB, 8.983219e-03_EB, 1.071517e-02_EB, 1.270153e-03_EB, &
   2.238178e-03_EB, &
   1.746456e-02_EB, 4.672224e-02_EB, 3.759388e-02_EB, 4.616045e-02_EB, 5.001133e-02_EB, 1.722111e-02_EB, &
   2.491507e-02_EB, &
   1.081581e-01_EB, 1.027192e-01_EB, 8.800407e-02_EB, 8.251442e-02_EB, 6.268006e-02_EB, 4.814096e-02_EB, &
   7.543431e-02_EB, &
   1.324171e-01_EB, 1.186396e-01_EB, 8.522107e-02_EB, 8.615581e-02_EB, 7.161614e-02_EB, 4.907251e-02_EB, &
   4.189169e-02_EB, &
   1.360592e-01_EB, 1.511535e-01_EB, 9.963712e-02_EB, 9.989906e-02_EB, 7.841100e-02_EB, 6.989163e-02_EB, &
   1.135539e-01_EB, &
   1.798230e-01_EB, 2.192535e-01_EB, 1.208341e-01_EB, 1.243637e-01_EB, 2.183718e-01_EB, 7.565057e-02_EB, &
   8.186008e-02_EB, &
   2.869951e-01_EB, 2.102309e-01_EB, 1.285038e-01_EB, 1.353365e-01_EB, 1.198975e-01_EB, 1.097888e-01_EB, &
   8.920775e-02_EB, &
   4.131559e-01_EB, 3.718172e-01_EB, 1.977670e-01_EB, 1.999998e-01_EB, 2.661506e-01_EB, 1.630822e-01_EB, &
   1.430789e-01_EB/),(/7,8/))

sd2_c7h16(1:7,9:16) = RESHAPE((/ &  ! 2750-2925 cm-1
   4.814798e-01_EB, 4.371742e-01_EB, 2.946385e-01_EB, 2.539099e-01_EB, 3.004489e-01_EB, 2.639523e-01_EB, &
   2.106653e-01_EB, &
   4.711680e-01_EB, 4.231342e-01_EB, 2.816821e-01_EB, 2.662044e-01_EB, 2.714629e-01_EB, 4.380265e-01_EB, &
   2.595817e-01_EB, &
   6.539776e-01_EB, 7.016029e-01_EB, 5.178696e-01_EB, 4.535956e-01_EB, 5.002499e-01_EB, 7.197349e-01_EB, &
   3.199467e-01_EB, &
   1.241588e+00_EB, 1.255753e+00_EB, 9.891572e-01_EB, 8.923575e-01_EB, 8.902290e-01_EB, 1.114933e+00_EB, &
   4.631429e-01_EB, &
   5.884384e+00_EB, 4.563318e+00_EB, 3.951471e+00_EB, 3.339596e+00_EB, 2.714647e+00_EB, 2.641774e+00_EB, &
   7.520423e-01_EB, &
   1.531155e+01_EB, 1.096024e+01_EB, 9.162135e+00_EB, 7.803908e+00_EB, 5.562181e+00_EB, 4.589694e+00_EB, &
   1.056916e+00_EB, &
   1.317309e+01_EB, 1.037752e+01_EB, 9.226521e+00_EB, 8.166911e+00_EB, 6.406625e+00_EB, 6.017800e+00_EB, &
   1.580795e+00_EB, &
   2.423199e+01_EB, 1.841550e+01_EB, 1.601599e+01_EB, 1.391139e+01_EB, 1.034058e+01_EB, 9.135221e+00_EB, &
   2.166977e+00_EB/),(/7,8/))

sd2_c7h16(1:7,17:24) = RESHAPE((/ &  ! 2950-3125 cm-1
   2.666381e+01_EB, 2.072769e+01_EB, 1.799865e+01_EB, 1.588750e+01_EB, 1.152743e+01_EB, 9.524580e+00_EB, &
   2.012903e+00_EB, &
   3.139322e+01_EB, 2.055468e+01_EB, 1.643220e+01_EB, 1.409534e+01_EB, 8.965984e+00_EB, 6.302236e+00_EB, &
   1.518462e+00_EB, &
   2.594752e+00_EB, 2.921881e+00_EB, 2.424761e+00_EB, 2.367343e+00_EB, 1.941979e+00_EB, 1.810548e+00_EB, &
   1.163951e+00_EB, &
   9.124058e-01_EB, 1.242167e+00_EB, 8.462169e-01_EB, 7.637253e-01_EB, 7.930042e-01_EB, 7.315673e-01_EB, &
   9.222189e-01_EB, &
   5.694335e-01_EB, 1.040431e+00_EB, 4.731323e-01_EB, 4.570560e-01_EB, 5.207545e-01_EB, 4.518627e-01_EB, &
   7.258027e-01_EB, &
   4.615311e-01_EB, 7.262374e-01_EB, 2.975859e-01_EB, 3.257008e-01_EB, 3.595477e-01_EB, 3.278019e-01_EB, &
   5.488059e-01_EB, &
   4.537172e-01_EB, 5.187525e-01_EB, 2.166416e-01_EB, 2.523836e-01_EB, 3.336731e-01_EB, 2.793265e-01_EB, &
   4.844761e-01_EB, &
   1.173983e-01_EB, 4.043441e-01_EB, 1.764942e-01_EB, 2.034332e-01_EB, 2.260464e-01_EB, 2.126733e-01_EB, &
   5.112998e-01_EB/),(/7,8/))

sd2_c7h16(1:7,25:30) = RESHAPE((/ &  ! 3150-3275 cm-1
   1.745347e-01_EB, 3.869948e-01_EB, 1.606330e-01_EB, 1.691991e-01_EB, 1.773258e-01_EB, 1.545312e-01_EB, &
   3.754934e-01_EB, &
   2.125656e-01_EB, 3.845411e-01_EB, 1.519506e-01_EB, 1.522530e-01_EB, 1.848037e-01_EB, 1.013866e-01_EB, &
   2.632396e-01_EB, &
   2.543277e-01_EB, 2.891570e-01_EB, 9.975616e-02_EB, 1.088793e-01_EB, 1.322502e-01_EB, 5.343256e-02_EB, &
   1.951975e-01_EB, &
   6.412029e-02_EB, 2.180048e-01_EB, 9.631661e-02_EB, 5.927599e-02_EB, 1.626811e-02_EB, 2.207700e-02_EB, &
   1.199629e-01_EB, &
   4.940691e-03_EB, 1.520018e-02_EB, 8.071023e-03_EB, 1.120400e-02_EB, 5.605272e-03_EB, 1.455284e-03_EB, &
   4.956089e-02_EB, &
   9.565465e-04_EB, 9.572222e-04_EB, 9.564413e-04_EB, 9.564850e-04_EB, 9.561277e-04_EB, 9.585852e-04_EB, &
   1.922450e-03_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_c7h16(n_temp_c7h16,30)) 

! band #2: 2550 cm-1 - 3275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 5.5204 % 

! print fine structure array gamma_d 

gammad2_c7h16(1:7,1:8) = RESHAPE((/ &  ! 2550-2725 cm-1
   9.146079e-04_EB, 1.220358e-02_EB, 3.937261e-03_EB, 7.187786e-03_EB, 5.227505e-03_EB, 9.355109e-04_EB, &
   1.819838e-03_EB, &
   9.874420e-03_EB, 2.084439e-01_EB, 2.173266e-01_EB, 8.182000e-03_EB, 2.973184e-02_EB, 9.482167e-03_EB, &
   5.745211e-01_EB, &
   9.786100e-03_EB, 1.858210e-01_EB, 1.533293e-02_EB, 5.349239e-02_EB, 2.414250e-02_EB, 4.205505e-04_EB, &
   3.283538e-02_EB, &
   1.016087e-02_EB, 5.218556e-02_EB, 6.697325e-02_EB, 2.379606e-01_EB, 8.937135e-03_EB, 9.990945e-03_EB, &
   5.966306e-01_EB, &
   1.783704e-02_EB, 2.299629e-02_EB, 4.295102e-02_EB, 6.377060e-02_EB, 3.169911e-02_EB, 5.529323e-03_EB, &
   7.456241e-01_EB, &
   9.757832e-03_EB, 5.468131e-03_EB, 4.867751e-02_EB, 1.172578e-01_EB, 4.899078e-04_EB, 2.599541e-02_EB, &
   1.519840e+00_EB, &
   1.701122e-03_EB, 1.906085e-02_EB, 8.916232e-02_EB, 3.655990e-01_EB, 7.351637e-03_EB, 3.626665e-02_EB, &
   1.919410e+00_EB, &
   3.999808e-03_EB, 7.531740e-03_EB, 1.072796e-01_EB, 1.003155e+00_EB, 1.614637e-03_EB, 7.721701e-02_EB, &
   1.224303e+00_EB/),(/7,8/))

gammad2_c7h16(1:7,9:16) = RESHAPE((/ &  ! 2750-2925 cm-1
   4.880113e-03_EB, 8.565604e-03_EB, 1.367064e-02_EB, 3.997591e-01_EB, 4.442547e-03_EB, 1.915330e-02_EB, &
   5.051128e+00_EB, &
   1.825295e-03_EB, 7.512748e-03_EB, 1.729303e-02_EB, 1.104825e+00_EB, 1.586989e-02_EB, 1.138387e-02_EB, &
   1.724797e+00_EB, &
   8.938939e-03_EB, 1.035703e-02_EB, 1.615494e-02_EB, 9.389559e-02_EB, 1.525842e-02_EB, 1.850825e-02_EB, &
   3.698894e+00_EB, &
   2.478387e-02_EB, 2.301877e-02_EB, 4.297098e-02_EB, 1.503277e-01_EB, 6.014925e-02_EB, 8.742386e-02_EB, &
   2.622943e+00_EB, &
   1.601922e-01_EB, 1.114843e-01_EB, 1.112220e-01_EB, 1.804939e-01_EB, 2.339949e-01_EB, 1.584145e-01_EB, &
   9.606857e+00_EB, &
   5.154427e-01_EB, 3.320188e-01_EB, 3.025580e-01_EB, 3.557579e-01_EB, 6.640807e-01_EB, 3.241466e-01_EB, &
   2.268534e+01_EB, &
   4.187811e-01_EB, 3.414132e-01_EB, 3.225035e-01_EB, 3.833688e-01_EB, 7.058557e-01_EB, 4.035161e-01_EB, &
   4.330011e+01_EB, &
   8.053331e-01_EB, 5.054746e-01_EB, 4.771400e-01_EB, 5.477855e-01_EB, 1.088287e+00_EB, 5.418580e-01_EB, &
   4.326237e+01_EB/),(/7,8/))

gammad2_c7h16(1:7,17:24) = RESHAPE((/ &  ! 2950-3125 cm-1
   9.236854e-01_EB, 6.215496e-01_EB, 5.912471e-01_EB, 6.204509e-01_EB, 1.340069e+00_EB, 6.140174e-01_EB, &
   4.330342e+01_EB, &
   5.662789e-01_EB, 5.083725e-01_EB, 4.767821e-01_EB, 4.500428e-01_EB, 1.099144e+00_EB, 4.936755e-01_EB, &
   2.497692e+01_EB, &
   5.143360e-02_EB, 5.233965e-02_EB, 8.660531e-02_EB, 8.976548e-02_EB, 1.810249e-01_EB, 2.371358e-01_EB, &
   4.306803e+01_EB, &
   1.542965e-02_EB, 1.268218e-02_EB, 2.772033e-02_EB, 9.250133e-02_EB, 2.649493e-02_EB, 3.172831e-01_EB, &
   8.950148e+00_EB, &
   8.562522e-03_EB, 4.977221e-03_EB, 3.610411e-02_EB, 1.439318e+00_EB, 1.049615e-02_EB, 6.108043e-01_EB, &
   4.534910e+00_EB, &
   1.597126e-03_EB, 3.635349e-03_EB, 5.555783e-02_EB, 3.967871e+00_EB, 6.887030e-03_EB, 5.613768e-01_EB, &
   3.059256e+00_EB, &
   6.799236e-04_EB, 3.384843e-03_EB, 1.010555e-01_EB, 2.884777e+00_EB, 2.380478e-03_EB, 1.681744e-01_EB, &
   2.528805e+00_EB, &
   8.299446e-03_EB, 2.940292e-03_EB, 9.141446e-02_EB, 3.099637e+00_EB, 1.826256e-03_EB, 9.931857e-02_EB, &
   3.001536e+00_EB/),(/7,8/))

gammad2_c7h16(1:7,25:30) = RESHAPE((/ &  ! 3150-3275 cm-1
   1.022601e-02_EB, 2.439564e-03_EB, 2.129847e-02_EB, 1.571134e+00_EB, 1.678148e-03_EB, 9.784188e-02_EB, &
   7.913187e+00_EB, &
   2.072404e-02_EB, 2.738722e-03_EB, 2.883252e-02_EB, 5.459085e-01_EB, 6.513406e-04_EB, 1.967974e-01_EB, &
   1.236417e+00_EB, &
   3.185761e-03_EB, 2.226331e-03_EB, 3.340747e-02_EB, 2.954516e+00_EB, 3.127038e-04_EB, 3.897705e-02_EB, &
   5.973433e-01_EB, &
   3.201661e-03_EB, 5.002671e-04_EB, 5.027477e-04_EB, 5.059749e-01_EB, 2.581852e-03_EB, 5.653921e-03_EB, &
   4.902323e-01_EB, &
   3.091588e-03_EB, 1.287212e-03_EB, 4.025758e-03_EB, 1.683050e-02_EB, 2.725481e-03_EB, 9.810405e-04_EB, &
   1.476595e-02_EB, &
   7.632406e-04_EB, 7.638486e-04_EB, 7.631088e-04_EB, 7.630990e-04_EB, 7.628205e-04_EB, 7.652755e-04_EB, &
   1.536272e-03_EB/),(/7,6/))

!-------------------------methanol data-------------------


! there are 4 bands for methanol

! band #1: 825 cm-1 - 1150 cm-1 
! band #2: 1125 cm-1 - 1700 cm-1 
! band #3: 2600 cm-1 - 3225 cm-1 
! band #4: 3525 cm-1 - 3850 cm-1 

ALLOCATE(sd_ch3oh_temp(n_temp_ch3oh)) 

! initialize bands wavenumber bounds for methanol ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_ch3oh(n_band_ch3oh,3)) 
ALLOCATE(be_ch3oh(n_band_ch3oh)) 

om_bnd_ch3oh = RESHAPE((/ &
   825._EB, 1125._EB, 2600._EB, 3525._EB, &
   1150._EB, 1700._EB, 3225._EB, 3850._EB, &
   5._EB, 25._EB, 25._EB, 25._EB/),(/n_band_ch3oh,3/)) 

sd_ch3oh_temp = (/ &
   293._EB, 396._EB, 443._EB, 483._EB, 570._EB, 804._EB,&
   1000._EB/)

be_ch3oh = (/ &
   1.000_EB, 1.000_EB, 1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_ch3oh(n_temp_ch3oh,58)) 

! band #1: 825 cm-1 - 1150 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 3.9209 % 

sd1_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 825-860 cm-1
   7.075310e-04_EB, 8.159230e-04_EB, 6.305365e-04_EB, 6.771661e-04_EB, 6.382997e-04_EB, 6.733725e-04_EB, &
   5.317197e-04_EB, &
   1.027754e-03_EB, 7.803944e-04_EB, 9.690581e-04_EB, 8.327490e-04_EB, 1.014832e-03_EB, 1.009424e-03_EB, &
   4.548089e-02_EB, &
   9.923660e-04_EB, 7.762778e-04_EB, 8.630137e-04_EB, 1.137131e-03_EB, 8.918616e-04_EB, 1.053434e-03_EB, &
   2.195153e+00_EB, &
   9.759840e-04_EB, 1.083038e-03_EB, 9.878635e-04_EB, 9.020455e-04_EB, 9.799016e-04_EB, 6.911637e-04_EB, &
   6.514857e-01_EB, &
   1.006816e-03_EB, 9.248700e-04_EB, 8.933603e-04_EB, 9.409542e-04_EB, 7.796524e-04_EB, 1.033193e-03_EB, &
   3.104940e-01_EB, &
   1.262552e-03_EB, 9.296738e-04_EB, 1.134121e-03_EB, 1.041003e-03_EB, 8.506538e-04_EB, 1.025871e-03_EB, &
   4.257321e-01_EB, &
   9.314553e-04_EB, 9.719139e-04_EB, 1.050566e-03_EB, 1.107681e-03_EB, 7.432552e-04_EB, 7.341790e-03_EB, &
   2.976279e-01_EB, &
   9.166184e-04_EB, 9.949250e-04_EB, 1.142691e-03_EB, 1.003991e-03_EB, 9.570837e-04_EB, 2.474105e-02_EB, &
   3.966025e-02_EB/),(/7,8/))

sd1_ch3oh(1:7,9:16) = RESHAPE((/ &  ! 865-900 cm-1
   8.873778e-04_EB, 1.014874e-03_EB, 8.954091e-04_EB, 8.658071e-04_EB, 8.266735e-04_EB, 3.960677e-02_EB, &
   4.795769e-01_EB, &
   9.318468e-04_EB, 1.081375e-03_EB, 8.779975e-04_EB, 6.417478e-04_EB, 8.380587e-04_EB, 5.064673e-02_EB, &
   2.069845e+00_EB, &
   9.172948e-04_EB, 6.246642e-04_EB, 8.418469e-04_EB, 1.002610e-03_EB, 8.748320e-04_EB, 6.309985e-02_EB, &
   2.535762e+00_EB, &
   5.635661e-04_EB, 9.022735e-04_EB, 7.536768e-04_EB, 1.384075e-02_EB, 4.169450e-03_EB, 6.860592e-02_EB, &
   4.255594e+00_EB, &
   5.238187e-04_EB, 9.020580e-04_EB, 8.975650e-04_EB, 2.989369e-02_EB, 2.964629e-03_EB, 5.639540e-01_EB, &
   4.081000e+00_EB, &
   4.272194e-03_EB, 1.285085e-03_EB, 1.101707e-03_EB, 3.363383e-02_EB, 2.223112e-03_EB, 7.068202e-01_EB, &
   5.035905e+00_EB, &
   7.391996e-03_EB, 1.149033e-03_EB, 4.043676e-03_EB, 3.024886e-02_EB, 1.237039e-02_EB, 2.917510e-01_EB, &
   5.127708e+00_EB, &
   1.769445e-02_EB, 1.162378e-03_EB, 1.101607e-02_EB, 3.265720e-02_EB, 2.700163e-02_EB, 2.303736e-01_EB, &
   7.211748e+00_EB/),(/7,8/))

sd1_ch3oh(1:7,17:24) = RESHAPE((/ &  ! 905-940 cm-1
   1.041542e-02_EB, 1.114301e-02_EB, 1.311655e-02_EB, 4.064275e-02_EB, 7.266914e-02_EB, 2.423119e-01_EB, &
   4.957147e+00_EB, &
   2.050418e-02_EB, 1.592394e-02_EB, 3.294495e-02_EB, 5.520246e-02_EB, 7.496393e-02_EB, 2.368645e-01_EB, &
   1.130981e+01_EB, &
   2.109122e-02_EB, 1.222826e-02_EB, 4.727193e-02_EB, 8.561265e-02_EB, 1.062168e-01_EB, 1.049801e+00_EB, &
   7.996949e+00_EB, &
   3.772088e-02_EB, 8.854601e-03_EB, 6.435677e-02_EB, 7.648012e-02_EB, 1.545208e-01_EB, 2.887097e-01_EB, &
   3.363412e+00_EB, &
   5.889149e-02_EB, 2.614746e-02_EB, 8.847350e-02_EB, 1.405710e-01_EB, 1.772912e-01_EB, 3.718949e-01_EB, &
   9.837331e+00_EB, &
   6.654570e-02_EB, 5.749401e-02_EB, 9.802181e-02_EB, 1.688979e-01_EB, 2.243950e-01_EB, 1.220047e+00_EB, &
   1.655983e+01_EB, &
   7.945784e-02_EB, 9.272308e-02_EB, 1.367127e-01_EB, 1.949308e-01_EB, 2.725708e-01_EB, 6.381513e-01_EB, &
   1.228651e+01_EB, &
   9.776093e-02_EB, 1.303068e-01_EB, 1.815588e-01_EB, 2.578179e-01_EB, 3.469963e-01_EB, 5.560945e-01_EB, &
   1.183943e+01_EB/),(/7,8/))

sd1_ch3oh(1:7,25:32) = RESHAPE((/ &  ! 945-980 cm-1
   1.006010e-01_EB, 1.730636e-01_EB, 2.270729e-01_EB, 3.315982e-01_EB, 4.573468e-01_EB, 6.708389e-01_EB, &
   1.024291e+01_EB, &
   1.194413e-01_EB, 2.369771e-01_EB, 3.106016e-01_EB, 4.274250e-01_EB, 5.708632e-01_EB, 7.929006e-01_EB, &
   1.330083e+01_EB, &
   1.736727e-01_EB, 3.344439e-01_EB, 4.261052e-01_EB, 5.527598e-01_EB, 7.561373e-01_EB, 9.600217e-01_EB, &
   1.114027e+01_EB, &
   2.600298e-01_EB, 4.700460e-01_EB, 5.820182e-01_EB, 7.302241e-01_EB, 9.234681e-01_EB, 1.070772e+00_EB, &
   1.188692e+01_EB, &
   3.671640e-01_EB, 6.491391e-01_EB, 7.852073e-01_EB, 9.389176e-01_EB, 1.126578e+00_EB, 1.172724e+00_EB, &
   6.310257e+00_EB, &
   5.738182e-01_EB, 9.168206e-01_EB, 1.052392e+00_EB, 1.223815e+00_EB, 1.363750e+00_EB, 1.356380e+00_EB, &
   1.276687e+01_EB, &
   8.847952e-01_EB, 1.262063e+00_EB, 1.370575e+00_EB, 1.545198e+00_EB, 1.612648e+00_EB, 1.449724e+00_EB, &
   1.203918e+01_EB, &
   1.380983e+00_EB, 1.728933e+00_EB, 1.776913e+00_EB, 1.952279e+00_EB, 1.903225e+00_EB, 1.560147e+00_EB, &
   1.409457e+01_EB/),(/7,8/))

sd1_ch3oh(1:7,33:40) = RESHAPE((/ &  ! 985-1020 cm-1
   1.989792e+00_EB, 2.216523e+00_EB, 2.173252e+00_EB, 2.302944e+00_EB, 2.154186e+00_EB, 1.615229e+00_EB, &
   1.614576e+01_EB, &
   2.901685e+00_EB, 2.807597e+00_EB, 2.638654e+00_EB, 2.691349e+00_EB, 2.424882e+00_EB, 1.732381e+00_EB, &
   1.516364e+01_EB, &
   3.935938e+00_EB, 3.363705e+00_EB, 3.049578e+00_EB, 3.000715e+00_EB, 2.568481e+00_EB, 1.745273e+00_EB, &
   2.252866e+01_EB, &
   5.148871e+00_EB, 3.896474e+00_EB, 3.364794e+00_EB, 3.239344e+00_EB, 2.682511e+00_EB, 1.748386e+00_EB, &
   2.390649e+01_EB, &
   6.170452e+00_EB, 4.190528e+00_EB, 3.476275e+00_EB, 3.307799e+00_EB, 2.673466e+00_EB, 1.690532e+00_EB, &
   2.170835e+01_EB, &
   6.905606e+00_EB, 4.214646e+00_EB, 3.441984e+00_EB, 3.180517e+00_EB, 2.517520e+00_EB, 1.596762e+00_EB, &
   1.932637e+01_EB, &
   6.940706e+00_EB, 3.956671e+00_EB, 3.247771e+00_EB, 2.994245e+00_EB, 2.445520e+00_EB, 1.673228e+00_EB, &
   1.306371e+01_EB, &
   5.958695e+00_EB, 3.454990e+00_EB, 2.948441e+00_EB, 2.807961e+00_EB, 2.435395e+00_EB, 1.873623e+00_EB, &
   1.081489e+01_EB/),(/7,8/))

sd1_ch3oh(1:7,41:48) = RESHAPE((/ &  ! 1025-1060 cm-1
   4.201331e+00_EB, 2.701090e+00_EB, 2.459137e+00_EB, 2.429934e+00_EB, 2.270440e+00_EB, 1.776380e+00_EB, &
   1.211452e+01_EB, &
   6.021324e+00_EB, 4.632103e+00_EB, 4.201019e+00_EB, 4.002494e+00_EB, 3.407819e+00_EB, 2.266345e+00_EB, &
   1.495652e+01_EB, &
   9.188156e+00_EB, 5.163112e+00_EB, 4.909468e+00_EB, 4.147644e+00_EB, 3.320954e+00_EB, 1.941572e+00_EB, &
   1.668271e+01_EB, &
   3.550622e+00_EB, 2.256569e+00_EB, 2.077046e+00_EB, 1.997234e+00_EB, 1.831325e+00_EB, 1.549349e+00_EB, &
   1.016637e+01_EB, &
   6.406905e+00_EB, 3.553412e+00_EB, 2.981786e+00_EB, 2.818896e+00_EB, 2.393274e+00_EB, 1.845333e+00_EB, &
   1.267998e+01_EB, &
   9.119940e+00_EB, 5.130999e+00_EB, 4.171838e+00_EB, 3.898993e+00_EB, 3.183156e+00_EB, 2.173550e+00_EB, &
   1.068448e+01_EB, &
   1.001487e+01_EB, 6.103507e+00_EB, 4.932645e+00_EB, 4.654735e+00_EB, 3.723719e+00_EB, 2.430534e+00_EB, &
   1.889930e+01_EB, &
   8.945286e+00_EB, 6.111634e+00_EB, 5.090176e+00_EB, 4.850566e+00_EB, 3.927607e+00_EB, 2.511476e+00_EB, &
   2.721090e+01_EB/),(/7,8/))

sd1_ch3oh(1:7,49:56) = RESHAPE((/ &  ! 1065-1100 cm-1
   6.548008e+00_EB, 5.234485e+00_EB, 4.593759e+00_EB, 4.456866e+00_EB, 3.759323e+00_EB, 2.389545e+00_EB, &
   1.930330e+01_EB, &
   4.016006e+00_EB, 3.880275e+00_EB, 3.590672e+00_EB, 3.603798e+00_EB, 3.180636e+00_EB, 2.124012e+00_EB, &
   7.098842e+00_EB, &
   2.039641e+00_EB, 2.438251e+00_EB, 2.437403e+00_EB, 2.528162e+00_EB, 2.390431e+00_EB, 1.727051e+00_EB, &
   1.564504e+01_EB, &
   8.890902e-01_EB, 1.306925e+00_EB, 1.439104e+00_EB, 1.573364e+00_EB, 1.582749e+00_EB, 1.325349e+00_EB, &
   1.309975e+01_EB, &
   4.289045e-01_EB, 6.181922e-01_EB, 7.402675e-01_EB, 8.358344e-01_EB, 9.479778e-01_EB, 8.718487e-01_EB, &
   8.578315e+00_EB, &
   3.110189e-01_EB, 3.181081e-01_EB, 3.757707e-01_EB, 4.355864e-01_EB, 4.933576e-01_EB, 9.265380e-01_EB, &
   5.124942e+00_EB, &
   3.078577e-01_EB, 2.309721e-01_EB, 2.504640e-01_EB, 2.713709e-01_EB, 2.640710e-01_EB, 6.768124e-01_EB, &
   7.096896e+00_EB, &
   2.629835e-01_EB, 1.842053e-01_EB, 1.854752e-01_EB, 2.041121e-01_EB, 2.026328e-01_EB, 4.463207e-01_EB, &
   1.030077e+01_EB/),(/7,8/))

sd1_ch3oh(1:7,57:58) = RESHAPE((/ &  ! 1125-1150 cm-1
   1.066160e-01_EB, 8.237972e-02_EB, 7.650034e-02_EB, 7.750927e-02_EB, 8.650473e-02_EB, 5.002688e-01_EB, &
   2.961022e+00_EB, &
   4.699306e-04_EB, 4.696984e-04_EB, 4.697130e-04_EB, 4.702606e-04_EB, 4.699096e-04_EB, 4.682572e-04_EB, &
   4.691749e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_ch3oh(n_temp_ch3oh,58)) 

! band #1: 825 cm-1 - 1150 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 3.9209 % 

! print fine structure array gamma_d 

gammad1_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 825-860 cm-1
   6.320478e-04_EB, 6.861763e-04_EB, 4.976022e-04_EB, 5.212687e-04_EB, 4.951307e-04_EB, 4.722485e-04_EB, &
   4.126358e-04_EB, &
   5.857091e-04_EB, 4.986862e-04_EB, 7.181338e-04_EB, 6.000370e-04_EB, 8.426048e-04_EB, 6.854139e-04_EB, &
   9.797634e-05_EB, &
   6.687840e-04_EB, 5.631883e-04_EB, 6.493401e-04_EB, 8.974743e-04_EB, 8.364450e-04_EB, 9.128100e-04_EB, &
   2.228170e-05_EB, &
   7.015973e-04_EB, 8.509882e-04_EB, 6.749386e-04_EB, 7.104894e-04_EB, 9.244817e-04_EB, 4.993607e-04_EB, &
   3.618598e-05_EB, &
   6.879541e-04_EB, 6.793889e-04_EB, 6.064846e-04_EB, 7.494890e-04_EB, 6.021521e-04_EB, 8.041830e-04_EB, &
   4.227068e-05_EB, &
   1.115437e-03_EB, 6.123108e-04_EB, 9.950372e-04_EB, 7.526989e-04_EB, 5.762193e-04_EB, 6.564205e-04_EB, &
   4.006551e-05_EB, &
   7.510632e-04_EB, 6.557116e-04_EB, 7.646716e-04_EB, 8.360125e-04_EB, 5.785516e-04_EB, 7.166127e-02_EB, &
   3.307520e-05_EB, &
   7.077134e-04_EB, 6.379859e-04_EB, 8.630039e-04_EB, 8.226252e-04_EB, 6.815535e-04_EB, 2.157192e-01_EB, &
   1.186508e-04_EB/),(/7,8/))

gammad1_ch3oh(1:7,9:16) = RESHAPE((/ &  ! 865-900 cm-1
   6.878550e-04_EB, 7.626337e-04_EB, 7.231979e-04_EB, 6.559557e-04_EB, 6.188317e-04_EB, 1.528514e-01_EB, &
   3.379362e-05_EB, &
   8.329430e-04_EB, 7.234964e-04_EB, 7.086063e-04_EB, 4.799785e-04_EB, 5.996684e-04_EB, 1.058608e-01_EB, &
   2.576308e-05_EB, &
   7.201976e-04_EB, 4.879678e-04_EB, 7.514895e-04_EB, 9.296249e-04_EB, 7.627646e-04_EB, 1.357754e-02_EB, &
   4.112898e-05_EB, &
   4.468231e-04_EB, 7.152109e-04_EB, 5.892145e-04_EB, 6.268495e-03_EB, 6.890662e-04_EB, 8.344058e-03_EB, &
   3.782398e-05_EB, &
   4.115617e-04_EB, 6.608084e-04_EB, 6.730864e-04_EB, 1.653705e-02_EB, 7.764525e-04_EB, 2.258557e-04_EB, &
   7.122008e-05_EB, &
   2.243172e-03_EB, 1.292195e-03_EB, 1.270788e-03_EB, 1.990852e-02_EB, 5.353682e-02_EB, 3.379605e-04_EB, &
   1.094851e-04_EB, &
   2.795631e-03_EB, 8.919100e-04_EB, 5.452501e-03_EB, 5.460032e-02_EB, 4.954888e-01_EB, 2.416852e-03_EB, &
   8.087006e-05_EB, &
   3.021916e-03_EB, 8.815026e-04_EB, 2.682624e-03_EB, 1.363428e-01_EB, 7.040860e-01_EB, 9.532477e-03_EB, &
   4.720845e-05_EB/),(/7,8/))

gammad1_ch3oh(1:7,17:24) = RESHAPE((/ &  ! 905-940 cm-1
   1.302330e-02_EB, 5.102234e-03_EB, 2.422490e-01_EB, 1.768788e-01_EB, 2.586447e-01_EB, 3.526812e-01_EB, &
   4.129581e-04_EB, &
   2.544888e-02_EB, 1.060450e-02_EB, 1.381583e-02_EB, 4.518627e-01_EB, 5.206079e-01_EB, 8.391062e-03_EB, &
   2.632629e-04_EB, &
   2.533298e-01_EB, 6.433475e-03_EB, 2.725983e-02_EB, 2.678852e-01_EB, 2.275531e-01_EB, 1.261497e-03_EB, &
   2.704299e-04_EB, &
   3.122114e-02_EB, 8.346136e-03_EB, 2.487987e-01_EB, 9.374202e-02_EB, 6.604262e-02_EB, 2.488402e-02_EB, &
   2.562359e-04_EB, &
   3.870602e-02_EB, 2.158011e-01_EB, 1.910186e-01_EB, 2.065914e-01_EB, 8.461225e-01_EB, 2.529266e-01_EB, &
   3.603466e-04_EB, &
   7.517862e-02_EB, 2.869129e-01_EB, 5.031265e-01_EB, 3.707385e-01_EB, 3.378211e-01_EB, 7.697419e-03_EB, &
   3.347542e-04_EB, &
   1.159008e-01_EB, 5.649415e-01_EB, 2.185566e-01_EB, 2.313598e+00_EB, 1.481866e+00_EB, 3.664479e-02_EB, &
   3.612492e-04_EB, &
   2.101148e-01_EB, 1.362737e+00_EB, 1.681817e-01_EB, 2.403762e+00_EB, 5.052362e+00_EB, 4.829425e+01_EB, &
   7.256697e-04_EB/),(/7,8/))

gammad1_ch3oh(1:7,25:32) = RESHAPE((/ &  ! 945-980 cm-1
   2.154223e-01_EB, 1.063134e+00_EB, 2.745270e+00_EB, 6.048114e+00_EB, 1.191546e+01_EB, 3.613074e+01_EB, &
   9.758013e-04_EB, &
   3.415582e-01_EB, 5.670554e+00_EB, 1.626968e+00_EB, 8.377516e+00_EB, 2.646977e+01_EB, 4.827207e+01_EB, &
   6.374449e-04_EB, &
   7.233287e-01_EB, 1.170538e+01_EB, 2.699023e+00_EB, 1.811784e+01_EB, 4.833090e+01_EB, 4.829379e+01_EB, &
   6.697140e-04_EB, &
   4.410713e+00_EB, 2.163944e+01_EB, 1.204215e+01_EB, 6.230210e+01_EB, 5.735624e+01_EB, 4.829428e+01_EB, &
   6.621854e-04_EB, &
   6.410135e+00_EB, 5.657130e+01_EB, 3.102339e+01_EB, 6.230874e+01_EB, 5.735677e+01_EB, 4.829401e+01_EB, &
   1.509172e-03_EB, &
   3.756246e+01_EB, 6.880041e+01_EB, 4.129381e+01_EB, 6.230370e+01_EB, 5.735689e+01_EB, 4.829405e+01_EB, &
   8.487264e-04_EB, &
   6.111221e+01_EB, 6.881385e+01_EB, 6.471263e+01_EB, 6.230875e+01_EB, 5.735697e+01_EB, 4.829401e+01_EB, &
   1.318189e-03_EB, &
   7.999989e+01_EB, 6.881377e+01_EB, 6.495080e+01_EB, 6.230820e+01_EB, 5.735697e+01_EB, 4.828988e+01_EB, &
   1.354435e-03_EB/),(/7,8/))

gammad1_ch3oh(1:7,33:40) = RESHAPE((/ &  ! 985-1020 cm-1
   7.999998e+01_EB, 6.881384e+01_EB, 6.494110e+01_EB, 6.230890e+01_EB, 5.735615e+01_EB, 4.829425e+01_EB, &
   1.479087e-03_EB, &
   7.999920e+01_EB, 6.881390e+01_EB, 6.506107e+01_EB, 6.230863e+01_EB, 5.735699e+01_EB, 4.829154e+01_EB, &
   1.456376e-03_EB, &
   7.999917e+01_EB, 6.881386e+01_EB, 6.499928e+01_EB, 6.230892e+01_EB, 5.735688e+01_EB, 1.272502e+01_EB, &
   7.723150e-04_EB, &
   7.999980e+01_EB, 6.881301e+01_EB, 6.505184e+01_EB, 6.230892e+01_EB, 5.735698e+01_EB, 4.828738e+01_EB, &
   6.758482e-04_EB, &
   8.000000e+01_EB, 6.881338e+01_EB, 6.504268e+01_EB, 6.230892e+01_EB, 5.735686e+01_EB, 4.829415e+01_EB, &
   7.919076e-04_EB, &
   8.000000e+01_EB, 6.881339e+01_EB, 6.506098e+01_EB, 6.230892e+01_EB, 5.735694e+01_EB, 4.829379e+01_EB, &
   1.044482e-03_EB, &
   8.000000e+01_EB, 6.881323e+01_EB, 6.503974e+01_EB, 6.230732e+01_EB, 5.735696e+01_EB, 4.829424e+01_EB, &
   1.883098e-03_EB, &
   8.000000e+01_EB, 6.881386e+01_EB, 6.498651e+01_EB, 6.230722e+01_EB, 5.735624e+01_EB, 2.615034e+00_EB, &
   1.731620e-03_EB/),(/7,8/))

gammad1_ch3oh(1:7,41:48) = RESHAPE((/ &  ! 1025-1060 cm-1
   7.999918e+01_EB, 6.881390e+01_EB, 6.501143e+01_EB, 6.230770e+01_EB, 5.735573e+01_EB, 4.828762e+01_EB, &
   1.844304e-03_EB, &
   7.999986e+01_EB, 6.881374e+01_EB, 2.509078e+00_EB, 6.230892e+01_EB, 5.735699e+01_EB, 4.829418e+01_EB, &
   2.140274e-03_EB, &
   5.272622e-01_EB, 6.881390e+01_EB, 7.822208e-01_EB, 6.230892e+01_EB, 5.735699e+01_EB, 4.829426e+01_EB, &
   1.608061e-03_EB, &
   7.999988e+01_EB, 6.881380e+01_EB, 6.416330e+01_EB, 6.230736e+01_EB, 5.735695e+01_EB, 4.829395e+01_EB, &
   1.715694e-03_EB, &
   8.000000e+01_EB, 6.881387e+01_EB, 6.501663e+01_EB, 6.230833e+01_EB, 5.735699e+01_EB, 4.829374e+01_EB, &
   1.856662e-03_EB, &
   8.000000e+01_EB, 6.881386e+01_EB, 6.506116e+01_EB, 6.230710e+01_EB, 5.735699e+01_EB, 4.829358e+01_EB, &
   2.738523e-03_EB, &
   8.000000e+01_EB, 6.881391e+01_EB, 6.506116e+01_EB, 6.230892e+01_EB, 5.735699e+01_EB, 4.829418e+01_EB, &
   1.589317e-03_EB, &
   7.999991e+01_EB, 6.881391e+01_EB, 6.506118e+01_EB, 6.230892e+01_EB, 5.735618e+01_EB, 4.829397e+01_EB, &
   1.507561e-03_EB/),(/7,8/))

gammad1_ch3oh(1:7,49:56) = RESHAPE((/ &  ! 1065-1100 cm-1
   8.000000e+01_EB, 6.881387e+01_EB, 6.506116e+01_EB, 6.230827e+01_EB, 5.735699e+01_EB, 4.829427e+01_EB, &
   2.253348e-03_EB, &
   7.999938e+01_EB, 6.881318e+01_EB, 6.505314e+01_EB, 6.230887e+01_EB, 5.735699e+01_EB, 4.829173e+01_EB, &
   3.424844e-03_EB, &
   7.999997e+01_EB, 6.881353e+01_EB, 6.505405e+01_EB, 6.230890e+01_EB, 5.735697e+01_EB, 4.828966e+01_EB, &
   1.125303e-03_EB, &
   6.451623e+01_EB, 6.881383e+01_EB, 6.492757e+01_EB, 6.230886e+01_EB, 5.735691e+01_EB, 4.823782e+01_EB, &
   9.265505e-04_EB, &
   2.835561e+00_EB, 4.814647e+01_EB, 1.930069e+01_EB, 4.670295e+01_EB, 5.733661e+01_EB, 4.981992e-01_EB, &
   1.101043e-03_EB, &
   3.405886e-01_EB, 8.624550e+00_EB, 8.265512e-01_EB, 4.097800e+00_EB, 1.789592e+01_EB, 1.840689e-02_EB, &
   1.204450e-03_EB, &
   2.136135e+00_EB, 6.560155e+00_EB, 2.276106e-01_EB, 2.826746e-01_EB, 4.031936e+00_EB, 9.675596e-03_EB, &
   8.693678e-04_EB, &
   3.865526e-01_EB, 2.663625e+00_EB, 2.801253e-01_EB, 1.479540e-01_EB, 6.549643e-01_EB, 4.007539e-03_EB, &
   2.400175e-04_EB/),(/7,8/))

gammad1_ch3oh(1:7,57:58) = RESHAPE((/ &  ! 1125-1150 cm-1
   1.031691e-01_EB, 1.063157e+00_EB, 8.796742e-02_EB, 3.121383e-01_EB, 3.923037e-01_EB, 4.063955e-04_EB, &
   1.096427e-04_EB, &
   3.748583e-04_EB, 3.747168e-04_EB, 3.744971e-04_EB, 3.751403e-04_EB, 3.748851e-04_EB, 3.734005e-04_EB, &
   3.740585e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_ch3oh(n_temp_ch3oh,24)) 

! band #2: 1125 cm-1 - 1700 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 1.5713 % 

sd2_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 1125-1300 cm-1
   8.927831e-02_EB, 6.438189e-02_EB, 6.336484e-02_EB, 7.782879e-02_EB, 7.726760e-02_EB, 3.430637e-01_EB, &
   2.111028e+00_EB, &
   1.428237e-01_EB, 1.080083e-01_EB, 1.208600e-01_EB, 1.421780e-01_EB, 1.436580e-01_EB, 1.355584e+00_EB, &
   5.627932e+00_EB, &
   1.431595e-01_EB, 1.414695e-01_EB, 1.558110e-01_EB, 1.760891e-01_EB, 2.027903e-01_EB, 5.295318e-01_EB, &
   8.429621e+00_EB, &
   2.293295e-01_EB, 2.455796e-01_EB, 2.542758e-01_EB, 2.648440e-01_EB, 2.899987e-01_EB, 4.124789e-01_EB, &
   7.705341e+00_EB, &
   2.983198e-01_EB, 2.968924e-01_EB, 2.931528e-01_EB, 3.015322e-01_EB, 3.289150e-01_EB, 4.332046e-01_EB, &
   1.274789e+01_EB, &
   3.760820e-01_EB, 3.648052e-01_EB, 3.591136e-01_EB, 3.639314e-01_EB, 3.850232e-01_EB, 5.091509e-01_EB, &
   9.900445e+00_EB, &
   4.286786e-01_EB, 4.103485e-01_EB, 4.009699e-01_EB, 4.152912e-01_EB, 4.397207e-01_EB, 4.994214e-01_EB, &
   1.080146e+01_EB, &
   6.526914e-01_EB, 4.766714e-01_EB, 4.276416e-01_EB, 4.115988e-01_EB, 4.189427e-01_EB, 4.846849e-01_EB, &
   8.381233e+00_EB/),(/7,8/))

sd2_ch3oh(1:7,9:16) = RESHAPE((/ &  ! 1325-1500 cm-1
   8.926014e-01_EB, 5.725096e-01_EB, 4.936350e-01_EB, 4.461230e-01_EB, 4.220260e-01_EB, 4.240504e-01_EB, &
   6.735177e+00_EB, &
   9.802829e-01_EB, 6.469535e-01_EB, 5.500127e-01_EB, 4.937137e-01_EB, 4.581570e-01_EB, 4.189329e-01_EB, &
   5.750931e+00_EB, &
   7.911468e-01_EB, 5.793679e-01_EB, 5.044808e-01_EB, 4.592149e-01_EB, 4.299574e-01_EB, 3.354935e-01_EB, &
   5.762184e+00_EB, &
   6.069133e-01_EB, 4.730774e-01_EB, 4.244640e-01_EB, 3.791906e-01_EB, 3.723725e-01_EB, 3.634556e-01_EB, &
   5.212994e+00_EB, &
   5.069178e-01_EB, 3.597810e-01_EB, 3.086202e-01_EB, 2.782571e-01_EB, 2.798365e-01_EB, 3.419783e-01_EB, &
   5.269252e+00_EB, &
   4.541349e-01_EB, 3.336580e-01_EB, 2.970106e-01_EB, 2.657370e-01_EB, 2.766112e-01_EB, 3.361569e-01_EB, &
   5.737162e+00_EB, &
   4.615691e-01_EB, 3.258696e-01_EB, 2.757947e-01_EB, 2.454101e-01_EB, 2.523896e-01_EB, 2.077158e-01_EB, &
   3.866394e+00_EB, &
   3.328678e-01_EB, 2.037433e-01_EB, 1.740683e-01_EB, 1.606349e-01_EB, 1.789734e-01_EB, 3.064642e-01_EB, &
   2.819085e+00_EB/),(/7,8/))

sd2_ch3oh(1:7,17:24) = RESHAPE((/ &  ! 1525-1700 cm-1
   2.653588e-01_EB, 1.712954e-01_EB, 1.458324e-01_EB, 1.336402e-01_EB, 1.497627e-01_EB, 6.752522e-01_EB, &
   1.211105e+00_EB, &
   1.719421e-01_EB, 9.991170e-02_EB, 8.707576e-02_EB, 9.851840e-02_EB, 1.045626e-01_EB, 1.257519e-01_EB, &
   1.729420e-01_EB, &
   9.561175e-02_EB, 1.214100e-01_EB, 1.009499e-01_EB, 8.557814e-02_EB, 1.061461e-01_EB, 5.279653e-02_EB, &
   1.371839e+00_EB, &
   5.679817e-02_EB, 1.064996e-01_EB, 8.829344e-02_EB, 5.535303e-02_EB, 8.540729e-02_EB, 6.869857e-02_EB, &
   2.390774e+00_EB, &
   2.299715e-02_EB, 5.668754e-02_EB, 4.611447e-02_EB, 2.942043e-02_EB, 5.283754e-02_EB, 1.368461e-02_EB, &
   8.704598e-02_EB, &
   9.994737e-04_EB, 6.379671e-03_EB, 5.806672e-03_EB, 3.377512e-03_EB, 2.456039e-02_EB, 4.999823e-04_EB, &
   1.518390e-02_EB, &
   6.504118e-04_EB, 6.258215e-04_EB, 6.777009e-04_EB, 5.831397e-04_EB, 1.258943e-02_EB, 6.613745e-04_EB, &
   3.926253e-03_EB, &
   4.699306e-04_EB, 4.696984e-04_EB, 4.697130e-04_EB, 4.702606e-04_EB, 4.699096e-04_EB, 4.682572e-04_EB, &
   4.691749e-04_EB/),(/7,8/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_ch3oh(n_temp_ch3oh,24)) 

! band #2: 1125 cm-1 - 1700 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 1.5713 % 

! print fine structure array gamma_d 

gammad2_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 1125-1300 cm-1
   2.036126e-02_EB, 2.341570e-01_EB, 2.115157e-01_EB, 3.745991e-02_EB, 2.969993e-01_EB, 6.950924e-04_EB, &
   1.090294e-04_EB, &
   4.318091e-02_EB, 3.213276e+00_EB, 1.818315e-01_EB, 6.478423e-02_EB, 7.314703e-01_EB, 8.584914e-04_EB, &
   1.822575e-04_EB, &
   1.945666e-02_EB, 1.155771e+00_EB, 8.508445e-01_EB, 1.670850e-01_EB, 2.280568e+00_EB, 4.205261e-03_EB, &
   3.102895e-04_EB, &
   7.511583e-01_EB, 6.453461e+00_EB, 3.392847e+00_EB, 1.506522e+00_EB, 5.522265e+00_EB, 1.714255e-02_EB, &
   3.985647e-04_EB, &
   1.205702e+00_EB, 9.217844e+00_EB, 2.021655e+00_EB, 2.387126e+00_EB, 4.807045e+00_EB, 2.526665e-02_EB, &
   3.559097e-04_EB, &
   2.746003e+00_EB, 9.503086e+00_EB, 2.328040e+00_EB, 4.485330e+00_EB, 9.413318e+00_EB, 1.788133e-02_EB, &
   4.729684e-04_EB, &
   7.531540e+00_EB, 1.737151e+01_EB, 3.326148e+00_EB, 8.198222e+00_EB, 1.532668e+01_EB, 4.522148e-02_EB, &
   4.549597e-04_EB, &
   2.576303e+01_EB, 2.203775e+01_EB, 4.437035e+00_EB, 9.898502e+00_EB, 1.872281e+01_EB, 3.906498e-02_EB, &
   4.821874e-04_EB/),(/7,8/))

gammad2_ch3oh(1:7,9:16) = RESHAPE((/ &  ! 1325-1500 cm-1
   6.250399e+01_EB, 3.711813e+01_EB, 4.149099e+00_EB, 1.502754e+01_EB, 1.389076e+01_EB, 4.031216e-02_EB, &
   5.674173e-04_EB, &
   7.999962e+01_EB, 5.252798e+01_EB, 1.230226e+01_EB, 1.351493e+01_EB, 1.576973e+01_EB, 3.751308e-02_EB, &
   4.823587e-04_EB, &
   6.092331e+01_EB, 3.805490e+01_EB, 1.271376e+01_EB, 1.115539e+01_EB, 1.762731e+01_EB, 7.120929e-02_EB, &
   2.295392e-04_EB, &
   3.106485e+01_EB, 2.102259e+01_EB, 4.586676e+00_EB, 6.831440e+00_EB, 9.421920e+00_EB, 1.574921e-02_EB, &
   2.188143e-04_EB, &
   1.678441e+01_EB, 9.704951e+00_EB, 2.972980e+00_EB, 2.004021e+00_EB, 5.086223e+00_EB, 6.963213e-03_EB, &
   1.950359e-04_EB, &
   9.295170e+00_EB, 9.919157e+00_EB, 2.792442e+00_EB, 3.932078e+00_EB, 5.121134e+00_EB, 9.304463e-03_EB, &
   2.472872e-04_EB, &
   8.180857e+00_EB, 9.293708e+00_EB, 2.066943e+00_EB, 4.076114e+00_EB, 4.294751e+00_EB, 2.727838e-02_EB, &
   4.265301e-04_EB, &
   8.289335e-02_EB, 2.291116e+00_EB, 8.466376e-01_EB, 1.801527e-01_EB, 1.505300e+00_EB, 1.548374e-03_EB, &
   2.852002e-04_EB/),(/7,8/))

gammad2_ch3oh(1:7,17:24) = RESHAPE((/ &  ! 1525-1700 cm-1
   7.067783e-02_EB, 1.479749e+00_EB, 5.187742e-01_EB, 1.372559e-01_EB, 7.400371e-01_EB, 2.949015e-04_EB, &
   4.953459e-04_EB, &
   1.628627e-02_EB, 8.476648e-01_EB, 2.184428e-01_EB, 1.854959e-02_EB, 1.335169e-01_EB, 2.574359e-04_EB, &
   6.601786e-04_EB, &
   1.778534e-02_EB, 2.130514e-01_EB, 1.040976e-01_EB, 1.031289e-01_EB, 1.600152e+00_EB, 2.537685e-02_EB, &
   1.933608e-05_EB, &
   6.344549e-03_EB, 1.274538e-01_EB, 5.382005e-02_EB, 6.261896e-01_EB, 4.338573e-01_EB, 3.757499e-01_EB, &
   6.947226e-05_EB, &
   2.030029e-04_EB, 6.186341e-02_EB, 2.709645e-02_EB, 5.640946e-02_EB, 3.663756e-01_EB, 2.286239e-02_EB, &
   1.077139e-04_EB, &
   8.095102e-04_EB, 6.039758e-03_EB, 3.546569e-03_EB, 3.077717e-03_EB, 3.777927e-01_EB, 4.051839e-04_EB, &
   8.507547e-03_EB, &
   4.920711e-04_EB, 4.753092e-04_EB, 5.221505e-04_EB, 4.531404e-04_EB, 1.775124e-01_EB, 5.028014e-04_EB, &
   4.011702e-03_EB, &
   3.748583e-04_EB, 3.747168e-04_EB, 3.744971e-04_EB, 3.751403e-04_EB, 3.748851e-04_EB, 3.734005e-04_EB, &
   3.740585e-04_EB/),(/7,8/))

!---------------------------------------------------------------------------
ALLOCATE(sd3_ch3oh(n_temp_ch3oh,26)) 

! band #3: 2600 cm-1 - 3225 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.0044 % 

sd3_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 2600-2775 cm-1
   6.834504e-04_EB, 6.520912e-04_EB, 6.293446e-04_EB, 7.602714e-04_EB, 6.782192e-04_EB, 5.546971e-04_EB, &
   4.991362e-04_EB, &
   9.382971e-04_EB, 9.146731e-04_EB, 9.400615e-04_EB, 9.494733e-04_EB, 9.709058e-04_EB, 8.544557e-04_EB, &
   1.012754e-01_EB, &
   8.700561e-04_EB, 8.964927e-04_EB, 9.462079e-04_EB, 7.740041e-04_EB, 9.999124e-04_EB, 3.951283e-03_EB, &
   5.295951e-02_EB, &
   9.087410e-04_EB, 8.586469e-04_EB, 1.061772e-03_EB, 9.199452e-04_EB, 7.793155e-04_EB, 2.101256e-02_EB, &
   1.014346e-01_EB, &
   9.791713e-04_EB, 7.304428e-04_EB, 3.244007e-03_EB, 2.026606e-03_EB, 4.855695e-03_EB, 6.777374e-02_EB, &
   6.074970e-01_EB, &
   9.694186e-03_EB, 8.820685e-03_EB, 2.014998e-02_EB, 1.559350e-02_EB, 2.203169e-02_EB, 1.271102e-01_EB, &
   1.115356e+00_EB, &
   7.029245e-02_EB, 5.091559e-02_EB, 6.271267e-02_EB, 6.250553e-02_EB, 7.985290e-02_EB, 2.008186e-01_EB, &
   7.823174e-01_EB, &
   1.832465e-01_EB, 1.685849e-01_EB, 1.767557e-01_EB, 1.929714e-01_EB, 2.182618e-01_EB, 3.370225e-01_EB, &
   9.170478e-01_EB/),(/7,8/))

sd3_ch3oh(1:7,9:16) = RESHAPE((/ &  ! 2800-2975 cm-1
   6.431025e-01_EB, 5.822856e-01_EB, 5.579032e-01_EB, 5.709313e-01_EB, 5.423547e-01_EB, 5.399764e-01_EB, &
   1.804940e+00_EB, &
   1.616617e+00_EB, 1.106968e+00_EB, 9.574460e-01_EB, 9.001676e-01_EB, 7.635145e-01_EB, 6.273449e-01_EB, &
   2.500203e+00_EB, &
   2.054107e+00_EB, 1.363124e+00_EB, 1.172320e+00_EB, 1.104390e+00_EB, 9.352928e-01_EB, 7.059677e-01_EB, &
   1.705210e+00_EB, &
   2.373043e+00_EB, 1.746858e+00_EB, 1.529178e+00_EB, 1.473170e+00_EB, 1.254399e+00_EB, 9.322847e-01_EB, &
   2.050430e+00_EB, &
   2.371351e+00_EB, 1.769945e+00_EB, 1.583483e+00_EB, 1.552549e+00_EB, 1.370339e+00_EB, 1.037337e+00_EB, &
   2.789724e+00_EB, &
   3.792115e+00_EB, 2.516241e+00_EB, 2.137945e+00_EB, 2.022427e+00_EB, 1.666670e+00_EB, 1.107211e+00_EB, &
   3.100678e+00_EB, &
   4.157174e+00_EB, 2.644308e+00_EB, 2.220834e+00_EB, 2.075201e+00_EB, 1.691865e+00_EB, 1.086754e+00_EB, &
   2.915371e+00_EB, &
   4.342764e+00_EB, 2.822870e+00_EB, 2.375959e+00_EB, 2.221659e+00_EB, 1.802922e+00_EB, 1.138893e+00_EB, &
   2.903207e+00_EB/),(/7,8/))

sd3_ch3oh(1:7,17:24) = RESHAPE((/ &  ! 3000-3175 cm-1
   2.739997e+00_EB, 1.931251e+00_EB, 1.703000e+00_EB, 1.624103e+00_EB, 1.381073e+00_EB, 9.258508e-01_EB, &
   2.377676e+00_EB, &
   1.657562e+00_EB, 1.227931e+00_EB, 1.089977e+00_EB, 1.050288e+00_EB, 9.156831e-01_EB, 6.349324e-01_EB, &
   1.478204e+00_EB, &
   7.792319e-01_EB, 6.776318e-01_EB, 6.362424e-01_EB, 6.325745e-01_EB, 5.728339e-01_EB, 4.149136e-01_EB, &
   8.597464e-01_EB, &
   2.993662e-01_EB, 3.092377e-01_EB, 3.217485e-01_EB, 3.267839e-01_EB, 3.261281e-01_EB, 2.499769e-01_EB, &
   5.435448e-01_EB, &
   7.556116e-02_EB, 1.013165e-01_EB, 1.286472e-01_EB, 1.311502e-01_EB, 1.509531e-01_EB, 1.344009e-01_EB, &
   4.354789e-01_EB, &
   2.499115e-02_EB, 1.505850e-02_EB, 4.570610e-02_EB, 3.921605e-02_EB, 5.480372e-02_EB, 5.712297e-02_EB, &
   9.435909e-01_EB, &
   4.505936e-03_EB, 1.692846e-03_EB, 7.624509e-03_EB, 2.230769e-03_EB, 1.179331e-02_EB, 1.334128e-02_EB, &
   4.731345e-01_EB, &
   8.990606e-04_EB, 7.925739e-04_EB, 1.307211e-03_EB, 8.823890e-04_EB, 3.090415e-03_EB, 1.825819e-03_EB, &
   5.903177e-03_EB/),(/7,8/))

sd3_ch3oh(1:7,25:26) = RESHAPE((/ &  ! 3200-3225 cm-1
   6.930486e-04_EB, 6.167383e-04_EB, 6.581892e-04_EB, 7.470153e-04_EB, 6.327299e-04_EB, 5.451720e-04_EB, &
   6.008264e-04_EB, &
   4.699306e-04_EB, 4.696984e-04_EB, 4.697130e-04_EB, 4.702606e-04_EB, 4.699096e-04_EB, 4.682572e-04_EB, &
   4.691749e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad3_ch3oh(n_temp_ch3oh,26)) 

! band #3: 2600 cm-1 - 3225 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.0044 % 

! print fine structure array gamma_d 

gammad3_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 2600-2775 cm-1
   5.383708e-04_EB, 5.150011e-04_EB, 4.711666e-04_EB, 5.877341e-04_EB, 5.325580e-04_EB, 4.526624e-04_EB, &
   3.953361e-04_EB, &
   6.469173e-04_EB, 6.295204e-04_EB, 6.229213e-04_EB, 6.140709e-04_EB, 7.108644e-04_EB, 5.982461e-04_EB, &
   3.920520e-04_EB, &
   6.148513e-04_EB, 6.637986e-04_EB, 6.390155e-04_EB, 5.400884e-04_EB, 6.784448e-04_EB, 3.991532e-03_EB, &
   4.998359e-03_EB, &
   6.548999e-04_EB, 6.227472e-04_EB, 7.730036e-04_EB, 6.279630e-04_EB, 5.750964e-04_EB, 8.889816e-02_EB, &
   4.610782e-03_EB, &
   7.658461e-04_EB, 5.677730e-04_EB, 3.744908e-03_EB, 1.790178e-03_EB, 4.848659e-03_EB, 5.816371e-02_EB, &
   7.773551e-04_EB, &
   9.155583e-03_EB, 4.468866e-02_EB, 7.356051e-02_EB, 9.249021e-03_EB, 4.330986e-02_EB, 7.759327e-01_EB, &
   1.031984e-03_EB, &
   3.240618e-01_EB, 4.813829e-01_EB, 2.376949e-01_EB, 9.693334e-02_EB, 1.522188e-01_EB, 1.138999e+00_EB, &
   2.244386e-03_EB, &
   1.958993e+00_EB, 1.883288e+00_EB, 9.427313e-01_EB, 2.032231e+00_EB, 2.865213e+00_EB, 8.627391e+00_EB, &
   2.968586e-03_EB/),(/7,8/))

gammad3_ch3oh(1:7,9:16) = RESHAPE((/ &  ! 2800-2975 cm-1
   5.408634e+01_EB, 4.260803e+01_EB, 1.570059e+01_EB, 2.520800e+01_EB, 3.812191e+01_EB, 3.643507e+01_EB, &
   1.949224e-03_EB, &
   8.000000e+01_EB, 6.880853e+01_EB, 4.378312e+01_EB, 6.230866e+01_EB, 5.732026e+01_EB, 3.975357e+01_EB, &
   1.722003e-03_EB, &
   7.999962e+01_EB, 6.881387e+01_EB, 6.160988e+01_EB, 6.230878e+01_EB, 5.734000e+01_EB, 4.829424e+01_EB, &
   3.362841e-03_EB, &
   7.999999e+01_EB, 6.881386e+01_EB, 6.505766e+01_EB, 6.230889e+01_EB, 5.735690e+01_EB, 4.829419e+01_EB, &
   3.574523e-03_EB, &
   7.999997e+01_EB, 6.881383e+01_EB, 6.505113e+01_EB, 6.230881e+01_EB, 5.735699e+01_EB, 4.829428e+01_EB, &
   2.799128e-03_EB, &
   7.999848e+01_EB, 6.881323e+01_EB, 6.505412e+01_EB, 6.230885e+01_EB, 5.735695e+01_EB, 4.829421e+01_EB, &
   2.167350e-03_EB, &
   7.999924e+01_EB, 6.881280e+01_EB, 6.505343e+01_EB, 6.230882e+01_EB, 5.735689e+01_EB, 4.829429e+01_EB, &
   2.206560e-03_EB, &
   7.999955e+01_EB, 6.881390e+01_EB, 6.505410e+01_EB, 6.230862e+01_EB, 5.735674e+01_EB, 4.829421e+01_EB, &
   2.226025e-03_EB/),(/7,8/))

gammad3_ch3oh(1:7,17:24) = RESHAPE((/ &  ! 3000-3175 cm-1
   8.000000e+01_EB, 6.881380e+01_EB, 6.505082e+01_EB, 6.230881e+01_EB, 5.735619e+01_EB, 4.829420e+01_EB, &
   2.316390e-03_EB, &
   7.999999e+01_EB, 6.881373e+01_EB, 4.658312e+01_EB, 6.230870e+01_EB, 5.735696e+01_EB, 4.160676e+01_EB, &
   1.752164e-03_EB, &
   7.996301e+01_EB, 6.881042e+01_EB, 7.394547e+00_EB, 3.799419e+01_EB, 2.909775e+01_EB, 1.888125e+01_EB, &
   1.732017e-03_EB, &
   4.261947e+00_EB, 5.499787e+00_EB, 4.853050e-01_EB, 6.858640e+00_EB, 8.915265e+00_EB, 4.488725e+00_EB, &
   1.039419e-03_EB, &
   3.182069e-01_EB, 2.435677e-01_EB, 7.167320e-02_EB, 3.943245e-01_EB, 1.078816e+00_EB, 5.812817e-01_EB, &
   5.867531e-04_EB, &
   4.059798e-01_EB, 7.132746e-03_EB, 2.110640e-02_EB, 1.885784e-02_EB, 1.959525e+00_EB, 1.180753e-01_EB, &
   1.052399e-04_EB, &
   6.304409e-03_EB, 7.479753e-04_EB, 2.811545e-03_EB, 1.771085e-03_EB, 3.653146e-01_EB, 6.263784e-03_EB, &
   5.856904e-05_EB, &
   6.352137e-04_EB, 5.913538e-04_EB, 7.804488e-04_EB, 6.507595e-04_EB, 5.533275e-03_EB, 5.897126e-04_EB, &
   1.359107e-03_EB/),(/7,8/))

gammad3_ch3oh(1:7,25:26) = RESHAPE((/ &  ! 3200-3225 cm-1
   5.327466e-04_EB, 4.880598e-04_EB, 4.928840e-04_EB, 6.016314e-04_EB, 5.375037e-04_EB, 4.225641e-04_EB, &
   4.726310e-04_EB, &
   3.748583e-04_EB, 3.747168e-04_EB, 3.744971e-04_EB, 3.751403e-04_EB, 3.748851e-04_EB, 3.734005e-04_EB, &
   3.740585e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(sd4_ch3oh(n_temp_ch3oh,14)) 

! band #4: 3525 cm-1 - 3850 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.52389 % 

sd4_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 3525-3700 cm-1
   6.712614e-04_EB, 6.207569e-04_EB, 6.722895e-04_EB, 6.499733e-04_EB, 6.502738e-04_EB, 6.875865e-04_EB, &
   5.963773e-04_EB, &
   8.659465e-04_EB, 7.116916e-03_EB, 2.432786e-03_EB, 5.672617e-03_EB, 1.138304e-03_EB, 8.423409e-04_EB, &
   3.290015e-02_EB, &
   1.352414e-02_EB, 3.203665e-02_EB, 2.755663e-02_EB, 2.769632e-02_EB, 2.312063e-02_EB, 5.525230e-02_EB, &
   2.986065e-02_EB, &
   4.836204e-02_EB, 5.823870e-02_EB, 5.880972e-02_EB, 6.018123e-02_EB, 6.237734e-02_EB, 6.647353e-02_EB, &
   7.540357e-02_EB, &
   2.393220e-01_EB, 2.245383e-01_EB, 2.030768e-01_EB, 2.229738e-01_EB, 2.012059e-01_EB, 1.546924e-01_EB, &
   9.432293e-02_EB, &
   8.391763e-01_EB, 5.928915e-01_EB, 4.976034e-01_EB, 4.745560e-01_EB, 3.784750e-01_EB, 2.223088e-01_EB, &
   1.157200e-01_EB, &
   1.156927e+00_EB, 6.735305e-01_EB, 5.357276e-01_EB, 4.781313e-01_EB, 3.584284e-01_EB, 1.968791e-01_EB, &
   6.758499e-02_EB, &
   1.214378e+00_EB, 7.022996e-01_EB, 5.366277e-01_EB, 4.875065e-01_EB, 3.656436e-01_EB, 2.096662e-01_EB, &
   6.457647e-01_EB/),(/7,8/))

sd4_ch3oh(1:7,9:14) = RESHAPE((/ &  ! 3725-3850 cm-1
   7.394573e-01_EB, 5.605469e-01_EB, 4.645577e-01_EB, 4.379377e-01_EB, 3.556244e-01_EB, 2.353483e-01_EB, &
   1.093933e+00_EB, &
   1.627441e-01_EB, 1.282880e-01_EB, 1.200820e-01_EB, 1.269546e-01_EB, 1.078069e-01_EB, 1.197911e-01_EB, &
   5.594711e-02_EB, &
   5.446422e-02_EB, 6.975799e-02_EB, 5.570485e-02_EB, 2.768076e-02_EB, 2.935100e-02_EB, 4.703585e-02_EB, &
   5.054507e-01_EB, &
   5.503798e-03_EB, 3.001374e-02_EB, 1.548540e-02_EB, 1.507541e-03_EB, 3.052107e-03_EB, 1.675834e-02_EB, &
   3.138134e-02_EB, &
   6.972876e-04_EB, 3.117091e-03_EB, 2.511646e-03_EB, 6.126449e-04_EB, 6.508959e-04_EB, 6.472961e-04_EB, &
   6.127914e-04_EB, &
   4.699306e-04_EB, 4.696984e-04_EB, 4.697130e-04_EB, 4.702606e-04_EB, 4.699096e-04_EB, 4.682572e-04_EB, &
   4.691749e-04_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(gammad4_ch3oh(n_temp_ch3oh,14)) 

! band #4: 3525 cm-1 - 3850 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.52389 % 

! print fine structure array gamma_d 

gammad4_ch3oh(1:7,1:8) = RESHAPE((/ &  ! 3525-3700 cm-1
   5.122880e-04_EB, 4.870888e-04_EB, 5.202731e-04_EB, 4.858325e-04_EB, 4.894141e-04_EB, 5.234560e-04_EB, &
   4.739454e-04_EB, &
   6.149338e-04_EB, 1.383613e-01_EB, 1.894640e-03_EB, 1.652906e-01_EB, 1.023597e-03_EB, 6.016220e-04_EB, &
   4.058911e-03_EB, &
   3.510951e-01_EB, 4.231322e-02_EB, 2.094183e-02_EB, 1.756637e-01_EB, 8.897080e-02_EB, 1.613316e-04_EB, &
   2.603264e-03_EB, &
   4.348848e-01_EB, 1.997004e+00_EB, 9.906747e-02_EB, 6.264650e-02_EB, 3.368153e-02_EB, 5.344243e-03_EB, &
   1.024940e-02_EB, &
   6.050885e+00_EB, 5.037136e+00_EB, 1.610741e+00_EB, 9.196821e-01_EB, 3.334769e+00_EB, 1.286001e-01_EB, &
   1.852205e-02_EB, &
   7.999493e+01_EB, 4.970642e+01_EB, 1.283805e+01_EB, 1.192477e+01_EB, 8.469410e+00_EB, 7.793742e-01_EB, &
   2.887610e-03_EB, &
   7.999963e+01_EB, 6.880115e+01_EB, 1.309018e+01_EB, 1.179133e+01_EB, 1.455179e+01_EB, 1.751987e-01_EB, &
   5.542345e-03_EB, &
   7.999935e+01_EB, 6.881382e+01_EB, 1.685368e+01_EB, 1.646915e+01_EB, 1.087273e+01_EB, 1.258692e+00_EB, &
   2.305314e-04_EB/),(/7,8/))

gammad4_ch3oh(1:7,9:14) = RESHAPE((/ &  ! 3725-3850 cm-1
   6.925146e+01_EB, 3.835584e+01_EB, 7.883505e+00_EB, 1.488897e+01_EB, 1.369263e+01_EB, 1.664808e+00_EB, &
   3.295444e-04_EB, &
   5.002891e-01_EB, 7.416100e-01_EB, 2.838529e-01_EB, 3.020566e-01_EB, 1.523043e+00_EB, 1.133794e-03_EB, &
   4.486303e-03_EB, &
   6.823303e-03_EB, 4.678797e-01_EB, 9.198815e-02_EB, 1.770170e-01_EB, 1.483378e-01_EB, 1.869665e-01_EB, &
   2.436489e-05_EB, &
   4.419159e-03_EB, 3.220528e-01_EB, 2.310260e-01_EB, 1.151947e-03_EB, 2.324707e-03_EB, 6.947554e-03_EB, &
   1.285540e-04_EB, &
   5.281954e-04_EB, 2.447127e-03_EB, 2.223998e-03_EB, 4.730215e-04_EB, 5.179966e-04_EB, 4.804723e-04_EB, &
   4.227659e-04_EB, &
   3.748583e-04_EB, 3.747168e-04_EB, 3.744971e-04_EB, 3.751403e-04_EB, 3.748851e-04_EB, 3.734005e-04_EB, &
   3.740585e-04_EB/),(/7,6/))

!-------------------------mma data-------------------


! there are 6 bands for mma

! band #1: 750 cm-1 - 880 cm-1 
! band #2: 875 cm-1 - 1055 cm-1 
! band #3: 1050 cm-1 - 1275 cm-1 
! band #4: 1250 cm-1 - 1575 cm-1 
! band #5: 1550 cm-1 - 1975 cm-1 
! band #6: 2650 cm-1 - 3275 cm-1 

ALLOCATE(sd_c5h8o2_temp(n_temp_c5h8o2)) 

! initialize bands wavenumber bounds for mma ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_c5h8o2(n_band_c5h8o2,3)) 
ALLOCATE(be_c5h8o2(n_band_c5h8o2)) 

om_bnd_c5h8o2 = RESHAPE((/ &
   750._EB,  875._EB, 1050._EB, 1250._EB, 1550._EB, 2650._EB, &
   880._EB, 1055._EB, 1275._EB, 1575._EB, 1975._EB, 3275._EB, &
   5._EB, 5._EB, 5._EB, 25._EB, 25._EB, 25._EB/),(/n_band_c5h8o2,3/)) 

sd_c5h8o2_temp = (/ &
   297._EB, 396._EB, 441._EB, 483._EB, 597._EB, 803._EB,&
   1014._EB/)

be_c5h8o2 = (/ &
   1.000_EB, 1.000_EB, 1.000_EB, 1.000_EB, 1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_c5h8o2(n_temp_c5h8o2,27)) 

! band #1: 750 cm-1 - 880 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 8.0046 % 

sd1_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 750-785 cm-1
   6.713660e-04_EB, 5.825251e-04_EB, 7.916094e-04_EB, 6.543173e-04_EB, 6.421674e-04_EB, 6.087532e-04_EB, &
   7.107910e-04_EB, &
   9.848065e-04_EB, 9.607899e-04_EB, 7.994027e-04_EB, 1.090851e-03_EB, 1.087202e-03_EB, 7.734671e-04_EB, &
   9.284135e-04_EB, &
   9.110967e-04_EB, 7.183894e-04_EB, 8.010414e-04_EB, 1.063396e-03_EB, 1.000702e-03_EB, 6.657414e-04_EB, &
   6.127113e-04_EB, &
   8.267669e-04_EB, 9.315188e-04_EB, 8.290909e-04_EB, 7.969467e-04_EB, 8.913743e-04_EB, 6.855529e-04_EB, &
   7.046623e-04_EB, &
   1.042663e-03_EB, 1.043538e-03_EB, 8.872311e-04_EB, 8.694322e-04_EB, 9.329158e-04_EB, 9.772556e-04_EB, &
   7.458198e-04_EB, &
   1.006578e-03_EB, 7.209128e-04_EB, 1.019682e-03_EB, 7.712991e-04_EB, 3.187441e-02_EB, 5.248309e-02_EB, &
   8.138133e-04_EB, &
   7.585006e-04_EB, 8.285426e-04_EB, 7.850058e-04_EB, 7.514949e-04_EB, 2.435151e+00_EB, 2.382451e-01_EB, &
   7.995357e-04_EB, &
   5.346436e-03_EB, 1.695093e-03_EB, 1.005551e-02_EB, 1.807196e-02_EB, 4.893034e+00_EB, 1.072109e+00_EB, &
   8.263934e-04_EB/),(/7,8/))

sd1_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 790-825 cm-1
   7.889654e-02_EB, 6.932866e-02_EB, 6.943889e-02_EB, 1.657322e-01_EB, 4.820655e+00_EB, 1.595563e+00_EB, &
   6.604864e-04_EB, &
   2.859272e-01_EB, 2.060379e-01_EB, 1.784774e-01_EB, 2.655142e-01_EB, 5.069557e+00_EB, 1.587178e+00_EB, &
   8.491025e-04_EB, &
   6.960576e-01_EB, 3.989489e-01_EB, 3.440851e-01_EB, 4.569709e-01_EB, 4.118131e+00_EB, 7.958553e-01_EB, &
   7.688756e-04_EB, &
   1.250951e+00_EB, 6.520021e-01_EB, 6.089583e-01_EB, 6.456713e-01_EB, 2.261652e+00_EB, 7.908891e-01_EB, &
   8.388555e-04_EB, &
   2.208710e+00_EB, 1.110849e+00_EB, 1.035737e+00_EB, 1.049217e+00_EB, 1.084162e+00_EB, 7.448165e-01_EB, &
   1.121289e-03_EB, &
   3.745705e+00_EB, 1.437198e+00_EB, 1.241512e+00_EB, 1.206161e+00_EB, 8.198992e-01_EB, 7.316076e-01_EB, &
   7.825346e-04_EB, &
   3.469611e+00_EB, 1.274191e+00_EB, 1.123297e+00_EB, 1.063600e+00_EB, 7.020966e-01_EB, 8.989672e-01_EB, &
   7.836122e-04_EB, &
   3.387588e+00_EB, 1.219143e+00_EB, 1.110425e+00_EB, 9.818704e-01_EB, 6.253886e-01_EB, 1.450947e+00_EB, &
   8.923055e-04_EB/),(/7,8/))

sd1_c5h8o2(1:7,17:24) = RESHAPE((/ &  ! 830-865 cm-1
   2.650263e+00_EB, 9.577831e-01_EB, 8.753717e-01_EB, 8.295316e-01_EB, 6.934556e-01_EB, 1.431367e+00_EB, &
   8.664361e-04_EB, &
   1.730390e+00_EB, 6.490013e-01_EB, 6.546973e-01_EB, 5.810689e-01_EB, 2.534337e+00_EB, 7.747536e-01_EB, &
   7.994517e-04_EB, &
   1.002821e+00_EB, 3.608083e-01_EB, 3.892233e-01_EB, 3.747237e-01_EB, 6.259059e+00_EB, 1.624771e+00_EB, &
   7.025073e-04_EB, &
   4.114065e-01_EB, 1.202763e-01_EB, 1.380148e-01_EB, 2.829570e-01_EB, 7.146988e+00_EB, 1.774384e+00_EB, &
   9.127284e-04_EB, &
   1.274057e-01_EB, 9.521312e-04_EB, 7.626901e-04_EB, 5.309748e-02_EB, 6.530803e+00_EB, 3.950028e-01_EB, &
   8.260701e-04_EB, &
   3.101890e-02_EB, 9.960440e-04_EB, 1.450745e-03_EB, 8.889566e-04_EB, 4.823124e+00_EB, 5.607234e-02_EB, &
   1.289121e-03_EB, &
   1.815363e-02_EB, 6.716907e-04_EB, 1.134173e-03_EB, 7.672105e-04_EB, 2.305115e+00_EB, 1.215781e-02_EB, &
   9.963103e-03_EB, &
   1.671588e-02_EB, 9.534992e-04_EB, 7.842672e-04_EB, 1.050185e-03_EB, 8.279107e-04_EB, 9.473562e-04_EB, &
   1.520855e-02_EB/),(/7,8/))

sd1_c5h8o2(1:7,25:27) = RESHAPE((/ &  ! 870-880 cm-1
   3.388305e-02_EB, 8.968940e-04_EB, 9.267088e-04_EB, 1.004537e-03_EB, 9.822602e-04_EB, 8.968929e-04_EB, &
   3.122362e-02_EB, &
   2.663964e-02_EB, 7.688918e-04_EB, 6.407620e-04_EB, 5.722080e-04_EB, 9.499110e-04_EB, 6.115005e-04_EB, &
   3.351275e-02_EB, &
   4.683959e-04_EB, 4.682485e-04_EB, 4.694920e-04_EB, 4.693301e-04_EB, 4.691565e-04_EB, 4.675221e-04_EB, &
   4.659141e-04_EB/),(/7,3/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_c5h8o2(n_temp_c5h8o2,27)) 

! band #1: 750 cm-1 - 880 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 8.0046 % 

! print fine structure array gamma_d 

gammad1_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 750-785 cm-1
   4.881844e-04_EB, 4.346454e-04_EB, 5.958338e-04_EB, 4.571093e-04_EB, 4.822105e-04_EB, 4.530250e-04_EB, &
   5.750719e-04_EB, &
   6.727325e-04_EB, 6.267647e-04_EB, 5.477130e-04_EB, 7.657452e-04_EB, 1.077214e-03_EB, 5.499367e-04_EB, &
   6.460319e-04_EB, &
   5.864217e-04_EB, 5.127786e-04_EB, 5.824402e-04_EB, 8.990885e-04_EB, 8.469550e-04_EB, 5.075038e-04_EB, &
   5.043115e-04_EB, &
   6.108399e-04_EB, 7.428392e-04_EB, 5.525420e-04_EB, 5.867387e-04_EB, 7.082636e-04_EB, 5.889027e-04_EB, &
   4.893654e-04_EB, &
   7.434780e-04_EB, 6.504896e-04_EB, 6.589451e-04_EB, 6.802512e-04_EB, 7.417022e-04_EB, 7.262805e-04_EB, &
   5.511207e-04_EB, &
   6.356817e-04_EB, 5.339199e-04_EB, 7.935459e-04_EB, 4.886757e-04_EB, 2.583559e-03_EB, 3.792064e-03_EB, &
   6.024047e-04_EB, &
   5.106586e-04_EB, 5.822074e-04_EB, 5.456352e-04_EB, 5.195168e-04_EB, 1.249121e-04_EB, 2.441600e-03_EB, &
   6.419995e-04_EB, &
   8.911013e-04_EB, 1.366487e-03_EB, 7.184717e-03_EB, 1.764827e-01_EB, 2.574464e-04_EB, 1.205059e-03_EB, &
   5.689721e-04_EB/),(/7,8/))

gammad1_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 790-825 cm-1
   4.218735e-01_EB, 2.260530e-02_EB, 2.192833e-01_EB, 8.393035e-01_EB, 6.309253e-04_EB, 1.299978e-03_EB, &
   4.943829e-04_EB, &
   2.011501e+01_EB, 1.452133e-01_EB, 1.379684e+00_EB, 3.290149e+00_EB, 1.006802e-03_EB, 2.196075e-03_EB, &
   6.147171e-04_EB, &
   7.999647e+01_EB, 3.727151e-01_EB, 4.909546e+00_EB, 1.504204e-01_EB, 1.827877e-03_EB, 1.390294e-02_EB, &
   5.960009e-04_EB, &
   8.000000e+01_EB, 1.192733e+01_EB, 7.335305e+00_EB, 2.477027e+01_EB, 6.343012e-03_EB, 3.167933e-02_EB, &
   6.338454e-04_EB, &
   8.000000e+01_EB, 1.337438e+00_EB, 1.200630e+01_EB, 8.358805e-01_EB, 5.050388e-02_EB, 5.528025e-02_EB, &
   6.645315e-04_EB, &
   8.000000e+01_EB, 1.575943e+00_EB, 2.254763e+01_EB, 5.284652e-01_EB, 4.899995e-01_EB, 4.748211e-02_EB, &
   5.682371e-04_EB, &
   8.000000e+01_EB, 1.559461e+00_EB, 2.060048e+01_EB, 9.864427e-01_EB, 2.822989e+01_EB, 2.089522e-02_EB, &
   5.933048e-04_EB, &
   8.000000e+01_EB, 9.294253e-01_EB, 7.518592e-01_EB, 2.225028e+00_EB, 1.882489e+01_EB, 5.697081e-03_EB, &
   6.342437e-04_EB/),(/7,8/))

gammad1_c5h8o2(1:7,17:24) = RESHAPE((/ &  ! 830-865 cm-1
   8.000000e+01_EB, 5.595499e-01_EB, 7.229848e-01_EB, 4.663622e-01_EB, 4.280954e-02_EB, 3.407433e-03_EB, &
   7.293662e-04_EB, &
   7.999678e+01_EB, 5.297318e-01_EB, 1.408882e-01_EB, 5.990214e-01_EB, 2.104101e-03_EB, 5.806669e-03_EB, &
   5.977527e-04_EB, &
   8.000000e+01_EB, 2.382664e+00_EB, 5.857715e-02_EB, 2.596665e-01_EB, 5.172312e-04_EB, 1.035222e-03_EB, &
   5.002373e-04_EB, &
   5.199739e+01_EB, 2.075323e-01_EB, 1.580424e-02_EB, 9.101537e-03_EB, 3.465216e-04_EB, 4.105029e-04_EB, &
   6.238247e-04_EB, &
   2.523145e+00_EB, 7.869446e-04_EB, 6.466080e-04_EB, 5.984897e-03_EB, 2.719400e-04_EB, 9.547889e-04_EB, &
   6.744617e-04_EB, &
   3.144180e-01_EB, 6.967438e-04_EB, 8.338212e-04_EB, 6.247549e-04_EB, 2.369812e-04_EB, 2.160043e-01_EB, &
   1.146982e-03_EB, &
   3.556347e-01_EB, 5.790029e-04_EB, 9.280773e-04_EB, 6.140435e-04_EB, 1.166505e-04_EB, 4.285636e-03_EB, &
   1.647599e-01_EB, &
   9.611149e-01_EB, 6.993787e-04_EB, 6.536778e-04_EB, 6.470283e-04_EB, 6.300754e-04_EB, 6.508677e-04_EB, &
   8.394488e-01_EB/),(/7,8/))

gammad1_c5h8o2(1:7,25:27) = RESHAPE((/ &  ! 870-880 cm-1
   8.980663e-01_EB, 6.254884e-04_EB, 6.536611e-04_EB, 7.628907e-04_EB, 7.906692e-04_EB, 7.774676e-04_EB, &
   3.406991e+00_EB, &
   1.692425e-01_EB, 6.157700e-04_EB, 4.905940e-04_EB, 4.383728e-04_EB, 8.114566e-04_EB, 5.029333e-04_EB, &
   3.441873e-01_EB, &
   3.735791e-04_EB, 3.732776e-04_EB, 3.745026e-04_EB, 3.743404e-04_EB, 3.741743e-04_EB, 3.727339e-04_EB, &
   3.713571e-04_EB/),(/7,3/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_c5h8o2(n_temp_c5h8o2,37)) 

! band #2: 875 cm-1 - 1055 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 16.2578 % 

sd2_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 875-910 cm-1
   2.556429e-02_EB, 6.579116e-04_EB, 7.354099e-04_EB, 7.413989e-04_EB, 6.278952e-04_EB, 6.133453e-04_EB, &
   4.127227e-02_EB, &
   4.274772e-02_EB, 9.965553e-04_EB, 9.559177e-04_EB, 9.609587e-04_EB, 7.244017e-04_EB, 9.125642e-04_EB, &
   1.253746e-01_EB, &
   4.062606e-02_EB, 8.758771e-04_EB, 7.717414e-04_EB, 1.136983e-03_EB, 8.565555e-02_EB, 8.384111e-01_EB, &
   8.119652e-01_EB, &
   1.969402e-02_EB, 8.690397e-04_EB, 8.197383e-04_EB, 1.021204e-01_EB, 7.942200e-01_EB, 1.336575e+00_EB, &
   1.034711e+00_EB, &
   6.150252e-03_EB, 8.420963e-04_EB, 7.601258e-03_EB, 5.822673e-01_EB, 1.849210e+00_EB, 2.267047e-01_EB, &
   6.054709e-01_EB, &
   5.255431e-03_EB, 3.759946e-03_EB, 6.674130e-02_EB, 1.933685e-01_EB, 4.747476e+00_EB, 3.382479e-01_EB, &
   2.717558e-01_EB, &
   6.508425e-02_EB, 7.385905e-02_EB, 1.516138e-01_EB, 2.231123e-01_EB, 3.585642e+00_EB, 4.915855e-01_EB, &
   2.578738e-01_EB, &
   2.660815e-01_EB, 2.584587e-01_EB, 3.361944e-01_EB, 4.729659e-01_EB, 4.590394e+00_EB, 7.118124e-01_EB, &
   2.680924e-01_EB/),(/7,8/))

sd2_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 915-950 cm-1
   7.131926e-01_EB, 5.670265e-01_EB, 5.731113e-01_EB, 6.713483e-01_EB, 5.036802e+00_EB, 1.022991e+00_EB, &
   2.779647e-01_EB, &
   1.616994e+00_EB, 1.062612e+00_EB, 9.571218e-01_EB, 9.882622e-01_EB, 4.253259e+00_EB, 9.880129e-01_EB, &
   3.040041e-01_EB, &
   3.093637e+00_EB, 1.592723e+00_EB, 1.358023e+00_EB, 1.283095e+00_EB, 2.191786e+00_EB, 8.168738e-01_EB, &
   3.592899e-01_EB, &
   4.490720e+00_EB, 2.020538e+00_EB, 1.689196e+00_EB, 1.592450e+00_EB, 1.670002e+00_EB, 9.479666e-01_EB, &
   4.149322e-01_EB, &
   5.933235e+00_EB, 2.516838e+00_EB, 2.089816e+00_EB, 1.959977e+00_EB, 1.677729e+00_EB, 1.102576e+00_EB, &
   4.329370e-01_EB, &
   7.889794e+00_EB, 2.918784e+00_EB, 2.348617e+00_EB, 2.059360e+00_EB, 1.457291e+00_EB, 9.460046e-01_EB, &
   4.177430e-01_EB, &
   5.687576e+00_EB, 2.199276e+00_EB, 1.750801e+00_EB, 1.589524e+00_EB, 9.841957e-01_EB, 7.638453e-01_EB, &
   4.276039e-01_EB, &
   3.912971e+00_EB, 1.521382e+00_EB, 1.199146e+00_EB, 1.132591e+00_EB, 7.312452e-01_EB, 5.627116e-01_EB, &
   4.340964e-01_EB/),(/7,8/))

sd2_c5h8o2(1:7,17:24) = RESHAPE((/ &  ! 955-990 cm-1
   2.475675e+00_EB, 1.047082e+00_EB, 8.418699e-01_EB, 7.489753e-01_EB, 5.673221e-01_EB, 3.884845e-01_EB, &
   4.323566e-01_EB, &
   1.320393e+00_EB, 6.735843e-01_EB, 5.433291e-01_EB, 5.032234e-01_EB, 1.017197e+00_EB, 2.871344e-01_EB, &
   4.412980e-01_EB, &
   5.592797e-01_EB, 4.019616e-01_EB, 3.047795e-01_EB, 3.003325e-01_EB, 2.137201e+00_EB, 2.538634e-01_EB, &
   4.343141e-01_EB, &
   2.669833e-01_EB, 2.361810e-01_EB, 1.839138e-01_EB, 1.867551e-01_EB, 2.802026e+00_EB, 2.467341e-01_EB, &
   8.431414e-01_EB, &
   3.073555e-01_EB, 1.587138e-01_EB, 1.753897e-01_EB, 1.623549e-01_EB, 2.391695e+00_EB, 2.295551e-01_EB, &
   1.916065e+00_EB, &
   5.710281e-01_EB, 2.449444e-01_EB, 2.284392e-01_EB, 2.247888e-01_EB, 2.482041e+00_EB, 2.586464e-01_EB, &
   2.918096e+00_EB, &
   8.249745e-01_EB, 3.134360e-01_EB, 2.982082e-01_EB, 2.831453e-01_EB, 1.314234e+00_EB, 3.428023e-01_EB, &
   2.320833e+00_EB, &
   9.299117e-01_EB, 3.934637e-01_EB, 3.941485e-01_EB, 4.276205e-01_EB, 6.380472e-01_EB, 3.539805e-01_EB, &
   1.717216e+00_EB/),(/7,8/))

sd2_c5h8o2(1:7,25:32) = RESHAPE((/ &  ! 995-1030 cm-1
   1.174651e+00_EB, 5.076500e-01_EB, 5.123310e-01_EB, 5.072349e-01_EB, 6.017282e-01_EB, 4.040867e-01_EB, &
   9.071058e-01_EB, &
   1.220989e+00_EB, 6.322818e-01_EB, 6.487477e-01_EB, 6.488986e-01_EB, 5.856188e-01_EB, 4.402415e-01_EB, &
   6.313728e-01_EB, &
   1.342002e+00_EB, 7.952593e-01_EB, 7.945733e-01_EB, 8.003597e-01_EB, 7.389281e-01_EB, 4.724643e-01_EB, &
   1.867048e-01_EB, &
   2.059673e+00_EB, 1.016101e+00_EB, 9.573528e-01_EB, 9.061740e-01_EB, 6.905226e-01_EB, 5.789912e-01_EB, &
   1.181913e-01_EB, &
   2.708834e+00_EB, 1.261292e+00_EB, 1.056622e+00_EB, 9.704385e-01_EB, 7.470068e-01_EB, 3.592001e-01_EB, &
   1.004398e-01_EB, &
   3.415071e+00_EB, 1.258792e+00_EB, 1.017232e+00_EB, 8.840290e-01_EB, 7.155771e-01_EB, 2.368274e-01_EB, &
   1.095937e-01_EB, &
   3.371272e+00_EB, 1.241748e+00_EB, 9.412804e-01_EB, 7.974575e-01_EB, 8.067485e-01_EB, 1.865056e-01_EB, &
   1.178328e-01_EB, &
   2.426657e+00_EB, 8.931409e-01_EB, 6.775266e-01_EB, 5.586942e-01_EB, 6.967074e-01_EB, 1.154895e-01_EB, &
   1.127254e-01_EB/),(/7,8/))

sd2_c5h8o2(1:7,33:37) = RESHAPE((/ &  ! 1035-1055 cm-1
   1.800189e+00_EB, 5.884496e-01_EB, 4.478842e-01_EB, 3.611785e-01_EB, 6.015273e-01_EB, 2.659912e-02_EB, &
   1.137885e-01_EB, &
   6.470612e-01_EB, 2.347793e-01_EB, 1.941059e-01_EB, 1.585485e-01_EB, 1.218801e+00_EB, 6.925778e-04_EB, &
   6.825179e-01_EB, &
   1.671523e-01_EB, 4.958688e-02_EB, 7.155642e-02_EB, 4.212993e-02_EB, 4.481255e-02_EB, 8.744158e-04_EB, &
   6.166493e-01_EB, &
   4.210092e-02_EB, 5.872012e-03_EB, 1.202013e-02_EB, 8.673460e-03_EB, 7.247323e-04_EB, 7.563806e-04_EB, &
   1.496397e-01_EB, &
   4.683959e-04_EB, 4.682485e-04_EB, 4.694920e-04_EB, 4.693301e-04_EB, 4.691565e-04_EB, 4.675221e-04_EB, &
   4.659141e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_c5h8o2(n_temp_c5h8o2,37)) 

! band #2: 875 cm-1 - 1055 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 16.2578 % 

! print fine structure array gamma_d 

gammad2_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 875-910 cm-1
   2.823227e-01_EB, 4.610224e-04_EB, 4.906774e-04_EB, 5.417213e-04_EB, 5.214023e-04_EB, 4.619138e-04_EB, &
   2.996879e-01_EB, &
   1.589262e+01_EB, 7.112738e-04_EB, 6.863480e-04_EB, 6.954244e-04_EB, 5.636836e-04_EB, 6.909238e-04_EB, &
   1.362275e+00_EB, &
   1.255210e+01_EB, 6.773562e-04_EB, 5.735201e-04_EB, 8.272739e-04_EB, 1.724787e-04_EB, 1.661576e-04_EB, &
   1.138271e-03_EB, &
   7.051283e-01_EB, 6.891085e-04_EB, 7.024554e-04_EB, 1.678129e-04_EB, 1.524293e-04_EB, 5.547075e-04_EB, &
   1.418645e-03_EB, &
   7.121508e-02_EB, 6.384517e-04_EB, 1.642701e-03_EB, 1.694248e-04_EB, 1.470911e-04_EB, 6.320090e-02_EB, &
   4.163931e-03_EB, &
   2.508564e-03_EB, 4.326233e-03_EB, 1.052820e-02_EB, 5.592027e-03_EB, 2.571364e-04_EB, 3.075452e-02_EB, &
   7.237307e-02_EB, &
   3.527481e-01_EB, 6.803584e-02_EB, 4.286749e-02_EB, 3.844453e-02_EB, 5.858598e-04_EB, 2.259184e-02_EB, &
   3.898073e+00_EB, &
   1.571067e+01_EB, 6.225180e-02_EB, 6.038369e-02_EB, 2.741791e-02_EB, 1.292359e-03_EB, 2.289985e-02_EB, &
   1.760940e+01_EB/),(/7,8/))

gammad2_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 915-950 cm-1
   8.000000e+01_EB, 1.428430e-01_EB, 2.419731e+00_EB, 2.230224e-01_EB, 2.255022e-03_EB, 2.094570e-02_EB, &
   2.051410e+01_EB, &
   8.000000e+01_EB, 1.709830e-01_EB, 3.380268e+00_EB, 5.240909e-01_EB, 4.853055e-03_EB, 4.030274e-02_EB, &
   2.349415e+01_EB, &
   8.000000e+01_EB, 2.636763e-01_EB, 3.653768e+00_EB, 1.717666e+00_EB, 2.060690e-02_EB, 2.522187e-01_EB, &
   2.863894e+01_EB, &
   8.000000e+01_EB, 3.872301e-01_EB, 3.145524e+01_EB, 1.984056e+00_EB, 7.238493e-02_EB, 1.744457e-01_EB, &
   1.234446e+01_EB, &
   8.000000e+01_EB, 5.199554e-01_EB, 5.201416e+01_EB, 1.219270e+00_EB, 1.366230e-01_EB, 1.024375e-01_EB, &
   6.859851e+00_EB, &
   8.000000e+01_EB, 6.059595e-01_EB, 6.357106e+01_EB, 2.510435e+00_EB, 2.483816e-01_EB, 1.870400e-01_EB, &
   6.030901e+00_EB, &
   8.000000e+01_EB, 4.613948e-01_EB, 6.558012e+01_EB, 9.318435e-01_EB, 8.280329e-01_EB, 1.702565e-01_EB, &
   3.671889e+00_EB, &
   8.000000e+01_EB, 2.805617e-01_EB, 4.753473e+01_EB, 4.243871e-01_EB, 1.998842e-01_EB, 1.897322e-01_EB, &
   3.182135e-01_EB/),(/7,8/))

gammad2_c5h8o2(1:7,17:24) = RESHAPE((/ &  ! 955-990 cm-1
   8.000000e+01_EB, 1.917580e-01_EB, 3.160507e+01_EB, 3.013932e+00_EB, 8.610150e-02_EB, 2.680998e+00_EB, &
   1.225672e-01_EB, &
   8.000000e+01_EB, 8.181705e-02_EB, 1.393264e+01_EB, 4.846622e+00_EB, 5.364990e-03_EB, 7.410589e+00_EB, &
   7.715299e-02_EB, &
   7.999744e+01_EB, 2.071421e-02_EB, 2.646321e+00_EB, 1.308531e+00_EB, 1.066053e-03_EB, 5.993038e+00_EB, &
   7.459897e-02_EB, &
   1.600979e+01_EB, 6.164883e-03_EB, 3.149166e-01_EB, 1.578432e+00_EB, 6.172180e-04_EB, 3.880878e+00_EB, &
   7.684969e-03_EB, &
   2.587709e+01_EB, 3.007808e-02_EB, 1.472547e-01_EB, 1.384504e+00_EB, 6.759922e-04_EB, 4.013279e+00_EB, &
   2.292324e-03_EB, &
   7.999663e+01_EB, 3.784955e-02_EB, 4.354468e-01_EB, 2.323880e+00_EB, 8.249666e-04_EB, 2.529707e-01_EB, &
   1.448969e-03_EB, &
   8.000000e+01_EB, 6.887494e-02_EB, 3.383023e-01_EB, 3.534843e+00_EB, 2.269583e-03_EB, 7.251788e-02_EB, &
   1.620063e-03_EB, &
   8.000000e+01_EB, 1.396664e-01_EB, 2.500157e-01_EB, 1.154998e-01_EB, 1.324480e-02_EB, 1.665776e+00_EB, &
   1.569506e-03_EB/),(/7,8/))

gammad2_c5h8o2(1:7,25:32) = RESHAPE((/ &  ! 995-1030 cm-1
   7.999999e+01_EB, 2.622383e-01_EB, 3.680714e-01_EB, 3.431326e-01_EB, 3.214335e-02_EB, 3.028902e+00_EB, &
   2.505970e-03_EB, &
   8.000000e+01_EB, 2.931223e-01_EB, 4.763129e-01_EB, 5.158658e-01_EB, 1.564583e-01_EB, 4.068234e-01_EB, &
   3.470044e-03_EB, &
   8.000000e+01_EB, 2.591539e-01_EB, 9.760853e-01_EB, 3.725883e-01_EB, 6.472451e-02_EB, 7.215647e-02_EB, &
   3.438220e+00_EB, &
   8.000000e+01_EB, 2.874591e-01_EB, 8.902719e-01_EB, 3.997568e-01_EB, 1.273205e-01_EB, 2.106275e-02_EB, &
   1.404359e+00_EB, &
   8.000000e+01_EB, 2.627488e-01_EB, 2.860341e+01_EB, 4.718599e-01_EB, 6.177797e-02_EB, 7.129412e-02_EB, &
   1.080505e+00_EB, &
   8.000000e+01_EB, 3.044334e-01_EB, 1.893461e+01_EB, 7.512970e-01_EB, 3.525571e-02_EB, 8.428393e+00_EB, &
   1.655302e+00_EB, &
   8.000000e+01_EB, 2.042339e-01_EB, 2.401079e+01_EB, 1.842385e+00_EB, 1.678758e-02_EB, 1.160390e+01_EB, &
   7.709535e-02_EB, &
   8.000000e+01_EB, 1.532358e-01_EB, 2.398691e+01_EB, 7.959409e+00_EB, 8.237636e-03_EB, 1.840439e+00_EB, &
   1.574438e-01_EB/),(/7,8/))

gammad2_c5h8o2(1:7,33:37) = RESHAPE((/ &  ! 1035-1055 cm-1
   8.000000e+01_EB, 1.344383e-01_EB, 2.700694e+00_EB, 2.144500e+00_EB, 2.397408e-03_EB, 2.813670e-01_EB, &
   5.280750e-02_EB, &
   7.999478e+01_EB, 7.115998e-02_EB, 2.822221e-01_EB, 9.009493e-01_EB, 3.364759e-04_EB, 5.682134e-04_EB, &
   7.024865e-04_EB, &
   1.748846e+01_EB, 1.217636e-01_EB, 4.857488e-02_EB, 5.398096e-01_EB, 4.000099e-03_EB, 5.990465e-04_EB, &
   6.551947e-04_EB, &
   1.394167e+01_EB, 2.046399e-01_EB, 8.910555e-03_EB, 6.157030e-03_EB, 5.599763e-04_EB, 5.579403e-04_EB, &
   9.116727e-04_EB, &
   3.735791e-04_EB, 3.732776e-04_EB, 3.745026e-04_EB, 3.743404e-04_EB, 3.741743e-04_EB, 3.727339e-04_EB, &
   3.713571e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(sd3_c5h8o2(n_temp_c5h8o2,18)) 

! band #3: 1050 cm-1 - 1275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 38.3736 % 

sd3_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 1050-1085 cm-1
   2.283540e-02_EB, 4.922172e-04_EB, 2.560169e-03_EB, 3.717395e-03_EB, 5.679232e-04_EB, 5.948702e-04_EB, &
   9.070861e-02_EB, &
   1.330120e-02_EB, 4.682485e-04_EB, 6.438257e-03_EB, 1.249445e-02_EB, 1.158898e-03_EB, 8.034057e-04_EB, &
   3.528869e-01_EB, &
   5.312158e-04_EB, 4.682485e-04_EB, 3.462862e-03_EB, 4.802455e-03_EB, 1.325886e-03_EB, 9.860835e-04_EB, &
   5.966300e-01_EB, &
   6.730728e-04_EB, 5.383842e-04_EB, 1.541199e-02_EB, 1.216344e-03_EB, 8.711834e-04_EB, 2.095446e-02_EB, &
   2.344768e-01_EB, &
   4.842466e-03_EB, 5.362735e-03_EB, 2.344159e-02_EB, 1.759081e-02_EB, 9.888946e-04_EB, 5.883826e-02_EB, &
   1.524721e+00_EB, &
   1.962391e-02_EB, 1.842385e-02_EB, 4.296535e-02_EB, 3.232193e-02_EB, 8.684561e-04_EB, 1.277349e-01_EB, &
   1.441831e+00_EB, &
   3.957410e-02_EB, 3.609774e-02_EB, 7.027800e-02_EB, 6.897527e-02_EB, 1.805630e-02_EB, 1.769389e-01_EB, &
   2.219797e+00_EB, &
   6.033916e-02_EB, 5.989100e-02_EB, 1.135187e-01_EB, 1.032544e-01_EB, 9.054237e-02_EB, 2.674346e-01_EB, &
   1.717514e+00_EB/),(/7,8/))

sd3_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 1090-1225 cm-1
   1.000494e-01_EB, 9.414817e-02_EB, 1.444499e-01_EB, 1.465825e-01_EB, 1.695405e-01_EB, 3.353830e-01_EB, &
   6.307295e-01_EB, &
   1.689087e-01_EB, 1.393048e-01_EB, 1.815843e-01_EB, 1.798559e-01_EB, 3.961284e-01_EB, 3.870093e-01_EB, &
   5.212840e-01_EB, &
   2.316504e-01_EB, 1.905944e-01_EB, 2.231146e-01_EB, 2.373041e-01_EB, 6.108176e-01_EB, 5.517315e-01_EB, &
   5.514970e-01_EB, &
   7.470458e-01_EB, 7.308402e-01_EB, 9.133596e-01_EB, 1.083670e+00_EB, 1.643403e+00_EB, 2.654844e+00_EB, &
   2.192364e+00_EB, &
   9.479777e+00_EB, 8.554789e+00_EB, 7.931693e+00_EB, 8.343454e+00_EB, 7.636411e+00_EB, 6.546753e+00_EB, &
   3.098983e+00_EB, &
   2.594269e+01_EB, 1.307137e+01_EB, 1.100488e+01_EB, 1.045080e+01_EB, 8.049716e+00_EB, 4.871530e+00_EB, &
   1.930570e+00_EB, &
   1.934931e+01_EB, 7.224029e+00_EB, 5.871145e+00_EB, 5.292871e+00_EB, 3.621708e+00_EB, 1.992252e+00_EB, &
   7.777877e-01_EB, &
   4.140490e+00_EB, 2.139609e+00_EB, 1.821476e+00_EB, 1.658745e+00_EB, 1.493490e+00_EB, 7.032707e-01_EB, &
   2.294717e+00_EB/),(/7,8/))

sd3_c5h8o2(1:7,17:18) = RESHAPE((/ &  ! 1250-1275 cm-1
   5.271264e-01_EB, 2.912938e-01_EB, 2.917520e-01_EB, 2.643227e-01_EB, 5.514219e-01_EB, 1.728159e-01_EB, &
   9.935754e-01_EB, &
   4.683959e-04_EB, 4.682485e-04_EB, 4.694920e-04_EB, 4.693301e-04_EB, 4.691565e-04_EB, 4.675221e-04_EB, &
   4.659141e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad3_c5h8o2(n_temp_c5h8o2,18)) 

! band #3: 1050 cm-1 - 1275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 38.3736 % 

! print fine structure array gamma_d 

gammad3_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 1050-1085 cm-1
   6.752973e-01_EB, 3.965369e-04_EB, 1.852648e-03_EB, 6.404227e-03_EB, 4.395271e-04_EB, 4.775384e-04_EB, &
   2.331620e-03_EB, &
   2.308916e-01_EB, 3.732776e-04_EB, 1.932620e-03_EB, 2.619408e-03_EB, 7.975806e-04_EB, 6.318845e-04_EB, &
   1.556634e-03_EB, &
   4.424232e-04_EB, 3.732776e-04_EB, 2.729364e-03_EB, 1.540902e-03_EB, 9.180318e-04_EB, 7.423991e-04_EB, &
   8.482000e-04_EB, &
   5.798906e-04_EB, 4.109335e-04_EB, 4.684962e-03_EB, 1.748840e-03_EB, 5.924212e-04_EB, 8.917151e-02_EB, &
   3.547868e-03_EB, &
   2.445999e-01_EB, 1.709414e-03_EB, 8.476824e-03_EB, 1.070389e-02_EB, 6.275957e-04_EB, 4.118998e-01_EB, &
   5.365592e-04_EB, &
   2.903925e+00_EB, 1.089566e-02_EB, 2.453372e-02_EB, 1.655663e-01_EB, 5.991526e-04_EB, 4.446612e+00_EB, &
   5.688486e-04_EB, &
   1.701839e+01_EB, 1.480657e-01_EB, 7.478547e-02_EB, 3.425426e-01_EB, 6.291144e-01_EB, 9.461117e-01_EB, &
   4.469811e-04_EB, &
   6.073385e-01_EB, 1.511806e-01_EB, 1.786293e-01_EB, 1.777844e+00_EB, 9.776316e-01_EB, 7.402583e+00_EB, &
   1.041732e-03_EB/),(/7,8/))

gammad3_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 1090-1225 cm-1
   1.122342e+00_EB, 2.158159e-01_EB, 6.598921e-01_EB, 8.611218e-01_EB, 1.210188e-01_EB, 1.184061e+01_EB, &
   8.335214e-03_EB, &
   5.162273e+00_EB, 3.637833e-01_EB, 7.101749e-01_EB, 1.752143e+00_EB, 7.256519e-03_EB, 3.453023e+01_EB, &
   6.017312e-02_EB, &
   1.631154e+01_EB, 1.875513e-01_EB, 2.486087e+00_EB, 3.547150e+00_EB, 7.521178e-03_EB, 2.393089e+01_EB, &
   3.039081e+01_EB, &
   8.000000e+01_EB, 3.366139e-01_EB, 1.096627e+01_EB, 2.588432e+01_EB, 2.813109e-01_EB, 1.156638e+00_EB, &
   1.806175e-01_EB, &
   8.000000e+01_EB, 6.558841e-01_EB, 1.962907e+00_EB, 1.479647e+00_EB, 1.072565e+00_EB, 2.709515e+00_EB, &
   3.334761e-01_EB, &
   8.000000e+01_EB, 1.895544e+00_EB, 9.110170e+00_EB, 3.023831e+00_EB, 8.442824e-01_EB, 1.836944e+00_EB, &
   1.622069e-01_EB, &
   8.000000e+01_EB, 1.580597e+00_EB, 6.536920e+01_EB, 2.757392e+00_EB, 6.660127e-01_EB, 4.829505e+01_EB, &
   1.550614e-01_EB, &
   7.999934e+01_EB, 3.683494e-01_EB, 5.153029e+00_EB, 1.387986e+00_EB, 7.620403e-02_EB, 3.958904e+01_EB, &
   1.950315e-03_EB/),(/7,8/))

gammad3_c5h8o2(1:7,17:18) = RESHAPE((/ &  ! 1250-1275 cm-1
   7.999908e+01_EB, 9.980022e-02_EB, 3.152017e-01_EB, 1.039694e+00_EB, 3.550883e-03_EB, 1.085236e+00_EB, &
   8.696632e-04_EB, &
   3.735791e-04_EB, 3.732776e-04_EB, 3.745026e-04_EB, 3.743404e-04_EB, 3.741743e-04_EB, 3.727339e-04_EB, &
   3.713571e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(sd4_c5h8o2(n_temp_c5h8o2,14)) 

! band #4: 1250 cm-1 - 1575 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 26.5551 % 

sd4_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 1250-1425 cm-1
   2.651267e-01_EB, 1.759638e-01_EB, 1.769296e-01_EB, 1.699091e-01_EB, 2.046845e-01_EB, 2.341468e-01_EB, &
   1.105642e+00_EB, &
   1.022795e+00_EB, 9.266422e-01_EB, 1.017170e+00_EB, 1.118834e+00_EB, 1.488453e+00_EB, 1.472991e+00_EB, &
   1.015529e+00_EB, &
   9.008874e+00_EB, 4.815679e+00_EB, 4.305873e+00_EB, 4.164932e+00_EB, 3.844447e+00_EB, 2.320999e+00_EB, &
   1.064962e+00_EB, &
   1.317933e+01_EB, 5.235972e+00_EB, 4.254435e+00_EB, 3.735700e+00_EB, 2.736514e+00_EB, 1.279456e+00_EB, &
   5.251123e-01_EB, &
   3.335495e+00_EB, 1.313429e+00_EB, 1.062255e+00_EB, 8.454394e-01_EB, 1.978232e+00_EB, 2.277879e-01_EB, &
   3.643094e-02_EB, &
   1.208610e+00_EB, 5.576549e-01_EB, 5.119243e-01_EB, 4.375281e-01_EB, 2.260111e+00_EB, 2.451094e-01_EB, &
   9.103457e-02_EB, &
   1.929443e+00_EB, 8.888330e-01_EB, 8.250508e-01_EB, 7.384043e-01_EB, 2.462957e+00_EB, 6.121275e-01_EB, &
   1.297550e-01_EB, &
   2.558288e+00_EB, 1.354005e+00_EB, 1.257364e+00_EB, 1.190951e+00_EB, 1.431084e+00_EB, 8.234527e-01_EB, &
   1.945175e-01_EB/),(/7,8/))

sd4_c5h8o2(1:7,9:14) = RESHAPE((/ &  ! 1450-1575 cm-1
   7.269681e+00_EB, 2.804495e+00_EB, 2.300589e+00_EB, 2.018439e+00_EB, 1.754481e+00_EB, 8.665175e-01_EB, &
   1.851434e-01_EB, &
   3.164709e+00_EB, 1.268072e+00_EB, 1.073237e+00_EB, 9.062173e-01_EB, 1.779966e+00_EB, 4.279472e-01_EB, &
   6.749668e-02_EB, &
   4.203220e-01_EB, 2.293462e-01_EB, 2.285557e-01_EB, 1.794185e-01_EB, 1.212331e+00_EB, 7.133443e-02_EB, &
   9.727769e-04_EB, &
   3.205661e-01_EB, 1.369195e-01_EB, 1.279659e-01_EB, 9.148681e-02_EB, 2.184744e-01_EB, 3.056242e-03_EB, &
   9.133049e-04_EB, &
   1.379741e-01_EB, 6.041193e-02_EB, 5.056159e-02_EB, 3.540082e-02_EB, 4.510264e-03_EB, 3.305673e-03_EB, &
   6.976907e-04_EB, &
   4.683959e-04_EB, 4.682485e-04_EB, 4.694920e-04_EB, 4.693301e-04_EB, 4.691565e-04_EB, 4.675221e-04_EB, &
   4.659141e-04_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(gammad4_c5h8o2(n_temp_c5h8o2,14)) 

! band #4: 1250 cm-1 - 1575 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 26.5551 % 

! print fine structure array gamma_d 

gammad4_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 1250-1425 cm-1
   1.900794e+01_EB, 4.523774e-02_EB, 1.704620e-01_EB, 3.550652e-01_EB, 3.540509e-02_EB, 1.516513e+00_EB, &
   1.492599e-03_EB, &
   8.000000e+01_EB, 2.226090e-01_EB, 4.037264e+01_EB, 2.183364e+01_EB, 1.078975e-01_EB, 5.536881e-01_EB, &
   4.287697e-02_EB, &
   8.000000e+01_EB, 1.040937e+00_EB, 6.528460e+01_EB, 3.062031e+00_EB, 2.657658e-01_EB, 1.018869e+00_EB, &
   4.868239e-02_EB, &
   8.000000e+01_EB, 1.163412e+00_EB, 6.557988e+01_EB, 3.224408e+00_EB, 2.125073e-01_EB, 2.326661e-01_EB, &
   7.389632e-03_EB, &
   8.000000e+01_EB, 2.308729e-01_EB, 6.363109e-01_EB, 5.086314e+00_EB, 4.111352e-03_EB, 1.626691e-01_EB, &
   1.308121e+00_EB, &
   8.000000e+01_EB, 1.969078e-01_EB, 1.609102e+00_EB, 5.383790e+00_EB, 1.498115e-03_EB, 3.386863e-02_EB, &
   5.998714e-02_EB, &
   8.000000e+01_EB, 3.410180e-01_EB, 6.907342e+00_EB, 7.887927e+00_EB, 3.683657e-03_EB, 2.330242e-02_EB, &
   2.134694e+00_EB, &
   8.000000e+01_EB, 3.694747e-01_EB, 4.202343e+00_EB, 1.302414e+00_EB, 3.479140e-02_EB, 7.927698e-02_EB, &
   7.049147e+00_EB/),(/7,8/))

gammad4_c5h8o2(1:7,9:14) = RESHAPE((/ &  ! 1450-1575 cm-1
   8.000000e+01_EB, 7.141208e-01_EB, 6.536640e+01_EB, 2.342027e+00_EB, 7.034590e-02_EB, 1.139549e-01_EB, &
   5.068590e+00_EB, &
   8.000000e+01_EB, 3.439583e-01_EB, 2.158187e+00_EB, 1.101709e+00_EB, 6.060588e-03_EB, 2.185348e-02_EB, &
   3.183783e+00_EB, &
   5.006561e+01_EB, 3.299422e-01_EB, 1.190294e+00_EB, 8.974212e-01_EB, 3.345373e-04_EB, 2.520096e-02_EB, &
   1.139827e-03_EB, &
   2.293095e+01_EB, 1.421608e-01_EB, 7.208656e-01_EB, 9.099123e-01_EB, 6.014407e-05_EB, 2.971282e-03_EB, &
   6.193551e-04_EB, &
   2.442962e+00_EB, 4.484537e-02_EB, 1.690359e-01_EB, 1.259316e-01_EB, 9.393127e-04_EB, 1.423587e-03_EB, &
   5.673518e-04_EB, &
   3.735791e-04_EB, 3.732776e-04_EB, 3.745026e-04_EB, 3.743404e-04_EB, 3.741743e-04_EB, 3.727339e-04_EB, &
   3.713571e-04_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(sd5_c5h8o2(n_temp_c5h8o2,18)) 

! band #5: 1550 cm-1 - 1975 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 33.9665 % 

sd5_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 1550-1725 cm-1
   1.283781e-01_EB, 5.564827e-02_EB, 4.552881e-02_EB, 3.248586e-02_EB, 2.962603e-03_EB, 5.775666e-03_EB, &
   6.549796e-04_EB, &
   3.228302e-01_EB, 1.330461e-01_EB, 1.436585e-01_EB, 1.010144e-01_EB, 3.539658e+00_EB, 7.633235e-01_EB, &
   3.006478e-03_EB, &
   4.297863e-01_EB, 2.130911e-01_EB, 3.958264e-01_EB, 2.086835e-01_EB, 5.228719e+00_EB, 1.073793e+00_EB, &
   5.259787e-02_EB, &
   1.636030e+00_EB, 7.901740e-01_EB, 7.307309e-01_EB, 6.594717e-01_EB, 1.967717e+00_EB, 3.028161e-01_EB, &
   7.555184e-02_EB, &
   2.786162e+00_EB, 9.811787e-01_EB, 7.947490e-01_EB, 6.679210e-01_EB, 3.937599e-01_EB, 1.926412e-01_EB, &
   4.895888e-02_EB, &
   5.371774e-01_EB, 2.528496e-01_EB, 2.792431e-01_EB, 2.552566e-01_EB, 4.163329e-01_EB, 2.237360e-01_EB, &
   1.885545e-01_EB, &
   1.178231e+00_EB, 6.416768e-01_EB, 6.572166e-01_EB, 7.189543e-01_EB, 7.686419e-01_EB, 9.688893e-01_EB, &
   6.754928e-01_EB, &
   7.583351e+00_EB, 5.772869e+00_EB, 5.148166e+00_EB, 5.223398e+00_EB, 4.376104e+00_EB, 3.382685e+00_EB, &
   1.483712e+00_EB/),(/7,8/))

sd5_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 1750-1925 cm-1
   2.394789e+01_EB, 1.094859e+01_EB, 8.981291e+00_EB, 8.316621e+00_EB, 5.789605e+00_EB, 3.531464e+00_EB, &
   1.436978e+00_EB, &
   1.830194e+00_EB, 1.308594e+00_EB, 1.243723e+00_EB, 1.234798e+00_EB, 1.062593e+00_EB, 8.475205e-01_EB, &
   4.802389e-01_EB, &
   6.055532e-01_EB, 3.245427e-01_EB, 3.283205e-01_EB, 2.996324e-01_EB, 2.752568e-01_EB, 1.605966e-01_EB, &
   1.734273e-01_EB, &
   1.185795e-01_EB, 5.910661e-02_EB, 5.293524e-02_EB, 4.089505e-02_EB, 8.894114e-02_EB, 1.564737e-02_EB, &
   4.355749e-02_EB, &
   4.355125e-02_EB, 1.987367e-02_EB, 8.799089e-02_EB, 2.271459e-02_EB, 2.895065e-02_EB, 3.488458e-03_EB, &
   2.838444e-02_EB, &
   3.376642e-01_EB, 1.368154e-01_EB, 1.357689e-01_EB, 1.064263e-01_EB, 6.258461e-02_EB, 1.334385e-02_EB, &
   8.818944e-03_EB, &
   2.105506e-01_EB, 6.949000e-02_EB, 6.742734e-02_EB, 5.316113e-02_EB, 2.371958e-02_EB, 9.426557e-04_EB, &
   1.609069e-03_EB, &
   1.850030e-02_EB, 4.412077e-03_EB, 7.600489e-04_EB, 2.224295e-03_EB, 7.330402e-04_EB, 7.865066e-04_EB, &
   7.116352e-04_EB/),(/7,8/))

sd5_c5h8o2(1:7,17:18) = RESHAPE((/ &  ! 1950-1975 cm-1
   9.010576e-04_EB, 6.552187e-04_EB, 7.003811e-04_EB, 6.830920e-04_EB, 7.111960e-04_EB, 6.441936e-04_EB, &
   6.424584e-04_EB, &
   4.683959e-04_EB, 4.682485e-04_EB, 4.694920e-04_EB, 4.693301e-04_EB, 4.691565e-04_EB, 4.675221e-04_EB, &
   4.659141e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad5_c5h8o2(n_temp_c5h8o2,18)) 

! band #5: 1550 cm-1 - 1975 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 33.9665 % 

! print fine structure array gamma_d 

gammad5_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 1550-1725 cm-1
   2.001538e+00_EB, 1.043359e-01_EB, 1.744313e-01_EB, 5.342194e-02_EB, 1.360786e-03_EB, 2.592801e-03_EB, &
   5.075444e-04_EB, &
   1.725159e+01_EB, 1.172078e-01_EB, 4.594855e-02_EB, 2.031343e-01_EB, 1.111304e-04_EB, 2.683498e-04_EB, &
   1.034216e-01_EB, &
   4.985857e+01_EB, 1.394095e-01_EB, 7.496345e-03_EB, 1.380074e-01_EB, 3.724253e-04_EB, 1.238566e-03_EB, &
   7.511094e-01_EB, &
   7.999340e+01_EB, 3.572097e-01_EB, 1.449507e+00_EB, 3.819534e+00_EB, 3.386923e-03_EB, 2.420339e+00_EB, &
   1.689780e+00_EB, &
   8.000000e+01_EB, 4.078442e-01_EB, 2.746087e+01_EB, 5.496911e+00_EB, 7.681082e-02_EB, 9.142768e+00_EB, &
   9.256224e-01_EB, &
   8.000000e+01_EB, 1.120101e+00_EB, 2.154457e-01_EB, 2.496980e-01_EB, 5.035481e-03_EB, 4.567797e+00_EB, &
   1.613476e+01_EB, &
   8.000000e+01_EB, 2.505163e-01_EB, 7.289726e+00_EB, 3.688168e-01_EB, 1.262311e-01_EB, 4.852373e+01_EB, &
   4.327932e+01_EB, &
   8.000000e+01_EB, 4.925207e-01_EB, 1.508704e+00_EB, 1.040179e+00_EB, 6.046273e-01_EB, 2.007314e+01_EB, &
   3.883409e+01_EB/),(/7,8/))

gammad5_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 1750-1925 cm-1
   8.000000e+01_EB, 1.703895e+00_EB, 9.191951e+00_EB, 2.709061e+00_EB, 9.359203e-01_EB, 3.834274e+01_EB, &
   7.078759e-01_EB, &
   7.999983e+01_EB, 1.728338e-01_EB, 6.428302e-01_EB, 6.480871e-01_EB, 2.403467e-01_EB, 4.865268e+01_EB, &
   1.451594e+01_EB, &
   7.999999e+01_EB, 7.815232e-02_EB, 8.730432e-02_EB, 2.873412e-01_EB, 2.837574e-02_EB, 3.159028e+00_EB, &
   3.214034e+00_EB, &
   6.362515e-02_EB, 6.189950e-03_EB, 5.670266e-03_EB, 1.609798e-02_EB, 4.919947e-04_EB, 1.569334e-01_EB, &
   1.400792e-01_EB, &
   4.873130e-03_EB, 6.913205e-03_EB, 1.797921e-04_EB, 1.336251e-02_EB, 8.086060e-03_EB, 1.891584e-02_EB, &
   2.687108e-02_EB, &
   2.972716e+01_EB, 4.496649e-02_EB, 1.698626e-02_EB, 7.124485e-02_EB, 3.445376e-01_EB, 2.327743e-01_EB, &
   6.610266e-01_EB, &
   7.398272e+00_EB, 4.698677e-02_EB, 1.308793e-02_EB, 3.824701e-02_EB, 6.967306e-02_EB, 1.066862e-03_EB, &
   8.184279e-04_EB, &
   9.082357e-03_EB, 3.297959e-03_EB, 4.284375e-04_EB, 1.604779e-03_EB, 5.294588e-04_EB, 6.058801e-04_EB, &
   5.178529e-04_EB/),(/7,8/))

gammad5_c5h8o2(1:7,17:18) = RESHAPE((/ &  ! 1950-1975 cm-1
   7.219679e-04_EB, 5.013825e-04_EB, 4.812947e-04_EB, 5.404546e-04_EB, 5.475653e-04_EB, 4.887417e-04_EB, &
   4.941073e-04_EB, &
   3.735791e-04_EB, 3.732776e-04_EB, 3.745026e-04_EB, 3.743404e-04_EB, 3.741743e-04_EB, 3.727339e-04_EB, &
   3.713571e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(sd6_c5h8o2(n_temp_c5h8o2,26)) 

! band #6: 2650 cm-1 - 3275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 11.3243 % 

sd6_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 2650-2825 cm-1
   6.617679e-04_EB, 7.038166e-04_EB, 5.823242e-04_EB, 6.498762e-04_EB, 6.305493e-04_EB, 6.702982e-04_EB, &
   2.492250e-02_EB, &
   8.413829e-04_EB, 9.287068e-04_EB, 9.382552e-04_EB, 8.391021e-04_EB, 8.961598e-04_EB, 8.253163e-04_EB, &
   5.284118e-02_EB, &
   9.964645e-03_EB, 3.825227e-03_EB, 5.258108e-03_EB, 3.430243e-03_EB, 8.600200e-04_EB, 6.425559e-03_EB, &
   7.468932e-02_EB, &
   4.273055e-02_EB, 1.049000e-02_EB, 1.562321e-02_EB, 1.351817e-02_EB, 7.806216e-03_EB, 2.714191e-02_EB, &
   1.221449e-01_EB, &
   8.007875e-02_EB, 2.912543e-02_EB, 3.379663e-02_EB, 3.252406e-02_EB, 1.200427e-02_EB, 2.850065e-02_EB, &
   1.385424e-01_EB, &
   9.973924e-02_EB, 2.113446e-02_EB, 3.076355e-02_EB, 2.639593e-02_EB, 6.644964e-03_EB, 2.242511e-02_EB, &
   1.660444e-01_EB, &
   1.028910e-01_EB, 5.220431e-02_EB, 4.300251e-02_EB, 4.795173e-02_EB, 3.143119e-02_EB, 6.642097e-02_EB, &
   1.979488e-01_EB, &
   2.202342e-01_EB, 1.345745e-01_EB, 1.259318e-01_EB, 1.271074e-01_EB, 1.238737e-01_EB, 1.520918e-01_EB, &
   2.525918e-01_EB/),(/7,8/))

sd6_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 2850-3025 cm-1
   1.136213e+00_EB, 4.781891e-01_EB, 4.249774e-01_EB, 3.828941e-01_EB, 3.367838e-01_EB, 2.329863e-01_EB, &
   2.717789e-01_EB, &
   8.634086e-01_EB, 4.219181e-01_EB, 3.728310e-01_EB, 3.478836e-01_EB, 2.675079e-01_EB, 2.733947e-01_EB, &
   3.309081e-01_EB, &
   1.146984e+00_EB, 5.316939e-01_EB, 5.062337e-01_EB, 4.805370e-01_EB, 3.860685e-01_EB, 3.811075e-01_EB, &
   4.090385e-01_EB, &
   2.056434e+00_EB, 9.985803e-01_EB, 9.030141e-01_EB, 8.563341e-01_EB, 6.732570e-01_EB, 6.443301e-01_EB, &
   5.434919e-01_EB, &
   4.615588e+00_EB, 2.244759e+00_EB, 2.030504e+00_EB, 1.928204e+00_EB, 1.523456e+00_EB, 1.179250e+00_EB, &
   7.021791e-01_EB, &
   5.650492e+00_EB, 2.472471e+00_EB, 2.136749e+00_EB, 1.930516e+00_EB, 1.442591e+00_EB, 9.997534e-01_EB, &
   5.978827e-01_EB, &
   4.609191e+00_EB, 1.782318e+00_EB, 1.525555e+00_EB, 1.347833e+00_EB, 9.229266e-01_EB, 6.392344e-01_EB, &
   4.770802e-01_EB, &
   2.820790e+00_EB, 1.104460e+00_EB, 9.776057e-01_EB, 8.631595e-01_EB, 5.721892e-01_EB, 4.422352e-01_EB, &
   3.522752e-01_EB/),(/7,8/))

sd6_c5h8o2(1:7,17:24) = RESHAPE((/ &  ! 3050-3225 cm-1
   1.187353e+00_EB, 5.460448e-01_EB, 4.928400e-01_EB, 4.472267e-01_EB, 3.395637e-01_EB, 2.809475e-01_EB, &
   2.653517e-01_EB, &
   4.704535e-01_EB, 2.491073e-01_EB, 2.537560e-01_EB, 2.409873e-01_EB, 1.982473e-01_EB, 2.279718e-01_EB, &
   2.611720e-01_EB, &
   6.039772e-01_EB, 3.431934e-01_EB, 3.284617e-01_EB, 3.129996e-01_EB, 2.653921e-01_EB, 2.526170e-01_EB, &
   2.011828e-01_EB, &
   5.920773e-01_EB, 2.451875e-01_EB, 2.368206e-01_EB, 2.072703e-01_EB, 2.138315e-01_EB, 1.682738e-01_EB, &
   1.552495e-01_EB, &
   1.700742e-01_EB, 7.290921e-02_EB, 8.217440e-02_EB, 6.629481e-02_EB, 5.469629e-02_EB, 8.994315e-02_EB, &
   9.732156e-02_EB, &
   1.118992e-01_EB, 4.366084e-02_EB, 4.890966e-02_EB, 2.885185e-02_EB, 1.279711e-02_EB, 6.143421e-02_EB, &
   5.932504e-02_EB, &
   5.320012e-02_EB, 1.864948e-02_EB, 1.833215e-02_EB, 8.372681e-03_EB, 9.556242e-04_EB, 2.724051e-02_EB, &
   4.060600e-02_EB, &
   1.562049e-02_EB, 9.377158e-03_EB, 1.509016e-03_EB, 2.400877e-03_EB, 7.719203e-04_EB, 4.833588e-03_EB, &
   3.145233e-03_EB/),(/7,8/))

sd6_c5h8o2(1:7,25:26) = RESHAPE((/ &  ! 3250-3275 cm-1
   1.234745e-03_EB, 1.894526e-03_EB, 7.790315e-04_EB, 6.688577e-04_EB, 6.556824e-04_EB, 6.274349e-04_EB, &
   6.008281e-04_EB, &
   4.683959e-04_EB, 4.682485e-04_EB, 4.694920e-04_EB, 4.693301e-04_EB, 4.691565e-04_EB, 4.675221e-04_EB, &
   4.659141e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad6_c5h8o2(n_temp_c5h8o2,26)) 

! band #6: 2650 cm-1 - 3275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 11.3243 % 

! print fine structure array gamma_d 

gammad6_c5h8o2(1:7,1:8) = RESHAPE((/ &  ! 2650-2825 cm-1
   5.202784e-04_EB, 5.525656e-04_EB, 4.434412e-04_EB, 4.998255e-04_EB, 4.935419e-04_EB, 5.063850e-04_EB, &
   1.092982e-02_EB, &
   5.613238e-04_EB, 6.098207e-04_EB, 7.086604e-04_EB, 6.039142e-04_EB, 6.289973e-04_EB, 6.322323e-04_EB, &
   1.194772e+00_EB, &
   2.420843e-03_EB, 1.761218e-03_EB, 3.068546e-03_EB, 3.918268e-03_EB, 6.421563e-04_EB, 2.061317e-03_EB, &
   4.004519e-01_EB, &
   8.686816e-04_EB, 5.690975e-03_EB, 1.368498e-01_EB, 2.276203e-01_EB, 1.896994e-03_EB, 1.942151e-03_EB, &
   1.781156e+00_EB, &
   8.775885e-02_EB, 2.065817e-02_EB, 1.307324e-01_EB, 4.176628e-02_EB, 1.462109e-02_EB, 2.246227e-03_EB, &
   1.423947e+00_EB, &
   1.400226e-03_EB, 7.669502e-03_EB, 2.007661e-01_EB, 5.908664e-02_EB, 9.431682e-02_EB, 1.143177e-02_EB, &
   1.718195e+00_EB, &
   1.151669e-01_EB, 1.850972e-02_EB, 1.407721e-01_EB, 1.023565e-01_EB, 3.200614e-02_EB, 1.001680e-01_EB, &
   2.390609e+00_EB, &
   6.151837e+00_EB, 6.526566e-02_EB, 8.225336e-01_EB, 7.234620e-01_EB, 1.810297e-02_EB, 6.092961e-01_EB, &
   6.507145e+00_EB/),(/7,8/))

gammad6_c5h8o2(1:7,9:16) = RESHAPE((/ &  ! 2850-3025 cm-1
   7.999998e+01_EB, 4.402301e-01_EB, 8.184655e+00_EB, 7.049003e+00_EB, 2.465981e-02_EB, 4.666996e-01_EB, &
   1.050980e+01_EB, &
   8.000000e+01_EB, 1.344989e-01_EB, 5.238427e+00_EB, 2.779348e+00_EB, 5.844913e-02_EB, 2.117354e-01_EB, &
   4.746041e+00_EB, &
   7.999999e+01_EB, 1.706723e+00_EB, 1.717115e+01_EB, 7.059894e+00_EB, 8.912304e-02_EB, 2.663843e-01_EB, &
   1.132361e+01_EB, &
   7.999995e+01_EB, 4.086512e-01_EB, 4.839686e+01_EB, 2.308238e+01_EB, 2.068307e-01_EB, 8.659772e+00_EB, &
   1.089381e+00_EB, &
   8.000000e+01_EB, 7.844302e-01_EB, 6.546691e+01_EB, 3.816010e+00_EB, 3.501117e-01_EB, 4.860275e+01_EB, &
   2.192674e+01_EB, &
   8.000000e+01_EB, 7.277374e-01_EB, 6.563413e+01_EB, 6.128970e+01_EB, 2.795272e-01_EB, 4.863825e+01_EB, &
   1.715815e+01_EB, &
   8.000000e+01_EB, 7.372699e-01_EB, 6.564590e+01_EB, 5.533248e+01_EB, 2.910770e-01_EB, 9.688164e+00_EB, &
   2.950219e+01_EB, &
   8.000000e+01_EB, 5.982188e-01_EB, 6.554359e+01_EB, 3.587869e+01_EB, 3.437982e-01_EB, 2.068358e+00_EB, &
   2.249578e+01_EB/),(/7,8/))

gammad6_c5h8o2(1:7,17:24) = RESHAPE((/ &  ! 3050-3225 cm-1
   8.000000e+01_EB, 4.235551e-01_EB, 2.497231e+01_EB, 9.302322e+00_EB, 8.931722e-02_EB, 5.150630e+00_EB, &
   1.414719e+01_EB, &
   6.010851e+01_EB, 1.879884e+00_EB, 2.783897e+00_EB, 3.173126e+00_EB, 8.402729e-02_EB, 3.467785e+00_EB, &
   1.176859e+01_EB, &
   7.499988e+01_EB, 1.744287e-01_EB, 2.563845e+00_EB, 4.952657e+00_EB, 6.326135e-02_EB, 3.738165e+00_EB, &
   8.586505e+00_EB, &
   7.999518e+01_EB, 3.410346e-01_EB, 4.666935e-01_EB, 2.152000e+00_EB, 7.943992e-03_EB, 9.823071e-01_EB, &
   3.227657e+00_EB, &
   6.470047e-01_EB, 2.696803e-01_EB, 4.931105e-02_EB, 1.777548e-01_EB, 6.548798e-03_EB, 4.997383e-02_EB, &
   2.368120e+00_EB, &
   1.850193e-02_EB, 2.343943e-01_EB, 2.130779e-02_EB, 1.645282e-01_EB, 5.228779e-03_EB, 1.024445e-01_EB, &
   3.884470e-01_EB, &
   7.638521e-03_EB, 1.078876e-01_EB, 1.679201e-02_EB, 8.467687e-03_EB, 6.767255e-04_EB, 1.893862e-01_EB, &
   1.745921e-01_EB, &
   9.543385e-03_EB, 5.137643e-04_EB, 8.369605e-04_EB, 1.442443e-03_EB, 5.585763e-04_EB, 6.679236e-03_EB, &
   1.904481e-03_EB/),(/7,8/))

gammad6_c5h8o2(1:7,25:26) = RESHAPE((/ &  ! 3250-3275 cm-1
   4.938200e-04_EB, 5.333619e-04_EB, 5.374406e-04_EB, 5.112562e-04_EB, 5.068363e-04_EB, 4.776250e-04_EB, &
   4.658301e-04_EB, &
   3.735791e-04_EB, 3.732776e-04_EB, 3.745026e-04_EB, 3.743404e-04_EB, 3.741743e-04_EB, 3.727339e-04_EB, &
   3.713571e-04_EB/),(/7,2/))

!-------------------------propane data-------------------


! there are 2 bands for propane

! band #1: 1175 cm-1 - 1675 cm-1 
! band #2: 2550 cm-1 - 3375 cm-1 

ALLOCATE(sd_c3h8_temp(n_temp_c3h8)) 

! initialize bands wavenumber bounds for propane ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_c3h8(n_band_c3h8,3)) 
ALLOCATE(be_c3h8(n_band_c3h8)) 

om_bnd_c3h8 = RESHAPE((/ &
   1175._EB, 2550._EB, &
   1675._EB, 3375._EB, &
   25._EB, 25._EB/),(/n_band_c3h8,3/)) 

sd_c3h8_temp = (/ &
   295._EB, 396._EB, 435._EB, 513._EB, 578._EB, 790._EB,&
   1009._EB/)

be_c3h8 = (/ &
   1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_c3h8(n_temp_c3h8,21)) 

! band #1: 1175 cm-1 - 1675 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.76397 % 

sd1_c3h8(1:7,1:8) = RESHAPE((/ &  ! 1175-1350 cm-1
   2.856475e-02_EB, 1.979289e-02_EB, 1.707802e-02_EB, 1.759506e-02_EB, 5.992940e-04_EB, 6.107200e-04_EB, &
   6.336991e-04_EB, &
   1.484719e-02_EB, 1.365281e-02_EB, 1.474969e-02_EB, 1.061569e+00_EB, 1.945253e-03_EB, 1.757434e-03_EB, &
   4.662084e-03_EB, &
   6.963764e-04_EB, 8.500329e-04_EB, 1.107148e-02_EB, 4.507160e-01_EB, 1.264992e-03_EB, 2.619320e-03_EB, &
   7.816939e-01_EB, &
   7.226016e-04_EB, 9.075238e-04_EB, 1.125338e-02_EB, 2.194527e+00_EB, 3.638673e-03_EB, 7.038039e-03_EB, &
   8.347445e-03_EB, &
   7.167376e-04_EB, 2.529382e-03_EB, 1.471445e-02_EB, 5.689394e-01_EB, 5.129735e-03_EB, 3.199170e-01_EB, &
   4.714861e-02_EB, &
   1.973741e-02_EB, 3.150872e-02_EB, 4.659984e-02_EB, 9.243516e-01_EB, 5.980704e-02_EB, 5.046688e-01_EB, &
   2.447019e+00_EB, &
   1.692269e-01_EB, 1.375209e-01_EB, 1.214367e-01_EB, 4.508625e-01_EB, 8.448186e-02_EB, 1.539041e-01_EB, &
   2.379877e+00_EB, &
   3.080345e-01_EB, 2.355289e-01_EB, 2.170087e-01_EB, 2.766257e-01_EB, 1.745191e-01_EB, 3.225558e-01_EB, &
   2.527352e-01_EB/),(/7,8/))

sd1_c3h8(1:7,9:16) = RESHAPE((/ &  ! 1375-1550 cm-1
   5.345587e-01_EB, 3.413403e-01_EB, 2.898386e-01_EB, 4.388584e-01_EB, 1.966052e-01_EB, 2.312638e-01_EB, &
   1.022691e+00_EB, &
   4.993827e-01_EB, 3.759253e-01_EB, 3.463795e-01_EB, 4.015651e-01_EB, 2.676596e-01_EB, 2.708667e-01_EB, &
   5.754884e-01_EB, &
   3.709223e-01_EB, 3.014857e-01_EB, 3.021999e-01_EB, 9.358461e-01_EB, 2.624617e-01_EB, 2.519600e-01_EB, &
   3.362995e-01_EB, &
   8.886871e-01_EB, 5.747695e-01_EB, 5.132261e-01_EB, 6.330728e-01_EB, 3.344818e-01_EB, 2.389406e-01_EB, &
   1.497565e+00_EB, &
   1.197095e+00_EB, 6.989701e-01_EB, 6.042183e-01_EB, 8.699088e-01_EB, 3.814121e-01_EB, 2.349714e-01_EB, &
   2.153495e-01_EB, &
   6.712891e-01_EB, 3.823312e-01_EB, 3.697114e-01_EB, 8.758687e-01_EB, 2.396388e-01_EB, 1.894728e-01_EB, &
   1.039864e+00_EB, &
   1.822239e-01_EB, 1.242199e-01_EB, 1.389102e-01_EB, 1.652805e+00_EB, 1.094348e-01_EB, 9.737249e-02_EB, &
   6.800781e-02_EB, &
   1.683063e-02_EB, 9.880390e-03_EB, 3.283471e-02_EB, 4.035590e+00_EB, 3.286196e-02_EB, 7.512521e-02_EB, &
   1.698967e+00_EB/),(/7,8/))

sd1_c3h8(1:7,17:21) = RESHAPE((/ &  ! 1575-1675 cm-1
   4.515698e-03_EB, 1.676975e-02_EB, 9.512287e-03_EB, 1.896304e+00_EB, 1.105889e-02_EB, 3.093088e-01_EB, &
   1.288482e+00_EB, &
   9.765260e-04_EB, 6.840034e-01_EB, 1.034854e-02_EB, 3.287022e-03_EB, 7.493779e-03_EB, 2.811245e-01_EB, &
   4.924827e-02_EB, &
   6.957537e-04_EB, 2.871557e-03_EB, 6.957746e-04_EB, 8.150887e-04_EB, 1.138683e-03_EB, 4.575742e-03_EB, &
   1.445285e-02_EB, &
   5.346547e-04_EB, 6.489883e-04_EB, 5.962380e-04_EB, 5.749846e-04_EB, 5.575030e-04_EB, 6.786480e-04_EB, &
   6.619676e-04_EB, &
   4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, &
   4.634348e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_c3h8(n_temp_c3h8,21)) 

! band #1: 1175 cm-1 - 1675 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.76397 % 

! print fine structure array gamma_d 

gammad1_c3h8(1:7,1:8) = RESHAPE((/ &  ! 1175-1350 cm-1
   1.992031e-02_EB, 1.141762e-02_EB, 1.366155e-01_EB, 1.973038e-02_EB, 4.249794e-04_EB, 4.342393e-04_EB, &
   4.468799e-04_EB, &
   8.577876e-03_EB, 6.212098e-03_EB, 1.781483e-01_EB, 9.603555e-05_EB, 4.057086e-03_EB, 1.500925e-04_EB, &
   6.183171e-04_EB, &
   4.465589e-04_EB, 5.682785e-04_EB, 2.127940e-01_EB, 5.658019e-05_EB, 2.225962e-04_EB, 3.326833e-04_EB, &
   5.709512e-06_EB, &
   5.278412e-04_EB, 5.926586e-04_EB, 2.442104e-01_EB, 3.443200e-05_EB, 3.645364e-04_EB, 1.452036e-03_EB, &
   7.864828e-02_EB, &
   5.090058e-04_EB, 2.365350e-03_EB, 6.569843e-01_EB, 4.564938e-05_EB, 6.851756e-03_EB, 3.833493e-05_EB, &
   6.658254e-01_EB, &
   1.552597e-02_EB, 2.378523e-02_EB, 8.041370e-01_EB, 1.425241e-04_EB, 5.497853e-03_EB, 2.341884e-04_EB, &
   8.855023e-05_EB, &
   1.117621e-01_EB, 2.577861e-02_EB, 1.116130e+00_EB, 1.259624e-03_EB, 1.821509e-01_EB, 3.535192e-03_EB, &
   8.671143e-05_EB, &
   1.808233e-01_EB, 6.199510e-02_EB, 4.856790e+00_EB, 1.205230e-02_EB, 3.669786e-02_EB, 2.932383e-03_EB, &
   4.381239e-03_EB/),(/7,8/))

gammad1_c3h8(1:7,9:16) = RESHAPE((/ &  ! 1375-1550 cm-1
   2.669661e-01_EB, 1.045895e-01_EB, 1.300781e+00_EB, 9.033578e-03_EB, 2.323335e-01_EB, 1.276771e-02_EB, &
   8.417636e-04_EB, &
   3.073214e-01_EB, 1.882989e-01_EB, 1.392323e+01_EB, 2.514566e-02_EB, 1.753049e-01_EB, 1.668277e-02_EB, &
   1.564119e-03_EB, &
   1.834951e-01_EB, 1.904052e+00_EB, 1.048436e+01_EB, 3.971260e-03_EB, 8.256188e-02_EB, 1.231218e-02_EB, &
   4.278796e-03_EB, &
   4.067671e-01_EB, 3.098373e-01_EB, 1.683329e+01_EB, 2.321737e-02_EB, 4.469660e+00_EB, 9.982721e-02_EB, &
   5.564305e-04_EB, &
   4.351927e-01_EB, 3.169913e-01_EB, 3.797676e-01_EB, 1.482328e-02_EB, 1.277859e-01_EB, 9.112234e-02_EB, &
   7.464685e-03_EB, &
   3.429172e-01_EB, 1.903809e+01_EB, 1.993602e+01_EB, 4.962290e-03_EB, 2.386361e-01_EB, 1.458725e-02_EB, &
   3.066693e-04_EB, &
   1.761592e-01_EB, 6.137305e+00_EB, 2.041270e+00_EB, 3.861036e-04_EB, 1.018629e+00_EB, 1.406376e-02_EB, &
   1.289657e-02_EB, &
   2.874081e-02_EB, 4.620945e-01_EB, 1.299484e+00_EB, 3.273246e-05_EB, 8.688233e-02_EB, 4.030533e-04_EB, &
   2.700077e-05_EB/),(/7,8/))

gammad1_c3h8(1:7,17:21) = RESHAPE((/ &  ! 1575-1675 cm-1
   1.208990e-03_EB, 6.799487e-03_EB, 3.701547e-01_EB, 1.237888e-05_EB, 1.464532e-02_EB, 2.595773e-05_EB, &
   4.545433e-05_EB, &
   5.930652e-04_EB, 2.514086e-05_EB, 6.469106e-03_EB, 2.000295e-03_EB, 3.669583e-03_EB, 3.705261e-05_EB, &
   6.244904e-03_EB, &
   4.994954e-04_EB, 9.490064e-04_EB, 4.754384e-04_EB, 5.579883e-04_EB, 9.783514e-04_EB, 1.672579e-03_EB, &
   1.298467e-04_EB, &
   4.041578e-04_EB, 4.765211e-04_EB, 4.390729e-04_EB, 4.199440e-04_EB, 4.075633e-04_EB, 4.679761e-04_EB, &
   4.596989e-04_EB, &
   3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, &
   3.684034e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_c3h8(n_temp_c3h8,34)) 

! band #2: 2550 cm-1 - 3375 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.8902 % 

sd2_c3h8(1:7,1:8) = RESHAPE((/ &  ! 2550-2725 cm-1
   5.503672e-04_EB, 5.994891e-04_EB, 8.662153e-04_EB, 6.007139e-04_EB, 5.508955e-04_EB, 6.062568e-04_EB, &
   9.609866e-04_EB, &
   9.165832e-03_EB, 6.514927e-03_EB, 5.498817e-03_EB, 2.820723e-03_EB, 4.346969e-03_EB, 3.366148e-03_EB, &
   1.036676e-01_EB, &
   3.816653e-02_EB, 2.423984e-02_EB, 2.297132e-02_EB, 1.131018e-02_EB, 1.290251e-02_EB, 9.585606e-03_EB, &
   7.230381e-01_EB, &
   6.229937e-02_EB, 3.989116e-02_EB, 3.343177e-02_EB, 2.820251e-02_EB, 1.841367e-02_EB, 1.058259e-02_EB, &
   1.643570e-02_EB, &
   8.226034e-02_EB, 4.582033e-02_EB, 3.980542e-02_EB, 2.801891e-02_EB, 1.700839e-02_EB, 1.139598e-02_EB, &
   2.021613e-01_EB, &
   7.640614e-02_EB, 3.373610e-02_EB, 3.023935e-02_EB, 2.980953e-02_EB, 2.095503e-02_EB, 1.799163e-02_EB, &
   3.412770e-01_EB, &
   3.340838e-02_EB, 2.504949e-02_EB, 2.777922e-02_EB, 2.817541e-02_EB, 3.213869e-02_EB, 3.140060e-02_EB, &
   1.470802e+00_EB, &
   8.685682e-02_EB, 6.091073e-02_EB, 7.104290e-02_EB, 7.142576e-02_EB, 6.557492e-02_EB, 6.248472e-02_EB, &
   2.113951e+00_EB/),(/7,8/))

sd2_c3h8(1:7,9:16) = RESHAPE((/ &  ! 2750-2925 cm-1
   1.699556e-01_EB, 1.209358e-01_EB, 1.186218e-01_EB, 1.132869e-01_EB, 1.150490e-01_EB, 1.001152e-01_EB, &
   1.402520e-01_EB, &
   1.666365e-01_EB, 1.410150e-01_EB, 1.367199e-01_EB, 1.403357e-01_EB, 1.627667e-01_EB, 1.378065e-01_EB, &
   3.019538e-01_EB, &
   2.031290e-01_EB, 1.630134e-01_EB, 1.593693e-01_EB, 1.694187e-01_EB, 1.815236e-01_EB, 2.026379e-01_EB, &
   2.932861e-01_EB, &
   3.542177e-01_EB, 3.249511e-01_EB, 3.280019e-01_EB, 3.450738e-01_EB, 3.407162e-01_EB, 3.576253e-01_EB, &
   5.198445e-01_EB, &
   1.000630e+00_EB, 9.246602e-01_EB, 9.305942e-01_EB, 9.189938e-01_EB, 8.909863e-01_EB, 7.810609e-01_EB, &
   9.521582e-01_EB, &
   3.753627e+00_EB, 2.683351e+00_EB, 2.484390e+00_EB, 2.041213e+00_EB, 1.810422e+00_EB, 1.239899e+00_EB, &
   1.183325e+00_EB, &
   5.234499e+00_EB, 3.770775e+00_EB, 3.516223e+00_EB, 2.940932e+00_EB, 2.625379e+00_EB, 1.827482e+00_EB, &
   1.807152e+00_EB, &
   5.572566e+00_EB, 4.685164e+00_EB, 4.488098e+00_EB, 3.945955e+00_EB, 3.562082e+00_EB, 2.450846e+00_EB, &
   2.095864e+00_EB/),(/7,8/))

sd2_c3h8(1:7,17:24) = RESHAPE((/ &  ! 2950-3125 cm-1
   1.055304e+01_EB, 7.617962e+00_EB, 7.020237e+00_EB, 5.655513e+00_EB, 4.924307e+00_EB, 2.987908e+00_EB, &
   2.490259e+00_EB, &
   1.309229e+01_EB, 8.614162e+00_EB, 7.735378e+00_EB, 5.958736e+00_EB, 5.015764e+00_EB, 2.751758e+00_EB, &
   2.176924e+00_EB, &
   5.664769e+00_EB, 4.285205e+00_EB, 3.984661e+00_EB, 3.391957e+00_EB, 2.931808e+00_EB, 1.747475e+00_EB, &
   1.598403e+00_EB, &
   7.428057e-01_EB, 8.846817e-01_EB, 9.107382e-01_EB, 1.009605e+00_EB, 9.708428e-01_EB, 7.720345e-01_EB, &
   9.265547e-01_EB, &
   7.186847e-02_EB, 1.458443e-01_EB, 1.735809e-01_EB, 2.489900e-01_EB, 2.701093e-01_EB, 3.154699e-01_EB, &
   6.322075e-01_EB, &
   4.223055e-02_EB, 7.998122e-02_EB, 1.014124e-01_EB, 1.476562e-01_EB, 1.444812e-01_EB, 1.778681e-01_EB, &
   7.948678e-01_EB, &
   5.301805e-02_EB, 7.789503e-02_EB, 9.507427e-02_EB, 1.304625e-01_EB, 1.288331e-01_EB, 1.398066e-01_EB, &
   9.260605e-01_EB, &
   8.229182e-02_EB, 7.916035e-02_EB, 9.156865e-02_EB, 1.229671e-01_EB, 1.154257e-01_EB, 1.069746e-01_EB, &
   1.474222e+00_EB/),(/7,8/))

sd2_c3h8(1:7,25:32) = RESHAPE((/ &  ! 3150-3325 cm-1
   1.097023e-01_EB, 9.711792e-02_EB, 1.018166e-01_EB, 1.145716e-01_EB, 1.629561e-01_EB, 7.870503e-02_EB, &
   1.054101e+00_EB, &
   1.528366e-01_EB, 1.042414e-01_EB, 1.045964e-01_EB, 9.557668e-02_EB, 9.257045e-02_EB, 5.716098e-02_EB, &
   2.230393e+00_EB, &
   1.448776e-01_EB, 8.693362e-02_EB, 8.143001e-02_EB, 8.822165e-02_EB, 7.154741e-02_EB, 3.844782e-02_EB, &
   5.459571e+00_EB, &
   7.739671e-02_EB, 4.498684e-02_EB, 5.166537e-02_EB, 4.834891e-02_EB, 4.653514e-02_EB, 1.900636e-02_EB, &
   3.668670e+00_EB, &
   2.193118e-02_EB, 1.636816e-02_EB, 1.616765e-02_EB, 2.883862e-02_EB, 6.071064e-02_EB, 3.305292e-03_EB, &
   4.440979e+00_EB, &
   7.607103e-04_EB, 2.864772e-03_EB, 1.569824e-03_EB, 5.600370e-03_EB, 2.231541e-03_EB, 1.431965e-03_EB, &
   4.371189e+00_EB, &
   6.767298e-04_EB, 7.265903e-04_EB, 6.896754e-04_EB, 7.704584e-04_EB, 7.843659e-04_EB, 2.492767e-03_EB, &
   3.017731e+00_EB, &
   6.753685e-04_EB, 8.082467e-04_EB, 8.419705e-04_EB, 8.502049e-04_EB, 7.333137e-04_EB, 2.689703e-03_EB, &
   2.546579e-01_EB/),(/7,8/))

sd2_c3h8(1:7,33:34) = RESHAPE((/ &  ! 3350-3375 cm-1
   5.953124e-04_EB, 6.647555e-04_EB, 6.295086e-04_EB, 6.270639e-04_EB, 5.605951e-04_EB, 1.004627e-03_EB, &
   5.997674e-04_EB, &
   4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, 4.634348e-04_EB, &
   4.634348e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_c3h8(n_temp_c3h8,34)) 

! band #2: 2550 cm-1 - 3375 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.8902 % 

! print fine structure array gamma_d 

gammad2_c3h8(1:7,1:8) = RESHAPE((/ &  ! 2550-2725 cm-1
   4.136492e-04_EB, 4.360494e-04_EB, 7.953702e-04_EB, 4.660382e-04_EB, 4.005291e-04_EB, 4.511091e-04_EB, &
   6.072020e-04_EB, &
   2.759961e-03_EB, 1.062067e-01_EB, 1.215285e-02_EB, 2.640006e-03_EB, 3.656880e-03_EB, 1.938678e-03_EB, &
   5.487963e-05_EB, &
   1.584875e-02_EB, 2.833951e-01_EB, 3.264932e-01_EB, 1.739356e-02_EB, 5.773845e-03_EB, 4.395147e-02_EB, &
   1.322632e-05_EB, &
   1.193671e-01_EB, 2.854658e-01_EB, 4.644661e-02_EB, 2.765612e-02_EB, 1.053011e-02_EB, 2.667963e-01_EB, &
   5.817185e-03_EB, &
   2.623899e-02_EB, 1.875377e-01_EB, 5.240968e-02_EB, 2.355271e-02_EB, 7.559665e-02_EB, 2.293658e-01_EB, &
   7.216244e-05_EB, &
   5.846743e-03_EB, 2.105678e-01_EB, 2.701368e-01_EB, 2.713083e-02_EB, 1.035818e-02_EB, 3.749967e-02_EB, &
   8.197520e-05_EB, &
   1.507777e-02_EB, 2.373195e-01_EB, 2.803885e-01_EB, 2.381271e-02_EB, 2.374969e-03_EB, 3.007283e-01_EB, &
   3.941570e-05_EB, &
   5.460024e-02_EB, 2.148020e-01_EB, 1.102276e-02_EB, 2.048424e-02_EB, 1.381473e-02_EB, 2.265009e-01_EB, &
   1.113763e-04_EB/),(/7,8/))

gammad2_c3h8(1:7,9:16) = RESHAPE((/ &  ! 2750-2925 cm-1
   7.609446e-02_EB, 6.650037e-01_EB, 3.402647e-01_EB, 7.088071e-02_EB, 1.528056e-02_EB, 1.623423e+00_EB, &
   1.555020e-02_EB, &
   1.068152e-01_EB, 3.599354e-01_EB, 9.327433e-01_EB, 8.061623e-02_EB, 1.044707e-02_EB, 1.207247e+00_EB, &
   6.148015e-03_EB, &
   1.347777e-01_EB, 1.366370e-01_EB, 1.464407e+00_EB, 8.162167e-02_EB, 3.748354e-02_EB, 6.697765e+00_EB, &
   5.715348e-02_EB, &
   1.784972e-01_EB, 1.113695e-01_EB, 1.456847e-01_EB, 8.124496e-02_EB, 8.031435e-02_EB, 3.755111e+01_EB, &
   4.466462e-02_EB, &
   3.769734e-01_EB, 1.378926e-01_EB, 2.066963e-01_EB, 1.647811e-01_EB, 1.988552e-01_EB, 4.888631e+01_EB, &
   6.856421e-02_EB, &
   1.209828e+00_EB, 4.483110e-01_EB, 4.819971e-01_EB, 3.919982e-01_EB, 3.657716e-01_EB, 4.888633e+01_EB, &
   1.716140e-01_EB, &
   1.907010e+00_EB, 6.962287e-01_EB, 6.668651e-01_EB, 5.363387e-01_EB, 5.417758e-01_EB, 4.888614e+01_EB, &
   1.732383e-01_EB, &
   2.035309e+00_EB, 8.005331e-01_EB, 8.689003e-01_EB, 6.922331e-01_EB, 7.737400e-01_EB, 4.888562e+01_EB, &
   3.566562e-01_EB/),(/7,8/))

gammad2_c3h8(1:7,17:24) = RESHAPE((/ &  ! 2950-3125 cm-1
   3.521357e+00_EB, 1.169922e+00_EB, 1.245608e+00_EB, 9.879083e-01_EB, 9.805019e-01_EB, 4.888630e+01_EB, &
   2.805866e-01_EB, &
   8.194584e+00_EB, 1.715767e+00_EB, 1.407065e+00_EB, 1.043498e+00_EB, 9.327292e-01_EB, 4.888633e+01_EB, &
   2.729909e-01_EB, &
   9.096142e-01_EB, 9.519580e-01_EB, 6.722470e-01_EB, 5.099619e-01_EB, 4.841140e-01_EB, 4.888628e+01_EB, &
   2.097093e-01_EB, &
   1.758692e-01_EB, 8.005995e-01_EB, 3.028875e-01_EB, 1.771971e-01_EB, 1.723475e-01_EB, 4.888622e+01_EB, &
   9.485115e-02_EB, &
   1.552016e-02_EB, 1.106460e+00_EB, 3.508781e-01_EB, 8.578778e-02_EB, 7.481873e-02_EB, 2.362585e+01_EB, &
   2.190639e-02_EB, &
   8.141419e-03_EB, 2.936549e-01_EB, 1.051647e-01_EB, 3.583341e-02_EB, 8.451454e-02_EB, 4.015433e+00_EB, &
   4.672611e-03_EB, &
   6.825195e-03_EB, 2.860617e-01_EB, 9.805913e-02_EB, 3.356674e-02_EB, 3.562043e-02_EB, 1.354576e+00_EB, &
   2.504862e-03_EB, &
   1.731937e-02_EB, 8.103158e-02_EB, 3.656235e-02_EB, 1.650823e-02_EB, 1.716846e-02_EB, 1.104930e+00_EB, &
   1.025944e-03_EB/),(/7,8/))

gammad2_c3h8(1:7,25:32) = RESHAPE((/ &  ! 3150-3325 cm-1
   2.761892e-02_EB, 5.927558e-02_EB, 6.812671e-02_EB, 1.536547e-02_EB, 2.866389e-03_EB, 2.743296e-01_EB, &
   8.973476e-04_EB, &
   4.898929e-02_EB, 2.217375e-01_EB, 1.147533e-01_EB, 5.759265e-02_EB, 1.511087e-02_EB, 9.820712e-02_EB, &
   2.297806e-04_EB, &
   1.075239e-02_EB, 9.660209e-02_EB, 7.236142e-02_EB, 1.745200e-02_EB, 1.115524e-02_EB, 1.145421e-01_EB, &
   4.285123e-05_EB, &
   1.476159e-02_EB, 2.024331e-01_EB, 7.181537e-03_EB, 4.721618e-02_EB, 3.344976e-03_EB, 1.461756e+00_EB, &
   3.907613e-05_EB, &
   7.912117e-03_EB, 2.082529e-02_EB, 6.432685e-03_EB, 1.913236e-02_EB, 7.216503e-05_EB, 1.281409e-03_EB, &
   1.568408e-05_EB, &
   5.304482e-04_EB, 1.965226e-03_EB, 8.885309e-04_EB, 1.080202e-02_EB, 6.214399e-04_EB, 2.290412e-04_EB, &
   1.111330e-05_EB, &
   4.741388e-04_EB, 5.040738e-04_EB, 4.522965e-04_EB, 5.090295e-04_EB, 4.634938e-04_EB, 3.121426e-03_EB, &
   1.265450e-05_EB, &
   4.771043e-04_EB, 5.489571e-04_EB, 5.406306e-04_EB, 6.198226e-04_EB, 4.746614e-04_EB, 1.079570e-03_EB, &
   2.377068e-05_EB/),(/7,8/))

gammad2_c3h8(1:7,33:34) = RESHAPE((/ &  ! 3350-3375 cm-1
   4.399730e-04_EB, 5.039732e-04_EB, 4.232400e-04_EB, 4.772033e-04_EB, 4.129961e-04_EB, 5.016760e-04_EB, &
   4.444099e-04_EB, &
   3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, 3.684034e-04_EB, &
   3.684034e-04_EB/),(/7,2/))

!-------------------------propylene data-------------------


! there are 3 bands for propylene

! band #1: 775 cm-1 - 1150 cm-1 
! band #2: 1225 cm-1 - 1975 cm-1 
! band #3: 2650 cm-1 - 3275 cm-1 

ALLOCATE(sd_c3h6_temp(n_temp_c3h6)) 

! initialize bands wavenumber bounds for propylene ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_c3h6(n_band_c3h6,3)) 
ALLOCATE(be_c3h6(n_band_c3h6)) 

om_bnd_c3h6 = RESHAPE((/ &
   775._EB, 1225._EB, 2650._EB, &
   1150._EB, 1975._EB, 3275._EB, &
   5._EB, 25._EB, 25._EB/),(/n_band_c3h6,3/)) 

sd_c3h6_temp = (/ &
   296._EB, 390._EB, 444._EB, 491._EB, 594._EB, 793._EB,&
   1003._EB/)

be_c3h6 = (/ &
   1.000_EB, 1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_c3h6(n_temp_c3h6,68)) 

! band #1: 775 cm-1 - 1150 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.3448 % 

sd1_c3h6(1:7,1:8) = RESHAPE((/ &  ! 775-810 cm-1
   2.921820e-04_EB, 1.885331e-04_EB, 1.578479e-04_EB, 1.195361e-04_EB, 1.447279e-04_EB, 1.698242e-04_EB, &
   3.510553e-04_EB, &
   5.786383e-04_EB, 4.698237e-04_EB, 3.471753e-04_EB, 3.368251e-04_EB, 1.561044e-04_EB, 2.151079e-04_EB, &
   4.880960e-04_EB, &
   5.146434e-04_EB, 4.282047e-04_EB, 3.953252e-04_EB, 3.606827e-04_EB, 3.453987e-04_EB, 3.100075e-04_EB, &
   5.373971e-04_EB, &
   1.002033e-03_EB, 3.754947e-04_EB, 3.830171e-04_EB, 3.397481e-04_EB, 4.382293e-03_EB, 1.980601e-04_EB, &
   4.996406e-04_EB, &
   5.974180e-04_EB, 4.192863e-04_EB, 2.153249e-04_EB, 1.486372e-03_EB, 3.085406e-02_EB, 4.786688e-03_EB, &
   4.055803e-04_EB, &
   3.186779e-02_EB, 1.217059e-02_EB, 8.908478e-03_EB, 1.366154e-02_EB, 1.604611e-01_EB, 1.277554e-02_EB, &
   3.868758e-04_EB, &
   2.287860e-02_EB, 1.347562e-02_EB, 1.761538e-02_EB, 1.670141e-02_EB, 3.006989e-01_EB, 1.352943e-02_EB, &
   4.562671e-04_EB, &
   3.847401e-01_EB, 1.311464e-01_EB, 1.806964e-02_EB, 2.138006e-02_EB, 6.075345e-01_EB, 2.375582e-02_EB, &
   4.382283e-04_EB/),(/7,8/))

sd1_c3h6(1:7,9:16) = RESHAPE((/ &  ! 815-850 cm-1
   3.612608e-01_EB, 3.310333e-02_EB, 2.135533e-02_EB, 6.703345e-02_EB, 3.383135e-01_EB, 3.802107e-02_EB, &
   1.684642e-03_EB, &
   4.242283e-01_EB, 3.658435e-02_EB, 3.805267e-02_EB, 5.077876e-02_EB, 3.413346e-01_EB, 6.513409e-02_EB, &
   1.587704e-04_EB, &
   1.644787e-01_EB, 5.568292e-02_EB, 5.486572e-02_EB, 6.891789e-02_EB, 1.332364e-01_EB, 6.265865e-02_EB, &
   4.187924e-03_EB, &
   7.882443e-02_EB, 7.779061e-02_EB, 7.929456e-02_EB, 1.096667e-01_EB, 1.347277e-01_EB, 9.669708e-02_EB, &
   2.546892e-02_EB, &
   9.754850e-02_EB, 1.076154e-01_EB, 1.043330e-01_EB, 1.111524e-01_EB, 1.883759e-01_EB, 1.329948e-01_EB, &
   5.845340e-02_EB, &
   1.538045e-01_EB, 1.691789e-01_EB, 1.519359e-01_EB, 1.496250e-01_EB, 2.160998e-01_EB, 1.631159e-01_EB, &
   9.210403e-02_EB, &
   2.041607e-01_EB, 2.369334e-01_EB, 1.944714e-01_EB, 2.025469e-01_EB, 2.541597e-01_EB, 1.848773e-01_EB, &
   1.157566e-01_EB, &
   2.823440e-01_EB, 3.067145e-01_EB, 2.545703e-01_EB, 3.082214e-01_EB, 3.115871e-01_EB, 2.118905e-01_EB, &
   1.367771e-01_EB/),(/7,8/))

sd1_c3h6(1:7,17:24) = RESHAPE((/ &  ! 855-890 cm-1
   3.934069e-01_EB, 3.902770e-01_EB, 3.324288e-01_EB, 3.655227e-01_EB, 4.852182e-01_EB, 2.278497e-01_EB, &
   1.418258e-01_EB, &
   5.098671e-01_EB, 4.707600e-01_EB, 3.953366e-01_EB, 4.123636e-01_EB, 5.023609e-01_EB, 3.093041e-01_EB, &
   2.172109e-01_EB, &
   6.861408e-01_EB, 5.874598e-01_EB, 4.823945e-01_EB, 4.930154e-01_EB, 5.112207e-01_EB, 3.097694e-01_EB, &
   1.945893e-01_EB, &
   8.621571e-01_EB, 7.069132e-01_EB, 5.626431e-01_EB, 5.778790e-01_EB, 6.253640e-01_EB, 4.501915e-01_EB, &
   1.364634e-01_EB, &
   1.117672e+00_EB, 8.597340e-01_EB, 6.702957e-01_EB, 6.824590e-01_EB, 6.721777e-01_EB, 4.204596e-01_EB, &
   2.336953e-01_EB, &
   1.296785e+00_EB, 9.660400e-01_EB, 7.508644e-01_EB, 7.657156e-01_EB, 6.694906e-01_EB, 4.150445e-01_EB, &
   2.912663e-01_EB, &
   1.612024e+00_EB, 1.132199e+00_EB, 8.825466e-01_EB, 8.439805e-01_EB, 7.306116e-01_EB, 4.492608e-01_EB, &
   2.586394e-01_EB, &
   1.769679e+00_EB, 1.240089e+00_EB, 9.826761e-01_EB, 9.126471e-01_EB, 7.901846e-01_EB, 4.958428e-01_EB, &
   2.975731e-01_EB/),(/7,8/))

sd1_c3h6(1:7,25:32) = RESHAPE((/ &  ! 895-930 cm-1
   2.035152e+00_EB, 1.372164e+00_EB, 1.084943e+00_EB, 9.759957e-01_EB, 8.670922e-01_EB, 5.511628e-01_EB, &
   5.051868e-01_EB, &
   2.110903e+00_EB, 1.411810e+00_EB, 1.107883e+00_EB, 1.032360e+00_EB, 8.891386e-01_EB, 6.607676e-01_EB, &
   4.105057e-01_EB, &
   2.241734e+00_EB, 1.509671e+00_EB, 1.226689e+00_EB, 1.120333e+00_EB, 1.019009e+00_EB, 6.960237e-01_EB, &
   6.300682e-01_EB, &
   4.690927e+00_EB, 3.301112e+00_EB, 2.617543e+00_EB, 2.271022e+00_EB, 1.868783e+00_EB, 1.104123e+00_EB, &
   6.791524e-01_EB, &
   3.538563e+00_EB, 2.596313e+00_EB, 2.151089e+00_EB, 1.982625e+00_EB, 1.674390e+00_EB, 1.121363e+00_EB, &
   6.319589e-01_EB, &
   2.529174e+00_EB, 1.690171e+00_EB, 1.371477e+00_EB, 1.268574e+00_EB, 1.124475e+00_EB, 7.538665e-01_EB, &
   4.997968e-01_EB, &
   2.442382e+00_EB, 1.621890e+00_EB, 1.292151e+00_EB, 1.185851e+00_EB, 1.045526e+00_EB, 6.515910e-01_EB, &
   4.350142e-01_EB, &
   2.359003e+00_EB, 1.574882e+00_EB, 1.255634e+00_EB, 1.148501e+00_EB, 1.005751e+00_EB, 6.742696e-01_EB, &
   4.032512e-01_EB/),(/7,8/))

sd1_c3h6(1:7,33:40) = RESHAPE((/ &  ! 935-970 cm-1
   2.292826e+00_EB, 1.573504e+00_EB, 1.257992e+00_EB, 1.148404e+00_EB, 9.819016e-01_EB, 6.786328e-01_EB, &
   3.987199e-01_EB, &
   2.256277e+00_EB, 1.614009e+00_EB, 1.308043e+00_EB, 1.170636e+00_EB, 1.005170e+00_EB, 6.614059e-01_EB, &
   3.824164e-01_EB, &
   2.128388e+00_EB, 1.579098e+00_EB, 1.312822e+00_EB, 1.165846e+00_EB, 1.012764e+00_EB, 6.750425e-01_EB, &
   4.478263e-01_EB, &
   2.002676e+00_EB, 1.570187e+00_EB, 1.337531e+00_EB, 1.201187e+00_EB, 1.060851e+00_EB, 6.651784e-01_EB, &
   5.518473e-01_EB, &
   1.769536e+00_EB, 1.476601e+00_EB, 1.281675e+00_EB, 1.168120e+00_EB, 1.059359e+00_EB, 6.557750e-01_EB, &
   4.612030e-01_EB, &
   1.513007e+00_EB, 1.314704e+00_EB, 1.154175e+00_EB, 1.078702e+00_EB, 9.708818e-01_EB, 7.985200e-01_EB, &
   4.095132e-01_EB, &
   1.300819e+00_EB, 1.112939e+00_EB, 9.729632e-01_EB, 9.305722e-01_EB, 9.107075e-01_EB, 6.817647e-01_EB, &
   3.869006e-01_EB, &
   1.122524e+00_EB, 9.240820e-01_EB, 7.997622e-01_EB, 8.045075e-01_EB, 8.225817e-01_EB, 5.935068e-01_EB, &
   4.076619e-01_EB/),(/7,8/))

sd1_c3h6(1:7,41:48) = RESHAPE((/ &  ! 975-1010 cm-1
   1.090787e+00_EB, 8.705917e-01_EB, 7.392297e-01_EB, 7.356904e-01_EB, 7.567491e-01_EB, 5.944416e-01_EB, &
   3.882052e-01_EB, &
   1.097680e+00_EB, 8.498178e-01_EB, 7.146480e-01_EB, 7.049228e-01_EB, 7.300468e-01_EB, 5.868185e-01_EB, &
   4.235053e-01_EB, &
   1.175184e+00_EB, 9.551431e-01_EB, 8.117959e-01_EB, 7.958744e-01_EB, 7.668238e-01_EB, 5.942202e-01_EB, &
   4.784101e-01_EB, &
   2.200917e+00_EB, 1.489225e+00_EB, 1.189617e+00_EB, 1.045546e+00_EB, 8.888219e-01_EB, 6.345618e-01_EB, &
   3.753116e-01_EB, &
   1.026273e+00_EB, 8.138120e-01_EB, 7.081956e-01_EB, 6.489008e-01_EB, 5.962778e-01_EB, 4.346981e-01_EB, &
   3.144834e-01_EB, &
   8.545948e-01_EB, 6.501326e-01_EB, 5.807415e-01_EB, 5.228005e-01_EB, 5.069268e-01_EB, 3.587344e-01_EB, &
   2.535124e-01_EB, &
   8.182404e-01_EB, 6.020895e-01_EB, 5.272828e-01_EB, 4.843306e-01_EB, 4.757196e-01_EB, 3.642910e-01_EB, &
   2.463011e-01_EB, &
   7.867336e-01_EB, 5.627888e-01_EB, 4.859257e-01_EB, 4.529062e-01_EB, 4.373186e-01_EB, 2.999539e-01_EB, &
   2.195520e-01_EB/),(/7,8/))

sd1_c3h6(1:7,49:56) = RESHAPE((/ &  ! 1015-1050 cm-1
   7.370735e-01_EB, 5.261814e-01_EB, 4.365308e-01_EB, 4.196529e-01_EB, 3.794978e-01_EB, 3.100276e-01_EB, &
   2.316864e-01_EB, &
   6.545905e-01_EB, 4.754812e-01_EB, 3.879959e-01_EB, 3.680691e-01_EB, 3.486950e-01_EB, 2.537942e-01_EB, &
   2.192666e-01_EB, &
   5.636151e-01_EB, 4.220184e-01_EB, 3.368352e-01_EB, 3.294973e-01_EB, 3.480874e-01_EB, 2.278502e-01_EB, &
   1.435619e-01_EB, &
   4.794628e-01_EB, 3.716607e-01_EB, 2.986875e-01_EB, 2.920226e-01_EB, 2.849022e-01_EB, 2.098157e-01_EB, &
   1.086582e-01_EB, &
   3.987250e-01_EB, 3.273405e-01_EB, 2.706509e-01_EB, 2.488392e-01_EB, 2.796603e-01_EB, 1.937969e-01_EB, &
   1.318512e-01_EB, &
   3.618884e-01_EB, 3.043278e-01_EB, 2.752540e-01_EB, 2.348242e-01_EB, 2.449431e-01_EB, 1.840862e-01_EB, &
   1.365370e-01_EB, &
   4.327196e-01_EB, 3.118666e-01_EB, 2.824791e-01_EB, 2.285493e-01_EB, 2.297279e-01_EB, 1.660089e-01_EB, &
   6.776279e-02_EB, &
   2.158610e-01_EB, 1.711786e-01_EB, 1.890422e-01_EB, 1.508251e-01_EB, 1.710093e-01_EB, 1.189088e-01_EB, &
   3.722276e-01_EB/),(/7,8/))

sd1_c3h6(1:7,57:64) = RESHAPE((/ &  ! 1055-1090 cm-1
   1.808366e-01_EB, 1.285694e-01_EB, 1.347601e-01_EB, 1.159091e-01_EB, 1.328783e-01_EB, 9.806016e-02_EB, &
   4.385238e-02_EB, &
   1.688521e-01_EB, 1.098959e-01_EB, 9.863147e-02_EB, 1.007307e-01_EB, 1.218736e-01_EB, 1.096671e-01_EB, &
   3.048481e-02_EB, &
   1.509152e-01_EB, 9.784628e-02_EB, 8.122169e-02_EB, 7.571922e-02_EB, 8.082920e-02_EB, 6.716860e-02_EB, &
   6.805299e-01_EB, &
   1.236023e-01_EB, 8.420975e-02_EB, 6.449748e-02_EB, 5.721581e-02_EB, 8.699266e-02_EB, 4.794369e-02_EB, &
   5.275249e-02_EB, &
   1.005482e-01_EB, 7.333190e-02_EB, 5.246391e-02_EB, 4.854230e-02_EB, 8.581595e-02_EB, 3.857051e-02_EB, &
   1.047535e+00_EB, &
   8.140574e-02_EB, 6.139425e-02_EB, 4.191643e-02_EB, 3.330029e-02_EB, 6.128860e-02_EB, 3.771955e-02_EB, &
   5.547067e-01_EB, &
   5.876207e-02_EB, 4.748412e-02_EB, 3.403509e-02_EB, 2.491123e-02_EB, 1.741583e-01_EB, 3.685622e-02_EB, &
   1.173988e+00_EB, &
   5.067589e-02_EB, 3.249100e-02_EB, 3.614176e-02_EB, 1.695154e-02_EB, 5.607839e-02_EB, 4.126856e-02_EB, &
   4.414252e-02_EB/),(/7,8/))

sd1_c3h6(1:7,65:68) = RESHAPE((/ &  ! 1095-1150 cm-1
   3.475507e-02_EB, 1.803126e-02_EB, 2.233876e-01_EB, 1.137397e-02_EB, 1.233139e-02_EB, 2.244437e-02_EB, &
   3.055494e-04_EB, &
   8.633613e-02_EB, 9.384641e-03_EB, 3.751848e-01_EB, 5.532895e-03_EB, 6.183470e-03_EB, 7.005080e-03_EB, &
   5.030021e-04_EB, &
   4.004646e-04_EB, 1.996236e-04_EB, 3.332406e-04_EB, 1.767862e-04_EB, 1.467998e-04_EB, 1.391831e-04_EB, &
   3.497050e-04_EB, &
   1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, &
   1.027560e-04_EB/),(/7,4/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_c3h6(n_temp_c3h6,68)) 

! band #1: 775 cm-1 - 1150 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.3448 % 

! print fine structure array gamma_d 

gammad1_c3h6(1:7,1:8) = RESHAPE((/ &  ! 775-810 cm-1
   1.865203e-04_EB, 1.481280e-04_EB, 1.121426e-04_EB, 8.618315e-05_EB, 9.939553e-05_EB, 1.228308e-04_EB, &
   2.571832e-04_EB, &
   2.784672e-04_EB, 2.629293e-04_EB, 2.613280e-04_EB, 2.166322e-04_EB, 1.060275e-04_EB, 1.512432e-04_EB, &
   3.364746e-04_EB, &
   3.584490e-04_EB, 2.571836e-04_EB, 2.937264e-04_EB, 2.336679e-04_EB, 2.672244e-04_EB, 2.128761e-04_EB, &
   2.711168e-04_EB, &
   1.050069e-03_EB, 2.060006e-04_EB, 2.640878e-04_EB, 2.179849e-04_EB, 2.169826e-04_EB, 1.354118e-04_EB, &
   2.786454e-04_EB, &
   5.631752e-04_EB, 2.967683e-04_EB, 1.474684e-04_EB, 6.104105e-04_EB, 4.865310e-04_EB, 1.882053e-01_EB, &
   2.951885e-04_EB, &
   8.382309e-05_EB, 1.999192e-03_EB, 2.110707e-01_EB, 2.657106e-03_EB, 4.279171e-04_EB, 8.580741e-01_EB, &
   2.380492e-04_EB, &
   4.133686e-04_EB, 2.280390e-03_EB, 8.751086e-01_EB, 3.295157e-03_EB, 5.879263e-04_EB, 8.275870e-01_EB, &
   2.164440e-04_EB, &
   6.177217e-05_EB, 2.052275e-04_EB, 1.783767e+00_EB, 1.180004e-02_EB, 4.686960e-04_EB, 2.654668e+00_EB, &
   2.743441e-04_EB/),(/7,8/))

gammad1_c3h6(1:7,9:16) = RESHAPE((/ &  ! 815-850 cm-1
   1.362853e-04_EB, 9.486143e-03_EB, 2.766580e+00_EB, 2.523733e-03_EB, 1.516289e-03_EB, 8.166695e+00_EB, &
   1.693790e-01_EB, &
   2.370331e-04_EB, 1.279314e+00_EB, 1.307893e+01_EB, 1.368493e-01_EB, 2.662950e-03_EB, 9.028076e+00_EB, &
   1.501266e-04_EB, &
   1.605909e-03_EB, 3.888338e+00_EB, 3.291564e+01_EB, 2.430327e+00_EB, 1.840118e-02_EB, 4.830532e+00_EB, &
   1.460557e-01_EB, &
   1.751757e-02_EB, 7.388926e+00_EB, 6.531922e+01_EB, 9.190465e-02_EB, 7.887271e-02_EB, 7.702450e+00_EB, &
   3.065075e+00_EB, &
   4.005618e-02_EB, 5.793342e+00_EB, 6.531359e+01_EB, 3.481985e+01_EB, 7.477765e-02_EB, 1.306801e+01_EB, &
   3.031110e+01_EB, &
   4.670461e-02_EB, 1.511945e-01_EB, 6.531029e+01_EB, 5.058553e+01_EB, 1.162566e-01_EB, 1.593096e+01_EB, &
   4.345904e+01_EB, &
   6.963672e-02_EB, 1.180935e-01_EB, 6.530489e+01_EB, 8.198354e-01_EB, 1.298311e-01_EB, 4.887623e+01_EB, &
   4.345787e+01_EB, &
   1.032951e-01_EB, 1.651780e-01_EB, 6.527527e+01_EB, 1.583085e-01_EB, 1.439537e-01_EB, 4.887639e+01_EB, &
   4.345745e+01_EB/),(/7,8/))

gammad1_c3h6(1:7,17:24) = RESHAPE((/ &  ! 855-890 cm-1
   1.554114e-01_EB, 2.463349e-01_EB, 1.552737e+00_EB, 2.489268e-01_EB, 7.175075e-02_EB, 4.877618e+01_EB, &
   4.344856e+01_EB, &
   2.216480e-01_EB, 3.911085e-01_EB, 6.528132e+01_EB, 4.847826e-01_EB, 1.055938e-01_EB, 7.890568e-01_EB, &
   4.345834e+01_EB, &
   3.344937e-01_EB, 5.937512e-01_EB, 6.440865e+01_EB, 5.685496e-01_EB, 2.154384e-01_EB, 4.887601e+01_EB, &
   4.345865e+01_EB, &
   4.563788e-01_EB, 6.544545e-01_EB, 6.527427e+01_EB, 5.645359e-01_EB, 1.603139e-01_EB, 2.220849e-01_EB, &
   4.345957e+01_EB, &
   6.971662e-01_EB, 8.419195e-01_EB, 6.602209e+00_EB, 6.185587e-01_EB, 2.185881e-01_EB, 1.026546e+00_EB, &
   4.345911e+01_EB, &
   9.177302e-01_EB, 8.698788e-01_EB, 4.219718e+00_EB, 5.938638e-01_EB, 3.550496e-01_EB, 4.887400e+01_EB, &
   4.345951e+01_EB, &
   1.024886e+00_EB, 8.931182e-01_EB, 1.906388e+00_EB, 7.671322e-01_EB, 4.049625e-01_EB, 4.886302e+01_EB, &
   4.345956e+01_EB, &
   1.138590e+00_EB, 8.055441e-01_EB, 1.266173e+00_EB, 7.247405e-01_EB, 4.313090e-01_EB, 2.888925e+00_EB, &
   4.345947e+01_EB/),(/7,8/))

gammad1_c3h6(1:7,25:32) = RESHAPE((/ &  ! 895-930 cm-1
   1.258896e+00_EB, 8.879635e-01_EB, 1.183716e+00_EB, 8.697144e-01_EB, 3.994743e-01_EB, 1.020734e+00_EB, &
   1.469295e-01_EB, &
   1.489911e+00_EB, 9.449031e-01_EB, 1.364975e+00_EB, 7.600728e-01_EB, 4.690280e-01_EB, 5.532825e-01_EB, &
   2.320330e+01_EB, &
   1.791634e+00_EB, 1.068040e+00_EB, 1.344772e+00_EB, 1.025605e+00_EB, 5.570787e-01_EB, 4.887568e+01_EB, &
   1.866191e-01_EB, &
   6.582008e-01_EB, 6.151496e-01_EB, 7.320593e-01_EB, 7.674633e-01_EB, 5.984986e-01_EB, 1.256776e+00_EB, &
   1.505021e+00_EB, &
   8.496429e-01_EB, 6.647729e-01_EB, 8.138054e-01_EB, 6.996761e-01_EB, 6.263920e-01_EB, 7.887495e-01_EB, &
   4.345957e+01_EB, &
   2.585850e+00_EB, 1.312862e+00_EB, 1.669178e+00_EB, 1.005278e+00_EB, 5.484994e-01_EB, 1.645085e+00_EB, &
   4.301288e+01_EB, &
   2.304432e+00_EB, 1.337031e+00_EB, 1.903523e+00_EB, 1.080752e+00_EB, 4.772959e-01_EB, 4.887638e+01_EB, &
   3.954397e+01_EB, &
   2.284405e+00_EB, 1.252628e+00_EB, 1.761317e+00_EB, 1.081162e+00_EB, 5.126407e-01_EB, 9.460859e-01_EB, &
   4.345934e+01_EB/),(/7,8/))

gammad1_c3h6(1:7,33:40) = RESHAPE((/ &  ! 935-970 cm-1
   2.011991e+00_EB, 1.160491e+00_EB, 1.606206e+00_EB, 1.044964e+00_EB, 6.166918e-01_EB, 8.541470e-01_EB, &
   4.345947e+01_EB, &
   1.951894e+00_EB, 1.077050e+00_EB, 1.340221e+00_EB, 1.005762e+00_EB, 5.377001e-01_EB, 1.813695e+00_EB, &
   4.345915e+01_EB, &
   1.862497e+00_EB, 1.110545e+00_EB, 1.226215e+00_EB, 1.004332e+00_EB, 5.640741e-01_EB, 1.260546e+00_EB, &
   2.061631e+01_EB, &
   1.735085e+00_EB, 1.152792e+00_EB, 1.235271e+00_EB, 1.029092e+00_EB, 5.304049e-01_EB, 3.364825e+00_EB, &
   4.109795e-01_EB, &
   1.514333e+00_EB, 1.096526e+00_EB, 1.190673e+00_EB, 9.738878e-01_EB, 4.991815e-01_EB, 2.598609e+01_EB, &
   4.345955e+01_EB, &
   1.323596e+00_EB, 1.020536e+00_EB, 1.276590e+00_EB, 9.362680e-01_EB, 5.648663e-01_EB, 3.325554e-01_EB, &
   4.345900e+01_EB, &
   1.090160e+00_EB, 8.522171e-01_EB, 1.383531e+00_EB, 8.696999e-01_EB, 3.985636e-01_EB, 6.142456e-01_EB, &
   4.345952e+01_EB, &
   9.568718e-01_EB, 7.965634e-01_EB, 2.061656e+00_EB, 6.937459e-01_EB, 3.185588e-01_EB, 1.043529e+00_EB, &
   4.345956e+01_EB/),(/7,8/))

gammad1_c3h6(1:7,41:48) = RESHAPE((/ &  ! 975-1010 cm-1
   8.780028e-01_EB, 7.171719e-01_EB, 1.797413e+00_EB, 7.191605e-01_EB, 3.123087e-01_EB, 5.370041e-01_EB, &
   6.031057e-01_EB, &
   7.683322e-01_EB, 6.697888e-01_EB, 1.564974e+00_EB, 6.922381e-01_EB, 2.866901e-01_EB, 6.100381e-01_EB, &
   5.296325e-01_EB, &
   7.710645e-01_EB, 6.686977e-01_EB, 1.316990e+00_EB, 7.356423e-01_EB, 3.962548e-01_EB, 7.103715e-01_EB, &
   3.121685e-01_EB, &
   7.152139e-01_EB, 7.181107e-01_EB, 9.272079e-01_EB, 7.998523e-01_EB, 4.476511e-01_EB, 4.403214e-01_EB, &
   1.432304e+00_EB, &
   7.215480e-01_EB, 5.222809e-01_EB, 7.034660e-01_EB, 5.939143e-01_EB, 3.557110e-01_EB, 7.369583e+00_EB, &
   4.343647e+01_EB, &
   7.509675e-01_EB, 5.541920e-01_EB, 6.013325e-01_EB, 5.827632e-01_EB, 2.516345e-01_EB, 4.887235e+01_EB, &
   4.345657e+01_EB, &
   6.655689e-01_EB, 5.909646e-01_EB, 7.175681e-01_EB, 5.191324e-01_EB, 2.139999e-01_EB, 4.270810e-01_EB, &
   3.186474e+00_EB, &
   5.936443e-01_EB, 6.111940e-01_EB, 8.318915e-01_EB, 4.526249e-01_EB, 2.171504e-01_EB, 4.874889e+01_EB, &
   4.345936e+01_EB/),(/7,8/))

gammad1_c3h6(1:7,49:56) = RESHAPE((/ &  ! 1015-1050 cm-1
   5.108597e-01_EB, 5.449935e-01_EB, 1.467166e+00_EB, 4.175570e-01_EB, 2.800624e-01_EB, 2.807650e-01_EB, &
   1.226529e-01_EB, &
   4.396741e-01_EB, 4.212443e-01_EB, 1.991238e+00_EB, 4.853656e-01_EB, 1.994488e-01_EB, 3.606763e+00_EB, &
   1.195304e-01_EB, &
   3.521564e-01_EB, 3.419630e-01_EB, 9.696171e+00_EB, 3.944128e-01_EB, 1.086559e-01_EB, 4.887097e+01_EB, &
   4.345940e+01_EB, &
   2.826655e-01_EB, 2.596188e-01_EB, 1.720670e+00_EB, 3.339207e-01_EB, 1.507599e-01_EB, 1.947669e+00_EB, &
   4.345956e+01_EB, &
   2.103926e-01_EB, 1.918390e-01_EB, 5.187904e-01_EB, 5.571757e-01_EB, 8.852491e-02_EB, 1.040530e+00_EB, &
   4.345956e+01_EB, &
   2.190684e-01_EB, 2.166025e-01_EB, 2.770976e-01_EB, 6.790640e-01_EB, 1.398958e-01_EB, 4.887609e+01_EB, &
   1.273985e+01_EB, &
   2.372381e-01_EB, 2.099844e-01_EB, 1.691730e-01_EB, 3.370993e-01_EB, 9.039351e-02_EB, 5.560844e-01_EB, &
   3.729305e+01_EB, &
   1.338894e-01_EB, 2.270514e-01_EB, 8.598819e-02_EB, 2.237195e-01_EB, 4.915478e-02_EB, 3.404705e+01_EB, &
   2.634276e-03_EB/),(/7,8/))

gammad1_c3h6(1:7,57:64) = RESHAPE((/ &  ! 1055-1090 cm-1
   8.590987e-02_EB, 2.939130e-01_EB, 1.322061e-01_EB, 2.325564e-01_EB, 4.121415e-02_EB, 1.890644e+01_EB, &
   6.653137e+00_EB, &
   5.828702e-02_EB, 1.919845e-01_EB, 3.558396e+00_EB, 8.935135e-02_EB, 2.151974e-02_EB, 1.813839e-02_EB, &
   5.291639e+00_EB, &
   4.560095e-02_EB, 8.801874e-02_EB, 4.240838e+00_EB, 1.971226e-01_EB, 2.805033e-02_EB, 3.159379e-01_EB, &
   4.231868e-04_EB, &
   5.236543e-02_EB, 4.750499e-02_EB, 7.719299e+00_EB, 2.230632e+00_EB, 1.029862e-02_EB, 3.970081e+00_EB, &
   1.254647e-02_EB, &
   4.584112e-02_EB, 2.781221e-02_EB, 3.521327e+00_EB, 1.069332e-01_EB, 5.513817e-03_EB, 7.390687e-01_EB, &
   8.785510e-05_EB, &
   2.952585e-02_EB, 1.659069e-02_EB, 1.232312e+00_EB, 1.610279e+00_EB, 3.746132e-03_EB, 1.888989e-02_EB, &
   1.172808e-05_EB, &
   2.783898e-02_EB, 1.191524e-02_EB, 1.737513e-01_EB, 6.337914e-01_EB, 3.927326e-04_EB, 1.256543e-02_EB, &
   9.593758e-05_EB, &
   1.161773e-02_EB, 1.234249e-02_EB, 7.173288e-03_EB, 4.651250e-01_EB, 6.905694e-04_EB, 6.625808e-03_EB, &
   1.875822e-05_EB/),(/7,8/))

gammad1_c3h6(1:7,65:68) = RESHAPE((/ &  ! 1095-1150 cm-1
   7.371375e-03_EB, 3.892307e-02_EB, 2.066948e-04_EB, 2.148055e-01_EB, 5.123043e-03_EB, 8.826795e-03_EB, &
   1.519386e-04_EB, &
   1.685293e-04_EB, 4.052437e-02_EB, 5.766893e-05_EB, 1.099394e-01_EB, 1.056916e-02_EB, 3.703491e-03_EB, &
   4.534644e-04_EB, &
   2.576720e-04_EB, 1.502775e-04_EB, 2.394291e-04_EB, 1.106949e-04_EB, 1.070129e-04_EB, 1.035473e-04_EB, &
   2.424308e-04_EB, &
   8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, &
   8.047518e-05_EB/),(/7,4/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_c3h6(n_temp_c3h6,31)) 

! band #2: 1225 cm-1 - 1975 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.89939 % 

sd2_c3h6(1:7,1:8) = RESHAPE((/ &  ! 1225-1400 cm-1
   3.558545e-04_EB, 3.234771e-04_EB, 3.304950e-04_EB, 1.651299e-04_EB, 1.336215e-04_EB, 1.510735e-04_EB, &
   1.993574e-04_EB, &
   6.741050e-04_EB, 3.773009e-04_EB, 4.527377e-04_EB, 3.445857e-04_EB, 2.150432e-04_EB, 3.112911e-04_EB, &
   3.405452e-03_EB, &
   8.500802e-04_EB, 5.928536e-03_EB, 5.109505e-02_EB, 5.939970e-03_EB, 5.544562e-03_EB, 4.216340e-02_EB, &
   7.821472e-02_EB, &
   6.900284e-04_EB, 1.754929e-02_EB, 1.351477e-01_EB, 1.875618e-02_EB, 1.689271e-02_EB, 5.437608e-01_EB, &
   8.042160e-02_EB, &
   6.597844e-03_EB, 3.076774e-02_EB, 7.267186e-02_EB, 3.875256e-02_EB, 4.198566e-02_EB, 1.377841e-01_EB, &
   4.353146e-02_EB, &
   1.425388e-01_EB, 7.292994e-02_EB, 1.006956e-01_EB, 8.530846e-02_EB, 8.207286e-02_EB, 1.024217e-01_EB, &
   9.040341e-02_EB, &
   2.757070e-01_EB, 2.011514e-01_EB, 2.183300e-01_EB, 1.901067e-01_EB, 1.708424e-01_EB, 1.802992e-01_EB, &
   1.306124e-01_EB, &
   5.487648e-01_EB, 3.944423e-01_EB, 3.612309e-01_EB, 3.179025e-01_EB, 2.635343e-01_EB, 2.207712e-01_EB, &
   1.537228e-01_EB/),(/7,8/))

sd2_c3h6(1:7,9:16) = RESHAPE((/ &  ! 1425-1600 cm-1
   7.347914e-01_EB, 4.838524e-01_EB, 4.349597e-01_EB, 3.807440e-01_EB, 3.206784e-01_EB, 2.496395e-01_EB, &
   1.756546e-01_EB, &
   1.230079e+00_EB, 7.176913e-01_EB, 5.767768e-01_EB, 4.861231e-01_EB, 3.530715e-01_EB, 2.411707e-01_EB, &
   1.663700e-01_EB, &
   1.076204e+00_EB, 6.242017e-01_EB, 5.169408e-01_EB, 4.238453e-01_EB, 3.142849e-01_EB, 2.371514e-01_EB, &
   1.547102e-01_EB, &
   4.989020e-01_EB, 3.259007e-01_EB, 2.647985e-01_EB, 2.506597e-01_EB, 2.113725e-01_EB, 1.667162e-01_EB, &
   1.268622e-01_EB, &
   3.050060e-01_EB, 1.549283e-01_EB, 1.437705e-01_EB, 1.418842e-01_EB, 1.284408e-01_EB, 1.165257e-01_EB, &
   1.181179e-01_EB, &
   2.901797e-01_EB, 1.013263e-01_EB, 8.839726e-02_EB, 9.945842e-02_EB, 9.936424e-02_EB, 9.382087e-02_EB, &
   1.109026e-01_EB, &
   9.170595e-01_EB, 1.077723e-01_EB, 1.126792e-01_EB, 1.012963e-01_EB, 1.008830e-01_EB, 1.343878e-01_EB, &
   1.433742e-01_EB, &
   1.723540e+00_EB, 1.340775e-01_EB, 1.803447e-01_EB, 1.684967e-01_EB, 1.917971e-01_EB, 2.426751e-01_EB, &
   1.993900e-01_EB/),(/7,8/))

sd2_c3h6(1:7,17:24) = RESHAPE((/ &  ! 1625-1800 cm-1
   7.534624e-01_EB, 5.257745e-01_EB, 4.873627e-01_EB, 4.420277e-01_EB, 3.752945e-01_EB, 2.872981e-01_EB, &
   2.158678e-01_EB, &
   1.173360e+00_EB, 6.714563e-01_EB, 5.383806e-01_EB, 4.764335e-01_EB, 3.802518e-01_EB, 2.635373e-01_EB, &
   1.882440e-01_EB, &
   7.574636e-01_EB, 4.662984e-01_EB, 4.191013e-01_EB, 3.278277e-01_EB, 2.460639e-01_EB, 1.547873e-01_EB, &
   1.018707e-01_EB, &
   3.622467e-02_EB, 3.135298e-02_EB, 1.070575e-02_EB, 1.197498e-02_EB, 2.172120e-02_EB, 1.437764e-02_EB, &
   1.246169e-02_EB, &
   6.790247e-04_EB, 4.477004e-04_EB, 4.003496e-04_EB, 2.372342e-03_EB, 2.112353e-04_EB, 1.899291e-04_EB, &
   4.337606e-04_EB, &
   6.821738e-04_EB, 3.841181e-04_EB, 4.398296e-04_EB, 9.336168e-04_EB, 3.215112e-04_EB, 1.577836e-03_EB, &
   1.387750e-03_EB, &
   5.357513e-03_EB, 9.513971e-03_EB, 9.586573e-03_EB, 1.177708e-02_EB, 1.894358e-02_EB, 1.778693e-02_EB, &
   2.131572e-02_EB, &
   1.558459e-01_EB, 1.086936e-01_EB, 1.035246e-01_EB, 8.484628e-02_EB, 7.171917e-02_EB, 5.505700e-02_EB, &
   3.996848e-02_EB/),(/7,8/))

sd2_c3h6(1:7,25:31) = RESHAPE((/ &  ! 1825-1975 cm-1
   2.622036e-01_EB, 1.516085e-01_EB, 1.303312e-01_EB, 1.109721e-01_EB, 9.231520e-02_EB, 7.419268e-02_EB, &
   4.811877e-02_EB, &
   1.824093e-01_EB, 1.314846e-01_EB, 1.287875e-01_EB, 9.898525e-02_EB, 8.120577e-02_EB, 6.960379e-02_EB, &
   4.479529e-02_EB, &
   4.747023e-02_EB, 3.221598e-02_EB, 3.332202e-02_EB, 2.756103e-02_EB, 2.423984e-02_EB, 2.922304e-02_EB, &
   1.964140e-02_EB, &
   1.531867e-02_EB, 9.915829e-03_EB, 3.275042e-02_EB, 4.962254e-03_EB, 3.664153e-03_EB, 1.736865e-03_EB, &
   2.806612e-03_EB, &
   1.672196e-03_EB, 6.907827e-04_EB, 9.686969e-04_EB, 4.056228e-04_EB, 1.810253e-04_EB, 2.048415e-04_EB, &
   3.860056e-04_EB, &
   4.937814e-04_EB, 1.785822e-04_EB, 3.303675e-04_EB, 1.807431e-04_EB, 1.381312e-04_EB, 1.579684e-04_EB, &
   2.163472e-04_EB, &
   1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, &
   1.027560e-04_EB/),(/7,7/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_c3h6(n_temp_c3h6,31)) 

! band #2: 1225 cm-1 - 1975 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.89939 % 

! print fine structure array gamma_d 

gammad2_c3h6(1:7,1:8) = RESHAPE((/ &  ! 1225-1400 cm-1
   2.622715e-04_EB, 2.415333e-04_EB, 2.285633e-04_EB, 1.138182e-04_EB, 1.023912e-04_EB, 1.165998e-04_EB, &
   1.100001e-04_EB, &
   5.083865e-04_EB, 2.678443e-04_EB, 3.093172e-04_EB, 2.072846e-04_EB, 1.338274e-04_EB, 2.187271e-04_EB, &
   3.932176e-03_EB, &
   5.760161e-04_EB, 3.452147e-02_EB, 1.478839e-04_EB, 3.228596e-03_EB, 3.197521e-03_EB, 2.939038e-05_EB, &
   1.383337e-03_EB, &
   4.544833e-04_EB, 8.199799e-02_EB, 3.788188e-04_EB, 1.070417e-02_EB, 1.812632e-02_EB, 7.216282e-05_EB, &
   3.552942e-03_EB, &
   4.941401e-04_EB, 7.477210e-02_EB, 2.805114e-03_EB, 1.723863e-02_EB, 2.045615e-02_EB, 2.807456e-03_EB, &
   4.838154e+00_EB, &
   3.708619e-03_EB, 1.501171e-01_EB, 2.088233e-02_EB, 5.723047e-02_EB, 2.236513e-01_EB, 5.128140e-02_EB, &
   4.345874e+01_EB, &
   5.695155e-02_EB, 2.842757e-01_EB, 8.565045e-02_EB, 1.664632e-01_EB, 2.596688e-01_EB, 6.733220e-02_EB, &
   4.345955e+01_EB, &
   1.810625e-01_EB, 3.644158e-01_EB, 2.284041e-01_EB, 2.941351e-01_EB, 2.787627e-01_EB, 1.246260e-01_EB, &
   4.345946e+01_EB/),(/7,8/))

gammad2_c3h6(1:7,9:16) = RESHAPE((/ &  ! 1425-1600 cm-1
   2.618227e-01_EB, 4.657918e-01_EB, 2.858630e-01_EB, 3.534168e-01_EB, 2.920606e-01_EB, 2.314596e-01_EB, &
   4.345293e+01_EB, &
   4.735753e-01_EB, 5.151092e-01_EB, 4.261774e-01_EB, 3.947653e-01_EB, 3.773149e-01_EB, 1.853547e-01_EB, &
   4.343606e+01_EB, &
   4.430151e-01_EB, 5.018903e-01_EB, 3.534603e-01_EB, 4.029514e-01_EB, 3.684172e-01_EB, 1.453836e-01_EB, &
   4.345940e+01_EB, &
   1.212938e-01_EB, 3.173309e-01_EB, 1.209401e+00_EB, 3.624133e-01_EB, 2.461181e-01_EB, 2.905692e-01_EB, &
   4.345952e+01_EB, &
   1.084288e-02_EB, 1.935743e-01_EB, 3.826811e-01_EB, 2.209504e-01_EB, 3.384784e-01_EB, 5.030323e+00_EB, &
   4.314517e+00_EB, &
   3.502530e-03_EB, 4.235158e-01_EB, 4.492198e+01_EB, 4.178424e-01_EB, 3.215097e-01_EB, 1.976012e+01_EB, &
   3.306259e+01_EB, &
   9.195743e-04_EB, 5.344729e-02_EB, 3.856228e-02_EB, 8.090926e-02_EB, 1.993074e-01_EB, 8.752817e-02_EB, &
   4.345936e+01_EB, &
   5.991440e-04_EB, 1.174575e-01_EB, 4.116301e-02_EB, 1.163800e-01_EB, 2.048843e-01_EB, 1.170091e-01_EB, &
   4.345936e+01_EB/),(/7,8/))

gammad2_c3h6(1:7,17:24) = RESHAPE((/ &  ! 1625-1800 cm-1
   1.048318e-01_EB, 4.077170e-01_EB, 3.250763e-01_EB, 4.204265e-01_EB, 4.172624e-01_EB, 3.501140e-01_EB, &
   4.345754e+01_EB, &
   3.891245e-01_EB, 5.297782e-01_EB, 5.753653e-01_EB, 4.671660e-01_EB, 3.550334e-01_EB, 5.129497e-01_EB, &
   1.852430e+01_EB, &
   1.315391e-01_EB, 2.233168e-01_EB, 1.643737e-01_EB, 2.846889e-01_EB, 2.529524e-01_EB, 5.109478e-01_EB, &
   5.217994e+00_EB, &
   2.202200e-04_EB, 6.992095e-05_EB, 6.651429e-02_EB, 1.261036e-01_EB, 2.793432e-04_EB, 3.400065e-01_EB, &
   2.976967e-01_EB, &
   4.407006e-04_EB, 3.060631e-04_EB, 3.078173e-04_EB, 1.997468e-02_EB, 1.414429e-04_EB, 1.272939e-04_EB, &
   2.920596e-04_EB, &
   4.075042e-04_EB, 2.572692e-04_EB, 2.809507e-04_EB, 7.570268e-04_EB, 2.248097e-04_EB, 6.557878e-04_EB, &
   1.214879e-02_EB, &
   5.533543e-04_EB, 2.280793e-03_EB, 5.727351e-03_EB, 1.258039e-01_EB, 2.877495e-03_EB, 1.528211e-01_EB, &
   3.641746e-01_EB, &
   2.750324e-02_EB, 5.042429e-02_EB, 3.703191e-02_EB, 1.112856e-01_EB, 1.394191e-01_EB, 2.852003e-01_EB, &
   5.289596e+00_EB/),(/7,8/))

gammad2_c3h6(1:7,25:31) = RESHAPE((/ &  ! 1825-1975 cm-1
   1.249123e-01_EB, 1.287617e-01_EB, 8.856076e-02_EB, 1.063268e-01_EB, 6.381303e-02_EB, 3.721067e-02_EB, &
   1.078497e+01_EB, &
   9.224931e-02_EB, 7.783664e-02_EB, 4.282341e-02_EB, 1.077656e-01_EB, 9.668966e-02_EB, 3.557773e-02_EB, &
   4.459730e+00_EB, &
   1.094300e-02_EB, 1.484660e-02_EB, 1.474169e-02_EB, 2.868794e-02_EB, 2.363462e-02_EB, 7.563694e-03_EB, &
   9.480691e-01_EB, &
   3.143676e-03_EB, 1.423552e-03_EB, 2.803890e-04_EB, 1.197387e-02_EB, 1.810420e-03_EB, 9.517693e-04_EB, &
   1.253779e-03_EB, &
   9.700043e-04_EB, 5.310996e-04_EB, 3.238508e-04_EB, 2.462185e-04_EB, 1.221073e-04_EB, 1.352439e-04_EB, &
   2.806980e-04_EB, &
   3.482794e-04_EB, 1.171460e-04_EB, 2.319252e-04_EB, 1.036437e-04_EB, 1.041061e-04_EB, 1.117176e-04_EB, &
   1.472074e-04_EB, &
   8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, &
   8.047518e-05_EB/),(/7,7/))

!---------------------------------------------------------------------------
ALLOCATE(sd3_c3h6(n_temp_c3h6,26)) 

! band #3: 2650 cm-1 - 3275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.84732 % 

sd3_c3h6(1:7,1:8) = RESHAPE((/ &  ! 2650-2825 cm-1
   9.842396e-04_EB, 9.773581e-04_EB, 4.737803e-04_EB, 1.731125e-04_EB, 1.636141e-04_EB, 1.470685e-04_EB, &
   1.906099e-04_EB, &
   2.026987e-03_EB, 1.186057e-03_EB, 9.296802e-04_EB, 1.762232e-03_EB, 6.242290e-04_EB, 2.052893e-03_EB, &
   1.454807e-02_EB, &
   1.375163e-02_EB, 8.857038e-03_EB, 1.157316e-02_EB, 1.514073e-02_EB, 7.647375e-03_EB, 2.028955e-02_EB, &
   7.263667e-02_EB, &
   6.091687e-02_EB, 3.794720e-02_EB, 3.692778e-02_EB, 3.669497e-02_EB, 2.246508e-02_EB, 2.837078e-02_EB, &
   4.660372e-02_EB, &
   7.796068e-02_EB, 5.480320e-02_EB, 5.466178e-02_EB, 5.561425e-02_EB, 3.895555e-02_EB, 4.987873e-02_EB, &
   5.885558e-02_EB, &
   5.994425e-02_EB, 4.512610e-02_EB, 5.182436e-02_EB, 5.600091e-02_EB, 4.405629e-02_EB, 5.968023e-02_EB, &
   9.479119e-02_EB, &
   9.789297e-02_EB, 7.603902e-02_EB, 8.316338e-02_EB, 8.497910e-02_EB, 7.201492e-02_EB, 9.033387e-02_EB, &
   1.461816e-01_EB, &
   1.180654e-01_EB, 1.124261e-01_EB, 1.279063e-01_EB, 1.423319e-01_EB, 1.323416e-01_EB, 1.543846e-01_EB, &
   1.754013e-01_EB/),(/7,8/))

sd3_c3h6(1:7,9:16) = RESHAPE((/ &  ! 2850-3025 cm-1
   4.934323e-01_EB, 3.715076e-01_EB, 3.606179e-01_EB, 3.406474e-01_EB, 2.903159e-01_EB, 2.686381e-01_EB, &
   2.632060e-01_EB, &
   9.542159e-01_EB, 7.012988e-01_EB, 6.454962e-01_EB, 5.929394e-01_EB, 5.042463e-01_EB, 4.353109e-01_EB, &
   4.219635e-01_EB, &
   1.390515e+00_EB, 1.083207e+00_EB, 1.007094e+00_EB, 9.304473e-01_EB, 8.055580e-01_EB, 6.623870e-01_EB, &
   5.866959e-01_EB, &
   2.165889e+00_EB, 1.557403e+00_EB, 1.377848e+00_EB, 1.234664e+00_EB, 1.013037e+00_EB, 7.883431e-01_EB, &
   6.628056e-01_EB, &
   2.955798e+00_EB, 2.081509e+00_EB, 1.819270e+00_EB, 1.598535e+00_EB, 1.286353e+00_EB, 9.640864e-01_EB, &
   7.596691e-01_EB, &
   2.749633e+00_EB, 1.810398e+00_EB, 1.545370e+00_EB, 1.334588e+00_EB, 1.045362e+00_EB, 8.023541e-01_EB, &
   6.307507e-01_EB, &
   2.482580e+00_EB, 1.760628e+00_EB, 1.526249e+00_EB, 1.327667e+00_EB, 1.049659e+00_EB, 7.936791e-01_EB, &
   6.269325e-01_EB, &
   1.476540e+00_EB, 1.169856e+00_EB, 1.094611e+00_EB, 9.890179e-01_EB, 8.496338e-01_EB, 7.194488e-01_EB, &
   6.198264e-01_EB/),(/7,8/))

sd3_c3h6(1:7,17:24) = RESHAPE((/ &  ! 3050-3225 cm-1
   1.158768e+00_EB, 9.970946e-01_EB, 9.717818e-01_EB, 9.018919e-01_EB, 8.089673e-01_EB, 6.965923e-01_EB, &
   5.756752e-01_EB, &
   1.617699e+00_EB, 1.161452e+00_EB, 1.056867e+00_EB, 9.368936e-01_EB, 7.790857e-01_EB, 6.211223e-01_EB, &
   4.978033e-01_EB, &
   1.668409e+00_EB, 1.154363e+00_EB, 1.028890e+00_EB, 8.900880e-01_EB, 7.204607e-01_EB, 5.621992e-01_EB, &
   4.193038e-01_EB, &
   6.627374e-01_EB, 5.511257e-01_EB, 5.408762e-01_EB, 4.741726e-01_EB, 4.065990e-01_EB, 3.500514e-01_EB, &
   2.603275e-01_EB, &
   1.943149e-01_EB, 1.812675e-01_EB, 1.997641e-01_EB, 1.725390e-01_EB, 1.584383e-01_EB, 1.668560e-01_EB, &
   1.239099e-01_EB, &
   3.848101e-02_EB, 4.502827e-02_EB, 8.182965e-02_EB, 6.478508e-02_EB, 6.632017e-02_EB, 7.060523e-02_EB, &
   5.939407e-02_EB, &
   1.823763e-02_EB, 6.889882e-03_EB, 9.139773e-02_EB, 8.162446e-02_EB, 2.475465e-02_EB, 3.168406e-02_EB, &
   1.903134e-02_EB, &
   6.725369e-04_EB, 3.888602e-04_EB, 1.980271e-01_EB, 3.039167e-02_EB, 1.104407e-02_EB, 1.220321e-02_EB, &
   4.484167e-03_EB/),(/7,8/))

sd3_c3h6(1:7,25:26) = RESHAPE((/ &  ! 3250-3275 cm-1
   4.840794e-04_EB, 3.312870e-04_EB, 5.789111e-03_EB, 3.214095e-03_EB, 3.311342e-03_EB, 2.947009e-03_EB, &
   1.302455e-04_EB, &
   1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, 1.027560e-04_EB, &
   1.027560e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad3_c3h6(n_temp_c3h6,26)) 

! band #3: 2650 cm-1 - 3275 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.84732 % 

! print fine structure array gamma_d 

gammad3_c3h6(1:7,1:8) = RESHAPE((/ &  ! 2650-2825 cm-1
   3.045289e-04_EB, 1.015641e-03_EB, 3.935463e-04_EB, 1.154192e-04_EB, 1.194489e-04_EB, 1.102285e-04_EB, &
   1.137523e-04_EB, &
   5.888900e-04_EB, 1.065118e-03_EB, 5.931355e-04_EB, 8.588542e-04_EB, 5.730959e-04_EB, 5.277001e-04_EB, &
   2.327868e-03_EB, &
   3.553832e-03_EB, 5.225375e-02_EB, 5.508119e-03_EB, 3.764204e-03_EB, 6.360672e-02_EB, 6.490060e-04_EB, &
   1.126737e-03_EB, &
   3.234273e-02_EB, 6.202832e-02_EB, 1.838655e-02_EB, 1.292609e-02_EB, 1.730868e-01_EB, 7.897375e-03_EB, &
   1.175871e-02_EB, &
   5.722213e-02_EB, 8.665120e-02_EB, 3.332681e-02_EB, 1.901355e-02_EB, 3.938769e-02_EB, 7.412539e-03_EB, &
   2.819413e-02_EB, &
   3.022697e-02_EB, 1.129809e-01_EB, 2.531468e-02_EB, 1.970225e-02_EB, 5.471032e-02_EB, 2.097209e-02_EB, &
   2.390018e-02_EB, &
   5.826107e-02_EB, 1.114796e-01_EB, 3.749457e-02_EB, 3.774872e-02_EB, 1.009734e-01_EB, 4.970945e-02_EB, &
   3.868475e-02_EB, &
   7.083218e-02_EB, 1.686291e-01_EB, 7.156561e-02_EB, 6.619537e-02_EB, 2.202441e-01_EB, 1.169745e-01_EB, &
   3.834597e-01_EB/),(/7,8/))

gammad3_c3h6(1:7,9:16) = RESHAPE((/ &  ! 2850-3025 cm-1
   2.387749e-01_EB, 3.465833e-01_EB, 1.768091e-01_EB, 2.165351e-01_EB, 3.664688e-01_EB, 2.483037e-01_EB, &
   8.926665e-01_EB, &
   5.324178e-01_EB, 5.980304e-01_EB, 3.793473e-01_EB, 4.509406e-01_EB, 5.538922e-01_EB, 4.079465e-01_EB, &
   5.945021e-01_EB, &
   8.715520e-01_EB, 9.293503e-01_EB, 6.306967e-01_EB, 7.198743e-01_EB, 7.551117e-01_EB, 6.154785e-01_EB, &
   8.223511e-01_EB, &
   1.499464e+00_EB, 1.399873e+00_EB, 9.866448e-01_EB, 1.025786e+00_EB, 9.581477e-01_EB, 7.125299e-01_EB, &
   1.107012e+00_EB, &
   2.156257e+00_EB, 1.637737e+00_EB, 1.241655e+00_EB, 1.280306e+00_EB, 1.138186e+00_EB, 8.441688e-01_EB, &
   1.384899e+00_EB, &
   3.002594e+00_EB, 1.714061e+00_EB, 1.257507e+00_EB, 1.172490e+00_EB, 1.028409e+00_EB, 6.569541e-01_EB, &
   2.100100e+00_EB, &
   2.734771e+00_EB, 1.525371e+00_EB, 1.227801e+00_EB, 1.124537e+00_EB, 1.012834e+00_EB, 6.561761e-01_EB, &
   1.786964e+00_EB, &
   1.391166e+00_EB, 1.039743e+00_EB, 8.609837e-01_EB, 8.549996e-01_EB, 8.423802e-01_EB, 6.266077e-01_EB, &
   1.102257e+00_EB/),(/7,8/))

gammad3_c3h6(1:7,17:24) = RESHAPE((/ &  ! 3050-3225 cm-1
   8.645946e-01_EB, 9.364314e-01_EB, 7.448637e-01_EB, 7.668057e-01_EB, 7.630071e-01_EB, 5.882839e-01_EB, &
   1.479949e+00_EB, &
   1.230685e+00_EB, 1.046279e+00_EB, 7.526449e-01_EB, 7.500084e-01_EB, 7.007694e-01_EB, 5.329368e-01_EB, &
   1.354177e+00_EB, &
   1.587192e+00_EB, 9.549637e-01_EB, 7.345138e-01_EB, 6.867975e-01_EB, 5.981748e-01_EB, 4.326344e-01_EB, &
   4.822820e+00_EB, &
   6.944559e-01_EB, 3.901526e-01_EB, 3.458891e-01_EB, 3.287273e-01_EB, 2.870524e-01_EB, 2.152839e-01_EB, &
   4.345222e+01_EB, &
   1.550994e-01_EB, 1.096247e-01_EB, 1.314631e-01_EB, 1.210936e-01_EB, 1.139491e-01_EB, 5.956151e-02_EB, &
   3.313368e+01_EB, &
   2.202410e-02_EB, 2.065956e-02_EB, 1.775205e-02_EB, 2.386698e-02_EB, 3.888119e-02_EB, 2.898450e-02_EB, &
   2.058336e+00_EB, &
   1.065388e-04_EB, 1.920534e-03_EB, 6.158394e-04_EB, 4.540651e-04_EB, 1.753716e-02_EB, 7.274816e-03_EB, &
   1.007187e+00_EB, &
   4.230956e-04_EB, 2.525852e-04_EB, 8.709782e-05_EB, 2.954917e-04_EB, 5.762183e-03_EB, 3.267777e-03_EB, &
   1.887507e-02_EB/),(/7,8/))

gammad3_c3h6(1:7,25:26) = RESHAPE((/ &  ! 3250-3275 cm-1
   3.563014e-04_EB, 2.388188e-04_EB, 5.650871e-04_EB, 1.353722e-03_EB, 1.592126e-03_EB, 7.501643e-04_EB, &
   1.070038e-04_EB, &
   8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, 8.047518e-05_EB, &
   8.047518e-05_EB/),(/7,2/))


!-------------------------toluene data-------------------


! there are 5 bands for toluene

! band #1: 700 cm-1 - 805 cm-1 
! band #2: 975 cm-1 - 1175 cm-1 
! band #3: 1275 cm-1 - 1675 cm-1 
! band #4: 1650 cm-1 - 2075 cm-1 
! band #5: 2675 cm-1 - 3225 cm-1 

ALLOCATE(sd_c7h8_temp(n_temp_c7h8)) 

! initialize bands wavenumber bounds for toluene ABS(option coefficients.
! bands are organized by row: 
! 1st column is the lower bound band limit in cm-1. 
! 2nd column is the upper bound band limit in cm-1. 
! 3rd column is the stride between wavenumbers in band.
! IF 3rd column = 0., band is calculated, not tabulated.

ALLOCATE(om_bnd_c7h8(n_band_c7h8,3)) 
ALLOCATE(be_c7h8(n_band_c7h8)) 

om_bnd_c7h8 = RESHAPE((/ &
   700._EB,  975._EB, 1275._EB, 1650._EB, 2675._EB, &
   805._EB, 1175._EB, 1675._EB, 2075._EB, 3225._EB, &
   5._EB, 5._EB, 25._EB, 25._EB, 25._EB/),(/n_band_c7h8,3/)) 

sd_c7h8_temp = (/ &
   300._EB, 396._EB, 440._EB, 477._EB, 587._EB, 795._EB,&
   999._EB/)

be_c7h8 = (/ &
   1.000_EB, 1.000_EB, 1.000_EB, 1.000_EB, 1.000_EB/)

!---------------------------------------------------------------------------
ALLOCATE(sd1_c7h8(n_temp_c7h8,22)) 

! band #1: 700 cm-1 - 805 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 5.5814 % 

sd1_c7h8(1:7,1:8) = RESHAPE((/ &  ! 700-735 cm-1
   1.112466e+00_EB, 9.044923e-01_EB, 9.061370e-01_EB, 8.938794e-01_EB, 9.084825e-01_EB, 1.435389e+00_EB, &
   1.297156e+01_EB, &
   2.840815e+00_EB, 2.482407e+00_EB, 2.360535e+00_EB, 2.276458e+00_EB, 2.122272e+00_EB, 1.940697e+00_EB, &
   5.112554e+00_EB, &
   3.455206e+00_EB, 3.166147e+00_EB, 2.863080e+00_EB, 2.620674e+00_EB, 2.467633e+00_EB, 2.005290e+00_EB, &
   7.244894e-01_EB, &
   4.288216e+00_EB, 3.577673e+00_EB, 2.963877e+00_EB, 2.789847e+00_EB, 2.401611e+00_EB, 1.912466e+00_EB, &
   2.024019e+00_EB, &
   4.666515e+00_EB, 3.754488e+00_EB, 3.254065e+00_EB, 3.058494e+00_EB, 2.875391e+00_EB, 2.451298e+00_EB, &
   4.223163e+00_EB, &
   7.363590e+00_EB, 6.122925e+00_EB, 5.403265e+00_EB, 5.328257e+00_EB, 4.641149e+00_EB, 3.466727e+00_EB, &
   1.373953e+01_EB, &
   9.587057e+00_EB, 6.742653e+00_EB, 5.600793e+00_EB, 5.363171e+00_EB, 4.064656e+00_EB, 2.508085e+00_EB, &
   1.195299e+00_EB, &
   4.715829e+00_EB, 3.108857e+00_EB, 2.753587e+00_EB, 2.650163e+00_EB, 2.249724e+00_EB, 1.963447e+00_EB, &
   7.262488e+00_EB/),(/7,8/))

sd1_c7h8(1:7,9:16) = RESHAPE((/ &  ! 740-775 cm-1
   4.190702e+00_EB, 2.766122e+00_EB, 2.537781e+00_EB, 2.386878e+00_EB, 2.084718e+00_EB, 1.344226e+00_EB, &
   3.448390e+00_EB, &
   2.838874e+00_EB, 2.217681e+00_EB, 2.081071e+00_EB, 1.986804e+00_EB, 1.810042e+00_EB, 1.549486e+00_EB, &
   7.126841e-01_EB, &
   1.327883e+00_EB, 1.440906e+00_EB, 1.390571e+00_EB, 1.421333e+00_EB, 1.328347e+00_EB, 1.296316e+00_EB, &
   1.254896e+00_EB, &
   4.533944e-01_EB, 8.351739e-01_EB, 7.591511e-01_EB, 7.962784e-01_EB, 8.186310e-01_EB, 8.370456e-01_EB, &
   2.685634e+00_EB, &
   1.131517e-01_EB, 7.896344e-01_EB, 3.378503e-01_EB, 3.523472e-01_EB, 4.364146e-01_EB, 5.594990e-01_EB, &
   3.259981e+00_EB, &
   3.791615e-02_EB, 8.912653e-01_EB, 4.878342e-01_EB, 1.492598e-01_EB, 2.103530e-01_EB, 3.269232e-01_EB, &
   2.144533e+00_EB, &
   4.206761e-02_EB, 7.530056e-02_EB, 1.233550e-01_EB, 7.429326e-02_EB, 9.115926e-02_EB, 1.525590e-01_EB, &
   1.109247e+00_EB, &
   6.711521e-02_EB, 8.510326e-04_EB, 6.739778e-02_EB, 5.411094e-02_EB, 6.504658e-02_EB, 5.754676e-02_EB, &
   5.875940e-01_EB/),(/7,8/))

sd1_c7h8(1:7,17:22) = RESHAPE((/ &  ! 780-805 cm-1
   6.275862e-02_EB, 6.962504e-04_EB, 4.265822e-02_EB, 2.937821e-02_EB, 3.183374e-02_EB, 1.453078e-02_EB, &
   8.515825e-02_EB, &
   5.193232e-02_EB, 7.365615e-04_EB, 1.939785e-02_EB, 1.719349e-02_EB, 2.278960e-02_EB, 1.331369e-03_EB, &
   5.520836e-02_EB, &
   3.535671e-02_EB, 8.579228e-04_EB, 1.765940e-03_EB, 1.491081e-02_EB, 1.181497e-02_EB, 5.221427e-04_EB, &
   2.334611e-03_EB, &
   6.188348e-03_EB, 7.316561e-04_EB, 8.181678e-04_EB, 6.095900e-03_EB, 3.380075e-03_EB, 5.968569e-04_EB, &
   7.411591e-04_EB, &
   6.328121e-04_EB, 6.599274e-04_EB, 5.510258e-04_EB, 5.104659e-04_EB, 5.023646e-04_EB, 5.757667e-04_EB, &
   6.209242e-04_EB, &
   4.621323e-04_EB, 4.638592e-04_EB, 4.665860e-04_EB, 4.636519e-04_EB, 4.646393e-04_EB, 4.633936e-04_EB, &
   4.615859e-04_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(gammad1_c7h8(n_temp_c7h8,22)) 

! band #1: 700 cm-1 - 805 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 5.5814 % 

! print fine structure array gamma_d 

gammad1_c7h8(1:7,1:8) = RESHAPE((/ &  ! 700-735 cm-1
   7.999461e+01_EB, 6.960752e+01_EB, 3.022045e+01_EB, 6.344408e+01_EB, 5.719128e+01_EB, 3.598038e-02_EB, &
   1.819109e-03_EB, &
   7.999993e+01_EB, 6.962945e+01_EB, 6.605713e+01_EB, 6.344357e+01_EB, 5.719074e+01_EB, 3.344331e-01_EB, &
   1.001242e-02_EB, &
   8.000000e+01_EB, 6.922109e+01_EB, 4.454022e+00_EB, 6.344038e+01_EB, 5.719013e+01_EB, 9.751497e-02_EB, &
   4.383971e+01_EB, &
   8.000000e+01_EB, 1.649903e+01_EB, 6.605768e+01_EB, 6.344412e+01_EB, 5.719149e+01_EB, 6.787669e-01_EB, &
   1.717471e-01_EB, &
   8.000000e+01_EB, 5.534136e+00_EB, 6.605782e+01_EB, 6.344412e+01_EB, 5.719151e+01_EB, 4.914284e+01_EB, &
   6.407876e-02_EB, &
   8.000000e+01_EB, 6.962930e+01_EB, 6.605773e+01_EB, 8.838043e+00_EB, 5.719150e+01_EB, 2.254685e+00_EB, &
   1.178609e-02_EB, &
   8.000000e+01_EB, 6.963055e+01_EB, 6.605773e+01_EB, 3.115435e+00_EB, 5.719151e+01_EB, 4.914358e+01_EB, &
   4.383966e+01_EB, &
   8.000000e+01_EB, 6.963104e+01_EB, 6.605612e+01_EB, 6.343942e+01_EB, 5.719109e+01_EB, 1.671127e-01_EB, &
   6.618119e-03_EB/),(/7,8/))

gammad1_c7h8(1:7,9:16) = RESHAPE((/ &  ! 740-775 cm-1
   8.000000e+01_EB, 6.963105e+01_EB, 6.605690e+01_EB, 6.333100e+01_EB, 5.718653e+01_EB, 4.913467e+01_EB, &
   1.209833e-02_EB, &
   8.000000e+01_EB, 6.962982e+01_EB, 6.605767e+01_EB, 6.342460e+01_EB, 5.718540e+01_EB, 2.011652e-01_EB, &
   4.383963e+01_EB, &
   7.999980e+01_EB, 6.960747e+01_EB, 6.602987e+01_EB, 6.326339e+01_EB, 5.714448e+01_EB, 3.303355e-01_EB, &
   2.494003e-02_EB, &
   4.634778e+01_EB, 1.999423e-01_EB, 1.445536e+01_EB, 3.851197e+01_EB, 5.716296e+01_EB, 2.419503e+00_EB, &
   4.403455e-03_EB, &
   1.788269e+00_EB, 8.781616e-03_EB, 1.079139e-01_EB, 1.421403e+01_EB, 9.057111e+00_EB, 1.781522e+01_EB, &
   1.669472e-03_EB, &
   2.864129e-01_EB, 1.515581e-03_EB, 2.119890e-03_EB, 5.013687e+00_EB, 3.747774e+00_EB, 1.844443e+00_EB, &
   1.069786e-03_EB, &
   1.102908e+01_EB, 4.896633e-02_EB, 1.798393e-02_EB, 6.492346e-01_EB, 8.863603e-02_EB, 2.773305e-01_EB, &
   5.841296e-04_EB, &
   1.251724e+00_EB, 6.490441e-04_EB, 2.113061e-02_EB, 1.326897e+00_EB, 8.049318e-03_EB, 8.055391e-03_EB, &
   2.478640e-04_EB/),(/7,8/))

gammad1_c7h8(1:7,17:22) = RESHAPE((/ &  ! 780-805 cm-1
   6.819808e-01_EB, 5.093301e-04_EB, 2.666337e-02_EB, 2.101043e-01_EB, 1.980342e-02_EB, 1.679664e-03_EB, &
   1.906432e-02_EB, &
   4.691940e-01_EB, 5.429738e-04_EB, 1.319467e-02_EB, 4.975402e-03_EB, 2.433969e-02_EB, 1.884068e-03_EB, &
   5.887236e-03_EB, &
   1.875981e-01_EB, 6.598198e-04_EB, 1.786636e-03_EB, 1.408548e-02_EB, 2.119727e-01_EB, 4.335098e-04_EB, &
   1.937100e-03_EB, &
   1.442921e-02_EB, 5.711043e-04_EB, 5.971323e-04_EB, 7.973510e-03_EB, 3.341200e-03_EB, 4.314868e-04_EB, &
   5.200300e-04_EB, &
   5.156874e-04_EB, 5.125082e-04_EB, 4.318513e-04_EB, 3.992184e-04_EB, 3.936730e-04_EB, 4.154189e-04_EB, &
   4.736124e-04_EB, &
   3.682775e-04_EB, 3.700575e-04_EB, 3.724614e-04_EB, 3.699127e-04_EB, 3.709688e-04_EB, 3.696394e-04_EB, &
   3.680801e-04_EB/),(/7,6/))

!---------------------------------------------------------------------------
ALLOCATE(sd2_c7h8(n_temp_c7h8,29)) 

! band #2: 975 cm-1 - 1175 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.91355 % 

sd2_c7h8(1:7,1:8) = RESHAPE((/ &  ! 975-1010 cm-1
   5.806929e-04_EB, 5.942617e-04_EB, 6.726459e-04_EB, 7.052569e-04_EB, 6.695985e-04_EB, 5.159322e-04_EB, &
   5.551144e-04_EB, &
   6.785508e-04_EB, 8.056701e-04_EB, 8.454611e-04_EB, 7.340346e-04_EB, 7.540252e-04_EB, 6.231487e-04_EB, &
   1.161712e-02_EB, &
   7.608293e-04_EB, 5.613980e-02_EB, 8.330584e-04_EB, 7.050990e-04_EB, 7.281009e-04_EB, 1.902747e-03_EB, &
   2.699830e-02_EB, &
   4.497896e-02_EB, 9.748044e-02_EB, 5.769989e-04_EB, 8.307786e-04_EB, 1.048030e-03_EB, 3.071367e-02_EB, &
   5.746590e-02_EB, &
   4.185662e-02_EB, 5.807183e-02_EB, 5.845956e-03_EB, 4.738381e-03_EB, 1.373288e-02_EB, 1.722897e-02_EB, &
   6.537832e-02_EB, &
   2.338216e-02_EB, 1.675285e-02_EB, 2.368298e-02_EB, 3.338373e-02_EB, 1.128221e-02_EB, 2.346587e-02_EB, &
   3.565197e-02_EB, &
   1.379198e-02_EB, 9.132916e-04_EB, 3.690042e-02_EB, 2.275924e-02_EB, 2.434022e-02_EB, 4.242178e-02_EB, &
   6.940845e-02_EB, &
   2.629167e-02_EB, 2.433785e-03_EB, 7.441229e-02_EB, 7.329327e-02_EB, 8.154865e-02_EB, 8.520170e-02_EB, &
   8.926539e-02_EB/),(/7,8/))

sd2_c7h8(1:7,9:16) = RESHAPE((/ &  ! 1015-1050 cm-1
   8.643728e-02_EB, 8.397464e-02_EB, 1.634758e-01_EB, 1.741891e-01_EB, 1.714205e-01_EB, 1.480038e-01_EB, &
   1.360983e-01_EB, &
   3.004729e-01_EB, 2.793190e-01_EB, 2.877798e-01_EB, 2.821584e-01_EB, 2.456569e-01_EB, 1.463510e-01_EB, &
   1.826942e-01_EB, &
   5.278878e-01_EB, 3.693913e-01_EB, 2.973200e-01_EB, 2.717632e-01_EB, 2.438624e-01_EB, 1.379726e-01_EB, &
   2.050903e-01_EB, &
   5.845808e-01_EB, 4.794762e-01_EB, 3.584954e-01_EB, 3.427664e-01_EB, 2.786454e-01_EB, 1.226379e-01_EB, &
   2.097466e-01_EB, &
   4.729738e-01_EB, 3.788669e-01_EB, 2.664993e-01_EB, 2.564802e-01_EB, 2.264334e-01_EB, 1.535049e-01_EB, &
   1.682791e-01_EB, &
   5.530868e-01_EB, 4.490113e-01_EB, 3.310727e-01_EB, 3.288768e-01_EB, 3.019277e-01_EB, 2.043861e-01_EB, &
   1.497158e-01_EB, &
   4.129235e-01_EB, 3.251444e-01_EB, 2.714955e-01_EB, 2.657956e-01_EB, 2.487687e-01_EB, 1.862679e-01_EB, &
   1.223084e-01_EB, &
   2.614169e-01_EB, 1.802809e-01_EB, 1.975745e-01_EB, 1.888234e-01_EB, 2.189458e-01_EB, 1.348280e-01_EB, &
   9.276862e-02_EB/),(/7,8/))

sd2_c7h8(1:7,17:24) = RESHAPE((/ &  ! 1055-1090 cm-1
   2.730370e-01_EB, 1.767450e-01_EB, 2.218796e-01_EB, 2.155878e-01_EB, 2.223112e-01_EB, 1.469945e-01_EB, &
   1.150109e-01_EB, &
   2.816011e-01_EB, 2.052777e-01_EB, 2.522512e-01_EB, 2.342314e-01_EB, 2.280630e-01_EB, 2.544731e-01_EB, &
   1.975636e-01_EB, &
   3.260058e-01_EB, 2.887673e-01_EB, 2.910159e-01_EB, 2.801621e-01_EB, 2.791987e-01_EB, 2.761989e-01_EB, &
   6.486488e-01_EB, &
   4.240289e-01_EB, 3.848016e-01_EB, 3.433772e-01_EB, 3.147066e-01_EB, 2.873021e-01_EB, 2.027784e-01_EB, &
   6.801384e-01_EB, &
   5.531291e-01_EB, 4.341564e-01_EB, 3.312172e-01_EB, 3.045238e-01_EB, 2.620640e-01_EB, 1.519635e-01_EB, &
   9.321833e-01_EB, &
   4.365089e-01_EB, 3.538039e-01_EB, 2.730701e-01_EB, 2.397118e-01_EB, 2.227185e-01_EB, 2.009481e-01_EB, &
   1.174461e+00_EB, &
   4.683512e-01_EB, 3.712916e-01_EB, 2.862375e-01_EB, 2.540628e-01_EB, 2.485619e-01_EB, 2.060439e-01_EB, &
   9.951951e-01_EB, &
   4.698469e-01_EB, 3.518932e-01_EB, 2.825430e-01_EB, 2.624504e-01_EB, 2.213420e-01_EB, 1.625252e-01_EB, &
   5.760843e-01_EB/),(/7,8/))

sd2_c7h8(1:7,25:29) = RESHAPE((/ &  ! 1095-1175 cm-1
   3.649981e-01_EB, 2.736054e-01_EB, 2.408070e-01_EB, 2.323294e-01_EB, 1.878487e-01_EB, 1.252916e-01_EB, &
   4.029469e-02_EB, &
   2.681413e-01_EB, 1.701855e-01_EB, 1.674993e-01_EB, 1.511744e-01_EB, 1.166087e-01_EB, 7.290456e-02_EB, &
   1.474340e-01_EB, &
   4.990165e-02_EB, 6.292687e-02_EB, 4.665808e-02_EB, 3.377088e-02_EB, 2.693981e-02_EB, 9.990414e-03_EB, &
   9.979516e-03_EB, &
   6.130796e-04_EB, 1.597828e-03_EB, 5.178852e-03_EB, 6.236302e-04_EB, 6.358216e-04_EB, 6.358332e-04_EB, &
   5.729829e-04_EB, &
   4.621323e-04_EB, 4.638592e-04_EB, 4.665860e-04_EB, 4.636519e-04_EB, 4.646393e-04_EB, 4.633936e-04_EB, &
   4.615859e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(gammad2_c7h8(n_temp_c7h8,29)) 

! band #2: 975 cm-1 - 1175 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.91355 % 

! print fine structure array gamma_d 

gammad2_c7h8(1:7,1:8) = RESHAPE((/ &  ! 975-1010 cm-1
   4.305850e-04_EB, 4.663364e-04_EB, 5.289580e-04_EB, 5.402269e-04_EB, 5.024478e-04_EB, 4.070836e-04_EB, &
   4.471953e-04_EB, &
   5.202785e-04_EB, 6.465361e-04_EB, 6.452178e-04_EB, 5.771613e-04_EB, 5.746126e-04_EB, 4.735889e-04_EB, &
   1.573544e-04_EB, &
   5.872521e-04_EB, 7.042066e-03_EB, 5.999286e-04_EB, 5.331992e-04_EB, 6.220132e-04_EB, 6.364420e-02_EB, &
   1.528351e-04_EB, &
   5.511126e-01_EB, 1.497331e-02_EB, 4.594486e-04_EB, 5.729761e-04_EB, 7.713361e-04_EB, 7.130935e-03_EB, &
   1.427152e-04_EB, &
   2.212250e+00_EB, 1.945735e-01_EB, 2.212567e-01_EB, 9.922155e-03_EB, 9.234920e-02_EB, 1.013608e-02_EB, &
   4.272371e-04_EB, &
   4.405486e-01_EB, 7.436540e-01_EB, 1.920483e-01_EB, 2.962863e-03_EB, 3.448915e-01_EB, 3.486423e-01_EB, &
   1.313079e-01_EB, &
   2.410517e-01_EB, 1.039535e-03_EB, 1.254824e-01_EB, 1.067258e-01_EB, 2.942350e-01_EB, 5.806212e-01_EB, &
   1.774699e+00_EB, &
   4.514842e-01_EB, 1.759071e-01_EB, 2.902656e-01_EB, 2.267366e-01_EB, 8.401946e-02_EB, 7.681227e+00_EB, &
   1.330635e+00_EB/),(/7,8/))

gammad2_c7h8(1:7,9:16) = RESHAPE((/ &  ! 1015-1050 cm-1
   2.224412e+00_EB, 1.484769e+00_EB, 1.494944e+00_EB, 1.832759e+00_EB, 3.314911e+00_EB, 5.865543e+00_EB, &
   6.106454e-01_EB, &
   2.999242e+01_EB, 7.193297e+00_EB, 3.422464e+00_EB, 4.614503e+00_EB, 2.188441e+00_EB, 7.394814e+00_EB, &
   3.908543e-02_EB, &
   7.999990e+01_EB, 3.581705e+00_EB, 7.707850e+00_EB, 3.420693e+00_EB, 6.146487e+00_EB, 7.025651e+00_EB, &
   3.617788e-02_EB, &
   7.998352e+01_EB, 4.444656e-01_EB, 1.228523e+01_EB, 3.433637e+00_EB, 2.072175e+01_EB, 4.335871e+00_EB, &
   7.248922e-02_EB, &
   7.410363e+01_EB, 5.532136e-01_EB, 5.267901e+00_EB, 5.384235e+00_EB, 5.570719e+00_EB, 7.078253e+00_EB, &
   7.681535e+00_EB, &
   7.999988e+01_EB, 1.364684e+01_EB, 1.337888e+01_EB, 5.824947e+00_EB, 5.724261e+00_EB, 6.510348e+00_EB, &
   2.593089e+00_EB, &
   5.835221e+01_EB, 2.086568e+01_EB, 4.818524e+00_EB, 1.625989e+00_EB, 4.963013e+00_EB, 9.091313e+00_EB, &
   2.862616e+00_EB, &
   1.807887e+01_EB, 1.240441e+01_EB, 1.649639e+00_EB, 3.870645e+00_EB, 4.710413e+00_EB, 5.816990e+00_EB, &
   2.011503e-01_EB/),(/7,8/))

gammad2_c7h8(1:7,17:24) = RESHAPE((/ &  ! 1055-1090 cm-1
   2.242814e+01_EB, 1.324109e+01_EB, 2.772538e+00_EB, 2.132624e+00_EB, 2.967247e+00_EB, 5.544644e+00_EB, &
   1.492977e-02_EB, &
   2.726543e+01_EB, 8.755811e+00_EB, 4.268738e+00_EB, 8.268829e+00_EB, 3.471302e+00_EB, 7.655547e+00_EB, &
   4.678557e-03_EB, &
   3.856386e+01_EB, 1.415885e+01_EB, 4.621129e+00_EB, 7.774568e+00_EB, 5.265771e-01_EB, 7.296652e+00_EB, &
   8.956780e-04_EB, &
   6.779154e+01_EB, 1.689515e+01_EB, 8.227216e+00_EB, 7.685408e+00_EB, 6.275187e+00_EB, 4.111618e+00_EB, &
   1.120381e-03_EB, &
   7.998412e+01_EB, 5.217367e+00_EB, 1.270480e+01_EB, 7.011278e+00_EB, 9.533137e+00_EB, 6.557224e+00_EB, &
   7.187912e-04_EB, &
   6.476932e+01_EB, 4.324490e+00_EB, 7.006941e+00_EB, 3.059077e+00_EB, 2.551623e+00_EB, 6.781168e+00_EB, &
   4.850850e-04_EB, &
   7.999697e+01_EB, 1.683983e+01_EB, 7.523811e+00_EB, 3.934915e+00_EB, 3.190601e+00_EB, 5.987392e+00_EB, &
   4.640575e-04_EB, &
   7.298422e+01_EB, 1.746059e+01_EB, 8.445975e+00_EB, 4.651083e+00_EB, 5.086497e+00_EB, 1.446891e+00_EB, &
   2.261182e-04_EB/),(/7,8/))

gammad2_c7h8(1:7,25:29) = RESHAPE((/ &  ! 1095-1175 cm-1
   5.656308e+01_EB, 1.932458e+01_EB, 5.566997e+00_EB, 2.210155e-01_EB, 4.019081e+00_EB, 7.861722e-02_EB, &
   3.404900e-03_EB, &
   2.190597e+01_EB, 1.019795e+01_EB, 2.244564e+00_EB, 1.625625e+00_EB, 5.393117e-01_EB, 1.054102e+00_EB, &
   7.309670e-05_EB, &
   4.172119e-01_EB, 6.066095e-01_EB, 3.794528e-01_EB, 2.201505e-01_EB, 8.146470e-02_EB, 1.686354e-01_EB, &
   2.521455e-03_EB, &
   4.904118e-04_EB, 8.240680e-03_EB, 2.319123e-01_EB, 4.751613e-04_EB, 5.109468e-04_EB, 5.032463e-04_EB, &
   4.422597e-04_EB, &
   3.682775e-04_EB, 3.700575e-04_EB, 3.724614e-04_EB, 3.699127e-04_EB, 3.709688e-04_EB, 3.696394e-04_EB, &
   3.680801e-04_EB/),(/7,5/))

!---------------------------------------------------------------------------
ALLOCATE(sd3_c7h8(n_temp_c7h8,17)) 

! band #3: 1275 cm-1 - 1675 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 1.8457 % 

sd3_c7h8(1:7,1:8) = RESHAPE((/ &  ! 1275-1450 cm-1
   6.100102e-04_EB, 7.780830e-03_EB, 1.134938e-03_EB, 8.737574e-04_EB, 6.027841e-03_EB, 4.138198e-03_EB, &
   9.056224e-04_EB, &
   7.403110e-04_EB, 3.544033e-02_EB, 3.023568e-02_EB, 1.508266e-02_EB, 2.813528e-02_EB, 2.043206e-02_EB, &
   4.350956e-01_EB, &
   3.178283e-03_EB, 3.702754e-02_EB, 4.100162e-02_EB, 2.423149e-02_EB, 4.352025e-02_EB, 3.192771e-02_EB, &
   9.676683e-01_EB, &
   5.605844e-02_EB, 6.834484e-02_EB, 8.335619e-02_EB, 6.360984e-02_EB, 6.991888e-02_EB, 5.926399e-02_EB, &
   6.248318e-02_EB, &
   2.335707e-01_EB, 2.104431e-01_EB, 2.014754e-01_EB, 1.664392e-01_EB, 1.529464e-01_EB, 1.145999e-01_EB, &
   4.880863e-02_EB, &
   2.297034e-01_EB, 1.952774e-01_EB, 2.009723e-01_EB, 1.554485e-01_EB, 1.472329e-01_EB, 9.989047e-02_EB, &
   5.114395e-02_EB, &
   2.466205e-01_EB, 2.254143e-01_EB, 2.146406e-01_EB, 1.846656e-01_EB, 1.917496e-01_EB, 1.694292e-01_EB, &
   9.885420e-02_EB, &
   6.338228e-01_EB, 4.862753e-01_EB, 4.470592e-01_EB, 3.954809e-01_EB, 3.379362e-01_EB, 2.524978e-01_EB, &
   2.245216e-01_EB/),(/7,8/))

sd3_c7h8(1:7,9:16) = RESHAPE((/ &  ! 1475-1650 cm-1
   8.091914e-01_EB, 6.844493e-01_EB, 6.417303e-01_EB, 6.044514e-01_EB, 5.830416e-01_EB, 4.886091e-01_EB, &
   3.460225e-01_EB, &
   1.985404e+00_EB, 1.434197e+00_EB, 1.197360e+00_EB, 1.090591e+00_EB, 8.747838e-01_EB, 5.725550e-01_EB, &
   3.573484e-01_EB, &
   5.171385e-01_EB, 4.486683e-01_EB, 4.336080e-01_EB, 3.554049e-01_EB, 3.143224e-01_EB, 2.347278e-01_EB, &
   1.742309e-01_EB, &
   2.372746e-01_EB, 2.197212e-01_EB, 2.373528e-01_EB, 1.855960e-01_EB, 1.760600e-01_EB, 1.693906e-01_EB, &
   1.742322e-01_EB, &
   2.366393e-01_EB, 3.786876e-01_EB, 4.105672e-01_EB, 2.381817e-01_EB, 2.272899e-01_EB, 2.664032e-01_EB, &
   3.863583e-01_EB, &
   8.572620e-01_EB, 1.527668e+00_EB, 1.063149e+00_EB, 7.295884e-01_EB, 4.745380e-01_EB, 3.602757e-01_EB, &
   2.354867e-01_EB, &
   5.567486e-01_EB, 5.472333e-01_EB, 5.016325e-01_EB, 3.512578e-01_EB, 1.717825e-01_EB, 9.711190e-02_EB, &
   2.718347e-02_EB, &
   2.625054e-02_EB, 5.072059e-03_EB, 1.694459e-02_EB, 9.001577e-03_EB, 5.938464e-04_EB, 6.681563e-04_EB, &
   5.196216e-04_EB/),(/7,8/))

sd3_c7h8(1:7,17:17) = RESHAPE((/ &  ! 1675-1675 cm-1
   4.621323e-04_EB, 4.638592e-04_EB, 4.665860e-04_EB, 4.636519e-04_EB, 4.646393e-04_EB, 4.633936e-04_EB, &
   4.615859e-04_EB/),(/7,1/))

!---------------------------------------------------------------------------
ALLOCATE(gammad3_c7h8(n_temp_c7h8,17)) 

! band #3: 1275 cm-1 - 1675 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 1.8457 % 

! print fine structure array gamma_d 

gammad3_c7h8(1:7,1:8) = RESHAPE((/ &  ! 1275-1450 cm-1
   4.795815e-04_EB, 2.124259e-03_EB, 1.024725e-03_EB, 7.095131e-04_EB, 6.311264e-03_EB, 2.305633e-03_EB, &
   5.045893e-04_EB, &
   5.695615e-04_EB, 1.692360e-02_EB, 7.400058e-03_EB, 1.063887e-01_EB, 2.392876e-01_EB, 3.837883e-04_EB, &
   3.113013e-05_EB, &
   2.535772e-03_EB, 3.802603e-02_EB, 2.314277e-02_EB, 4.430523e-01_EB, 1.112554e-01_EB, 1.504970e-01_EB, &
   4.347453e-05_EB, &
   7.174141e-02_EB, 1.160640e+00_EB, 2.191895e-02_EB, 4.952819e-01_EB, 4.273983e-01_EB, 3.451409e-01_EB, &
   2.592647e-01_EB, &
   8.276943e+00_EB, 3.502859e-01_EB, 5.113046e-02_EB, 2.027660e+00_EB, 2.707192e+00_EB, 3.833564e-01_EB, &
   3.464319e-01_EB, &
   1.015475e+01_EB, 4.674500e+00_EB, 4.689087e-02_EB, 3.954451e+00_EB, 3.031215e+00_EB, 1.038644e+00_EB, &
   6.356775e-02_EB, &
   1.948113e+01_EB, 7.058967e+00_EB, 9.857392e-01_EB, 6.992797e+00_EB, 4.529865e+00_EB, 2.876236e+00_EB, &
   1.523489e+00_EB, &
   7.999993e+01_EB, 4.801475e+01_EB, 2.297044e+00_EB, 2.307491e+01_EB, 2.255518e+01_EB, 6.921021e+00_EB, &
   3.793481e+00_EB/),(/7,8/))

gammad3_c7h8(1:7,9:16) = RESHAPE((/ &  ! 1475-1650 cm-1
   7.999994e+01_EB, 3.335220e+01_EB, 5.984703e-01_EB, 2.853393e+01_EB, 3.626523e+01_EB, 2.277723e+01_EB, &
   2.603661e+01_EB, &
   7.999988e+01_EB, 6.963104e+01_EB, 6.603136e+01_EB, 6.344400e+01_EB, 5.716242e+01_EB, 4.913311e+01_EB, &
   1.650821e+01_EB, &
   7.999998e+01_EB, 6.963102e+01_EB, 2.177188e+01_EB, 3.299642e+01_EB, 1.878956e+01_EB, 8.669014e+00_EB, &
   7.333480e+00_EB, &
   5.251899e+01_EB, 2.739552e+01_EB, 6.379395e+00_EB, 1.192897e+01_EB, 1.481504e+01_EB, 6.745394e+00_EB, &
   2.823356e+00_EB, &
   1.643612e-01_EB, 5.735340e-03_EB, 4.062201e-03_EB, 3.334872e-02_EB, 1.282097e+00_EB, 8.493361e-02_EB, &
   6.376486e-03_EB, &
   1.487048e-01_EB, 1.035601e-02_EB, 1.158524e-02_EB, 3.774427e-02_EB, 1.523506e-01_EB, 2.811809e-02_EB, &
   7.547214e-03_EB, &
   7.999361e+01_EB, 1.936478e-02_EB, 9.051617e-03_EB, 1.944160e-02_EB, 2.206877e-01_EB, 2.134674e-02_EB, &
   3.778034e-01_EB, &
   1.287431e-02_EB, 7.240508e-04_EB, 6.507482e-03_EB, 3.950013e-03_EB, 4.716179e-04_EB, 5.110565e-04_EB, &
   4.179497e-04_EB/),(/7,8/))

gammad3_c7h8(1:7,17:17) = RESHAPE((/ &  ! 1675-1675 cm-1
   3.682775e-04_EB, 3.700575e-04_EB, 3.724614e-04_EB, 3.699127e-04_EB, 3.709688e-04_EB, 3.696394e-04_EB, &
   3.680801e-04_EB/),(/7,1/))

!---------------------------------------------------------------------------
ALLOCATE(sd4_c7h8(n_temp_c7h8,18)) 

! band #4: 1650 cm-1 - 2075 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.43921 % 

sd4_c7h8(1:7,1:8) = RESHAPE((/ &  ! 1650-1825 cm-1
   1.608510e-02_EB, 8.538057e-04_EB, 2.868555e-03_EB, 8.340494e-04_EB, 6.015596e-04_EB, 4.952974e-04_EB, &
   5.998542e-04_EB, &
   4.848768e-02_EB, 9.249328e-03_EB, 2.165734e-02_EB, 1.578699e-02_EB, 7.483581e-04_EB, 5.699315e-04_EB, &
   7.083490e-04_EB, &
   3.748462e-02_EB, 1.617457e-02_EB, 2.578284e-02_EB, 1.643850e-02_EB, 1.017398e-02_EB, 2.017300e-03_EB, &
   7.500522e-03_EB, &
   1.132410e-01_EB, 2.967893e-01_EB, 6.916754e-02_EB, 5.538010e-02_EB, 3.122037e-02_EB, 1.282020e-02_EB, &
   6.409946e-02_EB, &
   4.946037e-01_EB, 4.749997e-02_EB, 4.638123e-02_EB, 4.364247e-02_EB, 3.296949e-02_EB, 3.250858e-02_EB, &
   6.216037e-02_EB, &
   9.754174e-01_EB, 6.827422e-01_EB, 5.097573e-01_EB, 4.113609e-01_EB, 9.436322e-02_EB, 1.578733e-01_EB, &
   1.005549e-01_EB, &
   5.760301e-01_EB, 1.058899e+00_EB, 5.428139e-01_EB, 2.925165e-01_EB, 1.220623e-01_EB, 2.559121e-01_EB, &
   8.881943e-02_EB, &
   1.187295e+00_EB, 1.013370e+00_EB, 7.226648e-01_EB, 5.498312e-01_EB, 8.139112e-02_EB, 4.364089e-01_EB, &
   8.132406e-02_EB/),(/7,8/))

sd4_c7h8(1:7,9:16) = RESHAPE((/ &  ! 1850-2025 cm-1
   8.856836e-01_EB, 9.186222e-01_EB, 4.794847e-01_EB, 4.024356e-01_EB, 1.306238e-01_EB, 3.063448e-01_EB, &
   9.935498e-02_EB, &
   1.034528e+00_EB, 7.215443e-01_EB, 5.929663e-01_EB, 4.319927e-01_EB, 8.863470e-02_EB, 3.544881e-01_EB, &
   5.660195e-02_EB, &
   1.495500e+00_EB, 3.838684e-01_EB, 1.268890e-01_EB, 4.632815e-02_EB, 1.148509e-01_EB, 6.472883e-02_EB, &
   8.665695e-02_EB, &
   5.568983e-01_EB, 5.685082e-01_EB, 5.071619e-01_EB, 3.647283e-01_EB, 1.319697e-01_EB, 1.328525e-01_EB, &
   8.945221e-02_EB, &
   2.939080e-01_EB, 3.996266e-01_EB, 2.679096e-01_EB, 1.696599e-01_EB, 1.249259e-01_EB, 8.650031e-02_EB, &
   7.539025e-02_EB, &
   4.677446e-01_EB, 2.881514e-02_EB, 6.793664e-02_EB, 3.565764e-02_EB, 1.939151e-02_EB, 1.637475e-02_EB, &
   3.310412e-02_EB, &
   3.406776e-02_EB, 1.092503e-02_EB, 1.198880e-02_EB, 1.579425e-02_EB, 1.184388e-02_EB, 1.945188e-03_EB, &
   2.385167e-02_EB, &
   1.429157e-02_EB, 8.354021e-04_EB, 4.409521e-03_EB, 4.730578e-03_EB, 4.878561e-03_EB, 8.332528e-04_EB, &
   2.664558e-02_EB/),(/7,8/))

sd4_c7h8(1:7,17:18) = RESHAPE((/ &  ! 2050-2075 cm-1
   5.862965e-04_EB, 5.938092e-04_EB, 6.747537e-04_EB, 5.948772e-04_EB, 6.040692e-04_EB, 6.075352e-04_EB, &
   3.980874e-03_EB, &
   4.621323e-04_EB, 4.638592e-04_EB, 4.665860e-04_EB, 4.636519e-04_EB, 4.646393e-04_EB, 4.633936e-04_EB, &
   4.615859e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(gammad4_c7h8(n_temp_c7h8,18)) 

! band #4: 1650 cm-1 - 2075 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 0.43921 % 

! print fine structure array gamma_d 

gammad4_c7h8(1:7,1:8) = RESHAPE((/ &  ! 1650-1825 cm-1
   9.073750e-01_EB, 7.155535e-04_EB, 1.417709e-01_EB, 6.862278e-04_EB, 4.770609e-04_EB, 3.962002e-04_EB, &
   4.599493e-04_EB, &
   4.897835e-01_EB, 2.228040e-03_EB, 1.889228e-03_EB, 3.037127e-03_EB, 5.957388e-04_EB, 4.412357e-04_EB, &
   5.439496e-04_EB, &
   8.222146e-01_EB, 1.967261e-01_EB, 1.176890e-01_EB, 1.993414e-01_EB, 3.392040e-01_EB, 1.285307e+00_EB, &
   7.330343e-02_EB, &
   1.711715e-02_EB, 2.387489e-04_EB, 4.930892e-03_EB, 6.589700e-03_EB, 1.627007e-01_EB, 8.239901e-03_EB, &
   5.086131e-03_EB, &
   2.105689e-04_EB, 4.813970e-03_EB, 5.785146e-03_EB, 2.367207e-02_EB, 2.418442e-01_EB, 7.693889e-03_EB, &
   3.096368e-02_EB, &
   3.797437e-04_EB, 2.010270e-04_EB, 4.063113e-04_EB, 7.930910e-04_EB, 6.260656e-02_EB, 3.447434e-03_EB, &
   2.052893e+00_EB, &
   5.974715e-03_EB, 9.625448e-04_EB, 1.385440e-03_EB, 4.460370e-03_EB, 7.510642e-02_EB, 1.409188e-03_EB, &
   5.077965e-01_EB, &
   2.598721e-04_EB, 9.700799e-05_EB, 1.875798e-04_EB, 3.395949e-04_EB, 1.967132e-02_EB, 5.386934e-04_EB, &
   7.215293e-01_EB/),(/7,8/))

gammad4_c7h8(1:7,9:16) = RESHAPE((/ &  ! 1850-2025 cm-1
   1.840578e-03_EB, 8.312411e-04_EB, 1.279188e-03_EB, 2.105186e-03_EB, 3.505265e-02_EB, 1.459376e-03_EB, &
   1.277784e+00_EB, &
   1.258299e-03_EB, 4.510848e-04_EB, 4.294071e-04_EB, 8.136257e-04_EB, 9.696057e-03_EB, 2.146063e-04_EB, &
   6.860048e-01_EB, &
   4.671598e-05_EB, 4.670991e-05_EB, 7.143144e-05_EB, 4.655956e-03_EB, 4.081375e-04_EB, 5.438877e-03_EB, &
   8.906804e-02_EB, &
   1.496812e-03_EB, 6.530123e-04_EB, 7.078867e-04_EB, 1.683327e-03_EB, 3.846642e-02_EB, 1.605013e-02_EB, &
   5.460498e-01_EB, &
   6.920106e-02_EB, 4.084630e-03_EB, 4.488361e-03_EB, 4.217654e-02_EB, 7.816064e-02_EB, 2.053699e-02_EB, &
   9.131605e-03_EB, &
   2.787181e-04_EB, 4.069568e-03_EB, 6.638116e-04_EB, 5.502709e-03_EB, 9.104577e-03_EB, 6.898146e-03_EB, &
   1.209388e-01_EB, &
   4.515662e-03_EB, 2.851322e-03_EB, 2.748099e-03_EB, 5.552875e-03_EB, 6.615247e-03_EB, 9.856806e-04_EB, &
   1.115308e-02_EB, &
   3.124492e-03_EB, 5.314268e-04_EB, 2.263215e-03_EB, 1.662162e-03_EB, 2.480485e-03_EB, 6.346151e-04_EB, &
   2.031025e-03_EB/),(/7,8/))

gammad4_c7h8(1:7,17:18) = RESHAPE((/ &  ! 2050-2075 cm-1
   4.587128e-04_EB, 4.524884e-04_EB, 5.170999e-04_EB, 4.681410e-04_EB, 4.819639e-04_EB, 4.673877e-04_EB, &
   8.021840e-04_EB, &
   3.682775e-04_EB, 3.700575e-04_EB, 3.724614e-04_EB, 3.699127e-04_EB, 3.709688e-04_EB, 3.696394e-04_EB, &
   3.680801e-04_EB/),(/7,2/))

!---------------------------------------------------------------------------
ALLOCATE(sd5_c7h8(n_temp_c7h8,23)) 

! band #5: 2675 cm-1 - 3225 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.2868 % 

sd5_c7h8(1:7,1:8) = RESHAPE((/ &  ! 2675-2850 cm-1
   5.916860e-04_EB, 6.494496e-04_EB, 7.055877e-04_EB, 6.603478e-04_EB, 6.674805e-04_EB, 5.589319e-04_EB, &
   6.137597e-04_EB, &
   7.500569e-04_EB, 7.074450e-04_EB, 9.098169e-04_EB, 7.561953e-04_EB, 1.607193e-03_EB, 1.018018e-02_EB, &
   7.780318e-03_EB, &
   7.029259e-04_EB, 8.479757e-04_EB, 9.039741e-04_EB, 8.123582e-04_EB, 2.347565e-02_EB, 3.011753e-02_EB, &
   2.717306e-02_EB, &
   8.267165e-04_EB, 8.801454e-04_EB, 2.940361e-03_EB, 1.990166e-03_EB, 4.660083e-02_EB, 4.429277e-02_EB, &
   3.863177e-02_EB, &
   8.156025e-04_EB, 1.506400e-03_EB, 6.730104e-03_EB, 3.740574e-03_EB, 2.586356e-02_EB, 3.987065e-02_EB, &
   4.205663e-02_EB, &
   9.809286e-03_EB, 1.335514e-02_EB, 2.488364e-02_EB, 1.572642e-02_EB, 4.731877e-02_EB, 6.718605e-02_EB, &
   6.938131e-02_EB, &
   5.850971e-02_EB, 6.588341e-02_EB, 7.434265e-02_EB, 6.744072e-02_EB, 1.018290e-01_EB, 1.263017e-01_EB, &
   1.302993e-01_EB, &
   2.663902e-01_EB, 2.671000e-01_EB, 2.879771e-01_EB, 2.569076e-01_EB, 2.774278e-01_EB, 2.512113e-01_EB, &
   2.305380e-01_EB/),(/7,8/))

sd5_c7h8(1:7,9:16) = RESHAPE((/ &  ! 2875-3050 cm-1
   1.192316e+00_EB, 9.205767e-01_EB, 8.299791e-01_EB, 7.348094e-01_EB, 6.292028e-01_EB, 4.584459e-01_EB, &
   3.314507e-01_EB, &
   1.111204e+00_EB, 9.394978e-01_EB, 8.641911e-01_EB, 8.038659e-01_EB, 7.457709e-01_EB, 6.165272e-01_EB, &
   4.772364e-01_EB, &
   2.020089e+00_EB, 1.631718e+00_EB, 1.548727e+00_EB, 1.354242e+00_EB, 1.167657e+00_EB, 8.526752e-01_EB, &
   6.021600e-01_EB, &
   1.932741e+00_EB, 1.403246e+00_EB, 1.298233e+00_EB, 1.078566e+00_EB, 8.872613e-01_EB, 5.993978e-01_EB, &
   4.326418e-01_EB, &
   9.805520e-01_EB, 7.107386e-01_EB, 6.778887e-01_EB, 5.479328e-01_EB, 4.804311e-01_EB, 3.935294e-01_EB, &
   3.642029e-01_EB, &
   1.226852e+00_EB, 9.501327e-01_EB, 8.230054e-01_EB, 8.274332e-01_EB, 8.323208e-01_EB, 7.816587e-01_EB, &
   6.650895e-01_EB, &
   3.262768e+00_EB, 2.723650e+00_EB, 2.296188e+00_EB, 2.270140e+00_EB, 1.935929e+00_EB, 1.369446e+00_EB, &
   9.666500e-01_EB, &
   3.556303e+00_EB, 2.816675e+00_EB, 2.390948e+00_EB, 2.282893e+00_EB, 1.948012e+00_EB, 1.461546e+00_EB, &
   1.092304e+00_EB/),(/7,8/))

sd5_c7h8(1:7,17:23) = RESHAPE((/ &  ! 3075-3225 cm-1
   2.694239e+00_EB, 2.217181e+00_EB, 1.960010e+00_EB, 1.867801e+00_EB, 1.644870e+00_EB, 1.225906e+00_EB, &
   8.559281e-01_EB, &
   1.505365e+00_EB, 1.180006e+00_EB, 1.047132e+00_EB, 9.563041e-01_EB, 8.157730e-01_EB, 5.594937e-01_EB, &
   3.588645e-01_EB, &
   3.952204e-01_EB, 3.072414e-01_EB, 2.787226e-01_EB, 2.333774e-01_EB, 1.986791e-01_EB, 1.344281e-01_EB, &
   8.253592e-02_EB, &
   2.640387e-02_EB, 3.373476e-02_EB, 3.927558e-02_EB, 2.200042e-02_EB, 3.743437e-02_EB, 3.283105e-02_EB, &
   1.723376e-02_EB, &
   2.071719e-02_EB, 1.345537e-02_EB, 1.924395e-02_EB, 6.019999e-03_EB, 1.908321e-02_EB, 1.111107e-02_EB, &
   6.304850e-04_EB, &
   1.966078e-03_EB, 1.825102e-03_EB, 3.977009e-03_EB, 6.125614e-04_EB, 1.893366e-03_EB, 2.342482e-03_EB, &
   5.793976e-04_EB, &
   4.621323e-04_EB, 4.638592e-04_EB, 4.665860e-04_EB, 4.636519e-04_EB, 4.646393e-04_EB, 4.633936e-04_EB, &
   4.615859e-04_EB/),(/7,7/))

!---------------------------------------------------------------------------
ALLOCATE(gammad5_c7h8(n_temp_c7h8,23)) 

! band #5: 2675 cm-1 - 3225 cm-1 

! snb fit with malkmus model 

! error associated with malkmus fit: 
! on transmissivity: 2.2868 % 

! print fine structure array gamma_d 

gammad5_c7h8(1:7,1:8) = RESHAPE((/ &  ! 2675-2850 cm-1
   4.498679e-04_EB, 5.206905e-04_EB, 5.314804e-04_EB, 5.225909e-04_EB, 5.382880e-04_EB, 4.354686e-04_EB, &
   4.954396e-04_EB, &
   5.733660e-04_EB, 5.664218e-04_EB, 6.760255e-04_EB, 5.720172e-04_EB, 1.214045e-03_EB, 6.728005e-03_EB, &
   6.901498e-02_EB, &
   5.501549e-04_EB, 6.422611e-04_EB, 6.751496e-04_EB, 6.170879e-04_EB, 2.505333e-02_EB, 2.929665e-02_EB, &
   2.183072e-01_EB, &
   6.226316e-04_EB, 6.634878e-04_EB, 2.204441e-03_EB, 1.315450e-03_EB, 1.427511e-01_EB, 5.684139e-02_EB, &
   1.991618e-01_EB, &
   6.173941e-04_EB, 8.607660e-04_EB, 4.084909e-03_EB, 1.485988e-03_EB, 6.217446e-02_EB, 1.696125e-01_EB, &
   6.362403e-01_EB, &
   2.271432e-01_EB, 2.575315e-02_EB, 1.742774e-02_EB, 4.872962e-03_EB, 1.222993e-01_EB, 2.494679e-01_EB, &
   4.960213e-01_EB, &
   3.499358e-01_EB, 9.426016e-01_EB, 5.241002e-02_EB, 4.920380e-02_EB, 6.407093e-01_EB, 1.208700e+00_EB, &
   3.229169e+00_EB, &
   2.056851e+01_EB, 2.155787e+00_EB, 1.241796e-01_EB, 3.443147e-01_EB, 8.158180e+00_EB, 8.235706e+00_EB, &
   4.544885e+00_EB/),(/7,8/))

gammad5_c7h8(1:7,9:16) = RESHAPE((/ &  ! 2875-3050 cm-1
   7.999651e+01_EB, 4.892463e+01_EB, 2.269762e+00_EB, 1.699238e+01_EB, 3.883167e+01_EB, 2.126359e+01_EB, &
   1.128396e+01_EB, &
   7.999972e+01_EB, 6.309208e+01_EB, 3.234965e+00_EB, 1.718937e+01_EB, 4.942497e+01_EB, 3.563188e+01_EB, &
   2.569318e+01_EB, &
   7.999990e+01_EB, 6.940301e+01_EB, 2.588169e+00_EB, 6.006982e+00_EB, 5.718259e+01_EB, 4.914102e+01_EB, &
   3.825054e+01_EB, &
   7.999998e+01_EB, 6.945636e+01_EB, 5.656023e+01_EB, 4.602265e+01_EB, 5.718979e+01_EB, 4.365360e+01_EB, &
   2.767624e+01_EB, &
   7.999997e+01_EB, 6.310533e+01_EB, 2.112740e+01_EB, 2.715276e+01_EB, 2.826351e+01_EB, 1.692961e+01_EB, &
   1.913465e+01_EB, &
   7.999883e+01_EB, 6.945732e+01_EB, 3.458543e+00_EB, 4.892432e+01_EB, 5.719058e+01_EB, 4.910638e+01_EB, &
   3.932336e+01_EB, &
   7.999636e+01_EB, 6.932335e+01_EB, 2.959489e+01_EB, 2.214939e+01_EB, 5.718273e+01_EB, 4.914071e+01_EB, &
   4.383529e+01_EB, &
   7.999930e+01_EB, 6.963014e+01_EB, 6.599550e+01_EB, 5.266039e+01_EB, 5.718071e+01_EB, 4.910109e+01_EB, &
   4.383604e+01_EB/),(/7,8/))

gammad5_c7h8(1:7,17:23) = RESHAPE((/ &  ! 3075-3225 cm-1
   7.999949e+01_EB, 6.961299e+01_EB, 6.604063e+01_EB, 2.849107e+01_EB, 5.719148e+01_EB, 4.914326e+01_EB, &
   4.383622e+01_EB, &
   7.999993e+01_EB, 6.944137e+01_EB, 6.603883e+01_EB, 3.834261e+01_EB, 5.713287e+01_EB, 3.198295e+01_EB, &
   1.533847e+01_EB, &
   5.103942e+01_EB, 7.534768e+00_EB, 6.064864e+00_EB, 1.229512e+00_EB, 4.034341e+00_EB, 1.526492e+00_EB, &
   9.982477e-01_EB, &
   4.966937e-01_EB, 3.630636e-01_EB, 1.127309e-01_EB, 1.271562e-02_EB, 2.831161e-01_EB, 2.892147e-01_EB, &
   8.723984e-02_EB, &
   2.678214e-01_EB, 1.122623e-01_EB, 2.684864e-02_EB, 1.687467e-03_EB, 8.582904e-03_EB, 1.126356e-01_EB, &
   4.870184e-04_EB, &
   5.319443e-03_EB, 2.133031e-03_EB, 2.657648e-03_EB, 4.684001e-04_EB, 1.112851e-03_EB, 5.919668e-03_EB, &
   4.552670e-04_EB, &
   3.682775e-04_EB, 3.700575e-04_EB, 3.724614e-04_EB, 3.699127e-04_EB, 3.709688e-04_EB, 3.696394e-04_EB, &
   3.680801e-04_EB/),(/7,7/))

!--------------------------------old radcal data----------------------------------------
!  initialize sd array

! temp,k= 300     600      1000      1500      2000      2500       

sd(1:6,1:8) = RESHAPE ((/  & ! 50-225
.950e+00_EB, .103e+00_EB, .420e-01_EB, .114e-01_EB, .450e-02_EB, .300e-02_EB, &
.208e+01_EB, .365e+00_EB, .113e+00_EB, .375e-01_EB, .195e-01_EB, .134e-01_EB, &
.368e+01_EB, .990e+00_EB, .300e+00_EB, .104e+00_EB, .577e-01_EB, .365e-01_EB,  & 
.650e+01_EB, .201e+01_EB, .650e+00_EB, .214e+00_EB, .128e+00_EB, .845e-01_EB, &   
.825e+01_EB, .325e+01_EB, .121e+01_EB, .415e+00_EB, .260e+00_EB, .168e+00_EB, &  
.870e+01_EB, .452e+01_EB, .189e+01_EB, .765e+00_EB, .450e+00_EB, .289e+00_EB, &  
.810e+01_EB, .540e+01_EB, .261e+01_EB, .126e+01_EB, .695e+00_EB, .460e+00_EB, &  
.682e+01_EB, .600e+01_EB, .337e+01_EB, .179e+01_EB, .101e+01_EB, .679e+00_EB/), &
(/6,8/))
sd(1:6,9:16) = RESHAPE ((/ &  ! 250-425
.493e+01_EB, .622e+01_EB, .407e+01_EB, .230e+01_EB, .135e+01_EB, .935e+00_EB,  &   
.316e+01_EB, .592e+01_EB, .456e+01_EB, .281e+01_EB, .172e+01_EB, .122e+01_EB,  &  
.199e+01_EB, .528e+01_EB, .479e+01_EB, .328e+01_EB, .213e+01_EB, .149e+01_EB,  &   
.113e+01_EB, .450e+01_EB, .484e+01_EB, .361e+01_EB, .249e+01_EB, .179e+01_EB,  &   
.585e+00_EB, .370e+01_EB, .471e+01_EB, .383e+01_EB, .284e+01_EB, .208e+01_EB, &    
.293e+00_EB, .289e+01_EB, .443e+01_EB, .394e+01_EB, .312e+01_EB, .237e+01_EB, &    
.138e+00_EB, .205e+01_EB, .400e+01_EB, .396e+01_EB, .330e+01_EB, .260e+01_EB, &    
.620e-01_EB, .143e+01_EB, .347e+01_EB, .388e+01_EB, .341e+01_EB, .280e+01_EB/), &  
(/6,8/))
sd(1:6,17:24) = RESHAPE ((/ & ! 450-625
.255e-01_EB, .950e+00_EB, .292e+01_EB, .370e+01_EB, .345e+01_EB, .295e+01_EB,   &  
.940e-02_EB, .610e+00_EB, .236e+01_EB, .343e+01_EB, .342e+01_EB, .304e+01_EB,  &   
.340e-02_EB, .386e+00_EB, .188e+01_EB, .310e+01_EB, .334e+01_EB, .309e+01_EB,   &  
.105e-02_EB, .236e+00_EB, .145e+01_EB, .274e+01_EB, .319e+01_EB, .307e+01_EB,   &  
.350e-03_EB, .144e+00_EB, .110e+01_EB, .238e+01_EB, .300e+01_EB, .301e+01_EB,   & 
.126e-03_EB, .820e-01_EB, .818e+00_EB, .204e+01_EB, .276e+01_EB, .289e+01_EB,   &  
.430e-04_EB, .445e-01_EB, .598e+00_EB, .174e+01_EB, .248e+01_EB, .275e+01_EB,  &  
.150e-04_EB, .242e-01_EB, .427e+00_EB, .145e+01_EB, .222e+01_EB, .260e+01_EB/), & 
(/6,8/))
sd(1:6,25:32) = RESHAPE ((/  &! 650-825
.510e-05_EB, .127e-01_EB, .294e+00_EB, .118e+01_EB, .195e+01_EB, .241e+01_EB,   &  
.170e-05_EB, .630e-02_EB, .200e+00_EB, .950e+00_EB, .169e+01_EB, .221e+01_EB,   &  
.570e-06_EB, .300e-02_EB, .134e+00_EB, .748e+00_EB, .146e+01_EB, .200e+01_EB,   &  
.195e-06_EB, .140e-02_EB, .902e-01_EB, .580e+00_EB, .124e+01_EB, .178e+01_EB,   &  
.680e-07_EB, .620e-03_EB, .590e-01_EB, .443e+00_EB, .103e+01_EB, .156e+01_EB,   &  
.385e-07_EB, .275e-03_EB, .450e-01_EB, .330e+00_EB, .845e+00_EB, .136e+01_EB,   &  
.670e-07_EB, .113e-03_EB, .355e-01_EB, .242e+00_EB, .695e+00_EB, .117e+01_EB,   &  
.113e-06_EB, .500e-04_EB, .289e-01_EB, .174e+00_EB, .560e+00_EB, .100e+01_EB/),  & 
(/6,8/))
sd(1:6,33:40) = RESHAPE ((/ & ! 850-1025
.195e-06_EB, .230e-04_EB, .245e-01_EB, .123e+00_EB, .450e+00_EB, .855e+00_EB,   &  
.328e-06_EB, .103e-04_EB, .214e-01_EB, .100e+00_EB, .357e+00_EB, .718e+00_EB,  &   
.560e-06_EB, .460e-05_EB, .189e-01_EB, .830e-01_EB, .278e+00_EB, .595e+00_EB,  &   
.950e-06_EB, .205e-05_EB, .174e-01_EB, .730e-01_EB, .239e+00_EB, .492e+00_EB,  &   
.160e-05_EB, .140e-05_EB, .166e-01_EB, .665e-01_EB, .211e+00_EB, .405e+00_EB, &    
.275e-05_EB, .350e-05_EB, .165e-01_EB, .630e-01_EB, .195e+00_EB, .352e+00_EB, &    
.470e-05_EB, .850e-05_EB, .167e-01_EB, .620e-01_EB, .190e+00_EB, .312e+00_EB, &   
.810e-05_EB, .215e-04_EB, .175e-01_EB, .630e-01_EB, .191e+00_EB, .289e+00_EB/), &  
(/6,8/))
sd(1:6,41:48) = RESHAPE ((/ & ! 1050-1225
.136e-04_EB, .570e-04_EB, .188e-01_EB, .675e-01_EB, .194e+00_EB, .281e+00_EB,   &  
.235e-04_EB, .150e-03_EB, .208e-01_EB, .745e-01_EB, .202e+00_EB, .283e+00_EB,   &  
.400e-04_EB, .380e-03_EB, .233e-01_EB, .865e-01_EB, .223e+00_EB, .314e+00_EB,  &   
.680e-04_EB, .950e-03_EB, .268e-01_EB, .122e+00_EB, .260e+00_EB, .380e+00_EB,  &   
.120e-03_EB, .245e-02_EB, .343e-01_EB, .176e+00_EB, .328e+00_EB, .461e+00_EB,  &   
.200e-03_EB, .620e-02_EB, .638e-01_EB, .251e+00_EB, .411e+00_EB, .511e+00_EB,  &   
.365e-03_EB, .140e-01_EB, .107e+00_EB, .330e+00_EB, .458e+00_EB, .542e+00_EB,  &   
.680e-03_EB, .330e-01_EB, .166e+00_EB, .405e+00_EB, .487e+00_EB, .571e+00_EB/),&   
(/6,8/))
sd(1:6,49:56) = RESHAPE ((/ & ! 1250-1425
.130e-02_EB, .635e-01_EB, .244e+00_EB, .459e+00_EB, .535e+00_EB, .557e+00_EB,  &   
.250e-02_EB, .123e+00_EB, .341e+00_EB, .477e+00_EB, .502e+00_EB, .562e+00_EB,  &   
.500e-02_EB, .212e+00_EB, .407e+00_EB, .547e+00_EB, .531e+00_EB, .514e+00_EB,  &   
.103e-01_EB, .285e+00_EB, .489e+00_EB, .592e+00_EB, .497e+00_EB, .486e+00_EB,  &   
.219e-01_EB, .328e+00_EB, .491e+00_EB, .558e+00_EB, .489e+00_EB, .485e+00_EB,  &   
.485e-01_EB, .345e+00_EB, .505e+00_EB, .521e+00_EB, .477e+00_EB, .484e+00_EB,  &   
.114e+00_EB, .361e+00_EB, .538e+00_EB, .563e+00_EB, .503e+00_EB, .502e+00_EB,  &   
.249e+00_EB, .460e+00_EB, .621e+00_EB, .624e+00_EB, .538e+00_EB, .538e+00_EB/), &  
(/6,8/))
sd(1:6,57:64) = RESHAPE ((/ & ! 1450-1625
.397e+00_EB, .569e+00_EB, .749e+00_EB, .768e+00_EB, .581e+00_EB, .565e+00_EB,  &   
.418e+00_EB, .627e+00_EB, .824e+00_EB, .849e+00_EB, .640e+00_EB, .594e+00_EB,  &   
.108e+01_EB, .125e+01_EB, .113e+01_EB, .940e+00_EB, .807e+00_EB, .663e+00_EB,  &   
.165e+01_EB, .155e+01_EB, .118e+01_EB, .670e+00_EB, .562e+00_EB, .483e+00_EB,  &   
.142e+01_EB, .675e+00_EB, .557e+00_EB, .349e+00_EB, .276e+00_EB, .263e+00_EB,  &   
.451e+00_EB, .202e+00_EB, .132e+00_EB, .118e+00_EB, .134e+00_EB, .156e+00_EB,  &   
.603e-01_EB, .538e-01_EB, .863e-01_EB, .112e+00_EB, .120e+00_EB, .125e+00_EB,  &   
.501e+00_EB, .252e+00_EB, .118e+00_EB, .112e+00_EB, .131e+00_EB, .140e+00_EB/), &  
(/6,8/))
sd(1:6,65:72) = RESHAPE ((/ & ! 1650-1825
.730e+00_EB, .430e+00_EB, .237e+00_EB, .191e+00_EB, .171e+00_EB, .170e+00_EB,   &  
.149e+01_EB, .506e+00_EB, .294e+00_EB, .238e+00_EB, .210e+00_EB, .201e+00_EB,  &   
.100e+01_EB, .553e+00_EB, .434e+00_EB, .340e+00_EB, .260e+00_EB, .220e+00_EB,   &  
.802e+00_EB, .658e+00_EB, .528e+00_EB, .411e+00_EB, .300e+00_EB, .240e+00_EB,    & 
.580e+00_EB, .527e+00_EB, .460e+00_EB, .378e+00_EB, .322e+00_EB, .283e+00_EB,  &   
.330e+00_EB, .403e+00_EB, .430e+00_EB, .356e+00_EB, .318e+00_EB, .270e+00_EB,  &   
.250e+00_EB, .393e+00_EB, .405e+00_EB, .342e+00_EB, .301e+00_EB, .275e+00_EB,  &   
.147e+00_EB, .249e+00_EB, .313e+00_EB, .318e+00_EB, .291e+00_EB, .268e+00_EB/), &  
(/6,8/))
sd(1:6,73:80) = RESHAPE ((/ & ! 1850-2025
.910e-01_EB, .252e+00_EB, .298e+00_EB, .295e+00_EB, .269e+00_EB, .253e+00_EB,   &  
.580e-01_EB, .158e+00_EB, .214e+00_EB, .244e+00_EB, .244e+00_EB, .245e+00_EB,   &  
.370e-01_EB, .113e+00_EB, .184e+00_EB, .218e+00_EB, .214e+00_EB, .218e+00_EB,   &  
.244e-01_EB, .118e+00_EB, .156e+00_EB, .188e+00_EB, .195e+00_EB, .200e+00_EB,   & 
.162e-01_EB, .606e-01_EB, .976e-01_EB, .141e+00_EB, .166e+00_EB, .179e+00_EB,   & 
.112e-01_EB, .425e-01_EB, .903e-01_EB, .133e+00_EB, .148e+00_EB, .156e+00_EB,   & 
.780e-02_EB, .400e-01_EB, .765e-01_EB, .112e+00_EB, .129e+00_EB, .137e+00_EB,   & 
.540e-02_EB, .352e-01_EB, .647e-01_EB, .876e-01_EB, .110e+00_EB, .118e+00_EB/),  &
(/6,8/))
sd(1:6,81:88) = RESHAPE ((/ & ! 2050-2225
.380e-02_EB, .252e-01_EB, .507e-01_EB, .705e-01_EB, .888e-01_EB, .100e+00_EB,  &   
.260e-02_EB, .179e-01_EB, .377e-01_EB, .546e-01_EB, .724e-01_EB, .828e-01_EB,  &   
.180e-02_EB, .123e-01_EB, .294e-01_EB, .443e-01_EB, .608e-01_EB, .686e-01_EB,  &   
.127e-02_EB, .850e-02_EB, .212e-01_EB, .378e-01_EB, .579e-01_EB, .640e-01_EB,   &  
.880e-03_EB, .680e-02_EB, .152e-01_EB, .275e-01_EB, .449e-01_EB, .521e-01_EB,  &   
.620e-02_EB, .400e-02_EB, .107e-01_EB, .214e-01_EB, .374e-01_EB, .453e-01_EB,  &   
.480e-03_EB, .298e-02_EB, .931e-02_EB, .189e-01_EB, .329e-01_EB, .403e-01_EB,  &   
.405e-03_EB, .175e-02_EB, .696e-02_EB, .152e-01_EB, .295e-01_EB, .365e-01_EB/), &  
(/6,8/))
sd(1:6,89:96) = RESHAPE ((/ & ! 2250-2425 
.321e-03_EB, .120e-02_EB, .452e-02_EB, .101e-01_EB, .252e-01_EB, .331e-01_EB, &    
.229e-03_EB, .721e-03_EB, .364e-02_EB, .930e-02_EB, .225e-01_EB, .305e-01_EB,  &   
.195e-03_EB, .544e-03_EB, .318e-02_EB, .750e-02_EB, .202e-01_EB, .284e-01_EB, &    
.154e-03_EB, .375e-03_EB, .185e-02_EB, .603e-02_EB, .175e-01_EB, .269e-01_EB, &    
.101e-03_EB, .263e-03_EB, .119e-02_EB, .480e-02_EB, .156e-01_EB, .253e-01_EB, &    
.852e-04_EB, .185e-03_EB, .909e-03_EB, .360e-02_EB, .133e-01_EB, .241e-01_EB,  &   
.763e-04_EB, .137e-03_EB, .711e-03_EB, .316e-02_EB, .122e-01_EB, .237e-01_EB,  &   
.615e-04_EB, .126e-03_EB, .610e-03_EB, .257e-02_EB, .101e-01_EB, .218e-01_EB/), &  
(/6,8/))
sd(1:6,97:104) = RESHAPE ((/ & ! 2450-2625 
.480e-04_EB, .113e-03_EB, .518e-03_EB, .201e-02_EB, .920e-02_EB, .200e-01_EB, &    
.372e-04_EB, .106e-03_EB, .435e-03_EB, .168e-02_EB, .785e-02_EB, .183e-01_EB,  &   
.355e-04_EB, .101e-03_EB, .376e-03_EB, .168e-02_EB, .669e-02_EB, .166e-01_EB, &    
.358e-04_EB, .990e-04_EB, .366e-03_EB, .167e-02_EB, .651e-02_EB, .156e-01_EB,  &   
.389e-04_EB, .102e-03_EB, .376e-03_EB, .167e-02_EB, .641e-02_EB, .152e-01_EB,  &   
.422e-04_EB, .106e-03_EB, .373e-03_EB, .168e-02_EB, .656e-02_EB, .150e-01_EB,  &   
.521e-04_EB, .111e-03_EB, .371e-03_EB, .170e-02_EB, .673e-02_EB, .152e-01_EB,  &   
.646e-04_EB, .121e-03_EB, .384e-03_EB, .179e-02_EB, .798e-02_EB, .179e-01_EB/), &  
(/6,8/)) 
sd(1:6,105:112) = RESHAPE ((/ & ! 2650-2825 
.742e-04_EB, .129e-03_EB, .479e-03_EB, .201e-02_EB, .788e-02_EB, .175e-01_EB,   &  
.953e-04_EB, .165e-03_EB, .544e-03_EB, .249e-02_EB, .945e-02_EB, .204e-01_EB,   &  
.101e-03_EB, .190e-03_EB, .761e-03_EB, .324e-02_EB, .106e-01_EB, .231e-01_EB,   &  
.147e-03_EB, .272e-03_EB, .892e-03_EB, .441e-02_EB, .125e-01_EB, .257e-01_EB,    & 
.195e-03_EB, .326e-03_EB, .100e-02_EB, .499e-02_EB, .147e-01_EB, .295e-01_EB,    & 
.261e-03_EB, .421e-03_EB, .145e-02_EB, .568e-02_EB, .161e-01_EB, .306e-01_EB,   &  
.305e-03_EB, .515e-03_EB, .195e-02_EB, .754e-02_EB, .185e-01_EB, .363e-01_EB,   &  
.362e-03_EB, .645e-03_EB, .237e-02_EB, .830e-02_EB, .205e-01_EB, .373e-01_EB/), &  
(/6,8/))
sd(1:6,113:120) = RESHAPE ((/ & ! 2850-3025 
.507e-03_EB, .850e-03_EB, .274e-02_EB, .888e-02_EB, .234e-01_EB, .431e-01_EB,     &
.799e-03_EB, .118e-02_EB, .322e-02_EB, .110e-01_EB, .262e-01_EB, .451e-01_EB,     &
.935e-03_EB, .160e-02_EB, .386e-02_EB, .126e-01_EB, .292e-01_EB, .530e-01_EB,     &
.108e-02_EB, .231e-02_EB, .451e-02_EB, .140e-01_EB, .306e-01_EB, .536e-01_EB,     &
.192e-02_EB, .271e-02_EB, .563e-02_EB, .159e-01_EB, .357e-01_EB, .629e-01_EB,     &
.263e-02_EB, .300e-02_EB, .625e-02_EB, .179e-01_EB, .385e-01_EB, .666e-01_EB,    & 
.295e-02_EB, .330e-02_EB, .701e-02_EB, .203e-01_EB, .460e-01_EB, .782e-01_EB,   &  
.310e-02_EB, .370e-02_EB, .846e-02_EB, .220e-01_EB, .519e-01_EB, .889e-01_EB/),  & 
(/6,8/))
sd(1:6,121:128) = RESHAPE ((/ & ! 3050-3225
.340e-02_EB, .400e-02_EB, .969e-02_EB, .279e-01_EB, .662e-01_EB, .109e+00_EB,     &
.730e-02_EB, .450e-02_EB, .111e-01_EB, .272e-01_EB, .676e-01_EB, .109e+00_EB,     &
.900e-02_EB, .480e-02_EB, .137e-01_EB, .372e-01_EB, .864e-01_EB, .133e+00_EB,     &
.100e-02_EB, .510e-02_EB, .162e-01_EB, .471e-01_EB, .100e+00_EB, .142e+00_EB,     &
.640e-03_EB, .550e-02_EB, .205e-01_EB, .530e-01_EB, .122e+00_EB, .168e+00_EB,     &
.160e-02_EB, .600e-02_EB, .247e-01_EB, .633e-01_EB, .135e+00_EB, .177e+00_EB,    & 
.330e-02_EB, .700e-02_EB, .283e-01_EB, .770e-01_EB, .153e+00_EB, .185e+00_EB,   &  
.410e-02_EB, .860e-02_EB, .376e-01_EB, .914e-01_EB, .166e+00_EB, .206e+00_EB/),&   
(/6,8/))
sd(1:6,129:136) = RESHAPE ((/ & ! 3250-3425
.410e-02_EB, .103e-01_EB, .514e-01_EB, .117e+00_EB, .194e+00_EB, .228e+00_EB,     &
.290e-02_EB, .129e-01_EB, .664e-01_EB, .147e+00_EB, .220e+00_EB, .254e+00_EB,     &
.220e-02_EB, .161e-01_EB, .834e-01_EB, .171e+00_EB, .237e+00_EB, .263e+00_EB,     &
.220e-02_EB, .212e-01_EB, .103e+00_EB, .201e+00_EB, .268e+00_EB, .283e+00_EB,     &
.250e-02_EB, .285e-01_EB, .135e+00_EB, .240e+00_EB, .295e+00_EB, .295e+00_EB,     &
.310e-02_EB, .385e-01_EB, .169e+00_EB, .272e+00_EB, .312e+00_EB, .301e+00_EB,     &
.420e-02_EB, .540e-01_EB, .214e+00_EB, .309e+00_EB, .329e+00_EB, .307e+00_EB,     &
.600e-02_EB, .770e-01_EB, .267e+00_EB, .343e+00_EB, .332e+00_EB, .314e+00_EB/),  & 
(/6,8/))
sd(1:6,137:144) = RESHAPE ((/ & ! 3450-3625
.940e-02_EB, .117e+00_EB, .333e+00_EB, .372e+00_EB, .344e+00_EB, .303e+00_EB,     &
.165e-01_EB, .173e+00_EB, .365e+00_EB, .385e+00_EB, .353e+00_EB, .300e+00_EB,     &
.360e-01_EB, .258e+00_EB, .438e+00_EB, .393e+00_EB, .315e+00_EB, .288e+00_EB,     &
.720e-01_EB, .375e+00_EB, .510e+00_EB, .409e+00_EB, .294e+00_EB, .271e+00_EB,     &
.133e+00_EB, .401e+00_EB, .499e+00_EB, .390e+00_EB, .281e+00_EB, .257e+00_EB,     &
.215e+00_EB, .500e+00_EB, .443e+00_EB, .341e+00_EB, .254e+00_EB, .230e+00_EB,     &
.318e+00_EB, .450e+00_EB, .346e+00_EB, .286e+00_EB, .245e+00_EB, .219e+00_EB,     &
.442e+00_EB, .400e+00_EB, .354e+00_EB, .279e+00_EB, .233e+00_EB, .216e+00_EB/),   &
(/6,8/))
sd(1:6,145:152) = RESHAPE ((/ & ! 3650-3825
.473e+00_EB, .405e+00_EB, .347e+00_EB, .281e+00_EB, .238e+00_EB, .219e+00_EB,     &
.568e+00_EB, .501e+00_EB, .423e+00_EB, .315e+00_EB, .243e+00_EB, .218e+00_EB,     &
.690e+00_EB, .708e+00_EB, .673e+00_EB, .432e+00_EB, .268e+00_EB, .189e+00_EB,     &
.617e+00_EB, .831e+00_EB, .566e+00_EB, .320e+00_EB, .194e+00_EB, .123e+00_EB,     &
.181e+01_EB, .520e+00_EB, .200e+00_EB, .131e+00_EB, .124e+00_EB, .107e+00_EB,     &
.136e+00_EB, .124e+00_EB, .120e+00_EB, .119e+00_EB, .115e+00_EB, .115e+00_EB,     &
.455e+00_EB, .298e+00_EB, .167e+00_EB, .129e+00_EB, .123e+00_EB, .112e+00_EB,     &
.760e+00_EB, .503e+00_EB, .242e+00_EB, .154e+00_EB, .129e+00_EB, .127e+00_EB/),   &
(/6,8/))
sd(1:6,153:160) = RESHAPE ((/ & ! 3850-4025
.836e+00_EB, .584e+00_EB, .277e+00_EB, .184e+00_EB, .161e+00_EB, .145e+00_EB,     &
.840e+00_EB, .728e+00_EB, .422e+00_EB, .236e+00_EB, .197e+00_EB, .167e+00_EB,     &
.505e+00_EB, .500e+00_EB, .379e+00_EB, .276e+00_EB, .227e+00_EB, .192e+00_EB,     &
.117e+00_EB, .400e+00_EB, .423e+00_EB, .315e+00_EB, .243e+00_EB, .202e+00_EB,     &
.460e-01_EB, .300e+00_EB, .358e+00_EB, .290e+00_EB, .230e+00_EB, .202e+00_EB,     &
.183e-01_EB, .205e+00_EB, .269e+00_EB, .235e+00_EB, .195e+00_EB, .192e+00_EB,     &
.730e-02_EB, .135e+00_EB, .186e+00_EB, .179e+00_EB, .159e+00_EB, .168e+00_EB,     &
.557e-02_EB, .790e-01_EB, .113e+00_EB, .124e+00_EB, .124e+00_EB, .134e+00_EB/),  & 
(/6,8/))
sd(1:6,161:168) = RESHAPE ((/ & ! 4050-4225
.283e-02_EB, .415e-01_EB, .662e-01_EB, .886e-01_EB, .103e+00_EB, .106e+00_EB,     &
.226e-02_EB, .197e-01_EB, .367e-01_EB, .594e-01_EB, .801e-01_EB, .879e-01_EB,     &
.155e-02_EB, .860e-02_EB, .211e-01_EB, .395e-01_EB, .503e-01_EB, .610e-01_EB,     &
.103e-02_EB, .521e-02_EB, .119e-01_EB, .246e-01_EB, .354e-01_EB, .480e-01_EB,     &
.821e-03_EB, .365e-02_EB, .759e-02_EB, .166e-01_EB, .258e-01_EB, .370e-01_EB,     &
.752e-03_EB, .183e-02_EB, .445e-02_EB, .100e-01_EB, .179e-01_EB, .268e-01_EB,     &
.429e-03_EB, .141e-02_EB, .354e-02_EB, .821e-02_EB, .142e-01_EB, .212e-01_EB,     &
.327e-03_EB, .902e-03_EB, .209e-02_EB, .588e-02_EB, .112e-01_EB, .172e-01_EB/),   &
(/6,8/))
sd(1:6,169:176) = RESHAPE ((/ & ! 4250-4425
.225e-03_EB, .685e-03_EB, .189e-02_EB, .512e-02_EB, .101e-01_EB, .164e-01_EB,     &
.186e-03_EB, .551e-03_EB, .156e-02_EB, .366e-02_EB, .812e-02_EB, .136e-01_EB,     &
.173e-03_EB, .472e-03_EB, .139e-02_EB, .306e-02_EB, .661e-02_EB, .115e-01_EB,     &
.138e-03_EB, .395e-03_EB, .110e-02_EB, .272e-02_EB, .587e-02_EB, .104e-01_EB,     &
.900e-04_EB, .270e-03_EB, .968e-03_EB, .222e-02_EB, .497e-02_EB, .921e-02_EB,     &
.752e-04_EB, .233e-03_EB, .744e-03_EB, .208e-02_EB, .466e-02_EB, .876e-02_EB,     &
.618e-04_EB, .175e-03_EB, .638e-03_EB, .185e-02_EB, .465e-02_EB, .914e-02_EB,     &
.504e-04_EB, .134e-03_EB, .499e-03_EB, .174e-02_EB, .455e-02_EB, .935e-02_EB/),   &
(/6,8/))
sd(1:6,177:184) = RESHAPE ((/ & ! 4450-4625
.375e-04_EB, .123e-03_EB, .485e-03_EB, .182e-02_EB, .456e-02_EB, .971e-02_EB,     &
.305e-04_EB, .892e-04_EB, .338e-03_EB, .134e-02_EB, .460e-02_EB, .104e-01_EB,     &
.257e-04_EB, .790e-04_EB, .329e-03_EB, .154e-02_EB, .477e-02_EB, .112e-01_EB,     &
.242e-04_EB, .740e-04_EB, .308e-03_EB, .135e-02_EB, .497e-02_EB, .122e-01_EB,     &
.215e-04_EB, .653e-04_EB, .282e-03_EB, .131e-02_EB, .521e-02_EB, .133e-01_EB,     &
.218e-04_EB, .660e-04_EB, .272e-03_EB, .152e-02_EB, .573e-02_EB, .148e-01_EB,     &
.215e-04_EB, .671e-04_EB, .268e-03_EB, .134e-02_EB, .607e-02_EB, .159e-01_EB,     &
.217e-04_EB, .695e-04_EB, .285e-03_EB, .161e-02_EB, .677e-02_EB, .173e-01_EB/),  & 
(/6,8/))
sd(1:6,185:192) = RESHAPE ((/ & ! 4650-4825
.219e-04_EB, .722e-04_EB, .297e-03_EB, .169e-02_EB, .783e-02_EB, .197e-01_EB,     &
.226e-04_EB, .771e-04_EB, .341e-03_EB, .236e-02_EB, .925e-02_EB, .226e-01_EB,     &
.250e-04_EB, .815e-04_EB, .387e-03_EB, .286e-02_EB, .106e-01_EB, .250e-01_EB,    &
.280e-04_EB, .845e-04_EB, .420e-03_EB, .357e-02_EB, .124e-01_EB, .276e-01_EB,     &
.351e-04_EB, .192e-03_EB, .470e-03_EB, .467e-02_EB, .166e-01_EB, .313e-01_EB,     &
.435e-04_EB, .200e-03_EB, .105e-02_EB, .566e-02_EB, .185e-01_EB, .341e-01_EB,     &
.522e-04_EB, .233e-03_EB, .129e-02_EB, .736e-02_EB, .229e-01_EB, .378e-01_EB,     &
.673e-04_EB, .306e-03_EB, .183e-02_EB, .982e-02_EB, .258e-01_EB, .404e-01_EB/),   &
(/6,8/))
sd(1:6,193:200) = RESHAPE ((/ & ! 4850-5025
.886e-04_EB, .399e-03_EB, .246e-02_EB, .128e-01_EB, .302e-01_EB, .430e-01_EB,     &
.113e-03_EB, .618e-03_EB, .346e-02_EB, .161e-01_EB, .358e-01_EB, .459e-01_EB,     &
.174e-03_EB, .825e-03_EB, .441e-02_EB, .200e-01_EB, .417e-01_EB, .493e-01_EB,     &
.265e-03_EB, .163e-02_EB, .777e-02_EB, .245e-01_EB, .450e-01_EB, .507e-01_EB,     &
.355e-03_EB, .200e-02_EB, .978e-02_EB, .317e-01_EB, .492e-01_EB, .527e-01_EB,     &
.538e-03_EB, .271e-02_EB, .167e-01_EB, .401e-01_EB, .503e-01_EB, .523e-01_EB,     &
.651e-03_EB, .301e-02_EB, .264e-01_EB, .467e-01_EB, .520e-01_EB, .526e-01_EB,     &
.987e-03_EB, .530e-02_EB, .321e-01_EB, .499e-01_EB, .523e-01_EB, .510e-01_EB/),   &
(/6,8/))
sd(1:6,201:208) = RESHAPE ((/ & ! 5050-5225
.135e-02_EB, .860e-02_EB, .389e-01_EB, .528e-01_EB, .513e-01_EB, .492e-01_EB,    &
.226e-02_EB, .130e-01_EB, .472e-01_EB, .559e-01_EB, .500e-01_EB, .469e-01_EB,     &
.431e-02_EB, .198e-01_EB, .526e-01_EB, .557e-01_EB, .480e-01_EB, .452e-01_EB,     &
.628e-02_EB, .282e-01_EB, .488e-01_EB, .495e-01_EB, .451e-01_EB, .430e-01_EB,     &
.900e-02_EB, .390e-01_EB, .471e-01_EB, .449e-01_EB, .430e-01_EB, .423e-01_EB,     &
.180e-01_EB, .462e-01_EB, .412e-01_EB, .391e-01_EB, .403e-01_EB, .415e-01_EB,     &
.348e-01_EB, .710e-01_EB, .402e-01_EB, .360e-01_EB, .384e-01_EB, .414e-01_EB,     &
.718e-01_EB, .590e-01_EB, .399e-01_EB, .360e-01_EB, .376e-01_EB, .420e-01_EB/),   &
(/6,8/))
sd(1:6,209:216) = RESHAPE ((/ & ! 5250-5425
.111e+00_EB, .368e-01_EB, .340e-01_EB, .369e-01_EB, .409e-01_EB, .454e-01_EB,     &
.329e-01_EB, .285e-01_EB, .365e-01_EB, .423e-01_EB, .461e-01_EB, .482e-01_EB,     &
.281e-01_EB, .270e-01_EB, .432e-01_EB, .505e-01_EB, .529e-01_EB, .511e-01_EB,     &
.121e+00_EB, .422e-01_EB, .589e-01_EB, .598e-01_EB, .572e-01_EB, .544e-01_EB,     &
.139e+00_EB, .105e+00_EB, .844e-01_EB, .687e-01_EB, .593e-01_EB, .560e-01_EB,     &
.774e-01_EB, .710e-01_EB, .683e-01_EB, .618e-01_EB, .556e-01_EB, .534e-01_EB,     &
.858e-01_EB, .483e-01_EB, .579e-01_EB, .547e-01_EB, .503e-01_EB, .495e-01_EB,     &
.985e-01_EB, .575e-01_EB, .589e-01_EB, .510e-01_EB, .451e-01_EB, .449e-01_EB/),   &
(/6,8/))
sd(1:6,217:224) = RESHAPE ((/ & ! 5450-5625
.996e-01_EB, .682e-01_EB, .539e-01_EB, .489e-01_EB, .454e-01_EB, .446e-01_EB,     &
.680e-01_EB, .680e-01_EB, .548e-01_EB, .495e-01_EB, .460e-01_EB, .458e-01_EB,     &
.325e-01_EB, .520e-01_EB, .515e-01_EB, .483e-01_EB, .449e-01_EB, .454e-01_EB,     &
.150e-01_EB, .350e-01_EB, .451e-01_EB, .464e-01_EB, .452e-01_EB, .449e-01_EB,     &
.620e-02_EB, .238e-01_EB, .369e-01_EB, .408e-01_EB, .414e-01_EB, .417e-01_EB,     &
.270e-02_EB, .158e-01_EB, .282e-01_EB, .339e-01_EB, .366e-01_EB, .384e-01_EB,     &
.113e-02_EB, .101e-01_EB, .203e-01_EB, .263e-01_EB, .303e-01_EB, .333e-01_EB,     &
.829e-03_EB, .590e-02_EB, .148e-01_EB, .206e-01_EB, .247e-01_EB, .295e-01_EB/),  & 
(/6,8/))
sd(1:6,225:232) = RESHAPE ((/ & ! 5650-5825
.365e-03_EB, .310e-02_EB, .969e-02_EB, .154e-01_EB, .203e-01_EB, .258e-01_EB,     &
.240e-03_EB, .130e-02_EB, .589e-02_EB, .112e-01_EB, .164e-01_EB, .222e-01_EB,     &
.158e-03_EB, .400e-03_EB, .417e-02_EB, .850e-02_EB, .134e-01_EB, .190e-01_EB,     &
.103e-03_EB, .262e-03_EB, .208e-02_EB, .594e-02_EB, .109e-01_EB, .162e-01_EB,     &
.741e-04_EB, .181e-03_EB, .142e-02_EB, .455e-02_EB, .907e-02_EB, .141e-01_EB,     &
.625e-04_EB, .135e-03_EB, .816e-03_EB, .316e-02_EB, .698e-02_EB, .121e-01_EB,     &
.499e-04_EB, .111e-03_EB, .624e-03_EB, .230e-02_EB, .551e-02_EB, .102e-01_EB,     &
.325e-04_EB, .677e-04_EB, .425e-03_EB, .124e-02_EB, .385e-02_EB, .818e-02_EB/),   &
(/6,8/))
sd(1:6,233:240) = RESHAPE ((/ & ! 5850-6025
.231e-04_EB, .563e-04_EB, .278e-03_EB, .986e-03_EB, .290e-02_EB, .672e-02_EB,    &
.165e-04_EB, .481e-04_EB, .247e-03_EB, .944e-03_EB, .253e-02_EB, .612e-02_EB,           &
.126e-04_EB, .432e-04_EB, .241e-03_EB, .886e-03_EB, .220e-02_EB, .582e-02_EB,          &
.118e-04_EB, .420e-04_EB, .235e-03_EB, .847e-03_EB, .209e-02_EB, .571e-02_EB,          &
.110e-04_EB, .408e-04_EB, .226e-03_EB, .812e-03_EB, .221e-02_EB, .604e-02_EB,         &
.101e-04_EB, .400e-04_EB, .213e-03_EB, .805e-03_EB, .239e-02_EB, .641e-02_EB,        &
.983e-05_EB, .395e-04_EB, .186e-03_EB, .801e-03_EB, .247e-02_EB, .691e-02_EB,       &
.979e-05_EB, .401e-04_EB, .193e-03_EB, .805e-03_EB, .260e-02_EB, .732e-02_EB/),    &
(/6,8/))
sd(1:6,241:248) = RESHAPE ((/ & ! 6050-6225
.976e-05_EB, .410e-04_EB, .201e-03_EB, .814e-03_EB, .285e-02_EB, .776e-02_EB,     &
.988e-05_EB, .420e-04_EB, .210e-03_EB, .832e-03_EB, .317e-02_EB, .842e-02_EB,    &
.991e-05_EB, .425e-04_EB, .219e-03_EB, .877e-03_EB, .340e-02_EB, .888e-02_EB,        &
.102e-04_EB, .435e-04_EB, .231e-03_EB, .937e-03_EB, .361e-02_EB, .929e-02_EB,       &
.110e-04_EB, .486e-04_EB, .244e-03_EB, .971e-03_EB, .402e-02_EB, .994e-02_EB,      &
.127e-04_EB, .579e-04_EB, .257e-03_EB, .111e-02_EB, .437e-02_EB, .104e-01_EB,     &
.131e-04_EB, .612e-04_EB, .277e-03_EB, .113e-02_EB, .465e-02_EB, .110e-01_EB,        &
.150e-04_EB, .783e-04_EB, .353e-03_EB, .116e-02_EB, .510e-02_EB, .116e-01_EB/),     &
(/6,8/))
sd(1:6,249:256) = RESHAPE ((/ & ! 6250-6425
.178e-04_EB, .922e-04_EB, .394e-03_EB, .157e-02_EB, .555e-02_EB, .123e-01_EB,      &
.203e-04_EB, .115e-03_EB, .481e-03_EB, .188e-02_EB, .601e-02_EB, .131e-01_EB,     &
.230e-04_EB, .145e-03_EB, .617e-03_EB, .183e-02_EB, .644e-02_EB, .139e-01_EB,      & 
.280e-04_EB, .187e-03_EB, .723e-03_EB, .202e-02_EB, .686e-02_EB, .146e-01_EB,     & 
.305e-04_EB, .209e-03_EB, .811e-03_EB, .243e-02_EB, .779e-02_EB, .157e-01_EB,     &
.455e-04_EB, .244e-03_EB, .935e-03_EB, .243e-02_EB, .844e-02_EB, .166e-01_EB,    &
.661e-04_EB, .320e-03_EB, .989e-03_EB, .288e-02_EB, .902e-02_EB, .173e-01_EB,       & 
.723e-04_EB, .397e-03_EB, .122e-02_EB, .359e-02_EB, .100e-01_EB, .184e-01_EB/),    &
(/6,8/))
sd(1:6,257:264) = RESHAPE ((/ & ! 6450-6625
.847e-04_EB, .481e-03_EB, .143e-02_EB, .429e-02_EB, .108e-01_EB, .192e-01_EB,     &
.103e-03_EB, .591e-03_EB, .174e-02_EB, .488e-02_EB, .116e-01_EB, .200e-01_EB,           &
.131e-03_EB, .703e-03_EB, .247e-02_EB, .549e-02_EB, .124e-01_EB, .205e-01_EB,          &
.165e-03_EB, .872e-03_EB, .265e-02_EB, .641e-02_EB, .131e-01_EB, .211e-01_EB,         &
.205e-03_EB, .110e-02_EB, .298e-02_EB, .749e-02_EB, .140e-01_EB, .218e-01_EB,         &
.253e-03_EB, .130e-02_EB, .346e-02_EB, .811e-02_EB, .150e-01_EB, .230e-01_EB,         &
.338e-03_EB, .150e-02_EB, .445e-02_EB, .890e-02_EB, .159e-01_EB, .237e-01_EB,        &
.437e-03_EB, .170e-02_EB, .491e-02_EB, .107e-01_EB, .170e-01_EB, .245e-01_EB/),     &
(/6,8/))
sd(1:6,265:272) = RESHAPE ((/ & ! 6650-6825
.581e-03_EB, .190e-02_EB, .537e-02_EB, .116e-01_EB, .179e-01_EB, .254e-01_EB,         &
.685e-03_EB, .220e-02_EB, .578e-02_EB, .128e-01_EB, .189e-01_EB, .263e-01_EB,        &
.900e-03_EB, .250e-02_EB, .649e-02_EB, .134e-01_EB, .195e-01_EB, .275e-01_EB,       &
.121e-02_EB, .280e-02_EB, .722e-02_EB, .142e-01_EB, .202e-01_EB, .281e-01_EB,      &
.152e-02_EB, .330e-02_EB, .813e-02_EB, .161e-01_EB, .212e-01_EB, .288e-01_EB,     &
.185e-02_EB, .370e-02_EB, .907e-02_EB, .168e-01_EB, .222e-01_EB, .292e-01_EB,    &
.220e-02_EB, .430e-02_EB, .929e-02_EB, .183e-01_EB, .233e-01_EB, .294e-01_EB,       &
.255e-02_EB, .500e-02_EB, .114e-01_EB, .195e-01_EB, .245e-01_EB, .289e-01_EB/),       &
(/6,8/))
sd(1:6,273:280) = RESHAPE ((/ & ! 6850-7025
.290e-02_EB, .580e-02_EB, .167e-01_EB, .215e-01_EB, .260e-01_EB, .291e-01_EB,        &
.320e-02_EB, .670e-02_EB, .208e-01_EB, .237e-01_EB, .274e-01_EB, .293e-01_EB,        &
.360e-02_EB, .880e-02_EB, .220e-01_EB, .253e-01_EB, .282e-01_EB, .300e-01_EB,       &
.400e-02_EB, .920e-02_EB, .238e-01_EB, .273e-01_EB, .290e-01_EB, .304e-01_EB,        &
.460e-02_EB, .108e-01_EB, .272e-01_EB, .279e-01_EB, .298e-01_EB, .310e-01_EB,       &
.530e-02_EB, .128e-01_EB, .304e-01_EB, .292e-01_EB, .297e-01_EB, .312e-01_EB,        &
.620e-02_EB, .152e-01_EB, .344e-01_EB, .303e-01_EB, .293e-01_EB, .310e-01_EB,       &
.760e-02_EB, .182e-01_EB, .341e-01_EB, .297e-01_EB, .290e-01_EB, .300e-01_EB/),    &
(/6,8/))
sd(1:6,281:288) = RESHAPE ((/ & ! 7050-7225
.980e-02_EB, .222e-01_EB, .398e-01_EB, .318e-01_EB, .291e-01_EB, .294e-01_EB,       & 
.132e-01_EB, .271e-01_EB, .402e-01_EB, .294e-01_EB, .274e-01_EB, .282e-01_EB,      & 
.190e-01_EB, .335e-01_EB, .421e-01_EB, .286e-01_EB, .262e-01_EB, .269e-01_EB,     & 
.240e-01_EB, .432e-01_EB, .431e-01_EB, .276e-01_EB, .245e-01_EB, .257e-01_EB,    & 
.288e-01_EB, .570e-01_EB, .458e-01_EB, .270e-01_EB, .228e-01_EB, .243e-01_EB,   &
.323e-01_EB, .740e-01_EB, .449e-01_EB, .261e-01_EB, .214e-01_EB, .221e-01_EB,        &
.570e-01_EB, .890e-01_EB, .435e-01_EB, .225e-01_EB, .199e-01_EB, .196e-01_EB,       &
.216e-01_EB, .680e-01_EB, .378e-01_EB, .239e-01_EB, .195e-01_EB, .192e-01_EB/),    &
(/6,8/))
sd(1:6,289:296) = RESHAPE ((/ & ! 7250-7425
.126e-01_EB, .475e-01_EB, .364e-01_EB, .238e-01_EB, .197e-01_EB, .192e-01_EB,      & 
.117e-01_EB, .369e-01_EB, .385e-01_EB, .249e-01_EB, .212e-01_EB, .204e-01_EB,     & 
.140e-01_EB, .370e-01_EB, .419e-01_EB, .272e-01_EB, .228e-01_EB, .213e-01_EB,    & 
.425e-01_EB, .418e-01_EB, .440e-01_EB, .280e-01_EB, .248e-01_EB, .229e-01_EB,   &
.640e-01_EB, .460e-01_EB, .427e-01_EB, .290e-01_EB, .263e-01_EB, .238e-01_EB,   &
.385e-01_EB, .385e-01_EB, .374e-01_EB, .259e-01_EB, .235e-01_EB, .224e-01_EB,        &
.182e-01_EB, .179e-01_EB, .282e-01_EB, .231e-01_EB, .211e-01_EB, .214e-01_EB,       &
.170e-01_EB, .810e-02_EB, .191e-01_EB, .175e-01_EB, .181e-01_EB, .194e-01_EB/),     & 
(/6,8/))
sd(1:6,297:304) = RESHAPE ((/ & ! 7450-7625
.161e-01_EB, .370e-02_EB, .105e-01_EB, .127e-01_EB, .152e-01_EB, .171e-01_EB,      & 
.145e-01_EB, .170e-02_EB, .554e-02_EB, .855e-02_EB, .113e-01_EB, .131e-01_EB,     & 
.175e-02_EB, .140e-02_EB, .385e-02_EB, .595e-02_EB, .803e-02_EB, .945e-02_EB,    &
.772e-03_EB, .751e-03_EB, .384e-02_EB, .575e-02_EB, .537e-02_EB, .594e-02_EB,       &
.491e-03_EB, .600e-03_EB, .301e-02_EB, .453e-02_EB, .380e-02_EB, .434e-02_EB,      &
.275e-03_EB, .410e-03_EB, .193e-02_EB, .366e-02_EB, .319e-02_EB, .332e-02_EB,     &
.185e-01_EB, .280e-03_EB, .131e-02_EB, .232e-02_EB, .247e-02_EB, .256e-02_EB,      &
.101e-03_EB, .160e-03_EB, .915e-03_EB, .150e-02_EB, .186e-02_EB, .197e-02_EB/),   &
(/6,8/))
sd(1:6,305:312) = RESHAPE ((/ & ! 7650-7825
.691e-04_EB, .110e-03_EB, .565e-03_EB, .114e-02_EB, .205e-02_EB, .192e-02_EB,        &  
.476e-04_EB, .750e-04_EB, .114e-02_EB, .124e-02_EB, .175e-02_EB, .187e-02_EB,         &
.305e-04_EB, .590e-04_EB, .529e-03_EB, .114e-02_EB, .160e-02_EB, .185e-02_EB,       &
.240e-04_EB, .480e-04_EB, .293e-03_EB, .842e-03_EB, .141e-02_EB, .184e-02_EB,       &
.170e-04_EB, .360e-04_EB, .122e-03_EB, .435e-03_EB, .124e-02_EB, .182e-02_EB,      &
.120e-04_EB, .240e-04_EB, .121e-03_EB, .435e-03_EB, .118e-02_EB, .187e-02_EB,     &
.810e-05_EB, .170e-04_EB, .103e-03_EB, .439e-03_EB, .126e-02_EB, .192e-02_EB,    &
.550e-05_EB, .120e-04_EB, .866e-04_EB, .367e-03_EB, .119e-02_EB, .193e-02_EB/),    &
(/6,8/))
sd(1:6,313:320) = RESHAPE ((/&  ! 7850-8025
.390e-05_EB, .900e-05_EB, .716e-04_EB, .351e-03_EB, .116e-02_EB, .194e-02_EB,        & 
.295e-05_EB, .830e-05_EB, .373e-04_EB, .254e-03_EB, .114e-02_EB, .196e-02_EB,       & 
.230e-05_EB, .800e-05_EB, .465e-04_EB, .298e-03_EB, .117e-02_EB, .201e-02_EB,      & 
.225e-05_EB, .820e-05_EB, .367e-04_EB, .252e-03_EB, .116e-02_EB, .205e-02_EB,     & 
.220e-05_EB, .840e-05_EB, .371e-04_EB, .268e-03_EB, .127e-02_EB, .211e-02_EB,    & 
.223e-05_EB, .920e-05_EB, .396e-04_EB, .273e-03_EB, .128e-02_EB, .216e-02_EB,   &
.235e-05_EB, .103e-04_EB, .415e-04_EB, .263e-03_EB, .121e-02_EB, .221e-02_EB,      &
.280e-05_EB, .125e-04_EB, .633e-04_EB, .363e-03_EB, .136e-02_EB, .231e-02_EB/),   &
(/6,8/))
sd(1:6,321:328) = RESHAPE ((/ & ! 8050-8225
.310e-05_EB, .150e-04_EB, .979e-04_EB, .492e-03_EB, .150e-02_EB, .241e-02_EB,           &
.370e-05_EB, .180e-04_EB, .120e-03_EB, .580e-03_EB, .167e-02_EB, .251e-02_EB,           &
.420e-05_EB, .200e-04_EB, .987e-04_EB, .509e-03_EB, .171e-02_EB, .257e-02_EB,           &
.510e-05_EB, .240e-04_EB, .134e-03_EB, .547e-03_EB, .173e-02_EB, .267e-02_EB,           &
.600e-05_EB, .270e-04_EB, .121e-03_EB, .534e-03_EB, .172e-02_EB, .274e-02_EB,           &
.720e-05_EB, .300e-04_EB, .204e-03_EB, .684e-03_EB, .184e-02_EB, .285e-02_EB,           &
.820e-05_EB, .330e-04_EB, .276e-03_EB, .819e-03_EB, .199e-02_EB, .297e-02_EB,           &
.100e-04_EB, .380e-04_EB, .317e-03_EB, .859e-03_EB, .214e-02_EB, .308e-02_EB/),         &
(/6,8/))
sd(1:6,329:336) = RESHAPE ((/ & ! 8250-8425
.125e-04_EB, .420e-04_EB, .240e-03_EB, .818e-03_EB, .220e-02_EB, .317e-02_EB,           &
.145e-04_EB, .500e-04_EB, .452e-03_EB, .109e-02_EB, .238e-02_EB, .293e-02_EB,           &
.175e-04_EB, .560e-04_EB, .301e-03_EB, .941e-03_EB, .243e-02_EB, .342e-02_EB,           &
.198e-04_EB, .630e-04_EB, .280e-03_EB, .107e-02_EB, .260e-02_EB, .353e-02_EB,           &
.230e-04_EB, .710e-04_EB, .276e-03_EB, .109e-02_EB, .272e-02_EB, .365e-02_EB,           &
.280e-04_EB, .830e-04_EB, .369e-03_EB, .127e-02_EB, .295e-02_EB, .377e-02_EB,           &
.330e-04_EB, .890e-04_EB, .430e-03_EB, .139e-02_EB, .306e-02_EB, .385e-02_EB,           &
.360e-04_EB, .950e-04_EB, .371e-03_EB, .135e-02_EB, .306e-02_EB, .384e-02_EB/),         &
(/6,8/))
sd(1:6,337:344) = RESHAPE ((/ & ! 8450-8625
.390e-04_EB, .980e-04_EB, .434e-03_EB, .147e-02_EB, .316e-02_EB, .385e-02_EB,           &
.400e-04_EB, .990e-04_EB, .397e-03_EB, .143e-02_EB, .318e-02_EB, .384e-02_EB,           &
.400e-04_EB, .980e-04_EB, .364e-03_EB, .141e-02_EB, .317e-02_EB, .381e-02_EB,           &
.390e-04_EB, .940e-04_EB, .390e-03_EB, .142e-02_EB, .314e-02_EB, .376e-02_EB,          &
.380e-04_EB, .900e-04_EB, .380e-03_EB, .145e-02_EB, .318e-02_EB, .375e-02_EB,         &
.380e-04_EB, .900e-04_EB, .380e-03_EB, .145e-02_EB, .318e-02_EB, .375e-02_EB,        &
.330e-04_EB, .750e-04_EB, .358e-03_EB, .138e-02_EB, .310e-02_EB, .372e-02_EB,       &
.270e-04_EB, .580e-04_EB, .382e-03_EB, .143e-02_EB, .315e-02_EB, .369e-02_EB/),    &
(/6,8/))
sd(1:6,345:352) = RESHAPE ((/ & ! 8650-8825
.240e-04_EB, .500e-04_EB, .343e-03_EB, .136e-02_EB, .306e-02_EB, .363e-02_EB,          &
.200e-04_EB, .450e-04_EB, .309e-03_EB, .134e-02_EB, .306e-02_EB, .359e-02_EB,         & 
.180e-04_EB, .400e-04_EB, .281e-03_EB, .127e-02_EB, .294e-02_EB, .341e-02_EB,        &
.170e-04_EB, .360e-04_EB, .276e-03_EB, .124e-02_EB, .290e-02_EB, .336e-02_EB,       & 
.160e-04_EB, .310e-04_EB, .272e-03_EB, .122e-02_EB, .283e-02_EB, .323e-02_EB,      &  
.140e-04_EB, .280e-04_EB, .241e-03_EB, .117e-02_EB, .273e-02_EB, .309e-02_EB,       &
.120e-04_EB, .250e-04_EB, .237e-03_EB, .115e-02_EB, .269e-02_EB, .297e-02_EB,      &
.100e-04_EB, .220e-04_EB, .218e-03_EB, .111e-02_EB, .259e-02_EB, .284e-02_EB/),   &
(/6,8/))
sd(1:6,353:360) = RESHAPE ((/ & ! 8850-9025
.920e-05_EB, .198e-04_EB, .206e-03_EB, .105e-02_EB, .246e-02_EB, .269e-02_EB,       &
.810e-05_EB, .170e-04_EB, .205e-03_EB, .100e-02_EB, .235e-02_EB, .257e-02_EB,       & 
.720e-05_EB, .160e-04_EB, .177e-03_EB, .921e-03_EB, .220e-02_EB, .245e-02_EB,      &
.650e-05_EB, .150e-04_EB, .172e-03_EB, .834e-03_EB, .205e-02_EB, .232e-02_EB,     &
.590e-05_EB, .130e-04_EB, .147e-03_EB, .735e-03_EB, .194e-02_EB, .218e-02_EB,    &
.510e-05_EB, .110e-04_EB, .120e-03_EB, .629e-03_EB, .177e-02_EB, .203e-02_EB,     &
.460e-05_EB, .950e-05_EB, .960e-04_EB, .513e-03_EB, .154e-02_EB, .180e-02_EB,      &
.420e-05_EB, .800e-05_EB, .578e-04_EB, .314e-03_EB, .123e-02_EB, .154e-02_EB/),     &
(/6,8/))
sd(1:6,361:368) = RESHAPE ((/ & ! 9050-9225
.380e-05_EB, .720e-05_EB, .529e-04_EB, .292e-03_EB, .114e-02_EB, .137e-02_EB,      &
.330e-05_EB, .660e-05_EB, .485e-04_EB, .269e-03_EB, .102e-02_EB, .122e-02_EB,      &
.290e-05_EB, .580e-05_EB, .430e-04_EB, .239e-03_EB, .896e-03_EB, .107e-02_EB,      &
.270e-05_EB, .520e-05_EB, .259e-04_EB, .193e-03_EB, .748e-03_EB, .944e-03_EB,     &
.240e-05_EB, .450e-05_EB, .316e-04_EB, .207e-03_EB, .671e-02_EB, .848e-03_EB,      &
.220e-05_EB, .400e-05_EB, .444e-05_EB, .602e-04_EB, .516e-03_EB, .750e-03_EB,     &
.190e-05_EB, .360e-05_EB, .324e-05_EB, .460e-04_EB, .439e-03_EB, .688e-03_EB,     & 
.170e-05_EB, .320e-05_EB, .180e-05_EB, .321e-04_EB, .384e-03_EB, .653e-03_EB/),  & 
(/6,8/))
sd(1:6,369:376) = RESHAPE ((/ & ! 9250-9300
.140e-05_EB, .280e-05_EB, .171e-05_EB, .344e-04_EB, .340e-03_EB, .616e-03_EB,   & 
.130e-05_EB, .250e-05_EB, .299e-05_EB, .600e-04_EB, .343e-03_EB, .619e-03_EB,  &  
.120e-05_EB, .220e-05_EB, .299e-05_EB, .600e-04_EB, .343e-03_EB, .619e-03_EB, &  
1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB,&
1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB/),&
(/6,8/))

! intialize gamma array

gamma = RESHAPE((/ &
!     line broadening PARAMETERs,gamma(i,j),
!     j=co2,h2o,ch4,co,o2,n2,self resonant.
!   i=     co2  h2o  ch4   co
!    j
.09_EB , .12_EB, .0_EB , .07_EB,&
.07_EB , .09_EB, .0_EB , .06_EB,&
!.0_EB  , .0_EB , .16_EB, .0_EB ,&
.0_EB  , .0_EB , .12_EB, .0_EB ,&
.06_EB , .10_EB, .0_EB , .06_EB,&
.055_EB, .04_EB, .0_EB , .05_EB,&
.07_EB , .09_EB, .0_EB , .06_EB, &
.01_EB , .44_EB, .0_EB , .0_EB /),(/4,7/))

! initialize sd15 array

!  the following are data for the 15.0 micron band of co2
!  temp, k=300     600      1200      1500  1800    2400    

sd15(1:6,1:8) = RESHAPE ((/  &! 500-535
.000e+00_EB, .000e+00_EB, .000e+00_EB, .105e-01_EB, .300e-01_EB, .880e-01_EB,          &
.000e+00_EB, .000e+00_EB, .000e+00_EB, .180e-01_EB, .490e-01_EB, .880e-01_EB,         & 
.000e+00_EB, .000e+00_EB, .000e+00_EB, .300e-01_EB, .540e-01_EB, .740e-01_EB,        &  
.000e+00_EB, .000e+00_EB, .000e+00_EB, .300e-01_EB, .560e-01_EB, .890e-01_EB,       &   
.000e+00_EB, .000e+00_EB, .000e+00_EB, .330e-01_EB, .690e-01_EB, .990e-01_EB,      &    
.000e+00_EB, .000e+00_EB, .880e-02_EB, .380e-01_EB, .720e-01_EB, .970e-01_EB,     &     
.000e+00_EB, .000e+00_EB, .110e-01_EB, .530e-01_EB, .950e-01_EB, .124e+00_EB,    &      
.000e+00_EB, .000e+00_EB, .285e-01_EB, .630e-01_EB, .990e-01_EB, .140e+00_EB/), &       
(/6,8/))
sd15(1:6,9:16) = RESHAPE ((/  &! 540-575
.000e+00_EB, .000e+00_EB, .330e-01_EB, .680e-01_EB, .103e+00_EB, .134e+00_EB,          &
.000e+00_EB, .000e+00_EB, .450e-01_EB, .920e-01_EB, .138e+00_EB, .176e+00_EB,          &
.000e+00_EB, .000e+00_EB, .490e-01_EB, .970e-01_EB, .148e+00_EB, .191e+00_EB,         & 
.000e+00_EB, .000e+00_EB, .490e-01_EB, .120e-01_EB, .188e+00_EB, .247e+00_EB,        &  
.000e+00_EB, .000e+00_EB, .480e-01_EB, .126e+00_EB, .201e+00_EB, .241e+00_EB,       &   
.000e+00_EB, .000e+00_EB, .820e-01_EB, .198e+00_EB, .270e+00_EB, .265e+00_EB,      &    
.000e+00_EB, .750e-02_EB, .690e-01_EB, .140e+00_EB, .225e+00_EB, .340e+00_EB,     &     
.000e+00_EB, .205e-01_EB, .820e-01_EB, .145e+00_EB, .236e+00_EB, .530e+00_EB/),  &      
(/6,8/))
sd15(1:6,17:24) = RESHAPE ((/ & ! 580-615
.000e+00_EB, .355e-01_EB, .117e+00_EB, .193e+00_EB, .295e+00_EB, .550e+00_EB,          &
.157e-01_EB, .520e-01_EB, .170e+00_EB, .235e+00_EB, .305e+00_EB, .410e+00_EB,          &
.150e-01_EB, .880e-01_EB, .270e+00_EB, .330e+00_EB, .440e+00_EB, .520e+00_EB,          &
.510e-01_EB, .130e+00_EB, .400e+00_EB, .530e+00_EB, .560e+00_EB, .540e+00_EB,         & 
.120e+00_EB, .165e+00_EB, .275e+00_EB, .320e+00_EB, .420e+00_EB, .560e+00_EB,        &  
.880e-01_EB, .190e+00_EB, .430e+00_EB, .540e+00_EB, .620e+00_EB, .680e+00_EB,       &   
.110e+00_EB, .350e+00_EB, .710e+00_EB, .760e+00_EB, .760e+00_EB, .690e+00_EB,      &    
.180e+00_EB, .470e+00_EB, .920e+00_EB, .970e+00_EB, .910e+00_EB, .670e+00_EB/),   &     
(/6,8/))
sd15(1:6,25:32) = RESHAPE ((/ & ! 620-655
.970e-01_EB, .265e+00_EB, .610e+00_EB, .720e+00_EB, .780e+00_EB, .730e+00_EB,          &
.175e+00_EB, .380e+00_EB, .720e+00_EB, .790e+00_EB, .830e+00_EB, .840e+00_EB,          &
.370e+00_EB, .640e+00_EB, .920e+00_EB, .960e+00_EB, .980e+00_EB, .940e+00_EB,          &
.590e+00_EB, .840e+00_EB, .107e+01_EB, .110e+01_EB, .111e+01_EB, .106e+01_EB,         & 
.940e+00_EB, .103e+01_EB, .115e+01_EB, .115e+01_EB, .115e+01_EB, .118e+01_EB,        &  
.196e+01_EB, .177e+01_EB, .146e+01_EB, .136e+01_EB, .132e+01_EB, .139e+01_EB,       &   
.345e+01_EB, .282e+01_EB, .198e+01_EB, .172e+01_EB, .156e+01_EB, .148e+01_EB,      &    
.282e+01_EB, .248e+01_EB, .200e+01_EB, .190e+01_EB, .186e+01_EB, .205e+01_EB/),   &     
(/6,8/))
sd15(1:6,33:40) = RESHAPE ((/  &! 660-695
.254e+01_EB, .234e+01_EB, .184e+01_EB, .176e+01_EB, .174e+01_EB, .203e+01_EB,          &
.142e+02_EB, .860e+01_EB, .370e+01_EB, .260e+01_EB, .196e+01_EB, .142e+01_EB,          &
.450e+01_EB, .570e+01_EB, .580e+01_EB, .520e+01_EB, .350e+01_EB, .420e+01_EB,         & 
.360e+01_EB, .310e+01_EB, .330e+01_EB, .290e+01_EB, .205e+01_EB, .200e+01_EB,        &  
.310e+01_EB, .260e+01_EB, .200e+01_EB, .196e+01_EB, .180e+01_EB, .210e+01_EB,       &   
.240e+01_EB, .250e+01_EB, .230e+01_EB, .220e+01_EB, .170e+01_EB, .194e+01_EB,      &    
.182e+01_EB, .200e+01_EB, .218e+01_EB, .205e+01_EB, .184e+01_EB, .130e+01_EB,     &     
.104e+01_EB, .135e+01_EB, .172e+01_EB, .172e+01_EB, .165e+01_EB, .130e+01_EB/),  &      
(/6,8/))
sd15(1:6,41:48) = RESHAPE ((/  &! 700-735
.550e+00_EB, .120e+01_EB, .143e+01_EB, .147e+01_EB, .148e+01_EB, .125e+01_EB,          &
.136e+01_EB, .128e+01_EB, .128e+01_EB, .135e+01_EB, .138e+01_EB, .134e+01_EB,          &
.210e+00_EB, .780e+00_EB, .127e+01_EB, .133e+01_EB, .137e+01_EB, .132e+01_EB,          &
.190e+00_EB, .780e+00_EB, .140e+01_EB, .146e+01_EB, .147e+01_EB, .142e+01_EB,         & 
.900e+00_EB, .106e+01_EB, .140e+01_EB, .150e+01_EB, .155e+01_EB, .134e+01_EB,        &  
.720e-01_EB, .300e+00_EB, .800e+00_EB, .100e+01_EB, .115e+01_EB, .126e+01_EB,       &   
.640e-01_EB, .210e+00_EB, .560e+00_EB, .720e+00_EB, .860e+00_EB, .102e+01_EB,      &    
.680e-01_EB, .210e+00_EB, .530e+00_EB, .670e+00_EB, .790e+00_EB, .101e+01_EB/),   &     
(/6,8/))
sd15(1:6,49:56) = RESHAPE ((/ & ! 740-775
.690e-01_EB, .210e+00_EB, .540e+00_EB, .690e+00_EB, .820e+00_EB, .910e+00_EB,          &
.330e-01_EB, .140e+00_EB, .390e+00_EB, .530e+00_EB, .690e+00_EB, .770e+00_EB,          &
.230e-01_EB, .780e-01_EB, .270e+00_EB, .410e+00_EB, .560e+00_EB, .890e+00_EB,          &
.300e-01_EB, .860e-01_EB, .280e+00_EB, .400e+00_EB, .520e+00_EB, .710e+00_EB,         & 
.175e-01_EB, .620e-01_EB, .225e+00_EB, .335e+00_EB, .450e+00_EB, .660e+00_EB,       &   
.105e-01_EB, .450e-01_EB, .180e+00_EB, .280e+00_EB, .380e+00_EB, .600e+00_EB,      &    
.450e-02_EB, .300e-01_EB, .148e+00_EB, .240e+00_EB, .345e+00_EB, .570e+00_EB,     &     
.000e+00_EB, .140e-01_EB, .124e+00_EB, .205e+00_EB, .285e+00_EB, .430e+00_EB/),  &      
(/6,8/))
sd15(1:6,57:64) = RESHAPE ((/ & ! 780-815
.000e+00_EB, .115e-01_EB, .110e+00_EB, .185e+00_EB, .260e+00_EB, .375e+00_EB,          &
.000e+00_EB, .135e-01_EB, .840e-01_EB, .140e+00_EB, .205e+00_EB, .335e+00_EB,         & 
.000e+00_EB, .430e-02_EB, .650e-01_EB, .120e+00_EB, .185e+00_EB, .325e+00_EB,        &  
.000e+00_EB, .000e+00_EB, .540e-01_EB, .115e+00_EB, .180e+00_EB, .315e+00_EB,       &   
.000e+00_EB, .000e+00_EB, .440e-01_EB, .950e-01_EB, .150e+00_EB, .270e+00_EB,      &    
.000e+00_EB, .000e+00_EB, .360e-01_EB, .790e-01_EB, .125e+00_EB, .205e+00_EB,     &     
.000e+00_EB, .000e+00_EB, .250e-01_EB, .650e-01_EB, .110e+00_EB, .178e+00_EB,    &      
.000e+00_EB, .000e+00_EB, .180e-01_EB, .620e-01_EB, .103e+00_EB, .153e+00_EB/), &       
(/6,8/))
sd15(1:6,65:72) = RESHAPE ((/ & ! 820-855
.000e+00_EB, .000e+00_EB, .320e-01_EB, .580e-01_EB, .860e-01_EB, .147e+00_EB,          &
.000e+00_EB, .000e+00_EB, .800e-02_EB, .510e-01_EB, .870e-01_EB, .134e+00_EB,          &
.000e+00_EB, .000e+00_EB, .600e-02_EB, .480e-01_EB, .830e-01_EB, .133e+00_EB,          &
.000e+00_EB, .000e+00_EB, .000e+00_EB, .430e-01_EB, .780e-01_EB, .118e+00_EB,          &
.000e+00_EB, .000e+00_EB, .000e+00_EB, .420e-01_EB, .700e-01_EB, .108e+00_EB,         & 
.000e+00_EB, .000e+00_EB, .000e+00_EB, .360e-01_EB, .640e-01_EB, .980e-01_EB,        &  
.000e+00_EB, .000e+00_EB, .000e+00_EB, .350e-01_EB, .610e-01_EB, .870e-01_EB,       &   
.000e+00_EB, .000e+00_EB, .000e+00_EB, .320e-01_EB, .580e-01_EB, .860e-01_EB/),    &    
(/6,8/))
sd15(1:6,73:80) = RESHAPE ((/ & ! 860-880
.000e+00_EB, .000e+00_EB, .000e+00_EB, .330e-01_EB, .560e-01_EB, .750e-01_EB,     &     
.000e+00_EB, .000e+00_EB, .000e+00_EB, .300e-01_EB, .530e-01_EB, .750e-01_EB,    &      
.000e+00_EB, .000e+00_EB, .000e+00_EB, .290e-01_EB, .530e-01_EB, .850e-01_EB,   &       
.000e+00_EB, .000e+00_EB, .000e+00_EB, .240e-01_EB, .470e-01_EB, .900e-01_EB,  &        
.000e+00_EB, .000e+00_EB, .000e+00_EB, .220e-01_EB, .450e-01_EB, .860e-01_EB, &         
.000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB,&
.000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB,&
.000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB, .000e+00_EB/), &
(/6,8/))

! initialize sd7 array

!   the following data are for the 7.7 micron band of ch4
! temp,k= 290 600   850

sd7(1:3,1:8) = RESHAPE ((/ &
0._EB, 0._EB, 0._EB,&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.22_EB,&
0.16_EB, 0.20_EB, 0.47_EB,&
0.34_EB, 0.34_EB, 0.62_EB,&
0.69_EB, 0.53_EB, 0.65_EB,&
1.27_EB, 0.88_EB, 1.09_EB,&
1.68_EB, 1.38_EB, 0.87_EB/),(/3,8/))
sd7(1:3,9:16) = RESHAPE ((/&
0.55_EB, 0.28_EB, 0.40_EB,&
1.25_EB, 0.86_EB, 0.93_EB,&
0.34_EB, 0.59_EB, 0.75_EB,&
0._EB, 0.13_EB, 0.25_EB,&
0._EB, 0._EB, 0.06_EB,&
0._EB, 0._EB, 0._EB,&
0._EB, 0._EB, 0._EB,&
0._EB, 0._EB, 0._EB/),(/3,8/))

! initialize sd3 array

!  the following data are for the 3.3 micron band of ch4
! temp, k= 290  600   850
sd3(1:3,1:8) = RESHAPE ((/&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.06_EB,&
0.03_EB, 0.03_EB, 0.09_EB,&
0.07_EB, 0.07_EB, 0.12_EB,&
0.09_EB, 0.09_EB, 0.12_EB,&
0.14_EB, 0.15_EB, 0.22_EB/),(/3,8/))
sd3(1:3,9:16) = RESHAPE ((/&
0.18_EB, 0.22_EB, 0.28_EB,&
0.24_EB, 0.31_EB, 0.37_EB,&
0.33_EB, 0.44_EB, 0.47_EB,&
0.45_EB, 0.50_EB, 0.53_EB,&
0.59_EB, 0.62_EB, 0.62_EB,&
0.74_EB, 0.70_EB, 0.68_EB,&
0.91_EB, 0.77_EB, 0.72_EB,&
1.00_EB, 0.81_EB, 0.75_EB/),(/3,8/))
sd3(1:3,17:24) = RESHAPE ((/&
1.03_EB, 0.84_EB, 0.78_EB,&
1.03_EB, 0.84_EB, 0.78_EB,&
1.00_EB, 0.81_EB, 0.75_EB,&
0.94_EB, 0.77_EB, 0.72_EB,&
0.72_EB, 0.68_EB, 0.68_EB,&
0.52_EB, 0.63_EB, 0.63_EB,&
0.33_EB, 0.50_EB, 0.56_EB,&
0.25_EB, 0.42_EB, 0.50_EB/),(/3,8/))
sd3(1:3,25:32) = RESHAPE ((/&
0.17_EB, 0.26_EB, 0.37_EB,&
0.08_EB, 0.18_EB, 0.31_EB,&
0.04_EB, 0.11_EB, 0.22_EB,&
0._EB, 0.06_EB, 0.16_EB,&
0._EB, 0.02_EB, 0.12_EB,&
0._EB, 0._EB, 0.06_EB,&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0._EB/),(/3,8/))
!------------------------------------------------------------------------------
END SUBROUTINE RCALLOC


!==============================================================================
SUBROUTINE RCDEALLOC2
!==============================================================================   
DEALLOCATE(incident_radiance)
DEALLOCATE(ttau)
DEALLOCATE(lambda)
DEALLOCATE(ab)
DEALLOCATE(wave_number)
   
END SUBROUTINE RCDEALLOC2

!==============================================================================
SUBROUTINE RCDEALLOC
!==============================================================================   

DEALLOCATE(gamma)
DEALLOCATE(sd15)
DEALLOCATE(sd)
DEALLOCATE(sd7)
DEALLOCATE(sd3)

!-------------------------deallocation methane variables-------------------
DEALLOCATE(sd_ch4_temp)
DEALLOCATE(om_bnd_ch4)
DEALLOCATE(sd1_ch4)
DEALLOCATE(sd2_ch4)

!-------------------------deallocation propane variables-------------------
DEALLOCATE(sd_c3h8_temp)
DEALLOCATE(be_c3h8)
DEALLOCATE(om_bnd_c3h8)
DEALLOCATE(sd1_c3h8)
DEALLOCATE(sd2_c3h8)
DEALLOCATE(gammad1_c3h8)
DEALLOCATE(gammad2_c3h8)

!-------------------------deallocation methanol variables-------------------
DEALLOCATE(sd_ch3oh_temp)
DEALLOCATE(be_ch3oh)
DEALLOCATE(om_bnd_ch3oh)
DEALLOCATE(sd1_ch3oh)
DEALLOCATE(sd2_ch3oh)
DEALLOCATE(sd3_ch3oh)
DEALLOCATE(sd4_ch3oh)
DEALLOCATE(gammad1_ch3oh)
DEALLOCATE(gammad2_ch3oh)
DEALLOCATE(gammad3_ch3oh)
DEALLOCATE(gammad4_ch3oh)

!-------------------------deallocation heptane variables-------------------
DEALLOCATE(sd_c7h16_temp)
DEALLOCATE(be_c7h16)
DEALLOCATE(om_bnd_c7h16)
DEALLOCATE(sd1_c7h16)
DEALLOCATE(sd2_c7h16)
DEALLOCATE(gammad1_c7h16)
DEALLOCATE(gammad2_c7h16)

!-------------------------deallocation toluene variables-------------------
DEALLOCATE(sd_c7h8_temp)
DEALLOCATE(be_c7h8)
DEALLOCATE(om_bnd_c7h8)
DEALLOCATE(sd1_c7h8)
DEALLOCATE(sd2_c7h8)
DEALLOCATE(sd3_c7h8)
DEALLOCATE(sd4_c7h8)
DEALLOCATE(sd5_c7h8)
DEALLOCATE(gammad1_c7h8)
DEALLOCATE(gammad2_c7h8)
DEALLOCATE(gammad3_c7h8)
DEALLOCATE(gammad4_c7h8)
DEALLOCATE(gammad5_c7h8)

!-------------------------deallocation propylene variables-------------------
DEALLOCATE(sd_c3h6_temp)
DEALLOCATE(be_c3h6)
DEALLOCATE(om_bnd_c3h6)
DEALLOCATE(sd1_c3h6)
DEALLOCATE(sd2_c3h6)
DEALLOCATE(sd3_c3h6)
DEALLOCATE(gammad1_c3h6)
DEALLOCATE(gammad2_c3h6)
DEALLOCATE(gammad3_c3h6)

!-------------------------deallocation mma variables-------------------
DEALLOCATE(sd_c5h8o2_temp)
DEALLOCATE(be_c5h8o2)
DEALLOCATE(om_bnd_c5h8o2)
DEALLOCATE(sd1_c5h8o2)
DEALLOCATE(sd2_c5h8o2)
DEALLOCATE(sd3_c5h8o2)
DEALLOCATE(sd4_c5h8o2)
DEALLOCATE(sd5_c5h8o2)
DEALLOCATE(sd6_c5h8o2)
DEALLOCATE(gammad1_c5h8o2)
DEALLOCATE(gammad2_c5h8o2)
DEALLOCATE(gammad3_c5h8o2)
DEALLOCATE(gammad4_c5h8o2)
DEALLOCATE(gammad5_c5h8o2)
DEALLOCATE(gammad6_c5h8o2)

!-------------------------deallocation ethane variables-------------------
DEALLOCATE(sd_c2h6_temp)
DEALLOCATE(be_c2h6)
DEALLOCATE(om_bnd_c2h6)
DEALLOCATE(sd1_c2h6)
DEALLOCATE(sd2_c2h6)
DEALLOCATE(sd3_c2h6)
DEALLOCATE(gammad1_c2h6)
DEALLOCATE(gammad2_c2h6)
DEALLOCATE(gammad3_c2h6)

!-------------------------deallocation ethylene variables-------------------
DEALLOCATE(sd_c2h4_temp)
DEALLOCATE(be_c2h4)
DEALLOCATE(om_bnd_c2h4)
DEALLOCATE(sd1_c2h4)
DEALLOCATE(sd2_c2h4)
DEALLOCATE(sd3_c2h4)
DEALLOCATE(sd4_c2h4)
DEALLOCATE(gammad1_c2h4)
DEALLOCATE(gammad2_c2h4)
DEALLOCATE(gammad3_c2h4)
DEALLOCATE(gammad4_c2h4)

!------------------------------------------------------------------------------   
END SUBROUTINE RCDEALLOC

END MODULE RADCAL_CALC


MODULE SPECDATA

USE PRECISION_PARAMETERS
IMPLICIT NONE
INTEGER ::  I,J
INTEGER, PARAMETER :: NWATERK=183
REAL(EB) :: CPLXREF_WATER(NWATERK,2)
INTEGER, PARAMETER :: NFUELK=94
REAL(EB) :: CPLXREF_FUEL(NFUELK,3)

! Hale, G.M. and Querry, M.R., "Optical Constants of Water in the 200-nm to 200-\mu m Wavelength Region", 
! Applied Optics, 12(3), 555-63 (1973)

DATA ((CPLXREF_WATER(I,J), J = 1,2), I = 1,10)  &
/ 1.0000000E-6_EB,  2.8647890E-6_EB,&
1.0200000E-6_EB,  2.1915636E-6_EB,&
1.0400000E-6_EB,  1.3241691E-6_EB,&
1.0600000E-6_EB,  1.0122254E-6_EB,&
1.0800000E-6_EB,  1.1172677E-6_EB,&
1.1000000E-6_EB,  1.4880987E-6_EB,&
1.1200000E-6_EB,  4.6345919E-6_EB,&
1.1400000E-6_EB,  5.9874090E-6_EB,&
1.1600000E-6_EB,  8.2155782E-6_EB,&
1.1800000E-6_EB,  9.7657473E-6_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =11,20)  &
/ 1.2000000E-6_EB,  9.9312684E-6_EB,&
1.2200000E-6_EB,  9.2230290E-6_EB,&
1.2400000E-6_EB,  8.6834937E-6_EB,&
1.2600000E-6_EB,  8.9238177E-6_EB,&
1.2800000E-6_EB,  9.9821980E-6_EB,&
1.3000000E-6_EB,  1.1483029E-5_EB,&
1.3200000E-6_EB,  1.1953809E-5_EB,&
1.3400000E-6_EB,  1.2614780E-5_EB,&
1.3600000E-6_EB,  2.9978533E-5_EB,&
1.3800000E-6_EB,  6.6988316E-5_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =21,30)  &
/ 1.4000000E-6_EB,  1.3803508E-4_EB,&
1.4200000E-6_EB,  2.4995602E-4_EB,&
1.4400000E-6_EB,  3.3002369E-4_EB,&
1.4600000E-6_EB,  3.2996003E-4_EB,&
1.4800000E-6_EB,  2.5003560E-4_EB,&
1.5000000E-6_EB,  2.0996516E-4_EB,&
1.5200000E-6_EB,  1.6994565E-4_EB,&
1.5400000E-6_EB,  1.4497583E-4_EB,&
1.5600000E-6_EB,  1.2004433E-4_EB,&
1.5800000E-6_EB,  9.9957388E-5_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =31,40)  &
/ 1.6000000E-6_EB,  8.5561825E-5_EB,&
1.6200000E-6_EB,  7.5028952E-5_EB,&
1.6400000E-6_EB,  6.4992643E-5_EB,&
1.6600000E-6_EB,  5.9972898E-5_EB,&
1.6800000E-6_EB,  6.0027012E-5_EB,&
1.7000000E-6_EB,  6.0065211E-5_EB,&
1.7200000E-6_EB,  6.9942231E-5_EB,&
1.7400000E-6_EB,  8.5017388E-5_EB,&
1.7600000E-6_EB,  1.0000023E-4_EB,&
1.7800000E-6_EB,  1.1501809E-4_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =41,50)  &
/ 1.8000000E-6_EB,  1.1502128E-4_EB,&
1.8200000E-6_EB,  1.3005838E-4_EB,&
1.8400000E-6_EB,  1.4993669E-4_EB,&
1.8600000E-6_EB,  2.1003200E-4_EB,&
1.8800000E-6_EB,  4.6497435E-4_EB,&
1.9000000E-6_EB,  1.0000183E-3_EB,&
1.9200000E-6_EB,  1.7500423E-3_EB,&
1.9400000E-6_EB,  1.8499391E-3_EB,&
1.9600000E-6_EB,  1.6500261E-3_EB,&
1.9800000E-6_EB,  1.4500559E-3_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =51,60)  &
/ 2.0000000E-6_EB,  1.1000790E-3_EB,&
2.0200000E-6_EB,  9.0001961E-4_EB,&
2.0400000E-6_EB,  7.3003417E-4_EB,&
2.0600000E-6_EB,  6.3998112E-4_EB,&
2.0800000E-6_EB,  5.2006742E-4_EB,&
2.1000000E-6_EB,  4.5003447E-4_EB,&
2.1200000E-6_EB,  4.0505888E-4_EB,&
2.1400000E-6_EB,  3.4995785E-4_EB,&
2.1600000E-6_EB,  3.2005422E-4_EB,&
2.1800000E-6_EB,  2.9994500E-4_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =61,70)  &
/ 2.2000000E-6_EB,  2.8904129E-4_EB,&
2.2200000E-6_EB,  2.8495578E-4_EB,&
2.2400000E-6_EB,  2.9500960E-4_EB,&
2.2600000E-6_EB,  3.1005293E-4_EB,&
2.2800000E-6_EB,  3.5997028E-4_EB,&
2.3000000E-6_EB,  4.0998313E-4_EB,&
2.3200000E-6_EB,  4.9496551E-4_EB,&
2.3400000E-6_EB,  5.9494505E-4_EB,&
2.3600000E-6_EB,  6.9994116E-4_EB,&
2.3800000E-6_EB,  8.2007768E-4_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =71,80)  &
/ 2.4000000E-6_EB,  9.5607557E-4_EB,&
2.4200000E-6_EB,  1.1500727E-3_EB,&
2.4400000E-6_EB,  1.2999617E-3_EB,&
2.4600000E-6_EB,  1.4999176E-3_EB,&
2.4800000E-6_EB,  1.6999912E-3_EB,&
2.5000000E-6_EB,  1.8000424E-3_EB,&
2.5200000E-6_EB,  2.0500716E-3_EB,&
2.5400000E-6_EB,  2.1999478E-3_EB,&
2.5600000E-6_EB,  2.3500946E-3_EB,&
2.5800000E-6_EB,  2.7000302E-3_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =81,90)  &
/ 2.6000000E-6_EB,  3.1699367E-3_EB,&
2.6500000E-6_EB,  6.7000889E-3_EB,&
2.7000000E-6_EB,  1.8999997E-2_EB,&
2.7500000E-6_EB,  5.9000926E-2_EB,&
2.8000000E-6_EB,  1.1500027E-1_EB,&
2.8500000E-6_EB,  1.8499960E-1_EB,&
2.9000000E-6_EB,  2.6769861E-1_EB,&
2.9500000E-6_EB,  2.9813700E-1_EB,&
3.0000000E-6_EB,  2.7215495E-1_EB,&
3.0500000E-6_EB,  2.4000020E-1_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =91,100)  &
/ 3.1000000E-6_EB,  1.9200142E-1_EB,&
3.1500000E-6_EB,  1.3500032E-1_EB,&
3.2000000E-6_EB,  9.2401540E-2_EB,&
3.2500000E-6_EB,  6.0999713E-2_EB,&
3.3000000E-6_EB,  3.6798931E-2_EB,&
3.3500000E-6_EB,  2.6099958E-2_EB,&
3.4000000E-6_EB,  1.9500046E-2_EB,&
3.4500000E-6_EB,  1.3199993E-2_EB,&
3.5000000E-6_EB,  9.4000888E-3_EB,&
3.6000000E-6_EB,  5.1500311E-3_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =101,110)  &
/ 3.7000000E-6_EB,  3.6000769E-3_EB,&
3.8000000E-6_EB,  3.4001225E-3_EB,&
3.9000000E-6_EB,  3.7999516E-3_EB,&
4.0000000E-6_EB,  4.5998962E-3_EB,&
4.1000000E-6_EB,  5.6118033E-3_EB,&
4.2000000E-6_EB,  6.8850428E-3_EB,&
4.3000000E-6_EB,  8.4519233E-3_EB,&
4.4000000E-6_EB,  1.0294142E-2_EB,&
4.5000000E-6_EB,  1.3392888E-2_EB,&
4.6000000E-6_EB,  1.4715466E-2_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =111,120)  &
/ 4.7000000E-6_EB,  1.5708593E-2_EB,&
4.8000000E-6_EB,  1.5011494E-2_EB,&
4.9000000E-6_EB,  1.3686529E-2_EB,&
5.0000000E-6_EB,  1.2414086E-2_EB,&
5.1000000E-6_EB,  1.1120156E-2_EB,&
5.2000000E-6_EB,  1.0096790E-2_EB,&
5.3000000E-6_EB,  9.7848459E-3_EB,&
5.4000000E-6_EB,  1.0313240E-2_EB,&
5.5000000E-6_EB,  1.1598416E-2_EB,&
5.6000000E-6_EB,  1.4215720E-2_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =121,130)  &
/ 5.7000000E-6_EB,  2.0320903E-2_EB,&
5.8000000E-6_EB,  3.3000777E-2_EB,&
5.9000000E-6_EB,  6.2209688E-2_EB,&
6.0000000E-6_EB,  1.0699987E-1_EB,&
6.1000000E-6_EB,  1.3101555E-1_EB,&
6.2000000E-6_EB,  8.8019050E-2_EB,&
6.3000000E-6_EB,  5.7002139E-2_EB,&
6.4000000E-6_EB,  4.4868962E-2_EB,&
6.5000000E-6_EB,  3.9207820E-2_EB,&
6.6000000E-6_EB,  3.5609327E-2_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =131,140)  &
/ 6.7000000E-6_EB,  3.3696285E-2_EB,&
6.8000000E-6_EB,  3.2684059E-2_EB,&
6.9000000E-6_EB,  3.2176355E-2_EB,&
7.0000000E-6_EB,  3.1974228E-2_EB,&
7.1000000E-6_EB,  3.1979003E-2_EB,&
7.2000000E-6_EB,  3.2085637E-2_EB,&
7.3000000E-6_EB,  3.2182721E-2_EB,&
7.4000000E-6_EB,  3.2388031E-2_EB,&
7.5000000E-6_EB,  3.2586975E-2_EB,&
7.6000000E-6_EB,  3.2779552E-2_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =141,150)  &
/ 7.7000000E-6_EB,  3.3088313E-2_EB,&
7.8000000E-6_EB,  3.3518031E-2_EB,&
7.9000000E-6_EB,  3.3884883E-2_EB,&
8.0000000E-6_EB,  3.4313806E-2_EB,&
8.2000000E-6_EB,  3.5106397E-2_EB,&
8.4000000E-6_EB,  3.6096341E-2_EB,&
8.6000000E-6_EB,  3.1001791E-2_EB,&
8.8000000E-6_EB,  3.8515496E-2_EB,&
9.0000000E-6_EB,  3.9892186E-2_EB,&
9.2000000E-6_EB,  4.1510792E-2_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =151,160)  &
/ 9.4000000E-6_EB,  4.3310835E-2_EB,&
9.6000000E-6_EB,  4.5378257E-2_EB,&
9.8000000E-6_EB,  4.7883356E-2_EB,&
1.0000000E-5_EB,  5.0770427E-2_EB,&
1.0500000E-5_EB,  6.6176625E-2_EB,&
1.1000000E-5_EB,  9.6813952E-2_EB,&
1.1500000E-5_EB,  1.4202987E-1_EB,&
1.2000000E-5_EB,  1.9900734E-1_EB,&
1.2500000E-5_EB,  2.5902467E-1_EB,&
1.3000000E-5_EB,  3.0497270E-1_EB / 
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =161,170)  &
/ 1.3500000E-5_EB,  3.4302267E-1_EB,&
1.4000000E-5_EB,  3.6998750E-1_EB,&
1.4500000E-5_EB,  3.8804760E-1_EB,&
1.5000000E-5_EB,  4.0202539E-1_EB,&
1.5500000E-5_EB,  4.1394609E-1_EB,&
1.6000000E-5_EB,  4.2195159E-1_EB,&
1.6500000E-5_EB,  4.2804722E-1_EB,&
1.7000000E-5_EB,  4.2897828E-1_EB,&
1.7500000E-5_EB,  4.2906183E-1_EB,&
1.8000000E-5_EB,  4.2599412E-1_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =171,180)  &
/ 1.8500000E-5_EB,  4.2104440E-1_EB,&
1.9000000E-5_EB,  4.1397792E-1_EB,&
1.9500000E-5_EB,  4.0407849E-1_EB,&
2.0000000E-5_EB,  3.9295355E-1_EB,&
3.0000000E-5_EB,  3.2801834E-1_EB,&
3.8000000E-5_EB,  3.6105890E-1_EB,&
5.0000000E-5_EB,  5.1407047E-1_EB,&
6.0000000E-5_EB,  5.8680428E-1_EB,&
7.0000000E-5_EB,  5.7598174E-1_EB,&
8.0000000E-5_EB,  5.4685638E-1_EB /
DATA ((CPLXREF_WATER(I,J), J = 1,2), I =181,NWATERK)  &
/ 9.0000000E-5_EB,  5.3571554E-1_EB,&
1.0000000E-4_EB,  5.3237328E-1_EB,&
2.0000000E-4_EB,  5.0452117E-1_EB /

! Heptane properties from
! L.A. Dombrovsky, S.S. Sazhin, S.V. Mikhalovsky, R. Wood, M.R. Heikal
! Spectral properties of diesel fuel PARTICLEs
! Fuel, Vol. 82, No. 1 (2003) pp. 15-22
DATA (CPLXREF_FUEL( 1,J), J=1,3) / 0.7_EB,   1.47_EB, 1.55E-07_EB /
DATA (CPLXREF_FUEL( 2,J), J=1,3) / 0.8_EB,   1.47_EB, 4.33E-07_EB /
DATA (CPLXREF_FUEL( 3,J), J=1,3) / 0.9_EB,   1.47_EB, 9.39E-07_EB /
DATA (CPLXREF_FUEL( 4,J), J=1,3) / 1.0_EB,   1.47_EB, 8.81E-07_EB /
DATA (CPLXREF_FUEL( 5,J), J=1,3) / 1.1_EB,   1.47_EB, 1.19E-06_EB /
DATA (CPLXREF_FUEL( 6,J), J=1,3) / 1.2_EB,   1.47_EB, 2.63E-06_EB /
DATA (CPLXREF_FUEL( 7,J), J=1,3) / 1.3_EB,   1.47_EB, 1.08E-05_EB /
DATA (CPLXREF_FUEL( 8,J), J=1,3) / 1.4_EB,   1.47_EB, 2.24E-05_EB /
DATA (CPLXREF_FUEL( 9,J), J=1,3) / 1.5_EB,   1.46_EB, 3.41E-05_EB /
DATA (CPLXREF_FUEL( 10,J), J=1,3) / 1.6_EB,  1.46_EB, 4.57E-05_EB /
DATA (CPLXREF_FUEL( 11,J), J=1,3) / 1.7_EB,  1.46_EB, 5.73E-05_EB /
DATA (CPLXREF_FUEL( 12,J), J=1,3) / 1.8_EB,  1.46_EB, 6.90E-05_EB /
DATA (CPLXREF_FUEL( 13,J), J=1,3) / 1.9_EB,  1.46_EB, 8.06E-05_EB /
DATA (CPLXREF_FUEL( 14,J), J=1,3) / 2.0_EB,  1.46_EB, 9.22E-05_EB /
DATA (CPLXREF_FUEL( 15,J), J=1,3) / 2.1_EB,  1.46_EB, 1.04E-04_EB /
DATA (CPLXREF_FUEL( 16,J), J=1,3) / 2.2_EB,  1.46_EB, 1.15E-04_EB /
DATA (CPLXREF_FUEL( 17,J), J=1,3) / 2.3_EB,  1.46_EB, 1.27E-04_EB /
DATA (CPLXREF_FUEL( 18,J), J=1,3) / 2.4_EB,  1.45_EB, 1.39E-04_EB /
DATA (CPLXREF_FUEL( 19,J), J=1,3) / 2.5_EB,  1.45_EB, 1.50E-04_EB /
DATA (CPLXREF_FUEL( 20,J), J=1,3) / 2.6_EB,  1.45_EB, 1.62E-04_EB /
DATA (CPLXREF_FUEL( 21,J), J=1,3) / 2.7_EB,  1.45_EB, 1.11E-04_EB /
DATA (CPLXREF_FUEL( 22,J), J=1,3) / 2.8_EB,  1.44_EB, 5.92E-05_EB /
DATA (CPLXREF_FUEL( 23,J), J=1,3) / 2.9_EB,  1.44_EB, 7.45E-05_EB /
DATA (CPLXREF_FUEL( 24,J), J=1,3) / 3.0_EB,  1.43_EB, 9.72E-05_EB /
DATA (CPLXREF_FUEL( 25,J), J=1,3) / 3.1_EB,  1.42_EB, 3.12E-04_EB /
DATA (CPLXREF_FUEL( 26,J), J=1,3) / 3.2_EB,  1.40_EB, 6.09E-04_EB /
DATA (CPLXREF_FUEL( 27,J), J=1,3) / 3.3_EB,  1.17_EB, 5.72E-02_EB /
DATA (CPLXREF_FUEL( 28,J), J=1,3) / 3.4_EB,  1.39_EB, 1.20E-01_EB /
DATA (CPLXREF_FUEL( 29,J), J=1,3) / 3.5_EB,  1.45_EB, 8.24E-02_EB /
DATA (CPLXREF_FUEL( 30,J), J=1,3) / 3.6_EB,  1.51_EB, 1.63E-03_EB /
DATA (CPLXREF_FUEL( 31,J), J=1,3) / 3.7_EB,  1.64_EB, 1.33E-03_EB /
DATA (CPLXREF_FUEL( 32,J), J=1,3) / 3.8_EB,  1.56_EB, 1.02E-03_EB /
DATA (CPLXREF_FUEL( 33,J), J=1,3) / 3.9_EB,  1.52_EB, 6.36E-04_EB /
DATA (CPLXREF_FUEL( 34,J), J=1,3) / 4.0_EB,  1.50_EB, 2.51E-04_EB /
DATA (CPLXREF_FUEL( 35,J), J=1,3) / 4.1_EB,  1.50_EB, 2.59E-04_EB /
DATA (CPLXREF_FUEL( 36,J), J=1,3) / 4.2_EB,  1.49_EB, 3.10E-04_EB /
DATA (CPLXREF_FUEL( 37,J), J=1,3) / 4.3_EB,  1.49_EB, 2.60E-04_EB /
DATA (CPLXREF_FUEL( 38,J), J=1,3) / 4.4_EB,  1.48_EB, 2.11E-04_EB /
DATA (CPLXREF_FUEL( 39,J), J=1,3) / 4.5_EB,  1.48_EB, 1.98E-04_EB /
DATA (CPLXREF_FUEL( 40,J), J=1,3) / 4.6_EB,  1.47_EB, 1.90E-04_EB /
DATA (CPLXREF_FUEL( 41,J), J=1,3) / 4.7_EB,  1.47_EB, 1.58E-04_EB /
DATA (CPLXREF_FUEL( 42,J), J=1,3) / 4.8_EB,  1.47_EB, 1.21E-04_EB /
DATA (CPLXREF_FUEL( 43,J), J=1,3) / 4.9_EB,  1.47_EB, 8.32E-05_EB /
DATA (CPLXREF_FUEL( 44,J), J=1,3) / 5.0_EB,  1.46_EB, 4.58E-05_EB /
DATA (CPLXREF_FUEL( 45,J), J=1,3) / 5.1_EB,  1.46_EB, 6.37E-05_EB /
DATA (CPLXREF_FUEL( 46,J), J=1,3) / 5.2_EB,  1.46_EB, 7.78E-05_EB /
DATA (CPLXREF_FUEL( 47,J), J=1,3) / 5.3_EB,  1.46_EB, 7.66E-05_EB /
DATA (CPLXREF_FUEL( 48,J), J=1,3) / 5.4_EB,  1.45_EB, 7.55E-05_EB /
DATA (CPLXREF_FUEL( 49,J), J=1,3) / 5.5_EB,  1.45_EB, 9.62E-05_EB /
DATA (CPLXREF_FUEL( 50,J), J=1,3) / 5.6_EB,  1.45_EB, 1.17E-04_EB /
DATA (CPLXREF_FUEL( 51,J), J=1,3) / 5.7_EB,  1.45_EB, 1.63E-04_EB /
DATA (CPLXREF_FUEL( 52,J), J=1,3) / 5.8_EB,  1.45_EB, 2.16E-04_EB /
DATA (CPLXREF_FUEL( 53,J), J=1,3) / 5.9_EB,  1.44_EB, 2.68E-04_EB /
DATA (CPLXREF_FUEL( 54,J), J=1,3) / 6.0_EB,  1.44_EB, 3.21E-04_EB /
DATA (CPLXREF_FUEL( 55,J), J=1,3) / 6.1_EB,  1.44_EB, 3.74E-04_EB /
DATA (CPLXREF_FUEL( 56,J), J=1,3) / 6.2_EB,  1.44_EB, 4.36E-04_EB /
DATA (CPLXREF_FUEL( 57,J), J=1,3) / 6.3_EB,  1.43_EB, 5.87E-04_EB /
DATA (CPLXREF_FUEL( 58,J), J=1,3) / 6.4_EB,  1.43_EB, 7.38E-04_EB /
DATA (CPLXREF_FUEL( 59,J), J=1,3) / 6.5_EB,  1.42_EB, 1.35E-03_EB /
DATA (CPLXREF_FUEL( 60,J), J=1,3) / 6.6_EB,  1.41_EB, 6.12E-03_EB /
DATA (CPLXREF_FUEL( 61,J), J=1,3) / 6.7_EB,  1.39_EB, 2.06E-02_EB /
DATA (CPLXREF_FUEL( 62,J), J=1,3) / 6.8_EB,  1.35_EB, 3.51E-02_EB /
DATA (CPLXREF_FUEL( 63,J), J=1,3) / 6.9_EB,  1.37_EB, 2.29E-02_EB /
DATA (CPLXREF_FUEL( 64,J), J=1,3) / 7.0_EB,  1.48_EB, 3.99E-03_EB /
DATA (CPLXREF_FUEL( 65,J), J=1,3) / 7.1_EB,  1.53_EB, 3.24E-03_EB /
DATA (CPLXREF_FUEL( 66,J), J=1,3) / 7.2_EB,  1.48_EB, 2.61E-03_EB /
DATA (CPLXREF_FUEL( 67,J), J=1,3) / 7.3_EB,  1.43_EB, 2.97E-03_EB /
DATA (CPLXREF_FUEL( 68,J), J=1,3) / 7.4_EB,  1.47_EB, 3.33E-03_EB /
DATA (CPLXREF_FUEL( 69,J), J=1,3) / 7.5_EB,  1.50_EB, 2.80E-03_EB /
DATA (CPLXREF_FUEL( 70,J), J=1,3) / 7.6_EB,  1.49_EB, 2.14E-03_EB /
DATA (CPLXREF_FUEL( 71,J), J=1,3) / 7.7_EB,  1.48_EB, 2.28E-03_EB /
DATA (CPLXREF_FUEL( 72,J), J=1,3) / 7.8_EB,  1.47_EB, 2.41E-03_EB /
DATA (CPLXREF_FUEL( 73,J), J=1,3) / 7.9_EB,  1.47_EB, 1.63E-03_EB /
DATA (CPLXREF_FUEL( 74,J), J=1,3) / 8.0_EB,  1.47_EB, 9.84E-04_EB /
DATA (CPLXREF_FUEL( 75,J), J=1,3) / 8.1_EB,  1.46_EB, 9.03E-04_EB /
DATA (CPLXREF_FUEL( 76,J), J=1,3) / 8.2_EB,  1.46_EB, 8.19E-04_EB /
DATA (CPLXREF_FUEL( 77,J), J=1,3) / 8.3_EB,  1.46_EB, 7.13E-04_EB /
DATA (CPLXREF_FUEL( 78,J), J=1,3) / 8.4_EB,  1.46_EB, 6.07E-04_EB /
DATA (CPLXREF_FUEL( 79,J), J=1,3) / 8.5_EB,  1.45_EB, 8.54E-04_EB /
DATA (CPLXREF_FUEL( 80,J), J=1,3) / 8.6_EB,  1.45_EB, 1.10E-03_EB /
DATA (CPLXREF_FUEL( 81,J), J=1,3) / 8.7_EB,  1.45_EB, 1.76E-03_EB /
DATA (CPLXREF_FUEL( 82,J), J=1,3) / 8.8_EB,  1.45_EB, 2.41E-03_EB /
DATA (CPLXREF_FUEL( 83,J), J=1,3) / 8.9_EB,  1.45_EB, 1.58E-03_EB /
DATA (CPLXREF_FUEL( 84,J), J=1,3) / 9.0_EB,  1.45_EB, 5.71E-04_EB /
DATA (CPLXREF_FUEL( 85,J), J=1,3) / 9.1_EB,  1.45_EB, 8.53E-04_EB /
DATA (CPLXREF_FUEL( 86,J), J=1,3) / 9.2_EB,  1.45_EB, 1.28E-03_EB /
DATA (CPLXREF_FUEL( 87,J), J=1,3) / 9.3_EB,  1.45_EB, 1.63E-03_EB /
DATA (CPLXREF_FUEL( 88,J), J=1,3) / 9.4_EB,  1.45_EB, 1.98E-03_EB /
DATA (CPLXREF_FUEL( 89,J), J=1,3) / 9.5_EB,  1.45_EB, 1.69E-03_EB /
DATA (CPLXREF_FUEL( 90,J), J=1,3) / 9.6_EB,  1.45_EB, 1.24E-03_EB /
DATA (CPLXREF_FUEL( 91,J), J=1,3) / 9.7_EB,  1.45_EB, 1.43E-03_EB /
DATA (CPLXREF_FUEL( 92,J), J=1,3) / 9.8_EB,  1.45_EB, 1.69E-03_EB /
DATA (CPLXREF_FUEL( 93,J), J=1,3) / 9.9_EB,  1.45_EB, 1.26E-03_EB /
DATA (CPLXREF_FUEL( 94,J), J=1,3) / 10.0_EB, 1.45_EB, 8.24E-04_EB /

END MODULE SPECDATA



MODULE MIEV

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS, ONLY: NUMBER_SPECTRAL_BANDS,NUMBER_RADIATION_ANGLES,LU_ERR
USE SPECDATA, ONLY: CPLXREF_WATER, NWATERK, CPLXREF_FUEL, NFUELK
USE RADCAL_CALC, ONLY: PLANCK
USE RADCONS
USE MEMORY_FUNCTIONS, ONLY : CHKMEMERR
USE MATH_FUNCTIONS, ONLY : INTERPOLATE1D

IMPLICIT NONE
REAL(EB), ALLOCATABLE :: RDMIE(:),  LMBDMIE(:),LMBDWGHT(:),REAL_REF_INDX(:),CMPLX_REF_INDX(:)
REAL(EB), ALLOCATABLE :: QSCA(:,:), QABS(:,:),  CHI_F(:,:)

PRIVATE 
PUBLIC MEAN_CROSS_SECTIONS

CONTAINS


SUBROUTINE MEAN_CROSS_SECTIONS(CLASS_NUMBER)

! This subroutine calculates the mean scattering and absorption
! coefficients for each radiation band and PARTICLE size group.

USE TYPES, ONLY: LAGRANGIAN_PARTICLE_CLASS_TYPE, LAGRANGIAN_PARTICLE_CLASS, TABLES_TYPE, TABLES
USE GLOBAL_CONSTANTS, ONLY: H2O_INDEX
USE MATH_FUNCTIONS, ONLY : INTERPOLATE1D
INTEGER  :: NSB,I,J,IBND,IZERO,NX,NLAMBDALOW(1),NLAMBDAHIGH(1),ND,CLASS_NUMBER,NL_DATA
REAL(EB) :: RMMAX,RMMIN,RDTMP,IB,IBSUM,AVAL,BVAL,ASUM,BSUM,B_WIEN
TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), POINTER :: LPC=>NULL()
TYPE (TABLES_TYPE),  POINTER :: TA=>NULL()
 
LPC => LAGRANGIAN_PARTICLE_CLASS(CLASS_NUMBER) 
NSB = NUMBER_SPECTRAL_BANDS

! Physical parameters
! minimum mean radius (m)
RMMIN = 0.5_EB*1.E-6*MIE_MINIMUM_DIAMETER
IF (RMMIN < EPSILON_EB) RMMIN = 0.5E-6_EB   

! maximum mean radius (m)
RMMAX = 0.5_EB*1.E-6*MIE_MAXIMUM_DIAMETER
IF (RMMAX < EPSILON_EB) THEN  
   RMMAX = 0.5_EB*LPC%DIAMETER
   ! Allow increase of the mean radius
   RMMAX = 1.5_EB*RMMAX
ENDIF

! Other constants

B_WIEN = 2.8977685E-3_EB
 
! Calculate parameters of the PARTICLE group lookup table
 
DGROUP_A = (LOG(RMMAX)-LOG(RMMIN))/(MIE_NDG-1)
DGROUP_B = LOG(RMMAX)-DGROUP_A*MIE_NDG

! Generate the PARTICLE radii for mie table (microns)
 
RDTMP = 0.5_EB*RMMIN
NX = 0
DO WHILE (RDTMP < 2._EB*RMMAX) 
   NX = NX + 1
   RDTMP = RDTMP + MIN(1.0E6_EB,0.2_EB*RDTMP**(1._EB))
ENDDO
NRDMIE = NX 

ALLOCATE(RDMIE(1:NRDMIE),STAT=IZERO)
CALL ChkMemErr('MIEV','RDMIE',IZERO)

RDTMP = 0.5_EB*RMMIN
RDMIE(1) = RDTMP
DO NX = 2, NRDMIE
   RDTMP = RDTMP + MIN(3.0E6_EB,0.2_EB*RDTMP**(1.0_EB))
   RDMIE(NX) = RDTMP
ENDDO

! Radiative properties

NLMBDMIE = NWATERK
IF (LPC%RADIATIVE_PROPERTY_INDEX>0) THEN
   TA => TABLES(LPC%RADIATIVE_PROPERTY_INDEX)
   NL_DATA = TA%NUMBER_ROWS
ELSEIF (LPC%Y_INDEX==H2O_INDEX) THEN
   NL_DATA = NLMBDMIE
ELSEIF (LPC%FUEL) THEN
   NL_DATA = NFUELK
ELSE
   NL_DATA = 1
ENDIF

!     Allocate arrays

ALLOCATE(QSCA(1:NRDMIE,1:NLMBDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','QSCA',IZERO)
ALLOCATE(QABS(1:NRDMIE,1:NLMBDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','QABS',IZERO)
ALLOCATE(CHI_F(1:NRDMIE,1:NLMBDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','CHI_F',IZERO)
!
ALLOCATE(LMBDMIE(1:NLMBDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','LMBDMIE',IZERO)
ALLOCATE(LMBDWGHT(1:NLMBDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','LMBDWGHT',IZERO)
ALLOCATE(REAL_REF_INDX(1:NLMBDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','REAL_REF_INDX',IZERO)
ALLOCATE(CMPLX_REF_INDX(1:NLMBDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','CMPLX_REF_INDX',IZERO)

! Radiative properties

LMBDMIE(1:NLMBDMIE) = CPLXREF_WATER(1:NLMBDMIE,1)
IF (LPC%RADIATIVE_PROPERTY_INDEX > 0) THEN
   DO NX = 1,NLMBDMIE
      CALL INTERPOLATE1D(TA%TABLE_DATA(:,1)*1.0E-6_EB,TA%TABLE_DATA(:,2),LMBDMIE(NX),REAL_REF_INDX(NX))
      CALL INTERPOLATE1D(TA%TABLE_DATA(:,1)*1.0E-6_EB,TA%TABLE_DATA(:,3),LMBDMIE(NX),CMPLX_REF_INDX(NX))
   ENDDO
ELSEIF (LPC%Y_INDEX==H2O_INDEX) THEN
   REAL_REF_INDX              = 1.33_EB
   CMPLX_REF_INDX(1:NLMBDMIE) = CPLXREF_WATER(1:NLMBDMIE,2)
ELSEIF (LPC%FUEL) THEN
   DO NX = 1,NLMBDMIE
      CALL INTERPOLATE1D(CPLXREF_FUEL(:,1)*1.0E-6_EB,CPLXREF_FUEL(:,2),LMBDMIE(NX),REAL_REF_INDX(NX))
      CALL INTERPOLATE1D(CPLXREF_FUEL(:,1)*1.0E-6_EB,CPLXREF_FUEL(:,3),LMBDMIE(NX),CMPLX_REF_INDX(NX))
   ENDDO
ELSE
   REAL_REF_INDX              = LPC%REAL_REFRACTIVE_INDEX
   CMPLX_REF_INDX             = LPC%COMPLEX_REFRACTIVE_INDEX
ENDIF

CALL MIE_SCATTERING

!     Generate integration weights for lambda

IF (NLMBDMIE == 1) THEN
   LMBDWGHT(1) = 1._EB
ELSE
   LMBDWGHT(1) = 0.5_EB*(LMBDMIE(2)-LMBDMIE(1))
   DO I = 2, NLMBDMIE-1
      LMBDWGHT(I) = 0.5_EB*(LMBDMIE(I+1) - LMBDMIE(I-1))
   ENDDO
   LMBDWGHT(NLMBDMIE) = 0.5_EB*(LMBDMIE(NLMBDMIE)-LMBDMIE(NLMBDMIE-1))
ENDIF


!     Loop over all radiation bands

BANDLOOP: DO IBND = 1,NSB
   !
   IF (NSB == 1) THEN
      NLAMBDALOW = 1
      NLAMBDAHIGH = NLMBDMIE
   ELSE
      NLAMBDALOW  = MINLOC(LMBDMIE,MASK=LMBDMIE>=WL_LOW(IBND)*1.E-6)
      NLAMBDAHIGH = MAXLOC(LMBDMIE,MASK=LMBDMIE<=WL_HIGH(IBND)*1.0E-6)
   ENDIF

   !     Loop over all PARTICLE size groups


   DRGROUPLOOP: DO ND = 1, MIE_NDG

      LPC%R50(ND) = EXP(DGROUP_A*REAL(ND,EB) + DGROUP_B)

      !     Loop over wavelengths

      IBSUM = 0._EB

      DO J = NLAMBDALOW(1),NLAMBDAHIGH(1)
         IB = PLANCK(RADTMP, LMBDMIE(J)*1.0E6_EB)
         IBSUM = IBSUM + LMBDWGHT(J)*IB

         ASUM = 0._EB
         BSUM = 0._EB

         ! Properties at d32

         CALL INTERPOLATE1D(RDMIE,QSCA(:,J),LPC%R50(ND),AVAL)
         CALL INTERPOLATE1D(RDMIE,CHI_F(:,J),LPC%R50(ND),BVAL)
         BVAL = (1._EB-BVAL)
         ASUM = AVAL*BVAL
         CALL INTERPOLATE1D(RDMIE,QABS(:,J),LPC%R50(ND),BVAL)
         BSUM = BVAL

         LPC%WQSCA(ND,IBND) = LPC%WQSCA(ND,IBND) + ASUM*LMBDWGHT(J)*IB
         LPC%WQABS(ND,IBND) = LPC%WQABS(ND,IBND) + BSUM*LMBDWGHT(J)*IB
      ENDDO

      !     Normalize with blackbody radiation

      LPC%WQSCA(ND,IBND)  = LPC%WQSCA(ND,IBND)/IBSUM
      LPC%WQABS(ND,IBND)  = LPC%WQABS(ND,IBND)/IBSUM

ENDDO DRGROUPLOOP
ENDDO BANDLOOP

DEALLOCATE(RDMIE)
DEALLOCATE(QSCA)
DEALLOCATE(QABS)
DEALLOCATE(CHI_F)
DEALLOCATE(LMBDMIE)
DEALLOCATE(LMBDWGHT)
DEALLOCATE(REAL_REF_INDX)
DEALLOCATE(CMPLX_REF_INDX)

END SUBROUTINE MEAN_CROSS_SECTIONS



SUBROUTINE MIE_SCATTERING
!
!     Calculates the scattering and absorption cross sections
!     and calculates forward scattering fraction by integrating 
!     the scattering phase function.
!
! ----------------------------------------------------------------------
! -----------  SPECIFICATIONS FOR SUBROUTINE  MIEV0  ---------------
! ----------------------------------------------------------------------
INTEGER   MOMDIM
PARAMETER  ( MOMDIM = 200)
LOGICAL   ANYANG, PERFCT, PRNT( 2 )
INTEGER   IPOLZN, NMOM
REAL(EB)  GQSC, MIMCUT, PMOM( 0:MOMDIM, 4 ), SPIKE, QE, QS
REAL(EB), ALLOCATABLE :: XXX(:)
COMPLEX(EB) SFORW, SBACK, TFORW( 2 ), TBACK( 2 ), CREFIN
COMPLEX(EB), ALLOCATABLE :: S1(:), S2(:) 

! --------------- LOCAL VARIABLES --------------------------------------

INTEGER   I, J, K, NX, IZERO
INTEGER   NLAMBDA, NRA, NMIEANG2
REAL(EB)  STMP, AIJ, FTMP, XX_MAX
REAL(EB)  MUMIN1, MUMIN2,THETALIM1, THETALIM2, mudloc
REAL(EB)  mu1, mu2, nu1, nu2, mud0loc, mudPiloc, mud1, mud2, dmud 
REAL(EB),ALLOCATABLE :: XMU1(:), XNU1(:),XMU2(:), &
          ANGLE1(:), ANGLE2(:), MUD(:), MUDX(:),PWGHT(:), &
          PHSFUN(:), PFOR(:,:), MUD0(:,:), MUDPI(:,:)
! ----------------------------------------------------------------------

NRA    = NUMBER_RADIATION_ANGLES

!     MIEV-code variables

MIMCUT = 1.E-6_EB
PERFCT = .FALSE.
ANYANG = .TRUE.
!      IPOLZN = +1234
IPOLZN = 0
NMOM   = 0
PRNT   = .FALSE.

!     Limit for XX

XX_MAX = 15000.0_EB

!     Integration limits

THETALIM1 = ACOS(1._EB - 2._EB/REAL(NRA))
MUMIN1 = COS(THETALIM1)
MUMIN2 = MUMIN1**2-(1._EB-MUMIN1**2)
THETALIM2 = ACOS(MUMIN2)
NMIEANG2 = NMIEANG*2

!     Allocate local arrays

ALLOCATE(XXX(1:NRDMIE),STAT=IZERO)
CALL ChkMemErr('INIT','XXX',IZERO)
ALLOCATE(S1(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','S1',IZERO)
ALLOCATE(S2(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','S2',IZERO)
ALLOCATE(XMU1(1:NMIEANG),STAT=IZERO)
CALL ChkMemErr('INIT','XMU1',IZERO)
ALLOCATE(XNU1(1:NMIEANG),STAT=IZERO)
CALL ChkMemErr('INIT','XNU1',IZERO)
ALLOCATE(XMU2(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','XMU2',IZERO)
ALLOCATE(MUD(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','MUD',IZERO)
ALLOCATE(MUDX(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','MUDX',IZERO)
ALLOCATE(PWGHT(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','PWGHT',IZERO)
ALLOCATE(ANGLE1(1:NMIEANG),STAT=IZERO)
CALL ChkMemErr('INIT','ANGLE1',IZERO)
ALLOCATE(ANGLE2(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','ANGLE2',IZERO)
ALLOCATE(PHSFUN(1:NMIEANG2),STAT=IZERO)
CALL ChkMemErr('INIT','PHSFUN',IZERO)
ALLOCATE(PFOR(1:NMIEANG,1:NMIEANG),STAT=IZERO)
CALL ChkMemErr('INIT','PFOR',IZERO)
ALLOCATE(MUD0(1:NMIEANG,1:NMIEANG),STAT=IZERO)
CALL ChkMemErr('INIT','MUD0',IZERO)
ALLOCATE(MUDPI(1:NMIEANG,1:NMIEANG),STAT=IZERO)
CALL ChkMemErr('INIT','MUDPI',IZERO)

!     Create solid angle integration arrays

DO I = 1,NMIEANG
   ANGLE1(I) = THETALIM1*REAL(I-1)/REAL(NMIEANG-1)
   ANGLE1(I) = ANGLE1(I)*(1._EB-0.99*REAL(NMIEANG-I)/REAL(NMIEANG))
ENDDO
DO I = 1,NMIEANG
   XMU1(NMIEANG-I+1) = COS(ANGLE1(I))
   XNU1(NMIEANG-I+1) = SIN(ANGLE1(I))
ENDDO

!    Create phase function ingtegration arrays

DO I = 1,NMIEANG2
   ANGLE2(I) = THETALIM2*REAL(I-1)/REAL(NMIEANG2-1)
   ANGLE2(I) = ANGLE2(I)*(1._EB-0.99*REAL(NMIEANG2-I)/REAL(NMIEANG2))
ENDDO
DO I = 1,NMIEANG2
   XMU2(NMIEANG2-I+1) = COS(ANGLE2(I))
ENDDO

!    Calculate phase function ingetration limits

DO J = 1,NMIEANG
   DO I = 1,NMIEANG
      MUD0(I,J)  = XMU1(I)*XMU1(J) + XNU1(I)*XNU1(J)
      MUDPI(I,J) = XMU1(I)*XMU1(J) - XNU1(I)*XNU1(J)
   ENDDO
ENDDO

!     Calculate phase function integration weights

mu1 = 0.7_EB
mu2 = 0.9_EB
nu1 = SQRT(1-mu1**2)
nu2 = SQRT(1-mu2**2)
mud0loc  = mu1*mu2 + nu1*nu2
mudPiloc = mu1*mu2 - nu1*nu2
mud1 = mudPiloc
mud2 = mud0loc
dmud = (mud2-mud1)/(NMIEANG2-1)
DO I = 1,NMIEANG2
   mud(I) = mud1+REAL(I-1)*dmud
ENDDO
mud(1)       = mud(1)       + 0.25_EB*dmud !empirical
mud(NMIEANG2) = mud(NMIEANG2) - 0.25_EB*dmud !empirical
DO I = 1,NMIEANG2
   MUDX(I) = (mud(I)-mud1)/(mud2-mud1)
ENDDO
DO I = 1,NMIEANG2
   PWGHT(I) = dmud/sqrt((nu1*nu2)**2-(mud(I)-mu1*mu2)**2)
ENDDO
PWGHT(2)         = 0.5_EB*PWGHT(2)
PWGHT(NMIEANG2-1) = 0.5_EB*PWGHT(NMIEANG2-1)

!     Loop over wavelength

LAMBDALOOP: DO NLAMBDA = 1, NLMBDMIE
!
   CREFIN = CMPLX( REAL_REF_INDX(NLAMBDA), CMPLX_REF_INDX(NLAMBDA) )

! Choose Perfectly reflecting sphere, if large real index is given.

   IF (REAL_REF_INDX(NLAMBDA) > 10._EB) PERFCT = .TRUE.

!     Loop over PARTICLE radius

   RADIUSLOOP: DO NX = 1, NRDMIE

      XXX(NX) = MIN(XX_MAX,2._EB*PI*RDMIE(NX)/LMBDMIE(NLAMBDA))
      CALL MIEV0( XXX(NX), CREFIN, PERFCT, MIMCUT, ANYANG,   &
                  NMIEANG2, XMU2, NMOM, IPOLZN, MOMDIM, PRNT,  &
                  QE, QS, GQSC, &
                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK, SPIKE )

      QSCA(NX,NLAMBDA) = QS
      QABS(NX,NLAMBDA) = QE-QS

!     Calculate single drop phase function

      IF (ABS(QS)>TWO_EPSILON_EB) THEN
         DO I = 1,NMIEANG2
            PHSFUN(I) = 2._EB*(abs(S1(I))**2 + abs(S2(I))**2 )
         ENDDO
         PHSFUN = PHSFUN/(QS*XXX(NX)**2)
      ELSE
         PHSFUN = 1.0_EB
      ENDIF

!     Calculate the innermost integral of the forward scattering fraction

      PFOR = 0._EB
      DO J = 1,NMIEANG
         DO I = J,NMIEANG
            IF (ABS(MUD0(I,J)-MUDPI(I,J))<=TWO_EPSILON_EB) THEN
               CALL INTERPOLATE1D(XMU2,PHSFUN,MUD0(I,J),FTMP)
               PFOR(I,J) = PI*FTMP
            ELSE
               mud1 = MUDPI(I,J) 
               mud2 = MUD0(I,J)
               STMP = 0.0_EB
               DO K = 1,NMIEANG2
                  mudloc = mud1+MUDX(K)*(mud2-mud1)
                  CALL INTERPOLATE1D(XMU2,PHSFUN,mudloc,FTMP)
                  STMP = STMP + PWGHT(K)*FTMP
               ENDDO
               PFOR(I,J) = STMP 
            ENDIF
         ENDDO
      ENDDO

!     Calculate the two outer integrals of the forward fraction

      STMP = 0._EB
      DO J = 1,NMIEANG-1
         DO I = J+1,NMIEANG-1
            AIJ = (XMU1(I+1)-XMU1(I))*(XMU1(J+1)-XMU1(J))/2._EB
            STMP = STMP +   (2._EB*PFOR(I,J)+PFOR(I+1,J)+PFOR(I,J+1)+2._EB*PFOR(I+1,J+1))*AIJ/3._EB
         ENDDO
      ENDDO
      DO I = 1,NMIEANG-1
         AIJ = ((XMU1(I+1)-XMU1(I))**2)/2._EB
         STMP = STMP + (PFOR(I,I)+PFOR(I+1,I)+PFOR(I+1,I+1))*AIJ/3._EB
      ENDDO
      CHI_F(NX,NLAMBDA) = 2._EB*STMP/(4._EB*PI/NRA)      
   ENDDO RADIUSLOOP
ENDDO LAMBDALOOP

DEALLOCATE(XXX)
DEALLOCATE(S1)
DEALLOCATE(S2)
DEALLOCATE(XMU1)
DEALLOCATE(XNU1)
DEALLOCATE(XMU2)
DEALLOCATE(MUD)
DEALLOCATE(MUDX)
DEALLOCATE(PWGHT)
DEALLOCATE(ANGLE1)
DEALLOCATE(ANGLE2)
DEALLOCATE(PHSFUN)
DEALLOCATE(PFOR)
DEALLOCATE(MUD0)
DEALLOCATE(MUDPI)

END SUBROUTINE MIE_SCATTERING

SUBROUTINE MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU, &
                  NMOM, IPOLZN, MOMDIM, PRNT, QEXT, QSCA, GQSC, &
                  PMOM, SFORW, SBACK, S1, S2, TFORW, TBACK, &
                  SPIKE )

!     Mie scattering for a single PARTICLE and wavelength.
!     Author:  Dr. Warren J. Wiscombe (wiscombe@climate.gsfc.nasa.gov)
!         NASA Goddard Space Flight Center
!         Code 913
!         Greenbelt, MD 20771

!     REFERENCES
!     ----------
!
!     (1) Wiscombe, W., 1979: Mie Scattering Calculations--Advances
!         in Technique And Fast, Vector-Speed Computer Codes,
!         Ncar Tech Note TN-140+STR, National Center For
!         Atmospheric Research, Boulder, Colorado (out of print
!         but an updated electronic version available)
!
!     (2) Wiscombe, W., 1980: Improved Mie scattering algorithms,
!         Appl. Opt. 19, 1505-1509
!

!    Computes Mie scattering and extinction efficiencies; asymmetry
!    factor;  forward- and backscatter amplitude;  scattering
!    amplitudes vs. scattering angle for incident polarization parallel
!    and perpendicular to the plane of scattering;
!    coefficients in the Legendre polynomial expansions of either the
!    unpolarized phase function or the polarized phase matrix;
!    some quantities needed in polarized radiative transfer;  and
!    information about whether or not a resonance has been hit.
!
!    Input and output variables are described in file MIEV.doc. 
!    Many statements are accompanied by comments referring to 
!    references in MIEV.doc, notably the NCAR Mie report which is now
!    available electronically and which is referred to using the
!    shorthand (Rn), meaning Eq. (n) of the report.

!    CALLING TREE:
!
!        MIEV0
!            TESTMI
!                TSTBAD
!                MIPRNT
!                ERRMSG
!            CKINMI
!                WRTBAD
!                WRTDIM
!                ERRMSG
!            SMALL1
!            SMALL2
!            ERRMSG
!            BIGA
!                CONFRA
!                    ERRMSG
!            LPCOEF
!                LPCO1T
!                LPCO2T
!                ERRMSG
!            MIPRNT
!
!   I N P U T   V A R I A B L E S
!   -----------------------------
!
!  ( Even if an input variable is not needed for a particular
!    application, make sure it has a legitimate value that can
!    be written out and read in -- no indefinites, etc. )
!
!  XX        Mie size parameter ( 2 * pi * radius / wavelength )
!
!  CREFIN    Complex refractive index ( imag part can be + or -,
!            but internally a negative imaginary index is assumed ).
!            If imag part is - ,  scattering amplitudes as in Van
!            de Hulst are returned;  if imag part is + , complex
!            conjugates of those scattering amplitudes are returned
!            (the latter is the convention in physics).
!            ** NOTE ** In the 'PERFECT' case, scattering amplitudes
!            in the Van de Hulst (Ref. 6 above) convention will
!            automatically be returned unless  Im(CREFIN)  is
!            positive;  otherwise, CREFIN plays no role.
!
!  PERFCT    TRUE, assume refractive index is infinite and use
!            special case formulas for Mie coefficients  'a'
!            and  'b'  ( see Kerker, M., The Scattering of
!            Light and Other Electromagnetic Radiation, p. 90 ).
!            This is sometimes called the 'totally reflecting',
!            sometimes the 'perfectly conducting' case.
!            ( see CREFIN for additional information )
!
!  MIMCUT    (positive) value below which imaginary refractive
!            index is regarded as zero (computation proceeds
!            faster for zero imaginary index)
!
!  ANYANG    TRUE, any angles whatsoever may be input through
!            XMU.  FALSE, the angles are monotone increasing
!            and mirror symmetric about 90 degrees (this option
!            is advantageous because the scattering amplitudes
!            S1,S2 for the angles between 90 and 180 degrees
!            are evaluable from symmetry relations, and hence
!            are obtained with little added computational cost.)
!
!  NUMANG    No. of angles at which scattering amplitudes
!            S1,S2 are to be evaluated  ( set = 0 to skip
!            calculation of S1,S2 ).  Make sure NUMANG does
!            not exceed the parameter MAXANG in the program.
!
!  XMU(N)    Cosines of angles ( N = 1 TO NUMANG ) at which S1,S2
!            are to be evaluated.  If ANYANG = FALSE, then
!
!             (a) the angles must be monotone increasing and
!                 mirror symmetric about 90 degrees (if 90-A is
!                 an angle, then 90+A must be also)
!
!             (b) if NUMANG is odd, 90 degrees must be among
!                 the angles
!
!  NMOM       Highest Legendre moment PMOM to calculate,
!             numbering from zero ( NMOM = 0 prevents
!             calculation of PMOM )
!
!  IPOLZN     POSITIVE, Compute Legendre moments PMOM for the
!                       Mueller matrix elements determined by the
!                       digits of IPOLZN, with 1 referring to M1,
!                       2 to M2, 3 to S21, and 4 to D21 (Ref. 3).
!                       E.g., if IPOLZN = 14 then only moments for
!                       M1 and D21 will be returned.
!
!             0,        Compute Legendre moments PMOM for the
!                       unpolarized unnormalized phase function.
!
!             NEGATIVE, Compute Legendre moments PMOM for the
!                       Sekera phase quantities determined by the
!                       digits of ABS(IPOLZN), with 1 referring to
!                       R1, 2 to R2, 3 to R3, and 4 to R4 (REF. 4).
!                       E.g., if IPOLZN = -14 then only moments for
!                       R1 and R4 will be returned.
!
!             ( NOT USED IF  NMOM = 0 )
!
!  MOMDIM     Determines first dimension of PMOM, which is dimensioned
!             internally as PMOM( 0:MOMDIM, * ) (second dimension must
!             be the larger of unity and the highest digit in
!             IPOLZN; if not, serious errors will occur).
!             Must be given a value, even if  NMOM = 0.  Minimum: 1.
!
!  PRT(L)     Print flags (LOGICAL).  L = 1  prints  S1,S2, their
!             squared absolute values, and degree of polarization,
!             provided NUMANG is non-zero.   L = 2  prints all
!             output variables other than  S1,S2.
!
!
! O U T P U T   V A R I A B L E S
! -------------------------------
!
!  QEXT      (REAL) extinction efficiency factor  ( Ref. 2, Eq. 1A )
!
!  QSCA      (REAL) scattering efficiency factor  ( Ref. 2, Eq. 1B )
!
!  GQSC      (REAL) asymmetry factor times scattering efficiency
!            ( Ref. 2, Eq. 1C )  ( allows calculation of radiation
!            pressure efficiency factor  QPR = QEXT - GQSC )
!
!  =====================================================================
!  ==== NOTE --  S1, S2, SFORW, SBACK, TFORW, AND TBACK are calculated
!  ====          internally for negative imaginary refractive index;
!  ====          for positive imaginary index, their complex conjugates
!  ====          are taken before they are returned, to correspond to
!  ====          customary usage in some parts of physics ( in parti-
!  ====          cular, in papers on CAM approximations to Mie theory ).
!  =====================================================================
!
!  S1(N),    (COMPLEX) Mie scattering amplitudes at angles specified
!  S2(N)     by XMU(N) ( N=1 to NUMANG )  ( Ref. 2, Eqs. 1d-e ).
!
!  SFORW     (COMPLEX) forward-scattering amplitude S1 at
!            0 degrees.  ( S2(0 deg) = S1(0 deg) )
!
!  SBACK     (COMPLEX) backscattering amplitude S1 at
!            180 degrees.   ( S2(180 deg) = - S1(180 deg) )
!
!  TFORW(I)  (COMPLEX) values of
!
!                I=1:  T1 = ( S2 - (MU)*S1 ) / ( 1 - MU**2 )
!                I=2:  T2 = ( S1 - (MU)*S2 ) / ( 1 - MU**2 )
!
!            At angle theta = 0 ( MU = COS(theta) = 1 ), where the
!            expressions on the right-hand side are indeterminate.
!            ( these quantities are required for doing polarized
!            radiative transfer (Ref. 4, Appendix). )
!  TBACK(I)  (COMPLEX) values of  T1 (for I=1) or  T2 (for I=2) at
!            angle  theta = 180 degrees ( MU = COS(theta) = - 1 ).
!
!  SPIKE     (REAL) magnitude of the smallest denominator of
!            either Mie coefficient (a-sub-n or b-sub-n),
!            taken over all terms in the Mie series past
!            N = size parameter XX.  Values of SPIKE below
!            about 0.3 signify a ripple spike, since these
!            spikes are produced by abnormally small denominators
!            in the Mie coefficients (normal denominators are of
!            order unity or higher).  Defaults to 1.0 when not
!            on a spike.  Does not identify all resonances
!            (we are still working on that).
!
! PMOM(M,NP) (REAL) moments  M = 0 to NMOM  of unnormalized NP-th
!            phase quantity PQ  ( moments with  M > 2*NTRM  are
!            zero, where  NTRM = no. terms in Mie series =
!            XX + 4*XX**1/3 + 1 ) :
!
!              PQ( MU, NP ) = sum( M=0 to infinity ) ( (2M+1)
!                                * PMOM( M,NP ) * P-sub-M( MU ) )
!
!            WHERE  MU = COS( scattering angle )
!                   P-sub-M = M-th Legendre polynomial
!
!            and the definition of 'PQ' is as follows:
!
!            IPOLZN>0:  PQ(MU,1) = CABS( S1(MU) )**2
!                          PQ(MU,2) = CABS( S2(MU) )**2
!                          PQ(MU,3) = RE( S1(MU)*CONJG( S2(MU) ) )
!                          PQ(MU,4) = - IM( S1(MU)*CONJG( S2(MU) ) )
!                          ( called M1, M2, S21, D21 in literature )
!
!            IPOLZN=0:  PQ(MU,1) = ( CABS(S1)**2 + CABS(S2)**2 ) / 2
!                       ( the unnormalized phase function )
!
!            IPOLZN<0:  PQ(MU,1) = CABS( T1(MU) )**2
!                          PQ(MU,2) = CABS( T2(MU) )**2
!                          PQ(MU,3) = RE( T1(MU)*CONJG( T2(MU) ) )
!                          PQ(MU,4) = - IM( T1(MU)*CONJG( T2(MU) ) )
!                          ( called R1, R2, R3, R4 in literature )
!
!            The sign of the 4th phase quantity is a source of
!            confusion.  It flips if the complex conjugates of
!            S1,S2  or  T1,T2  are used, as occurs when a
!            refractive index with positive imaginary part is
!            used (see discussion below).  The definition above
!            is consistent with a negative imaginary part.
!
!            ** WARNING **  Make sure the second dimension of PMOM
!            in the calling program is at least as large as the
!            absolute value of IPOLZN.
!
!            For small enough values of XX, or large enough values
!            of M,  PMOM  will tend to underflow.  Thus, it is
!            unwise to assume the values returned are non-zero and,
!            for example, to divide some quantity by them.
!


!      I N T E R N A L   V A R I A B L E S
!      -----------------------------------

!  AN,BN           Mie coefficients a-sub-n, b-sub-n ( Ref. 1, Eq. 16 )
!  ANM1,BNM1       Mie coefficients  a-sub-(n-1),
!                     b-sub-(n-1);  used in GQSC sum
!  ANP             Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
!  BNP             Coeffs. in S- expansion ( Ref. 2, p. 1507 )
!  ANPM            Coeffs. in S+ expansion ( Ref. 2, p. 1507 )
!                     when  MU  is replaced by  - MU
!  BNPM            Coeffs. in S- expansion ( Ref. 2, p. 1507 )
!                     when  MU  is replaced by  - MU
!  CALCMO(K)       TRUE, calculate moments for K-th phase quantity
!                     (derived from IPOLZN)
!  CBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
!                     ( COMPLEX version )
!  CDENAN,         (COMPLEX) denominators of An,Bn
!   CDENBN
!  CIOR            Complex index of refraction with negative
!                     imaginary part (Van de Hulst convention)
!  CIORIV          1 / cIoR
!  COEFF           ( 2N + 1 ) / ( N ( N + 1 ) )
!  CSUM1,2         temporary sum variables for TFORW, TBACK
!  FN              Floating point version of loop index for
!                     Mie series summation
!  LITA,LITB(N)    Mie coefficients An, Bn, saved in arrays for
!                     use in calculating Legendre moments PMOM
!  MAXTRM          Max. possible no. of terms in Mie series
!  MM              (-1)^(n+1), where n is Mie series sum index 
!  MIM             Magnitude of imaginary refractive index
!  MRE             Real part of refractive index
!  MAXANG          Max. possible value of input variable NUMANG
!  NANGD2          (NUMANG+1)/2 ( no. of angles in 0-90 deg; ANYANG=F )
!  NOABS           TRUE, sphere non-absorbing (determined by MIMCUT)
!  NP1DN           ( N + 1 ) / N
!  NPQUAN          Highest-numbered phase quantity for which moments are
!                     to be calculated (the largest digit in IPOLZN
!                     if  IPOLZN /= 0)
!  NTRM            No. of terms in Mie series
!  PASS1           TRUE on first entry, FALSE thereafter; for self-test
!  PIN(J)          Angular function pi-sub-n ( Ref. 2, Eq. 3 )
!                     at J-th angle
!  PINM1(J)        pi-sub-(n-1) ( see PIn ) at J-th angle
!  PSINM1          Ricatti-Bessel function psi-sub-(n-1), argument XX
!  PSIN            Ricatti-Bessel function psi-sub-n of argument XX
!                     ( Ref. 1, p. 11 ff. )
!  RBIGA(N)        Bessel function ratio A-sub-N (Ref. 2, Eq. 2)
!                     ( REAL version, for when imag refrac index = 0 )
!  RIORIV          1 / Mre
!  RN              1 / N
!  RTMP            (REAL) temporary variable
!  SP(J)           S+  for J-th angle  ( Ref. 2, p. 1507 )
!  SM(J)           S-  for J-TH angle  ( Ref. 2, p. 1507 )
!  SPS(J)          S+  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )
!  SMS(J)          S-  for (NUMANG+1-J)-th angle ( ANYANG=FALSE )
!  TAUN            Angular function tau-sub-n ( Ref. 2, Eq. 4 )
!                     at J-th angle
!  TCOEF           N ( N+1 ) ( 2N+1 ) (for summing TFORW,TBACK series)
!  TWONP1          2N + 1
!  YESANG          TRUE if scattering amplitudes are to be calculated
!  ZETNM1          Ricatti-Bessel function  zeta-sub-(n-1) of argument
!                     XX  ( Ref. 2, Eq. 17 )
!  ZETN            Ricatti-Bessel function  zeta-sub-n of argument XX
! ----------------------------------------------------------------------
!
IMPLICIT  NONE
!
! ----------------------------------------------------------------------
! --------  I / O SPECIFICATIONS FOR SUBROUTINE MIEV0  -----------------
! ----------------------------------------------------------------------
LOGICAL     ANYANG, PERFCT, PRNT(*)
INTEGER     IPOLZN, MOMDIM, NUMANG, NMOM
REAL(EB)    GQSC, MIMCUT, PMOM( 0:MOMDIM, * ), QEXT, QSCA, SPIKE,XMU(*), XX
COMPLEX(EB) CREFIN, SFORW, SBACK, S1(*), S2(*), TFORW(*), TBACK(*)
! ----------------------------------------------------------------------
!
!                                  ** NOTE --  MAXTRM = 10100  is neces-
!                                  ** sary to do some of the test probs,
!                                  ** but 1100 is sufficient for most
!                                  ** conceivable applications
!     .. Parameters ..

INTEGER     MAXANG
PARAMETER   ( MAXANG = 5000 )
INTEGER     MAXTRM
!      PARAMETER ( MAXTRM = 10100 )
PARAMETER   ( MAXTRM = 16000 ) ! works for FDS
REAL(EB)    ONETHR
PARAMETER   ( ONETHR = 1._EB / 3._EB )

!     .. Local Scalars ..

LOGICAL     NOABS, PASS1, YESANG, GT100
INTEGER     I, J, N, NANGD2, NPQUAN, NTRM
REAL(EB)    CHIN, CHINM1, COEFF, DENAN, DENBN, FN, MIM, MM, MRE, &
            NP1DN, PSIN, PSINM1, RATIO, RIORIV, RN, RTMP, TAUN, &
            TCOEF, TWONP1, XINV
COMPLEX(EB) AN, ANM1, ANP, ANPM, BN, BNM1, BNP, BNPM, CDENAN,&
            CDENBN, CIOR, CIORIV, CSUM1, CSUM2, ZET, ZETN, ZETNM1

!     .. Local Arrays ..

LOGICAL     CALCMO( 4 )
REAL(EB), ALLOCATABLE :: PIN(:), PINM1(:)
REAL(EB) RBIGA( MAXTRM )
COMPLEX(EB) CBIGA( MAXTRM ), LITA( MAXTRM ), LITB( MAXTRM )
COMPLEX(EB), ALLOCATABLE :: SM(:), SMS(:), SP(:), SPS(:)

SAVE      PASS1
DATA      PASS1 / .TRUE. /
!
ALLOCATE(PIN(1:NUMANG)) 
ALLOCATE(PINM1(1:NUMANG))
ALLOCATE(SM(1:NUMANG)) 
ALLOCATE(SMS(1:(NUMANG+1)/2)) 
ALLOCATE(SP(1:NUMANG)) 
ALLOCATE(SPS(1:(NUMANG+1)/2)) 

!                    ** Save some input variables and replace them
!                    ** with values needed to do the self-test

IF( PASS1 ) CALL TESTMI( .FALSE., XX, CREFIN, MIMCUT, PERFCT, &
                         ANYANG, NMOM, IPOLZN, NUMANG, XMU, QEXT, &
                         QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW, &
                         TBACK, PMOM, MOMDIM )
MAIN_MIEV: DO 
!                                        ** Check input and calculate
!                                        ** certain variables from input
CALL CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM, NMOM, &
             IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )

IF( PERFCT .AND. XX<=0.1_EB ) THEN
!                                            ** Use totally-reflecting
!                                            ** small-particle limit
   CALL SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK, &
                  S1, S2, TFORW, TBACK, LITA, LITB )
   NTRM = 2
ELSE
   NOABS = .TRUE.
   GT100 = .FALSE.

   IF( .NOT.PERFCT ) THEN
      CIOR = CREFIN
      IF( AIMAG(CIOR)>0.0_EB ) CIOR = CONJG( CIOR )
      MRE    = REAL( CIOR )
      MIM    = -AIMAG( CIOR )
      NOABS  = MIM<=MIMCUT
      CIORIV = 1.0_EB / CIOR
      RIORIV = 1.0_EB / MRE

      IF( XX*MAX( 1._EB, ABS(CIOR) )<=0.1_EB ) THEN
!                                    ** Use general-refractive-index
!                                    ** small-particle limit
         CALL SMALL2( XX, CIOR, MIM>MIMCUT, NUMANG, XMU, QEXT, &
                           QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW, &
                           TBACK, LITA, LITB )
         NTRM = 2
         GT100 = .TRUE.
      END IF
   END IF
   GT100IF: IF (.NOT. GT100) THEN

      NANGD2 = ( NUMANG + 1 ) / 2
      YESANG = NUMANG>0
   !                             ** Number of terms in Mie series; Eq R50
      IF( XX<=8.0_EB ) THEN
      NTRM = XX + 4._EB*XX**ONETHR + 1._EB
      ELSE IF( XX<4200._EB ) THEN
      NTRM = XX + 4.05*XX**ONETHR + 2._EB
      ELSE
      NTRM = XX + 4._EB*XX**ONETHR + 2._EB
      END IF
      IF( NTRM+1 > MAXTRM )CALL ERRMSG('MIEV0--PARAMETER MaxTrm TOO SMALL',.TRUE.)
   !                            ** Calculate logarithmic derivatives of
   !                            ** J-Bessel-fcn., A-sub-(1 to NTrm)
      IF( .NOT.PERFCT ) CALL BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )
   !                            ** Initialize Ricatti-Bessel functions
   !                            ** (psi,chi,zeta)-sub-(0,1) for upward
   !                            ** recurrence ( Eq. R19 )
      XINV   = 1.0_EB / XX
      PSINM1 = SIN( XX )
      CHINM1 = COS( XX )
      PSIN   = PSINM1*XINV - CHINM1
      CHIN   = CHINM1*XINV + PSINM1
      ZETNM1 = CMPLX( PSINM1, CHINM1 )
      ZETN   = CMPLX( PSIN, CHIN )
   !                                     ** Initialize previous coeffi-
   !                                     ** cients for GQSC series
      ANM1 = ( 0.0_EB, 0.0_EB )
      BNM1 = ( 0.0_EB, 0.0_EB )
   !                             ** Initialize angular function  pi
   !                             ** and sums for S+, S- ( Ref. 2, p. 1507 )
      IF( ANYANG ) THEN
      DO J = 1, NUMANG
   !                             ** Eq. R39
            PINM1( J ) = 0.0_EB
            PIN( J ) = 1.0_EB

               SP( J ) = ( 0.0_EB, 0.0_EB )
               SM( J ) = ( 0.0_EB, 0.0_EB )
         END DO
      ELSE
         DO J = 1, NANGD2
   !                          ** Eq. R39
            PINM1( J ) = 0.0_EB
            PIN( J ) = 1.0_EB
            SP( J ) = ( 0.0_EB, 0.0_EB )
            SM( J ) = ( 0.0_EB, 0.0_EB )
            SPS( J ) = ( 0.0_EB, 0.0_EB )
            SMS( J ) = ( 0.0_EB, 0.0_EB )
         END DO
      END IF
   !                       ** Initialize Mie sums for efficiencies, etc.
      QSCA  = 0.0_EB
      GQSC  = 0.0_EB
      SFORW = ( 0._EB, 0._EB )
      SBACK = ( 0._EB, 0._EB )
      CSUM1 = ( 0._EB, 0._EB )
      CSUM2 = ( 0._EB, 0._EB )
   !
   ! ---------  LOOP TO SUM MIE SERIES  -----------------------------------
      MM     = +1.0_EB
      SPIKE  = 1.0_EB
      DO N = 1, NTRM
   !                           ** Compute various numerical coefficients
         FN     = N
         RN     = 1.0_EB / FN
         NP1DN  = 1.0_EB + RN
         TWONP1 = 2*N + 1
         COEFF  = TWONP1 / ( FN * ( N + 1 ) )
         TCOEF  = TWONP1 * ( FN * ( N + 1 ) )
   !                          ** Calculate Mie series coefficients
         IF( PERFCT ) THEN
   !                                 ** Totally-reflecting case; Eq R/A.1,2
            AN = ( ( FN*XINV )*PSIN - PSINM1 ) / ( ( FN*XINV )*ZETN - ZETNM1 )
            BN = PSIN / ZETN
         ELSE IF( NOABS ) THEN
   !                                      ** No-absorption case; Eq (R16)
            CDENAN = ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            AN   = ( ( RIORIV*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 ) / CDENAN
            CDENBN = ( MRE*RBIGA(N) + ( FN*XINV ) ) * ZETN - ZETNM1
            BN   = ( ( MRE*RBIGA(N) + ( FN*XINV ) ) * PSIN - PSINM1 ) / CDENBN
         ELSE
   !                                       ** Absorptive case; Eq (R16)
            CDENAN = ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            CDENBN =   ( CIOR*CBIGA( N ) + ( FN*XINV ) )*ZETN - ZETNM1
            AN   = ( ( CIORIV*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 ) / CDENAN
            BN     = ( ( CIOR*CBIGA( N ) + ( FN*XINV ) )*PSIN - PSINM1 ) / CDENBN
   !                                         ** Eq (R7)
            QSCA   = QSCA + TWONP1*( SQ( AN ) + SQ( BN ) )
         END IF
   !                       ** Save Mie coefficients for PMOM calculation
         LITA( N ) = AN
         LITB( N ) = BN
   !
         IF( .NOT.PERFCT .AND. N>XX ) THEN
   !                                               ** Flag resonance spikes
            DENAN  = ABS( CDENAN )
            DENBN  = ABS( CDENBN )
   !                                                   ** Eq. R/B.9
            RATIO  = DENAN / DENBN
   !                                                   ** Eq. R/B.10
            IF( RATIO<=0.2_EB .OR. RATIO>=5.0_EB ) SPIKE = MIN( SPIKE, DENAN, DENBN )
         END IF
   !                                  ** Increment Mie sums for non-angle-
   !                                  ** dependent quantities
   !                                                   ** Eq. R/B.2
         SFORW = SFORW + TWONP1*( AN + BN )
   !                                                   ** Eq. R/B.5,6
         CSUM1 = CSUM1 + TCOEF *( AN - BN )
   !                                                   ** Eq. R/B.1
         SBACK = SBACK + ( MM*TWONP1 )*( AN - BN )
   !                                                   ** Eq. R/B.7,8
         CSUM2 = CSUM2 + ( MM*TCOEF ) *( AN + BN )
   !                                         ** Eq (R8)
         GQSC  = GQSC  + (FN - RN) * REAL(ANM1 * CONJG(AN) + BNM1 * CONJG(BN) ) + COEFF * REAL(AN * CONJG(BN))
         IF( YESANG ) THEN
   !                                      ** Put Mie coefficients in form
   !                                      ** needed for computing S+, S-
   !                                      ** ( Eq R10 )
         ANP = COEFF*( AN + BN )
         BNP = COEFF*( AN - BN )
   !                                      ** Increment Mie sums for S+, S-
   !                                      ** while upward recursing
   !                                      ** angular functions pi and tau
         IF( ANYANG ) THEN
   !                                         ** Arbitrary angles
   !                                              ** vectorizable loop
            DO J = 1, NUMANG
   !                                                 ** Eq. (R37b)
               RTMP = ( XMU(J) * PIN(J) ) - PINM1( J )
   !                                                 ** Eq. (R38b)
               TAUN   = FN * RTMP - PINM1( J )
   !                                                   ** Eq (R10)
               SP( J ) = SP( J ) + ANP * ( PIN( J ) + TAUN )
               SM( J ) = SM( J ) + BNP * ( PIN( J ) - TAUN )
               PINM1( J ) = PIN( J )
   !                                                 ** Eq. R37c
               PIN( J ) = ( XMU( J ) * PIN( J ) ) + NP1DN * RTMP
            END DO
         ELSE
   !                                  ** Angles symmetric about 90 degrees
            ANPM = MM*ANP
            BNPM = MM*BNP
   !                                          ** vectorizable loop
            DO J = 1, NANGD2
   !                                                 ** Eq. (R37b)
               RTMP = ( XMU(J) * PIN(J) ) - PINM1( J )
   !                                                 ** Eq. (R38b)
               TAUN = FN * RTMP - PINM1( J )
   !                                                 ** Eq (R10,12)
               SP ( J ) = SP ( J ) + ANP * ( PIN( J ) + TAUN )
               SMS( J ) = SMS( J ) + BNPM *( PIN( J ) + TAUN )
               SM ( J ) = SM ( J ) + BNP * ( PIN( J ) - TAUN )
               SPS( J ) = SPS( J ) + ANPM *( PIN( J ) - TAUN )
               PINM1( J ) = PIN( J )
         !                                                 ** Eq. R37c
               PIN( J ) = ( XMU(J) * PIN(J) ) + NP1DN * RTMP
            END DO
         END IF
         END IF
   !                          ** Update relevant quantities for next
   !                          ** pass through loop
         MM   = - MM
         ANM1 = AN
         BNM1 = BN
   !                           ** Upward recurrence for Ricatti-Bessel
   !                           ** functions ( Eq. R17 )
         ZET    = ( TWONP1*XINV ) * ZETN - ZETNM1
         ZETNM1 = ZETN
         ZETN   = ZET
         PSINM1 = PSIN
         PSIN   = REAL( ZETN )
      END DO
   ! ---------- END LOOP TO SUM MIE SERIES --------------------------------
   !
   !                                         ** Eq (R6)
      QEXT = 2._EB / XX**2*REAL( SFORW )
      IF( PERFCT .OR. NOABS ) THEN
         QSCA = QEXT
      ELSE
         QSCA = 2._EB/ XX**2 * QSCA
      END IF
      GQSC   = 4._EB/ XX**2 * GQSC
      SFORW  = 0.5_EB*SFORW
      SBACK  = 0.5_EB*SBACK
      TFORW( 1 ) =  0.5_EB*SFORW - 0.125_EB*CSUM1
      TFORW( 2 ) =  0.5_EB*SFORW + 0.125_EB*CSUM1
      TBACK( 1 ) = -0.5_EB*SBACK + 0.125_EB*CSUM2
      TBACK( 2 ) =  0.5_EB*SBACK + 0.125_EB*CSUM2
      IF( YESANG ) THEN
   !                                ** Recover scattering amplitudes
   !                                ** from S+, S- ( Eq (R11) )
         IF( ANYANG ) THEN
   !                                         ** vectorizable loop
            DO J = 1, NUMANG
   !                                                  ** Eq (R11)
               S1( J ) = 0.5_EB*( SP( J ) + SM( J ) )
               S2( J ) = 0.5_EB*( SP( J ) - SM( J ) )
            END DO
         ELSE
   !                                        ** vectorizable loop
            DO J = 1, NANGD2
   !                                                   ** Eq (R11)
               S1( J ) = 0.5_EB*( SP( J ) + SM( J ) )
               S2( J ) = 0.5_EB*( SP( J ) - SM( J ) )
            END DO
   !                                        ** vectorizable loop
            DO  J = 1, NANGD2
               S1( NUMANG + 1 - J ) = 0.5_EB*( SPS( J ) + SMS( J ) )
               S2( NUMANG + 1 - J ) = 0.5_EB*( SPS( J ) - SMS( J ) )
            END DO
         END IF
      END IF
   END IF GT100IF
   !                                 ** Calculate Legendre moments
END IF
IF( NMOM>0 ) CALL LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO,NPQUAN, LITA, LITB, PMOM )
IF( AIMAG( CREFIN )>0.0_EB ) THEN
!                                         ** Take complex conjugates
!                                         ** of scattering amplitudes
   SFORW = CONJG( SFORW )
   SBACK = CONJG( SBACK )
   DO I = 1, 2
      TFORW( I ) = CONJG( TFORW( I ) )
      TBACK( I ) = CONJG( TBACK( I ) )
   END DO
   DO J = 1, NUMANG
      S1( J ) = CONJG( S1( J ) )
      S2( J ) = CONJG( S2( J ) )
   END DO
END IF

IF( PASS1 ) THEN
!                           ** Compare test case results with
!                           ** correct answers and abort if bad;
!                           ** otherwise restore user input and proceed
   CALL TESTMI( .TRUE., XX, CREFIN, MIMCUT, PERFCT, ANYANG, NMOM, &
                  IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, &
                  SBACK, S1, S2, TFORW, TBACK, PMOM, MOMDIM )
   PASS1  = .FALSE.
   CYCLE MAIN_MIEV
END IF
EXIT MAIN_MIEV
END DO MAIN_MIEV

IF( PRNT( 1 ) .OR. PRNT( 2 ) )  &
  CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT, &
               QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM, &
               SFORW, SBACK, TFORW, TBACK, S1, S2 )
!
DEALLOCATE(PIN)
DEALLOCATE(PINM1)
DEALLOCATE(SM) 
DEALLOCATE(SMS)
DEALLOCATE(SP) 
DEALLOCATE(SPS) 

RETURN

END SUBROUTINE MIEV0

SUBROUTINE CKINMI( NUMANG, MAXANG, XX, PERFCT, CREFIN, MOMDIM, NMOM, IPOLZN, ANYANG, XMU, CALCMO, NPQUAN )

!        Check for bad input to MIEV0 and calculate CALCMO, NPQUAN

!     Routines called :  ERRMSG, WRTBAD, WRTDIM

!     .. Scalar Arguments ..

LOGICAL     ANYANG, PERFCT
INTEGER     IPOLZN, MAXANG, MOMDIM, NMOM, NPQUAN, NUMANG
REAL(EB)    XX
COMPLEX(EB) CREFIN
!     ..
!     .. Array Arguments ..

LOGICAL     CALCMO( * )
REAL(EB)    XMU( * )
!     ..
!     .. Local Scalars ..

CHARACTER(4) :: STRING
LOGICAL     INPERR
INTEGER     I, IP, J, L

INPERR = .FALSE.

IF( NUMANG>MAXANG ) INPERR = WRTDIM( 'MaxAng', NUMANG )
IF( NUMANG<0 ) INPERR = WRTBAD( 'NUMANG' )

IF( XX<0._EB ) INPERR = WRTBAD( 'XX' )

IF( .NOT.PERFCT .AND. REAL( CREFIN )<=0._EB )INPERR = WRTBAD( 'CREFIN' )

IF( MOMDIM<0 ) INPERR = WRTBAD( 'MOMDIM' )


IF( NMOM/=0 ) THEN

   IF( NMOM<0 .OR. NMOM>MOMDIM ) INPERR = WRTBAD( 'NMOM' )
   IF( ABS( IPOLZN )>4444 ) INPERR = WRTBAD( 'IPOLZN' )

   NPQUAN = 0

   DO L = 1, 4
      CALCMO( L ) = .FALSE.
   END DO

   IF( IPOLZN/=0 ) THEN
   !                                 ** Parse out IPOLZN into its digits
   !                                 ** to find which phase quantities are
   !                                 ** to have their moments calculated
      WRITE( STRING, '(I4)' ) ABS( IPOLZN )
      DO J = 1, 4
         IP = ICHAR( STRING( J:J ) ) - ICHAR( '0' )
         IF( IP>=1 .AND. IP<=4 ) CALCMO( IP ) = .TRUE.
         IF( IP==0 .OR. ( IP>=5 .AND. IP<=9 ) ) INPERR = WRTBAD( 'IPOLZN' )
         NPQUAN = MAX( NPQUAN, IP )
      END DO
   END IF
END IF

IF( ANYANG ) THEN
!                                ** Allow for slight imperfections in
!                                ** computation of cosine
   DO I = 1, NUMANG
      IF( XMU( I )<-1.00001_EB .OR. XMU( I )>1.00001_EB ) INPERR = WRTBAD( 'XMU' )
   END DO
   ELSE
   DO I = 1, ( NUMANG + 1 ) / 2
      IF( XMU( I )<-0.00001_EB .OR. XMU( I )>1.00001_EB ) INPERR = WRTBAD( 'XMU' )
   END DO
END IF

IF( INPERR ) CALL ERRMSG( 'MIEV0--Input error(S).  Aborting...', .TRUE. )
IF( XX>20000.0_EB .OR. REAL( CREFIN )>10.0_EB .OR. ABS( AIMAG( CREFIN ) )>10.0_EB ) &
    CALL ERRMSG( 'MIEV0--XX or CREFIN outside tested range',.FALSE.)
RETURN
END SUBROUTINE CKINMI

SUBROUTINE LPCOEF( NTRM, NMOM, IPOLZN, MOMDIM, CALCMO, NPQUAN, A, B, PMOM )
!         Calculate Legendre polynomial expansion coefficients (also
!         called moments) for phase quantities ( Ref. 5 formulation )

!     INPUT:  NTRM                    Number terms in Mie series
!             NMOM, IPOLZN, MOMDIM    MIEV0 arguments
!             CALCMO                  Flags calculated from IPOLZN
!             NPQUAN                  Defined in MIEV0
!             A, B                    Mie series coefficients
!
!     OUTPUT: PMOM                   Legendre moments (MIEV0 argument)
!
!     Routines called :  ERRMSG, LPCO1T, LPCO2T
!
!     *** NOTES ***
!
!         (1)  Eqs. 2-5 are in error in Dave, Appl. Opt. 9,
!         1888 (1970).  Eq. 2 refers to M1, not M2;  eq. 3 refers to
!         M2, not M1.  In eqs. 4 and 5, the subscripts on the second
!         term in square brackets should be interchanged.
!
!         (2)  The general-case logic in this subroutine works correctly
!         in the two-term Mie series case, but subroutine LPCO2T
!         is called instead, for speed.
!
!         (3)  Subroutine  LPCO1T, to do the one-term case, is never
!         called within the context of MIEV0, but is included for
!         complete generality.
!
!         (4)  Some improvement in speed is obtainable by combining the
!         310- and 410-loops, if moments for both the third and fourth
!         phase quantities are desired, because the third phase quantity
!         is the real part of a complex series, while the fourth phase
!         quantity is the imaginary part of that very same series.  But
!         most users are not interested in the fourth phase quantity,
!         which is related to circular polarization, so the present
!         scheme is usually more efficient.
!
!
!           ** Definitions of local variables ***
!      AM(M)       Numerical coefficients  a-sub-m-super-l
!                     in Dave, Eqs. 1-15, as simplified in Ref. 5.
!      BI(I)       Numerical coefficients  b-sub-i-super-l
!                     in Dave, Eqs. 1-15, as simplified in Ref. 5.
!      BIDEL(I)    1/2 Bi(I) times factor capital-del in Dave
!      CM,DM()     Arrays C and D in Dave, Eqs. 16-17 (Mueller form),
!                     calculated using recurrence derived in Ref. 5
!      CS,DS()     Arrays C and D in Ref. 4, Eqs. A5-A6 (Sekera form),
!                     calculated using recurrence derived in Ref. 5
!      C,D()       Either CM,DM or CS,DS, depending on IPOLZN
!      EVENL       True for even-numbered moments;  false otherwise
!      IDEL        1 + little-del  in Dave
!      MAXTRM      Max. no. of terms in Mie series
!      MAXMOM      Max. no. of non-zero moments
!      NUMMOM      Number of non-zero moments
!      RECIP(K)    1 / K

IMPLICIT  NONE

!     .. Parameters ..

INTEGER     MAXTRM, MAXMOM, MXMOM2, MAXRCP
PARAMETER ( MAXTRM=1001, MAXMOM = 2*MAXTRM, MXMOM2 = MAXMOM / 2, MAXRCP = 4*MAXTRM + 2 )
!     ..
!     .. Scalar Arguments ..

INTEGER   IPOLZN, MOMDIM, NMOM, NPQUAN, NTRM
!     ..
!     .. Array Arguments ..

LOGICAL     CALCMO( * )
REAL(EB)    PMOM( 0:MOMDIM, * )
COMPLEX(EB) A( * ), B( * )
!     ..
!     .. Local Scalars ..

LOGICAL     EVENL, PASS1
INTEGER     I, IDEL, IMAX, J, K, L, LD2, M, MMAX, NUMMOM
REAL(EB)    SUM
!     ..
!     .. Local Arrays ..

REAL(EB)    AM( 0:MAXTRM ), BI( 0:MXMOM2 ), BIDEL( 0:MXMOM2 ), RECIP( MAXRCP )
COMPLEX(EB) C( MAXTRM ), CM( MAXTRM ), CS( MAXTRM ), D( MAXTRM ), DM( MAXTRM ), DS( MAXTRM )
!     ..
!     .. Equivalences ..

EQUIVALENCE ( C, CM ), ( D, DM )
!     ..
SAVE      PASS1, RECIP
DATA      PASS1 / .TRUE. /

IF( PASS1 ) THEN
   DO  K = 1, MAXRCP
      RECIP( K ) = 1.0_EB / K
   END DO
   PASS1  = .FALSE.
END IF

DO J = 1, MAX( 1, NPQUAN )
   DO L = 0, NMOM
      PMOM( L, J ) = 0.0_EB
   END DO
END DO

IF( NTRM==1 ) THEN
   CALL LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )
   RETURN
ELSE IF( NTRM==2 ) THEN
   CALL LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )
   RETURN
END IF

IF( NTRM + 2>MAXTRM )CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)
!                                     ** Calculate Mueller C, D arrays
CM( NTRM + 2 ) = ( 0._EB, 0._EB )
DM( NTRM + 2 ) = ( 0._EB, 0._EB )
CM( NTRM + 1 ) = ( 1._EB - RECIP( NTRM+1 ) ) * B( NTRM )
DM( NTRM + 1 ) = ( 1._EB - RECIP( NTRM+1 ) ) * A( NTRM )
CM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * A( NTRM ) + ( 1._EB - RECIP( NTRM ) )*B( NTRM-1 )
DM( NTRM ) = ( RECIP( NTRM ) + RECIP( NTRM+1 ) ) * B( NTRM ) + ( 1._EB - RECIP( NTRM ) )*A( NTRM-1 )
DO K = NTRM-1, 2, -1
   CM( K ) = CM( K+2 ) - ( 1._EB + RECIP(K+1) ) * B( K+1 )+ ( RECIP(K) + RECIP(K+1) ) * A( K ) &
                        + ( 1._EB - RECIP(K) ) * B( K-1 )
   DM( K ) = DM( K+2 ) - ( 1._EB + RECIP(K+1) ) * A( K+1 )  + ( RECIP(K) + RECIP(K+1) ) * B( K ) &
                        + ( 1._EB - RECIP(K) ) * A( K-1 )
END DO

CM( 1 ) = CM( 3 ) + 1.5_EB * ( A( 1 ) - B( 2 ) )
DM( 1 ) = DM( 3 ) + 1.5_EB * ( B( 1 ) - A( 2 ) )

IF( IPOLZN>=0 ) THEN
   DO K = 1, NTRM + 2
      C( K ) = ( 2*K - 1 ) * CM( K )
      D( K ) = ( 2*K - 1 ) * DM( K )
   END DO
ELSE
!                                    ** Compute Sekera C and D arrays
   CS( NTRM + 2 ) = ( 0._EB, 0._EB )
   DS( NTRM + 2 ) = ( 0._EB, 0._EB )
   CS( NTRM + 1 ) = ( 0._EB, 0._EB )
   DS( NTRM + 1 ) = ( 0._EB, 0._EB )

   DO K = NTRM, 1, -1
      CS( K ) = CS( K+2 ) + ( 2*K + 1 ) * ( CM( K+1 ) - B( K ) )
      DS( K ) = DS( K+2 ) + ( 2*K + 1 ) * ( DM( K+1 ) - A( K ) )
   END DO

   DO K = 1, NTRM + 2
      C( K ) = ( 2*K - 1 ) * CS( K )
      D( K ) = ( 2*K - 1 ) * DS( K )
   END DO

END IF

IF( IPOLZN<0 ) NUMMOM = MIN( NMOM, 2*NTRM - 2 )
IF( IPOLZN>=0 ) NUMMOM = MIN( NMOM, 2*NTRM )

IF( NUMMOM>MAXMOM )  CALL ERRMSG('LPCoef--PARAMETER MaxTrm too small',.TRUE.)
!
!                          ** Loop over moments
L240: DO L = 0, NUMMOM
   LD2 = L / 2
   EVENL  = MOD( L, 2 )==0
!                                    ** Calculate numerical coefficients
!                                    ** a-sub-m and b-sub-i in Dave
!                                    ** double-sums for moments
   IF( L==0 ) THEN
      IDEL = 1
      DO M = 0, NTRM
         AM( M ) = 2.0_EB * RECIP( 2*M + 1 )
      END DO
      BI( 0 ) = 1.0_EB
   ELSE IF( EVENL ) THEN
      IDEL = 1
      DO M = LD2, NTRM
         AM( M ) = ( 1._EB + RECIP( 2*M - L + 1 ) ) * AM( M )
      END DO
      DO I = 0, LD2 - 1
         BI( I ) = ( 1._EB - RECIP( L - 2*I ) ) * BI( I )
      END DO
      BI( LD2 ) = ( 2._EB - RECIP( L ) ) * BI( LD2 - 1 )
   ELSE
   IDEL = 2
      DO M = LD2, NTRM
         AM( M ) = ( 1._EB - RECIP( 2*M + L + 2 ) ) * AM( M )
      END DO

      DO I = 0, LD2
         BI( I ) = ( 1._EB - RECIP( L + 2*I + 1 ) ) * BI( I )
      END DO
   END IF
!                                     ** Establish upper limits for sums
!                                     ** and incorporate factor capital-
!                                     ** del into b-sub-i
   MMAX = NTRM - IDEL
   IF( IPOLZN>=0 ) MMAX = MMAX + 1
   IMAX = MIN( LD2, MMAX - LD2 )
   IF( IMAX<0 ) EXIT L240
   DO I = 0, IMAX
      BIDEL( I ) = BI( I )
   END DO
   IF( EVENL ) BIDEL( 0 ) = 0.5_EB*BIDEL( 0 )
!                                    ** Perform double sums just for
!                                    ** phase quantities desired by user
   IF( IPOLZN==0 ) THEN
      DO I = 0, IMAX
!                                           ** vectorizable loop
         SUM = 0.0_EB
         DO M = LD2, MMAX - I
            SUM = SUM + AM(M) * ( REAL(C(M-I+1) * CONJG(C(M+I+IDEL)))+ REAL(D(M-I+1) * CONJG(D(M+I+IDEL))))
         END DO
         PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM
      END DO
      PMOM( L, 1 ) = 0.5_EB*PMOM( L, 1 )
      CYCLE L240
   END IF

   IF( CALCMO( 1 ) ) THEN
      DO I = 0, IMAX
         SUM = 0.0_EB
!                                          ** vectorizable loop
         DO M = LD2, MMAX - I
            SUM = SUM + AM( M ) * REAL( C(M-I+1) * CONJG( C(M+I+IDEL) ) )
         END DO
         PMOM( L, 1 ) = PMOM( L, 1 ) + BIDEL( I ) * SUM
      END DO
   END IF

   IF( CALCMO( 2 ) ) THEN
      DO I = 0, IMAX
         SUM = 0.0_EB
!                                           ** vectorizable loop
         DO M = LD2, MMAX - I
            SUM = SUM + AM( M ) * REAL( D(M-I+1) * CONJG( D(M+I+IDEL) ) )
         END DO
         PMOM( L, 2 ) = PMOM( L, 2 ) + BIDEL( I ) * SUM
      END DO
   END IF

   IF( CALCMO( 3 ) ) THEN
      DO I = 0, IMAX
         SUM = 0.0_EB
!                                           ** vectorizable loop
         DO M = LD2, MMAX - I
            SUM = SUM + AM(M) *( REAL(C(M-I+1) * CONJG(D(M+I+IDEL))) + REAL(C(M+I+IDEL) * CONJG(D(M-I+1))))
         END DO
         PMOM( L, 3 ) = PMOM( L, 3 ) + BIDEL( I ) * SUM
      END DO
      PMOM( L, 3 ) = 0.5_EB*PMOM( L, 3 )
   END IF

   IF( CALCMO( 4 ) ) THEN
      DO I = 0, IMAX
         SUM= 0.0_EB
!                                         ** vectorizable loop
         DO M = LD2, MMAX - I
            SUM = SUM + AM(M) *(AIMAG(C(M-I+1) * CONJG(D(M+I+IDEL)))+ AIMAG(C(M+I+IDEL) * CONJG(D(M-I+1))))
         END DO
         PMOM( L, 4 ) = PMOM( L, 4 ) + BIDEL( I ) * SUM
      END DO

      PMOM( L, 4 ) = - 0.5_EB * PMOM( L, 4 )
   END IF

END DO L240

RETURN
END SUBROUTINE LPCOEF

SUBROUTINE LPCO1T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )
!
!         Calculate Legendre polynomial expansion coefficients (also
!         called moments) for phase quantities in special case where
!         no. terms in Mie series = 1
!
!        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
!                CALCMO                   Flags calculated from IPOLZN
!                A(1), B(1)               Mie series coefficients
!
!        OUTPUT: PMOM                     Legendre moments

!     .. Scalar Arguments ..

INTEGER     IPOLZN, MOMDIM, NMOM
!     ..
!     .. Array Arguments ..

LOGICAL     CALCMO( * )
REAL(EB)    PMOM( 0:MOMDIM, * )
COMPLEX(EB) A( * ), B( * )
!     ..
!     .. Local Scalars ..

INTEGER     L, NUMMOM
REAL(EB)    A1SQ, B1SQ
COMPLEX(EB) A1B1C

A1SQ   = SQ( A( 1 ) )
B1SQ   = SQ( B( 1 ) )
A1B1C  = A( 1 ) * CONJG( B( 1 ) )


IF( IPOLZN<0 ) THEN
   IF( CALCMO( 1 ) ) PMOM( 0, 1 ) = 2.25_EB*B1SQ
   IF( CALCMO( 2 ) ) PMOM( 0, 2 ) = 2.25_EB*A1SQ
   IF( CALCMO( 3 ) ) PMOM( 0, 3 ) = 2.25_EB*REAL( A1B1C )
   IF( CALCMO( 4 ) ) PMOM( 0, 4 ) = 2.25_EB*AIMAG( A1B1C )
ELSE
   NUMMOM = MIN( NMOM, 2 )
   !                             ** Loop over moments
   L10: DO L = 0, NUMMOM
      IF( IPOLZN==0 ) THEN
         IF( L==0 ) PMOM( L, 1 ) = 1.5_EB*( A1SQ + B1SQ )
         IF( L==1 ) PMOM( L, 1 ) = 1.5_EB*REAL( A1B1C )
         IF( L==2 ) PMOM( L, 1 ) = 0.15_EB*( A1SQ + B1SQ )
         CYCLE L10
      END IF

      IF( CALCMO( 1 ) ) THEN
         IF( L==0 ) PMOM( L, 1 ) = 2.25_EB*( A1SQ + B1SQ / 3.)
         IF( L==1 ) PMOM( L, 1 ) = 1.5_EB*REAL( A1B1C )
         IF( L==2 ) PMOM( L, 1 ) = 0.3_EB*B1SQ
      END IF

      IF( CALCMO( 2 ) ) THEN
         IF( L==0 ) PMOM( L, 2 ) = 2.25_EB*( B1SQ + A1SQ / 3. )
         IF( L==1 ) PMOM( L, 2 ) = 1.5_EB*REAL( A1B1C )
         IF( L==2 ) PMOM( L, 2 ) = 0.3_EB*A1SQ
      END IF

      IF( CALCMO( 3 ) ) THEN
         IF( L==0 ) PMOM( L, 3 ) = 3.0_EB*REAL( A1B1C )
         IF( L==1 ) PMOM( L, 3 ) = 0.75_EB*( A1SQ + B1SQ )
         IF( L==2 ) PMOM( L, 3 ) = 0.3_EB*REAL( A1B1C )
      END IF

      IF( CALCMO( 4 ) ) THEN
         IF( L==0 ) PMOM( L, 4 ) = -1.5_EB*AIMAG( A1B1C )
         IF( L==1 ) PMOM( L, 4 ) = 0.0_EB
         IF( L==2 ) PMOM( L, 4 ) = 0.3_EB*AIMAG( A1B1C )
      END IF

   END DO L10
END IF

RETURN
END SUBROUTINE LPCO1T

SUBROUTINE LPCO2T( NMOM, IPOLZN, MOMDIM, CALCMO, A, B, PMOM )
!         Calculate Legendre polynomial expansion coefficients (also
!         called moments) for phase quantities in special case where
!         no. terms in Mie series = 2

!        INPUT:  NMOM, IPOLZN, MOMDIM     MIEV0 arguments
!                CALCMO                   Flags calculated from IPOLZN
!                A(1-2), B(1-2)           Mie series coefficients
!
!        OUTPUT: PMOM                     Legendre moments


!     .. Scalar Arguments ..

INTEGER     IPOLZN, MOMDIM, NMOM
!     ..
!     .. Array Arguments ..

LOGICAL     CALCMO( * )
REAL(EB)    PMOM( 0:MOMDIM, * )
COMPLEX(EB) A( * ), B( * )
!     ..
!     .. Local Scalars ..

INTEGER     L, NUMMOM
REAL(EB)    A2SQ, B2SQ, PM1, PM2
COMPLEX(EB) A2C, B2C, CA, CAC, CAT, CB, CBC, CBT, CG, CH

CA   = 3.*A( 1 ) - 5.*B( 2 )
CAT  = 3.*B( 1 ) - 5.*A( 2 )
CAC  = CONJG( CA )
A2SQ = SQ( A( 2 ) )
B2SQ = SQ( B( 2 ) )
A2C  = CONJG( A( 2 ) )
B2C  = CONJG( B( 2 ) )

IF( IPOLZN<0 ) THEN
!                                   ** Loop over Sekera moments
   NUMMOM = MIN( NMOM, 2 )
   DO L = 0, NUMMOM
      IF( CALCMO( 1 ) ) THEN
         IF( L==0 ) PMOM( L, 1 ) = 0.25_EB * ( SQ( CAT ) + (100._EB/3._EB)* B2SQ )
         IF( L==1 ) PMOM( L, 1 ) = (5._EB/3._EB)*REAL( CAT*B2C )
         IF( L==2 ) PMOM( L, 1 ) = (10._EB/3._EB)*B2SQ
      END IF

      IF( CALCMO( 2 ) ) THEN
         IF( L==0 ) PMOM( L, 2 ) = 0.25_EB * ( SQ( CA ) + (100._EB/3._EB) * A2SQ )
         IF( L==1 ) PMOM( L, 2 ) = (5._EB/3._EB)*REAL( CA*A2C )
         IF( L==2 ) PMOM( L, 2 ) = (10._EB/3._EB)*A2SQ
      END IF

      IF( CALCMO( 3 ) ) THEN
         IF( L==0 ) PMOM( L, 3 ) = 0.25_EB * REAL( CAT * CAC  + (100._EB/3._EB) * B(2) * A2C )
         IF( L==1 ) PMOM( L, 3 ) = 5._EB/6._EB* REAL( B(2)*CAC + CAT*A2C )
         IF( L==2 ) PMOM( L, 3 ) = 10._EB/3._EB* REAL( B(2)*A2C )
      END IF

      IF( CALCMO( 4 ) ) THEN
         IF( L==0 ) PMOM( L, 4 ) = -0.25_EB * AIMAG( CAT * CAC + (100._EB/3._EB)* B(2) * A2C )
         IF( L==1 ) PMOM( L, 4 ) = -5._EB/ 6._EB*  AIMAG( B(2)*CAC + CAT*A2C )
         IF( L==2 ) PMOM( L, 4 ) = -10._EB/ 3._EB* AIMAG( B(2)*A2C )
      END IF

   END DO

ELSE

   CB  = 3._EB*B( 1 ) + 5._EB*A( 2 )
   CBT = 3._EB*A( 1 ) + 5._EB*B( 2 )
   CBC = CONJG( CB )
   CG  = ( CBC*CBT + 10._EB*( CAC*A( 2 ) + B2C*CAT ) ) / 3._EB
   CH  = 2._EB*( CBC*A( 2 ) + B2C*CBT )

   !                               ** Loop over Mueller moments
   NUMMOM = MIN( NMOM, 4 )

   L20: DO L = 0, NUMMOM

      IF( IPOLZN==0 .OR. CALCMO( 1 ) ) THEN
         IF( L==0 ) PM1 = 0.25_EB*SQ( CA ) + SQ( CB ) / 12._EB + (5._EB/3._EB)*REAL( CA*B2C ) + 5._EB*B2SQ
         IF( L==1 ) PM1 = REAL( CB * ( CAC / 6._EB+ B2C ) )
         IF( L==2 ) PM1 = SQ( CB ) / 30._EB+ (20._EB/7._EB)*B2SQ + (2._EB/3._EB)*REAL( CA*B2C )
         IF( L==3 ) PM1 = (2._EB/7.) * REAL( CB*B2C )
         IF( L==4 ) PM1 = (40._EB/63._EB) * B2SQ
         IF( CALCMO( 1 ) ) PMOM( L, 1 ) = PM1
      END IF

      IF( IPOLZN==0 .OR. CALCMO( 2 ) ) THEN
         IF( L==0 ) PM2 = 0.25_EB*SQ( CAT ) + SQ( CBT ) / 12._EB + ( 5._EB/ 3._EB) * REAL( CAT*A2C ) + 5._EB*A2SQ
         IF( L==1 ) PM2 = REAL( CBT * ( CONJG( CAT ) / 6._EB+ A2C ) )
         IF( L==2 ) PM2 = SQ( CBT ) / 30._EB + ( 20._EB/7._EB) * A2SQ + ( 2._EB/3._EB) * REAL( CAT*A2C )
         IF( L==3 ) PM2 = (2._EB/7._EB) * REAL( CBT*A2C )
         IF( L==4 ) PM2 = (40._EB/63._EB) * A2SQ
         IF( CALCMO( 2 ) ) PMOM( L, 2 ) = PM2
      END IF

      IF( IPOLZN==0 ) THEN
         PMOM( L, 1 ) = 0.5_EB*( PM1 + PM2 )
         CYCLE L20
      END IF

      IF( CALCMO( 3 ) ) THEN
         IF( L==0 ) PMOM( L, 3 ) = 0.25_EB * REAL( CAC*CAT + CG + 20._EB* B2C * A(2) )
         IF( L==1 ) PMOM( L, 3 ) = REAL( CAC*CBT + CBC*CAT + 3._EB*CH ) / 12.
         IF( L==2 ) PMOM( L, 3 ) = 0.1_EB * REAL( CG + (200._EB/7._EB) * B2C * A(2) )
         IF( L==3 ) PMOM( L, 3 ) = REAL( CH ) / 14._EB
         IF( L==4 ) PMOM( L, 3 ) = 40._EB/63._EB* REAL( B2C*A(2) )
      END IF

      IF( CALCMO( 4 ) ) THEN
         IF( L==0 ) PMOM( L, 4 ) = 0.25_EB * AIMAG( CAC*CAT + CG + 20._EB* B2C * A(2) )
         IF( L==1 ) PMOM( L, 4 ) = AIMAG( CAC*CBT + CBC*CAT + 3._EB*CH ) / 12._EB
         IF( L==2 ) PMOM( L, 4 ) = 0.1_EB * AIMAG( CG + (200._EB/7._EB) * B2C * A(2) )
         IF( L==3 ) PMOM( L, 4 ) = AIMAG( CH ) / 14._EB
         IF( L==4 ) PMOM( L, 4 ) = 40._EB/63._EB* AIMAG( B2C*A(2) )
      END IF

   END DO L20

END IF

RETURN
END SUBROUTINE LPCO2T

SUBROUTINE BIGA( CIOR, XX, NTRM, NOABS, YESANG, RBIGA, CBIGA )

!        Calculate logarithmic derivatives of J-Bessel-function
!
!     Input :  CIOR, XX, NTRM, NOABS, YESANG  (defined in MIEV0)
!
!    Output :  RBIGA or CBIGA  (defined in MIEV0)
!
!    Routines called :  CONFRA
!
!
!    INTERNAL VARIABLES :
!
!       CONFRA     Value of Lentz continued fraction for CBIGA(NTRM),
!                     used to initialize downward recurrence
!
!       DOWN       = True, use down-recurrence.  False, do not.
!
!       F1,F2,F3   Arithmetic statement functions used in determining
!                     whether to use up-  or down-recurrence
!                     ( Ref. 2, Eqs. 6-8 )
!
!       MRE        Real refractive index
!       MIM        Imaginary refractive index
!
!       REZINV     1 / ( MRE * XX ); temporary variable for recurrence
!       ZINV       1 / ( CIOR * XX ); temporary variable for recurrence
!
!     .. Scalar Arguments ..

LOGICAL     NOABS, YESANG
INTEGER     NTRM
REAL(EB)    XX
COMPLEX(EB) CIOR

!     .. Array Arguments ..
REAL(EB)    RBIGA( * )
COMPLEX(EB) CBIGA( * )

!     .. Local Scalars ..
LOGICAL     DOWN
INTEGER     N
REAL(EB)    MIM, MRE, REZINV, RTMP
COMPLEX(EB) CTMP, ZINV

!                                  ** Decide whether BigA can be
!                                  ** calculated by up-recurrence
MRE = REAL( CIOR )
MIM = ABS( AIMAG( CIOR ) )
IF( MRE<1.0_EB .OR. MRE>10.0_EB .OR. MIM>10.0_EB ) THEN
   DOWN = .TRUE.
ELSE IF( YESANG ) THEN
   DOWN = .TRUE.
!                                                    ** Eq. R48
   IF( MIM*XX < F2( MRE ) ) DOWN = .FALSE.
ELSE
   DOWN = .TRUE.
!                                                    ** Eq. R48
   IF( MIM*XX < F1( MRE ) ) DOWN = .FALSE.
END IF

ZINV   = 1.0_EB / ( CIOR*XX )
REZINV = 1.0_EB / ( MRE*XX )


IF( DOWN ) THEN
!                          ** Compute initial high-order BigA using
!                          ** Lentz method ( Ref. 1, pp. 17-20 )
   CTMP = CONFRA( NTRM, ZINV )
!                                   *** Downward recurrence for BigA
   IF( NOABS ) THEN
!                                        ** No-absorption case; Eq (R23)
      RBIGA( NTRM ) = REAL( CTMP )
      DO N = NTRM, 2, -1
         RBIGA( N - 1 ) = ( N*REZINV ) - 1.0_EB / ( ( N*REZINV ) + RBIGA( N ) )
      END DO
   ELSE
!                                         ** Absorptive case; Eq (R23)
      CBIGA( NTRM ) = CTMP
      DO N = NTRM, 2, -1
         CBIGA( N-1 ) = (N*ZINV) - 1.0_EB / ( (N*ZINV) + CBIGA( N ) )
   END DO
   END IF

ELSE
!                            *** Upward recurrence for BigA
   IF( NOABS ) THEN
!                                  ** No-absorption case; Eq (R20,21)
      RTMP = SIN( MRE*XX )
      RBIGA( 1 ) = - REZINV + RTMP / ( RTMP*REZINV - COS( MRE*XX ) )

      DO N = 2, NTRM
         RBIGA( N ) = -( N*REZINV ) + 1.0_EB / ( ( N*REZINV ) - RBIGA( N - 1 ) )
      END DO

   ELSE
!                                     ** Absorptive case; Eq (R20,22)
      CTMP = EXP( - (0._EB,2._EB)*CIOR*XX )
      CBIGA( 1 ) = - ZINV + (1._EB-CTMP) /  ( ZINV * (1._EB-CTMP) - (0._EB,1._EB)*(1._EB+CTMP) )
      DO N = 2, NTRM
         CBIGA( N ) = - (N*ZINV) + 1.0_EB / ((N*ZINV) - CBIGA( N-1 ))
      END DO

END IF

END IF

CONTAINS

REAL(EB) FUNCTION F1(MRE)
REAL(EB), INTENT(IN) :: MRE
F1 = -8.0_EB + MRE**2*( 26.22_EB + MRE*( -0.4474_EB + MRE**3*( 0.00204_EB - 0.000175_EB*MRE ) ) )
END FUNCTION F1

REAL(EB) FUNCTION F2(MRE)  ! Eq. R47b
REAL(EB), INTENT(IN) :: MRE
F2 = 3.9_EB + MRE*( -10.8_EB + 13.78_EB*MRE )
END FUNCTION F2

END SUBROUTINE BIGA

COMPLEX(EB) FUNCTION CONFRA( N, ZINV )

!         Compute Bessel function ratio A-sub-N from its
!         continued fraction using Lentz method
!
!         ZINV = Reciprocal of argument of A
!
!
!    I N T E R N A L    V A R I A B L E S
!    ------------------------------------
!
!    CAK      Term in continued fraction expansion of A (Eq. R25)
!    CAPT     Factor used in Lentz iteration for A (Eq. R27)
!    CNUMER   Numerator   in capT  ( Eq. R28A )
!    CDENOM   Denominator in capT  ( Eq. R28B )
!    CDTD     Product of two successive denominators of capT factors
!                 ( Eq. R34C )
!    CNTN     Product of two successive numerators of capT factors
!                 ( Eq. R34B )
!    EPS1     Ill-conditioning criterion
!    EPS2     Convergence criterion
!    KK       Subscript k of cAk  ( Eq. R25B )
!    KOUNT    Iteration counter ( used to prevent infinite looping )
!    MAXIT    Max. allowed no. of iterations
!    MM       + 1  and - 1, alternately
! --------------------------------------------------------------------

!     .. Scalar Arguments ..

INTEGER     N
COMPLEX(EB) ZINV

!     .. Local Scalars ..

INTEGER     KK, KOUNT, MAXIT, MM
REAL(EB)    EPS1, EPS2
COMPLEX(EB) CAK, CAPT, CDENOM, CDTD, CNTN, CNUMER
!    
DATA      EPS1 / 1.E-2_EB / , EPS2 / 1.E-8_EB /
DATA      MAXIT / 10000 /
!                                 ** Eq. R25a
CONFRA = ( N + 1 ) * ZINV
MM     = - 1
KK     = 2*N + 3
!                                 ** Eq. R25b, k=2
CAK    = ( MM*KK ) * ZINV
CDENOM = CAK
CNUMER = CDENOM + 1.0_EB / CONFRA
KOUNT  = 1

L10: DO
   KOUNT = KOUNT + 1

   IF( KOUNT>MAXIT ) CALL ERRMSG('ConFra--Iteration failed to converge',.TRUE.)

   MM  = - MM
   KK  = KK + 2
   !                                 ** Eq. R25b
   CAK = ( MM*KK ) * ZINV
   !                                          ** Eq. R32
   IF( ABS( CNUMER / CAK )<=EPS1 .OR. ABS( CDENOM / CAK )<=EPS1 ) THEN

   !                                  ** Ill-conditioned case -- stride
   !                                  ** two terms instead of one

   !                                       ** Eq. R34
      CNTN   = CAK * CNUMER + 1.0_EB
      CDTD   = CAK * CDENOM + 1.0_EB
   !                                           ** Eq. R33
      CONFRA = ( CNTN / CDTD ) * CONFRA

      MM  = - MM
      KK  = KK + 2
   !                                 ** Eq. R25b
      CAK = ( MM*KK ) * ZINV
   !                                      ** Eq. R35
      CNUMER = CAK + CNUMER / CNTN
      CDENOM = CAK + CDENOM / CDTD
      KOUNT  = KOUNT + 1
      CYCLE L10

   ELSE
   !                           *** Well-conditioned case

   !                                  ** Eq. R27
      CAPT   = CNUMER / CDENOM
   !                                  ** Eq. R26
      CONFRA = CAPT * CONFRA
   !                                  ** Check for convergence; Eq. R31

      IF (      ABS( REAL (CAPT) - 1.0_EB )>=EPS2 .OR. ABS( AIMAG(CAPT) )      >=EPS2 )  THEN
   !                                        ** Eq. R30
         CNUMER = CAK + 1.0_EB / CNUMER
         CDENOM = CAK + 1.0_EB / CDENOM
         CYCLE L10
      END IF
   END IF
   EXIT L10
END DO L10

RETURN
END FUNCTION CONFRA

SUBROUTINE MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT, QSCA, GQSC, NMOM, IPOLZN, MOMDIM, &
                   CALCMO, PMOM, SFORW, SBACK, TFORW, TBACK, S1, S2 )

!         Print scattering quantities of a single particle

!     .. Scalar Arguments ..
LOGICAL     PERFCT
INTEGER     IPOLZN, MOMDIM, NMOM, NUMANG
REAL(EB)    GQSC, QEXT, QSCA, XX
COMPLEX(EB) CREFIN, SBACK, SFORW

!     .. Array Arguments ..
LOGICAL     CALCMO( * ), PRNT( * )
REAL(EB)    PMOM( 0:MOMDIM, * ), XMU( * )
COMPLEX(EB) S1( * ), S2( * ), TBACK( * ), TFORW( * )

!     .. Local Scalars ..
CHARACTER   FMAT*22
INTEGER     I, J, M
REAL(EB)    FNORM, I1, I2

!     .. Intrinsic Functions ..
!     INTRINSIC AIMAG, CONJG, REAL


IF( PERFCT ) WRITE(LU_ERR, '(''1'',10X,A,1P,E11.4)' )  'Perfectly Conducting Case, size parameter =', XX

IF( .NOT.PERFCT ) WRITE(LU_ERR, '(''1'',10X,3(A,1P,E11.4))' ) &
    'Refractive Index:  Real ', REAL( CREFIN ), '  Imag ', &
    AIMAG( CREFIN ), ',   Size Parameter =', XX


IF( PRNT( 1 ) .AND. NUMANG>0 ) THEN

   WRITE(LU_ERR, '(/,A)' ) &
         '    cos(angle)  ------- S1 ---------  ------- S2 ---------' &
         // '  --- S1*conjg(S2) ---   i1=S1**2   i2=S2**2  (i1+i2)/2' &
         // '  DEG POLZN'

   DO I = 1, NUMANG
      I1 = REAL( S1( I ) )**2 + AIMAG( S1( I ) )**2
      I2 = REAL( S2( I ) )**2 + AIMAG( S2( I ) )**2
      WRITE(LU_ERR, '( I4, F10.6, 1P,10E11.3 )'   ) I, XMU(I), S1(I), S2(I), &
         S1(I)*CONJG(S2(I)),I1, I2, 0.5_EB*(I1+I2), (I2-I1)/(I2+I1)
   END DO
END IF


IF( PRNT( 2 ) ) THEN

   WRITE (LU_ERR, '(/,A,9X,A,17X,A,17X,A,/,(0P,F7.2, 1P,6E12.3) )' ) &
            '  Angle', 'S-sub-1', 'T-sub-1', 'T-sub-2', &
                  0.0_EB,     SFORW,    TFORW(1),  TFORW(2), &
               180._EB,     SBACK,    TBACK(1),  TBACK(2)
   WRITE (LU_ERR, '(/,4(A,1P,E11.4))' ) &
            ' Efficiency Factors,  extinction:', QEXT, &
                                 '   scattering:', QSCA, &
                                 '   absorption:', QEXT-QSCA, &
                              '   rad. pressure:', QEXT-GQSC
   IF( NMOM>0 ) THEN
      WRITE(LU_ERR, '(/,A)' ) ' Normalized moments of :'
      IF( IPOLZN==0 ) WRITE(LU_ERR, '(''+'',27X,A)' ) 'Phase Fcn'
      IF( IPOLZN>0 ) WRITE(LU_ERR, '(''+'',33X,A)' )  'M1           M2          S21          D21'

      IF( IPOLZN<0 ) WRITE(LU_ERR, '(''+'',33X,A)' )  'R1           R2           R3           R4'
      FNORM = 4._EB/ ( XX**2 * QSCA )
      DO M = 0, NMOM
         WRITE(LU_ERR, '(A,I4)' ) '      Moment no.', M
         DO J = 1, 4
            IF( CALCMO( J ) ) THEN
               WRITE( FMAT, '(A,I2,A)' ) '( ''+'', T', 24+(J-1)*13, ', 1P,E13.4 )'
               WRITE(LU_ERR, FMAT ) FNORM * PMOM( M, J )
            END IF
         END DO
      END DO
   END IF
END IF

RETURN
END SUBROUTINE MIPRNT

SUBROUTINE SMALL1( XX, NUMANG, XMU, QEXT, QSCA, GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK, A, B )

!       Small-particle limit of Mie quantities in totally reflecting
!       limit ( Mie series truncated after 2 terms )
!
!        A,B       First two Mie coefficients, with numerator and
!                  denominator expanded in powers of XX ( a factor
!                  of XX**3 is missing but is restored before return
!                  to calling program )  ( Ref. 2, p. 1508 )

!     .. Parameters ..
REAL(EB)  TWOTHR, FIVTHR, FIVNIN
PARAMETER ( TWOTHR = 2._EB/3._EB, FIVTHR = 5._EB/3._EB, FIVNIN = 5._EB/9._EB )

!     .. Scalar Arguments ..
INTEGER     NUMANG
REAL(EB)    GQSC, QEXT, QSCA, XX
COMPLEX(EB) SBACK, SFORW

!     .. Array Arguments ..
REAL(EB)    XMU( * )
COMPLEX(EB) A(*), B( * ), S1( * ), S2( * ), TBACK( * ), TFORW( * )

!     .. Local Scalars ..
INTEGER     J
REAL(EB)    RTMP

!                                                       ** Eq. R/A.5
A( 1 ) = CMPLX( 0._EB, TWOTHR*( 1._EB - 0.2_EB*XX**2 ) ) / CMPLX( 1._EB - 0.5_EB*XX**2, TWOTHR*XX**3 )
!                                                      ** Eq. R/A.6
B( 1 ) = CMPLX( 0._EB, - ( 1._EB - 0.1_EB*XX**2 ) / 3._EB) / CMPLX( 1._EB + 0.5_EB*XX**2, - XX**3 / 3._EB)
!                                                       ** Eq. R/A.7,8
A( 2 ) = CMPLX( 0._EB,   XX**2 / 30._EB)
B( 2 ) = CMPLX( 0._EB, - XX**2 / 45._EB)
!                                                       ** Eq. R/A.9
QSCA = 6._EB* XX**4 *( SQ( A(1) ) + SQ( B(1) ) + FIVTHR*( SQ( A(2) ) + SQ( B(2) ) ) )
QEXT = QSCA
!                                                       ** Eq. R/A.10
GQSC = 6._EB* XX**4 *REAL( A(1)*CONJG( A(2) + B(1) ) + ( B(1) + FIVNIN*A(2) )*CONJG( B(2) ) )

RTMP   = 1.5_EB * XX**3
SFORW  = RTMP*( A(1) + B(1) + FIVTHR*( A(2) + B(2) ) )
SBACK  = RTMP*( A(1) - B(1) - FIVTHR*( A(2) - B(2) ) )
TFORW( 1 ) = RTMP*( B(1) + FIVTHR*( 2._EB*B(2) - A(2) ) )
TFORW( 2 ) = RTMP*( A(1) + FIVTHR*( 2._EB*A(2) - B(2) ) )
TBACK( 1 ) = RTMP*( B(1) - FIVTHR*( 2._EB*B(2) + A(2) ) )
TBACK( 2 ) = RTMP*( A(1) - FIVTHR*( 2._EB*A(2) + B(2) ) )

DO J = 1, NUMANG
!                                                    ** Eq. R/A.11,12
   S1( J ) = RTMP*( A(1) + B(1)*XMU( J ) + FIVTHR*( A(2)*XMU( J ) +  B(2)*( 2._EB*XMU( J )**2 - 1._EB) ) )
   S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) + FIVTHR*( B(2)*XMU( J ) +  A(2)*( 2._EB*XMU( J )**2 - 1._EB) ) )
END DO
!                                     ** Recover actual Mie coefficients
A( 1 ) = XX**3 * A(1)
A( 2 ) = XX**3 * A(2)
B( 1 ) = XX**3 * B(1)
B( 2 ) = XX**3 * B(2)

RETURN
END SUBROUTINE SMALL1

SUBROUTINE SMALL2( XX, CIOR, CALCQE, NUMANG, XMU, QEXT, QSCA,GQSC, SFORW, SBACK, S1, S2, TFORW, TBACK,A, B )

!       Small-particle limit of Mie quantities for general refractive
!       index ( Mie series truncated after 2 terms )
!
!        A,B       First two Mie coefficients, with numerator and
!                  denominator expanded in powers of XX ( a factor
!                  of XX**3 is missing but is restored before return
!                  to calling program )

!        CIORSQ    Square of refractive index

!     .. Parameters ..
REAL(EB)    TWOTHR, FIVTHR
PARAMETER   ( TWOTHR = 2._EB/3._EB, FIVTHR = 5._EB/3._EB)

!     .. Scalar Arguments ..
LOGICAL     CALCQE
INTEGER     NUMANG
REAL(EB)    GQSC, QEXT, QSCA, XX
COMPLEX(EB) CIOR, SBACK, SFORW

!     .. Array Arguments ..
REAL(EB)    XMU( * )
COMPLEX(EB) A(*), B(*), S1(*), S2(*), TBACK(*), TFORW( * )

!     .. Local Scalars ..
INTEGER   J
REAL(EB)  RTMP
COMPLEX(EB)   CIORSQ, CTMP


CIORSQ = CIOR**2
CTMP   = CMPLX( 0._EB, TWOTHR )*( CIORSQ - 1.0_EB )

!                                           ** Eq. R42a
A( 1 ) = CTMP*( 1._EB- 0.1_EB*XX**2 +   ( CIORSQ / 350._EB + 1._EB/280._EB)*XX**4 ) / &
          ( CIORSQ + 2._EB+ ( 1._EB- 0.7_EB*CIORSQ )*XX**2 -  &
          ( CIORSQ**2 / 175._EB- 0.275_EB*CIORSQ + 0.25_EB )*XX**4 + &
          XX**3 * CTMP * ( 1._EB- 0.1_EB*XX**2 ) )

!                                           ** Eq. R42b
B( 1 ) = ( XX**2 / 30._EB )*CTMP*( 1._EB+  ( CIORSQ / 35._EB - 1._EB/ 14._EB)*XX**2 ) /  &
         ( 1._EB- ( CIORSQ / 15._EB - 1._EB/6._EB)*XX**2 )

!                                           ** Eq. R42c

A( 2 ) = ( 0.1_EB*XX**2 )*CTMP*( 1._EB- XX**2 / 14._EB ) /  ( 2._EB*CIORSQ + 3._EB- &
         ( CIORSQ / 7._EB- 0.5_EB ) * XX**2 )

!                                           ** Eq. R40a

QSCA = (6._EB*XX**4) * ( SQ( A(1) ) + SQ( B(1) ) + FIVTHR * SQ( A(2) ) )

!                                           ** Eq. R40b
QEXT = QSCA
IF( CALCQE ) QEXT = 6._EB*XX * REAL( A(1) + B(1) + FIVTHR*A(2) )

!                                           ** Eq. R40c

GQSC = (6._EB*XX**4) * REAL( A(1)*CONJG( A(2) + B(1) ) )

RTMP   = 1.5_EB * XX**3
SFORW  = RTMP*( A(1) + B(1) + FIVTHR*A(2) )
SBACK  = RTMP*( A(1) - B(1) - FIVTHR*A(2) )
TFORW( 1 ) = RTMP*( B(1) - FIVTHR*A(2) )
TFORW( 2 ) = RTMP*( A(1) + 2._EB*FIVTHR*A(2) )
TBACK( 1 ) = TFORW(1)
TBACK( 2 ) = RTMP*( A(1) - 2._EB*FIVTHR*A(2) )


DO J = 1, NUMANG
!                                      ** Eq. R40d,e

   S1( J ) = RTMP*( A(1) + ( B(1) + FIVTHR*A(2) )*XMU( J ) )
   S2( J ) = RTMP*( B(1) + A(1)*XMU( J ) + FIVTHR*A(2)*( 2._EB*XMU( J )**2 - 1._EB) )
END DO

!                                     ** Recover actual Mie coefficients
A( 1 ) = XX**3 * A(1)
A( 2 ) = XX**3 * A(2)
B( 1 ) = XX**3 * B(1)
B( 2 ) = ( 0._EB, 0._EB)

RETURN
END SUBROUTINE SMALL2

SUBROUTINE TESTMI( COMPAR, XX, CREFIN, MIMCUT, PERFCT, ANYANG, &
                   NMOM, IPOLZN, NUMANG, XMU, QEXT, QSCA, GQSC, &
                   SFORW, SBACK, S1, S2, TFORW, TBACK, PMOM, &
                   MOMDIM )

!         Set up to run test case when  COMPAR = False;  when  = True,
!         compare Mie code test case results with correct answers
!         and abort if even one result is inaccurate.

!         The test case is :  Mie size parameter = 10
!                             refractive index   = 1.5 - 0.1 i
!                             scattering angle = 140 degrees
!                             1 Sekera moment

!         Results for this case may be found among the test cases
!         at the end of reference (1).

!         *** NOTE *** When running on some computers, esp. in single
!         precision, the Accur criterion below may have to be relaxed.
!         However, if Accur must be set larger than 10**-3 for some
!         size parameters, your computer is probably not accurate
!         enough to do Mie computations for those size parameters.

!     Routines called :  ERRMSG, MIPRNT, TSTBAD

!     .. Scalar Arguments ..

LOGICAL     ANYANG, COMPAR, PERFCT
INTEGER     IPOLZN, MOMDIM, NMOM, NUMANG
REAL(EB)    GQSC, MIMCUT, QEXT, QSCA, XX
COMPLEX(EB) CREFIN, SBACK, SFORW

!     .. Array Arguments ..

REAL(EB)    PMOM( 0:MOMDIM, * ), XMU( * )
COMPLEX(EB) S1( * ), S2( * ), TBACK( * ), TFORW( * )

!     .. Local Scalars ..

LOGICAL     ANYSAV, OK, PERSAV
INTEGER     IPOSAV, M, N, NMOSAV, NUMSAV
REAL(EB)    MIMSAV, TESTGQ, TESTQE, TESTQS,XMUSAV, XXSAV
COMPLEX(EB) CRESAV, TESTS1, TESTS2, TESTSB, TESTSF

!     .. Local Arrays ..

LOGICAL     CALCMO( 4 ), PRNT( 2 )
REAL(EB)    TESTPM( 0:1 )
COMPLEX(EB) TESTTB( 2 ), TESTTF( 2 )


SAVE      XXSAV, CRESAV, MIMSAV, PERSAV, ANYSAV, NMOSAV, IPOSAV,NUMSAV, XMUSAV

DATA      TESTQE / 2.459791_EB /,&
          TESTQS / 1.235144_EB /,&
          TESTGQ / 1.139235_EB /,&
          TESTSF / ( 61.49476_EB, -3.177994_EB ) /,&
          TESTSB / ( 1.493434_EB,  0.2963657_EB ) /,&
          TESTS1 / ( -0.1548380_EB, -1.128972_EB ) /,&
          TESTS2 / ( 0.05669755_EB, 0.5425681_EB ) /,&
          TESTTF / ( 12.95238_EB, -136.6436_EB ),&
                   ( 48.54238_EB, 133.4656_EB ) /,&
          TESTTB / ( 41.88414_EB, -15.57833_EB ),&
                   ( 43.37758_EB, -15.28196_EB ) /,&
          TESTPM / 227.1975_EB, 183.6898_EB /

IF( .NOT.COMPAR ) THEN
!                                   ** Save certain user input values
   XXSAV  = XX
   CRESAV = CREFIN
   MIMSAV = MIMCUT
   PERSAV = PERFCT
   ANYSAV = ANYANG
   NMOSAV = NMOM
   IPOSAV = IPOLZN
   NUMSAV = NUMANG
   XMUSAV = XMU( 1 )
!                                   ** Reset input values for test case
   XX     = 10.0_EB
   CREFIN = ( 1.5_EB, -0.1_EB )
   MIMCUT = 0.0_EB
   PERFCT = .FALSE.
   ANYANG = .TRUE.
   NMOM   = 1
   IPOLZN = -1
   NUMANG = 1
   XMU( 1 ) = -0.7660444_EB

ELSE
!                                    ** Compare test case results with
!                                    ** correct answers and abort if bad
   OK = .TRUE.

   IF( WRONG( QEXT,TESTQE ) ) OK = TSTBAD( 'QEXT', ABS( ( QEXT - TESTQE ) / TESTQE ) )

   IF( WRONG( QSCA,TESTQS ) )OK = TSTBAD( 'QSCA', ABS( ( QSCA - TESTQS ) / TESTQS ) )

   IF( WRONG( GQSC,TESTGQ ) ) OK = TSTBAD( 'GQSC', ABS( ( GQSC - TESTGQ ) / TESTGQ ) )

   IF( WRONG( REAL( SFORW ),REAL( TESTSF ) ) .OR. WRONG( AIMAG( SFORW ),AIMAG( TESTSF ) ) ) &
         OK = TSTBAD( 'SFORW', ABS( ( SFORW - TESTSF ) / TESTSF ) )

   IF( WRONG( REAL( SBACK ),REAL( TESTSB ) ) .OR. WRONG( AIMAG( SBACK ),AIMAG( TESTSB ) ) ) &
         OK = TSTBAD( 'SBACK', ABS( ( SBACK - TESTSB ) / TESTSB ) )

   IF( WRONG( REAL( S1(1) ),REAL( TESTS1 ) ) .OR. WRONG( AIMAG( S1(1) ),AIMAG( TESTS1 ) ) ) &
         OK = TSTBAD( 'S1', ABS( ( S1(1) - TESTS1 ) / TESTS1 ) )

   IF( WRONG( REAL( S2(1) ),REAL( TESTS2 ) ) .OR. WRONG( AIMAG( S2(1) ),AIMAG( TESTS2 ) ) ) &
         OK = TSTBAD( 'S2', ABS( ( S2(1) - TESTS2 ) / TESTS2 ) )

   DO N = 1, 2
      IF( WRONG( REAL( TFORW(N) ),REAL( TESTTF(N) ) ) .OR. &
               WRONG( AIMAG( TFORW(N) ), &
               AIMAG( TESTTF(N) ) ) ) OK = TSTBAD( 'TFORW', &
               ABS( ( TFORW(N) - TESTTF(N) ) / TESTTF(N) ) ) 

      IF( WRONG( REAL( TBACK(N) ),REAL( TESTTB(N) ) ) .OR. &
               WRONG( AIMAG( TBACK(N) ), &
               AIMAG( TESTTB(N) ) ) ) OK = TSTBAD( 'TBACK', &
               ABS( ( TBACK(N) - TESTTB(N) ) / TESTTB(N) ) )
   END DO

   DO M = 0, 1
      IF ( WRONG( PMOM(M,1), TESTPM(M) ) ) OK =  TSTBAD( 'PMOM', ABS( (PMOM(M,1)-TESTPM(M)) / TESTPM(M) ) )
   END DO


   IF( .NOT.OK ) THEN
      PRNT( 1 ) = .TRUE.
      PRNT( 2 ) = .TRUE.
      CALCMO( 1 ) = .TRUE.
      CALCMO( 2 ) = .FALSE.
      CALCMO( 3 ) = .FALSE.
      CALCMO( 4 ) = .FALSE.
      CALL MIPRNT( PRNT, XX, PERFCT, CREFIN, NUMANG, XMU, QEXT, &
                        QSCA, GQSC, NMOM, IPOLZN, MOMDIM, CALCMO, PMOM, &
                        SFORW, SBACK, TFORW, TBACK, S1, S2 )
      CALL ERRMSG( 'MIEV0 -- Self-test failed',.TRUE.)
   END IF
!                                       ** Restore user input values
   XX     = XXSAV
   CREFIN = CRESAV
   MIMCUT = MIMSAV
   PERFCT = PERSAV
   ANYANG = ANYSAV
   NMOM   = NMOSAV
   IPOLZN = IPOSAV
   NUMANG = NUMSAV
   XMU( 1 ) = XMUSAV

END IF

CONTAINS

LOGICAL FUNCTION WRONG(CALC,EXACT)
REAL(EB) ACCUR
REAL(EB), INTENT(IN) :: CALC,EXACT
DATA  ACCUR / 1.E-4_EB /
WRONG = ABS( ( CALC - EXACT ) / EXACT )>ACCUR
END FUNCTION WRONG

END SUBROUTINE TESTMI


SUBROUTINE ErrMsg( MESSAG, FATAL )
! Print out a warning or error message;  abort if error after making symbolic dump (machine-specific)
USE COMP_FUNCTIONS, ONLY: SHUTDOWN
CHARACTER :: MESSAG*(*)
LOGICAL ::  FATAL,MSGLIM
INTEGER ::  MAXMSG, NUMMSG
SAVE      MAXMSG, NUMMSG, MSGLIM
DATA      NUMMSG / 0 /,  MAXMSG / 100 /,  MSGLIM /.FALSE./
IF (FATAL) CALL SHUTDOWN(MESSAG)
NUMMSG = NUMMSG + 1
IF( MSGLIM ) RETURN
IF( NUMMSG<=MAXMSG ) THEN
   WRITE(LU_ERR, '(/,2A,/)' ) ' ****** WARNING *****  ', MESSAG
ELSE
   WRITE(LU_ERR, 9000 )
   MSGLIM = .True.
END IF
RETURN
9000 FORMAT( / , / , ' ****** TOO MANY WARNING MESSAGES --  ','They will no longer be printed *******', / , / )
END SUBROUTINE ErrMsg


LOGICAL FUNCTION WrtBad( VarNam )
!          Write names of erroneous variables and return 'TRUE'
!      INPUT :   VarNam = Name of erroneous variable to be written
!                         ( CHARACTER, any length )
CHARACTER VarNam*(*)
INTEGER   MAXMSG, NUMMSG
SAVE      NUMMSG, MAXMSG
DATA      NUMMSG / 0 /, MAXMSG / 50 /
WrtBad = .TRUE.
NUMMSG = NUMMSG + 1
WRITE(LU_ERR, '(3A)' ) ' ****  Input variable  ', VarNam,'  in error  ****'
IF( NUMMSG==MAXMSG ) CALL ErrMsg( 'Too many input errors.  Aborting...',.TRUE.)
RETURN
END FUNCTION WrtBad


LOGICAL FUNCTION WrtDim( DimNam, Minval )

! Write name of too-small symbolic dimension and the value it should be increased to;  return 'TRUE'

!      INPUT :  DimNam = Name of symbolic dimension which is too small
!                        ( CHARACTER, any length )
!               Minval = Value to which that dimension should be
!                        increased (at least)
CHARACTER :: DimNam*(*)
INTEGER ::  Minval

WRITE(LU_ERR, '(3A,I7)' ) ' ****  Symbolic dimension  ', DimNam,'  should be increased to at least ', Minval
WrtDim = .TRUE.
RETURN
END FUNCTION WrtDim


LOGICAL FUNCTION TstBad( VarNam, RelErr )
! Write name (VarNam) of variable failing self-test and its percent error from the correct value; return 'FALSE'
CHARACTER :: VarNam*(*)
REAL(EB) ::  RelErr
TstBad = .FALSE.
WRITE(LU_ERR, '(/,3A,1P,E11.2,A)' ) ' *** Output variable ', VarNam, ' differed by ', &
     100.*RelErr,' per cent from correct value.  Self-test failed.'
RETURN
END FUNCTION TstBad
 

REAL(EB) FUNCTION SQ(CTMP)
COMPLEX(EB), INTENT(IN) :: CTMP
SQ = REAL( CTMP )**2 + AIMAG( CTMP )**2
END FUNCTION SQ
 
END MODULE MIEV
