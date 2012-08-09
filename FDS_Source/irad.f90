!------------------------------------------------------------------------------
! This file contains 4 modules
! 1 - RADCONS
! 2 - RADCALV
! 3 - SPECDATA
! 4 - MIEV
!------------------------------------------------------------------------------


MODULE RADCONS
!------------------------------------------------------------------------------
! Radiation Parameters

USE PRECISION_PARAMETERS
IMPLICIT NONE

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: DLN
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: BBFRAC, WL_LOW, WL_HIGH
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: DLX, DLY, DLZ, DLB, RSA

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: DLM
INTEGER, ALLOCATABLE, DIMENSION(:)   :: NRP

REAL(EB) :: RADTMP, PATH_LENGTH, RADIATIVE_FRACTION
REAL(EB) :: DGROUP_A, DGROUP_B, WEIGH_CYL
REAL(EB) :: DPHI0, FOUR_SIGMA, RPI_SIGMA, LTSTEP, RTMPMAX, RTMPMIN

INTEGER :: TIME_STEP_INCREMENT,NMIEANG
INTEGER :: NRDMIE, NLMBDMIE, NDG = 50
INTEGER :: NRT,NCO,UIIDIM,NLAMBDAT,NKAPPAT,NKAPPAZ

LOGICAL :: WIDE_BAND_MODEL, CH4_BANDS

INTEGER :: N_RADCAL_SPECIES,RADCAL_SPECIES_INDEX(11),N_KAPPA_T=44,N_KAPPA_Y=50
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: Z2RADCAL_SPECIES
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: RADCAL_SPECIES2KAPPA
CHARACTER(30) :: RADCAL_SPECIES(11)='null'


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
!     NDG       Number of PARTICLE radii in WQABS and WQSCA arrays
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



MODULE RADCALV

!------------------------------------------------------------------------------
! VARIABLES
! OMMIN : (REAL) MINIMUM WAVE NUMBER IN SPECTRUM, CM-1
! OMMAX : (REAL) MAXIMUM WAVE NUMBER IN SPECTRUM, CM-1
! NOM   : (INTEGER) NUMBER OF WAVELENGTH INTERVALS
! GC    : (REAL) COLLISION BROADENED HALF-WIDTH
! AMBDA : (REAL) WAVELENGTH, IN MICROMETER 


!------------------------------------------------------------------------------
! Module wrapper for RadCal subroutine

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS, ONLY: AL2O3,NEW_ABSORPTION
IMPLICIT NONE

PRIVATE
CHARACTER(255), PARAMETER :: iradid='$Id$'
CHARACTER(255), PARAMETER :: iradrev='$Revision$'
CHARACTER(255), PARAMETER :: iraddate='$Date$'

PUBLIC OMMAX, OMMIN, DD, SPECIE, SVF, PLANCK, P, RCT, RCALLOC,                 &
       INIT_RADCAL, RADCAL, RCDEALLOC, GET_REV_irad

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: GAMMA, SD15, SD, SD7, SD3,             &
                                         SD1_CH4,     SD2_CH4,     OM_BND_CH4,  & ! FOR METHANE NEW DATA
                                         SD1_C3H8,    SD2_C3H8,    OM_BND_C3H8, & ! FOR PROPANE NEW DATA
                                         SD1_C7H16,   SD2_C7H16,   OM_BND_C7H16,& ! FOR HEPTANE NEW DATA
                                         SD1_CH3OH,   SD2_CH3OH,   OM_BND_CH3OH,& ! FOR METHANOL NEW DATA
                                         SD3_CH3OH,                             & 
                                         SD1_C7H8,    SD2_C7H8,    OM_BND_C7H8, & ! FOR TOLUENE
                                         SD3_C7H8,    SD4_C7H8,                 &
                                         SD1_C3H6,    SD2_C3H6,    OM_BND_C3H6, & ! FOR PROPYLENE NEW DATA
                                         SD3_C3H6,                              & 
                                         SD1_MMA,     SD2_MMA,     OM_BND_MMA      ! FOR MMA NEW DATA      


REAL(EB), ALLOCATABLE, DIMENSION(:)   :: SPECIE, QW, TTAU, XT, AB,             &
                                         AMBDA, ATOT, BCNT, P, UUU, GC, X

REAL(EB), ALLOCATABLE, DIMENSION(:)   :: SD_CH4_TEMP,   SD_C3H8_TEMP, SD_C7H16_TEMP, &
                                         SD_CH3OH_TEMP, SD_C7H8_TEMP, SD_C3H6_TEMP,  &
                                         SD_MMA_TEMP

REAL(EB) :: OMMIN, OMMAX, TWALL, RCT, AC, AD, DD, XPART, TAU,                  &
            SVF, TAUS, XTOT, XSTAR

INTEGER :: N_TEMP_CH4,   N_BAND_CH4,   &
           N_TEMP_C3H8,  N_BAND_C3H8,  &
           N_TEMP_C7H16, N_BAND_C7H16, & 
           N_TEMP_CH3OH, N_BAND_CH3OH, &
           N_TEMP_C7H8,  N_BAND_C7H8,  &
           N_TEMP_C3H6,  N_BAND_C3H6,  &
           N_TEMP_MMA,   N_BAND_MMA

INTEGER :: NOM

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


!==============================================================================
SUBROUTINE INIT_RADCAL
!==============================================================================
! 
 USE RADCONS, ONLY: RADTMP
!------------------------------------------------------------------------------

 TWALL = RADTMP

 IF(OMMAX<1100._EB) THEN 
    NOM=INT((OMMAX-OMMIN)/5._EB)
 ELSEIF(OMMIN>5000._EB) THEN
    NOM=INT((OMMAX-OMMIN)/50._EB)
 ELSEIF(OMMIN<1100._EB.AND.OMMAX>5000._EB) THEN
    NOM=INT((1100._EB-OMMIN)/5._EB)+INT((5000._EB-1100._EB)/25._EB) +INT((OMMAX-5000._EB)/50._EB)
 ELSEIF(OMMIN<1100._EB) THEN
    NOM=INT((1100._EB-OMMIN)/5._EB)+INT((OMMAX-1100._EB)/25._EB)
 ELSEIF(OMMAX>5000._EB) THEN
    NOM=INT((5000._EB-OMMIN)/25._EB)+INT((OMMAX-5000._EB)/50._EB)
 ELSE
    NOM=INT((OMMAX-OMMIN)/25._EB)
 ENDIF
!------------------------------------------------------------------------------
END SUBROUTINE INIT_RADCAL


!==============================================================================
SUBROUTINE RADCAL(AMEAN,AP0,RADCAL_ID)
!==============================================================================
! OUTPUT
! AMEAN : (REAL) EFFECTIVE ABSORPTION COEFFICIENT IN UNIFORM MEDIUM
! AP0   : (REAL) PLANCK MEAN ABSORPTION COEFFCIENT IN UNIFORM MEDIUM
!
! LOCAL
! DOM   : (REAL) INCREMENT OF WAVENUMBER
! OMEGA : (REAL) WAVENUMBER
! UUU   : (REAL) OPTICAL THICKNESS OF THE ith SPECIES,  UNITS: cmSTP
!
! RECALL: CO2: 5 BANDS (4 Modeled, 1 Tabulated)
!         H2O: 5 BANDS (5 Tabulated)
!         CO: 1 BANDS (Modeled)
!        old CH4: 3 BANDS (1 Modeled, 2 Tabulated)
!------------------------------------------------------------------------------
USE RADCONS, ONLY:RPI_SIGMA

! VARIABLES PASSED OUTPUT
REAL (EB), INTENT(OUT) :: AMEAN, AP0

! LOCAL VARIABLES

REAL(EB) :: DOM, ABGAS, PTOT, TEMP, UK, XD, YD, XX, ENN, ARG, ARGNEW, RSL, RSS,  &
            ABLONG, ABSHRT, ABIL, ABIS, OMEGA, WL, DAMBDA, SDWEAK, GDINV,        &
            GDDINV, YC, Y, AIWALL, XC, AOM, Q, LTERM, AZORCT, RCT4

INTEGER  :: I, II, KK, NM, N, MM, KMAX, KMIN

CHARACTER(30),INTENT(IN):: RADCAL_ID
! [NOTE: THE TOTAL INTENSITY CALCULATED IS THAT WHICH LEAVES INTERVAL J=1.
! P(I,J) IS PARTIAL PRESSURE, ATM, OF SPECIES I IN  INTERVAL J.
! I=1,2,3,4,5, OR 6 IMPLIES SPECIES IS CO2, H2O, CH4, CO, O2, OR N2, RESP.]

DOM    = 5.0_EB
OMEGA  = OMMIN-DOM
NM     = NOM-1
AZORCT = 273._EB/RCT
RCT4   = RCT**4

! COMPUTE TOTAL PRESSURE
PTOT   = SUM(P(1:6))

! LOOP 1000 COMPUTES EACH SPECTRAL CONTRIBUTION, LOOP OVER EACH WAVE NUMBER INTERVALS

L1000: DO KK=1, NOM

   OMEGA=OMEGA+DOM
   IF(OMEGA>1100._EB) OMEGA = OMEGA+20._EB
   IF(OMEGA>5000._EB) OMEGA = OMEGA+25._EB

   AMBDA(KK) = 10000._EB/OMEGA  ! COMPUTE LAMBDA (WAVELENGTH) IN microMETER
   ABGAS     = 0._EB
 
! LOOP 200 COMPUTES THE CONTRIBUTION OF EACH SPECIES TO TAU
! IF SPECIE(I) IS SET TO 0., THAT PARTICULAR RADIATING SPECIES IS NOT PRESENT.  THE SPECIES CONSIDERED ARE
!          I   SPECIES
!          1     CO2
!          2     H2O
!          3     FUEL (CH4, C3H8, C7H16, CH3OH, ...)
!          4     CO
!          5     PARTICULATES

   L200: DO I=1,4

      IF(SPECIE(I) <= ZERO_P) CYCLE L200
 
! LOOP 100 IS FOR EACH ELEMENT ALONG PATH
! (CALCULATION PROCEEDS IN ACCORDANCE WITH THE SLG MODEL, TABLE 5-18, IN NASA SP-3080._EB)
! MODEL BASED ON RANDOM BAND MODEL

      IF(KK<=1) THEN
         UUU(I) = AZORCT*P(I)*100.*DD ! COMPUTE OPTICAL THICKNESS FOR THE ith SPECIES

!----------------------------------------------------------------------------------------
! THIS CAN BE TRANSFORMED AS A FUNCTION
! Compute collisional broadening half-width at half height, EQ 5-34 NASA REPORT SP-3080
         GC(I) = 0._EB

! INCLUDE NON-RESONANT FOREIGN AND SELF-BROADENING COLLISIONS
         DO II=1,6
            GC(I)=GC(I)+GAMMA(I,II)*P(II)*SQRT(AZORCT)
         ENDDO

! INCLUDE RESONANT SELF-BROADENING COLLISIONS
         GC(I)=GC(I)+GAMMA(I,7)*P(I)*AZORCT     
!----------------------------------------------------------------------------------------
      ENDIF

      IF(P(I)<=ZERO_P) THEN
         XSTAR = 1.E-34_EB
         AC    = 1._EB
         AD    = 1._EB
         X(I)  = XSTAR
      ELSE
         TEMP  = RCT
         SELECT CASE (I)
            CASE (1)
               CALL CO2(OMEGA,TEMP,GC(1),SDWEAK,GDINV,GDDINV)
            CASE (2)
               CALL H2O(OMEGA,TEMP,GC(2),SDWEAK,GDINV,GDDINV)
            CASE (3)
               IF (NEW_ABSORPTION) THEN
                  SELECT CASE(RADCAL_ID)
                     CASE('METHANE')
                        CALL CH4_NEW(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                     CASE('PROPANE')
                        CALL C3H8(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                     CASE('N-HEPTANE')
                        CALL C7H16(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                     CASE('METHANOL')
                        CALL CH3OH(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                     CASE('TOLUENE')
                        CALL C7H8(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                     CASE('PROPYLENE')
                        CALL C3H6(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                     CASE('MMA')
                        CALL MMA(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                     CASE DEFAULT
                        CALL CH4(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
                  END SELECT                 
               ELSE
                  CALL CH4(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
               ENDIF
            CASE (4)
               CALL CO(OMEGA,TEMP,GC(4),SDWEAK,GDINV,GDDINV)
         END SELECT

         UK    = SDWEAK*UUU(I)
         XSTAR = UK+1.E-34_EB
         ABGAS = UK/DD+ABGAS
! RECALL: AD = Doppler broadened fine structure parameter: (line_half_width/line_spacing)
! RECALL: AC = Collision broadened fine structure parameter: (line_half_width/line_spacing)
         AD    = GDDINV
         AC    = GDINV

         IF(XSTAR>=1.E-6_EB) THEN

            XD = 1.7_EB*AD*SQRT(DLOG(1._EB+(XSTAR/(1.7_EB*AD))**2)) ! SP-3080 EQ. 3-53 Doppler broadened
                                                                    ! optical depth
            YD = 1._EB-(XD/XSTAR)**2
            XC = XSTAR/SQRT(1._EB+0.25_EB*XSTAR/AC)   ! SP-3080 eq. 3-52 Random band models
                                                      ! optical depth
 
! THE FOLLOWING LOOP COMPUTES THE OPTICAL THICKNESS, XC, FOR METHANE USING 
! THE GODSON EQUATION AND AN APPROXIMATION TO THE LADENBERG-REICHE
! FUNCTION AS RECOMMENDED BY BROSMER AND TIEN (JQSRT 33,P 521), EQ. 2, PAGE 523 
! ELSASSER MODEL
! THE ERROR FUNCTION IS FOUND FROM ITS SERIES EXPANSION.
! 
            IF((I==3).AND.(XC<=10).AND.(RADCAL_ID=='METHANE')) THEN
               AOM=XC
               XX = .5_EB*SQRTPI*XC

               IF(XX<=3._EB) THEN
                  ENN=1._EB

                  DO N=1,30
                     ENN = ENN*REAL(N,EB)  ! ENN = N!
                     MM  = 2*N+1
                     ARG = 1.128379_EB*(-1._EB)**N*((.88622693_EB*XC)**MM)/(REAL(MM,EB)*ENN)
! THIS CORRESPONDS TO ERR(SQRT(PI)/2*XC) = ERR(XX)
                     ARGNEW=ARG+AOM
!     IF(ABS(ARG/ARGNEW)<.000001)N=30
                     AOM=ARGNEW
                  ENDDO

               ELSE
                  AOM=1._EB-EXP(-XX**2)/(SQRTPI*XX)
               ENDIF

               IF (AOM>=1._EB) AOM=.9999999_EB
               XC=-DLOG(1._EB-AOM)
            ENDIF

            YC=1._EB-(XC/XSTAR)**2
            Y=MAX(1._EB/YC**2+1._EB/YD**2-1._EB,1._EB)
            X(I)=XSTAR*(SQRT(1._EB-(Y**(-.5_EB))))
         ELSE
            X(I)=XSTAR
         ENDIF

      ENDIF

   END DO L200

! DETERMINE OPTICAL DEPTH OF SOOT

   IF (SPECIE(5)<=ZERO_P) THEN
      XPART=0._EB
   ELSE
      CALL POD(OMEGA)
   ENDIF
   AB(KK)=ABGAS+XPART/DD

! Evaluate the combined spectral transmettance and radiance

   XTOT = 0._EB
   DO I=1,4
      IF(SPECIE(I)<=ZERO_P) X(I)=0._EB
      XTOT=X(I)+XTOT
   ENDDO

   XTOT=XTOT+XPART

   IF(XTOT>=99._EB) THEN
      TAU=0._EB
   ELSE
      TAU=EXP(-XTOT)
   ENDIF

   XT(KK)  = XTOT
   TTAU(KK)= TAU
   QW(KK)  = -(TAU-1._EB)*PLANCK(RCT,AMBDA(KK))+TTAU(KK)*PLANCK(TWALL,AMBDA(KK))

ENDDO L1000
     
! INTEGRATE THE RADIANCE OVER THE SPECTRUM

Q=QW(1)*(AMBDA(1)-AMBDA(2))

DO KK=2,NM
   Q=Q+QW(KK)*(AMBDA(KK-1)-AMBDA(KK+1))/2._EB
ENDDO

Q=Q+QW(NOM)*(AMBDA(NOM-1)-AMBDA(NOM))

! DETERMINE SOOT RADIANCE FOR SHORT AND LONG WAVELENGTHS.
     
RSL=0._EB
RSS=0._EB
ABLONG=0._EB
ABSHRT=0._EB
ABIL=0._EB
ABIS=0._EB

IF(.NOT.(SPECIE(5)<=ZERO_P .AND. TWALL<=ZERO_P)) THEN
   KMAX=INT(OMMIN)
   DO KK=5,KMAX,5
      OMEGA=FLOAT(KK)
      WL=10000._EB/OMEGA
      DAMBDA=10000._EB/(OMEGA-2.5_EB)-10000._EB/(OMEGA+2.5_EB)
      CALL POD(OMEGA)
   
      IF(XPART>=33._EB) THEN
         TAUS=0._EB
      ELSE
         TAUS=EXP(-XPART)
      ENDIF
      RSL=RSL-(TAUS-1._EB)*PLANCK(RCT,WL)*DAMBDA
      ABLONG=ABLONG+XPART/DD*PLANCK(RCT,WL)*DAMBDA/RPI_SIGMA/RCT4
      ABIL=ABIL+XPART/DD*PLANCK(TWALL,WL)*DAMBDA/RPI_SIGMA/(TWALL+.000001_EB)**4
      RSL=RSL+TAUS*PLANCK(TWALL,WL)*DAMBDA
   ENDDO
   KMIN=INT(OMMAX)
   DO KK=KMIN,25000,100
      OMEGA=FLOAT(KK)
      WL=10000._EB/OMEGA
      DAMBDA=10000._EB/(OMEGA-50._EB)-10000._EB/(OMEGA+50._EB)
      CALL POD(OMEGA)
      IF(XPART>=33._EB) THEN
         TAUS=0._EB
      ELSE
         TAUS=EXP(-XPART)
      ENDIF
      RSS=RSS-(TAUS-1._EB)*PLANCK(RCT,WL)*DAMBDA
      ABSHRT=ABSHRT+XPART/DD*PLANCK(RCT,WL)*DAMBDA/RPI_SIGMA/RCT4
      ABIS=ABIS+XPART/DD*PLANCK(TWALL,WL)*DAMBDA/RPI_SIGMA/(TWALL+.000001_EB)**4
      RSS=RSS+TAUS*PLANCK(TWALL,WL)*DAMBDA
   ENDDO
ENDIF

Q=Q+RSS+RSL

! THE FOLLOWING SECTION COMPUTES THE MEAN ABSORPTION COEFFICIENTS IF THE SYSTEM IS HOMOGENEOUS (IE., NPT=1).

NM=NOM-1
AIWALL=AB(1)*(AMBDA(1)-AMBDA(2))/2._EB*PLANCK(TWALL,AMBDA(1))
AP0=AB(1)*(AMBDA(1)-AMBDA(2))/2._EB*PLANCK(RCT,AMBDA(1))

DO KK=2,NM
   AIWALL=AIWALL+AB(KK)*(AMBDA(KK-1)-AMBDA(KK+1))/2._EB *PLANCK(TWALL,AMBDA(KK))
   AP0=AP0+AB(KK)*(AMBDA(KK-1)-AMBDA(KK+1))/2._EB*PLANCK(RCT,AMBDA(KK))
ENDDO

AP0=(AP0+AB(NOM)*(AMBDA(NM)-AMBDA(NOM))/2._EB *PLANCK(RCT,AMBDA(NOM)))/RPI_SIGMA/RCT4

IF(ABS(TWALL-RCT)<=SPACING(TWALL) .OR. TWALL<=ZERO_P) THEN
   AIWALL=AP0
   LTERM = MAX(1E-20_EB,(Q/RPI_SIGMA-RCT4)/(-RCT4))
   AMEAN=-1._EB/DD*DLOG(LTERM)
ELSE
   AIWALL=(AIWALL+AB(NOM)*(AMBDA(NM)-AMBDA(NOM))/2._EB*PLANCK(TWALL,AMBDA(NOM)))/RPI_SIGMA/TWALL**4
   LTERM = MAX(1E-20_EB,(Q/RPI_SIGMA-RCT4)/(TWALL**4-RCT4))   
   AMEAN=-1._EB/DD*DLOG(LTERM)
ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE RADCAL


!==============================================================================
SUBROUTINE CO2(OMEGA,TEMP,GC1,SDWEAK,GDINV,GDDINV)
!==============================================================================

INTEGER I,J,K,L
REAL(EB) OMEGA,TEMP,GC1,SDWEAK,GDINV,GDDINV,AA,BB,CC,QQ,EE,FF,GG, &
         SMINUS,SPLUS,SDSTRG,GD,OM1,OM2,OM3,T0,Q2,BE,COM1, &
         COM2,COM3,X13,X23,X33,XBAR,OM12,ALPHA,OMPRIM,V3,GAM, &
         OMVV3,DELTA,V,OMVBAR,F1,F2,UNFLO1,UNFLO2,UNFLO3,TEST, &
         VBAR1,OMA,OMB,TTEMP,TT,T1,TW,WW,TEMP1,TEMP2,TEMP3,WM, &
         DINV,A,B,D,G,W1,DINV1,DINV2,DINV3,Q2OT,T0OT,Q2OT0

IF(OMEGA>5725._EB) THEN
   SDWEAK = 0._EB
   GDINV  = 1._EB
   GDDINV = 1._EB
ELSE
   WM = 44._EB
! COMPUTE DOPPLER BROADENING HALF-WIDTHS SP-3080 EQ: 5-35
   GD = 5.94E-6_EB*OMEGA*SQRT(TEMP/(273._EB*WM))

   IF(OMEGA>4550._EB) THEN
! CONTRIBUTION TO 2.0 MICRON BAND FROM:
! (000)-(041)
! (000)-(121),
! (000)-(201) TRANSITIONS.

      OM1 = 1354.91_EB ! (100)
      OM2 = 673.0_EB   ! (010)
      OM3 = 2396.49_EB ! (001)

! BAND CENTER 

      BCNT(1) = 4860.5_EB ! (000)-(041)
      BCNT(2) = 4983.5_EB ! (000)-(121)
      BCNT(3) = 5109.0_EB ! (000)-(201)

      T0   = 300._EB
      Q2   = 1.4388_EB
      T0OT = T0/TEMP
      Q2OT = -Q2/TEMP

      BE   = 0.391635_EB ! ROTATIONAL CONSTANT

! COMPUTE THE ASSOCIATED TRANSITION WAVENUMBERS

      COM1 = 4._EB*OM2+OM3      ! (000)-(041)
      COM2 = OM1+2._EB*OM2+OM3  ! (000)-(121)
      COM3 = 2._EB*OM1+OM3      ! (000)-(201)


! COMPUTE THE INTEGRATED INTENSITY OF THE BAND
! 0.272: INTEGRATED INTENSITY (IN CM-2-ATM-1) OF BAND BCNT(1) = 4860.5_EB ((041)-(000)) AT T=300K
! 1.01 : INTEGRATED INTENSITY (IN CM-2-ATM-1) OF BAND BCNT(2) = 4983.5_EB ((121)-(000)) AT T=300K
! 0.426: INTEGRATED INTENSITY (IN CM-2-ATM-1) OF BAND BCNT(3) = 5109  _EB ((201)-(000)) AT T=300K
! VALUES FROM PENNER AND VARANASI JQSRT VOL 4 pp 799-806


      ATOT(1) = 0.272_EB*T0OT*(1._EB-EXP(Q2OT*COM1))/(1._EB-EXP(Q2OT*OM2))**4/(1._EB-EXP(Q2OT*OM3))
!  NOTE: [(1._EB-EXP(Q2OT*OM2))**4*(1._EB-EXP(Q2OT*OM3))] IS THE ROVIBRATIONAL PARTITION FUNCTION
!  EQ 2 IN PENNER AND VARANASI JQSRT VOL 4 pp 799-806

      ATOT(2) = 1.01_EB*T0OT*(1._EB-EXP(Q2OT*COM2))/(1._EB-EXP(Q2OT*OM1))/(1._EB-EXP(Q2OT*OM2))**2/ &
              (1._EB-EXP(Q2OT*OM3))

      ATOT(3) = 0.426_EB*T0OT*(1._EB-EXP(Q2OT*COM3))/(1._EB-EXP(Q2OT*OM1))**2/(1._EB-EXP(Q2OT*OM3))

      SDWEAK  = 0.0_EB

! COMPUTE THE LINE STRENGTH USING EQ 11-44 FROM PENNER IN QUANTITATIVE MOLECULAR SPECTROMETRY, 1955
! SDWEAK = SUM(ATOT*PROPABILITY) 
! JUST OVERLAPPING SPECTRAL LINE
      DO K=1,3
         SDWEAK = SDWEAK+ATOT(K)*(-Q2OT)/(4._EB*BE)*ABS(OMEGA-BCNT(K))*EXP(Q2OT/(4._EB*BE)*(OMEGA-BCNT(K))**2)
      ENDDO

      DINV   = 1._EB/(4._EB*BE)
      GDINV  = GC1*DINV
      GDDINV = GD*DINV
!***EXPRESS S/D AT STP, AS IS IN NASA SP-3080
      SDWEAK=SDWEAK*TEMP/273._EB
   ELSEIF((OMEGA<=4550._EB).AND.(OMEGA>3800._EB)) THEN
      SDWEAK=0._EB
      GDINV=1._EB
      GDDINV=1._EB
   ELSEIF((OMEGA<=3800._EB).AND.(OMEGA>3050._EB)) THEN
      B=.391635_EB
      A=.0030875_EB
      X13=-19.37_EB
      X23=-12.53_EB
      X33=-12.63_EB
      OM1=1354.91_EB
      OM2=673._EB
      OM3=2396.49_EB
      T0=300._EB
      T0OT = T0/TEMP
      Q2=1.4388_EB
      Q2OT = -Q2/TEMP
      Q2OT0 = -Q2/T0
      XBAR=.5_EB*(.5_EB*X13+X23)
      OM12=.5_EB*(.5_EB*OM1+OM2)
      SDWEAK=0._EB
      SDSTRG=0._EB
      IF(OMEGA<=2395._EB) THEN
         ALPHA=2700._EB
         OMPRIM=OM3
         AA=ALPHA*B*Q2/(A*(1._EB-EXP(OM3*Q2OT0))*(1._EB-EXP(OM12*Q2OT0))**3*(1._EB+EXP(OM12*Q2OT0))*(1._EB-EXP(OMPRIM*Q2OT0)))
         BB=(1._EB-EXP(Q2OT*OMEGA))*(1._EB-EXP(Q2OT*OM3))*(1._EB-EXP(OM12*Q2OT))**3*(1._EB+EXP(OM12*Q2OT)) &
           *(1._EB-EXP(Q2OT*OMPRIM))
         CC=AA*BB*OMEGA*T0/TEMP**2
         L202: DO J=1,20
            V=FLOAT(J-1)
            IF(J/2*2==J)G=(V+1._EB)*(V+3._EB)/4._EB
            IF(J/2*2/=J)G=(V+2._EB)*(V+2._EB)/4._EB
            L201: DO K=1,10
               V3=FLOAT(K-1)
               QQ=(V3+1._EB)*G*EXP(-(V3*OM3+V*OM12)*Q2OT)
               GAM=B-A*(V3+1._EB)
               OMVV3=OM3+.5_EB*X13+X23+2._EB*X33+XBAR*V+2._EB*X33*V3
               DELTA=A*(OMEGA-OMVV3)
               IF(GAM*GAM<=DELTA) CYCLE L202
               D=2._EB*SQRT(GAM*GAM-DELTA)
               OMVBAR=OMVV3*(1._EB-EXP(-OMVV3*Q2OT))
               F1=GAM-D/2
               F2=GAM+D/2._EB
               EE=Q2*GAM/(A*A*TEMP)
               UNFLO1=EE*DELTA*(1._EB+.5_EB*A/GAM)
               IF(UNFLO1<=-78.) CYCLE L202
               UNFLO2=EE*2._EB*GAM*F1
               IF(UNFLO2>=78.) CYCLE L202
               FF=EXP(EE*DELTA*(1._EB+.5_EB*A/GAM))
               SMINUS=CC*QQ/OMVBAR*ABS(F1)*FF*EXP(-EE*2._EB*GAM*F1)
               UNFLO3=EE*2._EB*GAM*F2
               IF(UNFLO3>=78.) THEN
                  SPLUS=0.
               ELSE
                  SPLUS=CC*QQ/OMVBAR*ABS(F2)*FF*EXP(-EE*2._EB*GAM*F2)
               ENDIF
               GG=SDWEAK
               SDWEAK=(SMINUS+SPLUS)/D+SDWEAK
               TEST=(SDWEAK-GG)/SDWEAK
               IF(TEST<.0001) CYCLE L202
               SDSTRG=SQRT(.5_EB*G)*(SQRT(SMINUS)+SQRT(SPLUS))/D+SDSTRG
            ENDDO L201
         ENDDO L202
         IF(SDWEAK<=ZERO_P) THEN
            SDWEAK=0._EB
            GDINV=1._EB
            GDDINV=1._EB
         ELSE
            DINV=SDSTRG*SDSTRG/SDWEAK
            GDINV=GC1*DINV
            GDDINV=GD*DINV
!***  EXPRESS S/D AT STP, AS IS K IN NASA SP-3080
            SDWEAK=SDWEAK*TEMP/273.
         ENDIF
      ELSE
!CALCULATE ABSORPTION COEF. AND LINE SPACING PARAMETER FOR 2.7 MICRON BAND
         L=1
!CONTRIBUTION TO 2.7 MICRON BAND FROM (000)-(021) AND (010)-(031) TRANS.
         ALPHA=28.5_EB
         OMPRIM=2._EB*OM2+OM3
         L120: DO
         AA=ALPHA*B*Q2/(A*(1._EB-EXP(OM3*Q2OT0))*(1._EB-EXP(OM12*Q2OT0))**3*(1._EB+EXP(OM12*Q2OT0))*(1._EB-EXP(OMPRIM*Q2OT0)))
         BB=(1._EB-EXP(Q2OT*OMEGA))*(1._EB-EXP(Q2OT*OM3))* (1._EB-EXP(OM12*Q2OT))**3*(1._EB+EXP(OM12*Q2OT)) &
           *(1._EB-EXP(Q2OT*OMPRIM))
         CC=AA*BB*OMEGA*T0/TEMP**2
         L102: DO J=1,20
         V=FLOAT(J-1)
         IF(J/2*2==J)G=(V+1._EB)*(V+3._EB)/4._EB
         IF(J/2*2/=J)G=(V+2._EB)*(V+2._EB)/4._EB
         VBAR1=-1._EB+(V+3._EB)*(V+4._EB)/(V+2.)/6._EB
         IF(J/2*2==J)VBAR1=-1._EB+(V+5._EB)/6._EB
         L101: DO K=1,10
         V3=FLOAT(K-1)
         QQ=(V3+1)*G*EXP((V3*OM3+V*OM12)*Q2OT)*(VBAR1+1._EB)
         GAM=B-A*(V3+1._EB)
         IF(L==2) THEN
            OMVV3=3728._EB-5._EB*V-47._EB*V3
            IF(V<=ZERO_P) OMVV3=3715._EB-47._EB*V3
         ELSE
            OMVV3=3598._EB-18._EB*V-47._EB*V3
            IF(V<=ZERO_P) OMVV3=3613._EB-47._EB*V3
         ENDIF
         DELTA=A*(OMEGA-OMVV3)
         IF(GAM*GAM<=DELTA) CYCLE L102
         D=2._EB*SQRT(GAM*GAM-DELTA)
         OMVBAR=OMVV3*(1._EB-EXP(OMVV3*Q2OT))
         F1=GAM-D/2._EB
         F2=GAM+D/2._EB
         EE=Q2*GAM/(A*A*TEMP)
         UNFLO1=EE*DELTA*(1._EB+.5_EB*A/GAM)
         IF(UNFLO1<=-78._EB) CYCLE L102
         UNFLO2=EE*2._EB*GAM*F1
         IF(UNFLO2>=78._EB) CYCLE L102
         FF=EXP(EE*DELTA*(1._EB+.5_EB*A/GAM))
         SMINUS=CC*QQ/OMVBAR*ABS(F1)*FF*EXP(-EE*2._EB*GAM*F1)
         UNFLO3=EE*2._EB*GAM*F2
         IF(UNFLO3>=78._EB) THEN
            SPLUS=0._EB
         ELSE
            SPLUS=CC*QQ/OMVBAR*ABS(F2)*FF*EXP(-EE*2._EB*GAM*F2)
         ENDIF
         GG=SDWEAK
         SDWEAK=(SMINUS+SPLUS)/D+SDWEAK
         TEST=(SDWEAK-GG)/SDWEAK
         IF(TEST<.0001_EB) CYCLE L102
         SDSTRG=SQRT(.5_EB*G)*(SQRT(SMINUS)+SQRT(SPLUS))/D+SDSTRG
         ENDDO L101
         ENDDO L102
         IF(L==2) EXIT L120
!CONTRIBUTION TO 2.7 MICRON BAND FROM (000)-(101) AND (010)-(111) TRANS.
         ALPHA=42.3_EB
         OMPRIM=OM1+OM3
         L=2
      ENDDO L120
!CALCULATE ABSORPTION COEF AND LINE SPACING PARAMETER FOR 4.3 MICRON BAND
      IF(SDWEAK<=ZERO_P) THEN
         SDWEAK=0._EB
         GDINV=1._EB
         GDDINV=1._EB
      ELSE
         DINV=SDSTRG*SDSTRG/SDWEAK
         GDINV=GC1*DINV
         GDDINV=GD*DINV
!***EXPRESS S/D AT STP, AS IS K IN NASA SP-3080
         SDWEAK=SDWEAK*TEMP/273.
      ENDIF
   ENDIF
ELSEIF((OMEGA<=3050._EB).AND.(OMEGA>2474.)) THEN
   SDWEAK=0._EB
   GDINV=1._EB
   GDDINV=1._EB
ELSEIF((OMEGA<=2474.).AND.(OMEGA>1975.)) THEN
   B=.391635_EB
   A=.0030875_EB
   X13=-19.37_EB
   X23=-12.53_EB
   X33=-12.63_EB
   OM1=1354.91_EB
   OM2=673._EB
   OM3=2396.49_EB
   T0=300._EB
   Q2=1.4388_EB
   XBAR=.5_EB*(.5_EB*X13+X23)
   OM12=.5_EB*(.5_EB*OM1+OM2)
   SDWEAK=0._EB
   SDSTRG=0._EB
   IF(OMEGA<=2395._EB) THEN
         ALPHA=2700._EB
         OMPRIM=OM3
         AA=ALPHA*B*Q2/(A*(1._EB-EXP(-OM3*Q2/T0))*(1._EB-EXP(-OM12*Q2 /T0))**3*(1._EB+EXP(-OM12*Q2/T0))*(1._EB-EXP(-OMPRIM*Q2/T0)))
         BB=(1._EB-EXP(-Q2*OMEGA/TEMP))*(1._EB-EXP(-Q2*OM3/TEMP))*(1._EB-EXP(-OM12*Q2/TEMP))**3*(1._EB+EXP(-OM12*Q2/TEMP)) &
           *(1._EB-EXP(-Q2*OMPRIM/TEMP))
         CC=AA*BB*OMEGA/TEMP*T0/TEMP
         L202A: DO J=1,20
            V=FLOAT(J-1)
            IF(J/2*2==J)G=(V+1._EB)*(V+3._EB)/4._EB
            IF(J/2*2/=J)G=(V+2._EB)*(V+2._EB)/4._EB
            L201A: DO K=1,10
               V3=FLOAT(K-1)
               QQ=(V3+1._EB)*G*EXP(-(V3*OM3+V*OM12)*Q2/TEMP)
               GAM=B-A*(V3+1._EB)
               OMVV3=OM3+.5_EB*X13+X23+2._EB*X33+XBAR*V+2._EB*X33*V3
               DELTA=A*(OMEGA-OMVV3)
               IF(GAM*GAM<=DELTA) CYCLE L202A
               D=2._EB*SQRT(GAM*GAM-DELTA)
               OMVBAR=OMVV3*(1._EB-EXP(-OMVV3*Q2/TEMP))
               F1=GAM-D/2._EB
               F2=GAM+D/2._EB
               EE=Q2*GAM/(A*A*TEMP)
               UNFLO1=EE*DELTA*(1._EB+.5_EB*A/GAM)
               IF(UNFLO1<=-78.) CYCLE L202A
               UNFLO2=EE*2._EB*GAM*F1
               IF(UNFLO2>=78.) CYCLE L202A
               FF=EXP(EE*DELTA*(1._EB+.5_EB*A/GAM))
               SMINUS=CC*QQ/OMVBAR*ABS(F1)*FF*EXP(-EE*2._EB*GAM*F1)
               UNFLO3=EE*2._EB*GAM*F2
               IF(UNFLO3>=78._EB) THEN
                  SPLUS=0._EB
                  ELSE
                  SPLUS=CC*QQ/OMVBAR*ABS(F2)*FF*EXP(-EE*2._EB*GAM*F2)
                  ENDIF
               GG=SDWEAK
               SDWEAK=(SMINUS+SPLUS)/D+SDWEAK
               TEST=(SDWEAK-GG)/SDWEAK
               IF(TEST<.0001_EB) CYCLE L202A
               SDSTRG=SQRT(0.5_EB*G)*(SQRT(SMINUS)+SQRT(SPLUS))/D+SDSTRG
            ENDDO L201A
         ENDDO L202A
         IF(SDWEAK<=ZERO_P) THEN
            SDWEAK=0._EB
            GDINV=1._EB
            GDDINV=1._EB
         ELSE
            DINV=SDSTRG*SDSTRG/SDWEAK
            GDINV=GC1*DINV
            GDDINV=GD*DINV
!***  EXPRESS S/D AT STP, AS IS K IN NASA SP-3080
            SDWEAK=SDWEAK*TEMP/273._EB
         ENDIF
      ELSE
!CALCULATE ABSORPTION COEF. AND LINE SPACING PARAMETER FOR 2.7 MICRON BAND
         L=1
!CONTRIBUTION TO 2.7 MICRON BAND FROM (000)-(021) AND (010)-(031) TRANS.
         ALPHA=28.5_EB
         OMPRIM=2._EB*OM2+OM3
         L120A: DO
         AA=ALPHA*B*Q2/(A*(1._EB-EXP(-OM3*Q2/T0))*(1._EB-EXP(-OM12*Q2 /T0))**3*(1._EB+EXP(-OM12*Q2/T0))*(1._EB-EXP(-OMPRIM*Q2/T0)))
         BB=(1._EB-EXP(-Q2*OMEGA/TEMP))*(1._EB-EXP(-Q2*OM3/TEMP))* (1._EB-EXP(-OM12*Q2/TEMP))**3*(1._EB+EXP(-OM12*Q2/TEMP)) &
           *(1._EB-EXP(-Q2*OMPRIM/TEMP))
         CC=AA*BB*OMEGA/TEMP*T0/TEMP
         L102A: DO J=1,20
            V=FLOAT(J-1)
            IF(J/2*2==J)G=(V+1._EB)*(V+3._EB)/4._EB
            IF(J/2*2/=J)G=(V+2._EB)*(V+2._EB)/4._EB
            VBAR1=-1._EB+(V+3._EB)*(V+4._EB)/(V+2._EB)/6._EB
            IF(J/2*2==J)VBAR1=-1._EB+(V+5._EB)/6._EB
            L101A: DO K=1,10
               V3=FLOAT(K-1)
               QQ=(V3+1)*G*EXP(-(V3*OM3+V*OM12)*Q2/TEMP)*(VBAR1+1._EB)
               GAM=B-A*(V3+1._EB)
               IF(L==2) THEN
                  OMVV3=3728._EB-5._EB*V-47._EB*V3
                  IF(V<=ZERO_P)OMVV3=3715._EB-47._EB*V3
                  ELSE
                  OMVV3=3598._EB-18._EB*V-47._EB*V3
                  IF(V<=ZERO_P)OMVV3=3613._EB-47._EB*V3
               ENDIF
               DELTA=A*(OMEGA-OMVV3)
               IF(GAM*GAM<=DELTA) CYCLE L102A
               D=2._EB*SQRT(GAM*GAM-DELTA)
               OMVBAR=OMVV3*(1._EB-EXP(-OMVV3*Q2/TEMP))
               F1=GAM-D/2._EB
               F2=GAM+D/2._EB
               EE=Q2*GAM/(A*A*TEMP)
               UNFLO1=EE*DELTA*(1._EB+.5_EB*A/GAM)
               IF(UNFLO1<=-78._EB) CYCLE L102A
               UNFLO2=EE*2._EB*GAM*F1
               IF(UNFLO2>=78._EB) CYCLE L102A
               FF=EXP(EE*DELTA*(1._EB+.5_EB*A/GAM))
               SMINUS=CC*QQ/OMVBAR*ABS(F1)*FF*EXP(-EE*2._EB*GAM*F1)
               UNFLO3=EE*2._EB*GAM*F2
               IF(UNFLO3>=78._EB) THEN
                  SPLUS=0._EB
                  ELSE
                  SPLUS=CC*QQ/OMVBAR*ABS(F2)*FF*EXP(-EE*2._EB*GAM*F2)
               ENDIF
               GG=SDWEAK
               SDWEAK=(SMINUS+SPLUS)/D+SDWEAK
               TEST=(SDWEAK-GG)/SDWEAK
               IF(TEST<.0001) CYCLE L102A
               SDSTRG=SQRT(.5_EB*G)*(SQRT(SMINUS)+SQRT(SPLUS))/D+SDSTRG
            ENDDO L101A
         ENDDO L102A
         IF(L==2) EXIT L120A
!CONTRIBUTION TO 2.7 MICRON BAND FROM (000)-(101) AND (010)-(111) TRANS.
         ALPHA=42.3_EB
         OMPRIM=OM1+OM3
         L=2
      ENDDO L120A
!CALCULATE ABSORPTION COEF AND LINE SPACING PARAMETER FOR 4.3 MICRON BAND
      IF(SDWEAK<=ZERO_P) THEN
         SDWEAK=0._EB
         GDINV=1._EB
         GDDINV=1._EB
         ELSE
         DINV=SDSTRG*SDSTRG/SDWEAK
         GDINV=GC1*DINV
         GDDINV=GD*DINV
!***EXPRESS S/D AT STP, AS IS K IN NASA SP-3080
         SDWEAK=SDWEAK*TEMP/273._EB
         ENDIF
      ENDIF
   ELSEIF((OMEGA<=1975._EB).AND.(OMEGA>1100._EB)) THEN
      SDWEAK=0._EB
      GDINV=1._EB
      GDDINV=1._EB
   ELSEIF((OMEGA<=1100._EB).AND.(OMEGA>880._EB)) THEN
!CONTRIBUTION TO 10.0 MICRON BAND FROM (100)-(001) AND (020)-(001) TRANS.
      OM1=1354.91_EB
      OM2=673._EB
      OM3=2396.49_EB
      Q2=1.4388_EB
      BCNT(1)=960.8_EB
      BCNT(2)=1063.6_EB
      OMA=OM3
      OMB=(OM1+2._EB*OM2)/2._EB
      T0=300._EB
      ATOT(1)=0.0219_EB
      ATOT(2)=0.0532_EB
      BE=0.391635_EB
      DO K=1,2
         ATOT(K)=T0/TEMP*ATOT(K)*EXP(Q2*OMB*(1._EB/T0-1._EB/TEMP)) *(1._EB-EXP(-Q2*(OMA-OMB)/TEMP))/(1._EB-EXP(-Q2*OMA/TEMP)) &
                   /(1._EB-EXP(-OMB*Q2/TEMP))
      ENDDO
      SDWEAK=0._EB
      DO I=1,2
         SDWEAK=SDWEAK+ATOT(I)*Q2/(4._EB*BE*TEMP)*ABS(OMEGA-BCNT(I))  *EXP(-Q2/(4._EB*BE*TEMP)*(OMEGA-BCNT(I))**2)
      ENDDO
      DINV=1._EB/4._EB/BE
      GDINV=GC1*DINV
      GDDINV=GD*DINV
!***EXPRESS S/D AT STP, AS IS IN NASA SP-3080
      SDWEAK=SDWEAK*TEMP/273.
   ELSEIF((OMEGA<=880._EB).AND.(OMEGA>500._EB))  THEN
!CONTRIBUTION TO 15.0 MICRON BAND FROM (000)-(010) TRANS.
      TTEMP=TEMP
      J=(OMEGA-495._EB)/5._EB
      W1=495._EB+5._EB*REAL(J,EB)
      WW=(OMEGA-W1)/5
      IF(TEMP>=2400._EB) TEMP = 2399.99_EB
      IF(TEMP < 300._EB) TEMP =  300.00_EB
      I = TEMP/300._EB
      SELECT CASE(I)
         CASE(3)
            I=2
            TT=(TEMP-600._EB)/600._EB
         CASE(6:7) 
            I=5
            TT=(TEMP-1800._EB)/600._EB
         CASE DEFAULT
            T1=REAL(I,EB)*300._EB
            TT=(TEMP-T1)/300._EB
            IF (I>=4) I=I-1     
      END SELECT
      TW=TT*WW
      SDWEAK=SD15(I,J)*(1._EB-TT-WW+TW)+SD15(I+1,J)*(TT-TW)  +SD15(I,J+1)*(WW-TW)+SD15(I+1,J+1)*TW
      IF(SDWEAK<=ZERO_P) THEN
         SDWEAK=0._EB
         GDINV=1._EB
         GDDINV=1._EB
         ELSE
!CALCULATE LINE SPACING PARAMETER FOR 15.0 MICRON BAND
         DINV1=1.2_EB
         DINV2=8.0_EB
         DINV3=30.0_EB
         TEMP1=300.0_EB
         TEMP2=550.0_EB
         TEMP3=830.0_EB
         DINV=DINV1*(TEMP-TEMP2)*(TEMP-TEMP3)/(TEMP1-TEMP2)  /(TEMP1-TEMP3)+DINV2*(TEMP-TEMP1)*(TEMP-TEMP3)&
                  /(TEMP2-TEMP1)/(TEMP2-TEMP3)+DINV3*(TEMP-TEMP1)  *(TEMP-TEMP2)/(TEMP3-TEMP1)/(TEMP3-TEMP2)
         GDINV=GC1*DINV
         GDDINV=GD*DINV
         ENDIF
         TEMP = TTEMP  ! Line added by Jason Floyd, Aug 30, 2002
   ELSE
      SDWEAK=0._EB
      GDINV=1._EB
      GDDINV=1._EB
      ENDIF
   ENDIF
!------------------------------------------------------------------------------
END SUBROUTINE CO2


!==============================================================================
SUBROUTINE H2O(OMEGA,TEMP,GC2,SDWEAK,GDINV,GDDINV)
!==============================================================================
INTEGER I,J
REAL(EB) OMEGA,TEMP,GC2,SDWEAK,GDINV,GDDINV,WM,W1,WW,T1,TT,TW, D,B,DINV,TTEMP,GD

IF (OMEGA>=9300..OR.OMEGA<50._EB) THEN
   SDWEAK=0._EB
   GDINV=1._EB
   GDDINV=1._EB
ELSE
   WM   = 18._EB
   GD   = 5.94E-6_EB*OMEGA*SQRT(TEMP/(273._EB*WM))
   J    = (OMEGA-25._EB)/25._EB
   TTEMP= TEMP

   IF(TEMP>=2500._EB) TEMP = 2499.99_EB
   IF(TEMP<300._EB)   TEMP = 300._EB

   I = TEMP/500._EB +1

   IF(I==2.AND.TEMP<600._EB) I=1

   W1 = 25._EB+25._EB*FLOAT(J)
   WW = (OMEGA-W1)/25._EB

   IF(I>2) THEN
      T1=FLOAT(I-1)*500._EB
      TT=(TEMP-T1)/500._EB
   ELSE
      IF(I==1) TT=(TEMP-300._EB)/300._EB
      IF(I==2) TT=(TEMP-600._EB)/400._EB
   ENDIF

   TW     = TT*WW
   SDWEAK = SD(I,J)*(1._EB-TT-WW+TW)+SD(I+1,J)*(TT-TW)+SD(I,J+1) *(WW-TW)+SD(I+1,J+1)*TW
   D      = -2.294_EB+.3004E-02_EB*TEMP-.366E-06_EB*TEMP**2
   B      = SIN(.0036_EB*OMEGA-8.043_EB)
   DINV   = EXP(.7941_EB*B+D)
!     DINV=EXP(0.00106*TEMP-1.21)
   GDINV  = GC2*DINV
   GDDINV = GD*DINV
   TEMP   = TTEMP
ENDIF
!------------------------------------------------------------------------------
END SUBROUTINE H2O


!==============================================================================
SUBROUTINE CO(OMEGA,TEMP,GC4,SDWEAK,GDINV,GDDINV)
!==============================================================================
INTEGER J
REAL(EB) OMEGA,TEMP,GC4,SDWEAK,GDINV,GDDINV,AA,BB,CC,QQ,EE,FF,GG, &
         SMINUS,SPLUS,SDSTRG,B,ALPHA,A,OME,WX,WY,OMPRIM,T0, &
         Q2,WM,GD,V,GAM,OMV,DELTA,D,OMVBAR,F1,F2,TEST,DINV,Q2OT,TOAZ

IF(OMEGA<1600._EB .OR. OMEGA>2400._EB) THEN
   SDWEAK = 0._EB
   GDINV  = 1._EB
   GDDINV = 1._EB
ELSE
   B     = 1.93139_EB
   ALPHA = 260._EB
   A     = .017485_EB
   OME   = 2170.21_EB
   WX    = 13.461_EB
   WY    = .0308_EB
   OMPRIM= OME-2._EB*WX+3.25_EB*WY

   T0 = 300._EB
   Q2 = 1.4388_EB

   TOAZ = TEMP/273._EB
   Q2OT = Q2/TEMP

   WM     = 28._EB
! DOPPLER BROADENING HALF-WIDTH
   GD     = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM)
   SDWEAK = 1.E-99_EB
   SDSTRG = 1.E-99_EB

   AA = ALPHA*B*Q2/(A*(1._EB-EXP(-OMPRIM*Q2/T0))**2)
   BB = (1._EB-EXP(-OMEGA*Q2OT))*(1._EB-EXP(-OMPRIM*Q2OT))**2
   CC = AA*BB*OMEGA*T0/TEMP**2

   L101: DO J=1,20
      V    = FLOAT(J-1)
      QQ   = (V+1._EB)*EXP(-V*OME*Q2OT)
      GAM  = B-A*(V+1._EB)
      OMV  = OME-2._EB*(V+1._EB)*WX+(3._EB*(V+1._EB)*(V+1._EB)+.25_EB)*WY
      DELTA= A*(OMEGA-OMV)

      IF(GAM**2<=DELTA) EXIT L101

      D      = 2._EB*SQRT(GAM*GAM-DELTA)
      OMVBAR = OMV*(1._EB-EXP(-OMV*Q2OT))
      F1     = GAM-0.5_EB*D
      F2     = GAM+0.5_EB*D
      EE     = Q2*GAM/(A*A*TEMP)
      
      FF     = EXP(EE*DELTA*(1._EB+.5_EB*A/GAM))
      SMINUS = CC*QQ/OMVBAR*ABS(F1)*FF*EXP(-EE*2._EB*GAM*F1)
      SPLUS  = CC*QQ/OMVBAR*ABS(F2)*FF*EXP(-EE*2._EB*GAM*F2)
      GG     = SDWEAK
      SDWEAK = (SMINUS+SPLUS)/D+SDWEAK
      TEST   = (SDWEAK-GG)/SDWEAK
      
      IF(TEST<.0001_EB) EXIT L101

      SDSTRG=(SQRT(SMINUS)+SQRT(SPLUS))/D+SDSTRG
   ENDDO L101

   DINV   = SDSTRG*SDSTRG/SDWEAK
   GDINV  = GC4*DINV
   GDDINV = GD*DINV
!***EXPRESS S/D AT STP, AS IS K IN NASA SP-3080
   SDWEAK = SDWEAK*TOAZ
ENDIF
!------------------------------------------------------------------------------
END SUBROUTINE CO


!==============================================================================
SUBROUTINE POD(OMEGA)
!==============================================================================
! POD CALCULATES PARTICLE OPTICAL DEPTH, XPART, OF THE VOLUME FRACTION OF SOOT PARTICLES IN GAS CLOUD.  RIN AND RIK ARE
! THE REAL AND IMAGINARY PARTS OF THE INDEX OF REFRACTION. THE PARTICLES ARE ASSUMED TO BE IN THE RAYLEIGH LIMIT.

REAL(EB) OMEGA,ABCO,FF_FAC,FF,LAMBDA!,RIN,RIK

LAMBDA=10000._EB/OMEGA
 
! ABSORPTION COEF. IS BASED UPON MEASUREMENTS OF WIDMANN AND MULHOLLAND

IF (AL2O3) THEN
   IF (RCT > 2570._EB) THEN
      FF_FAC = 0.017_EB
   ELSEIF (RCT < 500._EB) THEN
      FF_FAC = 5.E-6_EB
   ELSE
      FF_FAC = 0.00073_EB     
   ENDIF
ELSE
   FF_FAC = 8.9_EB
ENDIF

FF=FF_FAC/LAMBDA
ABCO=FF*SVF*1.E6_EB
XPART=ABCO*DD
!------------------------------------------------------------------------------
END SUBROUTINE POD


!==============================================================================
SUBROUTINE CH4(OMEGA,TEMP,PCH4,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
INTEGER I,J
REAL(EB) OMEGA,TEMP,PCH4,PTOT,GC3,SDWEAK,GDINV,GDDINV,BE,Q2, &
      WM,GD,OM1,OM2,OM3,OM4,COM1,COM2,COM3,COM4,DINV,PE,W1,SDB, &
      SDA,SDC,Q2OT,AZOT,TOAZ

IF(OMEGA>5000._EB .OR. OMEGA<1125._EB) THEN

   SDWEAK=0.0_EB
   GDINV=1._EB
   GDDINV=1._EB

ELSE

   BE=5.2412_EB
   Q2=1.4388_EB
   WM=16._EB
   Q2OT = -Q2/TEMP
   AZOT = 273._EB/TEMP
   TOAZ = TEMP/273._EB

   GD=5.94E-6_EB*OMEGA*SQRT(TOAZ)/4._EB

   IF(OMEGA>=3400._EB) THEN

! CONTRIBUTION TO 2.4 MICRON BAND FROM (0000)-(0110), (0000)-(0011),
! (0000)-(1001), AND (0000)-(0102) TRANS.  THE INTEGRATED BAND INTENSITIES
! OF VINCENT-GEISSE (ANNALES DE PHYSIQUE SER.12, V. 10, 1955) HAVE
! BEEN MULTIPLIED BY A FACTOR OF 4 AND THE LINE SPACING IS THAT
! OF V4 FROM GRAY AND PENNER (JQSRT V. 5, 1965).

      OM1=2914.2_EB
      OM2=1526.0_EB
      OM3=3020.3_EB
      OM4=1306.2_EB

      BCNT(1)=4123.0_EB
      BCNT(2)=4216.3_EB
      BCNT(3)=4313.2_EB
      BCNT(4)=4546.0_EB

      COM1=OM2+2._EB*OM4
      COM2=OM1+OM4
      COM3=OM3+OM4
      COM4=OM2+OM3

      ATOT(1)=0.64_EB*AZOT*(1._EB-EXP(Q2OT*COM1))/((1._EB-EXP(Q2OT*OM2))*(1._EB-EXP(Q2OT*OM4))**2)
      ATOT(2)=17.6_EB*AZOT*(1._EB-EXP(Q2OT*COM2))/((1._EB-EXP(Q2OT*OM1))*(1._EB-EXP(Q2OT*OM4)))
      ATOT(3)=14.8_EB*AZOT*(1._EB-EXP(Q2OT*COM3))/((1._EB-EXP(Q2OT*OM3))*(1._EB-EXP(Q2OT*OM4)))
      ATOT(4)=5.04_EB*AZOT*(1._EB-EXP(Q2OT*COM4))/((1._EB-EXP(Q2OT*OM2))*(1._EB-EXP(Q2OT*OM3)))

      DINV=1._EB/5.74_EB
      GDINV=GC3*DINV
      GDDINV=GD*DINV
      SDWEAK=0.0_EB

      DO I=1,4
         SDWEAK=SDWEAK+2._EB*(OMEGA-BCNT(I))**2*(-Q2OT*BE)*SQRT(-Q2OT*BE)*ATOT(I)/SQRTPI*DINV**3*EXP(Q2OT*BE*DINV**2 &
         *(OMEGA-BCNT(I))**2)
      ENDDO
      SDWEAK=SDWEAK*TOAZ
   ELSE
      PE=PTOT+.3_EB*PCH4

      IF(OMEGA>=2625._EB) THEN
! CONTRIBUTION TO 3.3 MICRON BAND FROM (0000)-(0010) TRANS.
! REFER TO BROSMER AND TIEN, JQSRT V. 33, P. 521

         GDINV  = .00734_EB*PE*SQRT(AZOT)*EXP(1.02_EB*(TOAZ-1._EB))
         GDDINV = GD/9.4_EB

         J  = (OMEGA-2600._EB)/25._EB
         W1 = 2600._EB+25._EB*FLOAT(J)
         SDB= SD3(2,J)+(OMEGA-W1)/25._EB*(SD3(2,J+1)-SD3(2,J))

         IF(TEMP>600._EB) THEN
            SDC=SD3(3,J)+(OMEGA-W1)/25._EB*(SD3(3,J+1)-SD3(3,J))
            SDWEAK=SDB+(TEMP-600._EB)/250._EB*(SDC-SDB)
            IF(SDWEAK<0._EB)SDWEAK=0._EB
         ELSE
            SDA=SD3(1,J)+(OMEGA-W1)/25._EB*(SD3(1,J+1)-SD3(1,J))
            SDWEAK=SDA+(TEMP-290._EB)/310._EB*(SDB-SDA)
            IF(SDWEAK<0._EB)SDWEAK=0._EB
         ENDIF

      ELSEIF(OMEGA>1450._EB) THEN

         SDWEAK=0.0_EB
         GDINV=1._EB
         GDDINV=1._EB

      ELSE
! CONTRIBUTION TO 7.7 MICRON BAND FROM (0000)-(0001) TRANS.
! REFER TO BROSMER AND TIEN, JQSRT V. 33, P. 521.
         GDINV  = .0243_EB*PE*(TOAZ)**.8_EB
         GDDINV = GD/5.1_EB

         J   = (OMEGA-1100._EB)/25._EB
         W1  = 1100._EB+25._EB*FLOAT(J)
         SDB = SD7(2,J)+(OMEGA-W1)/25._EB*(SD7(2,J+1)-SD7(2,J))

         IF(TEMP>600._EB) THEN
            SDC = SD7(3,J)+(OMEGA-W1)/25._EB*(SD7(3,J+1)-SD7(3,J))
            SDWEAK = SDB+(TEMP-600._EB)/250._EB*(SDC-SDB)
            IF(SDWEAK<0._EB)SDWEAK=0._EB
         ELSE
            SDA = SD7(1,J)+(OMEGA-W1)/25._EB*(SD7(1,J+1)-SD7(1,J))
            SDWEAK = SDA+(TEMP-290._EB)/310._EB*(SDB-SDA)
            IF(SDWEAK<0._EB)SDWEAK=0._EB
         ENDIF
      ENDIF
   ENDIF
ENDIF
!------------------------------------------------------------------------------
END SUBROUTINE CH4


!==============================================================================
SUBROUTINE CH4_NEW(OMEGA,TEMP,PCH4,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
! COMPUTE METHANE OPTICAL PROPERTIES USING SINGLE LINE GROUP
!
! BAND #1: 1150 cm-1 - 1600 cm-1 
! BAND #2: 2700 cm-1 - 3250 cm-1 
! BAND #3: 3400 cm-1 - 5000 cm-1
! 
! VARIABLES PASSED IN
! OMEGA: (REAL) WAVENUMBER IN CM-1
! TEMP : (REAL) TEMPERATURE IN K
! PCH4 : (REAL) PARTIAL PRESSURE OF METHANE IN ATM
! PTOT : (REAL) TOTAL PRESSURE IN ATM
! GC3  : (REAL) collisional broadening half-width at half heigh
! 
! VRAIBLES PASSED OUT
! SDWEAK : (REAL) SPECTRAL ABSORPTION COEFFICIENT
! GDINV  : (REAL) LINE WIDTH TO LINE SPACING RATIO
! GDDINV : (REAL) LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING
!
! LOCALS
! PRESSURE_EFFECTIVE : (REAL) EFFECTIVE PRESSURE, (FORMELY PE)
! BE_CH4             : (REAL) ROTATIONAL CONSTANT FOR CH4 [CM^-1]
! DINV_CH4           : INVERSE LINE SPACING [CM]

REAL(EB), INTENT(IN)  :: OMEGA, TEMP, PCH4, PTOT, GC3
REAL(EB), INTENT(OUT) :: SDWEAK, GDINV, GDDINV

REAL(EB) :: GD, PRESSURE_EFFECTIVE, Q2OT, AZOT, TOAZ, FACT1
REAL(EB) :: DINV_CH4

REAL(EB), PARAMETER :: Q2     = 1.4388_EB ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
REAL(EB), PARAMETER :: WM_CH4 = 16.0425_EB

!LINE SPACING D2_CH4=5.74 CM^-1 FOR THE V4 FUNDAMENTAL OF METHANE GRAY & PENNER,612
!FOOTNOTE AT BOTTOM
REAL(EB), PARAMETER :: BE_CH4 = 5.248_EB   ! ROTATIONAL CONSTANTS [CM^-1]

REAL(EB), DIMENSION(4), PARAMETER :: &
                 OM_CH4  = (/2914.2_EB, 1526.0_EB, 3020.3_EB, 1306.2_EB/),                                               &
                 COM_CH4 = (/1526.0_EB+2._EB*1306.2_EB, 2914.2_EB+1306.2_EB, 3020.3_EB+1306.2_EB, 1526.0_EB+3020.3_EB/), &
                 S2_CH4  = (/0.64_EB,17.6_EB,14.8_EB,5.04_EB/)

INTEGER I

! INITIALIZE OUTPUT

SDWEAK = 0.0_EB  !SPECTRAL ABSORPTION COEFFICIENT. 
GDINV  = 1._EB   !LINE WIDTH TO LINE SPACING RATIO: FINE STRUCTURE PARAMETER
GDDINV = 1._EB   !LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING FINE STRUCTURE PARAMETER

! INITIALIZE SOME COEFFICIENTS
Q2OT = -Q2/TEMP       
AZOT = 273._EB/TEMP
TOAZ = TEMP/273._EB

! COMPUTE DOPLLER HALF-WIDTH. EQ 5-35
GD   = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM_CH4) !DOPPLER HALF WIDTH [CM^-1]. NASA,222

! COMPUTE EFFECTIVE PRESSURE. BROSMER & TIEN,524 EQ. 7 (SPECIFIC TO METHANE)
PRESSURE_EFFECTIVE = PTOT+.3_EB*PCH4

!------------------------------------------------------------------------------
! COMPUTED PROPERTIES 2.4 MICRON BAND, BAND #3: 3400 cm-1 - 5000 cm-1
IF((OM_BND_CH4(N_BAND_CH4,1)<=OMEGA).AND.(OMEGA<=OM_BND_CH4(N_BAND_CH4,2))) THEN

! CONTRIBUTION OF THE V1+V4 COMBINAISON BAND (2.4 MICRON)
! CONTRIBUTION TO 2.4 MICRON BAND FROM
! (0000)-(0102)
! (0000)-(1001)
! (0000)-(0011)
! (0000)-(0110) TRANSITIONS
! THE INTEGRATED BAND INTENSITIES (S2_CH4)
! OF VINCENT-GEISSE (ANNALES DE PHYSIQUE SER.12, V. 10, 1955) HAVE
! BEEN MULTIPLIED BY A FACTOR OF 4 
! THE LINE SPACING (1/DINV_CH4) IS THAT
! OF V4 FROM GRAY AND PENNER (JQSRT V. 5, 1965).                 

   ATOT(1) = S2_CH4(1)*AZOT*(1._EB-EXP(Q2OT*COM_CH4(1)))/((1._EB-EXP(Q2OT*OM_CH4(2)))*(1._EB-EXP(Q2OT*OM_CH4(4)))**2)
   ATOT(2) = S2_CH4(2)*AZOT*(1._EB-EXP(Q2OT*COM_CH4(2)))/((1._EB-EXP(Q2OT*OM_CH4(1)))*(1._EB-EXP(Q2OT*OM_CH4(4))))
   ATOT(3) = S2_CH4(3)*AZOT*(1._EB-EXP(Q2OT*COM_CH4(3)))/((1._EB-EXP(Q2OT*OM_CH4(3)))*(1._EB-EXP(Q2OT*OM_CH4(4))))
   ATOT(4) = S2_CH4(4)*AZOT*(1._EB-EXP(Q2OT*COM_CH4(4)))/((1._EB-EXP(Q2OT*OM_CH4(2)))*(1._EB-EXP(Q2OT*OM_CH4(3))))

   DINV_CH4 = 1._EB/5.74_EB
   FACT1  = Q2OT*BE_CH4*DINV_CH4**2

   DO I=1,4
      SDWEAK = SDWEAK+(OMEGA-COM_CH4(I))**2*ATOT(I)*EXP(FACT1*(OMEGA-COM_CH4(I))**2)
   ENDDO

! EQ 6 IN GRAY AND PENNER JQSRT, VOL 5, PAGE 611-620, 1965
   SDWEAK=SDWEAK*2._EB*(-Q2OT*BE_CH4)*SQRT(-Q2OT*BE_CH4)/SQRTPI*DINV_CH4**3

!***EXPRESS S/D AT STANDARD TEMPERATURE AND PRESSURE, AS IS IN NASA SP-3080
   SDWEAK = SDWEAK*TOAZ
   GDINV  = GC3*DINV_CH4
   GDDINV = GD*DINV_CH4 

!------------------------------------------------------------------------------
! TABULATED PROPERTIES
ELSE IF((OM_BND_CH4(2,1)<=OMEGA).AND.(OMEGA<OM_BND_CH4(2,2))) THEN
      
!------------------------------------------------------------------------------
! CONTRIBUTION TO 3.3 MICRON BAND, BAND #2: 2700 cm-1 - 3250 cm-1
! (0000)-(0010) TRANSITION (V3 FUNDAMENTAL)
! 9.4 IS THE AVERAGE LINE POSITION FOR THE 3.3 MICRON BAND (CM^-1)
! SOURCE: BROSMER & TIEN, JQSRT V. 33, P. 525 (SPECIFIC TO METHANE)

   DINV_CH4 = 1._EB/9.4_EB

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH4_TEMP,OM_BND_CH4(2,:),SD2_CH4)
! LINE SHAPE PARAMETER BROSMER & TIEN, JQSRT V. 33, P. 525 EQ. 10 (SPECIFIC TO METHANE)
   GDINV  = .00734_EB*PRESSURE_EFFECTIVE*SQRT(AZOT)*EXP(1.02_EB*(TOAZ-1._EB)) 
   GDDINV = GD*DINV_CH4

ELSE IF((OM_BND_CH4(1,1)<=OMEGA).AND.(OMEGA<OM_BND_CH4(1,2))) THEN
!------------------------------------------------------------------------------
! CONTRIBUTION TO 7.7 MICRON BAND, BAND #1: 1150 cm-1 - 1600 cm-1 
! (0000)-(0001) TRANSITION (V4 FUNDAMENTAL)
! 5.1 (CM-1) IS THE AVERAGE LINE POSITION FOR THE 7.7 MICRON BAND 
! SOURCE: BROSMER & TIEN, JQSRT V. 33, P. 525 (SPECIFIC TO METHANE)

   DINV_CH4 = 1._EB/5.1_EB

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH4_TEMP,OM_BND_CH4(1,:),SD1_CH4)
! LINE SHAPE PARAMETER BROSMER & TIEN, JQSRT V. 33, 525 EQ. 11 (SPECIFIC TO METHANE)
   GDINV  = .0243_EB*PRESSURE_EFFECTIVE*(TOAZ)**.8_EB
   GDDINV = GD*DINV_CH4

ENDIF 

!------------------------------------------------------------------------------
END SUBROUTINE CH4_NEW

!==============================================================================
SUBROUTINE C3H8(OMEGA,TEMP,PC3H8,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
! COMPUTE PROPANE OPTICAL PROPERTIES USING SINGLE LINE GROUP
! 
! ABSORPTION BANDS: 3.3 MICRON BAND
!                   7 MICRON BAND
!
! VARIABLES PASSED IN
! OMEGA : (REAL) WAVENUMBER IN CM-1
! TEMP  : (REAL) TEMPERATURE IN K
! PC3H8 : (REAL) PARTIAL PRESSURE OF PROPANE
! PTOT  : (REAL) TOTAL PRESSURE
! GC3   : (REAL) collisional broadening half-width at half heigh
! 
! VRAIBLES PASSED OUT
! SDWEAK : (REAL) SPECTRAL ABSORPTION COEFFICIENT
! GDINV  : (REAL) LINE WIDTH TO LINE SPACING RATIO,  FINE STRUCTURE PARAMETER
! GDDINV : (REAL) LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING,
!                 FINE STRUCTURE PARAMETER
!
! LOCALS
! PRESSURE_EFFECTIVE : (REAL) EFFECTIVE PRESSURE, (FORMELY PE)
! DINV_C3H8          : INVERSE LINE SPACING [CM] FOR PROPANE

REAL(EB), INTENT(IN)  :: OMEGA, TEMP, PC3H8, PTOT, GC3
REAL(EB), INTENT(OUT) :: SDWEAK, GDINV, GDDINV

REAL(EB) :: Q2OT, AZOT, TOAZ
REAL(EB) :: GD, PRESSURE_EFFECTIVE
REAL(EB) :: DINV_C3H8 

REAL(EB), PARAMETER :: Q2      = 1.4388_EB  ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
REAL(EB), PARAMETER :: WM_C3H8 = 44.0956_EB ! NIST WEBBOOK DATA

! COMPUTE THERMAL COEFFICIENTS

Q2OT = -Q2/TEMP       
AZOT = 273._EB/TEMP
TOAZ = TEMP/273._EB

! COMPUTE DOPPLER BROADENING HALF WIDTH GD
GD   = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM_C3H8) !DOPPLER HALF WIDTH [CM^-1]. NASA,222

! LINE SPACING
! PROPANE IS A PROLATE SYMMETRIC TOP. THE LINE SPACING IS APPROXIMATED TO
! 2*ROTATIONAL_CONSTANT B. RECALL THAT ROTATIONAL_CONSTANT B = ROTATIONAL_CONSTANT C
! ROTATIONAL_CONSTANT_B OBTAINED FROM CCCBDV WEBSITE, USING EXPERIMENTAL VALUES  
! PROPANE BELONGS TO GROUP C2v
DINV_C3H8 = 1._EB/(2._EB*0.28173_EB)

! SET INITIAL VALUES TO SDWEAK, GDINV, GDDINV - VALUES ARE RETURNED IF OMEGA IS NOT
! WITHIN ABSOPRTION BANDS
SDWEAK = 0.0_EB
GDINV  = 1._EB 
GDDINV = 1._EB 

! COMPUTE PRESSURE EFFECTIVE. FROM BROSMER & TIEN,524 EQ. 7 (SPECIFIC TO METHANE)
! NEED TO FIND PROPANE EFFECTIVE PRESSURE FORMULATION FOR PROPANE

PRESSURE_EFFECTIVE = PTOT+0.3_EB*PC3H8 ! (SPECIFIC TO METHANE)

! IF OMEGA IS WITHIN ABSOPRTION BAND, THEN COMPUTE SDWEAK FROM TABULATED DATA
! PROPANE HAS ONLY TWO BANDS
! 7 MICRON AND 3.3 MICRON


IF ((OM_BND_C3H8(1,1)<=OMEGA).AND.(OMEGA<OM_BND_C3H8(1,2))) THEN
! BAND #1: 1300 cm-1 - 1600 cm-1  7 MICRON BAND
   
   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C3H8_TEMP,OM_BND_C3H8(1,:),SD1_C3H8)

!! LINE SHAPE PARAMETER BROSMER & TIEN, JQSRT V. 33, 525 EQ. 11 (SPECIFIC TO METHANE)
!  GDINV  = .0243_EB*PRESSURE_EFFECTIVE*(TOAZ)**.8_EB

   GDINV  = GC3*DINV_C3H8
   GDDINV = GD*DINV_C3H8

ELSE IF ((OM_BND_C3H8(2,1)<=OMEGA).AND.(OMEGA<OM_BND_C3H8(2,2))) THEN
! BAND #2: 2700 cm-1 - 3200 cm-1  3.3 MICRON BAND
   
   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C3H8_TEMP,OM_BND_C3H8(2,:),SD2_C3H8)

!! LINE SHAPE PARAMETER BROSMER & TIEN, JQSRT V. 33, P. 525 EQ. 10 (SPECIFIC TO METHANE)
!  GDINV  = .00734_EB*PRESSURE_EFFECTIVE*SQRT(AZOT)*EXP(1.02_EB*(TOAZ-1._EB)) 

   GDINV  = GC3*DINV_C3H8
   GDDINV = GD*DINV_C3H8
ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE C3H8

!==============================================================================
SUBROUTINE C7H16(OMEGA,TEMP,PC7H16,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
! COMPUTE HEPTANE OPTICAL PROPERTIES USING SINGLE LINE GROUP
! 
! ABSORPTION BANDS: 3.3 MICRON BAND
!                   7 MICRON BAND
!
! VARIABLES PASSED IN
! OMEGA  : (REAL) WAVENUMBER IN CM-1
! TEMP   : (REAL) TEMPERATURE IN K
! PC7H16 : (REAL) PARTIAL PRESSURE OF HEPTANE
! PTOT   : (REAL) TOTAL PRESSURE
! GC3    : (REAL) collisional broadening half-width at half heigh
! 
! VRAIBLES PASSED OUT
! SDWEAK : (REAL) SPECTRAL ABSORPTION COEFFICIENT
! GDINV  : (REAL) LINE WIDTH TO LINE SPACING RATIO,  FINE STRUCTURE PARAMETER
! GDDINV : (REAL) LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING,
!                 FINE STRUCTURE PARAMETER
!
! LOCALS
! PRESSURE_EFFECTIVE : (REAL) EFFECTIVE PRESSURE, (FORMELY PE)
! DINV_C7H16         : INVERSE LINE SPACING [CM] FOR HEPTANE

REAL(EB), INTENT(IN)  :: OMEGA, TEMP, PC7H16, PTOT, GC3
REAL(EB), INTENT(OUT) :: SDWEAK, GDINV, GDDINV

REAL(EB) :: Q2OT, AZOT, TOAZ
REAL(EB) :: GD, PRESSURE_EFFECTIVE
REAL(EB) :: DINV_C7H16 

REAL(EB), PARAMETER :: Q2       = 1.4388_EB ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
REAL(EB), PARAMETER :: WM_C7H16 = 100.2019_EB ! NIST WEBBOOK DATA

! COMPUTE THERMAL COEFFICIENTS

Q2OT = -Q2/TEMP       
AZOT = 273._EB/TEMP
TOAZ = TEMP/273._EB

! COMPUTE DOPPLER BROADENING HALF WIDTH GD
GD   = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM_C7H16) !DOPPLER HALF WIDTH [CM^-1]. NASA,222

! COMPUTE AVERAGE LINE SPACING
! n_HEPTANE IS A PROLATE SYMMETRIC TOP. THE LINE SPACING IS APPROXIMATED TO
! 2*ROTATIONAL_CONSTANT B. RECALL THAT ROTATIONAL_CONSTANT B = ROTATIONAL_CONSTANT C
! ROTATIONAL_CONSTANT_B OBTAINED FROM CCCBDV WEBSITE, USING EXPERIMENTAL VALUES  
! n-Heptane BELONGS TO GROUP C2v
DINV_C7H16=1._EB/(2.0_EB*0.024_EB)

! SET INITIAL VALUES TO SDWEAK, GDINV, GDDINV - VALUES ARE RETURNED IF OMEGA IS NOT
! WITHIN ABSOPRTION BANDS
SDWEAK = 0.0_EB
GDINV  = 1._EB 
GDDINV = 1._EB 

! COMPUTE PRESSURE EFFECTIVE. FROM BROSMER & TIEN,524 EQ. 7 (SPECIFIC TO METHANE)
! NEED TO FIND HEPTANE EFFECTIVE PRESSURE FORMULATION FOR HEPTANE
PRESSURE_EFFECTIVE = PTOT+0.3_EB*PC7H16 ! (SPECIFIC TO METHANE)

! IF OMEGA IS WITHIN ABSOPRTION BAND, THEN COMPUTE SDWEAK FROM TABULATED DATA
! HEPTANE HAS ONLY TWO BANDS
! 7 MICRON AND 3.3 MICRON


IF ((OM_BND_C7H16(1,1)<=OMEGA).AND.(OMEGA<OM_BND_C7H16(1,2))) THEN
! BAND #1: 1300 cm-1 - 1550 cm-1 7 MICRON BAND 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H16_TEMP,OM_BND_C7H16(1,:),SD1_C7H16)
   GDINV  = GC3*DINV_C7H16
   GDDINV = GD*DINV_C7H16

ELSE IF ((OM_BND_C7H16(2,1)<=OMEGA).AND.(OMEGA<OM_BND_C7H16(2,2))) THEN
!  BAND #2: 2700 cm-1 - 3200 cm-1  3.3 MICRON BAND

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H16_TEMP,OM_BND_C7H16(2,:),SD2_C7H16)
   GDINV  = GC3*DINV_C7H16
   GDDINV = GD*DINV_C7H16

ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE C7H16

!==============================================================================
SUBROUTINE CH3OH(OMEGA,TEMP,PCH3OH,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
! COMPUTE METHANOL OPTICAL PROPERTIES USING SINGLE LINE GROUP
! 
! THERE ARE 3 BANDS FOR Methanol
!
! BAND #1: 2700 cm-1 - 3200 cm-1 
! BAND #2: 3600 cm-1 - 3800 cm-1 
! BAND #3: 900 cm-1 - 1600 cm-1 
!
! VARIABLES PASSED IN
! OMEGA  : (REAL) WAVENUMBER IN CM-1
! TEMP   : (REAL) TEMPERATURE IN K
! PCH3OH : (REAL) PARTIAL PRESSURE OF METHANOL
! PTOT   : (REAL) TOTAL PRESSURE
! GC3    : (REAL) collisional broadening half-width at half heigh
! 
! VRAIBLES PASSED OUT
! SDWEAK : (REAL) SPECTRAL ABSORPTION COEFFICIENT
! GDINV  : (REAL) LINE WIDTH TO LINE SPACING RATIO,  FINE STRUCTURE PARAMETER
! GDDINV : (REAL) LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING,
!                 FINE STRUCTURE PARAMETER
!
! LOCALS
! PRESSURE_EFFECTIVE : (REAL) EFFECTIVE PRESSURE, (FORMELY PE)
! DINV_CH3OH         : INVERSE LINE SPACING [CM] FOR METHANOL

REAL(EB), INTENT(IN)  :: OMEGA, TEMP, PCH3OH, PTOT, GC3
REAL(EB), INTENT(OUT) :: SDWEAK, GDINV, GDDINV

REAL(EB) :: Q2OT, AZOT, TOAZ
REAL(EB) :: GD, PRESSURE_EFFECTIVE
REAL(EB) :: DINV_CH3OH 

REAL(EB), PARAMETER :: Q2       = 1.4388_EB  ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
REAL(EB), PARAMETER :: WM_CH3OH = 32.0419_EB ! NIST WEBBOOK DATA

! COMPUTE THERMAL COEFFICIENTS

Q2OT = -Q2/TEMP       
AZOT = 273._EB/TEMP
TOAZ = TEMP/273._EB

! COMPUTE DOPPLER BROADENING HALF WIDTH GD
GD   = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM_CH3OH) !DOPPLER HALF WIDTH [CM^-1]. NASA,222

! COMPUTE AVERAGE LINE SPACING
! METHANOL BELONGS Cs POINT GROUP
DINV_CH3OH=1._EB

! SET INITIAL VALUES TO SDWEAK, GDINV, GDDINV - VALUES ARE RETURNED IF OMEGA IS NOT
! WITHIN ABSOPRTION BANDS
SDWEAK = 0.0_EB
GDINV  = 1._EB 
GDDINV = 1._EB 

PRESSURE_EFFECTIVE = PTOT+0.3*PCH3OH

! IF OMEGA IS WITHIN ABSOPRTION BAND, THEN COMPUTE SDWEAK FROM TABULATED DATA

! THERE ARE 3 BANDS FOR Methanol
! BAND #1: 2700 cm-1 - 3200 cm-1, 3.5  MICRON BAND,  DATA IN SD1_CH3OH
! BAND #2: 3600 cm-1 - 3800 cm-1, 2.9  MICRON BAND,  DATA IN SD2_CH3OH
! BAND #3:  900 cm-1 - 1600 cm-1, 10.0 & 8.0 MICRON BANDS, DATA IN SD3_CH3OH

IF ((OM_BND_CH3OH(1,1)<=OMEGA).AND.(OMEGA<OM_BND_CH3OH(1,2))) THEN
! FIRST BAND: 2700 cm-1 - 3200 cm-1, 3.5  MICRON BAND

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH3OH_TEMP,OM_BND_CH3OH(1,:),SD1_CH3OH)
   GDINV  = GC3*DINV_CH3OH
   GDDINV = GD*DINV_CH3OH

ELSE IF ((OM_BND_CH3OH(2,1)<=OMEGA).AND.(OMEGA<OM_BND_CH3OH(2,2))) THEN
! SECOND BAND: 3600 cm-1 - 3800 cm-1, 2.9  MICRON BAND

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH3OH_TEMP,OM_BND_CH3OH(2,:),SD2_CH3OH)
   GDINV  = GC3*DINV_CH3OH
   GDDINV = GD*DINV_CH3OH

ELSE IF ((OM_BND_CH3OH(3,1)<=OMEGA).AND.(OMEGA<OM_BND_CH3OH(3,2))) THEN
! THIRD BAND:  900 cm-1 - 1600 cm-1, 10.0 & 8.0 MICRON BANDS

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH3OH_TEMP,OM_BND_CH3OH(3,:),SD3_CH3OH)
   GDINV  = GC3*DINV_CH3OH
   GDDINV = GD*DINV_CH3OH

ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE CH3OH


!==============================================================================
SUBROUTINE C7H8(OMEGA,TEMP,PC7H8,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
! COMPUTE TOLUENE OPTICAL PROPERTIES USING SINGLE LINE GROUP
! 
! ABSORPTION BANDS: 
! THERE ARE 4 BANDS FOR Toluene
! BAND #1: 1000 cm-1 - 1150 cm-1 
! BAND #2: 1300 cm-1 - 1936 cm-1 
! BAND #3: 2700 cm-1 - 3200 cm-1 
! BAND #4: 700 cm-1 - 800 cm-1 
!
! VARIABLES PASSED IN
! OMEGA  : (REAL) WAVENUMBER IN CM-1
! TEMP   : (REAL) TEMPERATURE IN K
! PC7H8 : (REAL) PARTIAL PRESSURE OF TOLUENE
! PTOT   : (REAL) TOTAL PRESSURE
! GC3    : (REAL) collisional broadening half-width at half heigh
! 
! VRAIBLES PASSED OUT
! SDWEAK : (REAL) SPECTRAL ABSORPTION COEFFICIENT
! GDINV  : (REAL) LINE WIDTH TO LINE SPACING RATIO, FINE STRUCTURE PARAMETER
! GDDINV : (REAL) LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING,
!                 FINE STRUCTURE PARAMETER
!
! LOCALS
! PRESSURE_EFFECTIVE : (REAL) EFFECTIVE PRESSURE, (FORMELY PE)
! DINV_C7H8         : INVERSE LINE SPACING [CM] FOR TOLUENE

REAL(EB), INTENT(IN)  :: OMEGA, TEMP, PC7H8, PTOT, GC3
REAL(EB), INTENT(OUT) :: SDWEAK, GDINV, GDDINV

REAL(EB) :: Q2OT, AZOT, TOAZ
REAL(EB) :: GD, PRESSURE_EFFECTIVE
REAL(EB) :: DINV_C7H8 

REAL(EB), PARAMETER :: Q2       = 1.4388_EB  ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
REAL(EB), PARAMETER :: WM_C7H8  = 92.1384_EB ! NIST WEBBOOK DATA

! COMPUTE THERMAL COEFFICIENTS

Q2OT = -Q2/TEMP       
AZOT = 273._EB/TEMP
TOAZ = TEMP/273._EB

! COMPUTE DOPPLER BROADENING HALF WIDTH GD
GD   = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM_C7H8) !DOPPLER HALF WIDTH [CM^-1]. NASA,222

! COMPUTE AVERAGE LINE SPACING
! 
DINV_C7H8=1._EB

! SET INITIAL VALUES TO SDWEAK, GDINV, GDDINV - VALUES ARE RETURNED IF OMEGA IS NOT
! WITHIN ABSOPRTION BANDS
SDWEAK = 0.0_EB
GDINV  = 1._EB 
GDDINV = 1._EB 

! COMPUTE PRESSURE EFFECTIVE. 
PRESSURE_EFFECTIVE = PTOT+0.3_EB*PC7H8 ! (SPECIFIC TO METHANE)

! IF OMEGA IS WITHIN ABSOPRTION BAND, THEN COMPUTE SDWEAK FROM TABULATED DATA
! TOLUENE HAS FOUR BANDS
! BAND #1: 1000 cm-1 - 1150 cm-1 
! BAND #2: 1300 cm-1 - 1936 cm-1 
! BAND #3: 2700 cm-1 - 3200 cm-1 
! BAND #4: 700 cm-1 - 800 cm-1 

IF ((OM_BND_C7H8(1,1)<=OMEGA).AND.(OMEGA<OM_BND_C7H8(1,2))) THEN
! BAND #1: 1000 cm-1 - 1150 cm-1

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H8_TEMP,OM_BND_C7H8(1,:),SD1_C7H8)
   GDINV  = GC3*DINV_C7H8
   GDDINV = GD*DINV_C7H8

ELSE IF ((OM_BND_C7H8(2,1)<=OMEGA).AND.(OMEGA<OM_BND_C7H8(2,2))) THEN
! BAND #2: 1300 cm-1 - 1936 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H8_TEMP,OM_BND_C7H8(2,:),SD2_C7H8)
   GDINV  = GC3*DINV_C7H8
   GDDINV = GD*DINV_C7H8

ELSE IF ((OM_BND_C7H8(3,1)<=OMEGA).AND.(OMEGA<OM_BND_C7H8(3,2))) THEN
! BAND #3: 2700 cm-1 - 3200 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H8_TEMP,OM_BND_C7H8(3,:),SD3_C7H8)
   GDINV  = GC3*DINV_C7H8
   GDDINV = GD*DINV_C7H8

ELSE IF ((OM_BND_C7H8(4,1)<=OMEGA).AND.(OMEGA<OM_BND_C7H8(4,2))) THEN
! BAND #4: 700 cm-1 - 800 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H8_TEMP,OM_BND_C7H8(4,:),SD4_C7H8)
   GDINV  = GC3*DINV_C7H8
   GDDINV = GD*DINV_C7H8

ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE C7H8

!==============================================================================
SUBROUTINE C3H6(OMEGA,TEMP,PC3H6,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
! COMPUTE PROPYLENE OPTICAL PROPERTIES USING SINGLE LINE GROUP
! 
! ABSORPTION BANDS: 
! THERE ARE 3 BANDS FOR Propylene
!
! BAND #1: 1250 cm-1 - 1950 cm-1 
! BAND #2: 2700 cm-1 - 3200 cm-1 
! BAND #3: 700 cm-1 - 1150 cm-1 
!
! VARIABLES PASSED IN
! OMEGA  : (REAL) WAVENUMBER IN CM-1
! TEMP   : (REAL) TEMPERATURE IN K
! PC3H6 : (REAL) PARTIAL PRESSURE OF PROPYLENE
! PTOT   : (REAL) TOTAL PRESSURE
! GC3    : (REAL) collisional broadening half-width at half heigh
! 
! VRAIBLES PASSED OUT
! SDWEAK : (REAL) SPECTRAL ABSORPTION COEFFICIENT
! GDINV  : (REAL) LINE WIDTH TO LINE SPACING RATIO, FINE STRUCTURE PARAMETER
! GDDINV : (REAL) LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING,
!                 FINE STRUCTURE PARAMETER
!
! LOCALS
! PRESSURE_EFFECTIVE : (REAL) EFFECTIVE PRESSURE, (FORMELY PE)
! DINV_C3H6         : INVERSE LINE SPACING [CM] FOR PROPYLENE

REAL(EB), INTENT(IN)  :: OMEGA, TEMP, PC3H6, PTOT, GC3
REAL(EB), INTENT(OUT) :: SDWEAK, GDINV, GDDINV

REAL(EB) :: Q2OT, AZOT, TOAZ
REAL(EB) :: GD, PRESSURE_EFFECTIVE
REAL(EB) :: DINV_C3H6 

REAL(EB), PARAMETER :: Q2       = 1.4388_EB  ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
REAL(EB), PARAMETER :: WM_C3H6  = 42.0797_EB ! NIST WEBBOOK DATA

! COMPUTE THERMAL COEFFICIENTS

Q2OT = -Q2/TEMP       
AZOT = 273._EB/TEMP
TOAZ = TEMP/273._EB

! COMPUTE DOPPLER BROADENING HALF WIDTH GD
GD   = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM_C3H6) !DOPPLER HALF WIDTH [CM^-1]. NASA,222

! COMPUTE AVERAGE LINE SPACING
! 
DINV_C3H6 = 1._EB

! SET INITIAL VALUES TO SDWEAK, GDINV, GDDINV - VALUES ARE RETURNED IF OMEGA IS NOT
! WITHIN ABSOPRTION BANDS
SDWEAK = 0.0_EB
GDINV  = 1._EB 
GDDINV = 1._EB 

! COMPUTE PRESSURE EFFECTIVE. 
PRESSURE_EFFECTIVE = PTOT+0.3_EB*PC3H6 ! (SPECIFIC TO METHANE)

! IF OMEGA IS WITHIN ABSOPRTION BAND, THEN COMPUTE SDWEAK FROM TABULATED DATA
! PROPYLENE HAS THREE BANDS
! BAND #1: 1250 cm-1 - 1950 cm-1 
! BAND #2: 2700 cm-1 - 3200 cm-1 
! BAND #3: 700 cm-1 - 1150 cm-1 

IF ((OM_BND_C3H6(1,1)<=OMEGA).AND.(OMEGA<OM_BND_C3H6(1,2))) THEN
! BAND #1: 1250 cm-1 - 1950 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C3H6_TEMP,OM_BND_C3H6(1,:),SD1_C3H6)
   GDINV  = GC3*DINV_C3H6
   GDDINV = GD*DINV_C3H6

ELSE IF ((OM_BND_C3H6(2,1)<=OMEGA).AND.(OMEGA<OM_BND_C3H6(2,2))) THEN
! BAND #2: 2700 cm-1 - 3200 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C3H6_TEMP,OM_BND_C3H6(2,:),SD2_C3H6)
   GDINV  = GC3*DINV_C3H6
   GDDINV = GD*DINV_C3H6

ELSE IF ((OM_BND_C3H6(3,1)<=OMEGA).AND.(OMEGA<OM_BND_C3H6(3,2))) THEN
! BAND #3: 700 cm-1 - 1150 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C3H6_TEMP,OM_BND_C3H6(3,:),SD3_C3H6)
   GDINV  = GC3*DINV_C3H6
   GDDINV = GD*DINV_C3H6

ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE C3H6


!==============================================================================
SUBROUTINE MMA(OMEGA,TEMP,PMMA,PTOT,GC3,SDWEAK,GDINV,GDDINV)
!==============================================================================
! COMPUTE MMA OPTICAL PROPERTIES USING SINGLE LINE GROUP
! 
! ABSORPTION BANDS: 
! THERE ARE 2 BANDS FOR MMA
!
! BAND #1: 2700 cm-1 - 3200 cm-1 
! BAND #2: 750 cm-1 - 1950 cm-1 
!
! VARIABLES PASSED IN
! OMEGA  : (REAL) WAVENUMBER IN CM-1
! TEMP   : (REAL) TEMPERATURE IN K
! PMMA : (REAL) PARTIAL PRESSURE OF MMA
! PTOT   : (REAL) TOTAL PRESSURE
! GC3    : (REAL) collisional broadening half-width at half heigh
! 
! VRAIBLES PASSED OUT
! SDWEAK : (REAL) SPECTRAL ABSORPTION COEFFICIENT
! GDINV  : (REAL) LINE WIDTH TO LINE SPACING RATIO, FINE STRUCTURE PARAMETER
! GDDINV : (REAL) LINE WIDTH TO LINE SPACING RATIO FOR DOPPLER BORADENING,
!                 FINE STRUCTURE PARAMETER
!
! LOCALS
! PRESSURE_EFFECTIVE : (REAL) EFFECTIVE PRESSURE, (FORMELY PE)
! DINV_MMA         : INVERSE LINE SPACING [CM] FOR MMA

REAL(EB), INTENT(IN)  :: OMEGA, TEMP, PMMA, PTOT, GC3
REAL(EB), INTENT(OUT) :: SDWEAK, GDINV, GDDINV

REAL(EB) :: Q2OT, AZOT, TOAZ
REAL(EB) :: GD, PRESSURE_EFFECTIVE
REAL(EB) :: DINV_MMA 

REAL(EB), PARAMETER :: Q2       = 1.4388_EB   ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
REAL(EB), PARAMETER :: WM_MMA   = 100.1158_EB ! NIST WEBBOOK DATA

! COMPUTE THERMAL COEFFICIENTS

Q2OT = -Q2/TEMP       
AZOT = 273._EB/TEMP
TOAZ = TEMP/273._EB

! COMPUTE DOPPLER BROADENING HALF WIDTH GD
GD   = 5.94E-6_EB*OMEGA*SQRT(TOAZ/WM_MMA) !DOPPLER HALF WIDTH [CM^-1]. NASA,222

! COMPUTE AVERAGE LINE SPACING
! 
DINV_MMA = 1._EB

! SET INITIAL VALUES TO SDWEAK, GDINV, GDDINV - VALUES ARE RETURNED IF OMEGA IS NOT
! WITHIN ABSOPRTION BANDS
SDWEAK = 0.0_EB
GDINV  = 1._EB 
GDDINV = 1._EB 

! COMPUTE PRESSURE EFFECTIVE. 
PRESSURE_EFFECTIVE = PTOT+0.3_EB*PMMA ! (SPECIFIC TO METHANE)

! IF OMEGA IS WITHIN ABSOPRTION BAND, THEN COMPUTE SDWEAK FROM TABULATED DATA
! MMA HAS TWO BANDS
! BAND #1: 2700 cm-1 - 3200 cm-1 
! BAND #2: 750 cm-1 - 1950 cm-1  

IF ((OM_BND_MMA(1,1)<=OMEGA).AND.(OMEGA<OM_BND_MMA(1,2))) THEN
! BAND #1: 2700 cm-1 - 3200 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_MMA_TEMP,OM_BND_MMA(1,:),SD1_MMA)
   GDINV  = GC3*DINV_MMA
   GDDINV = GD*DINV_MMA

ELSE IF ((OM_BND_MMA(2,1)<=OMEGA).AND.(OMEGA<OM_BND_MMA(2,2))) THEN
! BAND #2: 750 cm-1 - 1950 cm-1 

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_MMA_TEMP,OM_BND_MMA(2,:),SD2_MMA)
   GDINV  = GC3*DINV_MMA
   GDDINV = GD*DINV_MMA

ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE MMA


!==============================================================================
REAL(EB) FUNCTION GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_TEMP,BOUNDS,ABSORPTION_DATA)
!==============================================================================
! RETURNS INTERPOLATED VALUES OF SPECTRAL ABSORPTION
! NOTE ASSUMED SHAPE
!
! Input:
! OMEGA  : WAVENUMBER, UNIT = 1/cm
! TEMP   : GAS TEMPERATURE, UNIT = KELVIN
! SD_TEMP: VECTOR OF TEMPERATURE MEASUREMENTS, DIMENSION(*) UNIT = KELVIN
! BOUNDS : VECTOR OF WAVENUMBER BOUNDS AND SPACING FOR A GIVEN BAND, DIMENSIONS(3)
! ABSORPTION_DATA: ARRAY OF SPECIES SPECTRAL MEASUREMENTS, DIMENSIONS(*,*)
!                 NOTE = FIRST DIMENSION EQUALS TO DIMENSION OF SD_TEMP
!
! ACTION PERFORMED: BILINEAR INTERPOLATION OF ABSORPTION_DATA
!------------------------------------------------------------------------------
! Variables passed in

 REAL(EB), INTENT(IN) :: OMEGA
 REAL(EB), INTENT(IN) :: TEMP
 REAL(EB), DIMENSION(:), INTENT(IN)   :: SD_TEMP, BOUNDS
 REAL(EB), DIMENSION(:,:), INTENT(IN) :: ABSORPTION_DATA 

 REAL(EB) :: DELTA_T
 REAL(EB) :: DELTA_W
 REAL(EB) :: W1
 REAL(EB) :: OMEGA_MIN, OMEGA_MAX, DELTA_OMEGA
 REAL(EB) :: TTEMP
 REAL(EB) :: INTERPOLATED_VALUE

 INTEGER :: N_TEMP, N_BOUNDS, N_ABS_DATA1, N_ABS_DATA2
 INTEGER :: I_OMEGA, I_TEMP

! CVL NEED TO LINK MODULE THAT DEFINES ERRMSG
! PART OF THE CODE COMMENTED FOR NOW
!  CHARACTER(255) MESSAG
!  LOGICAL FATAL

 N_TEMP      = UBOUND(SD_TEMP,1)  ! GET NUMBER OF TEMPERATURE MEASUREMENTS
 N_BOUNDS    = UBOUND(BOUNDS,1)   ! GET NUMBER OF ELEMENTS IN BOUNDS. SHOULD BE 3
 N_ABS_DATA1 = UBOUND(ABSORPTION_DATA,1) ! GET NUMBER OF TEMPERATURE DATA
 N_ABS_DATA2 = UBOUND(ABSORPTION_DATA,2) ! GET NUMBER OF WAVENUMBER DATA


! ! CVL NEED TO LINK MODULE THAT DEFINES ERRMSG
! ! PART OF THE CODE COMMENTED FOR NOW
! !-----------------------------TEST DIMENSIONS OF BOUNDS AND ABSORPTION_DATA----
! ! N_BOUNDS SHOULD BE 3
! 
!  IF (N_BOUNDS/=3) THEN 
!     MESSAG='ERROR IN RADCAL. N_BOUNDS SHOULD BE EQUAL TO 3'
!     FATAL = .TRUE.
!     CALL ERRMSG( trim(MESSAG), FATAL )
!  ENDIF
! 
! ! N_ABS_DATA1 SHOULD BE EQUAL TO N_TEMP
! 
!  IF (N_ABS_DATA1/=N_TEMP) THEN 
!     MESSAG='ERROR IN RADCAL. FIRST DIMENSION OF ABSORPTION_DATA SHOULD BE EQUAL TO N_TEMP'
!     FATAL = .TRUE.
!     CALL ERRMSG( trim(MESSAG), FATAL )
!  ENDIF

!------------------------------------------------------------------------------
 OMEGA_MIN   = BOUNDS(1)
 OMEGA_MAX   = BOUNDS(2)
 DELTA_OMEGA = BOUNDS(3)
 
! BOUNDS TTEMP SO THAT MIN(SD_TEMP) <= TTEMP <= MAX(SD_TEMP)
! IF TTEMP < MIN(SD_TEMP) THEN TTEMP := MIN(SD_TEMP)
! IF MAX(SD_TEMP) < TTEMP THEN TTEMP := MAX(SD_TEMP)
 TTEMP       = MIN(MAXVAL(SD_TEMP),MAX(TEMP,MINVAL(SD_TEMP)))  

 IF (TTEMP <= MINVAL(SD_TEMP)) THEN
    I_TEMP  = 1
    DELTA_T = 0._EB
 ELSE IF (TTEMP >= MAXVAL(SD_TEMP)) THEN
    I_TEMP  = N_TEMP-1
    DELTA_T = 1.0_EB
 ELSE
! FIND INDEX I_TEMP SUCH THAT SD_TEMP(I_TEMP)<TTEMP<=SD_TEMP(I_TEMP+1)
    DO I_TEMP = 1, N_TEMP-1
       IF ((SD_TEMP(I_TEMP)<TTEMP).AND.(TTEMP<=SD_TEMP(I_TEMP+1))) EXIT
    ENDDO

    DELTA_T = (TTEMP-SD_TEMP(I_TEMP))/(SD_TEMP(I_TEMP+1)-SD_TEMP(I_TEMP))
 ENDIF

 IF ((OMEGA_MIN<=OMEGA).AND.(OMEGA<=OMEGA_MAX)) THEN 

    IF(OMEGA<OMEGA_MAX) THEN
! FIND I_OMEGA SUCH THAT OMEGA(I_OMEGA)<=OMEGA<OMEGA(I_OMEGA+1)
       I_OMEGA = FLOOR((OMEGA-OMEGA_MIN)/DELTA_OMEGA)+1
       I_OMEGA = MIN(I_OMEGA,N_ABS_DATA2-1)
       W1      = OMEGA_MIN+DELTA_OMEGA*(REAL(I_OMEGA-1,EB))
       DELTA_W = (OMEGA-W1)/DELTA_OMEGA
    ELSE
! CASE I_OMEGA = OMEGA_MAX => ASSUME THEN THAT I_OMEGA = N_ABS_DATA2 -1 
       I_OMEGA = N_ABS_DATA2 - 1
       W1      = OMEGA_MIN+DELTA_OMEGA*(REAL(I_OMEGA-1,EB))
       DELTA_W = (OMEGA-W1)/DELTA_OMEGA
    ENDIF

! PERFORM BILINEAR INTERPOLATION TO OBTAIN VALUES OF SDWEAK
    INTERPOLATED_VALUE = (1._EB-DELTA_W)*((1._EB-DELTA_T)*ABSORPTION_DATA(I_TEMP,I_OMEGA)+ & 
                         ABSORPTION_DATA(I_TEMP+1,I_OMEGA)*DELTA_T) +                      &
                         DELTA_W*(ABSORPTION_DATA(I_TEMP,I_OMEGA+1)*(1._EB-DELTA_T)+       &
                         ABSORPTION_DATA(I_TEMP+1,I_OMEGA+1)*DELTA_T)

 ELSE
    INTERPOLATED_VALUE = 0._EB
 ENDIF

! ENSURE POSITIVITY
 GET_SPECTRAL_ABSORPTION = MAX(0._EB,INTERPOLATED_VALUE)
!------------------------------------------------------------------------------
 END FUNCTION GET_SPECTRAL_ABSORPTION


!==============================================================================
REAL(EB) FUNCTION PLANCK(TEMP,LAMBDA)
!==============================================================================
! COMPUTES BLACKBODY FUNCTION IN UNITS OF W/M-2/MICRON/SR
! Input:
! TEMP: TEMPERATURE, UNIT = KELVIN
! LAMBDA: WAVELENGTH, UNIT = MicroMETERS
!------------------------------------------------------------------------------
! Variables passed in
 REAL (EB), INTENT(IN):: TEMP, LAMBDA

! local variables
 REAL (EB), PARAMETER :: Q1 = 1.19088E8_EB ! Q1 = 2*SPEED_OF_LIGHT^2*PLANCK_CNST
 REAL (EB), PARAMETER :: Q2 = 14388._EB    ! Q2 = SPEED_OF_LIGHT*PLANCK_CNS/BOLTZMANN
 REAL (EB) :: C   ! C = LAMBDA * TEMP

 IF(ABS(TEMP)<ZERO_P) THEN
    PLANCK = 0._EB
 ELSE
    C = TEMP * LAMBDA
    IF (Q2/C > 38._EB) THEN
       PLANCK = 0._EB
    ELSE
       PLANCK=Q1*(LAMBDA**(-5))/(EXP(Q2/C)-1._EB)
   END IF
 ENDIF
!------------------------------------------------------------------------------
END FUNCTION PLANCK

!==============================================================================
SUBROUTINE RCALLOC
!==============================================================================

ALLOCATE(P(6))
ALLOCATE(SPECIE(5))
ALLOCATE(QW(600))
ALLOCATE(TTAU(600))
ALLOCATE(XT(600))
ALLOCATE(X(4))
ALLOCATE(GC(4))
ALLOCATE(AMBDA(600))
ALLOCATE(UUU(4))
ALLOCATE(AB(600))
ALLOCATE(ATOT(4))
ALLOCATE(BCNT(4))
ALLOCATE(GAMMA(4,7))
ALLOCATE(SD15(6,80))
ALLOCATE(SD(6,376))
ALLOCATE(SD7(3,16))
ALLOCATE(SD3(3,32))

!-------------------------Methane DATA-------------------


! THERE ARE 2 BANDS FOR Methane

! BAND #1: 1150 cm-1 - 1600 cm-1 
! BAND #2: 2700 cm-1 - 3250 cm-1 
! BAND #3: 3400 cm-1 - 5000 cm-1

N_TEMP_CH4 = 23
N_BAND_CH4 = 3

ALLOCATE(SD_CH4_TEMP(N_TEMP_CH4)) 

! INITIALIZE BANDS WAVENUMBER BOUNDS FOR Methane ABSOPTION COEFFICIENTS.
! BANDS ARE ORGANIZED BY ROW: 
! 1ST COLUMN IS THE LOWER BOUND BAND LIMIT IN cm-1. 
! 2ND COLUMN IS THE UPPER BOUND BAND LIMIT IN cm-1. 
! 3RD COLUMN IS THE STRIDE BETWEEN WAVENUMBERS IN BAND.
! IF 3RD COLUMN = 0., THEN THE BAND IS CALCULATED AND NOT TABULATED.

ALLOCATE(OM_BND_CH4(N_BAND_CH4,3)) 

OM_BND_CH4 = RESHAPE((/ &
    1150._EB, 2700._EB, 3400._EB,&
    1600._EB, 3250._EB, 5000._EB,&
      25._EB,   25._EB,    0._EB/),(/N_BAND_CH4,3/)) 

SD_CH4_TEMP = (/ &
    300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
    600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
    900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
    1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(SD1_CH4(N_TEMP_CH4,19)) 

! BAND #1: 1150 cm-1 - 1600 cm-1 

SD1_CH4(1:23,1:8) = RESHAPE((/ &  ! 1150-1325 cm-1
    2.674955E-03_EB, 4.122616E-03_EB, 6.116504E-03_EB, 8.394889E-03_EB, 1.066588E-02_EB, 1.270929E-02_EB, &
    1.431368E-02_EB, 1.543500E-02_EB, 1.607798E-02_EB, 1.623394E-02_EB, 1.594103E-02_EB, 1.537005E-02_EB, &
    1.456109E-02_EB, 1.356591E-02_EB, 1.253062E-02_EB, 1.145338E-02_EB, 1.037156E-02_EB, 9.319748E-03_EB, &
    8.296535E-03_EB, 7.359273E-03_EB, 6.512254E-03_EB, 5.764075E-03_EB, 5.034367E-03_EB,  &
    1.056047E-02_EB, 1.689522E-02_EB, 2.382293E-02_EB, 3.013931E-02_EB, 3.518013E-02_EB, 3.873081E-02_EB, &
    4.056505E-02_EB, 4.099451E-02_EB, 4.031689E-02_EB, 3.867720E-02_EB, 3.630281E-02_EB, 3.360423E-02_EB, &
    3.072697E-02_EB, 2.769650E-02_EB, 2.480529E-02_EB, 2.203561E-02_EB, 1.945492E-02_EB, 1.710571E-02_EB, &
    1.490918E-02_EB, 1.297063E-02_EB, 1.127359E-02_EB, 9.806206E-03_EB, 8.433200E-03_EB,  &
    6.535368E-02_EB, 8.634735E-02_EB, 1.035285E-01_EB, 1.149338E-01_EB, 1.211259E-01_EB, 1.231160E-01_EB, &
    1.212403E-01_EB, 1.167139E-01_EB, 1.104665E-01_EB, 1.027695E-01_EB, 9.409726E-02_EB, 8.537322E-02_EB, &
    7.679556E-02_EB, 6.836183E-02_EB, 6.060003E-02_EB, 5.339354E-02_EB, 4.679438E-02_EB, 4.089038E-02_EB, &
    3.545789E-02_EB, 3.072468E-02_EB, 2.660506E-02_EB, 2.306236E-02_EB, 1.979016E-02_EB,  &
    2.797310E-01_EB, 2.878774E-01_EB, 2.844637E-01_EB, 2.718418E-01_EB, 2.544158E-01_EB, 2.351542E-01_EB, &
    2.144282E-01_EB, 1.940107E-01_EB, 1.744123E-01_EB, 1.554283E-01_EB, 1.374171E-01_EB, 1.210226E-01_EB, &
    1.060935E-01_EB, 9.238361E-02_EB, 8.033186E-02_EB, 6.960028E-02_EB, 6.015950E-02_EB, 5.186871E-02_EB, &
    4.448682E-02_EB, 3.815557E-02_EB, 3.272683E-02_EB, 2.816055E-02_EB, 2.399102E-02_EB,  &
    6.055294E-01_EB, 5.085203E-01_EB, 4.314546E-01_EB, 3.661128E-01_EB, 3.120083E-01_EB, 2.675996E-01_EB, &
    2.295433E-01_EB, 1.976798E-01_EB, 1.706085E-01_EB, 1.470613E-01_EB, 1.263463E-01_EB, 1.086762E-01_EB, &
    9.339644E-02_EB, 7.991483E-02_EB, 6.850499E-02_EB, 5.859531E-02_EB, 5.007218E-02_EB, 4.275497E-02_EB, &
    3.636644E-02_EB, 3.094626E-02_EB, 2.637845E-02_EB, 2.255529E-02_EB, 1.911152E-02_EB,  &
    7.291248E-01_EB, 6.247580E-01_EB, 5.567927E-01_EB, 5.047634E-01_EB, 4.622435E-01_EB, 4.253591E-01_EB, &
    3.901550E-01_EB, 3.566718E-01_EB, 3.248608E-01_EB, 2.936112E-01_EB, 2.628808E-01_EB, 2.344160E-01_EB, &
    2.078641E-01_EB, 1.829068E-01_EB, 1.606407E-01_EB, 1.403181E-01_EB, 1.221762E-01_EB, 1.060864E-01_EB, &
    9.154534E-02_EB, 7.892764E-02_EB, 6.808643E-02_EB, 5.882322E-02_EB, 5.031768E-02_EB,  &
    1.503577E+00_EB, 1.216438E+00_EB, 1.009787E+00_EB, 8.445781E-01_EB, 7.110515E-01_EB, 6.021463E-01_EB, &
    5.099874E-01_EB, 4.327438E-01_EB, 3.680691E-01_EB, 3.125303E-01_EB, 2.644163E-01_EB, 2.242054E-01_EB, &
    1.897708E-01_EB, 1.602130E-01_EB, 1.354708E-01_EB, 1.143770E-01_EB, 9.654761E-02_EB, 8.152145E-02_EB, &
    6.854917E-02_EB, 5.774320E-02_EB, 4.875464E-02_EB, 4.129146E-02_EB, 3.467888E-02_EB,  &
    1.113787E+00_EB, 8.666106E-01_EB, 6.937600E-01_EB, 5.636603E-01_EB, 4.651973E-01_EB, 3.899245E-01_EB, &
    3.293625E-01_EB, 2.804438E-01_EB, 2.403481E-01_EB, 2.062563E-01_EB, 1.767199E-01_EB, 1.518103E-01_EB, &
    1.304026E-01_EB, 1.116263E-01_EB, 9.571162E-02_EB, 8.192030E-02_EB, 7.003835E-02_EB, 5.986022E-02_EB, &
    5.092725E-02_EB, 4.337782E-02_EB, 3.700459E-02_EB, 3.164349E-02_EB, 2.682395E-02_EB/),(/23,8/))

SD1_CH4(1:23,9:16) = RESHAPE((/ &  ! 1350-1525 cm-1
    5.859637E-01_EB, 5.741640E-01_EB, 5.446284E-01_EB, 5.016955E-01_EB, 4.533548E-01_EB, 4.052375E-01_EB, &
    3.578427E-01_EB, 3.138042E-01_EB, 2.741405E-01_EB, 2.379703E-01_EB, 2.050667E-01_EB, 1.765003E-01_EB, &
    1.515652E-01_EB, 1.294551E-01_EB, 1.106049E-01_EB, 9.430003E-02_EB, 8.028460E-02_EB, 6.831086E-02_EB, &
    5.788843E-02_EB, 4.905327E-02_EB, 4.167658E-02_EB, 3.549325E-02_EB, 2.997796E-02_EB,  &
    5.143640E-02_EB, 7.394302E-02_EB, 9.480931E-02_EB, 1.108135E-01_EB, 1.211319E-01_EB, 1.262273E-01_EB, &
    1.262909E-01_EB, 1.227527E-01_EB, 1.167145E-01_EB, 1.087709E-01_EB, 9.943737E-02_EB, 9.002538E-02_EB, &
    8.060527E-02_EB, 7.137067E-02_EB, 6.299760E-02_EB, 5.512164E-02_EB, 4.801455E-02_EB, 4.169764E-02_EB, &
    3.593587E-02_EB, 3.095498E-02_EB, 2.664278E-02_EB, 2.295572E-02_EB, 1.959682E-02_EB,  &
    7.616134E-03_EB, 9.957002E-03_EB, 1.281279E-02_EB, 1.589717E-02_EB, 1.894916E-02_EB, 2.174268E-02_EB, &
    2.395920E-02_EB, 2.555557E-02_EB, 2.648904E-02_EB, 2.671274E-02_EB, 2.629477E-02_EB, 2.540283E-02_EB, &
    2.418621E-02_EB, 2.265546E-02_EB, 2.101542E-02_EB, 1.928221E-02_EB, 1.752831E-02_EB, 1.582992E-02_EB, &
    1.416926E-02_EB, 1.261250E-02_EB, 1.121041E-02_EB, 9.951538E-03_EB, 8.728757E-03_EB,  &
    1.277402E-02_EB, 1.327177E-02_EB, 1.322108E-02_EB, 1.267471E-02_EB, 1.186190E-02_EB, 1.093562E-02_EB, &
    9.920972E-03_EB, 8.920419E-03_EB, 7.971129E-03_EB, 7.057028E-03_EB, 6.197372E-03_EB, 5.429023E-03_EB, &
    4.732635E-03_EB, 4.092017E-03_EB, 3.548360E-03_EB, 3.062048E-03_EB, 2.633809E-03_EB, 2.263255E-03_EB, &
    1.935254E-03_EB, 1.656795E-03_EB, 1.416023E-03_EB, 1.215289E-03_EB, 1.031774E-03_EB,  &
    8.014464E-03_EB, 7.194752E-03_EB, 6.553713E-03_EB, 6.000398E-03_EB, 5.531684E-03_EB, 5.121393E-03_EB, &
    4.722456E-03_EB, 4.343916E-03_EB, 3.983544E-03_EB, 3.621247E-03_EB, 3.268422E-03_EB, 2.936068E-03_EB, &
    2.617151E-03_EB, 2.322363E-03_EB, 2.051315E-03_EB, 1.804840E-03_EB, 1.578471E-03_EB, 1.378963E-03_EB, &
    1.197596E-03_EB, 1.036958E-03_EB, 9.007021E-04_EB, 7.820103E-04_EB, 6.693842E-04_EB,  &
    2.792665E-03_EB, 2.575454E-03_EB, 2.537125E-03_EB, 2.578233E-03_EB, 2.648059E-03_EB, 2.713453E-03_EB, &
    2.741224E-03_EB, 2.722787E-03_EB, 2.662395E-03_EB, 2.559402E-03_EB, 2.415338E-03_EB, 2.254521E-03_EB, &
    2.077537E-03_EB, 1.891657E-03_EB, 1.708933E-03_EB, 1.533948E-03_EB, 1.367874E-03_EB, 1.212162E-03_EB, &
    1.064882E-03_EB, 9.339250E-04_EB, 8.169687E-04_EB, 7.180516E-04_EB, 6.220503E-04_EB,  &
    8.239408E-04_EB, 8.207289E-04_EB, 8.769232E-04_EB, 9.471812E-04_EB, 1.013068E-03_EB, 1.072041E-03_EB, &
    1.109555E-03_EB, 1.125480E-03_EB, 1.120373E-03_EB, 1.094883E-03_EB, 1.043855E-03_EB, 9.887558E-04_EB, &
    9.221581E-04_EB, 8.517148E-04_EB, 7.786200E-04_EB, 7.086512E-04_EB, 6.363394E-04_EB, 5.709486E-04_EB, &
    5.037963E-04_EB, 4.478632E-04_EB, 3.933607E-04_EB, 3.461225E-04_EB, 3.036065E-04_EB,  &
    1.787674E-02_EB, 1.651301E-02_EB, 1.538618E-02_EB, 1.427021E-02_EB, 1.321135E-02_EB, 1.220556E-02_EB, &
    1.120017E-02_EB, 1.022271E-02_EB, 9.293408E-03_EB, 8.381487E-03_EB, 7.491250E-03_EB, 6.673691E-03_EB, &
    5.914957E-03_EB, 5.201979E-03_EB, 4.566356E-03_EB, 3.988961E-03_EB, 3.472458E-03_EB, 3.017625E-03_EB, &
    2.605493E-03_EB, 2.248445E-03_EB, 1.941774E-03_EB, 1.676885E-03_EB, 1.438512E-03_EB/),(/23,8/))

SD1_CH4(1:23,17:19) = RESHAPE((/ &  ! 1550-1600 cm-1
    3.842192E-03_EB, 4.909077E-03_EB, 5.876585E-03_EB, 6.637230E-03_EB, 7.159498E-03_EB, 7.449619E-03_EB, &
    7.502474E-03_EB, 7.374704E-03_EB, 7.105819E-03_EB, 6.720456E-03_EB, 6.241450E-03_EB, 5.729337E-03_EB, &
    5.213056E-03_EB, 4.685851E-03_EB, 4.183214E-03_EB, 3.714659E-03_EB, 3.276286E-03_EB, 2.879051E-03_EB, &
    2.508044E-03_EB, 2.182116E-03_EB, 1.899894E-03_EB, 1.653926E-03_EB, 1.423668E-03_EB,  &
    1.587677E-03_EB, 1.446538E-03_EB, 1.552901E-03_EB, 1.804439E-03_EB, 2.128742E-03_EB, 2.458565E-03_EB, &
    2.743898E-03_EB, 2.955974E-03_EB, 3.088835E-03_EB, 3.134690E-03_EB, 3.095698E-03_EB, 3.002805E-03_EB, &
    2.862955E-03_EB, 2.683299E-03_EB, 2.486550E-03_EB, 2.282038E-03_EB, 2.073847E-03_EB, 1.873143E-03_EB, &
    1.674361E-03_EB, 1.489429E-03_EB, 1.321804E-03_EB, 1.173554E-03_EB, 1.026702E-03_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,3/))

ALLOCATE(SD2_CH4(N_TEMP_CH4,23)) 

! BAND #2: 2700 cm-1 - 3250 cm-1 

SD2_CH4(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
    1.217530E-02_EB, 1.117112E-02_EB, 1.018335E-02_EB, 9.113159E-03_EB, 8.069237E-03_EB, 7.089143E-03_EB, &
    6.164984E-03_EB, 5.339743E-03_EB, 4.605161E-03_EB, 3.952813E-03_EB, 3.369126E-03_EB, 2.873535E-03_EB, &
    2.442004E-03_EB, 2.069306E-03_EB, 1.752468E-03_EB, 1.483929E-03_EB, 1.255931E-03_EB, 1.056565E-03_EB, &
    8.920522E-04_EB, 7.508758E-04_EB, 6.317838E-04_EB, 5.352971E-04_EB, 4.530464E-04_EB,  &
    2.513916E-02_EB, 2.114376E-02_EB, 1.780513E-02_EB, 1.491293E-02_EB, 1.247243E-02_EB, 1.043356E-02_EB, &
    8.698933E-03_EB, 7.256269E-03_EB, 6.066735E-03_EB, 5.057876E-03_EB, 4.209670E-03_EB, 3.508303E-03_EB, &
    2.922322E-03_EB, 2.428483E-03_EB, 2.028794E-03_EB, 1.689558E-03_EB, 1.409261E-03_EB, 1.176715E-03_EB, &
    9.788265E-04_EB, 8.167162E-04_EB, 6.823326E-04_EB, 5.730157E-04_EB, 4.786727E-04_EB,  &
    2.557782E-02_EB, 2.133345E-02_EB, 1.832454E-02_EB, 1.586992E-02_EB, 1.378920E-02_EB, 1.202910E-02_EB, &
    1.045539E-02_EB, 9.073223E-03_EB, 7.878108E-03_EB, 6.813401E-03_EB, 5.858711E-03_EB, 5.036798E-03_EB, &
    4.322566E-03_EB, 3.689893E-03_EB, 3.157259E-03_EB, 2.690401E-03_EB, 2.292743E-03_EB, 1.951361E-03_EB, &
    1.658405E-03_EB, 1.405787E-03_EB, 1.195120E-03_EB, 1.018165E-03_EB, 8.624014E-04_EB,  &
    3.267234E-02_EB, 2.794974E-02_EB, 2.414894E-02_EB, 2.075791E-02_EB, 1.780983E-02_EB, 1.531517E-02_EB, &
    1.312860E-02_EB, 1.127375E-02_EB, 9.704346E-03_EB, 8.352011E-03_EB, 7.159191E-03_EB, 6.165400E-03_EB, &
    5.304681E-03_EB, 4.552437E-03_EB, 3.918329E-03_EB, 3.371374E-03_EB, 2.897155E-03_EB, 2.490584E-03_EB, &
    2.139129E-03_EB, 1.834560E-03_EB, 1.579448E-03_EB, 1.362473E-03_EB, 1.164728E-03_EB,  &
    4.729257E-02_EB, 4.020146E-02_EB, 3.655117E-02_EB, 3.486375E-02_EB, 3.451183E-02_EB, 3.486976E-02_EB, &
    3.526911E-02_EB, 3.549943E-02_EB, 3.539419E-02_EB, 3.477334E-02_EB, 3.359574E-02_EB, 3.209631E-02_EB, &
    3.033091E-02_EB, 2.826035E-02_EB, 2.617425E-02_EB, 2.398793E-02_EB, 2.182750E-02_EB, 1.975692E-02_EB, &
    1.770980E-02_EB, 1.581216E-02_EB, 1.409330E-02_EB, 1.256203E-02_EB, 1.106706E-02_EB,  &
    3.967221E-02_EB, 3.764359E-02_EB, 3.949747E-02_EB, 4.277841E-02_EB, 4.623172E-02_EB, 4.903094E-02_EB, &
    5.058120E-02_EB, 5.097742E-02_EB, 5.029388E-02_EB, 4.863988E-02_EB, 4.605751E-02_EB, 4.308689E-02_EB, &
    3.983596E-02_EB, 3.638363E-02_EB, 3.298206E-02_EB, 2.967115E-02_EB, 2.650179E-02_EB, 2.356866E-02_EB, &
    2.077236E-02_EB, 1.827032E-02_EB, 1.605444E-02_EB, 1.411180E-02_EB, 1.226114E-02_EB,  &
    7.223868E-02_EB, 8.544183E-02_EB, 9.962661E-02_EB, 1.104863E-01_EB, 1.171420E-01_EB, 1.196787E-01_EB, &
    1.182118E-01_EB, 1.139658E-01_EB, 1.078438E-01_EB, 1.002471E-01_EB, 9.161471E-02_EB, 8.298437E-02_EB, &
    7.452297E-02_EB, 6.622584E-02_EB, 5.857678E-02_EB, 5.147962E-02_EB, 4.505849E-02_EB, 3.932179E-02_EB, &
    3.408714E-02_EB, 2.948172E-02_EB, 2.551286E-02_EB, 2.211846E-02_EB, 1.897647E-02_EB,  &
    1.759866E-01_EB, 2.015363E-01_EB, 2.180490E-01_EB, 2.232613E-01_EB, 2.203020E-01_EB, 2.118597E-01_EB, &
    1.990075E-01_EB, 1.841205E-01_EB, 1.683887E-01_EB, 1.521963E-01_EB, 1.360646E-01_EB, 1.209307E-01_EB, &
    1.069128E-01_EB, 9.382435E-02_EB, 8.206146E-02_EB, 7.155641E-02_EB, 6.219310E-02_EB, 5.392479E-02_EB, &
    4.648837E-02_EB, 4.008493E-02_EB, 3.458414E-02_EB, 2.988156E-02_EB, 2.558176E-02_EB/),(/23,8/))

SD2_CH4(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
    3.311709E-01_EB, 3.344762E-01_EB, 3.267663E-01_EB, 3.099531E-01_EB, 2.886161E-01_EB, 2.659881E-01_EB, &
    2.423456E-01_EB, 2.191913E-01_EB, 1.973904E-01_EB, 1.764196E-01_EB, 1.564157E-01_EB, 1.382986E-01_EB, &
    1.218590E-01_EB, 1.067421E-01_EB, 9.331013E-02_EB, 8.132822E-02_EB, 7.068471E-02_EB, 6.139407E-02_EB, &
    5.297295E-02_EB, 4.571219E-02_EB, 3.948284E-02_EB, 3.416733E-02_EB, 2.931231E-02_EB,  &
    8.235613E-01_EB, 7.238968E-01_EB, 6.324968E-01_EB, 5.466365E-01_EB, 4.712117E-01_EB, 4.058971E-01_EB, &
    3.486088E-01_EB, 2.993153E-01_EB, 2.575278E-01_EB, 2.209771E-01_EB, 1.890168E-01_EB, 1.618217E-01_EB, &
    1.384064E-01_EB, 1.180022E-01_EB, 1.008139E-01_EB, 8.601549E-02_EB, 7.332027E-02_EB, 6.253544E-02_EB, &
    5.308446E-02_EB, 4.515135E-02_EB, 3.847721E-02_EB, 3.288247E-02_EB, 2.786602E-02_EB,  &
    6.821008E-01_EB, 5.401363E-01_EB, 4.353503E-01_EB, 3.533680E-01_EB, 2.892836E-01_EB, 2.389427E-01_EB, &
    1.981961E-01_EB, 1.653239E-01_EB, 1.387030E-01_EB, 1.163771E-01_EB, 9.777404E-02_EB, 8.236369E-02_EB, &
    6.949908E-02_EB, 5.853014E-02_EB, 4.942064E-02_EB, 4.176497E-02_EB, 3.531768E-02_EB, 2.985698E-02_EB, &
    2.516886E-02_EB, 2.128060E-02_EB, 1.802439E-02_EB, 1.532867E-02_EB, 1.292454E-02_EB,  &
    5.823907E-01_EB, 4.331161E-01_EB, 3.353248E-01_EB, 2.655224E-01_EB, 2.144294E-01_EB, 1.761169E-01_EB, &
    1.472744E-01_EB, 1.241615E-01_EB, 1.057407E-01_EB, 9.045740E-02_EB, 7.745806E-02_EB, 6.667193E-02_EB, &
    5.744892E-02_EB, 4.942743E-02_EB, 4.263289E-02_EB, 3.676096E-02_EB, 3.166692E-02_EB, 2.731256E-02_EB, &
    2.343114E-02_EB, 2.014747E-02_EB, 1.734533E-02_EB, 1.497211E-02_EB, 1.281431E-02_EB,  &
    3.210821E+00_EB, 2.793045E+00_EB, 2.475891E+00_EB, 2.198508E+00_EB, 1.954455E+00_EB, 1.740419E+00_EB, &
    1.543274E+00_EB, 1.366589E+00_EB, 1.208668E+00_EB, 1.063917E+00_EB, 9.309610E-01_EB, 8.135837E-01_EB, &
    7.091172E-01_EB, 6.150707E-01_EB, 5.334009E-01_EB, 4.611724E-01_EB, 3.980142E-01_EB, 3.431784E-01_EB, &
    2.944666E-01_EB, 2.526486E-01_EB, 2.171496E-01_EB, 1.870189E-01_EB, 1.596229E-01_EB,  &
    3.517151E-01_EB, 2.541996E-01_EB, 1.936159E-01_EB, 1.517283E-01_EB, 1.221750E-01_EB, 1.004308E-01_EB, &
    8.351777E-02_EB, 7.016371E-02_EB, 5.949379E-02_EB, 5.053966E-02_EB, 4.298934E-02_EB, 3.667062E-02_EB, &
    3.130112E-02_EB, 2.669085E-02_EB, 2.279490E-02_EB, 1.946235E-02_EB, 1.660600E-02_EB, 1.417075E-02_EB, &
    1.204697E-02_EB, 1.025668E-02_EB, 8.756988E-03_EB, 7.490243E-03_EB, 6.351877E-03_EB,  &
    7.921358E-01_EB, 6.055366E-01_EB, 4.798618E-01_EB, 3.867414E-01_EB, 3.174199E-01_EB, 2.642038E-01_EB, &
    2.216893E-01_EB, 1.874171E-01_EB, 1.594862E-01_EB, 1.359448E-01_EB, 1.159382E-01_EB, 9.908354E-02_EB, &
    8.473452E-02_EB, 7.233655E-02_EB, 6.187191E-02_EB, 5.286391E-02_EB, 4.518233E-02_EB, 3.863453E-02_EB, &
    3.288079E-02_EB, 2.803528E-02_EB, 2.394608E-02_EB, 2.051761E-02_EB, 1.743587E-02_EB,  &
    1.334409E+00_EB, 1.095178E+00_EB, 9.120448E-01_EB, 7.631562E-01_EB, 6.425576E-01_EB, 5.444728E-01_EB, &
    4.621443E-01_EB, 3.938260E-01_EB, 3.368579E-01_EB, 2.880892E-01_EB, 2.459378E-01_EB, 2.104005E-01_EB, &
    1.799327E-01_EB, 1.535802E-01_EB, 1.314063E-01_EB, 1.121907E-01_EB, 9.581696E-02_EB, 8.185813E-02_EB, &
    6.964125E-02_EB, 5.934825E-02_EB, 5.067442E-02_EB, 4.339843E-02_EB, 3.686411E-02_EB/),(/23,8/))

SD2_CH4(1:23,17:23) = RESHAPE((/ &  ! 3100-3250 cm-1
    9.674073E-01_EB, 9.013091E-01_EB, 8.256467E-01_EB, 7.423397E-01_EB, 6.619262E-01_EB, 5.879801E-01_EB, &
    5.191373E-01_EB, 4.573467E-01_EB, 4.022801E-01_EB, 3.523929E-01_EB, 3.072101E-01_EB, 2.674592E-01_EB, &
    2.326179E-01_EB, 2.012004E-01_EB, 1.742107E-01_EB, 1.505433E-01_EB, 1.298551E-01_EB, 1.119451E-01_EB, &
    9.604141E-02_EB, 8.242372E-02_EB, 7.089244E-02_EB, 6.110245E-02_EB, 5.221546E-02_EB,  &
    3.658653E-01_EB, 4.184364E-01_EB, 4.441884E-01_EB, 4.457809E-01_EB, 4.319013E-01_EB, 4.088715E-01_EB, &
    3.793819E-01_EB, 3.474065E-01_EB, 3.151904E-01_EB, 2.830904E-01_EB, 2.516607E-01_EB, 2.228145E-01_EB, &
    1.964195E-01_EB, 1.718178E-01_EB, 1.502324E-01_EB, 1.307846E-01_EB, 1.135071E-01_EB, 9.842138E-02_EB, &
    8.481592E-02_EB, 7.309981E-02_EB, 6.307130E-02_EB, 5.451710E-02_EB, 4.669207E-02_EB,  &
    8.373191E-02_EB, 1.229238E-01_EB, 1.573970E-01_EB, 1.825905E-01_EB, 1.975753E-01_EB, 2.040052E-01_EB, &
    2.024706E-01_EB, 1.957505E-01_EB, 1.855755E-01_EB, 1.726368E-01_EB, 1.580202E-01_EB, 1.433211E-01_EB, &
    1.287520E-01_EB, 1.144879E-01_EB, 1.013973E-01_EB, 8.925648E-02_EB, 7.818868E-02_EB, 6.826146E-02_EB, &
    5.920229E-02_EB, 5.129713E-02_EB, 4.441206E-02_EB, 3.851923E-02_EB, 3.306979E-02_EB,  &
    1.457007E-02_EB, 2.736439E-02_EB, 4.340085E-02_EB, 6.002981E-02_EB, 7.506875E-02_EB, 8.727601E-02_EB, &
    9.581284E-02_EB, 1.007875E-01_EB, 1.027285E-01_EB, 1.017913E-01_EB, 9.839048E-02_EB, 9.365514E-02_EB, &
    8.789001E-02_EB, 8.123320E-02_EB, 7.449732E-02_EB, 6.766345E-02_EB, 6.100922E-02_EB, 5.470686E-02_EB, &
    4.865358E-02_EB, 4.309088E-02_EB, 3.810993E-02_EB, 3.369467E-02_EB, 2.945916E-02_EB,  &
    1.191095E-03_EB, 1.552045E-03_EB, 2.059789E-03_EB, 2.649547E-03_EB, 3.236595E-03_EB, 3.762561E-03_EB, &
    4.165025E-03_EB, 4.448429E-03_EB, 4.600555E-03_EB, 4.628224E-03_EB, 4.534948E-03_EB, 4.374710E-03_EB, &
    4.151710E-03_EB, 3.880527E-03_EB, 3.599880E-03_EB, 3.300161E-03_EB, 2.998249E-03_EB, 2.715486E-03_EB, &
    2.425276E-03_EB, 2.164309E-03_EB, 1.929373E-03_EB, 1.713087E-03_EB, 1.506010E-03_EB,  &
    4.538291E-04_EB, 4.599393E-04_EB, 4.587246E-04_EB, 4.526653E-04_EB, 4.375098E-04_EB, 4.167315E-04_EB, &
    3.918634E-04_EB, 3.644121E-04_EB, 3.335314E-04_EB, 3.072189E-04_EB, 2.753771E-04_EB, 2.453914E-04_EB, &
    2.183408E-04_EB, 1.968316E-04_EB, 1.704817E-04_EB, 1.507858E-04_EB, 1.319531E-04_EB, 1.134243E-04_EB, &
    9.801919E-05_EB, 8.900765E-05_EB, 7.481087E-05_EB, 6.511635E-05_EB, 5.390619E-05_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,7/))

!-------------------------Propane DATA-------------------


! THERE ARE 2 BANDS FOR Propane

! BAND #1: 1300 cm-1 - 1600 cm-1 
! BAND #2: 2700 cm-1 - 3200 cm-1 

N_TEMP_C3H8 = 23
N_BAND_C3H8 = 2

ALLOCATE(SD_C3H8_TEMP(N_TEMP_C3H8)) 

! INITIALIZE BANDS WAVENUMBER BOUNDS FOR Propane ABSOPTION COEFFICIENTS.
! BANDS ARE ORGANIZED BY ROW: 
! 1ST COLUMN IS THE LOWER BOUND BAND LIMIT IN cm-1. 
! 2ND COLUMN IS THE UPPER BOUND BAND LIMIT IN cm-1. 
! 3RD COLUMN IS THE STRIDE BETWEEN WAVENUMBERS IN BAND.
! IF 3RD COLUMN = 0., THEN THE BAND IS CALCULATED AND NOT TABULATED.

ALLOCATE(OM_BND_C3H8(N_BAND_C3H8,3)) 

OM_BND_C3H8 = RESHAPE((/ &
    1300._EB, 2700._EB, &
    1600._EB, 3200._EB, &
     25._EB, 25._EB/),(/N_BAND_C3H8,3/)) 

SD_C3H8_TEMP = (/ &
    300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
    600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
    900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
    1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(SD1_C3H8(N_TEMP_C3H8,13)) 

! BAND #1: 1300 cm-1 - 1600 cm-1 

SD1_C3H8(1:23,1:8) = RESHAPE((/ &  ! 1300-1475 cm-1
    8.203975E-02_EB, 1.025595E-01_EB, 1.192993E-01_EB, 1.320128E-01_EB, 1.409410E-01_EB, 1.465632E-01_EB, &
    1.494326E-01_EB, 1.500917E-01_EB, 1.490304E-01_EB, 1.466629E-01_EB, 1.433344E-01_EB, 1.393226E-01_EB, &
    1.348485E-01_EB, 1.300827E-01_EB, 1.251594E-01_EB, 1.201775E-01_EB, 1.152163E-01_EB, 1.103297E-01_EB, &
    1.055578E-01_EB, 1.009285E-01_EB, 9.646159E-02_EB, 9.216865E-02_EB, 8.805680E-02_EB,  &
    1.737969E-01_EB, 1.740925E-01_EB, 1.715941E-01_EB, 1.674367E-01_EB, 1.623398E-01_EB, 1.567598E-01_EB, &
    1.509871E-01_EB, 1.452056E-01_EB, 1.395301E-01_EB, 1.340326E-01_EB, 1.287505E-01_EB, 1.237091E-01_EB, &
    1.189163E-01_EB, 1.143701E-01_EB, 1.100666E-01_EB, 1.059965E-01_EB, 1.021512E-01_EB, 9.851748E-02_EB, &
    9.508307E-02_EB, 9.183684E-02_EB, 8.876651E-02_EB, 8.586167E-02_EB, 8.311109E-02_EB,  &
    3.827428E-01_EB, 3.550833E-01_EB, 3.295953E-01_EB, 3.063920E-01_EB, 2.853372E-01_EB, 2.662267E-01_EB, &
    2.488572E-01_EB, 2.330396E-01_EB, 2.186084E-01_EB, 2.054164E-01_EB, 1.933331E-01_EB, 1.822448E-01_EB, &
    1.720521E-01_EB, 1.626641E-01_EB, 1.540036E-01_EB, 1.459995E-01_EB, 1.385918E-01_EB, 1.317223E-01_EB, &
    1.253426E-01_EB, 1.194105E-01_EB, 1.138848E-01_EB, 1.087298E-01_EB, 1.039142E-01_EB,  &
    5.812751E-01_EB, 4.885960E-01_EB, 4.174531E-01_EB, 3.613475E-01_EB, 3.161210E-01_EB, 2.790101E-01_EB, &
    2.481120E-01_EB, 2.220658E-01_EB, 1.998851E-01_EB, 1.808271E-01_EB, 1.643250E-01_EB, 1.499385E-01_EB, &
    1.373216E-01_EB, 1.261939E-01_EB, 1.163316E-01_EB, 1.075527E-01_EB, 9.970488E-02_EB, 9.266249E-02_EB, &
    8.631973E-02_EB, 8.058960E-02_EB, 7.539561E-02_EB, 7.067384E-02_EB, 6.637024E-02_EB,  &
    2.897005E-01_EB, 2.987455E-01_EB, 3.046227E-01_EB, 3.081640E-01_EB, 3.099165E-01_EB, 3.102657E-01_EB, &
    3.095051E-01_EB, 3.078639E-01_EB, 3.055268E-01_EB, 3.026418E-01_EB, 2.993323E-01_EB, 2.956948E-01_EB, &
    2.918139E-01_EB, 2.877536E-01_EB, 2.835679E-01_EB, 2.793007E-01_EB, 2.749899E-01_EB, 2.706625E-01_EB, &
    2.663418E-01_EB, 2.620468E-01_EB, 2.577932E-01_EB, 2.535914E-01_EB, 2.494510E-01_EB,  &
    5.628545E-01_EB, 5.054634E-01_EB, 4.607537E-01_EB, 4.245596E-01_EB, 3.943723E-01_EB, 3.685998E-01_EB, &
    3.461927E-01_EB, 3.264330E-01_EB, 3.088095E-01_EB, 2.929498E-01_EB, 2.785722E-01_EB, 2.654593E-01_EB, &
    2.534400E-01_EB, 2.423749E-01_EB, 2.321503E-01_EB, 2.226736E-01_EB, 2.138627E-01_EB, 2.056504E-01_EB, &
    1.979785E-01_EB, 1.907958E-01_EB, 1.840576E-01_EB, 1.777269E-01_EB, 1.717658E-01_EB,  &
    1.112299E+00_EB, 8.967789E-01_EB, 7.436999E-01_EB, 6.298272E-01_EB, 5.420622E-01_EB, 4.725162E-01_EB, &
    4.161761E-01_EB, 3.697100E-01_EB, 3.308193E-01_EB, 2.978686E-01_EB, 2.696563E-01_EB, 2.452875E-01_EB, &
    2.240753E-01_EB, 2.054812E-01_EB, 1.890865E-01_EB, 1.745515E-01_EB, 1.616040E-01_EB, 1.500183E-01_EB, &
    1.396097E-01_EB, 1.302233E-01_EB, 1.217297E-01_EB, 1.140202E-01_EB, 1.070013E-01_EB,  &
    7.888871E-01_EB, 6.401361E-01_EB, 5.335745E-01_EB, 4.537499E-01_EB, 3.918660E-01_EB, 3.425807E-01_EB, &
    3.024797E-01_EB, 2.692780E-01_EB, 2.413942E-01_EB, 2.176968E-01_EB, 1.973512E-01_EB, 1.797354E-01_EB, &
    1.643639E-01_EB, 1.508658E-01_EB, 1.389407E-01_EB, 1.283507E-01_EB, 1.189044E-01_EB, 1.104390E-01_EB, &
    1.028240E-01_EB, 9.594866E-02_EB, 8.972154E-02_EB, 8.406339E-02_EB, 7.890709E-02_EB/),(/23,8/))

SD1_C3H8(1:23,9:13) = RESHAPE((/ &  ! 1500-1600 cm-1
    3.219275E-01_EB, 2.933374E-01_EB, 2.710503E-01_EB, 2.529531E-01_EB, 2.377785E-01_EB, 2.247315E-01_EB, &
    2.132922E-01_EB, 2.031085E-01_EB, 1.939346E-01_EB, 1.855919E-01_EB, 1.779489E-01_EB, 1.709040E-01_EB, &
    1.643816E-01_EB, 1.583143E-01_EB, 1.526539E-01_EB, 1.473576E-01_EB, 1.423886E-01_EB, 1.377171E-01_EB, &
    1.333163E-01_EB, 1.291619E-01_EB, 1.252349E-01_EB, 1.215191E-01_EB, 1.179941E-01_EB,  &
    6.321990E-02_EB, 7.120389E-02_EB, 8.040725E-02_EB, 9.063308E-02_EB, 1.015952E-01_EB, 1.129846E-01_EB, &
    1.245185E-01_EB, 1.359572E-01_EB, 1.471110E-01_EB, 1.578347E-01_EB, 1.680289E-01_EB, 1.776274E-01_EB, &
    1.865881E-01_EB, 1.948961E-01_EB, 2.025510E-01_EB, 2.095642E-01_EB, 2.159540E-01_EB, 2.217488E-01_EB, &
    2.269759E-01_EB, 2.316698E-01_EB, 2.358607E-01_EB, 2.395826E-01_EB, 2.428673E-01_EB,  &
    1.880192E-02_EB, 2.855895E-02_EB, 3.913282E-02_EB, 4.997389E-02_EB, 6.067279E-02_EB, 7.094818E-02_EB, &
    8.061431E-02_EB, 8.956418E-02_EB, 9.774888E-02_EB, 1.051567E-01_EB, 1.118019E-01_EB, 1.177188E-01_EB, &
    1.229500E-01_EB, 1.275434E-01_EB, 1.315492E-01_EB, 1.350196E-01_EB, 1.380002E-01_EB, 1.405395E-01_EB, &
    1.426802E-01_EB, 1.444622E-01_EB, 1.459215E-01_EB, 1.470909E-01_EB, 1.480031E-01_EB,  &
    2.609848E-02_EB, 3.410624E-02_EB, 4.129149E-02_EB, 4.753321E-02_EB, 5.283784E-02_EB, 5.726710E-02_EB, &
    6.091036E-02_EB, 6.385853E-02_EB, 6.620344E-02_EB, 6.802440E-02_EB, 6.939964E-02_EB, 7.039318E-02_EB, &
    7.106279E-02_EB, 7.145846E-02_EB, 7.162160E-02_EB, 7.159069E-02_EB, 7.139654E-02_EB, 7.106573E-02_EB, &
    7.062225E-02_EB, 7.008588E-02_EB, 6.947162E-02_EB, 6.879546E-02_EB, 6.806812E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

ALLOCATE(SD2_C3H8(N_TEMP_C3H8,21)) 

! BAND #2: 2700 cm-1 - 3200 cm-1 

SD2_C3H8(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
    4.263288E-02_EB, 4.640354E-02_EB, 4.788356E-02_EB, 4.785035E-02_EB, 4.687052E-02_EB, 4.532933E-02_EB, &
    4.347875E-02_EB, 4.148309E-02_EB, 3.944429E-02_EB, 3.742594E-02_EB, 3.546708E-02_EB, 3.359049E-02_EB, &
    3.180706E-02_EB, 3.012097E-02_EB, 2.853484E-02_EB, 2.704526E-02_EB, 2.564667E-02_EB, 2.433839E-02_EB, &
    2.311301E-02_EB, 2.196575E-02_EB, 2.089234E-02_EB, 1.988555E-02_EB, 1.894305E-02_EB,  &
    1.282935E-01_EB, 1.200918E-01_EB, 1.104850E-01_EB, 1.008318E-01_EB, 9.175400E-02_EB, 8.348541E-02_EB, &
    7.607517E-02_EB, 6.948444E-02_EB, 6.363894E-02_EB, 5.845564E-02_EB, 5.385369E-02_EB, 4.975870E-02_EB, &
    4.610552E-02_EB, 4.283445E-02_EB, 3.989880E-02_EB, 3.725538E-02_EB, 3.486722E-02_EB, 3.270395E-02_EB, &
    3.073845E-02_EB, 2.894746E-02_EB, 2.731293E-02_EB, 2.581496E-02_EB, 2.444106E-02_EB,  &
    1.616155E-01_EB, 1.510497E-01_EB, 1.415692E-01_EB, 1.331431E-01_EB, 1.256577E-01_EB, 1.189832E-01_EB, &
    1.130001E-01_EB, 1.076030E-01_EB, 1.027065E-01_EB, 9.823662E-02_EB, 9.413358E-02_EB, 9.034827E-02_EB, &
    8.684063E-02_EB, 8.357597E-02_EB, 8.052840E-02_EB, 7.767451E-02_EB, 7.499103E-02_EB, 7.246490E-02_EB, &
    7.007974E-02_EB, 6.782231E-02_EB, 6.568657E-02_EB, 6.365791E-02_EB, 6.172891E-02_EB,  &
    1.734183E-01_EB, 1.645080E-01_EB, 1.556924E-01_EB, 1.473756E-01_EB, 1.396814E-01_EB, 1.326192E-01_EB, &
    1.261514E-01_EB, 1.202206E-01_EB, 1.147685E-01_EB, 1.097411E-01_EB, 1.050907E-01_EB, 1.007719E-01_EB, &
    9.674969E-02_EB, 9.299234E-02_EB, 8.947384E-02_EB, 8.616778E-02_EB, 8.305758E-02_EB, 8.012406E-02_EB, &
    7.735210E-02_EB, 7.472914E-02_EB, 7.224318E-02_EB, 6.988420E-02_EB, 6.764195E-02_EB,  &
    2.293196E-01_EB, 2.292635E-01_EB, 2.268458E-01_EB, 2.231654E-01_EB, 2.188284E-01_EB, 2.141781E-01_EB, &
    2.094091E-01_EB, 2.046326E-01_EB, 1.999092E-01_EB, 1.952760E-01_EB, 1.907517E-01_EB, 1.863415E-01_EB, &
    1.820517E-01_EB, 1.778814E-01_EB, 1.738302E-01_EB, 1.698952E-01_EB, 1.660725E-01_EB, 1.623613E-01_EB, &
    1.587575E-01_EB, 1.552588E-01_EB, 1.518615E-01_EB, 1.485625E-01_EB, 1.453616E-01_EB,  &
    5.215213E-01_EB, 5.049996E-01_EB, 4.888627E-01_EB, 4.735790E-01_EB, 4.592761E-01_EB, 4.459396E-01_EB, &
    4.334924E-01_EB, 4.218444E-01_EB, 4.109017E-01_EB, 4.005781E-01_EB, 3.907984E-01_EB, 3.815006E-01_EB, &
    3.726289E-01_EB, 3.641388E-01_EB, 3.559919E-01_EB, 3.481598E-01_EB, 3.406142E-01_EB, 3.333338E-01_EB, &
    3.263011E-01_EB, 3.195006E-01_EB, 3.129175E-01_EB, 3.065412E-01_EB, 3.003622E-01_EB,  &
    2.334470E+00_EB, 2.035823E+00_EB, 1.805796E+00_EB, 1.624173E+00_EB, 1.477603E+00_EB, 1.357040E+00_EB, &
    1.256187E+00_EB, 1.170563E+00_EB, 1.096906E+00_EB, 1.032800E+00_EB, 9.764217E-01_EB, 9.263800E-01_EB, &
    8.815949E-01_EB, 8.412163E-01_EB, 8.045723E-01_EB, 7.711233E-01_EB, 7.404275E-01_EB, 7.121285E-01_EB, &
    6.859289E-01_EB, 6.615801E-01_EB, 6.388762E-01_EB, 6.176400E-01_EB, 5.977224E-01_EB,  &
    4.234659E+00_EB, 3.466039E+00_EB, 2.915309E+00_EB, 2.503486E+00_EB, 2.185140E+00_EB, 1.932405E+00_EB, &
    1.727314E+00_EB, 1.557799E+00_EB, 1.415486E+00_EB, 1.294394E+00_EB, 1.190156E+00_EB, 1.099511E+00_EB, &
    1.019985E+00_EB, 9.496646E-01_EB, 8.870523E-01_EB, 8.309597E-01_EB, 7.804321E-01_EB, 7.346940E-01_EB, &
    6.931057E-01_EB, 6.551434E-01_EB, 6.203635E-01_EB, 5.883968E-01_EB, 5.589247E-01_EB/),(/23,8/))

SD2_C3H8(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
    4.907292E+00_EB, 4.273262E+00_EB, 3.791498E+00_EB, 3.412411E+00_EB, 3.105847E+00_EB, 2.852384E+00_EB, &
    2.638943E+00_EB, 2.456391E+00_EB, 2.298167E+00_EB, 2.159444E+00_EB, 2.036593E+00_EB, 1.926843E+00_EB, &
    1.828042E+00_EB, 1.738498E+00_EB, 1.656859E+00_EB, 1.582038E+00_EB, 1.513147E+00_EB, 1.449454E+00_EB, &
    1.390352E+00_EB, 1.335331E+00_EB, 1.283957E+00_EB, 1.235863E+00_EB, 1.190733E+00_EB,  &
    6.937833E+00_EB, 6.025116E+00_EB, 5.299461E+00_EB, 4.712045E+00_EB, 4.228735E+00_EB, 3.825250E+00_EB, &
    3.484000E+00_EB, 3.192020E+00_EB, 2.939592E+00_EB, 2.719329E+00_EB, 2.525531E+00_EB, 2.353754E+00_EB, &
    2.200484E+00_EB, 2.062914E+00_EB, 1.938775E+00_EB, 1.826219E+00_EB, 1.723728E+00_EB, 1.630036E+00_EB, &
    1.544086E+00_EB, 1.464987E+00_EB, 1.391977E+00_EB, 1.324412E+00_EB, 1.261729E+00_EB,  &
    1.308369E+01_EB, 1.033997E+01_EB, 8.405621E+00_EB, 6.987186E+00_EB, 5.913710E+00_EB, 5.079934E+00_EB, &
    4.418088E+00_EB, 3.882899E+00_EB, 3.443181E+00_EB, 3.076866E+00_EB, 2.767987E+00_EB, 2.504747E+00_EB, &
    2.278280E+00_EB, 2.081811E+00_EB, 1.910088E+00_EB, 1.758986E+00_EB, 1.625227E+00_EB, 1.506174E+00_EB, &
    1.399688E+00_EB, 1.304014E+00_EB, 1.217701E+00_EB, 1.139544E+00_EB, 1.068524E+00_EB,  &
    9.306308E+00_EB, 7.481934E+00_EB, 6.149918E+00_EB, 5.146688E+00_EB, 4.371798E+00_EB, 3.760528E+00_EB, &
    3.269574E+00_EB, 2.869071E+00_EB, 2.537880E+00_EB, 2.260708E+00_EB, 2.026270E+00_EB, 1.826103E+00_EB, &
    1.653750E+00_EB, 1.504224E+00_EB, 1.373617E+00_EB, 1.258835E+00_EB, 1.157400E+00_EB, 1.067305E+00_EB, &
    9.869156E-01_EB, 9.148808E-01_EB, 8.500868E-01_EB, 7.915940E-01_EB, 7.386168E-01_EB,  &
    2.115655E+00_EB, 2.066169E+00_EB, 1.962833E+00_EB, 1.836226E+00_EB, 1.703241E+00_EB, 1.572807E+00_EB, &
    1.449343E+00_EB, 1.334758E+00_EB, 1.229600E+00_EB, 1.133693E+00_EB, 1.046512E+00_EB, 9.673725E-01_EB, &
    8.955495E-01_EB, 8.303322E-01_EB, 7.710519E-01_EB, 7.171005E-01_EB, 6.679244E-01_EB, 6.230305E-01_EB, &
    5.819797E-01_EB, 5.443812E-01_EB, 5.098886E-01_EB, 4.781946E-01_EB, 4.490253E-01_EB,  &
    2.097634E-01_EB, 2.788466E-01_EB, 3.339837E-01_EB, 3.744436E-01_EB, 4.018120E-01_EB, 4.184140E-01_EB, &
    4.265634E-01_EB, 4.282740E-01_EB, 4.251898E-01_EB, 4.186081E-01_EB, 4.095357E-01_EB, 3.987385E-01_EB, &
    3.867979E-01_EB, 3.741541E-01_EB, 3.611337E-01_EB, 3.479841E-01_EB, 3.348822E-01_EB, 3.219640E-01_EB, &
    3.093208E-01_EB, 2.970228E-01_EB, 2.851121E-01_EB, 2.736194E-01_EB, 2.625615E-01_EB,  &
    6.093969E-02_EB, 8.506670E-02_EB, 1.057939E-01_EB, 1.222662E-01_EB, 1.345747E-01_EB, 1.432152E-01_EB, &
    1.488044E-01_EB, 1.519406E-01_EB, 1.531493E-01_EB, 1.528669E-01_EB, 1.514489E-01_EB, 1.491776E-01_EB, &
    1.462790E-01_EB, 1.429270E-01_EB, 1.392589E-01_EB, 1.353832E-01_EB, 1.313814E-01_EB, 1.273175E-01_EB, &
    1.232436E-01_EB, 1.191945E-01_EB, 1.152004E-01_EB, 1.112824E-01_EB, 1.074570E-01_EB,  &
    5.631146E-02_EB, 7.529887E-02_EB, 9.023923E-02_EB, 1.009473E-01_EB, 1.079248E-01_EB, 1.118747E-01_EB, &
    1.134892E-01_EB, 1.133565E-01_EB, 1.119533E-01_EB, 1.096477E-01_EB, 1.067217E-01_EB, 1.033861E-01_EB, &
    9.979894E-02_EB, 9.607594E-02_EB, 9.230299E-02_EB, 8.854075E-02_EB, 8.483530E-02_EB, 8.121631E-02_EB, &
    7.770611E-02_EB, 7.431770E-02_EB, 7.105960E-02_EB, 6.793766E-02_EB, 6.495208E-02_EB/),(/23,8/))

SD2_C3H8(1:23,17:21) = RESHAPE((/ &  ! 3100-3200 cm-1
    6.902145E-02_EB, 8.265347E-02_EB, 9.170932E-02_EB, 9.694038E-02_EB, 9.924906E-02_EB, 9.944236E-02_EB, &
    9.816219E-02_EB, 9.589395E-02_EB, 9.298693E-02_EB, 8.969575E-02_EB, 8.619570E-02_EB, 8.261261E-02_EB, &
    7.902878E-02_EB, 7.550604E-02_EB, 7.207904E-02_EB, 6.877272E-02_EB, 6.560020E-02_EB, 6.257083E-02_EB, &
    5.968714E-02_EB, 5.694706E-02_EB, 5.434865E-02_EB, 5.188782E-02_EB, 4.955802E-02_EB,  &
    9.643939E-02_EB, 1.016765E-01_EB, 1.019854E-01_EB, 9.936005E-02_EB, 9.512153E-02_EB, 9.009745E-02_EB, &
    8.478784E-02_EB, 7.948604E-02_EB, 7.435671E-02_EB, 6.948794E-02_EB, 6.491906E-02_EB, 6.066257E-02_EB, &
    5.671572E-02_EB, 5.306350E-02_EB, 4.969151E-02_EB, 4.657763E-02_EB, 4.370371E-02_EB, 4.104902E-02_EB, &
    3.859660E-02_EB, 3.632998E-02_EB, 3.423007E-02_EB, 3.228547E-02_EB, 3.048175E-02_EB,  &
    1.275255E-01_EB, 1.209418E-01_EB, 1.114686E-01_EB, 1.012405E-01_EB, 9.129826E-02_EB, 8.208997E-02_EB, &
    7.377494E-02_EB, 6.636074E-02_EB, 5.979550E-02_EB, 5.399532E-02_EB, 4.887543E-02_EB, 4.435102E-02_EB, &
    4.034604E-02_EB, 3.679147E-02_EB, 3.363081E-02_EB, 3.081204E-02_EB, 2.829140E-02_EB, 2.603255E-02_EB, &
    2.400099E-02_EB, 2.217151E-02_EB, 2.051857E-02_EB, 1.902142E-02_EB, 1.766427E-02_EB,  &
    1.409411E-01_EB, 1.258528E-01_EB, 1.108148E-01_EB, 9.712610E-02_EB, 8.513238E-02_EB, 7.480130E-02_EB, &
    6.595768E-02_EB, 5.839678E-02_EB, 5.192055E-02_EB, 4.635519E-02_EB, 4.155373E-02_EB, 3.739176E-02_EB, &
    3.376836E-02_EB, 3.059929E-02_EB, 2.781557E-02_EB, 2.535991E-02_EB, 2.318645E-02_EB, 2.125369E-02_EB, &
    1.953074E-02_EB, 1.798800E-02_EB, 1.660297E-02_EB, 1.535635E-02_EB, 1.423042E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

!-------------------------Heptane DATA-------------------


! THERE ARE 2 BANDS FOR Heptane

! BAND #1: 1300 cm-1 - 1550 cm-1 
! BAND #2: 2700 cm-1 - 3200 cm-1 

N_TEMP_C7H16 = 23
N_BAND_C7H16 = 2

ALLOCATE(SD_C7H16_TEMP(N_TEMP_C7H16)) 

! INITIALIZE BANDS WAVENUMBER BOUNDS FOR Heptane ABSOPTION COEFFICIENTS.
! BANDS ARE ORGANIZED BY ROW: 
! 1ST COLUMN IS THE LOWER BOUND BAND LIMIT IN cm-1. 
! 2ND COLUMN IS THE UPPER BOUND BAND LIMIT IN cm-1. 
! 3RD COLUMN IS THE STRIDE BETWEEN WAVENUMBERS IN BAND.
! IF 3RD COLUMN = 0., THEN THE BAND IS CALCULATED AND NOT TABULATED.

ALLOCATE(OM_BND_C7H16(N_BAND_C7H16,3)) 

OM_BND_C7H16 = RESHAPE((/ &
    1300._EB, 2700._EB, &
    1550._EB, 3200._EB, &
     25._EB, 25._EB/),(/N_BAND_C7H16,3/)) 

SD_C7H16_TEMP = (/ &
    300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
    600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
    900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
    1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(SD1_C7H16(N_TEMP_C7H16,11)) 

! BAND #1: 1300 cm-1 - 1550 cm-1 

SD1_C7H16(1:23,1:8) = RESHAPE((/ &  ! 1300-1475 cm-1
    1.059314E-01_EB, 1.211256E-01_EB, 1.335611E-01_EB, 1.436086E-01_EB, 1.516184E-01_EB, 1.579036E-01_EB, &
    1.627367E-01_EB, 1.663525E-01_EB, 1.689510E-01_EB, 1.707013E-01_EB, 1.717464E-01_EB, 1.722066E-01_EB, &
    1.721829E-01_EB, 1.717598E-01_EB, 1.710087E-01_EB, 1.699887E-01_EB, 1.687496E-01_EB, 1.673334E-01_EB, &
    1.657750E-01_EB, 1.641038E-01_EB, 1.623443E-01_EB, 1.605174E-01_EB, 1.586401E-01_EB,  &
    2.201113E-01_EB, 2.384142E-01_EB, 2.508856E-01_EB, 2.589997E-01_EB, 2.638285E-01_EB, 2.661649E-01_EB, &
    2.666050E-01_EB, 2.656031E-01_EB, 2.635092E-01_EB, 2.605941E-01_EB, 2.570686E-01_EB, 2.530979E-01_EB, &
    2.488103E-01_EB, 2.443071E-01_EB, 2.396668E-01_EB, 2.349507E-01_EB, 2.302071E-01_EB, 2.254731E-01_EB, &
    2.207771E-01_EB, 2.161411E-01_EB, 2.115815E-01_EB, 2.071106E-01_EB, 2.027370E-01_EB,  &
    7.465515E-01_EB, 6.967534E-01_EB, 6.561642E-01_EB, 6.216858E-01_EB, 5.915058E-01_EB, 5.645151E-01_EB, &
    5.400067E-01_EB, 5.175113E-01_EB, 4.967046E-01_EB, 4.773522E-01_EB, 4.592786E-01_EB, 4.423460E-01_EB, &
    4.264433E-01_EB, 4.114774E-01_EB, 3.973699E-01_EB, 3.840516E-01_EB, 3.714620E-01_EB, 3.595468E-01_EB, &
    3.482575E-01_EB, 3.375501E-01_EB, 3.273843E-01_EB, 3.177237E-01_EB, 3.085346E-01_EB,  &
    1.763409E+00_EB, 1.472903E+00_EB, 1.259466E+00_EB, 1.095869E+00_EB, 9.663397E-01_EB, 8.611613E-01_EB, &
    7.740307E-01_EB, 7.006806E-01_EB, 6.381136E-01_EB, 5.841583E-01_EB, 5.371982E-01_EB, 4.960023E-01_EB, &
    4.596146E-01_EB, 4.272802E-01_EB, 3.983938E-01_EB, 3.724640E-01_EB, 3.490874E-01_EB, 3.279299E-01_EB, &
    3.087116E-01_EB, 2.911971E-01_EB, 2.751861E-01_EB, 2.605081E-01_EB, 2.470161E-01_EB,  &
    1.910890E-01_EB, 2.183836E-01_EB, 2.430932E-01_EB, 2.651323E-01_EB, 2.845390E-01_EB, 3.014334E-01_EB, &
    3.159854E-01_EB, 3.283905E-01_EB, 3.388526E-01_EB, 3.475730E-01_EB, 3.547434E-01_EB, 3.605424E-01_EB, &
    3.651328E-01_EB, 3.686619E-01_EB, 3.712618E-01_EB, 3.730501E-01_EB, 3.741306E-01_EB, 3.745954E-01_EB, &
    3.745253E-01_EB, 3.739917E-01_EB, 3.730569E-01_EB, 3.717759E-01_EB, 3.701966E-01_EB,  &
    5.910628E-01_EB, 5.998443E-01_EB, 6.075257E-01_EB, 6.135952E-01_EB, 6.178881E-01_EB, 6.204261E-01_EB, &
    6.213296E-01_EB, 6.207646E-01_EB, 6.189136E-01_EB, 6.159572E-01_EB, 6.120656E-01_EB, 6.073931E-01_EB, &
    6.020773E-01_EB, 5.962388E-01_EB, 5.899821E-01_EB, 5.833970E-01_EB, 5.765601E-01_EB, 5.695369E-01_EB, &
    5.623825E-01_EB, 5.551432E-01_EB, 5.478584E-01_EB, 5.405604E-01_EB, 5.332765E-01_EB,  &
    3.396756E+00_EB, 2.764209E+00_EB, 2.294017E+00_EB, 1.934584E+00_EB, 1.653246E+00_EB, 1.428619E+00_EB, &
    1.246251E+00_EB, 1.096081E+00_EB, 9.709160E-01_EB, 8.654919E-01_EB, 7.758770E-01_EB, 6.990803E-01_EB, &
    6.327913E-01_EB, 5.751989E-01_EB, 5.248669E-01_EB, 4.806433E-01_EB, 4.415948E-01_EB, 4.069594E-01_EB, &
    3.761087E-01_EB, 3.485220E-01_EB, 3.237640E-01_EB, 3.014687E-01_EB, 2.813272E-01_EB,  &
    1.314273E+00_EB, 1.111468E+00_EB, 9.573636E-01_EB, 8.372661E-01_EB, 7.415287E-01_EB, 6.636950E-01_EB, &
    5.993381E-01_EB, 5.453461E-01_EB, 4.994788E-01_EB, 4.600878E-01_EB, 4.259355E-01_EB, 3.960766E-01_EB, &
    3.697765E-01_EB, 3.464569E-01_EB, 3.256563E-01_EB, 3.070020E-01_EB, 2.901902E-01_EB, 2.749709E-01_EB, &
    2.611365E-01_EB, 2.485129E-01_EB, 2.369540E-01_EB, 2.263355E-01_EB, 2.165510E-01_EB/),(/23,8/))

SD1_C7H16(1:23,9:11) = RESHAPE((/ &  ! 1500-1550 cm-1
    1.015621E-01_EB, 1.170190E-01_EB, 1.311723E-01_EB, 1.438600E-01_EB, 1.550615E-01_EB, 1.648312E-01_EB, &
    1.732628E-01_EB, 1.804675E-01_EB, 1.865629E-01_EB, 1.916644E-01_EB, 1.958814E-01_EB, 1.993160E-01_EB, &
    2.020609E-01_EB, 2.042000E-01_EB, 2.058082E-01_EB, 2.069521E-01_EB, 2.076909E-01_EB, 2.080763E-01_EB, &
    2.081547E-01_EB, 2.079664E-01_EB, 2.075467E-01_EB, 2.069267E-01_EB, 2.061341E-01_EB,  &
    6.046052E-02_EB, 6.366702E-02_EB, 6.721068E-02_EB, 7.074050E-02_EB, 7.410892E-02_EB, 7.725069E-02_EB, &
    8.013762E-02_EB, 8.276042E-02_EB, 8.512013E-02_EB, 8.722434E-02_EB, 8.908435E-02_EB, 9.071382E-02_EB, &
    9.212777E-02_EB, 9.334151E-02_EB, 9.437059E-02_EB, 9.522996E-02_EB, 9.593414E-02_EB, 9.649686E-02_EB, &
    9.693098E-02_EB, 9.724847E-02_EB, 9.746045E-02_EB, 9.757711E-02_EB, 9.760795E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,3/))

ALLOCATE(SD2_C7H16(N_TEMP_C7H16,21)) 

! BAND #2: 2700 cm-1 - 3200 cm-1 

SD2_C7H16(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
    1.165783E-01_EB, 1.229222E-01_EB, 1.270066E-01_EB, 1.295172E-01_EB, 1.309151E-01_EB, 1.315130E-01_EB, &
    1.315250E-01_EB, 1.311009E-01_EB, 1.303473E-01_EB, 1.293411E-01_EB, 1.281401E-01_EB, 1.267879E-01_EB, &
    1.253180E-01_EB, 1.237573E-01_EB, 1.221270E-01_EB, 1.204446E-01_EB, 1.187242E-01_EB, 1.169778E-01_EB, &
    1.152151E-01_EB, 1.134445E-01_EB, 1.116728E-01_EB, 1.099058E-01_EB, 1.081485E-01_EB,  &
    2.771574E-01_EB, 2.570380E-01_EB, 2.402856E-01_EB, 2.260838E-01_EB, 2.138541E-01_EB, 2.031760E-01_EB, &
    1.937374E-01_EB, 1.853027E-01_EB, 1.776916E-01_EB, 1.707648E-01_EB, 1.644128E-01_EB, 1.585497E-01_EB, &
    1.531067E-01_EB, 1.480286E-01_EB, 1.432706E-01_EB, 1.387961E-01_EB, 1.345746E-01_EB, 1.305807E-01_EB, &
    1.267930E-01_EB, 1.231934E-01_EB, 1.197662E-01_EB, 1.164980E-01_EB, 1.133766E-01_EB,  &
    1.888950E-01_EB, 1.952526E-01_EB, 1.999970E-01_EB, 2.036301E-01_EB, 2.064600E-01_EB, 2.086828E-01_EB, &
    2.104272E-01_EB, 2.117804E-01_EB, 2.128038E-01_EB, 2.135419E-01_EB, 2.140292E-01_EB, 2.142930E-01_EB, &
    2.143560E-01_EB, 2.142378E-01_EB, 2.139561E-01_EB, 2.135262E-01_EB, 2.129625E-01_EB, 2.122782E-01_EB, &
    2.114852E-01_EB, 2.105949E-01_EB, 2.096174E-01_EB, 2.085623E-01_EB, 2.074384E-01_EB,  &
    2.568613E-01_EB, 2.704250E-01_EB, 2.811145E-01_EB, 2.897314E-01_EB, 2.967954E-01_EB, 3.026533E-01_EB, &
    3.075431E-01_EB, 3.116323E-01_EB, 3.150429E-01_EB, 3.178662E-01_EB, 3.201734E-01_EB, 3.220213E-01_EB, &
    3.234578E-01_EB, 3.245233E-01_EB, 3.252533E-01_EB, 3.256792E-01_EB, 3.258297E-01_EB, 3.257303E-01_EB, &
    3.254045E-01_EB, 3.248741E-01_EB, 3.241584E-01_EB, 3.232760E-01_EB, 3.222432E-01_EB,  &
    5.247561E-01_EB, 5.334118E-01_EB, 5.398339E-01_EB, 5.447295E-01_EB, 5.485168E-01_EB, 5.514520E-01_EB, &
    5.536955E-01_EB, 5.553509E-01_EB, 5.564881E-01_EB, 5.571563E-01_EB, 5.573943E-01_EB, 5.572332E-01_EB, &
    5.567012E-01_EB, 5.558245E-01_EB, 5.546280E-01_EB, 5.531355E-01_EB, 5.513706E-01_EB, 5.493557E-01_EB, &
    5.471128E-01_EB, 5.446627E-01_EB, 5.420252E-01_EB, 5.392194E-01_EB, 5.362628E-01_EB,  &
    1.730313E+00_EB, 1.652531E+00_EB, 1.589848E+00_EB, 1.537847E+00_EB, 1.493684E+00_EB, 1.455424E+00_EB, &
    1.421693E+00_EB, 1.391486E+00_EB, 1.364054E+00_EB, 1.338827E+00_EB, 1.315372E+00_EB, 1.293353E+00_EB, &
    1.272513E+00_EB, 1.252650E+00_EB, 1.233609E+00_EB, 1.215269E+00_EB, 1.197535E+00_EB, 1.180333E+00_EB, &
    1.163606E+00_EB, 1.147307E+00_EB, 1.131401E+00_EB, 1.115859E+00_EB, 1.100659E+00_EB,  &
    8.615218E+00_EB, 7.282297E+00_EB, 6.276256E+00_EB, 5.493467E+00_EB, 4.869153E+00_EB, 4.360881E+00_EB, &
    3.939794E+00_EB, 3.585664E+00_EB, 3.283937E+00_EB, 3.023915E+00_EB, 2.797582E+00_EB, 2.598828E+00_EB, &
    2.422926E+00_EB, 2.266167E+00_EB, 2.125606E+00_EB, 1.998873E+00_EB, 1.884046E+00_EB, 1.779543E+00_EB, &
    1.684056E+00_EB, 1.596493E+00_EB, 1.515932E+00_EB, 1.441591E+00_EB, 1.372803E+00_EB,  &
    1.155626E+01_EB, 9.980610E+00_EB, 8.692382E+00_EB, 7.637410E+00_EB, 6.766963E+00_EB, 6.041884E+00_EB, &
    5.431828E+00_EB, 4.913515E+00_EB, 4.469075E+00_EB, 4.084711E+00_EB, 3.749691E+00_EB, 3.455585E+00_EB, &
    3.195718E+00_EB, 2.964745E+00_EB, 2.758356E+00_EB, 2.573035E+00_EB, 2.405897E+00_EB, 2.254549E+00_EB, &
    2.116995E+00_EB, 1.991556E+00_EB, 1.876811E+00_EB, 1.771548E+00_EB, 1.674728E+00_EB/),(/23,8/))

SD2_C7H16(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
    1.312382E+01_EB, 1.161977E+01_EB, 1.043014E+01_EB, 9.465560E+00_EB, 8.667236E+00_EB, 7.994989E+00_EB, &
    7.420459E+00_EB, 6.923098E+00_EB, 6.487683E+00_EB, 6.102729E+00_EB, 5.759426E+00_EB, 5.450920E+00_EB, &
    5.171808E+00_EB, 4.917783E+00_EB, 4.685368E+00_EB, 4.471730E+00_EB, 4.274536E+00_EB, 4.091847E+00_EB, &
    3.922039E+00_EB, 3.763736E+00_EB, 3.615767E+00_EB, 3.477125E+00_EB, 3.346941E+00_EB,  &
    2.115810E+01_EB, 1.827380E+01_EB, 1.597284E+01_EB, 1.411137E+01_EB, 1.258368E+01_EB, 1.131280E+01_EB, &
    1.024226E+01_EB, 9.330182E+00_EB, 8.545095E+00_EB, 7.863056E+00_EB, 7.265605E+00_EB, 6.738344E+00_EB, &
    6.269909E+00_EB, 5.851234E+00_EB, 5.475010E+00_EB, 5.135288E+00_EB, 4.827178E+00_EB, 4.546628E+00_EB, &
    4.290251E+00_EB, 4.055195E+00_EB, 3.839041E+00_EB, 3.639722E+00_EB, 3.455461E+00_EB,  &
    2.663189E+01_EB, 2.210694E+01_EB, 1.866232E+01_EB, 1.598051E+01_EB, 1.385104E+01_EB, 1.213070E+01_EB, &
    1.071961E+01_EB, 9.546548E+00_EB, 8.559696E+00_EB, 7.720674E+00_EB, 7.000583E+00_EB, 6.377341E+00_EB, &
    5.833839E+00_EB, 5.356656E+00_EB, 4.935137E+00_EB, 4.560735E+00_EB, 4.226520E+00_EB, 3.926826E+00_EB, &
    3.656976E+00_EB, 3.413077E+00_EB, 3.191871E+00_EB, 2.990605E+00_EB, 2.806941E+00_EB,  &
    8.943802E+00_EB, 7.800001E+00_EB, 6.853228E+00_EB, 6.065532E+00_EB, 5.405368E+00_EB, 4.847590E+00_EB, &
    4.372443E+00_EB, 3.964450E+00_EB, 3.611462E+00_EB, 3.303895E+00_EB, 3.034140E+00_EB, 2.796114E+00_EB, &
    2.584915E+00_EB, 2.396562E+00_EB, 2.227800E+00_EB, 2.075943E+00_EB, 1.938762E+00_EB, 1.814393E+00_EB, &
    1.701266E+00_EB, 1.598048E+00_EB, 1.503607E+00_EB, 1.416969E+00_EB, 1.337295E+00_EB,  &
    7.413487E-01_EB, 8.328178E-01_EB, 8.940677E-01_EB, 9.328753E-01_EB, 9.553320E-01_EB, 9.659550E-01_EB, &
    9.680242E-01_EB, 9.639073E-01_EB, 9.553161E-01_EB, 9.434951E-01_EB, 9.293561E-01_EB, 9.135727E-01_EB, &
    8.966474E-01_EB, 8.789582E-01_EB, 8.607912E-01_EB, 8.423651E-01_EB, 8.238472E-01_EB, 8.053662E-01_EB, &
    7.870213E-01_EB, 7.688885E-01_EB, 7.510261E-01_EB, 7.334782E-01_EB, 7.162776E-01_EB,  &
    4.447855E-01_EB, 4.761631E-01_EB, 4.981635E-01_EB, 5.135837E-01_EB, 5.242947E-01_EB, 5.315653E-01_EB, &
    5.362697E-01_EB, 5.390205E-01_EB, 5.402545E-01_EB, 5.402897E-01_EB, 5.393619E-01_EB, 5.376499E-01_EB, &
    5.352919E-01_EB, 5.323974E-01_EB, 5.290544E-01_EB, 5.253356E-01_EB, 5.213008E-01_EB, 5.170013E-01_EB, &
    5.124806E-01_EB, 5.077761E-01_EB, 5.029202E-01_EB, 4.979416E-01_EB, 4.928649E-01_EB,  &
    2.315471E-01_EB, 2.576360E-01_EB, 2.790490E-01_EB, 2.968662E-01_EB, 3.118735E-01_EB, 3.246434E-01_EB, &
    3.355972E-01_EB, 3.450489E-01_EB, 3.532361E-01_EB, 3.603416E-01_EB, 3.665083E-01_EB, 3.718506E-01_EB, &
    3.764623E-01_EB, 3.804211E-01_EB, 3.837930E-01_EB, 3.866349E-01_EB, 3.889966E-01_EB, 3.909220E-01_EB, &
    3.924504E-01_EB, 3.936170E-01_EB, 3.944539E-01_EB, 3.949902E-01_EB, 3.952524E-01_EB,  &
    1.444890E-01_EB, 1.683858E-01_EB, 1.889547E-01_EB, 2.067257E-01_EB, 2.221632E-01_EB, 2.356497E-01_EB, &
    2.474920E-01_EB, 2.579338E-01_EB, 2.671698E-01_EB, 2.753558E-01_EB, 2.826186E-01_EB, 2.890625E-01_EB, &
    2.947742E-01_EB, 2.998275E-01_EB, 3.042856E-01_EB, 3.082035E-01_EB, 3.116299E-01_EB, 3.146075E-01_EB, &
    3.171752E-01_EB, 3.193676E-01_EB, 3.212162E-01_EB, 3.227497E-01_EB, 3.239944E-01_EB/),(/23,8/))

SD2_C7H16(1:23,17:21) = RESHAPE((/ &  ! 3100-3200 cm-1
    9.033260E-02_EB, 1.105435E-01_EB, 1.286099E-01_EB, 1.446695E-01_EB, 1.589357E-01_EB, 1.716272E-01_EB, &
    1.829430E-01_EB, 1.930549E-01_EB, 2.021086E-01_EB, 2.102262E-01_EB, 2.175102E-01_EB, 2.240476E-01_EB, &
    2.299125E-01_EB, 2.351689E-01_EB, 2.398720E-01_EB, 2.440710E-01_EB, 2.478094E-01_EB, 2.511258E-01_EB, &
    2.540554E-01_EB, 2.566296E-01_EB, 2.588773E-01_EB, 2.608248E-01_EB, 2.624958E-01_EB,  &
    7.157048E-02_EB, 8.845027E-02_EB, 1.036743E-01_EB, 1.173022E-01_EB, 1.294767E-01_EB, 1.403590E-01_EB, &
    1.501016E-01_EB, 1.588403E-01_EB, 1.666915E-01_EB, 1.737542E-01_EB, 1.801124E-01_EB, 1.858377E-01_EB, &
    1.909916E-01_EB, 1.956271E-01_EB, 1.997908E-01_EB, 2.035235E-01_EB, 2.068619E-01_EB, 2.098384E-01_EB, &
    2.124827E-01_EB, 2.148213E-01_EB, 2.168786E-01_EB, 2.186767E-01_EB, 2.202362E-01_EB,  &
    8.822796E-02_EB, 9.804427E-02_EB, 1.059969E-01_EB, 1.124904E-01_EB, 1.178339E-01_EB, 1.222607E-01_EB, &
    1.259470E-01_EB, 1.290258E-01_EB, 1.315998E-01_EB, 1.337480E-01_EB, 1.355327E-01_EB, 1.370035E-01_EB, &
    1.382010E-01_EB, 1.391582E-01_EB, 1.399033E-01_EB, 1.404598E-01_EB, 1.408483E-01_EB, 1.410864E-01_EB, &
    1.411899E-01_EB, 1.411728E-01_EB, 1.410473E-01_EB, 1.408244E-01_EB, 1.405143E-01_EB,  &
    1.071646E-01_EB, 1.039738E-01_EB, 1.010031E-01_EB, 9.825636E-02_EB, 9.572028E-02_EB, 9.337487E-02_EB, &
    9.119844E-02_EB, 8.916983E-02_EB, 8.726975E-02_EB, 8.548107E-02_EB, 8.378892E-02_EB, 8.218057E-02_EB, &
    8.064545E-02_EB, 7.917453E-02_EB, 7.776042E-02_EB, 7.639694E-02_EB, 7.507907E-02_EB, 7.380249E-02_EB, &
    7.256381E-02_EB, 7.136007E-02_EB, 7.018892E-02_EB, 6.904833E-02_EB, 6.793654E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

!-------------------------Methanol DATA-------------------


! THERE ARE 3 BANDS FOR Methanol

! BAND #1: 2700 cm-1 - 3200 cm-1 
! BAND #2: 3600 cm-1 - 3800 cm-1 
! BAND #3: 900 cm-1 - 1600 cm-1 

N_TEMP_CH3OH = 23
N_BAND_CH3OH = 3

ALLOCATE(SD_CH3OH_TEMP(N_TEMP_CH3OH)) 

! INITIALIZE BANDS WAVENUMBER BOUNDS FOR Methanol ABSOPTION COEFFICIENTS.
! BANDS ARE ORGANIZED BY ROW: 
! 1ST COLUMN IS THE LOWER BOUND BAND LIMIT IN cm-1. 
! 2ND COLUMN IS THE UPPER BOUND BAND LIMIT IN cm-1. 
! 3RD COLUMN IS THE STRIDE BETWEEN WAVENUMBERS IN BAND.
! IF 3RD COLUMN = 0., THEN THE BAND IS CALCULATED AND NOT TABULATED.

ALLOCATE(OM_BND_CH3OH(N_BAND_CH3OH,3)) 

OM_BND_CH3OH = RESHAPE((/ &
    2700._EB, 3600._EB,  900._EB, &
    3200._EB, 3800._EB, 1600._EB, &
     25._EB, 25._EB, 25._EB/),(/N_BAND_CH3OH,3/)) 

SD_CH3OH_TEMP = (/ &
    300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
    600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
    900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
    1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(SD1_CH3OH(N_TEMP_CH3OH,21)) 

! BAND #1: 2700 cm-1 - 3200 cm-1 

SD1_CH3OH(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
    6.520051E-03_EB, 1.243656E-02_EB, 1.947415E-02_EB, 2.684155E-02_EB, 3.392758E-02_EB, 4.034511E-02_EB, &
    4.589468E-02_EB, 5.050938E-02_EB, 5.420601E-02_EB, 5.704902E-02_EB, 5.912630E-02_EB, 6.053375E-02_EB, &
    6.136636E-02_EB, 6.171313E-02_EB, 6.165472E-02_EB, 6.126252E-02_EB, 6.059871E-02_EB, 5.971689E-02_EB, &
    5.866272E-02_EB, 5.747485E-02_EB, 5.618581E-02_EB, 5.482279E-02_EB, 5.340843E-02_EB,  &
    1.557141E-02_EB, 2.712243E-02_EB, 4.009545E-02_EB, 5.316125E-02_EB, 6.537810E-02_EB, 7.619762E-02_EB, &
    8.537463E-02_EB, 9.286623E-02_EB, 9.875087E-02_EB, 1.031710E-01_EB, 1.062967E-01_EB, 1.083033E-01_EB, &
    1.093589E-01_EB, 1.096172E-01_EB, 1.092161E-01_EB, 1.082758E-01_EB, 1.069002E-01_EB, 1.051782E-01_EB, &
    1.031851E-01_EB, 1.009844E-01_EB, 9.862918E-02_EB, 9.616377E-02_EB, 9.362477E-02_EB,  &
    6.963432E-02_EB, 9.533441E-02_EB, 1.170516E-01_EB, 1.341145E-01_EB, 1.467519E-01_EB, 1.555542E-01_EB, &
    1.611947E-01_EB, 1.643048E-01_EB, 1.654296E-01_EB, 1.650198E-01_EB, 1.634394E-01_EB, 1.609787E-01_EB, &
    1.578671E-01_EB, 1.542854E-01_EB, 1.503750E-01_EB, 1.462468E-01_EB, 1.419870E-01_EB, 1.376624E-01_EB, &
    1.333250E-01_EB, 1.290147E-01_EB, 1.247614E-01_EB, 1.205877E-01_EB, 1.165105E-01_EB,  &
    2.786040E-01_EB, 3.234636E-01_EB, 3.525565E-01_EB, 3.693698E-01_EB, 3.771135E-01_EB, 3.783575E-01_EB, &
    3.750405E-01_EB, 3.685841E-01_EB, 3.600170E-01_EB, 3.500789E-01_EB, 3.393000E-01_EB, 3.280599E-01_EB, &
    3.166297E-01_EB, 3.052019E-01_EB, 2.939127E-01_EB, 2.828566E-01_EB, 2.720982E-01_EB, 2.616800E-01_EB, &
    2.516280E-01_EB, 2.419569E-01_EB, 2.326724E-01_EB, 2.237736E-01_EB, 2.152555E-01_EB,  &
    1.135988E+00_EB, 1.054777E+00_EB, 9.681041E-01_EB, 8.841707E-01_EB, 8.063652E-01_EB, 7.357579E-01_EB, &
    6.723391E-01_EB, 6.156292E-01_EB, 5.649786E-01_EB, 5.197096E-01_EB, 4.791811E-01_EB, 4.428136E-01_EB, &
    4.100964E-01_EB, 3.805843E-01_EB, 3.538919E-01_EB, 3.296863E-01_EB, 3.076804E-01_EB, 2.876254E-01_EB, &
    2.693060E-01_EB, 2.525353E-01_EB, 2.371499E-01_EB, 2.230073E-01_EB, 2.099827E-01_EB,  &
    1.733559E+00_EB, 1.480358E+00_EB, 1.285263E+00_EB, 1.131100E+00_EB, 1.006648E+00_EB, 9.043272E-01_EB, &
    8.188643E-01_EB, 7.464943E-01_EB, 6.844709E-01_EB, 6.307508E-01_EB, 5.837878E-01_EB, 5.423931E-01_EB, &
    5.056397E-01_EB, 4.727958E-01_EB, 4.432759E-01_EB, 4.166068E-01_EB, 3.924018E-01_EB, 3.703424E-01_EB, &
    3.501629E-01_EB, 3.316403E-01_EB, 3.145864E-01_EB, 2.988404E-01_EB, 2.842644E-01_EB,  &
    2.211436E+00_EB, 1.869525E+00_EB, 1.602463E+00_EB, 1.390638E+00_EB, 1.219944E+00_EB, 1.080319E+00_EB, &
    9.645271E-01_EB, 8.672978E-01_EB, 7.847365E-01_EB, 7.139231E-01_EB, 6.526375E-01_EB, 5.991687E-01_EB, &
    5.521817E-01_EB, 5.106224E-01_EB, 4.736489E-01_EB, 4.405817E-01_EB, 4.108675E-01_EB, 3.840509E-01_EB, &
    3.597543E-01_EB, 3.376621E-01_EB, 3.175083E-01_EB, 2.990672E-01_EB, 2.821467E-01_EB,  &
    1.954717E+00_EB, 1.797572E+00_EB, 1.655178E+00_EB, 1.528669E+00_EB, 1.416917E+00_EB, 1.318152E+00_EB, &
    1.230568E+00_EB, 1.152533E+00_EB, 1.082642E+00_EB, 1.019713E+00_EB, 9.627633E-01_EB, 9.109756E-01_EB, &
    8.636719E-01_EB, 8.202866E-01_EB, 7.803470E-01_EB, 7.434545E-01_EB, 7.092719E-01_EB, 6.775118E-01_EB, &
    6.479276E-01_EB, 6.203066E-01_EB, 5.944643E-01_EB, 5.702395E-01_EB, 5.474905E-01_EB/),(/23,8/))

SD1_CH3OH(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
    3.044471E+00_EB, 2.596990E+00_EB, 2.262836E+00_EB, 2.003917E+00_EB, 1.797420E+00_EB, 1.628854E+00_EB, &
    1.488576E+00_EB, 1.369929E+00_EB, 1.268176E+00_EB, 1.179860E+00_EB, 1.102405E+00_EB, 1.033857E+00_EB, &
    9.727041E-01_EB, 9.177657E-01_EB, 8.681035E-01_EB, 8.229638E-01_EB, 7.817344E-01_EB, 7.439126E-01_EB, &
    7.090818E-01_EB, 6.768933E-01_EB, 6.470530E-01_EB, 6.193111E-01_EB, 5.934532E-01_EB,  &
    3.779874E+00_EB, 3.098253E+00_EB, 2.608118E+00_EB, 2.240593E+00_EB, 1.955874E+00_EB, 1.729453E+00_EB, &
    1.545475E+00_EB, 1.393257E+00_EB, 1.265369E+00_EB, 1.156492E+00_EB, 1.062733E+00_EB, 9.811803E-01_EB, &
    9.096211E-01_EB, 8.463424E-01_EB, 7.900018E-01_EB, 7.395322E-01_EB, 6.940758E-01_EB, 6.529347E-01_EB, &
    6.155362E-01_EB, 5.814051E-01_EB, 5.501443E-01_EB, 5.214191E-01_EB, 4.949452E-01_EB,  &
    4.023285E+00_EB, 3.261083E+00_EB, 2.718464E+00_EB, 2.315174E+00_EB, 2.005228E+00_EB, 1.760531E+00_EB, &
    1.563031E+00_EB, 1.400644E+00_EB, 1.265007E+00_EB, 1.150171E+00_EB, 1.051796E+00_EB, 9.666545E-01_EB, &
    8.923004E-01_EB, 8.268490E-01_EB, 7.688272E-01_EB, 7.170684E-01_EB, 6.706375E-01_EB, 6.287761E-01_EB, &
    5.908638E-01_EB, 5.563871E-01_EB, 5.249183E-01_EB, 4.960977E-01_EB, 4.696209E-01_EB,  &
    3.332680E+00_EB, 2.816818E+00_EB, 2.431885E+00_EB, 2.134438E+00_EB, 1.898145E+00_EB, 1.706151E+00_EB, &
    1.547190E+00_EB, 1.413465E+00_EB, 1.299418E+00_EB, 1.200992E+00_EB, 1.115163E+00_EB, 1.039635E+00_EB, &
    9.726377E-01_EB, 9.127841E-01_EB, 8.589755E-01_EB, 8.103300E-01_EB, 7.661318E-01_EB, 7.257943E-01_EB, &
    6.888319E-01_EB, 6.548391E-01_EB, 6.234745E-01_EB, 5.944485E-01_EB, 5.675137E-01_EB,  &
    2.049144E+00_EB, 1.783895E+00_EB, 1.581092E+00_EB, 1.420837E+00_EB, 1.290872E+00_EB, 1.183225E+00_EB, &
    1.092483E+00_EB, 1.014841E+00_EB, 9.475525E-01_EB, 8.885837E-01_EB, 8.364018E-01_EB, 7.898304E-01_EB, &
    7.479535E-01_EB, 7.100480E-01_EB, 6.755359E-01_EB, 6.439504E-01_EB, 6.149097E-01_EB, 5.880990E-01_EB, &
    5.632567E-01_EB, 5.401623E-01_EB, 5.186296E-01_EB, 4.984996E-01_EB, 4.796354E-01_EB,  &
    1.106340E+00_EB, 1.014955E+00_EB, 9.398530E-01_EB, 8.769466E-01_EB, 8.233860E-01_EB, 7.771294E-01_EB, &
    7.366781E-01_EB, 7.009097E-01_EB, 6.689689E-01_EB, 6.401938E-01_EB, 6.140669E-01_EB, 5.901785E-01_EB, &
    5.682013E-01_EB, 5.478712E-01_EB, 5.289736E-01_EB, 5.113323E-01_EB, 4.948016E-01_EB, 4.792600E-01_EB, &
    4.646055E-01_EB, 4.507513E-01_EB, 4.376235E-01_EB, 4.251585E-01_EB, 4.133013E-01_EB,  &
    4.420023E-01_EB, 4.518956E-01_EB, 4.577579E-01_EB, 4.610461E-01_EB, 4.626276E-01_EB, 4.630366E-01_EB, &
    4.626111E-01_EB, 4.615700E-01_EB, 4.600581E-01_EB, 4.581733E-01_EB, 4.559843E-01_EB, 4.535404E-01_EB, &
    4.508791E-01_EB, 4.480298E-01_EB, 4.450163E-01_EB, 4.418593E-01_EB, 4.385768E-01_EB, 4.351848E-01_EB, &
    4.316977E-01_EB, 4.281293E-01_EB, 4.244918E-01_EB, 4.207964E-01_EB, 4.170537E-01_EB,  &
    1.747783E-01_EB, 1.869237E-01_EB, 1.965984E-01_EB, 2.044722E-01_EB, 2.109911E-01_EB, 2.164603E-01_EB, &
    2.210935E-01_EB, 2.250433E-01_EB, 2.284219E-01_EB, 2.313124E-01_EB, 2.337791E-01_EB, 2.358723E-01_EB, &
    2.376327E-01_EB, 2.390939E-01_EB, 2.402847E-01_EB, 2.412299E-01_EB, 2.419512E-01_EB, 2.424682E-01_EB, &
    2.427984E-01_EB, 2.429580E-01_EB, 2.429615E-01_EB, 2.428225E-01_EB, 2.425533E-01_EB/),(/23,8/))

SD1_CH3OH(1:23,17:21) = RESHAPE((/ &  ! 3100-3200 cm-1
    3.766400E-02_EB, 5.194706E-02_EB, 6.625533E-02_EB, 8.014988E-02_EB, 9.339578E-02_EB, 1.058819E-01_EB, &
    1.175689E-01_EB, 1.284580E-01_EB, 1.385725E-01_EB, 1.479462E-01_EB, 1.566179E-01_EB, 1.646284E-01_EB, &
    1.720176E-01_EB, 1.788244E-01_EB, 1.850859E-01_EB, 1.908374E-01_EB, 1.961121E-01_EB, 2.009414E-01_EB, &
    2.053548E-01_EB, 2.093800E-01_EB, 2.130429E-01_EB, 2.163682E-01_EB, 2.193784E-01_EB,  &
    1.454342E-02_EB, 2.067111E-02_EB, 2.688575E-02_EB, 3.296275E-02_EB, 3.877772E-02_EB, 4.426843E-02_EB, &
    4.940966E-02_EB, 5.419740E-02_EB, 5.863940E-02_EB, 6.274961E-02_EB, 6.654497E-02_EB, 7.004344E-02_EB, &
    7.326318E-02_EB, 7.622191E-02_EB, 7.893657E-02_EB, 8.142335E-02_EB, 8.369744E-02_EB, 8.577327E-02_EB, &
    8.766419E-02_EB, 8.938301E-02_EB, 9.094142E-02_EB, 9.235064E-02_EB, 9.362097E-02_EB,  &
    8.267643E-03_EB, 1.154703E-02_EB, 1.478935E-02_EB, 1.788560E-02_EB, 2.078293E-02_EB, 2.346120E-02_EB, &
    2.591840E-02_EB, 2.816195E-02_EB, 3.020379E-02_EB, 3.205759E-02_EB, 3.373742E-02_EB, 3.525683E-02_EB, &
    3.662865E-02_EB, 3.786480E-02_EB, 3.897626E-02_EB, 3.997312E-02_EB, 4.086467E-02_EB, 4.165940E-02_EB, &
    4.236514E-02_EB, 4.298904E-02_EB, 4.353770E-02_EB, 4.401717E-02_EB, 4.443306E-02_EB,  &
    6.233956E-03_EB, 8.676358E-03_EB, 1.107788E-02_EB, 1.335944E-02_EB, 1.548419E-02_EB, 1.743944E-02_EB, &
    1.922564E-02_EB, 2.084992E-02_EB, 2.232243E-02_EB, 2.365434E-02_EB, 2.485680E-02_EB, 2.594054E-02_EB, &
    2.691548E-02_EB, 2.779079E-02_EB, 2.857488E-02_EB, 2.927538E-02_EB, 2.989932E-02_EB, 3.045304E-02_EB, &
    3.094241E-02_EB, 3.137273E-02_EB, 3.174889E-02_EB, 3.207535E-02_EB, 3.235624E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

ALLOCATE(SD2_CH3OH(N_TEMP_CH3OH,9)) 

! BAND #2: 3600 cm-1 - 3800 cm-1 

SD2_CH3OH(1:23,1:8) = RESHAPE((/ &  ! 3600-3775 cm-1
    6.783104E-02_EB, 9.256382E-02_EB, 1.140839E-01_EB, 1.316075E-01_EB, 1.451522E-01_EB, 1.551139E-01_EB, &
    1.620218E-01_EB, 1.664123E-01_EB, 1.687725E-01_EB, 1.695216E-01_EB, 1.690074E-01_EB, 1.675149E-01_EB, &
    1.652733E-01_EB, 1.624665E-01_EB, 1.592413E-01_EB, 1.557137E-01_EB, 1.519765E-01_EB, 1.481025E-01_EB, &
    1.441495E-01_EB, 1.401629E-01_EB, 1.361782E-01_EB, 1.322232E-01_EB, 1.283190E-01_EB,  &
    3.983907E-01_EB, 3.926671E-01_EB, 3.826467E-01_EB, 3.705917E-01_EB, 3.577605E-01_EB, 3.448526E-01_EB, &
    3.322492E-01_EB, 3.201475E-01_EB, 3.086379E-01_EB, 2.977493E-01_EB, 2.874754E-01_EB, 2.777913E-01_EB, &
    2.686627E-01_EB, 2.600518E-01_EB, 2.519199E-01_EB, 2.442303E-01_EB, 2.369480E-01_EB, 2.300413E-01_EB, &
    2.234809E-01_EB, 2.172401E-01_EB, 2.112953E-01_EB, 2.056249E-01_EB, 2.002091E-01_EB,  &
    9.858871E-01_EB, 8.057830E-01_EB, 6.712477E-01_EB, 5.682437E-01_EB, 4.876653E-01_EB, 4.234373E-01_EB, &
    3.713961E-01_EB, 3.286178E-01_EB, 2.930030E-01_EB, 2.630145E-01_EB, 2.375065E-01_EB, 2.156116E-01_EB, &
    1.966636E-01_EB, 1.801440E-01_EB, 1.656447E-01_EB, 1.528409E-01_EB, 1.414713E-01_EB, 1.313240E-01_EB, &
    1.222256E-01_EB, 1.140325E-01_EB, 1.066262E-01_EB, 9.990666E-02_EB, 9.378992E-02_EB,  &
    1.008444E+00_EB, 7.974127E-01_EB, 6.496632E-01_EB, 5.417166E-01_EB, 4.601605E-01_EB, 3.968486E-01_EB, &
    3.465842E-01_EB, 3.059176E-01_EB, 2.724802E-01_EB, 2.446002E-01_EB, 2.210680E-01_EB, 2.009906E-01_EB, &
    1.836965E-01_EB, 1.686722E-01_EB, 1.555198E-01_EB, 1.439264E-01_EB, 1.336436E-01_EB, 1.244716E-01_EB, &
    1.162485E-01_EB, 1.088415E-01_EB, 1.021413E-01_EB, 9.605647E-02_EB, 9.051065E-02_EB,  &
    1.022169E+00_EB, 8.346469E-01_EB, 6.996219E-01_EB, 5.984777E-01_EB, 5.203170E-01_EB, 4.583768E-01_EB, &
    4.082580E-01_EB, 3.669867E-01_EB, 3.324875E-01_EB, 3.032720E-01_EB, 2.782482E-01_EB, 2.565988E-01_EB, &
    2.377012E-01_EB, 2.210740E-01_EB, 2.063398E-01_EB, 1.931987E-01_EB, 1.814106E-01_EB, 1.707804E-01_EB, &
    1.611485E-01_EB, 1.523832E-01_EB, 1.443747E-01_EB, 1.370313E-01_EB, 1.302751E-01_EB,  &
    2.632782E-01_EB, 2.607320E-01_EB, 2.572947E-01_EB, 2.535051E-01_EB, 2.496331E-01_EB, 2.458146E-01_EB, &
    2.421152E-01_EB, 2.385628E-01_EB, 2.351646E-01_EB, 2.319167E-01_EB, 2.288096E-01_EB, 2.258313E-01_EB, &
    2.229692E-01_EB, 2.202108E-01_EB, 2.175445E-01_EB, 2.149603E-01_EB, 2.124488E-01_EB, 2.100023E-01_EB, &
    2.076139E-01_EB, 2.052781E-01_EB, 2.029901E-01_EB, 2.007457E-01_EB, 1.985417E-01_EB,  &
    4.815390E-02_EB, 5.784827E-02_EB, 6.639819E-02_EB, 7.391181E-02_EB, 8.051877E-02_EB, 8.634311E-02_EB, &
    9.149478E-02_EB, 9.606700E-02_EB, 1.001381E-01_EB, 1.037730E-01_EB, 1.070256E-01_EB, 1.099406E-01_EB, &
    1.125555E-01_EB, 1.149020E-01_EB, 1.170067E-01_EB, 1.188926E-01_EB, 1.205797E-01_EB, 1.220850E-01_EB, &
    1.234242E-01_EB, 1.246106E-01_EB, 1.256565E-01_EB, 1.265729E-01_EB, 1.273696E-01_EB,  &
    1.792447E-02_EB, 2.373545E-02_EB, 2.921848E-02_EB, 3.427127E-02_EB, 3.886922E-02_EB, 4.302601E-02_EB, &
    4.677221E-02_EB, 5.014425E-02_EB, 5.317885E-02_EB, 5.591048E-02_EB, 5.837031E-02_EB, 6.058602E-02_EB, &
    6.258196E-02_EB, 6.437945E-02_EB, 6.599726E-02_EB, 6.745177E-02_EB, 6.875751E-02_EB, 6.992743E-02_EB, &
    7.097285E-02_EB, 7.190426E-02_EB, 7.273072E-02_EB, 7.346071E-02_EB, 7.410173E-02_EB/),(/23,8/))

SD2_CH3OH(1:23,9:9) = RESHAPE((/ &  ! 3800-3800 cm-1
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,1/))

ALLOCATE(SD3_CH3OH(N_TEMP_CH3OH,29)) 

! BAND #3: 900 cm-1 - 1600 cm-1 

SD3_CH3OH(1:23,1:8) = RESHAPE((/ &  ! 900-1075 cm-1
    4.645956E-03_EB, 1.418512E-02_EB, 3.248177E-02_EB, 6.138079E-02_EB, 1.013960E-01_EB, 1.519049E-01_EB, &
    2.115228E-01_EB, 2.784715E-01_EB, 3.508624E-01_EB, 4.268853E-01_EB, 5.049135E-01_EB, 5.835534E-01_EB, &
    6.616546E-01_EB, 7.382988E-01_EB, 8.127755E-01_EB, 8.845549E-01_EB, 9.532590E-01_EB, 1.018635E+00_EB, &
    1.080530E+00_EB, 1.138874E+00_EB, 1.193657E+00_EB, 1.244919E+00_EB, 1.292736E+00_EB,  &
    3.137794E-02_EB, 7.210298E-02_EB, 1.312594E-01_EB, 2.049635E-01_EB, 2.878652E-01_EB, 3.747190E-01_EB, &
    4.611764E-01_EB, 5.440240E-01_EB, 6.211177E-01_EB, 6.911955E-01_EB, 7.536651E-01_EB, 8.084111E-01_EB, &
    8.556391E-01_EB, 8.957567E-01_EB, 9.292866E-01_EB, 9.568064E-01_EB, 9.789080E-01_EB, 9.961696E-01_EB, &
    1.009140E+00_EB, 1.018331E+00_EB, 1.024210E+00_EB, 1.027201E+00_EB, 1.027685E+00_EB,  &
    3.109893E-01_EB, 4.812788E-01_EB, 6.536835E-01_EB, 8.135750E-01_EB, 9.525941E-01_EB, 1.067266E+00_EB, &
    1.157352E+00_EB, 1.224494E+00_EB, 1.271266E+00_EB, 1.300562E+00_EB, 1.315245E+00_EB, 1.317958E+00_EB, &
    1.311042E+00_EB, 1.296510E+00_EB, 1.276063E+00_EB, 1.251118E+00_EB, 1.222838E+00_EB, 1.192172E+00_EB, &
    1.159890E+00_EB, 1.126606E+00_EB, 1.092812E+00_EB, 1.058894E+00_EB, 1.025157E+00_EB,  &
    2.584314E+00_EB, 2.600417E+00_EB, 2.571456E+00_EB, 2.509950E+00_EB, 2.426612E+00_EB, 2.329850E+00_EB, &
    2.225905E+00_EB, 2.119235E+00_EB, 2.012918E+00_EB, 1.909008E+00_EB, 1.808814E+00_EB, 1.713120E+00_EB, &
    1.622340E+00_EB, 1.536635E+00_EB, 1.455995E+00_EB, 1.380299E+00_EB, 1.309357E+00_EB, 1.242934E+00_EB, &
    1.180776E+00_EB, 1.122621E+00_EB, 1.068207E+00_EB, 1.017278E+00_EB, 9.695923E-01_EB,  &
    5.746240E+00_EB, 4.687148E+00_EB, 3.915379E+00_EB, 3.329242E+00_EB, 2.870242E+00_EB, 2.502278E+00_EB, &
    2.201789E+00_EB, 1.952695E+00_EB, 1.743616E+00_EB, 1.566260E+00_EB, 1.414435E+00_EB, 1.283424E+00_EB, &
    1.169569E+00_EB, 1.069996E+00_EB, 9.824124E-01_EB, 9.049729E-01_EB, 8.361757E-01_EB, 7.747880E-01_EB, &
    7.197896E-01_EB, 6.703303E-01_EB, 6.256965E-01_EB, 5.852857E-01_EB, 5.485868E-01_EB,  &
    5.858720E+00_EB, 4.884128E+00_EB, 4.165995E+00_EB, 3.613983E+00_EB, 3.176190E+00_EB, 2.820619E+00_EB, &
    2.526390E+00_EB, 2.279242E+00_EB, 2.069058E+00_EB, 1.888440E+00_EB, 1.731842E+00_EB, 1.595013E+00_EB, &
    1.474639E+00_EB, 1.368095E+00_EB, 1.273277E+00_EB, 1.188476E+00_EB, 1.112291E+00_EB, 1.043563E+00_EB, &
    9.813267E-01_EB, 9.247701E-01_EB, 8.732066E-01_EB, 8.260522E-01_EB, 7.828073E-01_EB,  &
    6.713319E+00_EB, 5.768662E+00_EB, 5.042329E+00_EB, 4.459659E+00_EB, 3.978760E+00_EB, 3.573970E+00_EB, &
    3.228368E+00_EB, 2.930123E+00_EB, 2.670562E+00_EB, 2.443099E+00_EB, 2.242592E+00_EB, 2.064942E+00_EB, &
    1.906826E+00_EB, 1.765520E+00_EB, 1.638759E+00_EB, 1.524653E+00_EB, 1.421606E+00_EB, 1.328266E+00_EB, &
    1.243477E+00_EB, 1.166252E+00_EB, 1.095736E+00_EB, 1.031193E+00_EB, 9.719829E-01_EB,  &
    5.041451E-01_EB, 6.320024E-01_EB, 7.344743E-01_EB, 8.120779E-01_EB, 8.676687E-01_EB, 9.048069E-01_EB, &
    9.269851E-01_EB, 9.373104E-01_EB, 9.384098E-01_EB, 9.324392E-01_EB, 9.211357E-01_EB, 9.058814E-01_EB, &
    8.877668E-01_EB, 8.676460E-01_EB, 8.461848E-01_EB, 8.238989E-01_EB, 8.011855E-01_EB, 7.783490E-01_EB, &
    7.556201E-01_EB, 7.331723E-01_EB, 7.111336E-01_EB, 6.895973E-01_EB, 6.686290E-01_EB/),(/23,8/))

SD3_CH3OH(1:23,9:16) = RESHAPE((/ &  ! 1100-1275 cm-1
    1.570875E-01_EB, 1.740335E-01_EB, 1.855464E-01_EB, 1.928879E-01_EB, 1.970644E-01_EB, 1.988585E-01_EB, &
    1.988724E-01_EB, 1.975670E-01_EB, 1.952945E-01_EB, 1.923232E-01_EB, 1.888578E-01_EB, 1.850542E-01_EB, &
    1.810309E-01_EB, 1.768779E-01_EB, 1.726632E-01_EB, 1.684380E-01_EB, 1.642405E-01_EB, 1.600993E-01_EB, &
    1.560347E-01_EB, 1.520616E-01_EB, 1.481902E-01_EB, 1.444274E-01_EB, 1.407769E-01_EB,  &
    1.022431E-01_EB, 1.192492E-01_EB, 1.334763E-01_EB, 1.451979E-01_EB, 1.547220E-01_EB, 1.623515E-01_EB, &
    1.683657E-01_EB, 1.730124E-01_EB, 1.765081E-01_EB, 1.790383E-01_EB, 1.807614E-01_EB, 1.818107E-01_EB, &
    1.822996E-01_EB, 1.823226E-01_EB, 1.819595E-01_EB, 1.812774E-01_EB, 1.803325E-01_EB, 1.791720E-01_EB, &
    1.778356E-01_EB, 1.763569E-01_EB, 1.747636E-01_EB, 1.730797E-01_EB, 1.713250E-01_EB,  &
    7.332982E-02_EB, 9.818605E-02_EB, 1.218422E-01_EB, 1.435984E-01_EB, 1.631485E-01_EB, 1.804240E-01_EB, &
    1.954935E-01_EB, 2.084996E-01_EB, 2.196200E-01_EB, 2.290429E-01_EB, 2.369539E-01_EB, 2.435281E-01_EB, &
    2.489265E-01_EB, 2.532952E-01_EB, 2.567640E-01_EB, 2.594481E-01_EB, 2.614494E-01_EB, 2.628571E-01_EB, &
    2.637492E-01_EB, 2.641946E-01_EB, 2.642529E-01_EB, 2.639759E-01_EB, 2.634096E-01_EB,  &
    1.068325E-01_EB, 1.443289E-01_EB, 1.797892E-01_EB, 2.120736E-01_EB, 2.407192E-01_EB, 2.656619E-01_EB, &
    2.870594E-01_EB, 3.051829E-01_EB, 3.203508E-01_EB, 3.328907E-01_EB, 3.431186E-01_EB, 3.513271E-01_EB, &
    3.577816E-01_EB, 3.627184E-01_EB, 3.663462E-01_EB, 3.688477E-01_EB, 3.703823E-01_EB, 3.710883E-01_EB, &
    3.710857E-01_EB, 3.704784E-01_EB, 3.693561E-01_EB, 3.677964E-01_EB, 3.658665E-01_EB,  &
    2.021620E-01_EB, 2.477076E-01_EB, 2.803323E-01_EB, 3.014587E-01_EB, 3.132403E-01_EB, 3.178159E-01_EB, &
    3.170349E-01_EB, 3.123936E-01_EB, 3.050578E-01_EB, 2.959151E-01_EB, 2.856294E-01_EB, 2.746910E-01_EB, &
    2.634568E-01_EB, 2.521836E-01_EB, 2.410523E-01_EB, 2.301884E-01_EB, 2.196754E-01_EB, 2.095662E-01_EB, &
    1.998916E-01_EB, 1.906660E-01_EB, 1.818925E-01_EB, 1.735654E-01_EB, 1.656739E-01_EB,  &
    2.437475E-01_EB, 2.875029E-01_EB, 3.175354E-01_EB, 3.362489E-01_EB, 3.461446E-01_EB, 3.494020E-01_EB, &
    3.477872E-01_EB, 3.426782E-01_EB, 3.351286E-01_EB, 3.259322E-01_EB, 3.156817E-01_EB, 3.048152E-01_EB, &
    2.936537E-01_EB, 2.824305E-01_EB, 2.713127E-01_EB, 2.604174E-01_EB, 2.498250E-01_EB, 2.395889E-01_EB, &
    2.297412E-01_EB, 2.203000E-01_EB, 2.112721E-01_EB, 2.026561E-01_EB, 1.944455E-01_EB,  &
    3.164822E-01_EB, 3.637627E-01_EB, 3.931615E-01_EB, 4.086715E-01_EB, 4.139133E-01_EB, 4.118083E-01_EB, &
    4.045957E-01_EB, 3.939455E-01_EB, 3.810833E-01_EB, 3.668974E-01_EB, 3.520250E-01_EB, 3.369179E-01_EB, &
    3.218918E-01_EB, 3.071633E-01_EB, 2.928765E-01_EB, 2.791235E-01_EB, 2.659585E-01_EB, 2.534088E-01_EB, &
    2.414824E-01_EB, 2.301745E-01_EB, 2.194706E-01_EB, 2.093505E-01_EB, 1.997898E-01_EB,  &
    3.985343E-01_EB, 4.219812E-01_EB, 4.315127E-01_EB, 4.312739E-01_EB, 4.244043E-01_EB, 4.131882E-01_EB, &
    3.992551E-01_EB, 3.837524E-01_EB, 3.674789E-01_EB, 3.509849E-01_EB, 3.346426E-01_EB, 3.186983E-01_EB, &
    3.033094E-01_EB, 2.885705E-01_EB, 2.745322E-01_EB, 2.612148E-01_EB, 2.486176E-01_EB, 2.367259E-01_EB, &
    2.255164E-01_EB, 2.149595E-01_EB, 2.050228E-01_EB, 1.956725E-01_EB, 1.868739E-01_EB/),(/23,8/))

SD3_CH3OH(1:23,17:24) = RESHAPE((/ &  ! 1300-1475 cm-1
    6.899937E-01_EB, 6.029500E-01_EB, 5.333194E-01_EB, 4.761820E-01_EB, 4.283742E-01_EB, 3.877616E-01_EB, &
    3.528449E-01_EB, 3.225319E-01_EB, 2.960036E-01_EB, 2.726290E-01_EB, 2.519127E-01_EB, 2.334579E-01_EB, &
    2.169423E-01_EB, 2.021016E-01_EB, 1.887157E-01_EB, 1.766005E-01_EB, 1.656005E-01_EB, 1.555835E-01_EB, &
    1.464364E-01_EB, 1.380620E-01_EB, 1.303764E-01_EB, 1.233068E-01_EB, 1.167893E-01_EB,  &
    8.225446E-01_EB, 6.952057E-01_EB, 5.998490E-01_EB, 5.255387E-01_EB, 4.658527E-01_EB, 4.167829E-01_EB, &
    3.756973E-01_EB, 3.407914E-01_EB, 3.107827E-01_EB, 2.847299E-01_EB, 2.619239E-01_EB, 2.418187E-01_EB, &
    2.239853E-01_EB, 2.080816E-01_EB, 1.938304E-01_EB, 1.810049E-01_EB, 1.694169E-01_EB, 1.589099E-01_EB, &
    1.493516E-01_EB, 1.406299E-01_EB, 1.326492E-01_EB, 1.253272E-01_EB, 1.185931E-01_EB,  &
    7.609519E-01_EB, 6.572483E-01_EB, 5.775783E-01_EB, 5.141108E-01_EB, 4.621339E-01_EB, 4.186506E-01_EB, &
    3.816628E-01_EB, 3.497820E-01_EB, 3.220083E-01_EB, 2.975997E-01_EB, 2.759902E-01_EB, 2.567388E-01_EB, &
    2.394950E-01_EB, 2.239759E-01_EB, 2.099498E-01_EB, 1.972249E-01_EB, 1.856405E-01_EB, 1.750612E-01_EB, &
    1.653719E-01_EB, 1.564736E-01_EB, 1.482815E-01_EB, 1.407219E-01_EB, 1.337306E-01_EB,  &
    5.794263E-01_EB, 5.286912E-01_EB, 4.873533E-01_EB, 4.525736E-01_EB, 4.225957E-01_EB, 3.962844E-01_EB, &
    3.728756E-01_EB, 3.518341E-01_EB, 3.327711E-01_EB, 3.153935E-01_EB, 2.994740E-01_EB, 2.848314E-01_EB, &
    2.713171E-01_EB, 2.588080E-01_EB, 2.471997E-01_EB, 2.364025E-01_EB, 2.263391E-01_EB, 2.169421E-01_EB, &
    2.081515E-01_EB, 1.999151E-01_EB, 1.921859E-01_EB, 1.849222E-01_EB, 1.780864E-01_EB,  &
    4.243445E-01_EB, 3.932094E-01_EB, 3.675466E-01_EB, 3.456761E-01_EB, 3.265659E-01_EB, 3.095559E-01_EB, &
    2.942085E-01_EB, 2.802223E-01_EB, 2.673816E-01_EB, 2.555266E-01_EB, 2.445346E-01_EB, 2.343074E-01_EB, &
    2.247658E-01_EB, 2.158427E-01_EB, 2.074816E-01_EB, 1.996333E-01_EB, 1.922544E-01_EB, 1.853071E-01_EB, &
    1.787572E-01_EB, 1.725743E-01_EB, 1.667308E-01_EB, 1.612018E-01_EB, 1.559650E-01_EB,  &
    2.991172E-01_EB, 2.733217E-01_EB, 2.527760E-01_EB, 2.357921E-01_EB, 2.213461E-01_EB, 2.087877E-01_EB, &
    1.976867E-01_EB, 1.877478E-01_EB, 1.787612E-01_EB, 1.705725E-01_EB, 1.630652E-01_EB, 1.561480E-01_EB, &
    1.497484E-01_EB, 1.438071E-01_EB, 1.382751E-01_EB, 1.331105E-01_EB, 1.282779E-01_EB, 1.237469E-01_EB, &
    1.194906E-01_EB, 1.154852E-01_EB, 1.117103E-01_EB, 1.081469E-01_EB, 1.047789E-01_EB,  &
    3.841022E-01_EB, 3.513805E-01_EB, 3.251838E-01_EB, 3.034426E-01_EB, 2.848923E-01_EB, 2.687247E-01_EB, &
    2.544032E-01_EB, 2.415580E-01_EB, 2.299260E-01_EB, 2.193131E-01_EB, 2.095722E-01_EB, 2.005887E-01_EB, &
    1.922707E-01_EB, 1.845432E-01_EB, 1.773437E-01_EB, 1.706194E-01_EB, 1.643252E-01_EB, 1.584217E-01_EB, &
    1.528748E-01_EB, 1.476543E-01_EB, 1.427333E-01_EB, 1.380880E-01_EB, 1.336971E-01_EB,  &
    2.736364E-01_EB, 2.518386E-01_EB, 2.341757E-01_EB, 2.193566E-01_EB, 2.065874E-01_EB, 1.953585E-01_EB, &
    1.853307E-01_EB, 1.762704E-01_EB, 1.680112E-01_EB, 1.604306E-01_EB, 1.534353E-01_EB, 1.469527E-01_EB, &
    1.409238E-01_EB, 1.353010E-01_EB, 1.300437E-01_EB, 1.251179E-01_EB, 1.204936E-01_EB, 1.161450E-01_EB, &
    1.120494E-01_EB, 1.081866E-01_EB, 1.045383E-01_EB, 1.010884E-01_EB, 9.782232E-02_EB/),(/23,8/))

SD3_CH3OH(1:23,25:29) = RESHAPE((/ &  ! 1500-1600 cm-1
    1.476547E-01_EB, 1.417188E-01_EB, 1.366044E-01_EB, 1.320475E-01_EB, 1.278834E-01_EB, 1.240084E-01_EB, &
    1.203574E-01_EB, 1.168895E-01_EB, 1.135782E-01_EB, 1.104067E-01_EB, 1.073633E-01_EB, 1.044396E-01_EB, &
    1.016297E-01_EB, 9.892818E-02_EB, 9.633071E-02_EB, 9.383318E-02_EB, 9.143167E-02_EB, 8.912240E-02_EB, &
    8.690172E-02_EB, 8.476606E-02_EB, 8.271184E-02_EB, 8.073569E-02_EB, 7.883419E-02_EB,  &
    1.006138E-01_EB, 1.006161E-01_EB, 1.004414E-01_EB, 1.000995E-01_EB, 9.959951E-02_EB, 9.895281E-02_EB, &
    9.817455E-02_EB, 9.728122E-02_EB, 9.629013E-02_EB, 9.521811E-02_EB, 9.408087E-02_EB, 9.289241E-02_EB, &
    9.166540E-02_EB, 9.041096E-02_EB, 8.913833E-02_EB, 8.785556E-02_EB, 8.656941E-02_EB, 8.528550E-02_EB, &
    8.400853E-02_EB, 8.274235E-02_EB, 8.149005E-02_EB, 8.025408E-02_EB, 7.903649E-02_EB,  &
    4.557425E-02_EB, 5.203435E-02_EB, 5.843224E-02_EB, 6.456202E-02_EB, 7.029713E-02_EB, 7.556885E-02_EB, &
    8.034824E-02_EB, 8.463279E-02_EB, 8.843707E-02_EB, 9.178554E-02_EB, 9.470880E-02_EB, 9.723982E-02_EB, &
    9.941225E-02_EB, 1.012588E-01_EB, 1.028111E-01_EB, 1.040983E-01_EB, 1.051478E-01_EB, 1.059843E-01_EB, &
    1.066308E-01_EB, 1.071078E-01_EB, 1.074338E-01_EB, 1.076256E-01_EB, 1.076983E-01_EB,  &
    2.128939E-02_EB, 3.138345E-02_EB, 4.193721E-02_EB, 5.246279E-02_EB, 6.264071E-02_EB, 7.227369E-02_EB, &
    8.124992E-02_EB, 8.951601E-02_EB, 9.705773E-02_EB, 1.038867E-01_EB, 1.100310E-01_EB, 1.155279E-01_EB, &
    1.204203E-01_EB, 1.247531E-01_EB, 1.285711E-01_EB, 1.319184E-01_EB, 1.348363E-01_EB, 1.373642E-01_EB, &
    1.395383E-01_EB, 1.413925E-01_EB, 1.429574E-01_EB, 1.442613E-01_EB, 1.453299E-01_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

!-------------------------Toluene DATA-------------------


! THERE ARE 4 BANDS FOR Toluene

! BAND #1: 1000 cm-1 - 1150 cm-1 
! BAND #2: 1300 cm-1 - 1936 cm-1 
! BAND #3: 2700 cm-1 - 3200 cm-1 
! BAND #4: 700 cm-1 - 800 cm-1 

N_TEMP_C7H8 = 23
N_BAND_C7H8 = 4

ALLOCATE(SD_C7H8_TEMP(N_TEMP_C7H8)) 

! INITIALIZE BANDS WAVENUMBER BOUNDS FOR Toluene ABSOPTION COEFFICIENTS.
! BANDS ARE ORGANIZED BY ROW: 
! 1ST COLUMN IS THE LOWER BOUND BAND LIMIT IN cm-1. 
! 2ND COLUMN IS THE UPPER BOUND BAND LIMIT IN cm-1. 
! 3RD COLUMN IS THE STRIDE BETWEEN WAVENUMBERS IN BAND.
! IF 3RD COLUMN = 0., THEN THE BAND IS CALCULATED AND NOT TABULATED.

ALLOCATE(OM_BND_C7H8(N_BAND_C7H8,3)) 

OM_BND_C7H8 = RESHAPE((/ &
    1000._EB, 1300._EB, 2700._EB,  700._EB, &
    1150._EB, 1936._EB, 3200._EB,  800._EB, &
     25._EB, 25._EB, 25._EB, 25._EB/),(/N_BAND_C7H8,3/)) 

SD_C7H8_TEMP = (/ &
    300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
    600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
    900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
    1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(SD1_C7H8(N_TEMP_C7H8,7)) 

! BAND #1: 1000 cm-1 - 1150 cm-1 

SD1_C7H8(1:23,1:7) = RESHAPE((/ &  ! 1000-1150 cm-1
    1.037831E-01_EB, 1.244139E-01_EB, 1.458370E-01_EB, 1.668692E-01_EB, 1.864671E-01_EB, 2.039026E-01_EB, &
    2.187775E-01_EB, 2.309618E-01_EB, 2.405154E-01_EB, 2.476172E-01_EB, 2.525099E-01_EB, 2.554620E-01_EB, &
    2.567435E-01_EB, 2.566101E-01_EB, 2.552965E-01_EB, 2.530112E-01_EB, 2.499376E-01_EB, 2.462339E-01_EB, &
    2.420353E-01_EB, 2.374566E-01_EB, 2.325943E-01_EB, 2.275292E-01_EB, 2.223282E-01_EB,  &
    3.790980E-01_EB, 3.739208E-01_EB, 3.641351E-01_EB, 3.517367E-01_EB, 3.380077E-01_EB, 3.237647E-01_EB, &
    3.095221E-01_EB, 2.955978E-01_EB, 2.821806E-01_EB, 2.693755E-01_EB, 2.572329E-01_EB, 2.457675E-01_EB, &
    2.349717E-01_EB, 2.248241E-01_EB, 2.152952E-01_EB, 2.063514E-01_EB, 1.979570E-01_EB, 1.900764E-01_EB, &
    1.826746E-01_EB, 1.757181E-01_EB, 1.691753E-01_EB, 1.630164E-01_EB, 1.572136E-01_EB,  &
    2.460571E-01_EB, 2.679716E-01_EB, 2.841899E-01_EB, 2.958685E-01_EB, 3.039345E-01_EB, 3.091330E-01_EB, &
    3.120621E-01_EB, 3.132010E-01_EB, 3.129332E-01_EB, 3.115653E-01_EB, 3.093420E-01_EB, 3.064589E-01_EB, &
    3.030721E-01_EB, 2.993062E-01_EB, 2.952611E-01_EB, 2.910166E-01_EB, 2.866361E-01_EB, 2.821704E-01_EB, &
    2.776601E-01_EB, 2.731372E-01_EB, 2.686273E-01_EB, 2.641504E-01_EB, 2.597220E-01_EB,  &
    3.441742E-01_EB, 3.312444E-01_EB, 3.184065E-01_EB, 3.058737E-01_EB, 2.937592E-01_EB, 2.821275E-01_EB, &
    2.710143E-01_EB, 2.604353E-01_EB, 2.503916E-01_EB, 2.408741E-01_EB, 2.318667E-01_EB, 2.233487E-01_EB, &
    2.152971E-01_EB, 2.076866E-01_EB, 2.004924E-01_EB, 1.936894E-01_EB, 1.872532E-01_EB, 1.811607E-01_EB, &
    1.753898E-01_EB, 1.699196E-01_EB, 1.647306E-01_EB, 1.598044E-01_EB, 1.551239E-01_EB,  &
    1.176276E-01_EB, 1.201819E-01_EB, 1.218991E-01_EB, 1.228917E-01_EB, 1.232736E-01_EB, 1.231523E-01_EB, &
    1.226231E-01_EB, 1.217689E-01_EB, 1.206590E-01_EB, 1.193513E-01_EB, 1.178932E-01_EB, 1.163236E-01_EB, &
    1.146737E-01_EB, 1.129690E-01_EB, 1.112298E-01_EB, 1.094723E-01_EB, 1.077094E-01_EB, 1.059515E-01_EB, &
    1.042067E-01_EB, 1.024811E-01_EB, 1.007798E-01_EB, 9.910640E-02_EB, 9.746367E-02_EB,  &
    2.183023E-02_EB, 3.401372E-02_EB, 4.805533E-02_EB, 6.312733E-02_EB, 7.850550E-02_EB, 9.362912E-02_EB, &
    1.081018E-01_EB, 1.216673E-01_EB, 1.341786E-01_EB, 1.455690E-01_EB, 1.558284E-01_EB, 1.649843E-01_EB, &
    1.730883E-01_EB, 1.802058E-01_EB, 1.864094E-01_EB, 1.917744E-01_EB, 1.963746E-01_EB, 2.002816E-01_EB, &
    2.035624E-01_EB, 2.062798E-01_EB, 2.084916E-01_EB, 2.102507E-01_EB, 2.116054E-01_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,7/))

ALLOCATE(SD2_C7H8(N_TEMP_C7H8,26)) 

! BAND #2: 1300 cm-1 - 1936 cm-1 

SD2_C7H8(1:23,1:8) = RESHAPE((/ &  ! 1300-1475 cm-1
    1.692859E-02_EB, 3.320213E-02_EB, 5.336588E-02_EB, 7.530337E-02_EB, 9.716877E-02_EB, 1.176361E-01_EB, &
    1.358958E-01_EB, 1.515515E-01_EB, 1.645009E-01_EB, 1.748302E-01_EB, 1.827367E-01_EB, 1.884726E-01_EB, &
    1.923090E-01_EB, 1.945136E-01_EB, 1.953375E-01_EB, 1.950084E-01_EB, 1.937292E-01_EB, 1.916759E-01_EB, &
    1.890012E-01_EB, 1.858352E-01_EB, 1.822878E-01_EB, 1.784518E-01_EB, 1.744048E-01_EB,  &
    1.772093E-02_EB, 3.520701E-02_EB, 5.720689E-02_EB, 8.148831E-02_EB, 1.060309E-01_EB, 1.293307E-01_EB, &
    1.504260E-01_EB, 1.688046E-01_EB, 1.842846E-01_EB, 1.969038E-01_EB, 2.068347E-01_EB, 2.143228E-01_EB, &
    2.196448E-01_EB, 2.230819E-01_EB, 2.249032E-01_EB, 2.253577E-01_EB, 2.246689E-01_EB, 2.230351E-01_EB, &
    2.206288E-01_EB, 2.175991E-01_EB, 2.140736E-01_EB, 2.101609E-01_EB, 2.059525E-01_EB,  &
    4.926891E-02_EB, 7.927959E-02_EB, 1.105207E-01_EB, 1.399760E-01_EB, 1.658076E-01_EB, 1.871667E-01_EB, &
    2.038884E-01_EB, 2.162220E-01_EB, 2.246306E-01_EB, 2.296632E-01_EB, 2.318750E-01_EB, 2.317848E-01_EB, &
    2.298553E-01_EB, 2.264859E-01_EB, 2.220143E-01_EB, 2.167206E-01_EB, 2.108348E-01_EB, 2.045431E-01_EB, &
    1.979950E-01_EB, 1.913094E-01_EB, 1.845798E-01_EB, 1.778792E-01_EB, 1.712634E-01_EB,  &
    2.135190E-01_EB, 2.295950E-01_EB, 2.390418E-01_EB, 2.434827E-01_EB, 2.442819E-01_EB, 2.425036E-01_EB, &
    2.389476E-01_EB, 2.342027E-01_EB, 2.286969E-01_EB, 2.227380E-01_EB, 2.165458E-01_EB, 2.102759E-01_EB, &
    2.040363E-01_EB, 1.979020E-01_EB, 1.919224E-01_EB, 1.861298E-01_EB, 1.805433E-01_EB, 1.751731E-01_EB, &
    1.700230E-01_EB, 1.650919E-01_EB, 1.603760E-01_EB, 1.558691E-01_EB, 1.515637E-01_EB,  &
    1.301418E-01_EB, 1.470068E-01_EB, 1.610184E-01_EB, 1.725853E-01_EB, 1.820615E-01_EB, 1.897508E-01_EB, &
    1.959140E-01_EB, 2.007755E-01_EB, 2.045281E-01_EB, 2.073376E-01_EB, 2.093459E-01_EB, 2.106747E-01_EB, &
    2.114282E-01_EB, 2.116951E-01_EB, 2.115513E-01_EB, 2.110617E-01_EB, 2.102815E-01_EB, 2.092579E-01_EB, &
    2.080309E-01_EB, 2.066353E-01_EB, 2.051003E-01_EB, 2.034510E-01_EB, 2.017090E-01_EB,  &
    3.146415E-01_EB, 3.203372E-01_EB, 3.241099E-01_EB, 3.262934E-01_EB, 3.271492E-01_EB, 3.268962E-01_EB, &
    3.257218E-01_EB, 3.237890E-01_EB, 3.212381E-01_EB, 3.181896E-01_EB, 3.147461E-01_EB, 3.109945E-01_EB, &
    3.070073E-01_EB, 3.028458E-01_EB, 2.985599E-01_EB, 2.941921E-01_EB, 2.897764E-01_EB, 2.853411E-01_EB, &
    2.809094E-01_EB, 2.764999E-01_EB, 2.721276E-01_EB, 2.678045E-01_EB, 2.635403E-01_EB,  &
    5.746809E-01_EB, 5.376320E-01_EB, 5.068554E-01_EB, 4.804227E-01_EB, 4.571392E-01_EB, 4.362410E-01_EB, &
    4.172236E-01_EB, 3.997438E-01_EB, 3.835595E-01_EB, 3.684936E-01_EB, 3.544124E-01_EB, 3.412107E-01_EB, &
    3.288036E-01_EB, 3.171200E-01_EB, 3.060995E-01_EB, 2.956895E-01_EB, 2.858438E-01_EB, 2.765209E-01_EB, &
    2.676840E-01_EB, 2.592989E-01_EB, 2.513353E-01_EB, 2.437651E-01_EB, 2.365626E-01_EB,  &
    1.216680E+00_EB, 1.127921E+00_EB, 1.058738E+00_EB, 1.002359E+00_EB, 9.547658E-01_EB, 9.134679E-01_EB, &
    8.768605E-01_EB, 8.438767E-01_EB, 8.137853E-01_EB, 7.860700E-01_EB, 7.603557E-01_EB, 7.363618E-01_EB, &
    7.138728E-01_EB, 6.927186E-01_EB, 6.727620E-01_EB, 6.538894E-01_EB, 6.360054E-01_EB, 6.190284E-01_EB, &
    6.028873E-01_EB, 5.875201E-01_EB, 5.728714E-01_EB, 5.588918E-01_EB, 5.455365E-01_EB/),(/23,8/))

SD2_C7H8(1:23,9:16) = RESHAPE((/ &  ! 1500-1675 cm-1
    1.023400E+00_EB, 9.124216E-01_EB, 8.265950E-01_EB, 7.576231E-01_EB, 7.005116E-01_EB, 6.520972E-01_EB, &
    6.102884E-01_EB, 5.736500E-01_EB, 5.411635E-01_EB, 5.120838E-01_EB, 4.858505E-01_EB, 4.620323E-01_EB, &
    4.402892E-01_EB, 4.203481E-01_EB, 4.019868E-01_EB, 3.850202E-01_EB, 3.692941E-01_EB, 3.546767E-01_EB, &
    3.410560E-01_EB, 3.283348E-01_EB, 3.164286E-01_EB, 3.052635E-01_EB, 2.947744E-01_EB,  &
    2.357062E-01_EB, 2.424195E-01_EB, 2.478945E-01_EB, 2.522236E-01_EB, 2.555069E-01_EB, 2.578497E-01_EB, &
    2.593578E-01_EB, 2.601342E-01_EB, 2.602759E-01_EB, 2.598725E-01_EB, 2.590048E-01_EB, 2.577446E-01_EB, &
    2.561556E-01_EB, 2.542931E-01_EB, 2.522053E-01_EB, 2.499340E-01_EB, 2.475146E-01_EB, 2.449779E-01_EB, &
    2.423502E-01_EB, 2.396536E-01_EB, 2.369075E-01_EB, 2.341276E-01_EB, 2.313276E-01_EB,  &
    1.400331E-01_EB, 1.594231E-01_EB, 1.762783E-01_EB, 1.907722E-01_EB, 2.031260E-01_EB, 2.135669E-01_EB, &
    2.223125E-01_EB, 2.295643E-01_EB, 2.355059E-01_EB, 2.403024E-01_EB, 2.441009E-01_EB, 2.470320E-01_EB, &
    2.492106E-01_EB, 2.507375E-01_EB, 2.517010E-01_EB, 2.521781E-01_EB, 2.522363E-01_EB, 2.519341E-01_EB, &
    2.513222E-01_EB, 2.504454E-01_EB, 2.493420E-01_EB, 2.480458E-01_EB, 2.465858E-01_EB,  &
    2.613801E-01_EB, 2.971173E-01_EB, 3.273311E-01_EB, 3.527523E-01_EB, 3.740190E-01_EB, 3.916818E-01_EB, &
    4.062191E-01_EB, 4.180470E-01_EB, 4.275293E-01_EB, 4.349837E-01_EB, 4.406877E-01_EB, 4.448830E-01_EB, &
    4.477803E-01_EB, 4.495623E-01_EB, 4.503878E-01_EB, 4.503940E-01_EB, 4.497003E-01_EB, 4.484097E-01_EB, &
    4.466114E-01_EB, 4.443827E-01_EB, 4.417905E-01_EB, 4.388924E-01_EB, 4.357386E-01_EB,  &
    7.802304E-01_EB, 6.810847E-01_EB, 5.994172E-01_EB, 5.316692E-01_EB, 4.748939E-01_EB, 4.268126E-01_EB, &
    3.856938E-01_EB, 3.502182E-01_EB, 3.193725E-01_EB, 2.923675E-01_EB, 2.685805E-01_EB, 2.475144E-01_EB, &
    2.287668E-01_EB, 2.120089E-01_EB, 1.969698E-01_EB, 1.834231E-01_EB, 1.711797E-01_EB, 1.600793E-01_EB, &
    1.499859E-01_EB, 1.407829E-01_EB, 1.323705E-01_EB, 1.246619E-01_EB, 1.175823E-01_EB,  &
    9.933114E-02_EB, 1.010553E-01_EB, 1.027618E-01_EB, 1.044763E-01_EB, 1.061578E-01_EB, 1.077546E-01_EB, &
    1.092242E-01_EB, 1.105370E-01_EB, 1.116764E-01_EB, 1.126352E-01_EB, 1.134139E-01_EB, 1.140180E-01_EB, &
    1.144567E-01_EB, 1.147407E-01_EB, 1.148822E-01_EB, 1.148933E-01_EB, 1.147862E-01_EB, 1.145730E-01_EB, &
    1.142646E-01_EB, 1.138716E-01_EB, 1.134035E-01_EB, 1.128695E-01_EB, 1.122775E-01_EB,  &
    4.259418E-02_EB, 4.612365E-02_EB, 4.901499E-02_EB, 5.138251E-02_EB, 5.331256E-02_EB, 5.487248E-02_EB, &
    5.611660E-02_EB, 5.709005E-02_EB, 5.783099E-02_EB, 5.837207E-02_EB, 5.874148E-02_EB, 5.896357E-02_EB, &
    5.905943E-02_EB, 5.904734E-02_EB, 5.894315E-02_EB, 5.876054E-02_EB, 5.851139E-02_EB, 5.820596E-02_EB, &
    5.785316E-02_EB, 5.746063E-02_EB, 5.703504E-02_EB, 5.658212E-02_EB, 5.610684E-02_EB,  &
    4.557166E-02_EB, 4.915785E-02_EB, 5.209602E-02_EB, 5.450251E-02_EB, 5.646431E-02_EB, 5.804888E-02_EB, &
    5.931063E-02_EB, 6.029487E-02_EB, 6.104008E-02_EB, 6.157940E-02_EB, 6.194157E-02_EB, 6.215153E-02_EB, &
    6.223095E-02_EB, 6.219870E-02_EB, 6.207116E-02_EB, 6.186251E-02_EB, 6.158512E-02_EB, 6.124959E-02_EB, &
    6.086520E-02_EB, 6.043996E-02_EB, 5.998075E-02_EB, 5.949356E-02_EB, 5.898355E-02_EB/),(/23,8/))

SD2_C7H8(1:23,17:24) = RESHAPE((/ &  ! 1700-1875 cm-1
    7.586269E-03_EB, 1.711922E-02_EB, 3.043808E-02_EB, 4.632439E-02_EB, 6.337353E-02_EB, 8.035813E-02_EB, &
    9.636732E-02_EB, 1.108125E-01_EB, 1.233749E-01_EB, 1.339379E-01_EB, 1.425233E-01_EB, 1.492413E-01_EB, &
    1.542524E-01_EB, 1.577407E-01_EB, 1.598960E-01_EB, 1.609027E-01_EB, 1.609336E-01_EB, 1.601455E-01_EB, &
    1.586788E-01_EB, 1.566571E-01_EB, 1.541877E-01_EB, 1.513633E-01_EB, 1.482632E-01_EB,  &
    8.125833E-03_EB, 1.831305E-02_EB, 3.253283E-02_EB, 4.948498E-02_EB, 6.767454E-02_EB, 8.579670E-02_EB, &
    1.028833E-01_EB, 1.183089E-01_EB, 1.317340E-01_EB, 1.430333E-01_EB, 1.522294E-01_EB, 1.594378E-01_EB, &
    1.648285E-01_EB, 1.685959E-01_EB, 1.709415E-01_EB, 1.720607E-01_EB, 1.721370E-01_EB, 1.713372E-01_EB, &
    1.698104E-01_EB, 1.676883E-01_EB, 1.650853E-01_EB, 1.621000E-01_EB, 1.588174E-01_EB,  &
    1.979055E-02_EB, 3.419787E-02_EB, 5.154665E-02_EB, 7.044249E-02_EB, 8.954887E-02_EB, 1.078047E-01_EB, &
    1.244818E-01_EB, 1.391540E-01_EB, 1.516336E-01_EB, 1.619033E-01_EB, 1.700581E-01_EB, 1.762597E-01_EB, &
    1.807046E-01_EB, 1.836015E-01_EB, 1.851565E-01_EB, 1.855650E-01_EB, 1.850059E-01_EB, 1.836407E-01_EB, &
    1.816118E-01_EB, 1.790437E-01_EB, 1.760443E-01_EB, 1.727063E-01_EB, 1.691084E-01_EB,  &
    3.382598E-02_EB, 5.035882E-02_EB, 6.680258E-02_EB, 8.219769E-02_EB, 9.605987E-02_EB, 1.082064E-01_EB, &
    1.186306E-01_EB, 1.274215E-01_EB, 1.347146E-01_EB, 1.406641E-01_EB, 1.454254E-01_EB, 1.491471E-01_EB, &
    1.519661E-01_EB, 1.540061E-01_EB, 1.553770E-01_EB, 1.561758E-01_EB, 1.564872E-01_EB, 1.563849E-01_EB, &
    1.559328E-01_EB, 1.551863E-01_EB, 1.541929E-01_EB, 1.529938E-01_EB, 1.516241E-01_EB,  &
    2.397733E-02_EB, 4.086562E-02_EB, 6.055781E-02_EB, 8.177685E-02_EB, 1.034905E-01_EB, 1.249345E-01_EB, &
    1.455781E-01_EB, 1.650752E-01_EB, 1.832190E-01_EB, 1.999044E-01_EB, 2.150983E-01_EB, 2.288171E-01_EB, &
    2.411103E-01_EB, 2.520484E-01_EB, 2.617140E-01_EB, 2.701960E-01_EB, 2.775847E-01_EB, 2.839694E-01_EB, &
    2.894365E-01_EB, 2.940677E-01_EB, 2.979400E-01_EB, 3.011248E-01_EB, 3.036887E-01_EB,  &
    9.976805E-02_EB, 1.250605E-01_EB, 1.479614E-01_EB, 1.678470E-01_EB, 1.846035E-01_EB, 1.983821E-01_EB, &
    2.094543E-01_EB, 2.181335E-01_EB, 2.247329E-01_EB, 2.295463E-01_EB, 2.328394E-01_EB, 2.348475E-01_EB, &
    2.357755E-01_EB, 2.358008E-01_EB, 2.350757E-01_EB, 2.337302E-01_EB, 2.318752E-01_EB, 2.296048E-01_EB, &
    2.269988E-01_EB, 2.241247E-01_EB, 2.210394E-01_EB, 2.177911E-01_EB, 2.144203E-01_EB,  &
    1.141042E-01_EB, 1.239540E-01_EB, 1.330543E-01_EB, 1.414094E-01_EB, 1.490126E-01_EB, 1.558669E-01_EB, &
    1.619887E-01_EB, 1.674067E-01_EB, 1.721581E-01_EB, 1.762856E-01_EB, 1.798351E-01_EB, 1.828529E-01_EB, &
    1.853850E-01_EB, 1.874758E-01_EB, 1.891672E-01_EB, 1.904983E-01_EB, 1.915061E-01_EB, 1.922236E-01_EB, &
    1.926817E-01_EB, 1.929084E-01_EB, 1.929288E-01_EB, 1.927660E-01_EB, 1.924406E-01_EB,  &
    6.298505E-02_EB, 8.510611E-02_EB, 1.072605E-01_EB, 1.287323E-01_EB, 1.490945E-01_EB, 1.681118E-01_EB, &
    1.856727E-01_EB, 2.017464E-01_EB, 2.163536E-01_EB, 2.295473E-01_EB, 2.413998E-01_EB, 2.519942E-01_EB, &
    2.614187E-01_EB, 2.697626E-01_EB, 2.771126E-01_EB, 2.835528E-01_EB, 2.891625E-01_EB, 2.940159E-01_EB, &
    2.981823E-01_EB, 3.017252E-01_EB, 3.047034E-01_EB, 3.071704E-01_EB, 3.091754E-01_EB/),(/23,8/))

SD2_C7H8(1:23,25:26) = RESHAPE((/ &  ! 1900-1925 cm-1
    1.742141E-01_EB, 1.775203E-01_EB, 1.796168E-01_EB, 1.808072E-01_EB, 1.812946E-01_EB, 1.812215E-01_EB, &
    1.806919E-01_EB, 1.797874E-01_EB, 1.785732E-01_EB, 1.771035E-01_EB, 1.754243E-01_EB, 1.735747E-01_EB, &
    1.715881E-01_EB, 1.694934E-01_EB, 1.673150E-01_EB, 1.650741E-01_EB, 1.627886E-01_EB, 1.604738E-01_EB, &
    1.581424E-01_EB, 1.558055E-01_EB, 1.534721E-01_EB, 1.511499E-01_EB, 1.488452E-01_EB,  &
    4.428388E-02_EB, 4.541289E-02_EB, 4.628056E-02_EB, 4.694935E-02_EB, 4.745757E-02_EB, 4.783084E-02_EB, &
    4.808774E-02_EB, 4.824295E-02_EB, 4.830895E-02_EB, 4.829650E-02_EB, 4.821534E-02_EB, 4.807411E-02_EB, &
    4.788071E-02_EB, 4.764225E-02_EB, 4.736497E-02_EB, 4.705461E-02_EB, 4.671619E-02_EB, 4.635418E-02_EB, &
    4.597257E-02_EB, 4.557478E-02_EB, 4.516388E-02_EB, 4.474258E-02_EB, 4.431326E-02_EB/),(/23,2/))

ALLOCATE(SD3_C7H8(N_TEMP_C7H8,21)) 

! BAND #3: 2700 cm-1 - 3200 cm-1 

SD3_C7H8(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
    1.478412E-02_EB, 2.760122E-02_EB, 4.292571E-02_EB, 5.927613E-02_EB, 7.546969E-02_EB, 9.070428E-02_EB, &
    1.045076E-01_EB, 1.166519E-01_EB, 1.270745E-01_EB, 1.358168E-01_EB, 1.429811E-01_EB, 1.487010E-01_EB, &
    1.531230E-01_EB, 1.563949E-01_EB, 1.586582E-01_EB, 1.600454E-01_EB, 1.606771E-01_EB, 1.606618E-01_EB, &
    1.600966E-01_EB, 1.590668E-01_EB, 1.576477E-01_EB, 1.559048E-01_EB, 1.538952E-01_EB,  &
    3.908357E-02_EB, 5.729946E-02_EB, 7.597442E-02_EB, 9.384344E-02_EB, 1.101317E-01_EB, 1.244548E-01_EB, &
    1.366941E-01_EB, 1.468937E-01_EB, 1.551883E-01_EB, 1.617556E-01_EB, 1.667886E-01_EB, 1.704780E-01_EB, &
    1.730035E-01_EB, 1.745297E-01_EB, 1.752040E-01_EB, 1.751573E-01_EB, 1.745040E-01_EB, 1.733439E-01_EB, &
    1.717636E-01_EB, 1.698375E-01_EB, 1.676301E-01_EB, 1.651962E-01_EB, 1.625829E-01_EB,  &
    4.112356E-02_EB, 5.583260E-02_EB, 7.007568E-02_EB, 8.349028E-02_EB, 9.592829E-02_EB, 1.073558E-01_EB, &
    1.177970E-01_EB, 1.273027E-01_EB, 1.359345E-01_EB, 1.437565E-01_EB, 1.508317E-01_EB, 1.572192E-01_EB, &
    1.629746E-01_EB, 1.681493E-01_EB, 1.727904E-01_EB, 1.769415E-01_EB, 1.806431E-01_EB, 1.839318E-01_EB, &
    1.868416E-01_EB, 1.894040E-01_EB, 1.916475E-01_EB, 1.935986E-01_EB, 1.952815E-01_EB,  &
    3.022201E-02_EB, 4.424940E-02_EB, 5.893774E-02_EB, 7.368576E-02_EB, 8.811545E-02_EB, 1.020008E-01_EB, &
    1.152142E-01_EB, 1.276908E-01_EB, 1.394046E-01_EB, 1.503547E-01_EB, 1.605548E-01_EB, 1.700286E-01_EB, &
    1.788050E-01_EB, 1.869164E-01_EB, 1.943969E-01_EB, 2.012811E-01_EB, 2.076031E-01_EB, 2.133970E-01_EB, &
    2.186952E-01_EB, 2.235294E-01_EB, 2.279297E-01_EB, 2.319247E-01_EB, 2.355417E-01_EB,  &
    2.988847E-02_EB, 4.808938E-02_EB, 6.861285E-02_EB, 9.036903E-02_EB, 1.125471E-01_EB, 1.345782E-01_EB, &
    1.560786E-01_EB, 1.767978E-01_EB, 1.965795E-01_EB, 2.153330E-01_EB, 2.330130E-01_EB, 2.496052E-01_EB, &
    2.651180E-01_EB, 2.795737E-01_EB, 2.930054E-01_EB, 3.054525E-01_EB, 3.169587E-01_EB, 3.275702E-01_EB, &
    3.373338E-01_EB, 3.462967E-01_EB, 3.545053E-01_EB, 3.620049E-01_EB, 3.688390E-01_EB,  &
    7.300609E-02_EB, 1.058028E-01_EB, 1.382863E-01_EB, 1.689483E-01_EB, 1.970619E-01_EB, 2.223692E-01_EB, &
    2.448717E-01_EB, 2.647045E-01_EB, 2.820630E-01_EB, 2.971638E-01_EB, 3.102233E-01_EB, 3.214476E-01_EB, &
    3.310272E-01_EB, 3.391362E-01_EB, 3.459324E-01_EB, 3.515569E-01_EB, 3.561372E-01_EB, 3.597872E-01_EB, &
    3.626089E-01_EB, 3.646935E-01_EB, 3.661225E-01_EB, 3.669693E-01_EB, 3.672991E-01_EB,  &
    5.763692E-01_EB, 5.688459E-01_EB, 5.566362E-01_EB, 5.420274E-01_EB, 5.263228E-01_EB, 5.102726E-01_EB, &
    4.943062E-01_EB, 4.786651E-01_EB, 4.634798E-01_EB, 4.488151E-01_EB, 4.346967E-01_EB, 4.211283E-01_EB, &
    4.081010E-01_EB, 3.955986E-01_EB, 3.836017E-01_EB, 3.720893E-01_EB, 3.610404E-01_EB, 3.504335E-01_EB, &
    3.402491E-01_EB, 3.304672E-01_EB, 3.210695E-01_EB, 3.120383E-01_EB, 3.033567E-01_EB,  &
    1.206676E+00_EB, 1.069048E+00_EB, 9.625163E-01_EB, 8.773986E-01_EB, 8.076640E-01_EB, 7.493510E-01_EB, &
    6.997454E-01_EB, 6.569261E-01_EB, 6.194952E-01_EB, 5.864143E-01_EB, 5.568961E-01_EB, 5.303352E-01_EB, &
    5.062588E-01_EB, 4.842928E-01_EB, 4.641379E-01_EB, 4.455519E-01_EB, 4.283366E-01_EB, 4.123285E-01_EB, &
    3.973913E-01_EB, 3.834101E-01_EB, 3.702875E-01_EB, 3.579407E-01_EB, 3.462977E-01_EB/),(/23,8/))

SD3_C7H8(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
    1.100732E+00_EB, 1.064927E+00_EB, 1.035627E+00_EB, 1.010989E+00_EB, 9.898055E-01_EB, 9.712378E-01_EB, &
    9.546794E-01_EB, 9.396782E-01_EB, 9.258909E-01_EB, 9.130543E-01_EB, 9.009653E-01_EB, 8.894672E-01_EB, &
    8.784391E-01_EB, 8.677880E-01_EB, 8.574432E-01_EB, 8.473505E-01_EB, 8.374691E-01_EB, 8.277692E-01_EB, &
    8.182278E-01_EB, 8.088290E-01_EB, 7.995611E-01_EB, 7.904159E-01_EB, 7.813883E-01_EB,  &
    2.231706E+00_EB, 1.912652E+00_EB, 1.673768E+00_EB, 1.488171E+00_EB, 1.339758E+00_EB, 1.218290E+00_EB, &
    1.116949E+00_EB, 1.031021E+00_EB, 9.571499E-01_EB, 8.928820E-01_EB, 8.363870E-01_EB, 7.862736E-01_EB, &
    7.414671E-01_EB, 7.011249E-01_EB, 6.645779E-01_EB, 6.312884E-01_EB, 6.008189E-01_EB, 5.728103E-01_EB, &
    5.469647E-01_EB, 5.230325E-01_EB, 5.008035E-01_EB, 4.800983E-01_EB, 4.607637E-01_EB,  &
    1.078196E+00_EB, 9.387864E-01_EB, 8.327788E-01_EB, 7.493044E-01_EB, 6.817561E-01_EB, 6.258786E-01_EB, &
    5.788046E-01_EB, 5.385315E-01_EB, 5.036194E-01_EB, 4.730070E-01_EB, 4.458971E-01_EB, 4.216795E-01_EB, &
    3.998804E-01_EB, 3.801262E-01_EB, 3.621191E-01_EB, 3.456188E-01_EB, 3.304293E-01_EB, 3.163890E-01_EB, &
    3.033631E-01_EB, 2.912393E-01_EB, 2.799218E-01_EB, 2.693291E-01_EB, 2.593912E-01_EB,  &
    8.437348E-01_EB, 7.629131E-01_EB, 6.992342E-01_EB, 6.475115E-01_EB, 6.044849E-01_EB, 5.679907E-01_EB, &
    5.365302E-01_EB, 5.090309E-01_EB, 4.847048E-01_EB, 4.629598E-01_EB, 4.433439E-01_EB, 4.255066E-01_EB, &
    4.091724E-01_EB, 3.941221E-01_EB, 3.801801E-01_EB, 3.672033E-01_EB, 3.550750E-01_EB, 3.436979E-01_EB, &
    3.329918E-01_EB, 3.228881E-01_EB, 3.133295E-01_EB, 3.042662E-01_EB, 2.956559E-01_EB,  &
    1.655312E+00_EB, 1.597040E+00_EB, 1.545376E+00_EB, 1.499226E+00_EB, 1.457664E+00_EB, 1.419919E+00_EB, &
    1.385355E+00_EB, 1.353446E+00_EB, 1.323760E+00_EB, 1.295943E+00_EB, 1.269706E+00_EB, 1.244817E+00_EB, &
    1.221084E+00_EB, 1.198354E+00_EB, 1.176504E+00_EB, 1.155434E+00_EB, 1.135062E+00_EB, 1.115324E+00_EB, &
    1.096165E+00_EB, 1.077543E+00_EB, 1.059422E+00_EB, 1.041771E+00_EB, 1.024566E+00_EB,  &
    4.079883E+00_EB, 3.472018E+00_EB, 3.019929E+00_EB, 2.670772E+00_EB, 2.393071E+00_EB, 2.166922E+00_EB, &
    1.979128E+00_EB, 1.820604E+00_EB, 1.684899E+00_EB, 1.567313E+00_EB, 1.464348E+00_EB, 1.373352E+00_EB, &
    1.292279E+00_EB, 1.219527E+00_EB, 1.153828E+00_EB, 1.094165E+00_EB, 1.039710E+00_EB, 9.897859E-01_EB, &
    9.438328E-01_EB, 9.013822E-01_EB, 8.620396E-01_EB, 8.254703E-01_EB, 7.913883E-01_EB,  &
    2.374205E+00_EB, 2.197359E+00_EB, 2.053008E+00_EB, 1.932367E+00_EB, 1.829600E+00_EB, 1.740653E+00_EB, &
    1.662612E+00_EB, 1.593322E+00_EB, 1.531154E+00_EB, 1.474857E+00_EB, 1.423456E+00_EB, 1.376184E+00_EB, &
    1.332430E+00_EB, 1.291702E+00_EB, 1.253605E+00_EB, 1.217812E+00_EB, 1.184059E+00_EB, 1.152122E+00_EB, &
    1.121817E+00_EB, 1.092988E+00_EB, 1.065503E+00_EB, 1.039248E+00_EB, 1.014126E+00_EB,  &
    1.996169E+00_EB, 1.800315E+00_EB, 1.646394E+00_EB, 1.521669E+00_EB, 1.418149E+00_EB, 1.330541E+00_EB, &
    1.255185E+00_EB, 1.189465E+00_EB, 1.131460E+00_EB, 1.079726E+00_EB, 1.033161E+00_EB, 9.909102E-01_EB, &
    9.523003E-01_EB, 9.167962E-01_EB, 8.839679E-01_EB, 8.534662E-01_EB, 8.250043E-01_EB, 7.983457E-01_EB, &
    7.732927E-01_EB, 7.496792E-01_EB, 7.273640E-01_EB, 7.062268E-01_EB, 6.861638E-01_EB/),(/23,8/))

SD3_C7H8(1:23,17:21) = RESHAPE((/ &  ! 3100-3200 cm-1
    7.643893E-01_EB, 6.852582E-01_EB, 6.236262E-01_EB, 5.740686E-01_EB, 5.332134E-01_EB, 4.988450E-01_EB, &
    4.694430E-01_EB, 4.439271E-01_EB, 4.215087E-01_EB, 4.015982E-01_EB, 3.837473E-01_EB, 3.676090E-01_EB, &
    3.529113E-01_EB, 3.394389E-01_EB, 3.270185E-01_EB, 3.155098E-01_EB, 3.047980E-01_EB, 2.947886E-01_EB, &
    2.854026E-01_EB, 2.765737E-01_EB, 2.682460E-01_EB, 2.603715E-01_EB, 2.529089E-01_EB,  &
    1.305546E-01_EB, 1.285607E-01_EB, 1.268699E-01_EB, 1.254048E-01_EB, 1.241115E-01_EB, 1.229502E-01_EB, &
    1.218902E-01_EB, 1.209072E-01_EB, 1.199815E-01_EB, 1.190975E-01_EB, 1.182425E-01_EB, 1.174064E-01_EB, &
    1.165817E-01_EB, 1.157622E-01_EB, 1.149438E-01_EB, 1.141230E-01_EB, 1.132979E-01_EB, 1.124671E-01_EB, &
    1.116298E-01_EB, 1.107860E-01_EB, 1.099356E-01_EB, 1.090793E-01_EB, 1.082174E-01_EB,  &
    6.396032E-02_EB, 6.396678E-02_EB, 6.396719E-02_EB, 6.396303E-02_EB, 6.395420E-02_EB, 6.393951E-02_EB, &
    6.391705E-02_EB, 6.388461E-02_EB, 6.383988E-02_EB, 6.378072E-02_EB, 6.370534E-02_EB, 6.361232E-02_EB, &
    6.350063E-02_EB, 6.336969E-02_EB, 6.321926E-02_EB, 6.304944E-02_EB, 6.286059E-02_EB, 6.265328E-02_EB, &
    6.242828E-02_EB, 6.218649E-02_EB, 6.192887E-02_EB, 6.165644E-02_EB, 6.137029E-02_EB,  &
    4.958331E-02_EB, 5.240776E-02_EB, 5.466106E-02_EB, 5.649366E-02_EB, 5.800772E-02_EB, 5.927408E-02_EB, &
    6.034279E-02_EB, 6.124987E-02_EB, 6.202164E-02_EB, 6.267782E-02_EB, 6.323351E-02_EB, 6.370052E-02_EB, &
    6.408837E-02_EB, 6.440490E-02_EB, 6.465680E-02_EB, 6.484982E-02_EB, 6.498902E-02_EB, 6.507900E-02_EB, &
    6.512384E-02_EB, 6.512730E-02_EB, 6.509286E-02_EB, 6.502369E-02_EB, 6.492275E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

ALLOCATE(SD4_C7H8(N_TEMP_C7H8,5)) 

! BAND #4: 700 cm-1 - 800 cm-1 

SD4_C7H8(1:23,1:5) = RESHAPE((/ &  ! 700-800 cm-1
    3.436335E+00_EB, 3.206637E+00_EB, 2.996047E+00_EB, 2.803348E+00_EB, 2.627385E+00_EB, 2.466909E+00_EB, &
    2.320605E+00_EB, 2.187162E+00_EB, 2.065321E+00_EB, 1.953911E+00_EB, 1.851860E+00_EB, 1.758204E+00_EB, &
    1.672078E+00_EB, 1.592712E+00_EB, 1.519426E+00_EB, 1.451612E+00_EB, 1.388738E+00_EB, 1.330328E+00_EB, &
    1.275961E+00_EB, 1.225266E+00_EB, 1.177909E+00_EB, 1.133595E+00_EB, 1.092060E+00_EB,  &
    4.813064E+00_EB, 4.147262E+00_EB, 3.629948E+00_EB, 3.215234E+00_EB, 2.875187E+00_EB, 2.591536E+00_EB, &
    2.351670E+00_EB, 2.146532E+00_EB, 1.969409E+00_EB, 1.815204E+00_EB, 1.679974E+00_EB, 1.560611E+00_EB, &
    1.454640E+00_EB, 1.360060E+00_EB, 1.275239E+00_EB, 1.198834E+00_EB, 1.129731E+00_EB, 1.066995E+00_EB, &
    1.009842E+00_EB, 9.576050E-01_EB, 9.097170E-01_EB, 8.656911E-01_EB, 8.251083E-01_EB,  &
    2.178560E-01_EB, 2.675286E-01_EB, 3.130571E-01_EB, 3.536868E-01_EB, 3.892184E-01_EB, 4.197957E-01_EB, &
    4.457530E-01_EB, 4.675155E-01_EB, 4.855363E-01_EB, 5.002599E-01_EB, 5.121029E-01_EB, 5.214444E-01_EB, &
    5.286238E-01_EB, 5.339402E-01_EB, 5.376557E-01_EB, 5.399985E-01_EB, 5.411663E-01_EB, 5.413306E-01_EB, &
    5.406389E-01_EB, 5.392189E-01_EB, 5.371808E-01_EB, 5.346193E-01_EB, 5.316164E-01_EB,  &
    1.076836E-02_EB, 1.752852E-02_EB, 2.573567E-02_EB, 3.519001E-02_EB, 4.567683E-02_EB, 5.697273E-02_EB, &
    6.885432E-02_EB, 8.110841E-02_EB, 9.353976E-02_EB, 1.059771E-01_EB, 1.182750E-01_EB, 1.303146E-01_EB, &
    1.420018E-01_EB, 1.532648E-01_EB, 1.640516E-01_EB, 1.743266E-01_EB, 1.840679E-01_EB, 1.932649E-01_EB, &
    2.019162E-01_EB, 2.100270E-01_EB, 2.176081E-01_EB, 2.246742E-01_EB, 2.312430E-01_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

!-------------------------Propylene DATA-------------------


! THERE ARE 3 BANDS FOR Propylene

! BAND #1: 1250 cm-1 - 1950 cm-1 
! BAND #2: 2700 cm-1 - 3200 cm-1 
! BAND #3: 700 cm-1 - 1150 cm-1 

N_TEMP_C3H6 = 23
N_BAND_C3H6 = 3

ALLOCATE(SD_C3H6_TEMP(N_TEMP_C3H6)) 

! INITIALIZE BANDS WAVENUMBER BOUNDS FOR Propylene ABSOPTION COEFFICIENTS.
! BANDS ARE ORGANIZED BY ROW: 
! 1ST COLUMN IS THE LOWER BOUND BAND LIMIT IN cm-1. 
! 2ND COLUMN IS THE UPPER BOUND BAND LIMIT IN cm-1. 
! 3RD COLUMN IS THE STRIDE BETWEEN WAVENUMBERS IN BAND.
! IF 3RD COLUMN = 0., THEN THE BAND IS CALCULATED AND NOT TABULATED.

ALLOCATE(OM_BND_C3H6(N_BAND_C3H6,3)) 

OM_BND_C3H6 = RESHAPE((/ &
    1250._EB, 2700._EB,  700._EB, &
    1950._EB, 3200._EB, 1150._EB, &
     25._EB, 25._EB, 25._EB/),(/N_BAND_C3H6,3/)) 

SD_C3H6_TEMP = (/ &
    300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
    600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
    900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
    1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(SD1_C3H6(N_TEMP_C3H6,29)) 

! BAND #1: 1250 cm-1 - 1950 cm-1 

SD1_C3H6(1:23,1:8) = RESHAPE((/ &  ! 1250-1425 cm-1
    9.980634E-04_EB, 1.754119E-03_EB, 2.632367E-03_EB, 3.563637E-03_EB, 4.491205E-03_EB, 5.382035E-03_EB, &
    6.212940E-03_EB, 6.974658E-03_EB, 7.660558E-03_EB, 8.269868E-03_EB, 8.806189E-03_EB, 9.274394E-03_EB, &
    9.677995E-03_EB, 1.002468E-02_EB, 1.031826E-02_EB, 1.056498E-02_EB, 1.076900E-02_EB, 1.093483E-02_EB, &
    1.106676E-02_EB, 1.116983E-02_EB, 1.124640E-02_EB, 1.129853E-02_EB, 1.133066E-02_EB,  &
    7.161402E-03_EB, 1.001729E-02_EB, 1.258444E-02_EB, 1.471671E-02_EB, 1.637746E-02_EB, 1.759582E-02_EB, &
    1.842695E-02_EB, 1.893531E-02_EB, 1.918172E-02_EB, 1.922054E-02_EB, 1.910075E-02_EB, 1.885910E-02_EB, &
    1.852802E-02_EB, 1.813235E-02_EB, 1.769187E-02_EB, 1.722158E-02_EB, 1.673396E-02_EB, 1.623805E-02_EB, &
    1.573999E-02_EB, 1.524681E-02_EB, 1.476142E-02_EB, 1.428665E-02_EB, 1.382503E-02_EB,  &
    1.584584E-02_EB, 1.993794E-02_EB, 2.338498E-02_EB, 2.618794E-02_EB, 2.839993E-02_EB, 3.009841E-02_EB, &
    3.135858E-02_EB, 3.225147E-02_EB, 3.284167E-02_EB, 3.318345E-02_EB, 3.332270E-02_EB, 3.330033E-02_EB, &
    3.314615E-02_EB, 3.288735E-02_EB, 3.254700E-02_EB, 3.214140E-02_EB, 3.168688E-02_EB, 3.119487E-02_EB, &
    3.067402E-02_EB, 3.013564E-02_EB, 2.958450E-02_EB, 2.902450E-02_EB, 2.846168E-02_EB,  &
    2.586388E-02_EB, 3.348021E-02_EB, 4.025472E-02_EB, 4.608188E-02_EB, 5.097597E-02_EB, 5.500449E-02_EB, &
    5.825744E-02_EB, 6.083116E-02_EB, 6.282015E-02_EB, 6.430721E-02_EB, 6.537213E-02_EB, 6.607949E-02_EB, &
    6.648679E-02_EB, 6.664418E-02_EB, 6.659220E-02_EB, 6.636777E-02_EB, 6.600167E-02_EB, 6.552032E-02_EB, &
    6.494261E-02_EB, 6.428840E-02_EB, 6.357194E-02_EB, 6.280715E-02_EB, 6.200504E-02_EB,  &
    1.153633E-01_EB, 1.245526E-01_EB, 1.285002E-01_EB, 1.289000E-01_EB, 1.269717E-01_EB, 1.235632E-01_EB, &
    1.192526E-01_EB, 1.144392E-01_EB, 1.093883E-01_EB, 1.042797E-01_EB, 9.923347E-02_EB, 9.432412E-02_EB, &
    8.960104E-02_EB, 8.509117E-02_EB, 8.080763E-02_EB, 7.675598E-02_EB, 7.293147E-02_EB, 6.933082E-02_EB, &
    6.594450E-02_EB, 6.276318E-02_EB, 5.977554E-02_EB, 5.696758E-02_EB, 5.433254E-02_EB,  &
    3.460355E-01_EB, 3.205380E-01_EB, 2.960906E-01_EB, 2.734307E-01_EB, 2.527266E-01_EB, 2.339252E-01_EB, &
    2.168922E-01_EB, 2.014663E-01_EB, 1.874880E-01_EB, 1.748046E-01_EB, 1.632788E-01_EB, 1.527883E-01_EB, &
    1.432227E-01_EB, 1.344836E-01_EB, 1.264841E-01_EB, 1.191485E-01_EB, 1.124100E-01_EB, 1.062065E-01_EB, &
    1.004861E-01_EB, 9.520356E-02_EB, 9.031470E-02_EB, 8.578244E-02_EB, 8.157744E-02_EB,  &
    4.804312E-01_EB, 4.183837E-01_EB, 3.677948E-01_EB, 3.259852E-01_EB, 2.909746E-01_EB, 2.613148E-01_EB, &
    2.359310E-01_EB, 2.140177E-01_EB, 1.949566E-01_EB, 1.782671E-01_EB, 1.635708E-01_EB, 1.505600E-01_EB, &
    1.389897E-01_EB, 1.286554E-01_EB, 1.193894E-01_EB, 1.110524E-01_EB, 1.035252E-01_EB, 9.670962E-02_EB, &
    9.051819E-02_EB, 8.487945E-02_EB, 7.973241E-02_EB, 7.501998E-02_EB, 7.069744E-02_EB,  &
    9.835833E-01_EB, 7.945568E-01_EB, 6.600323E-01_EB, 5.598722E-01_EB, 4.826554E-01_EB, 4.214727E-01_EB, &
    3.719185E-01_EB, 3.310569E-01_EB, 2.968658E-01_EB, 2.678966E-01_EB, 2.430942E-01_EB, 2.216657E-01_EB, &
    2.030050E-01_EB, 1.866420E-01_EB, 1.722059E-01_EB, 1.593976E-01_EB, 1.479767E-01_EB, 1.377484E-01_EB, &
    1.285486E-01_EB, 1.202428E-01_EB, 1.127186E-01_EB, 1.058789E-01_EB, 9.964358E-02_EB/),(/23,8/))

SD1_C3H6(1:23,9:16) = RESHAPE((/ &  ! 1450-1625 cm-1
    8.795663E-01_EB, 6.709369E-01_EB, 5.300948E-01_EB, 4.299970E-01_EB, 3.559907E-01_EB, 2.995492E-01_EB, &
    2.554174E-01_EB, 2.202028E-01_EB, 1.916273E-01_EB, 1.681072E-01_EB, 1.485122E-01_EB, 1.320175E-01_EB, &
    1.180034E-01_EB, 1.060008E-01_EB, 9.564792E-02_EB, 8.665880E-02_EB, 7.880988E-02_EB, 7.191897E-02_EB, &
    6.583953E-02_EB, 6.045221E-02_EB, 5.565931E-02_EB, 5.137750E-02_EB, 4.753956E-02_EB,  &
    6.350346E-01_EB, 5.238615E-01_EB, 4.408428E-01_EB, 3.768285E-01_EB, 3.261766E-01_EB, 2.852455E-01_EB, &
    2.515976E-01_EB, 2.235395E-01_EB, 1.998629E-01_EB, 1.796807E-01_EB, 1.623252E-01_EB, 1.472877E-01_EB, &
    1.341718E-01_EB, 1.226631E-01_EB, 1.125105E-01_EB, 1.035104E-01_EB, 9.549831E-02_EB, 8.833628E-02_EB, &
    8.190969E-02_EB, 7.612322E-02_EB, 7.089524E-02_EB, 6.616162E-02_EB, 6.185914E-02_EB,  &
    2.202119E-01_EB, 2.060176E-01_EB, 1.928546E-01_EB, 1.808220E-01_EB, 1.698573E-01_EB, 1.598559E-01_EB, &
    1.507093E-01_EB, 1.423249E-01_EB, 1.346197E-01_EB, 1.275202E-01_EB, 1.209659E-01_EB, 1.149033E-01_EB, &
    1.092849E-01_EB, 1.040664E-01_EB, 9.921580E-02_EB, 9.469788E-02_EB, 9.048339E-02_EB, 8.654652E-02_EB, &
    8.286434E-02_EB, 7.941573E-02_EB, 7.618158E-02_EB, 7.314464E-02_EB, 7.029001E-02_EB,  &
    9.542899E-02_EB, 1.014545E-01_EB, 1.056655E-01_EB, 1.085461E-01_EB, 1.104218E-01_EB, 1.115195E-01_EB, &
    1.120100E-01_EB, 1.120195E-01_EB, 1.116476E-01_EB, 1.109741E-01_EB, 1.100630E-01_EB, 1.089650E-01_EB, &
    1.077220E-01_EB, 1.063667E-01_EB, 1.049296E-01_EB, 1.034309E-01_EB, 1.018902E-01_EB, 1.003227E-01_EB, &
    9.874044E-02_EB, 9.715155E-02_EB, 9.556725E-02_EB, 9.399139E-02_EB, 9.243064E-02_EB,  &
    8.237186E-02_EB, 8.493470E-02_EB, 8.661658E-02_EB, 8.765113E-02_EB, 8.819174E-02_EB, 8.834497E-02_EB, &
    8.818983E-02_EB, 8.778874E-02_EB, 8.718846E-02_EB, 8.643351E-02_EB, 8.555281E-02_EB, 8.457473E-02_EB, &
    8.352155E-02_EB, 8.241295E-02_EB, 8.126317E-02_EB, 8.008358E-02_EB, 7.888705E-02_EB, 7.768167E-02_EB, &
    7.647201E-02_EB, 7.526675E-02_EB, 7.406861E-02_EB, 7.288216E-02_EB, 7.170974E-02_EB,  &
    9.910066E-02_EB, 1.005414E-01_EB, 1.014820E-01_EB, 1.020368E-01_EB, 1.022723E-01_EB, 1.022407E-01_EB, &
    1.019817E-01_EB, 1.015300E-01_EB, 1.009146E-01_EB, 1.001654E-01_EB, 9.930285E-02_EB, 9.834962E-02_EB, &
    9.732226E-02_EB, 9.623703E-02_EB, 9.510777E-02_EB, 9.394528E-02_EB, 9.275824E-02_EB, 9.155552E-02_EB, &
    9.034433E-02_EB, 8.913050E-02_EB, 8.791746E-02_EB, 8.671085E-02_EB, 8.551261E-02_EB,  &
    2.220822E-01_EB, 2.363820E-01_EB, 2.472837E-01_EB, 2.555994E-01_EB, 2.618760E-01_EB, 2.665030E-01_EB, &
    2.697771E-01_EB, 2.719284E-01_EB, 2.731450E-01_EB, 2.735819E-01_EB, 2.733685E-01_EB, 2.726131E-01_EB, &
    2.714065E-01_EB, 2.698269E-01_EB, 2.679405E-01_EB, 2.658037E-01_EB, 2.634625E-01_EB, 2.609580E-01_EB, &
    2.583236E-01_EB, 2.555880E-01_EB, 2.527763E-01_EB, 2.499093E-01_EB, 2.470017E-01_EB,  &
    8.877792E-01_EB, 7.242881E-01_EB, 6.076915E-01_EB, 5.206825E-01_EB, 4.534192E-01_EB, 3.999423E-01_EB, &
    3.564505E-01_EB, 3.204177E-01_EB, 2.901037E-01_EB, 2.642720E-01_EB, 2.420167E-01_EB, 2.226659E-01_EB, &
    2.057039E-01_EB, 1.907288E-01_EB, 1.774279E-01_EB, 1.655489E-01_EB, 1.548861E-01_EB, 1.452726E-01_EB, &
    1.365700E-01_EB, 1.286619E-01_EB, 1.214527E-01_EB, 1.148579E-01_EB, 1.088103E-01_EB/),(/23,8/))

SD1_C3H6(1:23,17:24) = RESHAPE((/ &  ! 1650-1825 cm-1
    9.361306E-01_EB, 7.419497E-01_EB, 6.064115E-01_EB, 5.071834E-01_EB, 4.318192E-01_EB, 3.728930E-01_EB, &
    3.257369E-01_EB, 2.872711E-01_EB, 2.553952E-01_EB, 2.286306E-01_EB, 2.059008E-01_EB, 1.864103E-01_EB, &
    1.695581E-01_EB, 1.548774E-01_EB, 1.420047E-01_EB, 1.306521E-01_EB, 1.205873E-01_EB, 1.116199E-01_EB, &
    1.035983E-01_EB, 9.639160E-02_EB, 8.989452E-02_EB, 8.401681E-02_EB, 7.868237E-02_EB,  &
    8.105935E-02_EB, 8.277922E-02_EB, 8.203160E-02_EB, 7.990058E-02_EB, 7.703610E-02_EB, 7.382203E-02_EB, &
    7.048559E-02_EB, 6.716131E-02_EB, 6.392661E-02_EB, 6.082398E-02_EB, 5.787481E-02_EB, 5.508982E-02_EB, &
    5.246957E-02_EB, 5.000978E-02_EB, 4.770399E-02_EB, 4.554420E-02_EB, 4.352174E-02_EB, 4.162794E-02_EB, &
    3.985372E-02_EB, 3.819042E-02_EB, 3.663012E-02_EB, 3.516562E-02_EB, 3.378864E-02_EB,  &
    6.606969E-03_EB, 5.631867E-03_EB, 4.930957E-03_EB, 4.401271E-03_EB, 3.985800E-03_EB, 3.647403E-03_EB, &
    3.366811E-03_EB, 3.128637E-03_EB, 2.923731E-03_EB, 2.743062E-03_EB, 2.584643E-03_EB, 2.442215E-03_EB, &
    2.314885E-03_EB, 2.198876E-03_EB, 2.093795E-03_EB, 1.996850E-03_EB, 1.908256E-03_EB, 1.826013E-03_EB, &
    1.750425E-03_EB, 1.678223E-03_EB, 1.613855E-03_EB, 1.552073E-03_EB, 1.494762E-03_EB,  &
    1.655318E-03_EB, 1.670286E-03_EB, 1.702375E-03_EB, 1.741163E-03_EB, 1.781307E-03_EB, 1.819887E-03_EB, &
    1.855244E-03_EB, 1.886795E-03_EB, 1.912896E-03_EB, 1.934802E-03_EB, 1.952228E-03_EB, 1.965758E-03_EB, &
    1.975189E-03_EB, 1.981712E-03_EB, 1.984135E-03_EB, 1.985299E-03_EB, 1.982658E-03_EB, 1.978857E-03_EB, &
    1.972222E-03_EB, 1.963353E-03_EB, 1.954485E-03_EB, 1.942791E-03_EB, 1.931686E-03_EB,  &
    3.809286E-05_EB, 1.672711E-04_EB, 5.096017E-04_EB, 1.209734E-03_EB, 2.416914E-03_EB, 4.250851E-03_EB, &
    6.793251E-03_EB, 1.008504E-02_EB, 1.412649E-02_EB, 1.888457E-02_EB, 2.430221E-02_EB, 3.031022E-02_EB, &
    3.682278E-02_EB, 4.375954E-02_EB, 5.103488E-02_EB, 5.856822E-02_EB, 6.628476E-02_EB, 7.411525E-02_EB, &
    8.200016E-02_EB, 8.988563E-02_EB, 9.772530E-02_EB, 1.054765E-01_EB, 1.131089E-01_EB,  &
    2.269251E-02_EB, 2.711702E-02_EB, 3.090471E-02_EB, 3.430491E-02_EB, 3.751968E-02_EB, 4.068651E-02_EB, &
    4.387267E-02_EB, 4.709767E-02_EB, 5.034279E-02_EB, 5.357460E-02_EB, 5.675325E-02_EB, 5.983588E-02_EB, &
    6.278745E-02_EB, 6.557738E-02_EB, 6.818631E-02_EB, 7.059691E-02_EB, 7.280086E-02_EB, 7.479465E-02_EB, &
    7.657930E-02_EB, 7.815813E-02_EB, 7.953963E-02_EB, 8.072988E-02_EB, 8.174037E-02_EB,  &
    1.877832E-01_EB, 1.568816E-01_EB, 1.334951E-01_EB, 1.152428E-01_EB, 1.006439E-01_EB, 8.873226E-02_EB, &
    7.885008E-02_EB, 7.053997E-02_EB, 6.347027E-02_EB, 5.739733E-02_EB, 5.213662E-02_EB, 4.754521E-02_EB, &
    4.351293E-02_EB, 3.995209E-02_EB, 3.679156E-02_EB, 3.397455E-02_EB, 3.145225E-02_EB, 2.918472E-02_EB, &
    2.713990E-02_EB, 2.529129E-02_EB, 2.361318E-02_EB, 2.208783E-02_EB, 2.069470E-02_EB,  &
    2.116948E-01_EB, 1.734868E-01_EB, 1.456017E-01_EB, 1.244878E-01_EB, 1.080169E-01_EB, 9.485045E-02_EB, &
    8.410990E-02_EB, 7.519998E-02_EB, 6.770203E-02_EB, 6.131704E-02_EB, 5.582198E-02_EB, 5.105453E-02_EB, &
    4.688190E-02_EB, 4.320635E-02_EB, 3.995113E-02_EB, 3.705173E-02_EB, 3.445669E-02_EB, 3.212219E-02_EB, &
    3.001713E-02_EB, 2.811005E-02_EB, 2.637679E-02_EB, 2.479795E-02_EB, 2.335289E-02_EB/),(/23,8/))

SD1_C3H6(1:23,25:29) = RESHAPE((/ &  ! 1850-1950 cm-1
    5.774795E-02_EB, 5.535830E-02_EB, 5.311844E-02_EB, 5.104058E-02_EB, 4.910913E-02_EB, 4.730734E-02_EB, &
    4.561875E-02_EB, 4.402894E-02_EB, 4.252760E-02_EB, 4.110349E-02_EB, 3.975237E-02_EB, 3.846563E-02_EB, &
    3.724047E-02_EB, 3.607211E-02_EB, 3.495615E-02_EB, 3.389086E-02_EB, 3.287204E-02_EB, 3.189863E-02_EB, &
    3.096847E-02_EB, 3.007834E-02_EB, 2.922456E-02_EB, 2.840810E-02_EB, 2.762612E-02_EB,  &
    1.557832E-02_EB, 1.482703E-02_EB, 1.412488E-02_EB, 1.347507E-02_EB, 1.287650E-02_EB, 1.232304E-02_EB, &
    1.180849E-02_EB, 1.132941E-02_EB, 1.088181E-02_EB, 1.046182E-02_EB, 1.006627E-02_EB, 9.693858E-03_EB, &
    9.341518E-03_EB, 9.008538E-03_EB, 8.693044E-03_EB, 8.393639E-03_EB, 8.110524E-03_EB, 7.840521E-03_EB, &
    7.584428E-03_EB, 7.340945E-03_EB, 7.108491E-03_EB, 6.886869E-03_EB, 6.676465E-03_EB,  &
    4.009069E-03_EB, 3.867844E-03_EB, 3.719408E-03_EB, 3.572717E-03_EB, 3.429429E-03_EB, 3.293552E-03_EB, &
    3.163506E-03_EB, 3.041057E-03_EB, 2.924358E-03_EB, 2.813501E-03_EB, 2.708005E-03_EB, 2.607862E-03_EB, &
    2.513173E-03_EB, 2.423652E-03_EB, 2.338312E-03_EB, 2.257064E-03_EB, 2.179522E-03_EB, 2.105190E-03_EB, &
    2.036022E-03_EB, 1.968803E-03_EB, 1.905969E-03_EB, 1.844691E-03_EB, 1.787309E-03_EB,  &
    2.353264E-03_EB, 2.559145E-03_EB, 2.727021E-03_EB, 2.866746E-03_EB, 2.981627E-03_EB, 3.076238E-03_EB, &
    3.154380E-03_EB, 3.217329E-03_EB, 3.268287E-03_EB, 3.307261E-03_EB, 3.336106E-03_EB, 3.357731E-03_EB, &
    3.371471E-03_EB, 3.378485E-03_EB, 3.382093E-03_EB, 3.379363E-03_EB, 3.372051E-03_EB, 3.362410E-03_EB, &
    3.349056E-03_EB, 3.333176E-03_EB, 3.315344E-03_EB, 3.295172E-03_EB, 3.273345E-03_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

ALLOCATE(SD2_C3H6(N_TEMP_C3H6,21)) 

! BAND #2: 2700 cm-1 - 3200 cm-1 

SD2_C3H6(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
    2.587896E-02_EB, 2.398045E-02_EB, 2.225585E-02_EB, 2.071310E-02_EB, 1.934140E-02_EB, 1.812515E-02_EB, &
    1.704045E-02_EB, 1.607259E-02_EB, 1.520148E-02_EB, 1.441710E-02_EB, 1.370374E-02_EB, 1.305275E-02_EB, &
    1.245778E-02_EB, 1.190948E-02_EB, 1.140298E-02_EB, 1.093380E-02_EB, 1.049834E-02_EB, 1.009162E-02_EB, &
    9.711401E-03_EB, 9.356418E-03_EB, 9.021219E-03_EB, 8.706386E-03_EB, 8.410937E-03_EB,  &
    6.548104E-02_EB, 5.426426E-02_EB, 4.616171E-02_EB, 4.005642E-02_EB, 3.530258E-02_EB, 3.150123E-02_EB, &
    2.839636E-02_EB, 2.581125E-02_EB, 2.362728E-02_EB, 2.175695E-02_EB, 2.013686E-02_EB, 1.871913E-02_EB, &
    1.746789E-02_EB, 1.635465E-02_EB, 1.535843E-02_EB, 1.446116E-02_EB, 1.364833E-02_EB, 1.290872E-02_EB, &
    1.223298E-02_EB, 1.161437E-02_EB, 1.104275E-02_EB, 1.051693E-02_EB, 1.002865E-02_EB,  &
    5.918988E-02_EB, 5.447961E-02_EB, 5.082039E-02_EB, 4.787680E-02_EB, 4.544455E-02_EB, 4.339035E-02_EB, &
    4.162033E-02_EB, 4.007265E-02_EB, 3.869783E-02_EB, 3.746413E-02_EB, 3.634175E-02_EB, 3.531413E-02_EB, &
    3.436317E-02_EB, 3.347795E-02_EB, 3.265048E-02_EB, 3.187093E-02_EB, 3.113365E-02_EB, 3.043319E-02_EB, &
    2.976711E-02_EB, 2.913036E-02_EB, 2.852194E-02_EB, 2.793788E-02_EB, 2.737622E-02_EB,  &
    6.341989E-02_EB, 6.089936E-02_EB, 5.862855E-02_EB, 5.658255E-02_EB, 5.473108E-02_EB, 5.304571E-02_EB, &
    5.150083E-02_EB, 5.007577E-02_EB, 4.875018E-02_EB, 4.750984E-02_EB, 4.634208E-02_EB, 4.523707E-02_EB, &
    4.418721E-02_EB, 4.318556E-02_EB, 4.222611E-02_EB, 4.130534E-02_EB, 4.042012E-02_EB, 3.956687E-02_EB, &
    3.874304E-02_EB, 3.794737E-02_EB, 3.717742E-02_EB, 3.643222E-02_EB, 3.570982E-02_EB,  &
    1.009020E-01_EB, 9.791175E-02_EB, 9.500413E-02_EB, 9.225220E-02_EB, 8.967341E-02_EB, 8.726523E-02_EB, &
    8.501309E-02_EB, 8.290059E-02_EB, 8.090779E-02_EB, 7.902456E-02_EB, 7.723414E-02_EB, 7.552729E-02_EB, &
    7.389357E-02_EB, 7.232500E-02_EB, 7.081595E-02_EB, 6.936026E-02_EB, 6.795385E-02_EB, 6.659390E-02_EB, &
    6.527670E-02_EB, 6.399982E-02_EB, 6.276093E-02_EB, 6.155905E-02_EB, 6.039067E-02_EB,  &
    2.066071E-01_EB, 2.061851E-01_EB, 2.052859E-01_EB, 2.041108E-01_EB, 2.027714E-01_EB, 2.013317E-01_EB, &
    1.998276E-01_EB, 1.982756E-01_EB, 1.966849E-01_EB, 1.950602E-01_EB, 1.934043E-01_EB, 1.917198E-01_EB, &
    1.900068E-01_EB, 1.882673E-01_EB, 1.865047E-01_EB, 1.847196E-01_EB, 1.829158E-01_EB, 1.810974E-01_EB, &
    1.792650E-01_EB, 1.774251E-01_EB, 1.755781E-01_EB, 1.737277E-01_EB, 1.718777E-01_EB,  &
    6.195013E-01_EB, 5.436008E-01_EB, 4.857281E-01_EB, 4.400313E-01_EB, 4.029495E-01_EB, 3.721839E-01_EB, &
    3.461864E-01_EB, 3.238739E-01_EB, 3.044664E-01_EB, 2.873895E-01_EB, 2.722115E-01_EB, 2.586037E-01_EB, &
    2.463086E-01_EB, 2.351259E-01_EB, 2.248937E-01_EB, 2.154836E-01_EB, 2.067885E-01_EB, 1.987233E-01_EB, &
    1.912140E-01_EB, 1.841990E-01_EB, 1.776297E-01_EB, 1.714615E-01_EB, 1.656531E-01_EB,  &
    1.082018E+00_EB, 9.388312E-01_EB, 8.305876E-01_EB, 7.457663E-01_EB, 6.774072E-01_EB, 6.210524E-01_EB, &
    5.737154E-01_EB, 5.333192E-01_EB, 4.983732E-01_EB, 4.677887E-01_EB, 4.407465E-01_EB, 4.166237E-01_EB, &
    3.949373E-01_EB, 3.753073E-01_EB, 3.574317E-01_EB, 3.410673E-01_EB, 3.260174E-01_EB, 3.121158E-01_EB, &
    2.992308E-01_EB, 2.872455E-01_EB, 2.760644E-01_EB, 2.656063E-01_EB, 2.558015E-01_EB/),(/23,8/))

SD2_C3H6(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
    1.685544E+00_EB, 1.454186E+00_EB, 1.280686E+00_EB, 1.145634E+00_EB, 1.037413E+00_EB, 9.486406E-01_EB, &
    8.743988E-01_EB, 8.112883E-01_EB, 7.568895E-01_EB, 7.094301E-01_EB, 6.675909E-01_EB, 6.303676E-01_EB, &
    5.969859E-01_EB, 5.668365E-01_EB, 5.394395E-01_EB, 5.144048E-01_EB, 4.914202E-01_EB, 4.702240E-01_EB, &
    4.506044E-01_EB, 4.323793E-01_EB, 4.153985E-01_EB, 3.995342E-01_EB, 3.846745E-01_EB,  &
    2.184683E+00_EB, 1.834252E+00_EB, 1.576727E+00_EB, 1.379896E+00_EB, 1.224783E+00_EB, 1.099504E+00_EB, &
    9.962479E-01_EB, 9.096830E-01_EB, 8.360487E-01_EB, 7.726255E-01_EB, 7.174028E-01_EB, 6.688607E-01_EB, &
    6.258370E-01_EB, 5.874227E-01_EB, 5.529022E-01_EB, 5.217017E-01_EB, 4.933586E-01_EB, 4.674924E-01_EB, &
    4.437904E-01_EB, 4.219904E-01_EB, 4.018756E-01_EB, 3.832571E-01_EB, 3.659770E-01_EB,  &
    2.696430E+00_EB, 2.206555E+00_EB, 1.853844E+00_EB, 1.589268E+00_EB, 1.384351E+00_EB, 1.221497E+00_EB, &
    1.089290E+00_EB, 9.800282E-01_EB, 8.883436E-01_EB, 8.103932E-01_EB, 7.433590E-01_EB, 6.851369E-01_EB, &
    6.341236E-01_EB, 5.890774E-01_EB, 5.490303E-01_EB, 5.132084E-01_EB, 4.809893E-01_EB, 4.518705E-01_EB, &
    4.254354E-01_EB, 4.013439E-01_EB, 3.793064E-01_EB, 3.590836E-01_EB, 3.404670E-01_EB,  &
    2.478186E+00_EB, 1.991151E+00_EB, 1.647750E+00_EB, 1.394590E+00_EB, 1.201383E+00_EB, 1.049787E+00_EB, &
    9.280994E-01_EB, 8.285417E-01_EB, 7.457584E-01_EB, 6.759617E-01_EB, 6.163997E-01_EB, 5.650325E-01_EB, &
    5.203230E-01_EB, 4.810890E-01_EB, 4.464094E-01_EB, 4.155584E-01_EB, 3.879543E-01_EB, 3.631277E-01_EB, &
    3.406928E-01_EB, 3.203367E-01_EB, 3.017938E-01_EB, 2.848433E-01_EB, 2.693022E-01_EB,  &
    1.868585E+00_EB, 1.575934E+00_EB, 1.360199E+00_EB, 1.194879E+00_EB, 1.064308E+00_EB, 9.586422E-01_EB, &
    8.713982E-01_EB, 7.981372E-01_EB, 7.357249E-01_EB, 6.818895E-01_EB, 6.349468E-01_EB, 5.936258E-01_EB, &
    5.569482E-01_EB, 5.241537E-01_EB, 4.946381E-01_EB, 4.679200E-01_EB, 4.436116E-01_EB, 4.213917E-01_EB, &
    4.009980E-01_EB, 3.822088E-01_EB, 3.648421E-01_EB, 3.487394E-01_EB, 3.337691E-01_EB,  &
    1.152911E+00_EB, 1.045788E+00_EB, 9.607596E-01_EB, 8.913222E-01_EB, 8.333256E-01_EB, 7.839823E-01_EB, &
    7.413436E-01_EB, 7.040041E-01_EB, 6.709249E-01_EB, 6.413223E-01_EB, 6.145927E-01_EB, 5.902689E-01_EB, &
    5.679794E-01_EB, 5.474300E-01_EB, 5.283852E-01_EB, 5.106492E-01_EB, 4.940653E-01_EB, 4.785030E-01_EB, &
    4.638495E-01_EB, 4.500169E-01_EB, 4.369225E-01_EB, 4.245022E-01_EB, 4.126975E-01_EB,  &
    1.248446E+00_EB, 1.114445E+00_EB, 1.010445E+00_EB, 9.271033E-01_EB, 8.586166E-01_EB, 8.011753E-01_EB, &
    7.521682E-01_EB, 7.097475E-01_EB, 6.725628E-01_EB, 6.396112E-01_EB, 6.101272E-01_EB, 5.835255E-01_EB, &
    5.593432E-01_EB, 5.372165E-01_EB, 5.168522E-01_EB, 4.980156E-01_EB, 4.805104E-01_EB, 4.641810E-01_EB, &
    4.488911E-01_EB, 4.345311E-01_EB, 4.210063E-01_EB, 4.082364E-01_EB, 3.961514E-01_EB,  &
    1.528788E+00_EB, 1.256299E+00_EB, 1.060062E+00_EB, 9.127049E-01_EB, 7.984009E-01_EB, 7.073946E-01_EB, &
    6.333702E-01_EB, 5.720691E-01_EB, 5.205228E-01_EB, 4.766063E-01_EB, 4.387596E-01_EB, 4.058168E-01_EB, &
    3.768887E-01_EB, 3.512895E-01_EB, 3.284781E-01_EB, 3.080256E-01_EB, 2.895877E-01_EB, 2.728856E-01_EB, &
    2.576847E-01_EB, 2.437986E-01_EB, 2.310656E-01_EB, 2.193530E-01_EB, 2.085443E-01_EB/),(/23,8/))

SD2_C3H6(1:23,17:21) = RESHAPE((/ &  ! 3100-3200 cm-1
    1.203135E+00_EB, 1.010739E+00_EB, 8.699824E-01_EB, 7.627625E-01_EB, 6.784870E-01_EB, 6.105595E-01_EB, &
    5.546666E-01_EB, 5.078694E-01_EB, 4.681049E-01_EB, 4.338834E-01_EB, 4.041041E-01_EB, 3.779368E-01_EB, &
    3.547442E-01_EB, 3.340357E-01_EB, 3.154184E-01_EB, 2.985822E-01_EB, 2.832754E-01_EB, 2.692929E-01_EB, &
    2.564647E-01_EB, 2.446505E-01_EB, 2.337312E-01_EB, 2.236097E-01_EB, 2.142000E-01_EB,  &
    3.451842E-01_EB, 3.128765E-01_EB, 2.873847E-01_EB, 2.666566E-01_EB, 2.493968E-01_EB, 2.347499E-01_EB, &
    2.221175E-01_EB, 2.110739E-01_EB, 2.013050E-01_EB, 1.925735E-01_EB, 1.846985E-01_EB, 1.775379E-01_EB, &
    1.709819E-01_EB, 1.649412E-01_EB, 1.593451E-01_EB, 1.541353E-01_EB, 1.492649E-01_EB, 1.446930E-01_EB, &
    1.403896E-01_EB, 1.363246E-01_EB, 1.324758E-01_EB, 1.288239E-01_EB, 1.253501E-01_EB,  &
    1.093789E-01_EB, 1.045930E-01_EB, 1.006339E-01_EB, 9.728149E-02_EB, 9.438645E-02_EB, 9.184667E-02_EB, &
    8.958668E-02_EB, 8.755155E-02_EB, 8.569648E-02_EB, 8.399070E-02_EB, 8.240509E-02_EB, 8.091989E-02_EB, &
    7.951959E-02_EB, 7.819068E-02_EB, 7.692194E-02_EB, 7.570492E-02_EB, 7.453385E-02_EB, 7.340329E-02_EB, &
    7.230710E-02_EB, 7.124285E-02_EB, 7.020811E-02_EB, 6.920004E-02_EB, 6.821593E-02_EB,  &
    3.015706E-02_EB, 2.989298E-02_EB, 2.966172E-02_EB, 2.945668E-02_EB, 2.926926E-02_EB, 2.909860E-02_EB, &
    2.893866E-02_EB, 2.878661E-02_EB, 2.864098E-02_EB, 2.849770E-02_EB, 2.835530E-02_EB, 2.821241E-02_EB, &
    2.806845E-02_EB, 2.792147E-02_EB, 2.777137E-02_EB, 2.761796E-02_EB, 2.745989E-02_EB, 2.729859E-02_EB, &
    2.713428E-02_EB, 2.696461E-02_EB, 2.679222E-02_EB, 2.661603E-02_EB, 2.643731E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

ALLOCATE(SD3_C3H6(N_TEMP_C3H6,19)) 

! BAND #3: 700 cm-1 - 1150 cm-1 

SD3_C3H6(1:23,1:8) = RESHAPE((/ &  ! 700-875 cm-1
    1.076543E-04_EB, 5.863694E-04_EB, 2.101445E-03_EB, 5.666284E-03_EB, 1.250105E-02_EB, 2.380933E-02_EB, &
    4.059015E-02_EB, 6.352911E-02_EB, 9.296029E-02_EB, 1.288911E-01_EB, 1.710584E-01_EB, 2.189940E-01_EB, &
    2.720883E-01_EB, 3.296561E-01_EB, 3.909743E-01_EB, 4.553230E-01_EB, 5.220117E-01_EB, 5.903906E-01_EB, &
    6.598690E-01_EB, 7.299116E-01_EB, 8.000497E-01_EB, 8.698692E-01_EB, 9.390180E-01_EB,  &
    2.134918E-04_EB, 9.098435E-04_EB, 2.740934E-03_EB, 6.508258E-03_EB, 1.303446E-02_EB, 2.301361E-02_EB, &
    3.691789E-02_EB, 5.497439E-02_EB, 7.717564E-02_EB, 1.033181E-01_EB, 1.330668E-01_EB, 1.659930E-01_EB, &
    2.016203E-01_EB, 2.394563E-01_EB, 2.790184E-01_EB, 3.198518E-01_EB, 3.615373E-01_EB, 4.036951E-01_EB, &
    4.459941E-01_EB, 4.881431E-01_EB, 5.298956E-01_EB, 5.710411E-01_EB, 6.114069E-01_EB,  &
    4.134157E-04_EB, 1.560348E-03_EB, 4.273547E-03_EB, 9.427166E-03_EB, 1.783938E-02_EB, 3.014756E-02_EB, &
    4.674748E-02_EB, 6.778572E-02_EB, 9.318024E-02_EB, 1.226720E-01_EB, 1.558773E-01_EB, 1.923300E-01_EB, &
    2.315236E-01_EB, 2.729432E-01_EB, 3.160884E-01_EB, 3.604863E-01_EB, 4.057037E-01_EB, 4.513526E-01_EB, &
    4.970890E-01_EB, 5.426129E-01_EB, 5.876717E-01_EB, 6.320491E-01_EB, 6.755655E-01_EB,  &
    6.060969E-04_EB, 2.024227E-03_EB, 5.056621E-03_EB, 1.033987E-02_EB, 1.831455E-02_EB, 2.918194E-02_EB, &
    4.292187E-02_EB, 5.933911E-02_EB, 7.812567E-02_EB, 9.891197E-02_EB, 1.213029E-01_EB, 1.449114E-01_EB, &
    1.693757E-01_EB, 1.943651E-01_EB, 2.195929E-01_EB, 2.448092E-01_EB, 2.698066E-01_EB, 2.944108E-01_EB, &
    3.184846E-01_EB, 3.419188E-01_EB, 3.646249E-01_EB, 3.865430E-01_EB, 4.076246E-01_EB,  &
    3.681383E-03_EB, 9.882239E-03_EB, 2.085049E-02_EB, 3.724294E-02_EB, 5.900830E-02_EB, 8.555691E-02_EB, &
    1.159819E-01_EB, 1.492571E-01_EB, 1.843713E-01_EB, 2.204102E-01_EB, 2.566015E-01_EB, 2.923147E-01_EB, &
    3.270612E-01_EB, 3.604812E-01_EB, 3.923143E-01_EB, 4.223925E-01_EB, 4.506168E-01_EB, 4.769408E-01_EB, &
    5.013659E-01_EB, 5.239186E-01_EB, 5.446518E-01_EB, 5.636344E-01_EB, 5.809443E-01_EB,  &
    7.311227E-02_EB, 1.051986E-01_EB, 1.356778E-01_EB, 1.626784E-01_EB, 1.853575E-01_EB, 2.035540E-01_EB, &
    2.175042E-01_EB, 2.276462E-01_EB, 2.344917E-01_EB, 2.385620E-01_EB, 2.403409E-01_EB, 2.402593E-01_EB, &
    2.386902E-01_EB, 2.359489E-01_EB, 2.323013E-01_EB, 2.279628E-01_EB, 2.231156E-01_EB, 2.179022E-01_EB, &
    2.124423E-01_EB, 2.068285E-01_EB, 2.011353E-01_EB, 1.954262E-01_EB, 1.897430E-01_EB,  &
    5.047291E-01_EB, 5.062757E-01_EB, 5.003725E-01_EB, 4.892699E-01_EB, 4.747690E-01_EB, 4.582211E-01_EB, &
    4.405998E-01_EB, 4.225886E-01_EB, 4.046490E-01_EB, 3.870892E-01_EB, 3.701049E-01_EB, 3.538125E-01_EB, &
    3.382753E-01_EB, 3.235199E-01_EB, 3.095455E-01_EB, 2.963363E-01_EB, 2.838674E-01_EB, 2.721039E-01_EB, &
    2.610108E-01_EB, 2.505511E-01_EB, 2.406874E-01_EB, 2.313812E-01_EB, 2.225976E-01_EB,  &
    1.432416E+00_EB, 1.222431E+00_EB, 1.061034E+00_EB, 9.326660E-01_EB, 8.280165E-01_EB, 7.411094E-01_EB, &
    6.678807E-01_EB, 6.054502E-01_EB, 5.517040E-01_EB, 5.050445E-01_EB, 4.642437E-01_EB, 4.283347E-01_EB, &
    3.965501E-01_EB, 3.682698E-01_EB, 3.429893E-01_EB, 3.202919E-01_EB, 2.998327E-01_EB, 2.813241E-01_EB, &
    2.645211E-01_EB, 2.492202E-01_EB, 2.352428E-01_EB, 2.224407E-01_EB, 2.106839E-01_EB/),(/23,8/))

SD3_C3H6(1:23,9:16) = RESHAPE((/ &  ! 900-1075 cm-1
    2.625468E+00_EB, 2.194945E+00_EB, 1.869328E+00_EB, 1.614846E+00_EB, 1.411009E+00_EB, 1.244592E+00_EB, &
    1.106624E+00_EB, 9.907838E-01_EB, 8.924726E-01_EB, 8.082634E-01_EB, 7.355415E-01_EB, 6.722919E-01_EB, &
    6.169214E-01_EB, 5.681663E-01_EB, 5.250068E-01_EB, 4.866153E-01_EB, 4.523113E-01_EB, 4.215327E-01_EB, &
    3.938130E-01_EB, 3.687566E-01_EB, 3.460337E-01_EB, 3.253622E-01_EB, 3.065029E-01_EB,  &
    1.997010E+00_EB, 1.678672E+00_EB, 1.437948E+00_EB, 1.249234E+00_EB, 1.097350E+00_EB, 9.726505E-01_EB, &
    8.686572E-01_EB, 7.808292E-01_EB, 7.058671E-01_EB, 6.413086E-01_EB, 5.852754E-01_EB, 5.363079E-01_EB, &
    4.932512E-01_EB, 4.551825E-01_EB, 4.213558E-01_EB, 3.911608E-01_EB, 3.640940E-01_EB, 3.397355E-01_EB, &
    3.177366E-01_EB, 2.978019E-01_EB, 2.796800E-01_EB, 2.631587E-01_EB, 2.480553E-01_EB,  &
    1.283517E+00_EB, 1.178267E+00_EB, 1.088862E+00_EB, 1.011138E+00_EB, 9.425480E-01_EB, 8.813984E-01_EB, &
    8.264883E-01_EB, 7.769105E-01_EB, 7.319479E-01_EB, 6.910218E-01_EB, 6.536485E-01_EB, 6.194186E-01_EB, &
    5.879843E-01_EB, 5.590446E-01_EB, 5.323382E-01_EB, 5.076374E-01_EB, 4.847452E-01_EB, 4.634844E-01_EB, &
    4.437021E-01_EB, 4.252605E-01_EB, 4.080394E-01_EB, 3.919308E-01_EB, 3.768388E-01_EB,  &
    1.098432E+00_EB, 9.763396E-01_EB, 8.796099E-01_EB, 8.001406E-01_EB, 7.331764E-01_EB, 6.757110E-01_EB, &
    6.257276E-01_EB, 5.817940E-01_EB, 5.428554E-01_EB, 5.081080E-01_EB, 4.769205E-01_EB, 4.487874E-01_EB, &
    4.232975E-01_EB, 4.001124E-01_EB, 3.789465E-01_EB, 3.595616E-01_EB, 3.417564E-01_EB, 3.253547E-01_EB, &
    3.102076E-01_EB, 2.961855E-01_EB, 2.831762E-01_EB, 2.710790E-01_EB, 2.598089E-01_EB,  &
    6.280058E-01_EB, 5.618064E-01_EB, 5.079883E-01_EB, 4.629333E-01_EB, 4.244373E-01_EB, 3.910628E-01_EB, &
    3.618097E-01_EB, 3.359497E-01_EB, 3.129351E-01_EB, 2.923350E-01_EB, 2.738070E-01_EB, 2.570688E-01_EB, &
    2.418913E-01_EB, 2.280804E-01_EB, 2.154723E-01_EB, 2.039295E-01_EB, 1.933322E-01_EB, 1.835770E-01_EB, &
    1.745760E-01_EB, 1.662530E-01_EB, 1.585392E-01_EB, 1.513760E-01_EB, 1.447112E-01_EB,  &
    3.347405E-01_EB, 3.121675E-01_EB, 2.928571E-01_EB, 2.758632E-01_EB, 2.606378E-01_EB, 2.468392E-01_EB, &
    2.342398E-01_EB, 2.226755E-01_EB, 2.120207E-01_EB, 2.021749E-01_EB, 1.930543E-01_EB, 1.845874E-01_EB, &
    1.767133E-01_EB, 1.693757E-01_EB, 1.625279E-01_EB, 1.561259E-01_EB, 1.501327E-01_EB, 1.445118E-01_EB, &
    1.392359E-01_EB, 1.342736E-01_EB, 1.296032E-01_EB, 1.251983E-01_EB, 1.210412E-01_EB,  &
    1.066212E-01_EB, 1.054359E-01_EB, 1.038046E-01_EB, 1.018578E-01_EB, 9.969273E-02_EB, 9.738655E-02_EB, &
    9.500040E-02_EB, 9.257731E-02_EB, 9.015422E-02_EB, 8.775376E-02_EB, 8.539638E-02_EB, 8.309420E-02_EB, &
    8.085461E-02_EB, 7.868532E-02_EB, 7.658829E-02_EB, 7.456489E-02_EB, 7.261600E-02_EB, 7.074026E-02_EB, &
    6.893621E-02_EB, 6.720112E-02_EB, 6.553530E-02_EB, 6.393522E-02_EB, 6.239623E-02_EB,  &
    3.197637E-02_EB, 3.150557E-02_EB, 3.091345E-02_EB, 3.023879E-02_EB, 2.951277E-02_EB, 2.875654E-02_EB, &
    2.798627E-02_EB, 2.721424E-02_EB, 2.644932E-02_EB, 2.569764E-02_EB, 2.496409E-02_EB, 2.425274E-02_EB, &
    2.356371E-02_EB, 2.289854E-02_EB, 2.225919E-02_EB, 2.164312E-02_EB, 2.105239E-02_EB, 2.048523E-02_EB, &
    1.994116E-02_EB, 1.941882E-02_EB, 1.891869E-02_EB, 1.843854E-02_EB, 1.797807E-02_EB/),(/23,8/))

SD3_C3H6(1:23,17:19) = RESHAPE((/ &  ! 1100-1150 cm-1
    3.561877E-03_EB, 3.616052E-03_EB, 3.640603E-03_EB, 3.642157E-03_EB, 3.627829E-03_EB, 3.598687E-03_EB, &
    3.560864E-03_EB, 3.515054E-03_EB, 3.463698E-03_EB, 3.408822E-03_EB, 3.351226E-03_EB, 3.292650E-03_EB, &
    3.230572E-03_EB, 3.171707E-03_EB, 3.110991E-03_EB, 3.051734E-03_EB, 2.992286E-03_EB, 2.934199E-03_EB, &
    2.878552E-03_EB, 2.823195E-03_EB, 2.768721E-03_EB, 2.715898E-03_EB, 2.664728E-03_EB,  &
    6.878183E-03_EB, 7.337399E-03_EB, 7.680223E-03_EB, 7.931005E-03_EB, 8.106697E-03_EB, 8.222283E-03_EB, &
    8.289862E-03_EB, 8.317408E-03_EB, 8.314763E-03_EB, 8.286200E-03_EB, 8.238563E-03_EB, 8.175919E-03_EB, &
    8.101007E-03_EB, 8.016642E-03_EB, 7.923611E-03_EB, 7.826099E-03_EB, 7.724889E-03_EB, 7.622221E-03_EB, &
    7.515941E-03_EB, 7.410061E-03_EB, 7.303107E-03_EB, 7.197118E-03_EB, 7.090262E-03_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,3/))

!-------------------------MMA DATA-------------------


! THERE ARE 2 BANDS FOR MMA

! BAND #1: 2700 cm-1 - 3200 cm-1 
! BAND #2: 750 cm-1 - 1950 cm-1 

N_TEMP_MMA = 23
N_BAND_MMA = 2

ALLOCATE(SD_MMA_TEMP(N_TEMP_MMA)) 

! INITIALIZE BANDS WAVENUMBER BOUNDS FOR MMA ABSOPTION COEFFICIENTS.
! BANDS ARE ORGANIZED BY ROW: 
! 1ST COLUMN IS THE LOWER BOUND BAND LIMIT IN cm-1. 
! 2ND COLUMN IS THE UPPER BOUND BAND LIMIT IN cm-1. 
! 3RD COLUMN IS THE STRIDE BETWEEN WAVENUMBERS IN BAND.
! IF 3RD COLUMN = 0., THEN THE BAND IS CALCULATED AND NOT TABULATED.

ALLOCATE(OM_BND_MMA(N_BAND_MMA,3)) 

OM_BND_MMA = RESHAPE((/ &
    2700._EB,  750._EB, &
    3200._EB, 1950._EB, &
     25._EB, 25._EB/),(/N_BAND_MMA,3/)) 

SD_MMA_TEMP = (/ &
    300._EB, 350._EB, 400._EB, 450._EB, 500._EB, 550._EB,&
    600._EB, 650._EB, 700._EB, 750._EB, 800._EB, 850._EB,&
    900._EB, 950._EB, 1000._EB, 1050._EB, 1100._EB, 1150._EB,&
    1200._EB, 1250._EB, 1300._EB, 1350._EB, 1400._EB/)

ALLOCATE(SD1_MMA(N_TEMP_MMA,21)) 

! BAND #1: 2700 cm-1 - 3200 cm-1 

SD1_MMA(1:23,1:8) = RESHAPE((/ &  ! 2700-2875 cm-1
    4.834464E-03_EB, 8.147664E-03_EB, 1.174120E-02_EB, 1.529128E-02_EB, 1.859299E-02_EB, 2.153801E-02_EB, &
    2.408446E-02_EB, 2.623323E-02_EB, 2.800491E-02_EB, 2.943352E-02_EB, 3.055627E-02_EB, 3.141079E-02_EB, &
    3.203295E-02_EB, 3.245449E-02_EB, 3.270408E-02_EB, 3.280842E-02_EB, 3.279106E-02_EB, 3.267006E-02_EB, &
    3.246391E-02_EB, 3.218587E-02_EB, 3.184890E-02_EB, 3.146341E-02_EB, 3.104139E-02_EB,  &
    1.017707E-02_EB, 1.588161E-02_EB, 2.147095E-02_EB, 2.647555E-02_EB, 3.068371E-02_EB, 3.405663E-02_EB, &
    3.664340E-02_EB, 3.852804E-02_EB, 3.981774E-02_EB, 4.060724E-02_EB, 4.098664E-02_EB, 4.103567E-02_EB, &
    4.081890E-02_EB, 4.039274E-02_EB, 3.980417E-02_EB, 3.909034E-02_EB, 3.828292E-02_EB, 3.740566E-02_EB, &
    3.648100E-02_EB, 3.552474E-02_EB, 3.455009E-02_EB, 3.356936E-02_EB, 3.258812E-02_EB,  &
    1.581069E-02_EB, 2.200694E-02_EB, 2.751488E-02_EB, 3.212081E-02_EB, 3.581086E-02_EB, 3.866680E-02_EB, &
    4.080007E-02_EB, 4.232531E-02_EB, 4.335201E-02_EB, 4.397173E-02_EB, 4.425989E-02_EB, 4.428370E-02_EB, &
    4.409547E-02_EB, 4.373896E-02_EB, 4.324815E-02_EB, 4.265356E-02_EB, 4.197994E-02_EB, 4.124461E-02_EB, &
    4.046493E-02_EB, 3.965368E-02_EB, 3.881998E-02_EB, 3.797429E-02_EB, 3.712166E-02_EB,  &
    1.587352E-02_EB, 2.276450E-02_EB, 2.919974E-02_EB, 3.485746E-02_EB, 3.963860E-02_EB, 4.356413E-02_EB, &
    4.670764E-02_EB, 4.916542E-02_EB, 5.103207E-02_EB, 5.239891E-02_EB, 5.334389E-02_EB, 5.393764E-02_EB, &
    5.423766E-02_EB, 5.429636E-02_EB, 5.415281E-02_EB, 5.384355E-02_EB, 5.339928E-02_EB, 5.284668E-02_EB, &
    5.220397E-02_EB, 5.149279E-02_EB, 5.072384E-02_EB, 4.991292E-02_EB, 4.906889E-02_EB,  &
    4.265960E-02_EB, 5.683681E-02_EB, 6.844221E-02_EB, 7.729705E-02_EB, 8.365317E-02_EB, 8.790686E-02_EB, &
    9.046738E-02_EB, 9.169838E-02_EB, 9.190407E-02_EB, 9.132900E-02_EB, 9.016778E-02_EB, 8.857492E-02_EB, &
    8.666893E-02_EB, 8.454304E-02_EB, 8.226768E-02_EB, 7.989957E-02_EB, 7.748257E-02_EB, 7.504784E-02_EB, &
    7.262267E-02_EB, 7.022432E-02_EB, 6.786746E-02_EB, 6.556373E-02_EB, 6.331914E-02_EB,  &
    3.051354E-01_EB, 3.000208E-01_EB, 2.917036E-01_EB, 2.812331E-01_EB, 2.694998E-01_EB, 2.571722E-01_EB, &
    2.447142E-01_EB, 2.324305E-01_EB, 2.205136E-01_EB, 2.090773E-01_EB, 1.981824E-01_EB, 1.878581E-01_EB, &
    1.781063E-01_EB, 1.689171E-01_EB, 1.602702E-01_EB, 1.521443E-01_EB, 1.445093E-01_EB, 1.373397E-01_EB, &
    1.306064E-01_EB, 1.242811E-01_EB, 1.183391E-01_EB, 1.127540E-01_EB, 1.075030E-01_EB,  &
    6.482805E-01_EB, 5.553303E-01_EB, 4.831665E-01_EB, 4.259525E-01_EB, 3.797084E-01_EB, 3.416880E-01_EB, &
    3.099513E-01_EB, 2.830997E-01_EB, 2.601095E-01_EB, 2.402166E-01_EB, 2.228392E-01_EB, 2.075349E-01_EB, &
    1.939546E-01_EB, 1.818216E-01_EB, 1.709197E-01_EB, 1.610703E-01_EB, 1.521289E-01_EB, 1.439781E-01_EB, &
    1.365192E-01_EB, 1.296675E-01_EB, 1.233550E-01_EB, 1.175205E-01_EB, 1.121147E-01_EB,  &
    4.611878E-01_EB, 4.411914E-01_EB, 4.219957E-01_EB, 4.040890E-01_EB, 3.875657E-01_EB, 3.723650E-01_EB, &
    3.583690E-01_EB, 3.454445E-01_EB, 3.334641E-01_EB, 3.223134E-01_EB, 3.118921E-01_EB, 3.021136E-01_EB, &
    2.929079E-01_EB, 2.842106E-01_EB, 2.759725E-01_EB, 2.681482E-01_EB, 2.607003E-01_EB, 2.535965E-01_EB, &
    2.468112E-01_EB, 2.403186E-01_EB, 2.340995E-01_EB, 2.281357E-01_EB, 2.224085E-01_EB/),(/23,8/))

SD1_MMA(1:23,9:16) = RESHAPE((/ &  ! 2900-3075 cm-1
    8.374059E-01_EB, 7.671901E-01_EB, 7.086493E-01_EB, 6.591593E-01_EB, 6.167636E-01_EB, 5.799963E-01_EB, &
    5.477535E-01_EB, 5.191911E-01_EB, 4.936552E-01_EB, 4.706388E-01_EB, 4.497400E-01_EB, 4.306425E-01_EB, &
    4.130860E-01_EB, 3.968655E-01_EB, 3.818113E-01_EB, 3.677831E-01_EB, 3.546647E-01_EB, 3.423599E-01_EB, &
    3.307854E-01_EB, 3.198724E-01_EB, 3.095609E-01_EB, 2.997964E-01_EB, 2.905357E-01_EB,  &
    1.627006E+00_EB, 1.512521E+00_EB, 1.419392E+00_EB, 1.341677E+00_EB, 1.275493E+00_EB, 1.218162E+00_EB, &
    1.167779E+00_EB, 1.122941E+00_EB, 1.082593E+00_EB, 1.045935E+00_EB, 1.012342E+00_EB, 9.813266E-01_EB, &
    9.525025E-01_EB, 9.255628E-01_EB, 9.002594E-01_EB, 8.763902E-01_EB, 8.537937E-01_EB, 8.323318E-01_EB, &
    8.118933E-01_EB, 7.923800E-01_EB, 7.737167E-01_EB, 7.558322E-01_EB, 7.386668E-01_EB,  &
    3.819288E+00_EB, 3.231655E+00_EB, 2.798847E+00_EB, 2.467269E+00_EB, 2.205342E+00_EB, 1.993272E+00_EB, &
    1.818048E+00_EB, 1.670769E+00_EB, 1.545161E+00_EB, 1.436680E+00_EB, 1.341958E+00_EB, 1.258456E+00_EB, &
    1.184225E+00_EB, 1.117743E+00_EB, 1.057810E+00_EB, 1.003466E+00_EB, 9.539344E-01_EB, 9.085798E-01_EB, &
    8.668781E-01_EB, 8.283920E-01_EB, 7.927543E-01_EB, 7.596519E-01_EB, 7.288244E-01_EB,  &
    2.770726E+00_EB, 2.283527E+00_EB, 1.932937E+00_EB, 1.669696E+00_EB, 1.465413E+00_EB, 1.302641E+00_EB, &
    1.170100E+00_EB, 1.060195E+00_EB, 9.676397E-01_EB, 8.886537E-01_EB, 8.204631E-01_EB, 7.609985E-01_EB, &
    7.086822E-01_EB, 6.622935E-01_EB, 6.208760E-01_EB, 5.836709E-01_EB, 5.500651E-01_EB, 5.195599E-01_EB, &
    4.917490E-01_EB, 4.662920E-01_EB, 4.429071E-01_EB, 4.213555E-01_EB, 4.014332E-01_EB,  &
    2.052261E+00_EB, 1.675352E+00_EB, 1.405754E+00_EB, 1.204469E+00_EB, 1.049104E+00_EB, 9.259405E-01_EB, &
    8.261451E-01_EB, 7.437845E-01_EB, 6.747452E-01_EB, 6.160910E-01_EB, 5.656775E-01_EB, 5.219035E-01_EB, &
    4.835529E-01_EB, 4.496897E-01_EB, 4.195791E-01_EB, 3.926387E-01_EB, 3.683997E-01_EB, 3.464821E-01_EB, &
    3.265760E-01_EB, 3.084232E-01_EB, 2.918080E-01_EB, 2.765489E-01_EB, 2.624951E-01_EB,  &
    1.127518E+00_EB, 9.394817E-01_EB, 8.033166E-01_EB, 7.004966E-01_EB, 6.202906E-01_EB, 5.560734E-01_EB, &
    5.035428E-01_EB, 4.597918E-01_EB, 4.227894E-01_EB, 3.910805E-01_EB, 3.635924E-01_EB, 3.395224E-01_EB, &
    3.182610E-01_EB, 2.993304E-01_EB, 2.823589E-01_EB, 2.670495E-01_EB, 2.531653E-01_EB, 2.405087E-01_EB, &
    2.289216E-01_EB, 2.182702E-01_EB, 2.084460E-01_EB, 1.993520E-01_EB, 1.909123E-01_EB,  &
    3.604821E-01_EB, 3.446453E-01_EB, 3.313076E-01_EB, 3.198650E-01_EB, 3.098925E-01_EB, 3.010808E-01_EB, &
    2.931966E-01_EB, 2.860630E-01_EB, 2.795408E-01_EB, 2.735216E-01_EB, 2.679213E-01_EB, 2.626678E-01_EB, &
    2.577102E-01_EB, 2.530023E-01_EB, 2.485120E-01_EB, 2.442077E-01_EB, 2.400698E-01_EB, 2.360782E-01_EB, &
    2.322188E-01_EB, 2.284786E-01_EB, 2.248490E-01_EB, 2.213189E-01_EB, 2.178844E-01_EB,  &
    2.241879E-01_EB, 2.324233E-01_EB, 2.388385E-01_EB, 2.439662E-01_EB, 2.481421E-01_EB, 2.515909E-01_EB, &
    2.544638E-01_EB, 2.568632E-01_EB, 2.588646E-01_EB, 2.605208E-01_EB, 2.618737E-01_EB, 2.629530E-01_EB, &
    2.637859E-01_EB, 2.643932E-01_EB, 2.647939E-01_EB, 2.650036E-01_EB, 2.650397E-01_EB, 2.649137E-01_EB, &
    2.646398E-01_EB, 2.642286E-01_EB, 2.636918E-01_EB, 2.630417E-01_EB, 2.622865E-01_EB/),(/23,8/))

SD1_MMA(1:23,17:21) = RESHAPE((/ &  ! 3100-3200 cm-1
    4.282367E-01_EB, 3.941938E-01_EB, 3.678239E-01_EB, 3.466912E-01_EB, 3.293044E-01_EB, 3.146858E-01_EB, &
    3.021716E-01_EB, 2.912907E-01_EB, 2.817010E-01_EB, 2.731468E-01_EB, 2.654337E-01_EB, 2.584109E-01_EB, &
    2.519623E-01_EB, 2.459986E-01_EB, 2.404435E-01_EB, 2.352388E-01_EB, 2.303376E-01_EB, 2.257011E-01_EB, &
    2.212976E-01_EB, 2.171016E-01_EB, 2.130908E-01_EB, 2.092462E-01_EB, 2.055544E-01_EB,  &
    1.096342E-01_EB, 1.102533E-01_EB, 1.105792E-01_EB, 1.107234E-01_EB, 1.107520E-01_EB, 1.106984E-01_EB, &
    1.105864E-01_EB, 1.104256E-01_EB, 1.102221E-01_EB, 1.099798E-01_EB, 1.097014E-01_EB, 1.093880E-01_EB, &
    1.090401E-01_EB, 1.086583E-01_EB, 1.082439E-01_EB, 1.077979E-01_EB, 1.073213E-01_EB, 1.068171E-01_EB, &
    1.062859E-01_EB, 1.057294E-01_EB, 1.051514E-01_EB, 1.045509E-01_EB, 1.039317E-01_EB,  &
    4.689425E-02_EB, 4.941880E-02_EB, 5.137056E-02_EB, 5.291564E-02_EB, 5.416310E-02_EB, 5.518302E-02_EB, &
    5.602720E-02_EB, 5.672808E-02_EB, 5.731110E-02_EB, 5.779476E-02_EB, 5.819191E-02_EB, 5.851473E-02_EB, &
    5.876869E-02_EB, 5.896361E-02_EB, 5.910409E-02_EB, 5.919595E-02_EB, 5.924261E-02_EB, 5.924797E-02_EB, &
    5.921534E-02_EB, 5.914940E-02_EB, 5.905238E-02_EB, 5.892604E-02_EB, 5.877447E-02_EB,  &
    1.927148E-02_EB, 2.111334E-02_EB, 2.259913E-02_EB, 2.381615E-02_EB, 2.482722E-02_EB, 2.567715E-02_EB, &
    2.639654E-02_EB, 2.700886E-02_EB, 2.753380E-02_EB, 2.798352E-02_EB, 2.836768E-02_EB, 2.869669E-02_EB, &
    2.897447E-02_EB, 2.920918E-02_EB, 2.940425E-02_EB, 2.956267E-02_EB, 2.968905E-02_EB, 2.978562E-02_EB, &
    2.985598E-02_EB, 2.990189E-02_EB, 2.992675E-02_EB, 2.993009E-02_EB, 2.991569E-02_EB,  &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,5/))

ALLOCATE(SD2_MMA(N_TEMP_MMA,49)) 

! BAND #2: 750 cm-1 - 1950 cm-1 

SD2_MMA(1:23,1:8) = RESHAPE((/ &  ! 750-925 cm-1
    1.025137E-03_EB, 3.895639E-03_EB, 1.122501E-02_EB, 2.664816E-02_EB, 5.467756E-02_EB, 1.000990E-01_EB, &
    1.672995E-01_EB, 2.597500E-01_EB, 3.797101E-01_EB, 5.281678E-01_EB, 7.049270E-01_EB, 9.088045E-01_EB, &
    1.137852E+00_EB, 1.389593E+00_EB, 1.661222E+00_EB, 1.949778E+00_EB, 2.252281E+00_EB, 2.565828E+00_EB, &
    2.887671E+00_EB, 3.215257E+00_EB, 3.546257E+00_EB, 3.878580E+00_EB, 4.210370E+00_EB,  &
    4.751271E-02_EB, 7.616567E-02_EB, 1.126448E-01_EB, 1.575400E-01_EB, 2.119847E-01_EB, 2.776398E-01_EB, &
    3.564324E-01_EB, 4.502131E-01_EB, 5.604663E-01_EB, 6.881310E-01_EB, 8.335303E-01_EB, 9.963970E-01_EB, &
    1.175949E+00_EB, 1.370990E+00_EB, 1.580028E+00_EB, 1.801371E+00_EB, 2.033228E+00_EB, 2.273775E+00_EB, &
    2.521222E+00_EB, 2.773848E+00_EB, 3.030039E+00_EB, 3.288302E+00_EB, 3.547276E+00_EB,  &
    1.433914E+00_EB, 1.266683E+00_EB, 1.138412E+00_EB, 1.035775E+00_EB, 9.511150E-01_EB, 8.796924E-01_EB, &
    8.183951E-01_EB, 7.650744E-01_EB, 7.181857E-01_EB, 6.765815E-01_EB, 6.393883E-01_EB, 6.059223E-01_EB, &
    5.756412E-01_EB, 5.481068E-01_EB, 5.229564E-01_EB, 4.998961E-01_EB, 4.786749E-01_EB, 4.590826E-01_EB, &
    4.409405E-01_EB, 4.240955E-01_EB, 4.084127E-01_EB, 3.937797E-01_EB, 3.800947E-01_EB,  &
    7.982731E-01_EB, 6.921371E-01_EB, 6.102233E-01_EB, 5.447290E-01_EB, 4.910140E-01_EB, 4.461112E-01_EB, &
    4.080120E-01_EB, 3.752913E-01_EB, 3.469076E-01_EB, 3.220760E-01_EB, 3.001886E-01_EB, 2.807714E-01_EB, &
    2.634443E-01_EB, 2.479029E-01_EB, 2.338953E-01_EB, 2.212166E-01_EB, 2.096943E-01_EB, 1.991872E-01_EB, &
    1.895690E-01_EB, 1.807401E-01_EB, 1.726116E-01_EB, 1.651065E-01_EB, 1.581598E-01_EB,  &
    1.286717E-02_EB, 1.469682E-02_EB, 1.762203E-02_EB, 2.284794E-02_EB, 3.191942E-02_EB, 4.657260E-02_EB, &
    6.851330E-02_EB, 9.924421E-02_EB, 1.399155E-01_EB, 1.912727E-01_EB, 2.536483E-01_EB, 3.269857E-01_EB, &
    4.109064E-01_EB, 5.047776E-01_EB, 6.077700E-01_EB, 7.189365E-01_EB, 8.372551E-01_EB, 9.616748E-01_EB, &
    1.091156E+00_EB, 1.224697E+00_EB, 1.361339E+00_EB, 1.500194E+00_EB, 1.640448E+00_EB,  &
    6.157895E-03_EB, 1.230974E-02_EB, 2.114024E-02_EB, 3.245860E-02_EB, 4.586804E-02_EB, 6.088403E-02_EB, &
    7.702163E-02_EB, 9.383227E-02_EB, 1.109358E-01_EB, 1.280177E-01_EB, 1.448305E-01_EB, 1.611855E-01_EB, &
    1.769419E-01_EB, 1.920025E-01_EB, 2.063046E-01_EB, 2.198097E-01_EB, 2.325011E-01_EB, 2.443776E-01_EB, &
    2.554522E-01_EB, 2.657423E-01_EB, 2.752770E-01_EB, 2.840858E-01_EB, 2.922023E-01_EB,  &
    4.718063E-01_EB, 5.060694E-01_EB, 5.352013E-01_EB, 5.592650E-01_EB, 5.784582E-01_EB, 5.931329E-01_EB, &
    6.037368E-01_EB, 6.107592E-01_EB, 6.146914E-01_EB, 6.159966E-01_EB, 6.150922E-01_EB, 6.123560E-01_EB, &
    6.081126E-01_EB, 6.026439E-01_EB, 5.961898E-01_EB, 5.889540E-01_EB, 5.811085E-01_EB, 5.727969E-01_EB, &
    5.641395E-01_EB, 5.552368E-01_EB, 5.461737E-01_EB, 5.370157E-01_EB, 5.278213E-01_EB,  &
    3.249465E+00_EB, 2.623698E+00_EB, 2.174113E+00_EB, 1.836681E+00_EB, 1.575103E+00_EB, 1.367239E+00_EB, &
    1.198789E+00_EB, 1.060085E+00_EB, 9.443468E-01_EB, 8.466759E-01_EB, 7.634443E-01_EB, 6.919088E-01_EB, &
    6.299582E-01_EB, 5.759445E-01_EB, 5.285628E-01_EB, 4.867671E-01_EB, 4.497122E-01_EB, 4.167074E-01_EB, &
    3.871836E-01_EB, 3.606688E-01_EB, 3.367684E-01_EB, 3.151516E-01_EB, 2.955362E-01_EB/),(/23,8/))

SD2_MMA(1:23,9:16) = RESHAPE((/ &  ! 950-1125 cm-1
    7.284805E-01_EB, 6.295084E-01_EB, 5.557612E-01_EB, 4.983242E-01_EB, 4.521021E-01_EB, 4.139657E-01_EB, &
    3.818837E-01_EB, 3.544713E-01_EB, 3.307499E-01_EB, 3.100011E-01_EB, 2.916913E-01_EB, 2.754062E-01_EB, &
    2.608234E-01_EB, 2.476896E-01_EB, 2.357954E-01_EB, 2.249743E-01_EB, 2.150857E-01_EB, 2.060142E-01_EB, &
    1.976647E-01_EB, 1.899522E-01_EB, 1.828084E-01_EB, 1.761713E-01_EB, 1.699903E-01_EB,  &
    4.585157E-01_EB, 4.381171E-01_EB, 4.193556E-01_EB, 4.017948E-01_EB, 3.852291E-01_EB, 3.695629E-01_EB, &
    3.547397E-01_EB, 3.407205E-01_EB, 3.274710E-01_EB, 3.149555E-01_EB, 3.031402E-01_EB, 2.919859E-01_EB, &
    2.814554E-01_EB, 2.715118E-01_EB, 2.621204E-01_EB, 2.532430E-01_EB, 2.448498E-01_EB, 2.369050E-01_EB, &
    2.293817E-01_EB, 2.222517E-01_EB, 2.154881E-01_EB, 2.090667E-01_EB, 2.029670E-01_EB,  &
    1.416841E+00_EB, 1.194161E+00_EB, 1.031840E+00_EB, 9.076867E-01_EB, 8.092920E-01_EB, 7.291881E-01_EB, &
    6.625971E-01_EB, 6.063128E-01_EB, 5.580940E-01_EB, 5.163216E-01_EB, 4.797892E-01_EB, 4.475805E-01_EB, &
    4.189839E-01_EB, 3.934354E-01_EB, 3.704844E-01_EB, 3.497664E-01_EB, 3.309802E-01_EB, 3.138769E-01_EB, &
    2.982469E-01_EB, 2.839173E-01_EB, 2.707368E-01_EB, 2.585782E-01_EB, 2.473330E-01_EB,  &
    7.478840E-01_EB, 5.631802E-01_EB, 4.390387E-01_EB, 3.513186E-01_EB, 2.869487E-01_EB, 2.383047E-01_EB, &
    2.006680E-01_EB, 1.709760E-01_EB, 1.471694E-01_EB, 1.278136E-01_EB, 1.118862E-01_EB, 9.864048E-02_EB, &
    8.752013E-02_EB, 7.810634E-02_EB, 7.007487E-02_EB, 6.317642E-02_EB, 5.721233E-02_EB, 5.202592E-02_EB, &
    4.749435E-02_EB, 4.351230E-02_EB, 3.999819E-02_EB, 3.688344E-02_EB, 3.411123E-02_EB,  &
    2.193302E-03_EB, 4.196652E-03_EB, 7.819119E-03_EB, 1.370559E-02_EB, 2.225589E-02_EB, 3.355967E-02_EB, &
    4.742675E-02_EB, 6.348296E-02_EB, 8.126100E-02_EB, 1.002743E-01_EB, 1.200662E-01_EB, 1.402343E-01_EB, &
    1.604448E-01_EB, 1.804275E-01_EB, 1.999763E-01_EB, 2.189355E-01_EB, 2.371970E-01_EB, 2.546853E-01_EB, &
    2.713578E-01_EB, 2.871929E-01_EB, 3.021850E-01_EB, 3.163426E-01_EB, 3.296819E-01_EB,  &
    5.039379E-02_EB, 7.726045E-02_EB, 1.060089E-01_EB, 1.349347E-01_EB, 1.628608E-01_EB, 1.890516E-01_EB, &
    2.130964E-01_EB, 2.348188E-01_EB, 2.541910E-01_EB, 2.712860E-01_EB, 2.862309E-01_EB, 2.991833E-01_EB, &
    3.103135E-01_EB, 3.197969E-01_EB, 3.277958E-01_EB, 3.344680E-01_EB, 3.399591E-01_EB, 3.443997E-01_EB, &
    3.479107E-01_EB, 3.505991E-01_EB, 3.525635E-01_EB, 3.538874E-01_EB, 3.546480E-01_EB,  &
    1.801396E-01_EB, 2.721529E-01_EB, 3.698337E-01_EB, 4.677985E-01_EB, 5.623828E-01_EB, 6.513087E-01_EB, &
    7.333140E-01_EB, 8.078454E-01_EB, 8.748175E-01_EB, 9.344455E-01_EB, 9.871178E-01_EB, 1.033317E+00_EB, &
    1.073572E+00_EB, 1.108414E+00_EB, 1.138357E+00_EB, 1.163892E+00_EB, 1.185474E+00_EB, 1.203522E+00_EB, &
    1.218415E+00_EB, 1.230498E+00_EB, 1.240081E+00_EB, 1.247441E+00_EB, 1.252831E+00_EB,  &
    1.612314E+00_EB, 2.147516E+00_EB, 2.670779E+00_EB, 3.164487E+00_EB, 3.618659E+00_EB, 4.028670E+00_EB, &
    4.393461E+00_EB, 4.714209E+00_EB, 4.993394E+00_EB, 5.234182E+00_EB, 5.440019E+00_EB, 5.614370E+00_EB, &
    5.760573E+00_EB, 5.881755E+00_EB, 5.980788E+00_EB, 6.060277E+00_EB, 6.122565E+00_EB, 6.169742E+00_EB, &
    6.203668E+00_EB, 6.225990E+00_EB, 6.238163E+00_EB, 6.241472E+00_EB, 6.237051E+00_EB/),(/23,8/))

SD2_MMA(1:23,17:24) = RESHAPE((/ &  ! 1150-1325 cm-1
    1.884164E+01_EB, 1.610097E+01_EB, 1.406049E+01_EB, 1.247199E+01_EB, 1.119373E+01_EB, 1.013904E+01_EB, &
    9.251836E+00_EB, 8.494076E+00_EB, 7.838877E+00_EB, 7.266616E+00_EB, 6.762545E+00_EB, 6.315325E+00_EB, &
    5.916057E+00_EB, 5.557640E+00_EB, 5.234321E+00_EB, 4.941387E+00_EB, 4.674932E+00_EB, 4.431688E+00_EB, &
    4.208898E+00_EB, 4.004222E+00_EB, 3.815654E+00_EB, 3.641471E+00_EB, 3.480180E+00_EB,  &
    1.130865E+01_EB, 9.595553E+00_EB, 8.287192E+00_EB, 7.250012E+00_EB, 6.405758E+00_EB, 5.705003E+00_EB, &
    5.114619E+00_EB, 4.611325E+00_EB, 4.178133E+00_EB, 3.802265E+00_EB, 3.473877E+00_EB, 3.185235E+00_EB, &
    2.930167E+00_EB, 2.703677E+00_EB, 2.501681E+00_EB, 2.320806E+00_EB, 2.158243E+00_EB, 2.011635E+00_EB, &
    1.878990E+00_EB, 1.758622E+00_EB, 1.649083E+00_EB, 1.549137E+00_EB, 1.457715E+00_EB,  &
    7.487490E+00_EB, 5.801604E+00_EB, 4.649409E+00_EB, 3.820657E+00_EB, 3.201009E+00_EB, 2.723538E+00_EB, &
    2.346693E+00_EB, 2.043413E+00_EB, 1.795363E+00_EB, 1.589703E+00_EB, 1.417188E+00_EB, 1.271006E+00_EB, &
    1.146033E+00_EB, 1.038351E+00_EB, 9.449086E-01_EB, 8.633103E-01_EB, 7.916426E-01_EB, 7.283673E-01_EB, &
    6.722323E-01_EB, 6.222099E-01_EB, 5.774522E-01_EB, 5.372524E-01_EB, 5.010175E-01_EB,  &
    9.789231E-01_EB, 9.295129E-01_EB, 8.835128E-01_EB, 8.406546E-01_EB, 8.006386E-01_EB, 7.632108E-01_EB, &
    7.281650E-01_EB, 6.953268E-01_EB, 6.645441E-01_EB, 6.356760E-01_EB, 6.085946E-01_EB, 5.831770E-01_EB, &
    5.593085E-01_EB, 5.368818E-01_EB, 5.157934E-01_EB, 4.959505E-01_EB, 4.772638E-01_EB, 4.596531E-01_EB, &
    4.430407E-01_EB, 4.273553E-01_EB, 4.125331E-01_EB, 3.985137E-01_EB, 3.852415E-01_EB,  &
    3.294228E-01_EB, 3.896080E-01_EB, 4.454235E-01_EB, 4.959355E-01_EB, 5.408165E-01_EB, 5.801207E-01_EB, &
    6.141261E-01_EB, 6.432285E-01_EB, 6.678744E-01_EB, 6.885221E-01_EB, 7.056151E-01_EB, 7.195669E-01_EB, &
    7.307567E-01_EB, 7.395282E-01_EB, 7.461867E-01_EB, 7.510042E-01_EB, 7.542201E-01_EB, 7.560464E-01_EB, &
    7.566653E-01_EB, 7.562429E-01_EB, 7.549174E-01_EB, 7.528138E-01_EB, 7.500416E-01_EB,  &
    2.577277E+00_EB, 2.629910E+00_EB, 2.666932E+00_EB, 2.688311E+00_EB, 2.695713E+00_EB, 2.691278E+00_EB, &
    2.677143E+00_EB, 2.655236E+00_EB, 2.627222E+00_EB, 2.594499E+00_EB, 2.558227E+00_EB, 2.519350E+00_EB, &
    2.478638E+00_EB, 2.436707E+00_EB, 2.394056E+00_EB, 2.351079E+00_EB, 2.308088E+00_EB, 2.265331E+00_EB, &
    2.222999E+00_EB, 2.181242E+00_EB, 2.140173E+00_EB, 2.099876E+00_EB, 2.060415E+00_EB,  &
    7.100483E+00_EB, 6.127430E+00_EB, 5.363729E+00_EB, 4.748424E+00_EB, 4.242049E+00_EB, 3.818126E+00_EB, &
    3.458265E+00_EB, 3.149277E+00_EB, 2.881428E+00_EB, 2.647360E+00_EB, 2.441382E+00_EB, 2.259018E+00_EB, &
    2.096690E+00_EB, 1.951500E+00_EB, 1.821069E+00_EB, 1.703437E+00_EB, 1.596959E+00_EB, 1.500253E+00_EB, &
    1.412149E+00_EB, 1.331651E+00_EB, 1.257901E+00_EB, 1.190163E+00_EB, 1.127798E+00_EB,  &
    5.779449E+00_EB, 4.301252E+00_EB, 3.325323E+00_EB, 2.645172E+00_EB, 2.151166E+00_EB, 1.780562E+00_EB, &
    1.495240E+00_EB, 1.270881E+00_EB, 1.091332E+00_EB, 9.454943E-01_EB, 8.255233E-01_EB, 7.257350E-01_EB, &
    6.419272E-01_EB, 5.709336E-01_EB, 5.103316E-01_EB, 4.582428E-01_EB, 4.131876E-01_EB, 3.739956E-01_EB, &
    3.397215E-01_EB, 3.096070E-01_EB, 2.830246E-01_EB, 2.594661E-01_EB, 2.385058E-01_EB/),(/23,8/))

SD2_MMA(1:23,25:32) = RESHAPE((/ &  ! 1350-1525 cm-1
    5.683804E-01_EB, 5.324585E-01_EB, 5.009940E-01_EB, 4.729980E-01_EB, 4.477695E-01_EB, 4.248073E-01_EB, &
    4.037597E-01_EB, 3.843620E-01_EB, 3.664133E-01_EB, 3.497540E-01_EB, 3.342522E-01_EB, 3.197992E-01_EB, &
    3.062984E-01_EB, 2.936691E-01_EB, 2.818369E-01_EB, 2.707375E-01_EB, 2.603109E-01_EB, 2.505060E-01_EB, &
    2.412740E-01_EB, 2.325709E-01_EB, 2.243589E-01_EB, 2.166000E-01_EB, 2.092639E-01_EB,  &
    9.513542E-01_EB, 8.415325E-01_EB, 7.558512E-01_EB, 6.864731E-01_EB, 6.286952E-01_EB, 5.795317E-01_EB, &
    5.370002E-01_EB, 4.997294E-01_EB, 4.667321E-01_EB, 4.372780E-01_EB, 4.108070E-01_EB, 3.868805E-01_EB, &
    3.651530E-01_EB, 3.453364E-01_EB, 3.271982E-01_EB, 3.105409E-01_EB, 2.951988E-01_EB, 2.810300E-01_EB, &
    2.679131E-01_EB, 2.557412E-01_EB, 2.444247E-01_EB, 2.338805E-01_EB, 2.240395E-01_EB,  &
    1.116372E+00_EB, 1.037754E+00_EB, 9.742196E-01_EB, 9.208877E-01_EB, 8.748003E-01_EB, 8.340920E-01_EB, &
    7.975436E-01_EB, 7.643246E-01_EB, 7.338582E-01_EB, 7.057211E-01_EB, 6.795966E-01_EB, 6.552407E-01_EB, &
    6.324584E-01_EB, 6.110862E-01_EB, 5.909935E-01_EB, 5.720632E-01_EB, 5.541989E-01_EB, 5.373138E-01_EB, &
    5.213309E-01_EB, 5.061813E-01_EB, 4.918059E-01_EB, 4.781482E-01_EB, 4.651571E-01_EB,  &
    3.273312E+00_EB, 2.728178E+00_EB, 2.337179E+00_EB, 2.043288E+00_EB, 1.814236E+00_EB, 1.630543E+00_EB, &
    1.479819E+00_EB, 1.353826E+00_EB, 1.246883E+00_EB, 1.154942E+00_EB, 1.075042E+00_EB, 1.004963E+00_EB, &
    9.430036E-01_EB, 8.878429E-01_EB, 8.384283E-01_EB, 7.939204E-01_EB, 7.536341E-01_EB, 7.170062E-01_EB, &
    6.835711E-01_EB, 6.529370E-01_EB, 6.247751E-01_EB, 5.988042E-01_EB, 5.747847E-01_EB,  &
    3.376325E+00_EB, 2.617646E+00_EB, 2.097702E+00_EB, 1.723128E+00_EB, 1.442747E+00_EB, 1.226471E+00_EB, &
    1.055586E+00_EB, 9.178996E-01_EB, 8.051513E-01_EB, 7.115618E-01_EB, 6.329732E-01_EB, 5.663182E-01_EB, &
    5.092931E-01_EB, 4.601281E-01_EB, 4.174507E-01_EB, 3.801774E-01_EB, 3.474425E-01_EB, 3.185487E-01_EB, &
    2.929278E-01_EB, 2.701123E-01_EB, 2.497151E-01_EB, 2.314145E-01_EB, 2.149387E-01_EB,  &
    6.024239E-01_EB, 5.422734E-01_EB, 4.963008E-01_EB, 4.596326E-01_EB, 4.293937E-01_EB, 4.037909E-01_EB, &
    3.816595E-01_EB, 3.622124E-01_EB, 3.448982E-01_EB, 3.293225E-01_EB, 3.151924E-01_EB, 3.022813E-01_EB, &
    2.904192E-01_EB, 2.794667E-01_EB, 2.693132E-01_EB, 2.598660E-01_EB, 2.510485E-01_EB, 2.427962E-01_EB, &
    2.350552E-01_EB, 2.277763E-01_EB, 2.209197E-01_EB, 2.144463E-01_EB, 2.083276E-01_EB,  &
    1.479397E-01_EB, 1.507368E-01_EB, 1.525819E-01_EB, 1.536698E-01_EB, 1.541358E-01_EB, 1.540882E-01_EB, &
    1.536146E-01_EB, 1.527867E-01_EB, 1.516703E-01_EB, 1.503195E-01_EB, 1.487798E-01_EB, 1.470909E-01_EB, &
    1.452868E-01_EB, 1.433961E-01_EB, 1.414402E-01_EB, 1.394396E-01_EB, 1.374103E-01_EB, 1.353680E-01_EB, &
    1.333210E-01_EB, 1.312810E-01_EB, 1.292530E-01_EB, 1.272443E-01_EB, 1.252593E-01_EB,  &
    1.190071E-01_EB, 1.189246E-01_EB, 1.185859E-01_EB, 1.180199E-01_EB, 1.172507E-01_EB, 1.163020E-01_EB, &
    1.151983E-01_EB, 1.139632E-01_EB, 1.126212E-01_EB, 1.111955E-01_EB, 1.097021E-01_EB, 1.081595E-01_EB, &
    1.065828E-01_EB, 1.049858E-01_EB, 1.033759E-01_EB, 1.017635E-01_EB, 1.001582E-01_EB, 9.856306E-02_EB, &
    9.698340E-02_EB, 9.542390E-02_EB, 9.388753E-02_EB, 9.237608E-02_EB, 9.089244E-02_EB/),(/23,8/))

SD2_MMA(1:23,33:40) = RESHAPE((/ &  ! 1550-1725 cm-1
    1.071482E-01_EB, 1.106118E-01_EB, 1.124954E-01_EB, 1.132423E-01_EB, 1.131767E-01_EB, 1.125339E-01_EB, &
    1.114762E-01_EB, 1.101249E-01_EB, 1.085677E-01_EB, 1.068704E-01_EB, 1.050784E-01_EB, 1.032293E-01_EB, &
    1.013509E-01_EB, 9.946241E-02_EB, 9.757923E-02_EB, 9.571387E-02_EB, 9.387499E-02_EB, 9.206748E-02_EB, &
    9.029688E-02_EB, 8.856563E-02_EB, 8.687714E-02_EB, 8.523074E-02_EB, 8.362699E-02_EB,  &
    1.791350E-01_EB, 1.781579E-01_EB, 1.771148E-01_EB, 1.759549E-01_EB, 1.746470E-01_EB, 1.731814E-01_EB, &
    1.715609E-01_EB, 1.697959E-01_EB, 1.679061E-01_EB, 1.659094E-01_EB, 1.638262E-01_EB, 1.616761E-01_EB, &
    1.594742E-01_EB, 1.572384E-01_EB, 1.549818E-01_EB, 1.527150E-01_EB, 1.504497E-01_EB, 1.481937E-01_EB, &
    1.459536E-01_EB, 1.437370E-01_EB, 1.415456E-01_EB, 1.393859E-01_EB, 1.372605E-01_EB,  &
    3.881465E-01_EB, 3.949685E-01_EB, 3.997660E-01_EB, 4.029267E-01_EB, 4.047169E-01_EB, 4.053405E-01_EB, &
    4.049656E-01_EB, 4.037407E-01_EB, 4.017975E-01_EB, 3.992479E-01_EB, 3.961949E-01_EB, 3.927267E-01_EB, &
    3.889211E-01_EB, 3.848451E-01_EB, 3.805562E-01_EB, 3.761031E-01_EB, 3.715284E-01_EB, 3.668673E-01_EB, &
    3.621498E-01_EB, 3.574024E-01_EB, 3.526434E-01_EB, 3.478935E-01_EB, 3.431675E-01_EB,  &
    1.672640E+00_EB, 1.345702E+00_EB, 1.116168E+00_EB, 9.472291E-01_EB, 8.182603E-01_EB, 7.169090E-01_EB, &
    6.353739E-01_EB, 5.685020E-01_EB, 5.127719E-01_EB, 4.656940E-01_EB, 4.254665E-01_EB, 3.907467E-01_EB, &
    3.605209E-01_EB, 3.340073E-01_EB, 3.105906E-01_EB, 2.897868E-01_EB, 2.712013E-01_EB, 2.545174E-01_EB, &
    2.394731E-01_EB, 2.258505E-01_EB, 2.134696E-01_EB, 2.021779E-01_EB, 1.918458E-01_EB,  &
    6.948718E-01_EB, 5.621678E-01_EB, 4.720523E-01_EB, 4.074800E-01_EB, 3.592174E-01_EB, 3.218978E-01_EB, &
    2.922227E-01_EB, 2.680720E-01_EB, 2.480294E-01_EB, 2.311148E-01_EB, 2.166373E-01_EB, 2.040909E-01_EB, &
    1.931020E-01_EB, 1.833843E-01_EB, 1.747222E-01_EB, 1.669434E-01_EB, 1.599121E-01_EB, 1.535199E-01_EB, &
    1.476771E-01_EB, 1.423154E-01_EB, 1.373690E-01_EB, 1.327923E-01_EB, 1.285411E-01_EB,  &
    4.229856E-01_EB, 4.260594E-01_EB, 4.278430E-01_EB, 4.285515E-01_EB, 4.283233E-01_EB, 4.272674E-01_EB, &
    4.254826E-01_EB, 4.230616E-01_EB, 4.200917E-01_EB, 4.166526E-01_EB, 4.128223E-01_EB, 4.086681E-01_EB, &
    4.042497E-01_EB, 3.996231E-01_EB, 3.948354E-01_EB, 3.899268E-01_EB, 3.849338E-01_EB, 3.798855E-01_EB, &
    3.748073E-01_EB, 3.697233E-01_EB, 3.646499E-01_EB, 3.596016E-01_EB, 3.545935E-01_EB,  &
    1.276432E+00_EB, 1.471014E+00_EB, 1.635605E+00_EB, 1.774228E+00_EB, 1.890482E+00_EB, 1.987460E+00_EB, &
    2.067792E+00_EB, 2.133732E+00_EB, 2.187217E+00_EB, 2.229919E+00_EB, 2.263292E+00_EB, 2.288599E+00_EB, &
    2.306933E+00_EB, 2.319250E+00_EB, 2.326379E+00_EB, 2.329039E+00_EB, 2.327857E+00_EB, 2.323377E+00_EB, &
    2.316073E+00_EB, 2.306356E+00_EB, 2.294583E+00_EB, 2.281065E+00_EB, 2.266073E+00_EB,  &
    1.480365E+01_EB, 1.219102E+01_EB, 1.034180E+01_EB, 8.968958E+00_EB, 7.911019E+00_EB, 7.071040E+00_EB, &
    6.387735E+00_EB, 5.820678E+00_EB, 5.342225E+00_EB, 4.932874E+00_EB, 4.578498E+00_EB, 4.268611E+00_EB, &
    3.995270E+00_EB, 3.752338E+00_EB, 3.535006E+00_EB, 3.339436E+00_EB, 3.162537E+00_EB, 3.001779E+00_EB, &
    2.855079E+00_EB, 2.720698E+00_EB, 2.597177E+00_EB, 2.483275E+00_EB, 2.377940E+00_EB/),(/23,8/))

SD2_MMA(1:23,41:48) = RESHAPE((/ &  ! 1750-1925 cm-1
    7.194352E+00_EB, 6.108145E+00_EB, 5.307516E+00_EB, 4.690351E+00_EB, 4.197953E+00_EB, 3.794357E+00_EB, &
    3.456386E+00_EB, 3.168464E+00_EB, 2.919737E+00_EB, 2.702405E+00_EB, 2.510704E+00_EB, 2.340270E+00_EB, &
    2.187725E+00_EB, 2.050408E+00_EB, 1.926182E+00_EB, 1.813312E+00_EB, 1.710361E+00_EB, 1.616137E+00_EB, &
    1.529632E+00_EB, 1.449990E+00_EB, 1.376477E+00_EB, 1.308464E+00_EB, 1.245399E+00_EB,  &
    3.321285E-01_EB, 3.739563E-01_EB, 4.087610E-01_EB, 4.370667E-01_EB, 4.596328E-01_EB, 4.772445E-01_EB, &
    4.906296E-01_EB, 5.004317E-01_EB, 5.072119E-01_EB, 5.114539E-01_EB, 5.135705E-01_EB, 5.139118E-01_EB, &
    5.127805E-01_EB, 5.104273E-01_EB, 5.070710E-01_EB, 5.028903E-01_EB, 4.980427E-01_EB, 4.926572E-01_EB, &
    4.868442E-01_EB, 4.806960E-01_EB, 4.742911E-01_EB, 4.676961E-01_EB, 4.609669E-01_EB,  &
    1.852642E-01_EB, 1.718763E-01_EB, 1.607701E-01_EB, 1.513312E-01_EB, 1.431381E-01_EB, 1.359058E-01_EB, &
    1.294262E-01_EB, 1.235592E-01_EB, 1.181953E-01_EB, 1.132573E-01_EB, 1.086839E-01_EB, 1.044288E-01_EB, &
    1.004548E-01_EB, 9.673120E-02_EB, 9.323487E-02_EB, 8.994275E-02_EB, 8.683809E-02_EB, 8.390492E-02_EB, &
    8.112959E-02_EB, 7.850050E-02_EB, 7.600714E-02_EB, 7.363986E-02_EB, 7.138900E-02_EB,  &
    8.712487E-03_EB, 9.055414E-03_EB, 9.321856E-03_EB, 9.527903E-03_EB, 9.687670E-03_EB, 9.807490E-03_EB, &
    9.893895E-03_EB, 9.952344E-03_EB, 9.986048E-03_EB, 9.998023E-03_EB, 9.993729E-03_EB, 9.971896E-03_EB, &
    9.938079E-03_EB, 9.892966E-03_EB, 9.837815E-03_EB, 9.774575E-03_EB, 9.704130E-03_EB, 9.628520E-03_EB, &
    9.547458E-03_EB, 9.462492E-03_EB, 9.375194E-03_EB, 9.284290E-03_EB, 9.192702E-03_EB,  &
    8.012408E-02_EB, 7.702656E-02_EB, 7.451003E-02_EB, 7.237284E-02_EB, 7.049222E-02_EB, 6.878977E-02_EB, &
    6.721527E-02_EB, 6.573598E-02_EB, 6.432747E-02_EB, 6.297777E-02_EB, 6.167617E-02_EB, 6.041428E-02_EB, &
    5.919233E-02_EB, 5.800358E-02_EB, 5.684832E-02_EB, 5.572530E-02_EB, 5.463131E-02_EB, 5.356624E-02_EB, &
    5.253146E-02_EB, 5.152492E-02_EB, 5.054527E-02_EB, 4.959357E-02_EB, 4.866818E-02_EB,  &
    2.097905E-01_EB, 1.654856E-01_EB, 1.346542E-01_EB, 1.121809E-01_EB, 9.519814E-02_EB, 8.198461E-02_EB, &
    7.145622E-02_EB, 6.290142E-02_EB, 5.583508E-02_EB, 4.991753E-02_EB, 4.490405E-02_EB, 4.061273E-02_EB, &
    3.690619E-02_EB, 3.368104E-02_EB, 3.085511E-02_EB, 2.836385E-02_EB, 2.615560E-02_EB, 2.418992E-02_EB, &
    2.243165E-02_EB, 2.085160E-02_EB, 1.942730E-02_EB, 1.813859E-02_EB, 1.696989E-02_EB,  &
    1.935450E-02_EB, 1.711233E-02_EB, 1.544530E-02_EB, 1.415599E-02_EB, 1.312469E-02_EB, 1.227728E-02_EB, &
    1.156701E-02_EB, 1.095813E-02_EB, 1.042903E-02_EB, 9.961480E-03_EB, 9.544182E-03_EB, 9.168290E-03_EB, &
    8.826586E-03_EB, 8.514584E-03_EB, 8.226351E-03_EB, 7.960714E-03_EB, 7.712217E-03_EB, 7.481350E-03_EB, &
    7.264705E-03_EB, 7.061987E-03_EB, 6.869985E-03_EB, 6.687822E-03_EB, 6.517833E-03_EB,  &
    3.527490E-03_EB, 3.690590E-03_EB, 3.818609E-03_EB, 3.919835E-03_EB, 3.998656E-03_EB, 4.060714E-03_EB, &
    4.108260E-03_EB, 4.143528E-03_EB, 4.168572E-03_EB, 4.184353E-03_EB, 4.192236E-03_EB, 4.192141E-03_EB, &
    4.186784E-03_EB, 4.176649E-03_EB, 4.161739E-03_EB, 4.142638E-03_EB, 4.119845E-03_EB, 4.095094E-03_EB, &
    4.068399E-03_EB, 4.038290E-03_EB, 4.007600E-03_EB, 3.974765E-03_EB, 3.940956E-03_EB/),(/23,8/))

SD2_MMA(1:23,49:49) = RESHAPE((/ &  ! 1950-1950 cm-1
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, &
    0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB, 0.000000E+00_EB/),(/23,1/))

!--------------------------------OLD METHANE DATA----------------------------------------
!  Initialize SD Array

! TEMP,K= 300     600      1000      1500      2000      2500       

SD(1:6,1:8) = RESHAPE ((/  & ! 50-225
.950E+00_EB, .103E+00_EB, .420E-01_EB, .114E-01_EB, .450E-02_EB, .300E-02_EB, &
.208E+01_EB, .365E+00_EB, .113E+00_EB, .375E-01_EB, .195E-01_EB, .134E-01_EB, &
.368E+01_EB, .990E+00_EB, .300E+00_EB, .104E+00_EB, .577E-01_EB, .365E-01_EB,  & 
.650E+01_EB, .201E+01_EB, .650E+00_EB, .214E+00_EB, .128E+00_EB, .845E-01_EB, &   
.825E+01_EB, .325E+01_EB, .121E+01_EB, .415E+00_EB, .260E+00_EB, .168E+00_EB, &  
.870E+01_EB, .452E+01_EB, .189E+01_EB, .765E+00_EB, .450E+00_EB, .289E+00_EB, &  
.810E+01_EB, .540E+01_EB, .261E+01_EB, .126E+01_EB, .695E+00_EB, .460E+00_EB, &  
.682E+01_EB, .600E+01_EB, .337E+01_EB, .179E+01_EB, .101E+01_EB, .679E+00_EB/), &
(/6,8/))
SD(1:6,9:16) = RESHAPE ((/ &  ! 250-425
.493E+01_EB, .622E+01_EB, .407E+01_EB, .230E+01_EB, .135E+01_EB, .935E+00_EB,  &   
.316E+01_EB, .592E+01_EB, .456E+01_EB, .281E+01_EB, .172E+01_EB, .122E+01_EB,  &  
.199E+01_EB, .528E+01_EB, .479E+01_EB, .328E+01_EB, .213E+01_EB, .149E+01_EB,  &   
.113E+01_EB, .450E+01_EB, .484E+01_EB, .361E+01_EB, .249E+01_EB, .179E+01_EB,  &   
.585E+00_EB, .370E+01_EB, .471E+01_EB, .383E+01_EB, .284E+01_EB, .208E+01_EB, &    
.293E+00_EB, .289E+01_EB, .443E+01_EB, .394E+01_EB, .312E+01_EB, .237E+01_EB, &    
.138E+00_EB, .205E+01_EB, .400E+01_EB, .396E+01_EB, .330E+01_EB, .260E+01_EB, &    
.620E-01_EB, .143E+01_EB, .347E+01_EB, .388E+01_EB, .341E+01_EB, .280E+01_EB/), &  
(/6,8/))
SD(1:6,17:24) = RESHAPE ((/ & ! 450-625
.255E-01_EB, .950E+00_EB, .292E+01_EB, .370E+01_EB, .345E+01_EB, .295E+01_EB,   &  
.940E-02_EB, .610E+00_EB, .236E+01_EB, .343E+01_EB, .342E+01_EB, .304E+01_EB,  &   
.340E-02_EB, .386E+00_EB, .188E+01_EB, .310E+01_EB, .334E+01_EB, .309E+01_EB,   &  
.105E-02_EB, .236E+00_EB, .145E+01_EB, .274E+01_EB, .319E+01_EB, .307E+01_EB,   &  
.350E-03_EB, .144E+00_EB, .110E+01_EB, .238E+01_EB, .300E+01_EB, .301E+01_EB,   & 
.126E-03_EB, .820E-01_EB, .818E+00_EB, .204E+01_EB, .276E+01_EB, .289E+01_EB,   &  
.430E-04_EB, .445E-01_EB, .598E+00_EB, .174E+01_EB, .248E+01_EB, .275E+01_EB,  &  
.150E-04_EB, .242E-01_EB, .427E+00_EB, .145E+01_EB, .222E+01_EB, .260E+01_EB/), & 
(/6,8/))
SD(1:6,25:32) = RESHAPE ((/  &! 650-825
.510E-05_EB, .127E-01_EB, .294E+00_EB, .118E+01_EB, .195E+01_EB, .241E+01_EB,   &  
.170E-05_EB, .630E-02_EB, .200E+00_EB, .950E+00_EB, .169E+01_EB, .221E+01_EB,   &  
.570E-06_EB, .300E-02_EB, .134E+00_EB, .748E+00_EB, .146E+01_EB, .200E+01_EB,   &  
.195E-06_EB, .140E-02_EB, .902E-01_EB, .580E+00_EB, .124E+01_EB, .178E+01_EB,   &  
.680E-07_EB, .620E-03_EB, .590E-01_EB, .443E+00_EB, .103E+01_EB, .156E+01_EB,   &  
.385E-07_EB, .275E-03_EB, .450E-01_EB, .330E+00_EB, .845E+00_EB, .136E+01_EB,   &  
.670E-07_EB, .113E-03_EB, .355E-01_EB, .242E+00_EB, .695E+00_EB, .117E+01_EB,   &  
.113E-06_EB, .500E-04_EB, .289E-01_EB, .174E+00_EB, .560E+00_EB, .100E+01_EB/),  & 
(/6,8/))
SD(1:6,33:40) = RESHAPE ((/ & ! 850-1025
.195E-06_EB, .230E-04_EB, .245E-01_EB, .123E+00_EB, .450E+00_EB, .855E+00_EB,   &  
.328E-06_EB, .103E-04_EB, .214E-01_EB, .100E+00_EB, .357E+00_EB, .718E+00_EB,  &   
.560E-06_EB, .460E-05_EB, .189E-01_EB, .830E-01_EB, .278E+00_EB, .595E+00_EB,  &   
.950E-06_EB, .205E-05_EB, .174E-01_EB, .730E-01_EB, .239E+00_EB, .492E+00_EB,  &   
.160E-05_EB, .140E-05_EB, .166E-01_EB, .665E-01_EB, .211E+00_EB, .405E+00_EB, &    
.275E-05_EB, .350E-05_EB, .165E-01_EB, .630E-01_EB, .195E+00_EB, .352E+00_EB, &    
.470E-05_EB, .850E-05_EB, .167E-01_EB, .620E-01_EB, .190E+00_EB, .312E+00_EB, &   
.810E-05_EB, .215E-04_EB, .175E-01_EB, .630E-01_EB, .191E+00_EB, .289E+00_EB/), &  
(/6,8/))
SD(1:6,41:48) = RESHAPE ((/ & ! 1050-1225
.136E-04_EB, .570E-04_EB, .188E-01_EB, .675E-01_EB, .194E+00_EB, .281E+00_EB,   &  
.235E-04_EB, .150E-03_EB, .208E-01_EB, .745E-01_EB, .202E+00_EB, .283E+00_EB,   &  
.400E-04_EB, .380E-03_EB, .233E-01_EB, .865E-01_EB, .223E+00_EB, .314E+00_EB,  &   
.680E-04_EB, .950E-03_EB, .268E-01_EB, .122E+00_EB, .260E+00_EB, .380E+00_EB,  &   
.120E-03_EB, .245E-02_EB, .343E-01_EB, .176E+00_EB, .328E+00_EB, .461E+00_EB,  &   
.200E-03_EB, .620E-02_EB, .638E-01_EB, .251E+00_EB, .411E+00_EB, .511E+00_EB,  &   
.365E-03_EB, .140E-01_EB, .107E+00_EB, .330E+00_EB, .458E+00_EB, .542E+00_EB,  &   
.680E-03_EB, .330E-01_EB, .166E+00_EB, .405E+00_EB, .487E+00_EB, .571E+00_EB/),&   
(/6,8/))
SD(1:6,49:56) = RESHAPE ((/ & ! 1250-1425
.130E-02_EB, .635E-01_EB, .244E+00_EB, .459E+00_EB, .535E+00_EB, .557E+00_EB,  &   
.250E-02_EB, .123E+00_EB, .341E+00_EB, .477E+00_EB, .502E+00_EB, .562E+00_EB,  &   
.500E-02_EB, .212E+00_EB, .407E+00_EB, .547E+00_EB, .531E+00_EB, .514E+00_EB,  &   
.103E-01_EB, .285E+00_EB, .489E+00_EB, .592E+00_EB, .497E+00_EB, .486E+00_EB,  &   
.219E-01_EB, .328E+00_EB, .491E+00_EB, .558E+00_EB, .489E+00_EB, .485E+00_EB,  &   
.485E-01_EB, .345E+00_EB, .505E+00_EB, .521E+00_EB, .477E+00_EB, .484E+00_EB,  &   
.114E+00_EB, .361E+00_EB, .538E+00_EB, .563E+00_EB, .503E+00_EB, .502E+00_EB,  &   
.249E+00_EB, .460E+00_EB, .621E+00_EB, .624E+00_EB, .538E+00_EB, .538E+00_EB/), &  
(/6,8/))
SD(1:6,57:64) = RESHAPE ((/ & ! 1450-1625
.397E+00_EB, .569E+00_EB, .749E+00_EB, .768E+00_EB, .581E+00_EB, .565E+00_EB,  &   
.418E+00_EB, .627E+00_EB, .824E+00_EB, .849E+00_EB, .640E+00_EB, .594E+00_EB,  &   
.108E+01_EB, .125E+01_EB, .113E+01_EB, .940E+00_EB, .807E+00_EB, .663E+00_EB,  &   
.165E+01_EB, .155E+01_EB, .118E+01_EB, .670E+00_EB, .562E+00_EB, .483E+00_EB,  &   
.142E+01_EB, .675E+00_EB, .557E+00_EB, .349E+00_EB, .276E+00_EB, .263E+00_EB,  &   
.451E+00_EB, .202E+00_EB, .132E+00_EB, .118E+00_EB, .134E+00_EB, .156E+00_EB,  &   
.603E-01_EB, .538E-01_EB, .863E-01_EB, .112E+00_EB, .120E+00_EB, .125E+00_EB,  &   
.501E+00_EB, .252E+00_EB, .118E+00_EB, .112E+00_EB, .131E+00_EB, .140E+00_EB/), &  
(/6,8/))
SD(1:6,65:72) = RESHAPE ((/ & ! 1650-1825
.730E+00_EB, .430E+00_EB, .237E+00_EB, .191E+00_EB, .171E+00_EB, .170E+00_EB,   &  
.149E+01_EB, .506E+00_EB, .294E+00_EB, .238E+00_EB, .210E+00_EB, .201E+00_EB,  &   
.100E+01_EB, .553E+00_EB, .434E+00_EB, .340E+00_EB, .260E+00_EB, .220E+00_EB,   &  
.802E+00_EB, .658E+00_EB, .528E+00_EB, .411E+00_EB, .300E+00_EB, .240E+00_EB,    & 
.580E+00_EB, .527E+00_EB, .460E+00_EB, .378E+00_EB, .322E+00_EB, .283E+00_EB,  &   
.330E+00_EB, .403E+00_EB, .430E+00_EB, .356E+00_EB, .318E+00_EB, .270E+00_EB,  &   
.250E+00_EB, .393E+00_EB, .405E+00_EB, .342E+00_EB, .301E+00_EB, .275E+00_EB,  &   
.147E+00_EB, .249E+00_EB, .313E+00_EB, .318E+00_EB, .291E+00_EB, .268E+00_EB/), &  
(/6,8/))
SD(1:6,73:80) = RESHAPE ((/ & ! 1850-2025
.910E-01_EB, .252E+00_EB, .298E+00_EB, .295E+00_EB, .269E+00_EB, .253E+00_EB,   &  
.580E-01_EB, .158E+00_EB, .214E+00_EB, .244E+00_EB, .244E+00_EB, .245E+00_EB,   &  
.370E-01_EB, .113E+00_EB, .184E+00_EB, .218E+00_EB, .214E+00_EB, .218E+00_EB,   &  
.244E-01_EB, .118E+00_EB, .156E+00_EB, .188E+00_EB, .195E+00_EB, .200E+00_EB,   & 
.162E-01_EB, .606E-01_EB, .976E-01_EB, .141E+00_EB, .166E+00_EB, .179E+00_EB,   & 
.112E-01_EB, .425E-01_EB, .903E-01_EB, .133E+00_EB, .148E+00_EB, .156E+00_EB,   & 
.780E-02_EB, .400E-01_EB, .765E-01_EB, .112E+00_EB, .129E+00_EB, .137E+00_EB,   & 
.540E-02_EB, .352E-01_EB, .647E-01_EB, .876E-01_EB, .110E+00_EB, .118E+00_EB/),  &
(/6,8/))
SD(1:6,81:88) = RESHAPE ((/ & ! 2050-2225
.380E-02_EB, .252E-01_EB, .507E-01_EB, .705E-01_EB, .888E-01_EB, .100E+00_EB,  &   
.260E-02_EB, .179E-01_EB, .377E-01_EB, .546E-01_EB, .724E-01_EB, .828E-01_EB,  &   
.180E-02_EB, .123E-01_EB, .294E-01_EB, .443E-01_EB, .608E-01_EB, .686E-01_EB,  &   
.127E-02_EB, .850E-02_EB, .212E-01_EB, .378E-01_EB, .579E-01_EB, .640E-01_EB,   &  
.880E-03_EB, .680E-02_EB, .152E-01_EB, .275E-01_EB, .449E-01_EB, .521E-01_EB,  &   
.620E-02_EB, .400E-02_EB, .107E-01_EB, .214E-01_EB, .374E-01_EB, .453E-01_EB,  &   
.480E-03_EB, .298E-02_EB, .931E-02_EB, .189E-01_EB, .329E-01_EB, .403E-01_EB,  &   
.405E-03_EB, .175E-02_EB, .696E-02_EB, .152E-01_EB, .295E-01_EB, .365E-01_EB/), &  
(/6,8/))
SD(1:6,89:96) = RESHAPE ((/ & ! 2250-2425 
.321E-03_EB, .120E-02_EB, .452E-02_EB, .101E-01_EB, .252E-01_EB, .331E-01_EB, &    
.229E-03_EB, .721E-03_EB, .364E-02_EB, .930E-02_EB, .225E-01_EB, .305E-01_EB,  &   
.195E-03_EB, .544E-03_EB, .318E-02_EB, .750E-02_EB, .202E-01_EB, .284E-01_EB, &    
.154E-03_EB, .375E-03_EB, .185E-02_EB, .603E-02_EB, .175E-01_EB, .269E-01_EB, &    
.101E-03_EB, .263E-03_EB, .119E-02_EB, .480E-02_EB, .156E-01_EB, .253E-01_EB, &    
.852E-04_EB, .185E-03_EB, .909E-03_EB, .360E-02_EB, .133E-01_EB, .241E-01_EB,  &   
.763E-04_EB, .137E-03_EB, .711E-03_EB, .316E-02_EB, .122E-01_EB, .237E-01_EB,  &   
.615E-04_EB, .126E-03_EB, .610E-03_EB, .257E-02_EB, .101E-01_EB, .218E-01_EB/), &  
(/6,8/))
SD(1:6,97:104) = RESHAPE ((/ & ! 2450-2625 
.480E-04_EB, .113E-03_EB, .518E-03_EB, .201E-02_EB, .920E-02_EB, .200E-01_EB, &    
.372E-04_EB, .106E-03_EB, .435E-03_EB, .168E-02_EB, .785E-02_EB, .183E-01_EB,  &   
.355E-04_EB, .101E-03_EB, .376E-03_EB, .168E-02_EB, .669E-02_EB, .166E-01_EB, &    
.358E-04_EB, .990E-04_EB, .366E-03_EB, .167E-02_EB, .651E-02_EB, .156E-01_EB,  &   
.389E-04_EB, .102E-03_EB, .376E-03_EB, .167E-02_EB, .641E-02_EB, .152E-01_EB,  &   
.422E-04_EB, .106E-03_EB, .373E-03_EB, .168E-02_EB, .656E-02_EB, .150E-01_EB,  &   
.521E-04_EB, .111E-03_EB, .371E-03_EB, .170E-02_EB, .673E-02_EB, .152E-01_EB,  &   
.646E-04_EB, .121E-03_EB, .384E-03_EB, .179E-02_EB, .798E-02_EB, .179E-01_EB/), &  
(/6,8/)) 
SD(1:6,105:112) = RESHAPE ((/ & ! 2650-2825 
.742E-04_EB, .129E-03_EB, .479E-03_EB, .201E-02_EB, .788E-02_EB, .175E-01_EB,   &  
.953E-04_EB, .165E-03_EB, .544E-03_EB, .249E-02_EB, .945E-02_EB, .204E-01_EB,   &  
.101E-03_EB, .190E-03_EB, .761E-03_EB, .324E-02_EB, .106E-01_EB, .231E-01_EB,   &  
.147E-03_EB, .272E-03_EB, .892E-03_EB, .441E-02_EB, .125E-01_EB, .257E-01_EB,    & 
.195E-03_EB, .326E-03_EB, .100E-02_EB, .499E-02_EB, .147E-01_EB, .295E-01_EB,    & 
.261E-03_EB, .421E-03_EB, .145E-02_EB, .568E-02_EB, .161E-01_EB, .306E-01_EB,   &  
.305E-03_EB, .515E-03_EB, .195E-02_EB, .754E-02_EB, .185E-01_EB, .363E-01_EB,   &  
.362E-03_EB, .645E-03_EB, .237E-02_EB, .830E-02_EB, .205E-01_EB, .373E-01_EB/), &  
(/6,8/))
SD(1:6,113:120) = RESHAPE ((/ & ! 2850-3025 
.507E-03_EB, .850E-03_EB, .274E-02_EB, .888E-02_EB, .234E-01_EB, .431E-01_EB,     &
.799E-03_EB, .118E-02_EB, .322E-02_EB, .110E-01_EB, .262E-01_EB, .451E-01_EB,     &
.935E-03_EB, .160E-02_EB, .386E-02_EB, .126E-01_EB, .292E-01_EB, .530E-01_EB,     &
.108E-02_EB, .231E-02_EB, .451E-02_EB, .140E-01_EB, .306E-01_EB, .536E-01_EB,     &
.192E-02_EB, .271E-02_EB, .563E-02_EB, .159E-01_EB, .357E-01_EB, .629E-01_EB,     &
.263E-02_EB, .300E-02_EB, .625E-02_EB, .179E-01_EB, .385E-01_EB, .666E-01_EB,    & 
.295E-02_EB, .330E-02_EB, .701E-02_EB, .203E-01_EB, .460E-01_EB, .782E-01_EB,   &  
.310E-02_EB, .370E-02_EB, .846E-02_EB, .220E-01_EB, .519E-01_EB, .889E-01_EB/),  & 
(/6,8/))
SD(1:6,121:128) = RESHAPE ((/ & ! 3050-3225
.340E-02_EB, .400E-02_EB, .969E-02_EB, .279E-01_EB, .662E-01_EB, .109E+00_EB,     &
.730E-02_EB, .450E-02_EB, .111E-01_EB, .272E-01_EB, .676E-01_EB, .109E+00_EB,     &
.900E-02_EB, .480E-02_EB, .137E-01_EB, .372E-01_EB, .864E-01_EB, .133E+00_EB,     &
.100E-02_EB, .510E-02_EB, .162E-01_EB, .471E-01_EB, .100E+00_EB, .142E+00_EB,     &
.640E-03_EB, .550E-02_EB, .205E-01_EB, .530E-01_EB, .122E+00_EB, .168E+00_EB,     &
.160E-02_EB, .600E-02_EB, .247E-01_EB, .633E-01_EB, .135E+00_EB, .177E+00_EB,    & 
.330E-02_EB, .700E-02_EB, .283E-01_EB, .770E-01_EB, .153E+00_EB, .185E+00_EB,   &  
.410E-02_EB, .860E-02_EB, .376E-01_EB, .914E-01_EB, .166E+00_EB, .206E+00_EB/),&   
(/6,8/))
SD(1:6,129:136) = RESHAPE ((/ & ! 3250-3425
.410E-02_EB, .103E-01_EB, .514E-01_EB, .117E+00_EB, .194E+00_EB, .228E+00_EB,     &
.290E-02_EB, .129E-01_EB, .664E-01_EB, .147E+00_EB, .220E+00_EB, .254E+00_EB,     &
.220E-02_EB, .161E-01_EB, .834E-01_EB, .171E+00_EB, .237E+00_EB, .263E+00_EB,     &
.220E-02_EB, .212E-01_EB, .103E+00_EB, .201E+00_EB, .268E+00_EB, .283E+00_EB,     &
.250E-02_EB, .285E-01_EB, .135E+00_EB, .240E+00_EB, .295E+00_EB, .295E+00_EB,     &
.310E-02_EB, .385E-01_EB, .169E+00_EB, .272E+00_EB, .312E+00_EB, .301E+00_EB,     &
.420E-02_EB, .540E-01_EB, .214E+00_EB, .309E+00_EB, .329E+00_EB, .307E+00_EB,     &
.600E-02_EB, .770E-01_EB, .267E+00_EB, .343E+00_EB, .332E+00_EB, .314E+00_EB/),  & 
(/6,8/))
SD(1:6,137:144) = RESHAPE ((/ & ! 3450-3625
.940E-02_EB, .117E+00_EB, .333E+00_EB, .372E+00_EB, .344E+00_EB, .303E+00_EB,     &
.165E-01_EB, .173E+00_EB, .365E+00_EB, .385E+00_EB, .353E+00_EB, .300E+00_EB,     &
.360E-01_EB, .258E+00_EB, .438E+00_EB, .393E+00_EB, .315E+00_EB, .288E+00_EB,     &
.720E-01_EB, .375E+00_EB, .510E+00_EB, .409E+00_EB, .294E+00_EB, .271E+00_EB,     &
.133E+00_EB, .401E+00_EB, .499E+00_EB, .390E+00_EB, .281E+00_EB, .257E+00_EB,     &
.215E+00_EB, .500E+00_EB, .443E+00_EB, .341E+00_EB, .254E+00_EB, .230E+00_EB,     &
.318E+00_EB, .450E+00_EB, .346E+00_EB, .286E+00_EB, .245E+00_EB, .219E+00_EB,     &
.442E+00_EB, .400E+00_EB, .354E+00_EB, .279E+00_EB, .233E+00_EB, .216E+00_EB/),   &
(/6,8/))
SD(1:6,145:152) = RESHAPE ((/ & ! 3650-3825
.473E+00_EB, .405E+00_EB, .347E+00_EB, .281E+00_EB, .238E+00_EB, .219E+00_EB,     &
.568E+00_EB, .501E+00_EB, .423E+00_EB, .315E+00_EB, .243E+00_EB, .218E+00_EB,     &
.690E+00_EB, .708E+00_EB, .673E+00_EB, .432E+00_EB, .268E+00_EB, .189E+00_EB,     &
.617E+00_EB, .831E+00_EB, .566E+00_EB, .320E+00_EB, .194E+00_EB, .123E+00_EB,     &
.181E+01_EB, .520E+00_EB, .200E+00_EB, .131E+00_EB, .124E+00_EB, .107E+00_EB,     &
.136E+00_EB, .124E+00_EB, .120E+00_EB, .119E+00_EB, .115E+00_EB, .115E+00_EB,     &
.455E+00_EB, .298E+00_EB, .167E+00_EB, .129E+00_EB, .123E+00_EB, .112E+00_EB,     &
.760E+00_EB, .503E+00_EB, .242E+00_EB, .154E+00_EB, .129E+00_EB, .127E+00_EB/),   &
(/6,8/))
SD(1:6,153:160) = RESHAPE ((/ & ! 3850-4025
.836E+00_EB, .584E+00_EB, .277E+00_EB, .184E+00_EB, .161E+00_EB, .145E+00_EB,     &
.840E+00_EB, .728E+00_EB, .422E+00_EB, .236E+00_EB, .197E+00_EB, .167E+00_EB,     &
.505E+00_EB, .500E+00_EB, .379E+00_EB, .276E+00_EB, .227E+00_EB, .192E+00_EB,     &
.117E+00_EB, .400E+00_EB, .423E+00_EB, .315E+00_EB, .243E+00_EB, .202E+00_EB,     &
.460E-01_EB, .300E+00_EB, .358E+00_EB, .290E+00_EB, .230E+00_EB, .202E+00_EB,     &
.183E-01_EB, .205E+00_EB, .269E+00_EB, .235E+00_EB, .195E+00_EB, .192E+00_EB,     &
.730E-02_EB, .135E+00_EB, .186E+00_EB, .179E+00_EB, .159E+00_EB, .168E+00_EB,     &
.557E-02_EB, .790E-01_EB, .113E+00_EB, .124E+00_EB, .124E+00_EB, .134E+00_EB/),  & 
(/6,8/))
SD(1:6,161:168) = RESHAPE ((/ & ! 4050-4225
.283E-02_EB, .415E-01_EB, .662E-01_EB, .886E-01_EB, .103E+00_EB, .106E+00_EB,     &
.226E-02_EB, .197E-01_EB, .367E-01_EB, .594E-01_EB, .801E-01_EB, .879E-01_EB,     &
.155E-02_EB, .860E-02_EB, .211E-01_EB, .395E-01_EB, .503E-01_EB, .610E-01_EB,     &
.103E-02_EB, .521E-02_EB, .119E-01_EB, .246E-01_EB, .354E-01_EB, .480E-01_EB,     &
.821E-03_EB, .365E-02_EB, .759E-02_EB, .166E-01_EB, .258E-01_EB, .370E-01_EB,     &
.752E-03_EB, .183E-02_EB, .445E-02_EB, .100E-01_EB, .179E-01_EB, .268E-01_EB,     &
.429E-03_EB, .141E-02_EB, .354E-02_EB, .821E-02_EB, .142E-01_EB, .212E-01_EB,     &
.327E-03_EB, .902E-03_EB, .209E-02_EB, .588E-02_EB, .112E-01_EB, .172E-01_EB/),   &
(/6,8/))
SD(1:6,169:176) = RESHAPE ((/ & ! 4250-4425
.225E-03_EB, .685E-03_EB, .189E-02_EB, .512E-02_EB, .101E-01_EB, .164E-01_EB,     &
.186E-03_EB, .551E-03_EB, .156E-02_EB, .366E-02_EB, .812E-02_EB, .136E-01_EB,     &
.173E-03_EB, .472E-03_EB, .139E-02_EB, .306E-02_EB, .661E-02_EB, .115E-01_EB,     &
.138E-03_EB, .395E-03_EB, .110E-02_EB, .272E-02_EB, .587E-02_EB, .104E-01_EB,     &
.900E-04_EB, .270E-03_EB, .968E-03_EB, .222E-02_EB, .497E-02_EB, .921E-02_EB,     &
.752E-04_EB, .233E-03_EB, .744E-03_EB, .208E-02_EB, .466E-02_EB, .876E-02_EB,     &
.618E-04_EB, .175E-03_EB, .638E-03_EB, .185E-02_EB, .465E-02_EB, .914E-02_EB,     &
.504E-04_EB, .134E-03_EB, .499E-03_EB, .174E-02_EB, .455E-02_EB, .935E-02_EB/),   &
(/6,8/))
SD(1:6,177:184) = RESHAPE ((/ & ! 4450-4625
.375E-04_EB, .123E-03_EB, .485E-03_EB, .182E-02_EB, .456E-02_EB, .971E-02_EB,     &
.305E-04_EB, .892E-04_EB, .338E-03_EB, .134E-02_EB, .460E-02_EB, .104E-01_EB,     &
.257E-04_EB, .790E-04_EB, .329E-03_EB, .154E-02_EB, .477E-02_EB, .112E-01_EB,     &
.242E-04_EB, .740E-04_EB, .308E-03_EB, .135E-02_EB, .497E-02_EB, .122E-01_EB,     &
.215E-04_EB, .653E-04_EB, .282E-03_EB, .131E-02_EB, .521E-02_EB, .133E-01_EB,     &
.218E-04_EB, .660E-04_EB, .272E-03_EB, .152E-02_EB, .573E-02_EB, .148E-01_EB,     &
.215E-04_EB, .671E-04_EB, .268E-03_EB, .134E-02_EB, .607E-02_EB, .159E-01_EB,     &
.217E-04_EB, .695E-04_EB, .285E-03_EB, .161E-02_EB, .677E-02_EB, .173E-01_EB/),  & 
(/6,8/))
SD(1:6,185:192) = RESHAPE ((/ & ! 4650-4825
.219E-04_EB, .722E-04_EB, .297E-03_EB, .169E-02_EB, .783E-02_EB, .197E-01_EB,     &
.226E-04_EB, .771E-04_EB, .341E-03_EB, .236E-02_EB, .925E-02_EB, .226E-01_EB,     &
.250E-04_EB, .815E-04_EB, .387E-03_EB, .286E-02_EB, .106E-01_EB, .250E-01_EB,    &
.280E-04_EB, .845E-04_EB, .420E-03_EB, .357E-02_EB, .124E-01_EB, .276E-01_EB,     &
.351E-04_EB, .192E-03_EB, .470E-03_EB, .467E-02_EB, .166E-01_EB, .313E-01_EB,     &
.435E-04_EB, .200E-03_EB, .105E-02_EB, .566E-02_EB, .185E-01_EB, .341E-01_EB,     &
.522E-04_EB, .233E-03_EB, .129E-02_EB, .736E-02_EB, .229E-01_EB, .378E-01_EB,     &
.673E-04_EB, .306E-03_EB, .183E-02_EB, .982E-02_EB, .258E-01_EB, .404E-01_EB/),   &
(/6,8/))
SD(1:6,193:200) = RESHAPE ((/ & ! 4850-5025
.886E-04_EB, .399E-03_EB, .246E-02_EB, .128E-01_EB, .302E-01_EB, .430E-01_EB,     &
.113E-03_EB, .618E-03_EB, .346E-02_EB, .161E-01_EB, .358E-01_EB, .459E-01_EB,     &
.174E-03_EB, .825E-03_EB, .441E-02_EB, .200E-01_EB, .417E-01_EB, .493E-01_EB,     &
.265E-03_EB, .163E-02_EB, .777E-02_EB, .245E-01_EB, .450E-01_EB, .507E-01_EB,     &
.355E-03_EB, .200E-02_EB, .978E-02_EB, .317E-01_EB, .492E-01_EB, .527E-01_EB,     &
.538E-03_EB, .271E-02_EB, .167E-01_EB, .401E-01_EB, .503E-01_EB, .523E-01_EB,     &
.651E-03_EB, .301E-02_EB, .264E-01_EB, .467E-01_EB, .520E-01_EB, .526E-01_EB,     &
.987E-03_EB, .530E-02_EB, .321E-01_EB, .499E-01_EB, .523E-01_EB, .510E-01_EB/),   &
(/6,8/))
SD(1:6,201:208) = RESHAPE ((/ & ! 5050-5225
.135E-02_EB, .860E-02_EB, .389E-01_EB, .528E-01_EB, .513E-01_EB, .492E-01_EB,    &
.226E-02_EB, .130E-01_EB, .472E-01_EB, .559E-01_EB, .500E-01_EB, .469E-01_EB,     &
.431E-02_EB, .198E-01_EB, .526E-01_EB, .557E-01_EB, .480E-01_EB, .452E-01_EB,     &
.628E-02_EB, .282E-01_EB, .488E-01_EB, .495E-01_EB, .451E-01_EB, .430E-01_EB,     &
.900E-02_EB, .390E-01_EB, .471E-01_EB, .449E-01_EB, .430E-01_EB, .423E-01_EB,     &
.180E-01_EB, .462E-01_EB, .412E-01_EB, .391E-01_EB, .403E-01_EB, .415E-01_EB,     &
.348E-01_EB, .710E-01_EB, .402E-01_EB, .360E-01_EB, .384E-01_EB, .414E-01_EB,     &
.718E-01_EB, .590E-01_EB, .399E-01_EB, .360E-01_EB, .376E-01_EB, .420E-01_EB/),   &
(/6,8/))
SD(1:6,209:216) = RESHAPE ((/ & ! 5250-5425
.111E+00_EB, .368E-01_EB, .340E-01_EB, .369E-01_EB, .409E-01_EB, .454E-01_EB,     &
.329E-01_EB, .285E-01_EB, .365E-01_EB, .423E-01_EB, .461E-01_EB, .482E-01_EB,     &
.281E-01_EB, .270E-01_EB, .432E-01_EB, .505E-01_EB, .529E-01_EB, .511E-01_EB,     &
.121E+00_EB, .422E-01_EB, .589E-01_EB, .598E-01_EB, .572E-01_EB, .544E-01_EB,     &
.139E+00_EB, .105E+00_EB, .844E-01_EB, .687E-01_EB, .593E-01_EB, .560E-01_EB,     &
.774E-01_EB, .710E-01_EB, .683E-01_EB, .618E-01_EB, .556E-01_EB, .534E-01_EB,     &
.858E-01_EB, .483E-01_EB, .579E-01_EB, .547E-01_EB, .503E-01_EB, .495E-01_EB,     &
.985E-01_EB, .575E-01_EB, .589E-01_EB, .510E-01_EB, .451E-01_EB, .449E-01_EB/),   &
(/6,8/))
SD(1:6,217:224) = RESHAPE ((/ & ! 5450-5625
.996E-01_EB, .682E-01_EB, .539E-01_EB, .489E-01_EB, .454E-01_EB, .446E-01_EB,     &
.680E-01_EB, .680E-01_EB, .548E-01_EB, .495E-01_EB, .460E-01_EB, .458E-01_EB,     &
.325E-01_EB, .520E-01_EB, .515E-01_EB, .483E-01_EB, .449E-01_EB, .454E-01_EB,     &
.150E-01_EB, .350E-01_EB, .451E-01_EB, .464E-01_EB, .452E-01_EB, .449E-01_EB,     &
.620E-02_EB, .238E-01_EB, .369E-01_EB, .408E-01_EB, .414E-01_EB, .417E-01_EB,     &
.270E-02_EB, .158E-01_EB, .282E-01_EB, .339E-01_EB, .366E-01_EB, .384E-01_EB,     &
.113E-02_EB, .101E-01_EB, .203E-01_EB, .263E-01_EB, .303E-01_EB, .333E-01_EB,     &
.829E-03_EB, .590E-02_EB, .148E-01_EB, .206E-01_EB, .247E-01_EB, .295E-01_EB/),  & 
(/6,8/))
SD(1:6,225:232) = RESHAPE ((/ & ! 5650-5825
.365E-03_EB, .310E-02_EB, .969E-02_EB, .154E-01_EB, .203E-01_EB, .258E-01_EB,     &
.240E-03_EB, .130E-02_EB, .589E-02_EB, .112E-01_EB, .164E-01_EB, .222E-01_EB,     &
.158E-03_EB, .400E-03_EB, .417E-02_EB, .850E-02_EB, .134E-01_EB, .190E-01_EB,     &
.103E-03_EB, .262E-03_EB, .208E-02_EB, .594E-02_EB, .109E-01_EB, .162E-01_EB,     &
.741E-04_EB, .181E-03_EB, .142E-02_EB, .455E-02_EB, .907E-02_EB, .141E-01_EB,     &
.625E-04_EB, .135E-03_EB, .816E-03_EB, .316E-02_EB, .698E-02_EB, .121E-01_EB,     &
.499E-04_EB, .111E-03_EB, .624E-03_EB, .230E-02_EB, .551E-02_EB, .102E-01_EB,     &
.325E-04_EB, .677E-04_EB, .425E-03_EB, .124E-02_EB, .385E-02_EB, .818E-02_EB/),   &
(/6,8/))
SD(1:6,233:240) = RESHAPE ((/ & ! 5850-6025
.231E-04_EB, .563E-04_EB, .278E-03_EB, .986E-03_EB, .290E-02_EB, .672E-02_EB,    &
.165E-04_EB, .481E-04_EB, .247E-03_EB, .944E-03_EB, .253E-02_EB, .612E-02_EB,           &
.126E-04_EB, .432E-04_EB, .241E-03_EB, .886E-03_EB, .220E-02_EB, .582E-02_EB,          &
.118E-04_EB, .420E-04_EB, .235E-03_EB, .847E-03_EB, .209E-02_EB, .571E-02_EB,          &
.110E-04_EB, .408E-04_EB, .226E-03_EB, .812E-03_EB, .221E-02_EB, .604E-02_EB,         &
.101E-04_EB, .400E-04_EB, .213E-03_EB, .805E-03_EB, .239E-02_EB, .641E-02_EB,        &
.983E-05_EB, .395E-04_EB, .186E-03_EB, .801E-03_EB, .247E-02_EB, .691E-02_EB,       &
.979E-05_EB, .401E-04_EB, .193E-03_EB, .805E-03_EB, .260E-02_EB, .732E-02_EB/),    &
(/6,8/))
SD(1:6,241:248) = RESHAPE ((/ & ! 6050-6225
.976E-05_EB, .410E-04_EB, .201E-03_EB, .814E-03_EB, .285E-02_EB, .776E-02_EB,     &
.988E-05_EB, .420E-04_EB, .210E-03_EB, .832E-03_EB, .317E-02_EB, .842E-02_EB,    &
.991E-05_EB, .425E-04_EB, .219E-03_EB, .877E-03_EB, .340E-02_EB, .888E-02_EB,        &
.102E-04_EB, .435E-04_EB, .231E-03_EB, .937E-03_EB, .361E-02_EB, .929E-02_EB,       &
.110E-04_EB, .486E-04_EB, .244E-03_EB, .971E-03_EB, .402E-02_EB, .994E-02_EB,      &
.127E-04_EB, .579E-04_EB, .257E-03_EB, .111E-02_EB, .437E-02_EB, .104E-01_EB,     &
.131E-04_EB, .612E-04_EB, .277E-03_EB, .113E-02_EB, .465E-02_EB, .110E-01_EB,        &
.150E-04_EB, .783E-04_EB, .353E-03_EB, .116E-02_EB, .510E-02_EB, .116E-01_EB/),     &
(/6,8/))
SD(1:6,249:256) = RESHAPE ((/ & ! 6250-6425
.178E-04_EB, .922E-04_EB, .394E-03_EB, .157E-02_EB, .555E-02_EB, .123E-01_EB,      &
.203E-04_EB, .115E-03_EB, .481E-03_EB, .188E-02_EB, .601E-02_EB, .131E-01_EB,     &
.230E-04_EB, .145E-03_EB, .617E-03_EB, .183E-02_EB, .644E-02_EB, .139E-01_EB,      & 
.280E-04_EB, .187E-03_EB, .723E-03_EB, .202E-02_EB, .686E-02_EB, .146E-01_EB,     & 
.305E-04_EB, .209E-03_EB, .811E-03_EB, .243E-02_EB, .779E-02_EB, .157E-01_EB,     &
.455E-04_EB, .244E-03_EB, .935E-03_EB, .243E-02_EB, .844E-02_EB, .166E-01_EB,    &
.661E-04_EB, .320E-03_EB, .989E-03_EB, .288E-02_EB, .902E-02_EB, .173E-01_EB,       & 
.723E-04_EB, .397E-03_EB, .122E-02_EB, .359E-02_EB, .100E-01_EB, .184E-01_EB/),    &
(/6,8/))
SD(1:6,257:264) = RESHAPE ((/ & ! 6450-6625
.847E-04_EB, .481E-03_EB, .143E-02_EB, .429E-02_EB, .108E-01_EB, .192E-01_EB,     &
.103E-03_EB, .591E-03_EB, .174E-02_EB, .488E-02_EB, .116E-01_EB, .200E-01_EB,           &
.131E-03_EB, .703E-03_EB, .247E-02_EB, .549E-02_EB, .124E-01_EB, .205E-01_EB,          &
.165E-03_EB, .872E-03_EB, .265E-02_EB, .641E-02_EB, .131E-01_EB, .211E-01_EB,         &
.205E-03_EB, .110E-02_EB, .298E-02_EB, .749E-02_EB, .140E-01_EB, .218E-01_EB,         &
.253E-03_EB, .130E-02_EB, .346E-02_EB, .811E-02_EB, .150E-01_EB, .230E-01_EB,         &
.338E-03_EB, .150E-02_EB, .445E-02_EB, .890E-02_EB, .159E-01_EB, .237E-01_EB,        &
.437E-03_EB, .170E-02_EB, .491E-02_EB, .107E-01_EB, .170E-01_EB, .245E-01_EB/),     &
(/6,8/))
SD(1:6,265:272) = RESHAPE ((/ & ! 6650-6825
.581E-03_EB, .190E-02_EB, .537E-02_EB, .116E-01_EB, .179E-01_EB, .254E-01_EB,         &
.685E-03_EB, .220E-02_EB, .578E-02_EB, .128E-01_EB, .189E-01_EB, .263E-01_EB,        &
.900E-03_EB, .250E-02_EB, .649E-02_EB, .134E-01_EB, .195E-01_EB, .275E-01_EB,       &
.121E-02_EB, .280E-02_EB, .722E-02_EB, .142E-01_EB, .202E-01_EB, .281E-01_EB,      &
.152E-02_EB, .330E-02_EB, .813E-02_EB, .161E-01_EB, .212E-01_EB, .288E-01_EB,     &
.185E-02_EB, .370E-02_EB, .907E-02_EB, .168E-01_EB, .222E-01_EB, .292E-01_EB,    &
.220E-02_EB, .430E-02_EB, .929E-02_EB, .183E-01_EB, .233E-01_EB, .294E-01_EB,       &
.255E-02_EB, .500E-02_EB, .114E-01_EB, .195E-01_EB, .245E-01_EB, .289E-01_EB/),       &
(/6,8/))
SD(1:6,273:280) = RESHAPE ((/ & ! 6850-7025
.290E-02_EB, .580E-02_EB, .167E-01_EB, .215E-01_EB, .260E-01_EB, .291E-01_EB,        &
.320E-02_EB, .670E-02_EB, .208E-01_EB, .237E-01_EB, .274E-01_EB, .293E-01_EB,        &
.360E-02_EB, .880E-02_EB, .220E-01_EB, .253E-01_EB, .282E-01_EB, .300E-01_EB,       &
.400E-02_EB, .920E-02_EB, .238E-01_EB, .273E-01_EB, .290E-01_EB, .304E-01_EB,        &
.460E-02_EB, .108E-01_EB, .272E-01_EB, .279E-01_EB, .298E-01_EB, .310E-01_EB,       &
.530E-02_EB, .128E-01_EB, .304E-01_EB, .292E-01_EB, .297E-01_EB, .312E-01_EB,        &
.620E-02_EB, .152E-01_EB, .344E-01_EB, .303E-01_EB, .293E-01_EB, .310E-01_EB,       &
.760E-02_EB, .182E-01_EB, .341E-01_EB, .297E-01_EB, .290E-01_EB, .300E-01_EB/),    &
(/6,8/))
SD(1:6,281:288) = RESHAPE ((/ & ! 7050-7225
.980E-02_EB, .222E-01_EB, .398E-01_EB, .318E-01_EB, .291E-01_EB, .294E-01_EB,       & 
.132E-01_EB, .271E-01_EB, .402E-01_EB, .294E-01_EB, .274E-01_EB, .282E-01_EB,      & 
.190E-01_EB, .335E-01_EB, .421E-01_EB, .286E-01_EB, .262E-01_EB, .269E-01_EB,     & 
.240E-01_EB, .432E-01_EB, .431E-01_EB, .276E-01_EB, .245E-01_EB, .257E-01_EB,    & 
.288E-01_EB, .570E-01_EB, .458E-01_EB, .270E-01_EB, .228E-01_EB, .243E-01_EB,   &
.323E-01_EB, .740E-01_EB, .449E-01_EB, .261E-01_EB, .214E-01_EB, .221E-01_EB,        &
.570E-01_EB, .890E-01_EB, .435E-01_EB, .225E-01_EB, .199E-01_EB, .196E-01_EB,       &
.216E-01_EB, .680E-01_EB, .378E-01_EB, .239E-01_EB, .195E-01_EB, .192E-01_EB/),    &
(/6,8/))
SD(1:6,289:296) = RESHAPE ((/ & ! 7250-7425
.126E-01_EB, .475E-01_EB, .364E-01_EB, .238E-01_EB, .197E-01_EB, .192E-01_EB,      & 
.117E-01_EB, .369E-01_EB, .385E-01_EB, .249E-01_EB, .212E-01_EB, .204E-01_EB,     & 
.140E-01_EB, .370E-01_EB, .419E-01_EB, .272E-01_EB, .228E-01_EB, .213E-01_EB,    & 
.425E-01_EB, .418E-01_EB, .440E-01_EB, .280E-01_EB, .248E-01_EB, .229E-01_EB,   &
.640E-01_EB, .460E-01_EB, .427E-01_EB, .290E-01_EB, .263E-01_EB, .238E-01_EB,   &
.385E-01_EB, .385E-01_EB, .374E-01_EB, .259E-01_EB, .235E-01_EB, .224E-01_EB,        &
.182E-01_EB, .179E-01_EB, .282E-01_EB, .231E-01_EB, .211E-01_EB, .214E-01_EB,       &
.170E-01_EB, .810E-02_EB, .191E-01_EB, .175E-01_EB, .181E-01_EB, .194E-01_EB/),     & 
(/6,8/))
SD(1:6,297:304) = RESHAPE ((/ & ! 7450-7625
.161E-01_EB, .370E-02_EB, .105E-01_EB, .127E-01_EB, .152E-01_EB, .171E-01_EB,      & 
.145E-01_EB, .170E-02_EB, .554E-02_EB, .855E-02_EB, .113E-01_EB, .131E-01_EB,     & 
.175E-02_EB, .140E-02_EB, .385E-02_EB, .595E-02_EB, .803E-02_EB, .945E-02_EB,    &
.772E-03_EB, .751E-03_EB, .384E-02_EB, .575E-02_EB, .537E-02_EB, .594E-02_EB,       &
.491E-03_EB, .600E-03_EB, .301E-02_EB, .453E-02_EB, .380E-02_EB, .434E-02_EB,      &
.275E-03_EB, .410E-03_EB, .193E-02_EB, .366E-02_EB, .319E-02_EB, .332E-02_EB,     &
.185E-01_EB, .280E-03_EB, .131E-02_EB, .232E-02_EB, .247E-02_EB, .256E-02_EB,      &
.101E-03_EB, .160E-03_EB, .915E-03_EB, .150E-02_EB, .186E-02_EB, .197E-02_EB/),   &
(/6,8/))
SD(1:6,305:312) = RESHAPE ((/ & ! 7650-7825
.691E-04_EB, .110E-03_EB, .565E-03_EB, .114E-02_EB, .205E-02_EB, .192E-02_EB,        &  
.476E-04_EB, .750E-04_EB, .114E-02_EB, .124E-02_EB, .175E-02_EB, .187E-02_EB,         &
.305E-04_EB, .590E-04_EB, .529E-03_EB, .114E-02_EB, .160E-02_EB, .185E-02_EB,       &
.240E-04_EB, .480E-04_EB, .293E-03_EB, .842E-03_EB, .141E-02_EB, .184E-02_EB,       &
.170E-04_EB, .360E-04_EB, .122E-03_EB, .435E-03_EB, .124E-02_EB, .182E-02_EB,      &
.120E-04_EB, .240E-04_EB, .121E-03_EB, .435E-03_EB, .118E-02_EB, .187E-02_EB,     &
.810E-05_EB, .170E-04_EB, .103E-03_EB, .439E-03_EB, .126E-02_EB, .192E-02_EB,    &
.550E-05_EB, .120E-04_EB, .866E-04_EB, .367E-03_EB, .119E-02_EB, .193E-02_EB/),    &
(/6,8/))
SD(1:6,313:320) = RESHAPE ((/&  ! 7850-8025
.390E-05_EB, .900E-05_EB, .716E-04_EB, .351E-03_EB, .116E-02_EB, .194E-02_EB,        & 
.295E-05_EB, .830E-05_EB, .373E-04_EB, .254E-03_EB, .114E-02_EB, .196E-02_EB,       & 
.230E-05_EB, .800E-05_EB, .465E-04_EB, .298E-03_EB, .117E-02_EB, .201E-02_EB,      & 
.225E-05_EB, .820E-05_EB, .367E-04_EB, .252E-03_EB, .116E-02_EB, .205E-02_EB,     & 
.220E-05_EB, .840E-05_EB, .371E-04_EB, .268E-03_EB, .127E-02_EB, .211E-02_EB,    & 
.223E-05_EB, .920E-05_EB, .396E-04_EB, .273E-03_EB, .128E-02_EB, .216E-02_EB,   &
.235E-05_EB, .103E-04_EB, .415E-04_EB, .263E-03_EB, .121E-02_EB, .221E-02_EB,      &
.280E-05_EB, .125E-04_EB, .633E-04_EB, .363E-03_EB, .136E-02_EB, .231E-02_EB/),   &
(/6,8/))
SD(1:6,321:328) = RESHAPE ((/ & ! 8050-8225
.310E-05_EB, .150E-04_EB, .979E-04_EB, .492E-03_EB, .150E-02_EB, .241E-02_EB,           &
.370E-05_EB, .180E-04_EB, .120E-03_EB, .580E-03_EB, .167E-02_EB, .251E-02_EB,           &
.420E-05_EB, .200E-04_EB, .987E-04_EB, .509E-03_EB, .171E-02_EB, .257E-02_EB,           &
.510E-05_EB, .240E-04_EB, .134E-03_EB, .547E-03_EB, .173E-02_EB, .267E-02_EB,           &
.600E-05_EB, .270E-04_EB, .121E-03_EB, .534E-03_EB, .172E-02_EB, .274E-02_EB,           &
.720E-05_EB, .300E-04_EB, .204E-03_EB, .684E-03_EB, .184E-02_EB, .285E-02_EB,           &
.820E-05_EB, .330E-04_EB, .276E-03_EB, .819E-03_EB, .199E-02_EB, .297E-02_EB,           &
.100E-04_EB, .380E-04_EB, .317E-03_EB, .859E-03_EB, .214E-02_EB, .308E-02_EB/),         &
(/6,8/))
SD(1:6,329:336) = RESHAPE ((/ & ! 8250-8425
.125E-04_EB, .420E-04_EB, .240E-03_EB, .818E-03_EB, .220E-02_EB, .317E-02_EB,           &
.145E-04_EB, .500E-04_EB, .452E-03_EB, .109E-02_EB, .238E-02_EB, .293E-02_EB,           &
.175E-04_EB, .560E-04_EB, .301E-03_EB, .941E-03_EB, .243E-02_EB, .342E-02_EB,           &
.198E-04_EB, .630E-04_EB, .280E-03_EB, .107E-02_EB, .260E-02_EB, .353E-02_EB,           &
.230E-04_EB, .710E-04_EB, .276E-03_EB, .109E-02_EB, .272E-02_EB, .365E-02_EB,           &
.280E-04_EB, .830E-04_EB, .369E-03_EB, .127E-02_EB, .295E-02_EB, .377E-02_EB,           &
.330E-04_EB, .890E-04_EB, .430E-03_EB, .139E-02_EB, .306E-02_EB, .385E-02_EB,           &
.360E-04_EB, .950E-04_EB, .371E-03_EB, .135E-02_EB, .306E-02_EB, .384E-02_EB/),         &
(/6,8/))
SD(1:6,337:344) = RESHAPE ((/ & ! 8450-8625
.390E-04_EB, .980E-04_EB, .434E-03_EB, .147E-02_EB, .316E-02_EB, .385E-02_EB,           &
.400E-04_EB, .990E-04_EB, .397E-03_EB, .143E-02_EB, .318E-02_EB, .384E-02_EB,           &
.400E-04_EB, .980E-04_EB, .364E-03_EB, .141E-02_EB, .317E-02_EB, .381E-02_EB,           &
.390E-04_EB, .940E-04_EB, .390E-03_EB, .142E-02_EB, .314E-02_EB, .376E-02_EB,          &
.380E-04_EB, .900E-04_EB, .380E-03_EB, .145E-02_EB, .318E-02_EB, .375E-02_EB,         &
.380E-04_EB, .900E-04_EB, .380E-03_EB, .145E-02_EB, .318E-02_EB, .375E-02_EB,        &
.330E-04_EB, .750E-04_EB, .358E-03_EB, .138E-02_EB, .310E-02_EB, .372E-02_EB,       &
.270E-04_EB, .580E-04_EB, .382E-03_EB, .143E-02_EB, .315E-02_EB, .369E-02_EB/),    &
(/6,8/))
SD(1:6,345:352) = RESHAPE ((/ & ! 8650-8825
.240E-04_EB, .500E-04_EB, .343E-03_EB, .136E-02_EB, .306E-02_EB, .363E-02_EB,          &
.200E-04_EB, .450E-04_EB, .309E-03_EB, .134E-02_EB, .306E-02_EB, .359E-02_EB,         & 
.180E-04_EB, .400E-04_EB, .281E-03_EB, .127E-02_EB, .294E-02_EB, .341E-02_EB,        &
.170E-04_EB, .360E-04_EB, .276E-03_EB, .124E-02_EB, .290E-02_EB, .336E-02_EB,       & 
.160E-04_EB, .310E-04_EB, .272E-03_EB, .122E-02_EB, .283E-02_EB, .323E-02_EB,      &  
.140E-04_EB, .280E-04_EB, .241E-03_EB, .117E-02_EB, .273E-02_EB, .309E-02_EB,       &
.120E-04_EB, .250E-04_EB, .237E-03_EB, .115E-02_EB, .269E-02_EB, .297E-02_EB,      &
.100E-04_EB, .220E-04_EB, .218E-03_EB, .111E-02_EB, .259E-02_EB, .284E-02_EB/),   &
(/6,8/))
SD(1:6,353:360) = RESHAPE ((/ & ! 8850-9025
.920E-05_EB, .198E-04_EB, .206E-03_EB, .105E-02_EB, .246E-02_EB, .269E-02_EB,       &
.810E-05_EB, .170E-04_EB, .205E-03_EB, .100E-02_EB, .235E-02_EB, .257E-02_EB,       & 
.720E-05_EB, .160E-04_EB, .177E-03_EB, .921E-03_EB, .220E-02_EB, .245E-02_EB,      &
.650E-05_EB, .150E-04_EB, .172E-03_EB, .834E-03_EB, .205E-02_EB, .232E-02_EB,     &
.590E-05_EB, .130E-04_EB, .147E-03_EB, .735E-03_EB, .194E-02_EB, .218E-02_EB,    &
.510E-05_EB, .110E-04_EB, .120E-03_EB, .629E-03_EB, .177E-02_EB, .203E-02_EB,     &
.460E-05_EB, .950E-05_EB, .960E-04_EB, .513E-03_EB, .154E-02_EB, .180E-02_EB,      &
.420E-05_EB, .800E-05_EB, .578E-04_EB, .314E-03_EB, .123E-02_EB, .154E-02_EB/),     &
(/6,8/))
SD(1:6,361:368) = RESHAPE ((/ & ! 9050-9225
.380E-05_EB, .720E-05_EB, .529E-04_EB, .292E-03_EB, .114E-02_EB, .137E-02_EB,      &
.330E-05_EB, .660E-05_EB, .485E-04_EB, .269E-03_EB, .102E-02_EB, .122E-02_EB,      &
.290E-05_EB, .580E-05_EB, .430E-04_EB, .239E-03_EB, .896E-03_EB, .107E-02_EB,      &
.270E-05_EB, .520E-05_EB, .259E-04_EB, .193E-03_EB, .748E-03_EB, .944E-03_EB,     &
.240E-05_EB, .450E-05_EB, .316E-04_EB, .207E-03_EB, .671E-02_EB, .848E-03_EB,      &
.220E-05_EB, .400E-05_EB, .444E-05_EB, .602E-04_EB, .516E-03_EB, .750E-03_EB,     &
.190E-05_EB, .360E-05_EB, .324E-05_EB, .460E-04_EB, .439E-03_EB, .688E-03_EB,     & 
.170E-05_EB, .320E-05_EB, .180E-05_EB, .321E-04_EB, .384E-03_EB, .653E-03_EB/),  & 
(/6,8/))
SD(1:6,369:376) = RESHAPE ((/ & ! 9250-9300
.140E-05_EB, .280E-05_EB, .171E-05_EB, .344E-04_EB, .340E-03_EB, .616E-03_EB,   & 
.130E-05_EB, .250E-05_EB, .299E-05_EB, .600E-04_EB, .343E-03_EB, .619E-03_EB,  &  
.120E-05_EB, .220E-05_EB, .299E-05_EB, .600E-04_EB, .343E-03_EB, .619E-03_EB, &  
1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB,&
1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB, 1._EB/),&
(/6,8/))

! Intialize GAMMA Array

GAMMA = RESHAPE((/ &
!     LINE BROADENING PARAMETERS,GAMMA(I,J),
!     J=CO2,H2O,CH4,CO,O2,N2,SELF RESONANT.
!   I=     CO2  H2O  CH4   CO
!    J
.09_EB , .12_EB, .0_EB , .07_EB,&
.07_EB , .09_EB, .0_EB , .06_EB,&
!.0_EB  , .0_EB , .16_EB, .0_EB ,&
.0_EB  , .0_EB , .12_EB, .0_EB ,&
.06_EB , .10_EB, .0_EB , .06_EB,&
.055_EB, .04_EB, .0_EB , .05_EB,&
.07_EB , .09_EB, .0_EB , .06_EB, &
.01_EB , .44_EB, .0_EB , .0_EB /),(/4,7/))

! Initialize SD15 Array

!  THE FOLLOWING ARE DATA FOR THE 15.0 MICRON BAND OF CO2
!  TEMP, K=300     600      1200      1500  1800    2400    

SD15(1:6,1:8) = RESHAPE ((/  &! 500-535
.000E+00_EB, .000E+00_EB, .000E+00_EB, .105E-01_EB, .300E-01_EB, .880E-01_EB,          &
.000E+00_EB, .000E+00_EB, .000E+00_EB, .180E-01_EB, .490E-01_EB, .880E-01_EB,         & 
.000E+00_EB, .000E+00_EB, .000E+00_EB, .300E-01_EB, .540E-01_EB, .740E-01_EB,        &  
.000E+00_EB, .000E+00_EB, .000E+00_EB, .300E-01_EB, .560E-01_EB, .890E-01_EB,       &   
.000E+00_EB, .000E+00_EB, .000E+00_EB, .330E-01_EB, .690E-01_EB, .990E-01_EB,      &    
.000E+00_EB, .000E+00_EB, .880E-02_EB, .380E-01_EB, .720E-01_EB, .970E-01_EB,     &     
.000E+00_EB, .000E+00_EB, .110E-01_EB, .530E-01_EB, .950E-01_EB, .124E+00_EB,    &      
.000E+00_EB, .000E+00_EB, .285E-01_EB, .630E-01_EB, .990E-01_EB, .140E+00_EB/), &       
(/6,8/))
SD15(1:6,9:16) = RESHAPE ((/  &! 540-575
.000E+00_EB, .000E+00_EB, .330E-01_EB, .680E-01_EB, .103E+00_EB, .134E+00_EB,          &
.000E+00_EB, .000E+00_EB, .450E-01_EB, .920E-01_EB, .138E+00_EB, .176E+00_EB,          &
.000E+00_EB, .000E+00_EB, .490E-01_EB, .970E-01_EB, .148E+00_EB, .191E+00_EB,         & 
.000E+00_EB, .000E+00_EB, .490E-01_EB, .120E-01_EB, .188E+00_EB, .247E+00_EB,        &  
.000E+00_EB, .000E+00_EB, .480E-01_EB, .126E+00_EB, .201E+00_EB, .241E+00_EB,       &   
.000E+00_EB, .000E+00_EB, .820E-01_EB, .198E+00_EB, .270E+00_EB, .265E+00_EB,      &    
.000E+00_EB, .750E-02_EB, .690E-01_EB, .140E+00_EB, .225E+00_EB, .340E+00_EB,     &     
.000E+00_EB, .205E-01_EB, .820E-01_EB, .145E+00_EB, .236E+00_EB, .530E+00_EB/),  &      
(/6,8/))
SD15(1:6,17:24) = RESHAPE ((/ & ! 580-615
.000E+00_EB, .355E-01_EB, .117E+00_EB, .193E+00_EB, .295E+00_EB, .550E+00_EB,          &
.157E-01_EB, .520E-01_EB, .170E+00_EB, .235E+00_EB, .305E+00_EB, .410E+00_EB,          &
.150E-01_EB, .880E-01_EB, .270E+00_EB, .330E+00_EB, .440E+00_EB, .520E+00_EB,          &
.510E-01_EB, .130E+00_EB, .400E+00_EB, .530E+00_EB, .560E+00_EB, .540E+00_EB,         & 
.120E+00_EB, .165E+00_EB, .275E+00_EB, .320E+00_EB, .420E+00_EB, .560E+00_EB,        &  
.880E-01_EB, .190E+00_EB, .430E+00_EB, .540E+00_EB, .620E+00_EB, .680E+00_EB,       &   
.110E+00_EB, .350E+00_EB, .710E+00_EB, .760E+00_EB, .760E+00_EB, .690E+00_EB,      &    
.180E+00_EB, .470E+00_EB, .920E+00_EB, .970E+00_EB, .910E+00_EB, .670E+00_EB/),   &     
(/6,8/))
SD15(1:6,25:32) = RESHAPE ((/ & ! 620-655
.970E-01_EB, .265E+00_EB, .610E+00_EB, .720E+00_EB, .780E+00_EB, .730E+00_EB,          &
.175E+00_EB, .380E+00_EB, .720E+00_EB, .790E+00_EB, .830E+00_EB, .840E+00_EB,          &
.370E+00_EB, .640E+00_EB, .920E+00_EB, .960E+00_EB, .980E+00_EB, .940E+00_EB,          &
.590E+00_EB, .840E+00_EB, .107E+01_EB, .110E+01_EB, .111E+01_EB, .106E+01_EB,         & 
.940E+00_EB, .103E+01_EB, .115E+01_EB, .115E+01_EB, .115E+01_EB, .118E+01_EB,        &  
.196E+01_EB, .177E+01_EB, .146E+01_EB, .136E+01_EB, .132E+01_EB, .139E+01_EB,       &   
.345E+01_EB, .282E+01_EB, .198E+01_EB, .172E+01_EB, .156E+01_EB, .148E+01_EB,      &    
.282E+01_EB, .248E+01_EB, .200E+01_EB, .190E+01_EB, .186E+01_EB, .205E+01_EB/),   &     
(/6,8/))
SD15(1:6,33:40) = RESHAPE ((/  &! 660-695
.254E+01_EB, .234E+01_EB, .184E+01_EB, .176E+01_EB, .174E+01_EB, .203E+01_EB,          &
.142E+02_EB, .860E+01_EB, .370E+01_EB, .260E+01_EB, .196E+01_EB, .142E+01_EB,          &
.450E+01_EB, .570E+01_EB, .580E+01_EB, .520E+01_EB, .350E+01_EB, .420E+01_EB,         & 
.360E+01_EB, .310E+01_EB, .330E+01_EB, .290E+01_EB, .205E+01_EB, .200E+01_EB,        &  
.310E+01_EB, .260E+01_EB, .200E+01_EB, .196E+01_EB, .180E+01_EB, .210E+01_EB,       &   
.240E+01_EB, .250E+01_EB, .230E+01_EB, .220E+01_EB, .170E+01_EB, .194E+01_EB,      &    
.182E+01_EB, .200E+01_EB, .218E+01_EB, .205E+01_EB, .184E+01_EB, .130E+01_EB,     &     
.104E+01_EB, .135E+01_EB, .172E+01_EB, .172E+01_EB, .165E+01_EB, .130E+01_EB/),  &      
(/6,8/))
SD15(1:6,41:48) = RESHAPE ((/  &! 700-735
.550E+00_EB, .120E+01_EB, .143E+01_EB, .147E+01_EB, .148E+01_EB, .125E+01_EB,          &
.136E+01_EB, .128E+01_EB, .128E+01_EB, .135E+01_EB, .138E+01_EB, .134E+01_EB,          &
.210E+00_EB, .780E+00_EB, .127E+01_EB, .133E+01_EB, .137E+01_EB, .132E+01_EB,          &
.190E+00_EB, .780E+00_EB, .140E+01_EB, .146E+01_EB, .147E+01_EB, .142E+01_EB,         & 
.900E+00_EB, .106E+01_EB, .140E+01_EB, .150E+01_EB, .155E+01_EB, .134E+01_EB,        &  
.720E-01_EB, .300E+00_EB, .800E+00_EB, .100E+01_EB, .115E+01_EB, .126E+01_EB,       &   
.640E-01_EB, .210E+00_EB, .560E+00_EB, .720E+00_EB, .860E+00_EB, .102E+01_EB,      &    
.680E-01_EB, .210E+00_EB, .530E+00_EB, .670E+00_EB, .790E+00_EB, .101E+01_EB/),   &     
(/6,8/))
SD15(1:6,49:56) = RESHAPE ((/ & ! 740-775
.690E-01_EB, .210E+00_EB, .540E+00_EB, .690E+00_EB, .820E+00_EB, .910E+00_EB,          &
.330E-01_EB, .140E+00_EB, .390E+00_EB, .530E+00_EB, .690E+00_EB, .770E+00_EB,          &
.230E-01_EB, .780E-01_EB, .270E+00_EB, .410E+00_EB, .560E+00_EB, .890E+00_EB,          &
.300E-01_EB, .860E-01_EB, .280E+00_EB, .400E+00_EB, .520E+00_EB, .710E+00_EB,         & 
.175E-01_EB, .620E-01_EB, .225E+00_EB, .335E+00_EB, .450E+00_EB, .660E+00_EB,       &   
.105E-01_EB, .450E-01_EB, .180E+00_EB, .280E+00_EB, .380E+00_EB, .600E+00_EB,      &    
.450E-02_EB, .300E-01_EB, .148E+00_EB, .240E+00_EB, .345E+00_EB, .570E+00_EB,     &     
.000E+00_EB, .140E-01_EB, .124E+00_EB, .205E+00_EB, .285E+00_EB, .430E+00_EB/),  &      
(/6,8/))
SD15(1:6,57:64) = RESHAPE ((/ & ! 780-815
.000E+00_EB, .115E-01_EB, .110E+00_EB, .185E+00_EB, .260E+00_EB, .375E+00_EB,          &
.000E+00_EB, .135E-01_EB, .840E-01_EB, .140E+00_EB, .205E+00_EB, .335E+00_EB,         & 
.000E+00_EB, .430E-02_EB, .650E-01_EB, .120E+00_EB, .185E+00_EB, .325E+00_EB,        &  
.000E+00_EB, .000E+00_EB, .540E-01_EB, .115E+00_EB, .180E+00_EB, .315E+00_EB,       &   
.000E+00_EB, .000E+00_EB, .440E-01_EB, .950E-01_EB, .150E+00_EB, .270E+00_EB,      &    
.000E+00_EB, .000E+00_EB, .360E-01_EB, .790E-01_EB, .125E+00_EB, .205E+00_EB,     &     
.000E+00_EB, .000E+00_EB, .250E-01_EB, .650E-01_EB, .110E+00_EB, .178E+00_EB,    &      
.000E+00_EB, .000E+00_EB, .180E-01_EB, .620E-01_EB, .103E+00_EB, .153E+00_EB/), &       
(/6,8/))
SD15(1:6,65:72) = RESHAPE ((/ & ! 820-855
.000E+00_EB, .000E+00_EB, .320E-01_EB, .580E-01_EB, .860E-01_EB, .147E+00_EB,          &
.000E+00_EB, .000E+00_EB, .800E-02_EB, .510E-01_EB, .870E-01_EB, .134E+00_EB,          &
.000E+00_EB, .000E+00_EB, .600E-02_EB, .480E-01_EB, .830E-01_EB, .133E+00_EB,          &
.000E+00_EB, .000E+00_EB, .000E+00_EB, .430E-01_EB, .780E-01_EB, .118E+00_EB,          &
.000E+00_EB, .000E+00_EB, .000E+00_EB, .420E-01_EB, .700E-01_EB, .108E+00_EB,         & 
.000E+00_EB, .000E+00_EB, .000E+00_EB, .360E-01_EB, .640E-01_EB, .980E-01_EB,        &  
.000E+00_EB, .000E+00_EB, .000E+00_EB, .350E-01_EB, .610E-01_EB, .870E-01_EB,       &   
.000E+00_EB, .000E+00_EB, .000E+00_EB, .320E-01_EB, .580E-01_EB, .860E-01_EB/),    &    
(/6,8/))
SD15(1:6,73:80) = RESHAPE ((/ & ! 860-880
.000E+00_EB, .000E+00_EB, .000E+00_EB, .330E-01_EB, .560E-01_EB, .750E-01_EB,     &     
.000E+00_EB, .000E+00_EB, .000E+00_EB, .300E-01_EB, .530E-01_EB, .750E-01_EB,    &      
.000E+00_EB, .000E+00_EB, .000E+00_EB, .290E-01_EB, .530E-01_EB, .850E-01_EB,   &       
.000E+00_EB, .000E+00_EB, .000E+00_EB, .240E-01_EB, .470E-01_EB, .900E-01_EB,  &        
.000E+00_EB, .000E+00_EB, .000E+00_EB, .220E-01_EB, .450E-01_EB, .860E-01_EB, &         
.000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB,&
.000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB,&
.000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB, .000E+00_EB/), &
(/6,8/))

! Initialize SD7 Array

!   THE FOLLOWING DATA ARE FOR THE 7.7 MICRON BAND OF CH4
! TEMP,K= 290 600   850

SD7(1:3,1:8) = RESHAPE ((/ &
0._EB, 0._EB, 0._EB,&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.22_EB,&
0.16_EB, 0.20_EB, 0.47_EB,&
0.34_EB, 0.34_EB, 0.62_EB,&
0.69_EB, 0.53_EB, 0.65_EB,&
1.27_EB, 0.88_EB, 1.09_EB,&
1.68_EB, 1.38_EB, 0.87_EB/),(/3,8/))
SD7(1:3,9:16) = RESHAPE ((/&
0.55_EB, 0.28_EB, 0.40_EB,&
1.25_EB, 0.86_EB, 0.93_EB,&
0.34_EB, 0.59_EB, 0.75_EB,&
0._EB, 0.13_EB, 0.25_EB,&
0._EB, 0._EB, 0.06_EB,&
0._EB, 0._EB, 0._EB,&
0._EB, 0._EB, 0._EB,&
0._EB, 0._EB, 0._EB/),(/3,8/))

! Initialize SD3 Array

!  THE FOLLOWING DATA ARE FOR THE 3.3 MICRON BAND OF CH4
! TEMP, K= 290  600   850
SD3(1:3,1:8) = RESHAPE ((/&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.03_EB,&
0._EB, 0._EB, 0.06_EB,&
0.03_EB, 0.03_EB, 0.09_EB,&
0.07_EB, 0.07_EB, 0.12_EB,&
0.09_EB, 0.09_EB, 0.12_EB,&
0.14_EB, 0.15_EB, 0.22_EB/),(/3,8/))
SD3(1:3,9:16) = RESHAPE ((/&
0.18_EB, 0.22_EB, 0.28_EB,&
0.24_EB, 0.31_EB, 0.37_EB,&
0.33_EB, 0.44_EB, 0.47_EB,&
0.45_EB, 0.50_EB, 0.53_EB,&
0.59_EB, 0.62_EB, 0.62_EB,&
0.74_EB, 0.70_EB, 0.68_EB,&
0.91_EB, 0.77_EB, 0.72_EB,&
1.00_EB, 0.81_EB, 0.75_EB/),(/3,8/))
SD3(1:3,17:24) = RESHAPE ((/&
1.03_EB, 0.84_EB, 0.78_EB,&
1.03_EB, 0.84_EB, 0.78_EB,&
1.00_EB, 0.81_EB, 0.75_EB,&
0.94_EB, 0.77_EB, 0.72_EB,&
0.72_EB, 0.68_EB, 0.68_EB,&
0.52_EB, 0.63_EB, 0.63_EB,&
0.33_EB, 0.50_EB, 0.56_EB,&
0.25_EB, 0.42_EB, 0.50_EB/),(/3,8/))
SD3(1:3,25:32) = RESHAPE ((/&
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
SUBROUTINE RCDEALLOC
!==============================================================================   
DEALLOCATE(P)
DEALLOCATE(SPECIE)
DEALLOCATE(QW)
DEALLOCATE(TTAU)
DEALLOCATE(XT)
DEALLOCATE(X)
DEALLOCATE(GC)
DEALLOCATE(AMBDA)
DEALLOCATE(UUU)
DEALLOCATE(AB)
DEALLOCATE(ATOT)
DEALLOCATE(BCNT)
DEALLOCATE(GAMMA)
DEALLOCATE(SD15)
DEALLOCATE(SD)
DEALLOCATE(SD7)
DEALLOCATE(SD3)

!-------------------------DEALLOCATION Methane VARIABLES-------------------
DEALLOCATE(SD_CH4_TEMP)
DEALLOCATE(OM_BND_CH4)
DEALLOCATE(SD1_CH4)
DEALLOCATE(SD2_CH4)

!-------------------------DEALLOCATION Propane VARIABLES-------------------
DEALLOCATE(SD_C3H8_TEMP)
DEALLOCATE(OM_BND_C3H8)
DEALLOCATE(SD1_C3H8)
DEALLOCATE(SD2_C3H8)

!-------------------------DEALLOCATION Heptane VARIABLES-------------------
DEALLOCATE(SD_C7H16_TEMP)
DEALLOCATE(OM_BND_C7H16)
DEALLOCATE(SD1_C7H16)
DEALLOCATE(SD2_C7H16)

!-------------------------DEALLOCATION Methanol VARIABLES-------------------
DEALLOCATE(SD_CH3OH_TEMP)
DEALLOCATE(OM_BND_CH3OH)
DEALLOCATE(SD1_CH3OH)
DEALLOCATE(SD2_CH3OH)
DEALLOCATE(SD3_CH3OH)

!-------------------------DEALLOCATION Toluene VARIABLES-------------------
DEALLOCATE(SD_C7H8_TEMP)
DEALLOCATE(OM_BND_C7H8)
DEALLOCATE(SD1_C7H8)
DEALLOCATE(SD2_C7H8)
DEALLOCATE(SD3_C7H8)
DEALLOCATE(SD4_C7H8)

!-------------------------DEALLOCATION Propylene VARIABLES-------------------
DEALLOCATE(SD_C3H6_TEMP)
DEALLOCATE(OM_BND_C3H6)
DEALLOCATE(SD1_C3H6)
DEALLOCATE(SD2_C3H6)
DEALLOCATE(SD3_C3H6)

!-------------------------DEALLOCATION MMA VARIABLES-------------------
DEALLOCATE(SD_MMA_TEMP)
DEALLOCATE(OM_BND_MMA)
DEALLOCATE(SD1_MMA)
DEALLOCATE(SD2_MMA)

!------------------------------------------------------------------------------   
END SUBROUTINE RCDEALLOC

!------------------------------------------------------------------------------
END MODULE RADCALV



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
USE RADCALV, ONLY: PLANCK
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
 
! Physical parameters
 
RMMIN = 0.5E-6_EB   ! minimum mean radius (m)
B_WIEN = 2.8977685E-3_EB

NSB = NUMBER_SPECTRAL_BANDS
 
! Find the maximum mean PARTICLE radius
 
RMMAX = 0.0_EB
LPC => LAGRANGIAN_PARTICLE_CLASS(CLASS_NUMBER) 
IF (LPC%DIAMETER>2._EB*RMMAX) RMMAX = 0.5_EB*LPC%DIAMETER
 
! Allow increase of the mean radius
RMMAX = 1.5_EB*RMMAX
 
! Calculate parameters of the PARTICLE group lookup table
 
DGROUP_A = (LOG(RMMAX)-LOG(RMMIN))/(NDG-1)
DGROUP_B = LOG(RMMAX)-DGROUP_A*NDG

! Generate the PARTICLE radii for mie table (microns)
 
RDTMP = 0.2_EB
NX = 0
DO WHILE (RDTMP < 3._EB*RMMAX*1.0E6_EB) 
   NX = NX + 1
   RDTMP = RDTMP + MIN(3000._EB,0.2_EB*RDTMP**(1._EB))
ENDDO
NRDMIE = NX 

ALLOCATE(RDMIE(1:NRDMIE),STAT=IZERO)
CALL ChkMemErr('MIEV','RDMIE',IZERO)

RDTMP = 0.2_EB
RDMIE(1) = RDTMP
DO NX = 2, NRDMIE
   RDTMP = RDTMP + MIN(3000._EB,0.2_EB*RDTMP**(1.0_EB))
   RDMIE(NX) = RDTMP
ENDDO
RDMIE = RDMIE*1.0E-6_EB

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

LMBDMIE(1:NLMBDMIE)          = CPLXREF_WATER(1:NLMBDMIE,1)
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


   DRGROUPLOOP: DO ND = 1, NDG

      LPC%R50(ND) = EXP(DGROUP_A*REAL(ND,EB) + DGROUP_B)

      !     Loop over wavelengths

      IBSUM = 0._EB

      DO J = NLAMBDALOW(1),NLAMBDAHIGH(1)
         IB = PLANCK(RADTMP, LMBDMIE(J)*1.0E6_EB)
         IBSUM = IBSUM + LMBDWGHT(J)*IB

         ASUM = 0._EB
         BSUM = 0._EB

         !     Loop over PARTICLE size distribution

!         DO I = 0,NRDINT
!
!            !     Integrate effective scattering cross section 
!            !     = scattering cross section * (1-forward fraction)
!
!            CALL INTERPOLATE1D(RDMIE,QSCA(:,J),RDDIST(I),AVAL)
!            CALL INTERPOLATE1D(RDMIE,CHI_F(:,J),RDDIST(I),BVAL)
!            BVAL = (1._EB-BVAL)
!            AVAL = AVAL*BVAL*PI*RDDIST(I)**2
!            ASUM = ASUM + RDWGHT(I)*AVAL
!
!            !     Integrate absorption cross sections
!
!            CALL INTERPOLATE1D(RDMIE,QABS(:,J),RDDIST(I),BVAL)
!            BVAL = BVAL*PI*RDDIST(I)**2
!            BSUM = BSUM + RDWGHT(I)*BVAL
!         ENDDO
  
         ! Properties simply at d32 - Elizabeth's idea

         CALL INTERPOLATE1D(RDMIE,QSCA(:,J),LPC%R50(ND),AVAL)
         CALL INTERPOLATE1D(RDMIE,CHI_F(:,J),LPC%R50(ND),BVAL)
         BVAL = (1._EB-BVAL)
         ASUM = AVAL*BVAL
         CALL INTERPOLATE1D(RDMIE,QABS(:,J),LPC%R50(ND),BVAL)
         BSUM = BVAL
         ! End Elizabeth's 

         LPC%WQSCA(ND,IBND) = LPC%WQSCA(ND,IBND) + ASUM*LMBDWGHT(J)*IB
         LPC%WQABS(ND,IBND) = LPC%WQABS(ND,IBND) + BSUM*LMBDWGHT(J)*IB
      ENDDO

      !     Normalize with blackbody radiation

      LPC%WQSCA(ND,IBND)  = LPC%WQSCA(ND,IBND)/IBSUM
      LPC%WQABS(ND,IBND)  = LPC%WQABS(ND,IBND)/IBSUM

      !     Transform cross sections back to efficiency factors

!      LPC%WQSCA(ND,IBND)  = LPC%WQSCA(ND,IBND)/(PI*LPC%R50(ND)**2)
!      LPC%WQABS(ND,IBND)  = LPC%WQABS(ND,IBND)/(PI*LPC%R50(ND)**2)
!     For d32-based properties, no need to divide by drop area
!      LPC%WQSCA(ND,IBND)  = LPC%WQSCA(ND,IBND)/(PI*(LPC%DMN(ND)/2._EB)**2)
!      LPC%WQABS(ND,IBND)  = LPC%WQABS(ND,IBND)/(PI*(LPC%DMN(ND)/2._EB)**2)

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

      IF (ABS(QS)>ZERO_P) THEN
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
            IF (ABS(MUD0(I,J)-MUDPI(I,J))<=ZERO_P) THEN
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
