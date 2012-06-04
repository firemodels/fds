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
USE GLOBAL_CONSTANTS, ONLY: AL2O3, RADCAL_FUEL
IMPLICIT NONE

PRIVATE
CHARACTER(255), PARAMETER :: iradid='$Id$'
CHARACTER(255), PARAMETER :: iradrev='$Revision$'
CHARACTER(255), PARAMETER :: iraddate='$Date$'

PUBLIC OMMAX, OMMIN, DD, SPECIE, SVF, PLANCK, P, RCT, RCALLOC,                 &
       INIT_RADCAL, RADCAL, RCDEALLOC, GET_REV_irad

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: GAMMA, SD15, SD, SD7, SD3,             &
                                         SD3_CH4_NEW, SD7_CH4_NEW, OM_BND_CH4,  & ! FOR METHANE NEW DATA
                                         SD3_C3H8,    SD7_C3H8,    OM_BND_C3H8, & ! FOR PROPANE NEW DATA
                                         SD3_C7H16,   SD7_C7H16,   OM_BND_C7H16,& ! FOR HEPTANE NEW DATA
                                         SD3A_CH3OH,  SD3B_CH3OH,  OM_BND_CH3OH,& ! FOR METHANOL NEW DATA
                                         SD8_CH3OH,   SD10_CH3OH 

REAL(EB), ALLOCATABLE, DIMENSION(:)   :: SPECIE, QW, TTAU, XT, AB,             &
                                         AMBDA, ATOT, BCNT, P, UUU, GC, X

REAL(EB), ALLOCATABLE, DIMENSION(:)   :: SD_CH4_TEMP, SD_C3H8_TEMP, SD_C7H16_TEMP, &
                                         SD_CH3OH_TEMP

REAL(EB) :: OMMIN, OMMAX, TWALL, RCT, AC, AD, DD, XPART, TAU,                  &
            SVF, TAUS, XTOT, XSTAR

INTEGER :: N_TEMP_CH4,   N_BAND_CH4,   &
           N_TEMP_C3H8,  N_BAND_C3H8,  &
           N_TEMP_C7H16, N_BAND_C7H16, & 
           N_TEMP_CH3OH, N_BAND_CH3OH

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
SUBROUTINE RADCAL(AMEAN,AP0)
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
REAL (EB) :: AMEAN, AP0

! LOCAL VARIABLES

REAL(EB) :: DOM, ABGAS, PTOT, TEMP, UK, XD, YD, XX, ENN, ARG, ARGNEW, RSL, RSS,  &
            ABLONG, ABSHRT, ABIL, ABIS, OMEGA, WL, DAMBDA, SDWEAK, GDINV,        &
            GDDINV, YC, Y, AIWALL, XC, AOM, Q, LTERM, AZORCT, RCT4

INTEGER  :: I, II, KK, NM, N, MM, KMAX, KMIN

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
               IF (RADCAL_FUEL=='METHANE') THEN 
                  CALL CH4_NEW(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
               ELSE IF (RADCAL_FUEL=='PROPANE') THEN
                  CALL C3H8(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
               ELSE IF (RADCAL_FUEL=='N-HEPTANE') THEN
                  CALL C7H16(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
               ELSE IF (RADCAL_FUEL=='METHANOL') THEN
                  CALL CH3OH(OMEGA,TEMP,P(3),PTOT,GC(3),SDWEAK,GDINV,GDDINV)
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
            IF((I==3).AND.(XC<=10)) THEN
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
            X(I)=XSTAR*((1._EB-(Y**(-.5_EB)))**.5_EB)
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
               D=2._EB*(GAM*GAM-DELTA)**.5_EB
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
               SDSTRG=(.5_EB*G)**.5_EB*(SQRT(SMINUS)+SQRT(SPLUS))/D+SDSTRG
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
         D=2._EB*(GAM*GAM-DELTA)**.5_EB
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
               D=2._EB*(GAM*GAM-DELTA)**.5_EB
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
               SDSTRG=(.5_EB*G)**.5_EB*(SMINUS**.5+SPLUS**.5)/D+SDSTRG
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
               D=2._EB*(GAM*GAM-DELTA)**.5_EB
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
               SDSTRG=(.5_EB*G)**.5_EB*(SMINUS**.5+SPLUS**.5)/D+SDSTRG
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
   GD   = 5.94E-6_EB*OMEGA*(TEMP/(273._EB*WM))**.5_EB
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
         SDWEAK=SDWEAK+2._EB*(OMEGA-BCNT(I))**2*(-Q2OT*BE)**1.5_EB  *ATOT(I)/SQRTPI*DINV**3*EXP(Q2OT*BE*DINV**2 &
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
REAL(EB), PARAMETER :: WM_CH4 = 16._EB

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
! COMPUTED PROPERTIES 2.4 MICRON BAND
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
   SDWEAK=SDWEAK*2._EB*(-Q2OT*BE_CH4)**1.5_EB/SQRTPI*DINV_CH4**3

!***EXPRESS S/D AT STANDARD TEMPERATURE AND PRESSURE, AS IS IN NASA SP-3080
   SDWEAK = SDWEAK*TOAZ
   GDINV  = GC3*DINV_CH4
   GDDINV = GD*DINV_CH4 

!------------------------------------------------------------------------------
! TABULATED PROPERTIES
ELSE IF((OM_BND_CH4(2,1)<=OMEGA).AND.(OMEGA<OM_BND_CH4(2,2))) THEN
      
!------------------------------------------------------------------------------
! CONTRIBUTION TO 3.3 MICRON BAND FROM
! (0000)-(0010) TRANSITION (V3 FUNDAMENTAL)
! 9.4 IS THE AVERAGE LINE POSITION FOR THE 3.3 MICRON BAND (CM^-1)
! SOURCE: BROSMER & TIEN, JQSRT V. 33, P. 525 (SPECIFIC TO METHANE)

   DINV_CH4 = 1._EB/9.4_EB

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH4_TEMP,OM_BND_CH4(2,:),SD3_CH4_NEW)
! LINE SHAPE PARAMETER BROSMER & TIEN, JQSRT V. 33, P. 525 EQ. 10 (SPECIFIC TO METHANE)
   GDINV  = .00734_EB*PRESSURE_EFFECTIVE*SQRT(AZOT)*EXP(1.02_EB*(TOAZ-1._EB)) 
   GDDINV = GD*DINV_CH4

ELSE IF((OM_BND_CH4(1,1)<=OMEGA).AND.(OMEGA<OM_BND_CH4(1,2))) THEN
!------------------------------------------------------------------------------
! CONTRIBUTION TO 7.7 MICRON BAND FROM
! (0000)-(0001) TRANSITION (V4 FUNDAMENTAL)
! 5.1 (CM-1) IS THE AVERAGE LINE POSITION FOR THE 7.7 MICRON BAND 
! SOURCE: BROSMER & TIEN, JQSRT V. 33, P. 525 (SPECIFIC TO METHANE)

   DINV_CH4 = 1._EB/5.1_EB

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH4_TEMP,OM_BND_CH4(1,:),SD7_CH4_NEW)
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
! FIRST BAND: 7 MICRON BAND (1250 CM^-1 TO 1600 CM^-1)   
   
   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C3H8_TEMP,OM_BND_C3H8(1,:),SD7_C3H8)

!! LINE SHAPE PARAMETER BROSMER & TIEN, JQSRT V. 33, 525 EQ. 11 (SPECIFIC TO METHANE)
!  GDINV  = .0243_EB*PRESSURE_EFFECTIVE*(TOAZ)**.8_EB

   GDINV  = GC3*DINV_C3H8
   GDDINV = GD*DINV_C3H8

ELSE IF ((OM_BND_C3H8(2,1)<=OMEGA).AND.(OMEGA<OM_BND_C3H8(2,2))) THEN
! SECOND BAND: 3.3 MICRON BAND (2600 CM^-1 TO 3350 CM^-1)
   
   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C3H8_TEMP,OM_BND_C3H8(2,:),SD3_C3H8)

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
! FIRST BAND: 7 MICRON BAND (1250 CM^-1 TO 1600 CM^-1)   

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H16_TEMP,OM_BND_C7H16(1,:),SD7_C7H16)
   GDINV  = GC3*DINV_C7H16
   GDDINV = GD*DINV_C7H16

ELSE IF ((OM_BND_C7H16(2,1)<=OMEGA).AND.(OMEGA<OM_BND_C7H16(2,2))) THEN
! SECOND BAND: 3.3 MICRON BAND (2600 CM^-1 TO 3350 CM^-1)

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_C7H16_TEMP,OM_BND_C7H16(2,:),SD3_C7H16)
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
! ABSORPTION BANDS: 3.3 MICRON BAND
!                   7 MICRON BAND
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

! COMPUTE PRESSURE EFFECTIVE. FROM BROSMER & TIEN,524 EQ. 7 (SPECIFIC TO METHANE)
! NEED TO FIND METHANOL EFFECTIVE PRESSURE FORMULATION FOR METHANOL
PRESSURE_EFFECTIVE = PTOT+0.3_EB*PCH3OH ! (SPECIFIC TO METHANE)

! IF OMEGA IS WITHIN ABSOPRTION BAND, THEN COMPUTE SDWEAK FROM TABULATED DATA
! METHANOL HAS FOUR BANDS
! 10.0 MICRON BAND:  900-1150 cm-1, DATA IN SD10_CH3OH
! 8.0  MICRON BAND: 1150-1500 cm-1, DATA IN SD8_CH3OH
! 3.5  MICRON BAND: 2600-3300 cm-1, DATA IN SD3A_CH3OH
! 2.9  MICRON BAND: 3600-3800 cm-1, DATA IN SD3B_CH3OH

IF ((OM_BND_CH3OH(1,1)<=OMEGA).AND.(OMEGA<OM_BND_CH3OH(1,2))) THEN
! FIRST BAND: 10 MICRON BAND (900-1150 cm-1)   

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH3OH_TEMP,OM_BND_CH3OH(1,:),SD10_CH3OH)
   GDINV  = GC3*DINV_CH3OH
   GDDINV = GD*DINV_CH3OH

ELSE IF ((OM_BND_CH3OH(2,1)<=OMEGA).AND.(OMEGA<OM_BND_CH3OH(2,2))) THEN
! SECOND BAND: 8.0  MICRON BAND (1150-1500 cm-1)

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH3OH_TEMP,OM_BND_CH3OH(2,:),SD8_CH3OH)
   GDINV  = GC3*DINV_CH3OH
   GDDINV = GD*DINV_CH3OH

ELSE IF ((OM_BND_CH3OH(3,1)<=OMEGA).AND.(OMEGA<OM_BND_CH3OH(3,2))) THEN
! THIRD BAND: 3.5  MICRON BAND (2600-3300 cm-1)

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH3OH_TEMP,OM_BND_CH3OH(3,:),SD3A_CH3OH)
   GDINV  = GC3*DINV_CH3OH
   GDDINV = GD*DINV_CH3OH

ELSE IF ((OM_BND_CH3OH(4,1)<=OMEGA).AND.(OMEGA<OM_BND_CH3OH(4,2))) THEN
! FOURTH BAND: 2.9  MICRON BAND (3600-3800 cm-1)

   SDWEAK = GET_SPECTRAL_ABSORPTION(OMEGA,TEMP,SD_CH3OH_TEMP,OM_BND_CH3OH(4,:),SD3B_CH3OH)
   GDINV  = GC3*DINV_CH3OH
   GDDINV = GD*DINV_CH3OH

ENDIF

!------------------------------------------------------------------------------
END SUBROUTINE CH3OH

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

!------------------------------------------------------------------------------
! NEW METHANE DATA
N_TEMP_CH4 = 7
N_BAND_CH4 = 3

ALLOCATE(SD_CH4_TEMP(N_TEMP_CH4))
ALLOCATE(SD3_CH4_NEW(N_TEMP_CH4,33))
ALLOCATE(SD7_CH4_NEW(N_TEMP_CH4,17))
ALLOCATE(OM_BND_CH4(N_BAND_CH4,3))

! Initialize wavenumber limits and spacing for bands in the new CH4 absorption 
! 1st column is lower band limit.
! 2nd column being upper band limit.
! 3rd column being spacing between frequencies in band.  If 3rd column = 0., the band is calculated and not tabulated. 

! OM_BND_CH4 = RESHAPE((/ &
!     1100._EB,1500._EB,25._EB, &
!     2600._EB,3400._EB,25._EB, &
!     3400._EB,5000._EB,0._EB /),(/N_BAND_CH4,3/))

OM_BND_CH4 = RESHAPE((/ &
    1100._EB,2600._EB,3400._EB, &
    1500._EB,3400._EB,5000._EB, &
    25._EB  ,25._EB  ,0._EB /),(/N_BAND_CH4,3/))

SD_CH4_TEMP = (/300._EB, 400._EB, 450._EB, 500._EB, 600._EB, 800._EB, 1000._EB/)
! TEMP,K= 300     400      450      500      600      800      1000   Verify temperatures   

! Initialize SD3_CH4_NEW array
! THIS IS THE NEW CH4 DATA FOR 3.3 MICRON BAND
SD3_CH4_NEW(1:7,1:8) = RESHAPE ((/ &  ! 2600-2800 cm-1
   1.834E-03_EB, 8.687E-05_EB, 2.770E-03_EB, 1.697E-03_EB, 3.182E-04_EB, 7.942E-03_EB, 1.256E-02_EB,&
   3.084E-03_EB, 4.315E-03_EB, 1.037E-02_EB, 1.016E-02_EB, 5.816E-03_EB, 1.214E-02_EB, 1.746E-02_EB,&
   1.393E-03_EB, 2.722E-03_EB, 5.856E-03_EB, 6.409E-03_EB, 3.114E-03_EB, 1.141E-02_EB, 1.364E-02_EB,&
   5.957E-03_EB, 6.076E-03_EB, 1.116E-02_EB, 9.532E-03_EB, 4.563E-03_EB, 1.232E-02_EB, 1.453E-02_EB,&
   5.812E-03_EB, 4.492E-03_EB, 9.593E-03_EB, 6.866E-03_EB, 4.168E-03_EB, 1.206E-02_EB, 2.303E-02_EB,&
   1.833E-02_EB, 1.329E-02_EB, 1.730E-02_EB, 1.223E-02_EB, 8.811E-03_EB, 2.005E-02_EB, 3.239E-02_EB,&
   1.901E-02_EB, 1.550E-02_EB, 2.191E-02_EB, 1.824E-02_EB, 1.509E-02_EB, 2.740E-02_EB, 4.528E-02_EB,&
   2.541E-02_EB, 2.239E-02_EB, 3.033E-02_EB, 2.873E-02_EB, 2.535E-02_EB, 4.444E-02_EB, 6.472E-02_EB/),(/7,8/))
SD3_CH4_NEW(1:7,9:16) = RESHAPE ((/ &  ! 2800-3000 cm-1
   4.356E-02_EB, 3.951E-02_EB, 4.510E-02_EB, 4.426E-02_EB, 4.862E-02_EB, 7.205E-02_EB, 9.634E-02_EB,&
   2.998E-02_EB, 4.170E-02_EB, 5.256E-02_EB, 5.541E-02_EB, 6.694E-02_EB, 9.550E-02_EB, 1.224E-01_EB,&
   6.730E-02_EB, 1.122E-01_EB, 1.281E-01_EB, 1.339E-01_EB, 1.453E-01_EB, 1.622E-01_EB, 1.680E-01_EB,&
   1.555E-01_EB, 2.090E-01_EB, 2.114E-01_EB, 2.127E-01_EB, 2.064E-01_EB, 1.912E-01_EB, 1.853E-01_EB,&
   2.802E-01_EB, 2.943E-01_EB, 2.744E-01_EB, 2.607E-01_EB, 2.356E-01_EB, 1.965E-01_EB, 1.824E-01_EB,&
   5.762E-01_EB, 4.872E-01_EB, 4.239E-01_EB, 3.820E-01_EB, 3.138E-01_EB, 2.227E-01_EB, 1.909E-01_EB,&
   4.079E-01_EB, 3.157E-01_EB, 2.633E-01_EB, 2.331E-01_EB, 1.848E-01_EB, 1.367E-01_EB, 1.575E-01_EB,&
   3.553E-01_EB, 2.608E-01_EB, 2.312E-01_EB, 2.182E-01_EB, 2.250E-01_EB, 3.165E-01_EB, 4.357E-01_EB/),(/7,8/))
SD3_CH4_NEW(1:7,17:24) = RESHAPE ((/ &  ! 3000-3200 cm-1
   2.174E+00_EB, 2.001E+00_EB, 1.807E+00_EB, 1.684E+00_EB, 1.462E+00_EB, 1.083E+00_EB, 8.224E-01_EB,&
   2.236E-01_EB, 1.481E-01_EB, 1.202E-01_EB, 1.013E-01_EB, 7.980E-02_EB, 6.123E-02_EB, 7.703E-02_EB,&
   3.869E-01_EB, 2.746E-01_EB, 2.386E-01_EB, 2.147E-01_EB, 1.830E-01_EB, 1.531E-01_EB, 1.482E-01_EB,&
   5.404E-01_EB, 4.491E-01_EB, 4.020E-01_EB, 3.779E-01_EB, 3.332E-01_EB, 2.677E-01_EB, 2.442E-01_EB,&
   4.780E-01_EB, 4.504E-01_EB, 4.257E-01_EB, 4.040E-01_EB, 3.718E-01_EB, 3.222E-01_EB, 2.969E-01_EB,&
   2.537E-01_EB, 3.193E-01_EB, 3.247E-01_EB, 3.248E-01_EB, 3.197E-01_EB, 3.055E-01_EB, 2.960E-01_EB,&
   7.367E-02_EB, 1.464E-01_EB, 1.662E-01_EB, 1.781E-01_EB, 1.948E-01_EB, 2.213E-01_EB, 2.365E-01_EB,&
   1.035E-02_EB, 4.986E-02_EB, 6.685E-02_EB, 8.169E-02_EB, 1.011E-01_EB, 1.405E-01_EB, 1.757E-01_EB/),(/7,8/))
SD3_CH4_NEW(1:7,25:33) = RESHAPE ((/ &  ! 3200-3400 cm-1
   1.817E-03_EB, 1.287E-02_EB, 2.347E-02_EB, 2.940E-02_EB, 4.057E-02_EB, 7.610E-02_EB, 1.049E-01_EB,&
   1.248E-03_EB, 2.738E-03_EB, 1.478E-02_EB, 1.504E-02_EB, 1.683E-02_EB, 4.041E-02_EB, 5.690E-02_EB,&
   2.256E-03_EB, 1.031E-03_EB, 1.133E-02_EB, 7.837E-03_EB, 6.460E-03_EB, 1.969E-02_EB, 3.311E-02_EB,&
   1.852E-03_EB, 4.849E-04_EB, 8.952E-03_EB, 4.252E-03_EB, 2.205E-03_EB, 1.162E-02_EB, 2.115E-02_EB,&
   2.394E-03_EB, 7.241E-04_EB, 1.055E-02_EB, 3.946E-03_EB, 1.096E-03_EB, 1.046E-02_EB, 1.345E-02_EB,&
   1.856E-03_EB, 1.064E-03_EB, 8.816E-03_EB, 3.652E-03_EB, 4.306E-04_EB, 9.953E-03_EB, 7.618E-03_EB,&
   2.841E-03_EB, 2.559E-03_EB, 7.556E-03_EB, 4.026E-03_EB, 4.953E-04_EB, 8.845E-03_EB, 1.265E-02_EB,&
   3.600E-03_EB, 1.078E-03_EB, 8.019E-03_EB, 5.251E-03_EB, 1.139E-03_EB, 1.083E-02_EB, 8.573E-03_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB/),(/7,9/))

!Initialize SD7_CH4_NEW Array
! DATA FOR CH4 FOR THE 7.7 MICRON BAND
! TEMP,K= 300     400      450      500      600      800      1000        

SD7_CH4_NEW(1:7,1:8) = RESHAPE ((/ &  ! 1100-1275 cm-1
   0.000E+000_EB, 0.000E+000_EB, 3.832E-004_EB, 3.819E-003_EB, 3.633E-003_EB, 1.660E-002_EB, 3.956E-002_EB,&  !1100
   0.000E+000_EB, 0.000E+000_EB, 7.972E-004_EB, 6.138E-003_EB, 8.512E-003_EB, 3.221E-002_EB, 5.657E-002_EB,&  !1125
   1.826E-005_EB, 4.211E-004_EB, 7.895E-003_EB, 1.501E-002_EB, 2.355E-002_EB, 5.485E-002_EB, 8.598E-002_EB,&  !1150
   5.287E-003_EB, 1.833E-002_EB, 3.369E-002_EB, 4.756E-002_EB, 6.345E-002_EB, 1.004E-001_EB, 1.267E-001_EB,&  !1175
   5.990E-002_EB, 1.034E-001_EB, 1.200E-001_EB, 1.330E-001_EB, 1.478E-001_EB, 1.601E-001_EB, 1.547E-001_EB,&  !1200
   2.616E-001_EB, 2.743E-001_EB, 2.622E-001_EB, 2.553E-001_EB, 2.335E-001_EB, 2.067E-001_EB, 1.931E-001_EB,&  !1225
   5.434E-001_EB, 4.105E-001_EB, 3.412E-001_EB, 3.128E-001_EB, 2.649E-001_EB, 2.339E-001_EB, 2.248E-001_EB,&  !1250
   6.659E-001_EB, 5.648E-001_EB, 5.112E-001_EB, 4.982E-001_EB, 4.777E-001_EB, 4.637E-001_EB, 4.419E-001_EB/),(/7,8/))
SD7_CH4_NEW(1:7,9:17) = RESHAPE ((/ &  ! 1300-1500 cm-1
   1.137E+000_EB, 8.575E-001_EB, 7.112E-001_EB, 6.389E-001_EB, 5.057E-001_EB, 3.361E-001_EB, 2.645E-001_EB,& !1300
   8.309E-001_EB, 5.892E-001_EB, 4.893E-001_EB, 4.353E-001_EB, 3.536E-001_EB, 2.720E-001_EB, 2.415E-001_EB,& !1325
   4.924E-001_EB, 5.079E-001_EB, 4.564E-001_EB, 4.432E-001_EB, 4.013E-001_EB, 3.322E-001_EB, 3.012E-001_EB,& !1350
   3.538E-002_EB, 9.188E-002_EB, 1.045E-001_EB, 1.296E-001_EB, 1.444E-001_EB, 1.661E-001_EB, 1.873E-001_EB,& !1375
   1.597E-003_EB, 5.040E-003_EB, 1.315E-002_EB, 2.384E-002_EB, 3.093E-002_EB, 5.027E-002_EB, 7.346E-002_EB,& !1400
   6.513E-003_EB, 6.616E-003_EB, 1.117E-002_EB, 1.500E-002_EB, 1.240E-002_EB, 1.911E-002_EB, 2.972E-002_EB,& !1425
   8.471E-004_EB, 7.246E-004_EB, 1.479E-003_EB, 1.166E-002_EB, 3.081E-003_EB, 9.764E-003_EB, 2.664E-002_EB,& !1450
   2.373E-004_EB, 0.000E+000_EB, 2.518E-003_EB, 9.771E-003_EB, 3.017E-003_EB, 8.824E-003_EB, 1.870E-002_EB,& !1475
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB/),(/7,9/))

!---------------------------------PROPANE DATA-------------------------------------------

N_TEMP_C3H8 = 6
ALLOCATE(SD_C3H8_TEMP(N_TEMP_C3H8))   
SD_C3H8_TEMP = (/295._EB, 396._EB, 435._EB, 513._EB, 790._EB, 1009._EB/)
! TEMP CHECK OKAY.

! PROPANE DATA FOR 2 BANDS: 3.3 MICRON BAND AND 7 MICRONS
N_BAND_C3H8 = 2 
ALLOCATE(OM_BND_C3H8(N_BAND_C3H8,3))
OM_BND_C3H8 = RESHAPE((/ &
    1250._EB,2600._EB, &
    1600._EB,3350._EB, &
    25._EB  ,  25._EB /),(/N_BAND_C3H8,3/))

! 1ST COLUMN IS TEHE LOWER BAND LIMIT IN CM^-1
! 2ND COLUMN IS THE UPPER BAND LIMIT IN CM^-1
! 3RD COLUMN IS THE SPACING BETWEEN MEASUREMENT WAVENUMBERS. IF EQUALS TO 0, THEN
! THE BAND IS CALCULATED, NOT TABULATED

ALLOCATE(SD7_C3H8(N_TEMP_C3H8,15))
ALLOCATE(SD3_C3H8(N_TEMP_C3H8,31))

! INITIALIZE SD7_C3H8 array
! DATA FOR C3H8 FOR THE 7.7 MICRON BAND
! 6 TEMPERATURES (K): 295, 396, 435, 513, 790, 1009       
! 15 DIFFERENT WAVENUMBERS FOR THIS BAND

SD7_C3H8(1:6,1:8) = RESHAPE ((/ &  ! 1250-1450 CM-1
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.082E-001_EB, 1.442E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.179E-001_EB, 2.036E-001_EB,&
   3.019E-002_EB, 1.529E-001_EB, 8.475E-002_EB, 8.865E-002_EB, 1.695E-001_EB, 2.806E-001_EB,&
   1.042E-001_EB, 2.192E-001_EB, 1.335E-001_EB, 1.013E-001_EB, 2.036E-001_EB, 2.455E-001_EB,&
   3.069E-001_EB, 3.556E-001_EB, 2.682E-001_EB, 2.106E-001_EB, 2.563E-001_EB, 2.871E-001_EB,&
   5.349E-001_EB, 4.394E-001_EB, 3.351E-001_EB, 2.124E-001_EB, 2.748E-001_EB, 2.991E-001_EB,&
   2.163E-001_EB, 3.108E-001_EB, 2.484E-001_EB, 2.670E-001_EB, 3.040E-001_EB, 2.689E-001_EB,&
   4.949E-001_EB, 4.872E-001_EB, 3.975E-001_EB, 3.264E-001_EB, 3.011E-001_EB, 2.894E-001_EB/),(/6,8/))

SD7_C3H8(1:6,9:15) = RESHAPE ((/ &  ! 1450-1600 CM-1
   1.079E+000_EB, 7.161E-001_EB, 5.953E-001_EB, 4.989E-001_EB, 3.400E-001_EB, 2.718E-001_EB,&
   7.405E-001_EB, 5.134E-001_EB, 4.199E-001_EB, 3.283E-001_EB, 2.718E-001_EB, 2.327E-001_EB,&
   2.758E-001_EB, 2.173E-001_EB, 1.870E-001_EB, 2.524E-001_EB, 1.978E-001_EB, 1.861E-001_EB,&
   0.000E+000_EB, 5.167E-002_EB, 2.682E-002_EB, 1.212E-001_EB, 1.470E-001_EB, 1.539E-001_EB,&
   0.000E+000_EB, 1.754E-002_EB, 0.000E+000_EB, 5.852E-002_EB, 1.149E-001_EB, 1.568E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.366E-002_EB, 1.111E-001_EB, 1.198E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB/),(/6,7/))

! INITIALIZE SD3_C3H8 ARRAY
! DATA FOR C3H8 FOR THE 3.3 MICRON BAND
! 6 TEMPERATURES (K):  295, 396, 435, 513, 790, 1009
! 31 DIFFERENT WAVENUMBERS FOR THIS BAND

SD3_C3H8(1:6,1:8) = RESHAPE ((/ &  ! 2600-2800 cm-1
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.169E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.422E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 8.964E-002_EB, 0.000E+000_EB, 2.027E-001_EB, 1.685E-001_EB,&
   0.000E+000_EB, 2.243E-002_EB, 1.013E-001_EB, 0.000E+000_EB, 2.582E-001_EB, 1.861E-001_EB,&
   1.588E-002_EB, 1.013E-001_EB, 1.093E-001_EB, 8.050E-002_EB, 3.040E-001_EB, 2.225E-001_EB,&
   6.236E-002_EB, 1.179E-001_EB, 1.881E-001_EB, 1.013E-001_EB, 3.186E-001_EB, 2.270E-001_EB,&
   9.257E-002_EB, 1.812E-001_EB, 2.027E-001_EB, 1.013E-001_EB, 3.917E-001_EB, 3.166E-001_EB/),(/6,8/))

SD3_C3H8(1:6,9:16) = RESHAPE ((/ &  ! 2800-3000 cm-1
   1.160E-001_EB, 2.251E-001_EB, 2.534E-001_EB, 1.472E-001_EB, 4.852E-001_EB, 4.189E-001_EB,&
   4.093E-001_EB, 4.794E-001_EB, 5.145E-001_EB, 3.800E-001_EB, 7.619E-001_EB, 6.703E-001_EB,&
   2.245E+000_EB, 1.799E+000_EB, 1.782E+000_EB, 1.275E+000_EB, 1.330E+000_EB, 1.026E+000_EB,&
   4.231E+000_EB, 2.891E+000_EB, 2.724E+000_EB, 1.896E+000_EB, 1.794E+000_EB, 1.366E+000_EB,&
   4.850E+000_EB, 3.813E+000_EB, 3.685E+000_EB, 2.719E+000_EB, 2.523E+000_EB, 1.824E+000_EB,&
   6.917E+000_EB, 5.300E+000_EB, 5.105E+000_EB, 3.714E+000_EB, 3.090E+000_EB, 2.144E+000_EB,&
   1.330E+001_EB, 8.452E+000_EB, 7.668E+000_EB, 5.141E+000_EB, 3.477E+000_EB, 2.236E+000_EB,&
   9.395E+000_EB, 6.221E+000_EB, 5.621E+000_EB, 3.796E+000_EB, 2.666E+000_EB, 1.730E+000_EB/),(/6,8/))

SD3_C3H8(1:6,17:24) = RESHAPE ((/ &  ! 3000-3200 cm-1
   1.985E+000_EB, 2.010E+000_EB, 1.930E+000_EB, 1.489E+000_EB, 1.495E+000_EB, 1.173E+000_EB,&
   8.176E-002_EB, 3.438E-001_EB, 4.091E-001_EB, 3.429E-001_EB, 7.628E-001_EB, 6.469E-001_EB,&
   0.000E+000_EB, 1.013E-001_EB, 1.608E-001_EB, 1.013E-001_EB, 4.628E-001_EB, 4.521E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 1.013E-001_EB, 1.013E-001_EB, 4.053E-001_EB, 3.855E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 1.033E-001_EB, 9.938E-002_EB, 3.956E-001_EB, 3.322E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 1.325E-001_EB, 7.305E-002_EB, 3.283E-001_EB, 2.836E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 1.706E-001_EB, 2.439E-002_EB, 3.040E-001_EB, 2.397E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 1.413E-001_EB, 1.071E-002_EB, 3.040E-001_EB, 2.153E-001_EB/),(/6,8/))

SD3_C3H8(1:6,25:31) = RESHAPE ((/ &  ! 3200-3350 cm-1
   0.000E+000_EB, 0.000E+000_EB, 1.013E-001_EB, 0.000E+000_EB, 2.981E-001_EB, 1.938E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 2.553E-001_EB, 1.520E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.729E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.530E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.403E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 1.345E-001_EB,&
   0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB, 0.000E+000_EB/),(/6,7/))

!---------------------------------HEPTANE DATA-------------------------------------------

N_TEMP_C7H16 = 6
ALLOCATE(SD_C7H16_TEMP(N_TEMP_C7H16))   

SD_C7H16_TEMP = (/293._EB, 400._EB, 450._EB, 490._EB, 593._EB, 794._EB/)

! HEPTANE DATA FOR 2 BANDS: 3.3 MICRON BAND AND 7 MICRONS
N_BAND_C7H16 = 2 

!--TO DO: CJHECK BOUNDS OF WAVENUMBER

ALLOCATE(OM_BND_C7H16(N_BAND_C7H16,3))

OM_BND_C7H16 = RESHAPE((/ &
    1300._EB,2600._EB,&
    1600._EB,3350._EB,&
      25._EB,  25._EB/),(/N_BAND_C7H16,3/))

! 1ST COLUMN IS TEHE LOWER BAND LIMIT IN CM^-1
! 2ND COLUMN IS THE UPPER BAND LIMIT IN CM^-1
! 3RD COLUMN IS THE SPACING BETWEEN MEASUREMENT WAVENUMBERS. IF EQUALS TO 0, THEN
! THE BAND IS CALCULATED, NOT TABULATED

ALLOCATE(SD7_C7H16(N_TEMP_C7H16,13))
ALLOCATE(SD3_C7H16(N_TEMP_C7H16,31))

! INITIALIZE SD7_C7H16 array
! DATA FOR C7H16 FOR THE 7.7 MICRON BAND
! 6 TEMPERATURES (K): 293, 400, 450, 490, 593, 794       
! 13 DIFFERENT WAVENUMBERS FOR THIS BAND

SD7_C7H16(1:6,1:8) = RESHAPE ((/ &  ! 1300-1475 cm-1
   1.059E-01_EB, 1.336E-01_EB, 1.436E-01_EB, 1.516E-01_EB, 1.627E-01_EB, 1.717E-01_EB,&
   2.201E-01_EB, 2.509E-01_EB, 2.590E-01_EB, 2.638E-01_EB, 2.666E-01_EB, 2.571E-01_EB,&
   7.466E-01_EB, 6.562E-01_EB, 6.217E-01_EB, 5.915E-01_EB, 5.400E-01_EB, 4.593E-01_EB,&
   1.763E+00_EB, 1.259E+00_EB, 1.096E+00_EB, 9.663E-01_EB, 7.740E-01_EB, 5.372E-01_EB,&
   1.911E-01_EB, 2.431E-01_EB, 2.651E-01_EB, 2.845E-01_EB, 3.160E-01_EB, 3.547E-01_EB,&
   5.911E-01_EB, 6.075E-01_EB, 6.136E-01_EB, 6.179E-01_EB, 6.213E-01_EB, 6.121E-01_EB,&
   3.397E+00_EB, 2.294E+00_EB, 1.935E+00_EB, 1.653E+00_EB, 1.246E+00_EB, 7.759E-01_EB,&
   1.314E+00_EB, 9.574E-01_EB, 8.373E-01_EB, 7.415E-01_EB, 5.993E-01_EB, 4.259E-01_EB/),(/6,8/))
SD7_C7H16(1:6,9:13) = RESHAPE ((/ &  ! 1500-1600 cm-1
   1.016E-01_EB, 1.312E-01_EB, 1.439E-01_EB, 1.551E-01_EB, 1.733E-01_EB, 1.959E-01_EB,&
   6.047E-02_EB, 6.721E-02_EB, 7.074E-02_EB, 7.411E-02_EB, 8.014E-02_EB, 8.908E-02_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.00E+000_EB, 0.00E+000_EB, 0.00E+000_EB, 0.00E+000_EB/),(/6,5/))

! INITIALIZE SD3_C7H16 ARRAY
! DATA FOR C7H16 FOR THE 3.3 MICRON BAND
! 6 TEMPERATURES (K):  293, 400, 450, 490, 593, 794       
! 31 DIFFERENT WAVENUMBERS FOR THIS BAND

SD3_C7H16(1:6,1:8) = RESHAPE ((/ &  ! 2600-2800 cm-1
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   3.987E-09_EB, 7.440E-09_EB, 8.769E-09_EB, 9.797E-09_EB, 1.096E-08_EB, 1.116E-08_EB,&
   1.166E-01_EB, 1.270E-01_EB, 1.295E-01_EB, 1.309E-01_EB, 1.315E-01_EB, 1.281E-01_EB,&
   2.772E-01_EB, 2.403E-01_EB, 2.261E-01_EB, 2.139E-01_EB, 1.937E-01_EB, 1.644E-01_EB,&
   1.889E-01_EB, 2.000E-01_EB, 2.036E-01_EB, 2.065E-01_EB, 2.104E-01_EB, 2.140E-01_EB,&
   2.569E-01_EB, 2.811E-01_EB, 2.897E-01_EB, 2.968E-01_EB, 3.075E-01_EB, 3.202E-01_EB/),(/6,8/))
SD3_C7H16(1:6,9:16) = RESHAPE ((/ &  ! 2800-3000 cm-1
   5.248E-01_EB, 5.398E-01_EB, 5.447E-01_EB, 5.485E-01_EB, 5.537E-01_EB, 5.574E-01_EB,&
   1.730E+00_EB, 1.590E+00_EB, 1.538E+00_EB, 1.494E+00_EB, 1.422E+00_EB, 1.315E+00_EB,&
   8.615E+00_EB, 6.276E+00_EB, 5.493E+00_EB, 4.869E+00_EB, 3.940E+00_EB, 2.798E+00_EB,&
   1.156E+01_EB, 8.692E+00_EB, 7.637E+00_EB, 6.767E+00_EB, 5.432E+00_EB, 3.750E+00_EB,&
   1.312E+01_EB, 1.043E+01_EB, 9.466E+00_EB, 8.667E+00_EB, 7.420E+00_EB, 5.759E+00_EB,&
   2.116E+01_EB, 1.597E+01_EB, 1.411E+01_EB, 1.258E+01_EB, 1.024E+01_EB, 7.266E+00_EB,&
   2.663E+01_EB, 1.866E+01_EB, 1.598E+01_EB, 1.385E+01_EB, 1.072E+01_EB, 7.001E+00_EB,&
   8.944E+00_EB, 6.853E+00_EB, 6.066E+00_EB, 5.405E+00_EB, 4.372E+00_EB, 3.034E+00_EB/),(/6,8/))
SD3_C7H16(1:6,17:24) = RESHAPE ((/ &  ! 3000-3200 cm-1
   7.413E-01_EB, 8.941E-01_EB, 9.329E-01_EB, 9.553E-01_EB, 9.680E-01_EB, 9.294E-01_EB,&
   4.448E-01_EB, 4.982E-01_EB, 5.136E-01_EB, 5.243E-01_EB, 5.363E-01_EB, 5.394E-01_EB,&
   2.315E-01_EB, 2.790E-01_EB, 2.969E-01_EB, 3.119E-01_EB, 3.356E-01_EB, 3.665E-01_EB,&
   1.445E-01_EB, 1.890E-01_EB, 2.067E-01_EB, 2.222E-01_EB, 2.475E-01_EB, 2.826E-01_EB,&
   9.034E-02_EB, 1.286E-01_EB, 1.447E-01_EB, 1.589E-01_EB, 1.829E-01_EB, 2.175E-01_EB,&
   7.157E-02_EB, 1.037E-01_EB, 1.173E-01_EB, 1.295E-01_EB, 1.501E-01_EB, 1.801E-01_EB,&
   8.822E-02_EB, 1.060E-01_EB, 1.125E-01_EB, 1.178E-01_EB, 1.259E-01_EB, 1.355E-01_EB,&
   1.072E-01_EB, 1.010E-01_EB, 9.826E-02_EB, 9.572E-02_EB, 9.120E-02_EB, 8.379E-02_EB/),(/6,8/))
SD3_C7H16(1:6,25:31) = RESHAPE ((/ &  ! 3200-3350 cm-1
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB/),(/6,7/))

!---------------------------------METHANOL DATA-------------------------------------------

N_TEMP_CH3OH = 7
ALLOCATE(SD_CH3OH_TEMP(N_TEMP_C7H16))   

SD_CH3OH_TEMP = (/293._EB, 396._EB, 450._EB, 500._EB, 570._EB, 804._EB, 1000._EB/)

N_BAND_CH3OH = 4
ALLOCATE(OM_BND_CH3OH(N_BAND_CH3OH,3))

! Initialize wavenumber limits and spacing for bands in C3H8 absorption 
! 1st column is lower band limit.
! 2nd column being upper band limit.
! 3rd column being spacing between frequencies in band.  If 3rd column = 0., the band is calculated and not tabulated.

OM_BND_CH3OH = RESHAPE((/ &
     900._EB,1150._EB,2600._EB, 3600._EB,&
    1150._EB,1500._EB,3300._EB, 3800._EB,&
      25._EB,  25._EB,  25._EB,   25._EB/),(/N_BAND_CH3OH,3/))

ALLOCATE(SD10_CH3OH(N_TEMP_CH3OH,10)) ! DATA FOR CH3OH FOR THE 10.0 MICRON BAND
ALLOCATE(SD8_CH3OH(N_TEMP_CH3OH,15))  ! DATA FOR CH3OH FOR THE 8.0 MICRON BAND
ALLOCATE(SD3A_CH3OH(N_TEMP_CH3OH,28)) ! DATA FOR CH3OH FOR THE 3.5 MICRON BAND
ALLOCATE(SD3B_CH3OH(N_TEMP_CH3OH,9))  ! DATA FOR CH3OH FOR THE 2.9 MICRON BANDLOOP

! Initialize SD10_CH3OH array
! DATA FOR CH3OH FOR THE 10.0 MICRON BAND
! TEMP,K= 293       396      450      500     570     804       1000

SD10_CH3OH(1:7,1:8) = RESHAPE ((/ &  ! 900-1100 cm-1
   4.647E-003_EB, 3.248E-002_EB, 6.138E-002_EB, 1.014E-001_EB, 2.115E-001_EB, 5.049E-001_EB, 8.128E-001_EB,&
   3.139E-002_EB, 1.313E-001_EB, 2.050E-001_EB, 2.879E-001_EB, 4.612E-001_EB, 7.537E-001_EB, 9.293E-001_EB,&
   3.110E-001_EB, 6.537E-001_EB, 8.136E-001_EB, 9.526E-001_EB, 1.157E+000_EB, 1.315E+000_EB, 1.276E+000_EB,&
   2.584E+000_EB, 2.571E+000_EB, 2.510E+000_EB, 2.427E+000_EB, 2.226E+000_EB, 1.809E+000_EB, 1.456E+000_EB,&
   5.746E+000_EB, 3.915E+000_EB, 3.329E+000_EB, 2.870E+000_EB, 2.202E+000_EB, 1.414E+000_EB, 9.824E-001_EB,&
   5.859E+000_EB, 4.166E+000_EB, 3.614E+000_EB, 3.176E+000_EB, 2.526E+000_EB, 1.732E+000_EB, 1.273E+000_EB,&
   6.713E+000_EB, 5.042E+000_EB, 4.460E+000_EB, 3.979E+000_EB, 3.228E+000_EB, 2.243E+000_EB, 1.639E+000_EB,&
   5.041E-001_EB, 7.345E-001_EB, 8.121E-001_EB, 8.677E-001_EB, 9.270E-001_EB, 9.211E-001_EB, 8.462E-001_EB/),(/7,8/))
SD10_CH3OH(1:7,9:10) = RESHAPE ((/ &  ! 1100-1150 cm-1
   1.571E-001_EB, 1.855E-001_EB, 1.929E-001_EB, 1.971E-001_EB, 1.989E-001_EB, 1.889E-001_EB, 1.727E-001_EB,&
   1.022E-001_EB, 1.335E-001_EB, 1.452E-001_EB, 1.547E-001_EB, 1.684E-001_EB, 1.808E-001_EB, 1.820E-001_EB/),(/7,2/))

! Initialize SD8_CH3OH array
! DATA FOR CH3OH FOR THE 8.0 MICRON BAND
! TEMP,K= 293       396      450      500     570     804       1000

SD8_CH3OH(1:7,1:8) = RESHAPE ((/ &  ! 1150-1350 cm-1
   7.334E-02_EB, 1.218E-01_EB, 1.436E-01_EB, 1.631E-01_EB, 1.955E-01_EB, 2.370E-01_EB, 2.568E-01_EB,&
   1.068E-01_EB, 1.798E-01_EB, 2.121E-01_EB, 2.407E-01_EB, 2.871E-01_EB, 3.431E-01_EB, 3.663E-01_EB,&
   2.022E-01_EB, 2.803E-01_EB, 3.015E-01_EB, 3.132E-01_EB, 3.170E-01_EB, 2.856E-01_EB, 2.411E-01_EB,&
   2.437E-01_EB, 3.175E-01_EB, 3.362E-01_EB, 3.461E-01_EB, 3.478E-01_EB, 3.157E-01_EB, 2.713E-01_EB,&
   3.165E-01_EB, 3.932E-01_EB, 4.087E-01_EB, 4.139E-01_EB, 4.046E-01_EB, 3.520E-01_EB, 2.929E-01_EB,&
   3.985E-01_EB, 4.315E-01_EB, 4.313E-01_EB, 4.244E-01_EB, 3.993E-01_EB, 3.346E-01_EB, 2.745E-01_EB,&
   6.900E-01_EB, 5.333E-01_EB, 4.762E-01_EB, 4.284E-01_EB, 3.528E-01_EB, 2.519E-01_EB, 1.887E-01_EB,&
   8.225E-01_EB, 5.998E-01_EB, 5.255E-01_EB, 4.659E-01_EB, 3.757E-01_EB, 2.619E-01_EB, 1.938E-01_EB/),(/7,8/))
SD8_CH3OH(1:7,9:15) = RESHAPE ((/ &  ! 1350-1500 cm-1
   7.610E-01_EB, 5.776E-01_EB, 5.141E-01_EB, 4.621E-01_EB, 3.817E-01_EB, 2.760E-01_EB, 2.099E-01_EB,&
   5.794E-01_EB, 4.874E-01_EB, 4.526E-01_EB, 4.226E-01_EB, 3.729E-01_EB, 2.995E-01_EB, 2.472E-01_EB,&
   4.243E-01_EB, 3.675E-01_EB, 3.457E-01_EB, 3.266E-01_EB, 2.942E-01_EB, 2.445E-01_EB, 2.075E-01_EB,&
   2.991E-01_EB, 2.528E-01_EB, 2.358E-01_EB, 2.213E-01_EB, 1.977E-01_EB, 1.631E-01_EB, 1.383E-01_EB,&
   3.841E-01_EB, 3.252E-01_EB, 3.034E-01_EB, 2.849E-01_EB, 2.544E-01_EB, 2.096E-01_EB, 1.773E-01_EB,&
   2.736E-01_EB, 2.342E-01_EB, 2.194E-01_EB, 2.066E-01_EB, 1.853E-01_EB, 1.534E-01_EB, 1.300E-01_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB/),(/7,7/))
! Initialize SD3A_CH3OH array
! DATA FOR CH3OH FOR THE 3.5 MICRON BAND
! TEMP,K= 293       396      450      500     570     804       1000 
SD3A_CH3OH(1:7,1:8) = RESHAPE ((/ &  ! 2600-2800 cm-1
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   5.471E-09_EB, 1.639E-08_EB, 2.250E-08_EB, 2.852E-08_EB, 3.861E-08_EB, 4.966E-08_EB, 5.180E-08_EB,&
   6.520E-03_EB, 1.947E-02_EB, 2.684E-02_EB, 3.392E-02_EB, 4.590E-02_EB, 5.913E-02_EB, 6.165E-02_EB,&
   1.557E-02_EB, 4.010E-02_EB, 5.316E-02_EB, 6.538E-02_EB, 8.537E-02_EB, 1.063E-01_EB, 1.092E-01_EB,&
   6.963E-02_EB, 1.171E-01_EB, 1.341E-01_EB, 1.468E-01_EB, 1.612E-01_EB, 1.634E-01_EB, 1.504E-01_EB,&
   2.786E-01_EB, 3.526E-01_EB, 3.694E-01_EB, 3.771E-01_EB, 3.750E-01_EB, 3.393E-01_EB, 2.939E-01_EB/),(/7,8/))
SD3A_CH3OH(1:7,9:16) = RESHAPE ((/ &  ! 2800-3000 cm-1
   1.136E+00_EB, 9.681E-01_EB, 8.842E-01_EB, 8.064E-01_EB, 6.723E-01_EB, 4.792E-01_EB, 3.539E-01_EB,&
   1.734E+00_EB, 1.285E+00_EB, 1.131E+00_EB, 1.007E+00_EB, 8.189E-01_EB, 5.838E-01_EB, 4.433E-01_EB,&
   2.211E+00_EB, 1.602E+00_EB, 1.391E+00_EB, 1.220E+00_EB, 9.645E-01_EB, 6.526E-01_EB, 4.736E-01_EB,&
   1.955E+00_EB, 1.655E+00_EB, 1.529E+00_EB, 1.417E+00_EB, 1.231E+00_EB, 9.628E-01_EB, 7.803E-01_EB,&
   3.044E+00_EB, 2.263E+00_EB, 2.004E+00_EB, 1.797E+00_EB, 1.489E+00_EB, 1.102E+00_EB, 8.681E-01_EB,&
   3.780E+00_EB, 2.608E+00_EB, 2.241E+00_EB, 1.956E+00_EB, 1.545E+00_EB, 1.063E+00_EB, 7.900E-01_EB,&
   4.023E+00_EB, 2.718E+00_EB, 2.315E+00_EB, 2.005E+00_EB, 1.563E+00_EB, 1.052E+00_EB, 7.688E-01_EB,&
   3.333E+00_EB, 2.432E+00_EB, 2.134E+00_EB, 1.898E+00_EB, 1.547E+00_EB, 1.115E+00_EB, 8.590E-01_EB/),(/7,8/))
SD3A_CH3OH(1:7,17:24) = RESHAPE ((/ &  ! 3000-3200 cm-1
   2.049E+00_EB, 1.581E+00_EB, 1.421E+00_EB, 1.291E+00_EB, 1.092E+00_EB, 8.364E-01_EB, 6.755E-01_EB,&
   1.106E+00_EB, 9.399E-01_EB, 8.769E-01_EB, 8.234E-01_EB, 7.367E-01_EB, 6.141E-01_EB, 5.290E-01_EB,&
   4.420E-01_EB, 4.578E-01_EB, 4.610E-01_EB, 4.626E-01_EB, 4.626E-01_EB, 4.560E-01_EB, 4.450E-01_EB,&
   1.748E-01_EB, 1.966E-01_EB, 2.045E-01_EB, 2.110E-01_EB, 2.211E-01_EB, 2.338E-01_EB, 2.403E-01_EB,&
   3.766E-02_EB, 6.625E-02_EB, 8.016E-02_EB, 9.340E-02_EB, 1.176E-01_EB, 1.566E-01_EB, 1.851E-01_EB,&
   1.454E-02_EB, 2.689E-02_EB, 3.296E-02_EB, 3.878E-02_EB, 4.941E-02_EB, 6.655E-02_EB, 7.894E-02_EB,&
   8.268E-03_EB, 1.479E-02_EB, 1.788E-02_EB, 2.078E-02_EB, 2.592E-02_EB, 3.373E-02_EB, 3.898E-02_EB,&
   6.234E-03_EB, 1.108E-02_EB, 1.336E-02_EB, 1.548E-02_EB, 1.923E-02_EB, 2.486E-02_EB, 2.857E-02_EB/),(/7,8/))
SD3A_CH3OH(1:7,25:28) = RESHAPE ((/ &  ! 3200-3300 cm-1
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB/),(/7,4/))

! Initialize SD3B_CH3OH array
! DATA FOR CH3OH FOR THE 2.9 MICRON BAND
! TEMP,K= 293       396      450      500     570     804       1000
SD3B_CH3OH(1:7,1:9) = RESHAPE ((/ &  ! 3600-3800 cm-1
   6.783E-02_EB, 1.141E-01_EB, 1.316E-01_EB, 1.452E-01_EB, 1.620E-01_EB, 1.690E-01_EB, 1.592E-01_EB,&
   3.984E-01_EB, 3.826E-01_EB, 3.706E-01_EB, 3.578E-01_EB, 3.322E-01_EB, 2.875E-01_EB, 2.519E-01_EB,&
   9.859E-01_EB, 6.712E-01_EB, 5.682E-01_EB, 4.877E-01_EB, 3.714E-01_EB, 2.375E-01_EB, 1.656E-01_EB,&
   1.008E+00_EB, 6.497E-01_EB, 5.417E-01_EB, 4.602E-01_EB, 3.466E-01_EB, 2.211E-01_EB, 1.555E-01_EB,&
   1.022E+00_EB, 6.996E-01_EB, 5.985E-01_EB, 5.203E-01_EB, 4.083E-01_EB, 2.782E-01_EB, 2.063E-01_EB,&
   2.633E-01_EB, 2.573E-01_EB, 2.535E-01_EB, 2.496E-01_EB, 2.421E-01_EB, 2.288E-01_EB, 2.175E-01_EB,&
   4.815E-02_EB, 6.640E-02_EB, 7.391E-02_EB, 8.052E-02_EB, 9.149E-02_EB, 1.070E-01_EB, 1.170E-01_EB,&
   1.793E-02_EB, 2.922E-02_EB, 3.428E-02_EB, 3.887E-02_EB, 4.677E-02_EB, 5.837E-02_EB, 6.600E-02_EB,&
   0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB, 0.000E+00_EB/),(/7,9/))

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
.0_EB  , .0_EB , .16_EB, .0_EB ,&
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

! METHANE NEW DATA ARRAYS
DEALLOCATE(SD_CH4_TEMP)
DEALLOCATE(OM_BND_CH4)
DEALLOCATE(SD3_CH4_NEW)
DEALLOCATE(SD7_CH4_NEW)

! PROPANE NEW DATA ARRAYS
DEALLOCATE(SD_C3H8_TEMP)
DEALLOCATE(OM_BND_C3H8)
DEALLOCATE(SD3_C3H8)
DEALLOCATE(SD7_C3H8)

! HEPTANE NEW DATA ARRAYS
DEALLOCATE(SD_C7H16_TEMP)
DEALLOCATE(OM_BND_C7H16)
DEALLOCATE(SD3_C7H16)
DEALLOCATE(SD7_C7H16)

! METHANOL NEW DATA ARRAYS
DEALLOCATE(SD_CH3OH_TEMP)
DEALLOCATE(OM_BND_CH3OH)
DEALLOCATE(SD10_CH3OH)
DEALLOCATE(SD8_CH3OH)
DEALLOCATE(SD3A_CH3OH)
DEALLOCATE(SD3B_CH3OH)
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
