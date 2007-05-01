      MODULE RAD
C
C Radiation heat transfer
C
      USE PREC
      USE CONS
      USE VARS
      USE RADCONS
      USE RADCALV
      USE MIEV
C
      IMPLICIT NONE
      PRIVATE
C
C*******************************************************************
C
C     BBFRAC    Fraction of blackbody radiation
C     DLX       Mean X-component of the control angle ray vector
C     DLY       Mean Y-component of the control angle ray vector
C     DLZ       Mean Z-component of the control angle ray vector
C     DLB       Mean Bottom-component of rayn vector (cylindrical case)
C     DLM       Mirroring indexes
C     DLN       Wall normal matrix
C     DPHI0     Opening angle of the cylindrical domain
C     E_WALL    Wall emissivity
C     ILW       Radiation intensities on solid walls, mirrors
C     KAPZT     Array of combustion gas absorption coefficients
C     KAPWT     Array of water gas absorption coefficients
C     KWR       Array of droplet radii for Mie-calculations
C     NDG       Number of droplet radii in WQABS and WQSCA arrays
C     NLMBDMIE  Number of wave lengths in Mie calculations
C     NMIEANG   Number of angle bins in forward scattering integration
C     NRA       Total number of radiation control angles
C     NRDMIE    Number of droplet radii in Mie calculations
C     NRT       Number of radiation theta angles
C     NRP       Number of radiation phi angles on each theta band
C     NSB       Number of spectral bands 
C               (1=gray, 6=wide band, 9=wide band w. CH4)
C     PATH      Mean path length for the gray gas abs. coef.
C     PHIUP, PHILOW     Phi limits of the solid angle
C     RADTMP    Radiation temperature for absorption properties (Mie)
C     RSA       Array of solid angles
C     THETAUP,THETALOW  Theta limits of the solid angle
C     UII       Integrated intensity
C     UIID      Parts of UII
C        if WIDE_BAND_MODEL = TRUE, UIID contains the band specific intensity
C        else                       UIID contains the ANGLE_INCREMENTs of intensity
C     WL_LOW, WL_HIGH   Wavelength limits of the spectral bands
C     WQABS     Absorption efficiency factor array 
C     WQSCA     Scattering efficiency factor array
C     Y_CO2, Y_H2O species mass fractions for finite-rate-reac. absorption
C     Y_CO2_MAX Maximum value of CO2 mass fraction found in fires. (finite-rate-reac)
C     Y_SOOT_F  Soot mass fraction in the flame zone (finite-rate-reac)
C
C     ANGLE_INCREMENT         How many angles are skipped on each update
C     NUMBER_RADIATION_ANGLES Input for NRA
C     NUMBER_SPECTRAL_BANDS   Input for NSB
C     TIME_STEP_INCREMENT     How often is the radiation solver called
C
C*******************************************************************
C
      REAL(EB), ALLOCATABLE, DIMENSION(:,:)       :: DLN
      REAL(EB), ALLOCATABLE, DIMENSION(:) :: RSA, DLX, DLY, DLZ, DLB
      REAL(EB)  DPHI0, R4PI, W_AXI, FOUR_SIGMA, RPI_SIGMA
C
      INTEGER, ALLOCATABLE, DIMENSION(:,:)      :: DLM
      INTEGER, ALLOCATABLE, DIMENSION(:)        :: NRP
      INTEGER  NSB, NRA, NRT, TIME_STEP_INCREMENT,
     .         ANGLE_INCREMENT, UIIDIM
C
C     Spectral and absorption coeff. stuff
C
      LOGICAL  WIDE_BAND_MODEL, CH4_BANDS
      INTEGER  NLAMBDAT, NKAPPAT, NKAPPAZ
      REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: KAPZT
      REAL(EB), ALLOCATABLE, DIMENSION(:,:)   :: KAPWT
      REAL(EB) LTSTEP, PATH,Y_CO2_MAX
C
      PUBLIC INIT_RADIATION, RADIATION_FVM, NSB, NRA,ANGLE_INCREMENT,
     .      TIME_STEP_INCREMENT,UIIDIM,NRT,RSA,NRP,PATH
C
      CONTAINS
C
C
      SUBROUTINE INIT_RADIATION
C
      REAL(EB) THETAUP, THETALOW, PHIUP, PHILOW, F_THETA
      REAL(EB) PLANCK_C2, KSI, LT, ZZ, RCRHO,
     .     BBF, AP0, AMEAN,MTOT,
     .     XLENG,YLENG,ZLENG
      INTEGER N,I,J,K,IZERO,NN,NI,II,JJ,IOS,IIM,JJM,
     .     IBND,IYY
      REAL(EB) Y_H2O, Y_SOOT,Y_CO2, Y_SOOT_F
      REAL(EB) NU_CO2,NU_SOOT, FTMP
      NAMELIST /RADI/ TIME_STEP_INCREMENT,NUMBER_RADIATION_ANGLES,
     .                ANGLE_INCREMENT,KAPPA0,
     .                WIDE_BAND_MODEL, CH4_BANDS, PATH,
     .                NMIEANG, RADTMP
C
C     Defaults
C
      NUMBER_RADIATION_ANGLES = 100
      TIME_STEP_INCREMENT = 3
      IF (TWO_D) THEN
         NUMBER_RADIATION_ANGLES = 50
         TIME_STEP_INCREMENT = 2
      ENDIF
      ANGLE_INCREMENT = -1
C
      KAPPA0          = 0.
      RADTMP          = 900.
C
      WIDE_BAND_MODEL = .FALSE.
      CH4_BANDS       = .FALSE.
      NMIEANG         = 15
      PATH            = -1.0 ! calculate path based on the geometry
C
C     Read input file looking for RADI namelist group
C
      READ_LOOP: DO
      CALL CHECKREAD('RADI',LU5,IOS) ; IF (IOS.EQ.1) EXIT READ_LOOP
      READ(LU5,RADI,END=23,ERR=24,IOSTAT=IOS)
   24 IF (IOS.GT.0) THEN
         CALL SHUTDOWN(' ERROR: Problem with RADI line')
         ENDIF
      ENDDO READ_LOOP
   23 REWIND(LU5)
C
      RADTMP = RADTMP + TMPM
C
C     Determine the number of polar angles (theta)
C
      NRA = NUMBER_RADIATION_ANGLES
      IF (CYLINDRICAL) THEN
         NRT = NINT(SQRT(REAL(NRA)))
      ELSEIF (TWO_D) THEN
         NRT = 1
      ELSE
         NRT = 2*NINT(0.5*1.17*REAL(NRA)**(1./2.26))
      ENDIF      
C
      ALLOCATE(NRP(1:NRT),STAT=IZERO)
      CALL ChkMemErr('INIT','NRP',IZERO)
C
C     Determine number of azimuthal angles (phi)
C
      N = 0
      DO I=1,NRT
         IF (CYLINDRICAL) THEN
            NRP(I) = NINT(REAL(NRA)/(REAL(NRT)))
            ELSEIF (TWO_D) THEN
            NRP(I) = 4*NINT(0.25*REAL(NRA))
            ELSE
            THETALOW = PI*REAL(I-1)/REAL(NRT)
            THETAUP  = PI*REAL(I)/REAL(NRT)
            NRP(I) = NINT(0.5*REAL(NRA)*(COS(THETALOW)-COS(THETAUP)))
            NRP(I) = MAX(4,NRP(I))
            NRP(I) = 4*NINT(0.25*REAL(NRP(I)))
            ENDIF
         N = N + NRP(I)
      ENDDO
      NRA = N
      NUMBER_RADIATION_ANGLES = NRA
C
C     Set the opening angle of the cylindrical geometry
C     equal to the azimuthal angle
C
      IF (CYLINDRICAL) DPHI0 = PI/REAL(NRP(1))
C
C     Determine the number of arrays used for intermittent intensities
C
      IF (ANGLE_INCREMENT.LT.0) 
     .   ANGLE_INCREMENT = MAX(1,MIN(5,NUMBER_RADIATION_ANGLES/15))
      UIIDIM = ANGLE_INCREMENT
C
C     Set the number of spectral bands
C
      IF (WIDE_BAND_MODEL) THEN
         IF (CH4_BANDS) THEN
            NUMBER_SPECTRAL_BANDS = 9
         ELSE
            NUMBER_SPECTRAL_BANDS = 6
         ENDIF
         TIME_STEP_INCREMENT=MAX(1,TIME_STEP_INCREMENT*ANGLE_INCREMENT)
         ANGLE_INCREMENT = 1
         UIIDIM=NUMBER_SPECTRAL_BANDS
      ELSE
         NUMBER_SPECTRAL_BANDS = 1
      ENDIF
      NSB = NUMBER_SPECTRAL_BANDS
C
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
C
C     Determine mean direction normals and sweeping orders
C
      N = 0
      DO I=1,NRT
      DO J=1,NRP(I)
         N = N + 1
         THETALOW  = PI*REAL(I-1)/REAL(NRT)
         THETAUP   = PI*REAL(I)/REAL(NRT)
         F_THETA   = 0.5*(THETAUP-THETALOW
     .      - COS(THETAUP)*SIN(THETAUP) + COS(THETALOW)*SIN(THETALOW))
         IF (CYLINDRICAL) THEN
            PHILOW = PI*REAL(J-1)/REAL(NRP(I))
            PHIUP  = PI*REAL(J)/REAL(NRP(I))
         ELSEIF (TWO_D) THEN
            PHILOW = TWOPI*REAL(J-1)/REAL(NRP(I)) + PI/2.
            PHIUP  = TWOPI*REAL(J)/REAL(NRP(I))   + PI/2.
         ELSE
            PHILOW = TWOPI*REAL(J-1)/REAL(NRP(I))
            PHIUP  = TWOPI*REAL(J)/REAL(NRP(I))
         ENDIF
         RSA(N) = (PHIUP-PHILOW)*(COS(THETALOW)-COS(THETAUP))
         IF (CYLINDRICAL) THEN
            DLX(N) = 2.*SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW)) *F_THETA
            DLY(N) =  (-SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW))
     .                 +COS(DPHI0/2.)*(COS(PHILOW)-COS(PHIUP)))*F_THETA
            DLB(N) =  (-SIN(DPHI0/2.)*(SIN(PHIUP)-SIN(PHILOW))
     .                 -COS(DPHI0/2.)*(COS(PHILOW)-COS(PHIUP)))*F_THETA
            DLZ(N)    = 0.5*(PHIUP-PHILOW)
     .           * ((SIN(THETAUP))**2-(SIN(THETALOW))**2)
         ELSEIF (TWO_D) THEN
            DLX(N) = (SIN(PHIUP)-SIN(PHILOW))*F_THETA
            DLY(N) = 0._EB
            DLZ(N) = (COS(PHILOW)-COS(PHIUP))*F_THETA
         ELSE
            DLX(N) = (SIN(PHIUP)-SIN(PHILOW))*F_THETA
            DLY(N) = (COS(PHILOW)-COS(PHIUP))*F_THETA
            DLZ(N)    = 0.5*(PHIUP-PHILOW)
     .           * ((SIN(THETAUP))**2-(SIN(THETALOW))**2)
         ENDIF
      ENDDO
      ENDDO
C
C     Set (wall normal)*(angle vector) value
C
      DO N = 1,NRA
         DLN(-1,N) = -DLX(N)
         DLN( 1,N) =  DLX(N)
         DLN(-2,N) = -DLY(N)
         DLN( 2,N) =  DLY(N)
         DLN(-3,N) = -DLZ(N)
         DLN( 3,N) =  DLZ(N)
      ENDDO
C
C     Calculate mirroring matrix
C
      N = 0
      DO I=1,NRT
      DO J=1,NRP(I)
      N = N + 1
      DO K=1,3
C
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
C
         JJM = MODULO(JJM,NRP(I))
         IF (JJM.EQ.0) JJM = NRP(I)
C
         ELSE
C
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
C
         IIM = MODULO(IIM,NRT)
         JJM = MODULO(JJM,NRP(I))
         IF (IIM.EQ.0) IIM = NRT
         IF (JJM.EQ.0) JJM = NRP(I)
         ENDIF
C
         NN = 0
         DO II = 1,IIM
         DO JJ = 1,NRP(II)
            NN = NN + 1
            IF ((II.EQ.IIM).AND.(JJ.EQ.JJM)) NI = NN
         ENDDO
         ENDDO
         DLM(N,K) = NI
      ENDDO
      ENDDO
      ENDDO
C
C-----------------------------------------------------
C
C     Spectral information
C
C-----------------------------------------------------
      INIT_WIDE_BAND: IF (WIDE_BAND_MODEL) THEN
C
C     Fraction of blackbody emission in a vavelength interval
C
      PLANCK_C2 = 14387.69             ! micron.K
      NLAMBDAT = 4000
      LTSTEP   = 25.0           ! maximum LAMBDA*T = NLANBDAT*LTSTEP
      ALLOCATE(BBFRAC(0:NLAMBDAT),STAT=IZERO)
      CALL ChkMemErr('INIT','BBFRAC',IZERO)
      BBFRAC = 0.
      LT = 0.
      DO I = 1,NLAMBDAT
         LT =  LT + LTSTEP
         KSI = PLANCK_C2/LT
         DO J = 1,50
            BBFRAC(I) = BBFRAC(I) +  EXP(-KSI*REAL(J))/REAL(J)
     .      * (   KSI**3         + 3.*KSI**2/REAL(J) + 
     .         6.*KSI/REAL(J)**2 + 6./REAL(J)**3)
         ENDDO
      ENDDO
      BBFRAC =  BBFRAC * 15./PI**4
C
C     Define band limit wave lengths in micro meters
C
      ALLOCATE(WL_LOW(1:NSB),STAT=IZERO)
      CALL ChkMemErr('INIT','WL_LOW',IZERO)
      ALLOCATE(WL_HIGH(1:NSB),STAT=IZERO)
      CALL ChkMemErr('INIT','WL_HIGH',IZERO)
      IF (CH4_BANDS) THEN
         WL_LOW(1:NSB) =(/1.0, 2.63, 2.94, 3.57, 4.17, 4.6, 
     .                    7.0, 8.62, 10.0 /)
         WL_HIGH(1:NSB)=(/     2.63, 2.94, 3.57, 4.17, 4.6, 
     .                    7.0, 8.62, 10.0, 200. /) 
      ELSE
         WL_LOW(1:NSB) =(/1.0, 2.63, 2.94, 4.17, 4.6, 10.0 /)
         WL_HIGH(1:NSB)=(/     2.63, 2.94, 4.17, 4.6, 10.0, 200.0 /)
      ENDIF
c
      ENDIF INIT_WIDE_BAND
c
c
      MIXTURE_FRACTION_IF: IF (MIXTURE_FRACTION) THEN
C-----------------------------------------------------
C
C     Tables for gas phase absorption coefficient
C
C     **************
C
C     CONTROLLING PROGRAM FOR SUBROUTINE "RADCAL", A NARROW-BAND
C     MODEL FOR CALCULATING SPECTRAL INTENSITY (W/M-2/SR/MICRON) AND
C     SPECTRAL TRANSMITTANCE VERSUS WAVELENGTH (MICRONS) IN A NONISO-
C     THERMAL, VARIABLE COMPOSITION  MIXTURE OF CO2, H2O, CO, N2, O2,
C     CH4, AND SOOT.   FOR A HOMOGENEOUS PATH, THE PROGRAM ALSO COMPUTES
C     THE PLANCK-MEAN ABSORPTION COEF., AP0, THE INCIDENT-MEAN ABSORPTION
C     COEFFICIENT, AIWALL, AND THE EFFECTIVE-MEAN ABSORPTION COEFFICIENT,
C     AMEAN, ALL IN UNITS OF INVERSE METERS.
C
C     INPUT PARAMETERS:
C          NPT=NUMBER OF HOMOGENEOUS ELEMENTS
C          DD(J)=THICKNESS OF J TH ELEMENT, M
C          RCT(J)=TEMPERATURE OF J TH ELEMENT, K.
C          P(I,J)=PARTIAL PRESSURE OF GASEOUS COMPONENTS, kPa:
C                  I   GASEOUS SPECIES
C                  1        CO2
C                  2        H2O
C                  3        CH4
C                  4        CO
C                  5        O2
C                  6        N2
C          SVF(J)=SOOT VOLUME FRACTION OF J TH ELEMENT
C          OMMIN=MINIMUM WAVE NUMBER IN SPECTRUM, CM-1.
C          OMMAX=MAXIMUM WAVE NUMBER IN SPECTRUM, CM-1.
C
C
      CALL RCALLOC
C
C     20% of mean beam length Eq 8-51 J.P. Holman 7th Edition Heat Transfer.
C     Length = 3.6 Volume/Area
C
      XLENG = MESH(1)%XF-MESH(1)%XS
      YLENG = MESH(1)%YF-MESH(1)%YS
      ZLENG = MESH(1)%ZF-MESH(1)%ZS
      IF (PATH.LT.0.0) THEN  ! default was -1.0
         IF (TWO_D) THEN ! calculate based on the geometry
            PATH = MIN( 10._EB , 0.1*3.6*XLENG*ZLENG/(XLENG+ZLENG) )
         ELSE
            PATH = MIN( 10._EB , 0.1*3.6*XLENG*YLENG*ZLENG/
     .             (XLENG*YLENG+XLENG*ZLENG+YLENG*ZLENG) )
         ENDIF
      ENDIF
      DD(1) = MAX(PATH,1.0E-4_EB)
C
      NKAPPAT = 21
      NKAPPAZ = 40
C
      ALLOCATE(KAPZT(0:NKAPPAZ,0:NKAPPAT,1:NSB),STAT=IZERO)
      CALL ChkMemErr('INIT','KAPZT',IZERO)
C
C Water vapor absorption array
C
      IF (WATER_EVAPORATION) THEN
         ALLOCATE(KAPWT(0:NKAPPAT,1:NSB),STAT=IZERO)
         CALL ChkMemErr('INIT','KAPWT',IZERO)
         ENDIF
C
      DO IBND = 1,NSB
C
      IF (NSB.GT.1) THEN
         OMMIN = REAL(NINT(1.E4/WL_HIGH(IBND)),EB)
         OMMAX = REAL(NINT(1.E4/WL_LOW(IBND)),EB)
      ELSE
         OMMIN = 50.
         OMMAX = 10000.
      ENDIF
C
      CALL INIT_RADCAL
C
      DO K = 0,NKAPPAT
      RCT(1) = 300. + K*(2400.-300.)/NKAPPAT
C
      WATERIF: IF (WATER_EVAPORATION) THEN
         SPECIE(1) = 0.
         SPECIE(2) = 1.
         SPECIE(3) = 0.
         SPECIE(4) = 0.
         SPECIE(5) = 0.
         P(1,1) = 0.
         P(2,1) = 1.
         P(3,1) = 0.
         P(4,1) = 0.
         P(5,1) = 0.
         P(6,1) = 0.
         SVF(1) = 0.
         CALL RADCAL(AMEAN,AP0)
         IF (NSB.EQ.1 .AND. PATH.GT.0.) THEN
            KAPWT(K,IBND) = AMEAN
         ELSE
            IF (NSB.EQ.1) THEN
            BBF = 1.
            ELSE
            BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),RCT(1))
            ENDIF
            KAPWT(K,IBND) = AP0/BBF
         ENDIF
      ENDIF WATERIF
C
      DO J=0,NKAPPAZ
         ZZ = IYY2ZZ(J,Z_F,NKAPPAZ)
         IYY = MAX(0,NINT(ZZ*10000.))
         IYY = MIN(10000,IYY)
         RCRHO = MW_AVG(IYY)*PINF/(R0*RCT(1))
         MTOT = Y_STATE(IYY,1)/MW_FUEL + Y_STATE(IYY,2)/MW_O2+
     .          Y_STATE(IYY,3)/MW_N2   + Y_STATE(IYY,4)/MW_H2O+
     .          Y_STATE(IYY,5)/MW_CO2  + Y_STATE(IYY,6)/MW_CO
         SPECIE(1) = Y_STATE(IYY,5)
         SPECIE(2) = Y_STATE(IYY,4)
         SPECIE(3) = Y_STATE(IYY,1)
         SPECIE(4) = Y_STATE(IYY,6)
         SPECIE(5) = Y_STATE(IYY,8)*RCRHO/RHO_SOOT
         P(1,1) = Y_STATE(IYY,5)/MW_CO2/MTOT
         P(2,1) = Y_STATE(IYY,4)/MW_H2O/MTOT
         P(3,1) = Y_STATE(IYY,1)/MW_FUEL/MTOT
         P(4,1) = Y_STATE(IYY,6)/MW_CO/MTOT
         P(5,1) = Y_STATE(IYY,2)/MW_O2/MTOT
         P(6,1) = Y_STATE(IYY,3)/MW_N2/MTOT
         SVF(1) = Y_STATE(IYY,8)*RCRHO/RHO_SOOT
         CALL RADCAL(AMEAN,AP0)
         IF ((NSB .EQ. 1).AND.(PATH.GT.0.0)) THEN
            KAPZT(J,K,IBND) = MIN(AMEAN,AP0)
         ELSE
            IF (NSB.EQ.1) THEN
            BBF = 1.
            ELSE
            BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),RCT(1))
            ENDIF
            KAPZT(J,K,IBND) = AP0/BBF
         ENDIF
      ENDDO
      ENDDO
      ENDDO
      CALL RCDEALLOC
C
      ELSE MIXTURE_FRACTION_IF  ! Finite rate reaction
C
      CO2_EXISTS: IF (ICO2.GT.0 .AND. NUN(ICO2).NE.0.) THEN
C
         CALL RCALLOC
C
         DD(1) = 1.0E-4_EB
C
         NKAPPAT = 21
         NKAPPAZ = 40
C
         ALLOCATE(KAPZT(0:NKAPPAZ,0:NKAPPAT,1:NSB),STAT=IZERO)
         CALL ChkMemErr('INIT','KAPZT',IZERO)
C
         BAND_LOOP: DO IBND = 1,NSB
C
         IF (NSB.GT.1) THEN
            OMMIN = REAL(NINT(1.E4/WL_HIGH(IBND)),EB)
            OMMAX = REAL(NINT(1.E4/WL_LOW(IBND)),EB)
         ELSE
            OMMIN = 50.
            OMMAX = 10000.
         ENDIF
C
         CALL INIT_RADCAL
C
         Y_SOOT_F = 0.05
         Y_CO2_MAX = 0.2
C
         T_LOOP: DO K = 0,NKAPPAT
C
         RCT(1) = 300. + K*(2400.-300.)/NKAPPAT
C
         Z_LOOP: DO J=0,NKAPPAZ
C
C        At cold region, KAPPA is based on CO2 concentration 
C        and SOOT_YIELD
C
C        AT over 1000 K, increase SOOT mass fraction so that peak 
C        is Y_SOOT_F (typically about 0.08 for ethylene ??) at 1800 K.
C
            ZZ = J*Y_CO2_MAX/NKAPPAZ
            Y_CO2 = ZZ
            NU_SOOT = SOOT_YIELD*MW_FUEL/MW_C
            NU_CO2  = NUN(ICO2) - NU_SOOT
            Y_SOOT  = ZZ*NU_SOOT*MW_C/(NU_CO2*MW_CO2)
            IF (RCT(1) .GT. 1000. ) THEN ! flame region
               FTMP = MIN(1.0_EB,(RCT(1)-1000.)/800.)
               Y_SOOT  = (1.-FTMP)*Y_SOOT + FTMP*(ZZ/Y_CO2_MAX)*Y_SOOT_F
               Y_CO2 = (1.-FTMP)*Y_CO2
               ENDIF
            IF (IWATER.GT.0 .AND. NUN(IWATER).NE.0.) THEN
               Y_H2O = Y_CO2*NUN(IWATER)*MW_H2O/(NU_CO2*MW_CO2)
               ELSE
               Y_H2O = Y_CO2
               ENDIF
C
            RCRHO = (29. + 10.*(Y_CO2/Y_CO2_MAX))*PINF/(R0*RCT(1))
            MTOT = 0./MW_FUEL + 0./MW_O2+
     .          (1.-ZZ-Y_H2O)/MW_N2   + Y_H2O/MW_H2O+
     .          ZZ/MW_CO2  + 0./MW_CO
            SPECIE(1) = Y_CO2
            SPECIE(2) = Y_H2O
            SPECIE(3) = 0.
            SPECIE(4) = 0.
            SPECIE(5) = Y_SOOT*RCRHO/RHO_SOOT
            P(1,1) = Y_CO2/MW_CO2/MTOT
            P(2,1) = Y_H2O/MW_H2O/MTOT
            P(3,1) = 0./MW_FUEL/MTOT
            P(4,1) = 0./MW_CO/MTOT
            P(5,1) = 0./MW_O2/MTOT
            P(6,1) = (1.-Y_CO2-Y_H2O-Y_SOOT)/MW_N2/MTOT
            SVF(1) = Y_SOOT*RCRHO/RHO_SOOT
            CALL RADCAL(AMEAN,AP0)
            IF (NSB.EQ.1 .AND. PATH.GT.0.) THEN
               KAPZT(J,K,IBND) = MIN(AMEAN,AP0)
            ELSE
               IF (NSB.EQ.1) THEN
                  BBF = 1.
                  ELSE
                  BBF = BLACKBODY_FRACTION(WL_LOW(IBND),
     .                                     WL_HIGH(IBND),RCT(1))
                  ENDIF
               KAPZT(J,K,IBND) = AP0/BBF
            ENDIF
C
      ENDDO Z_LOOP
      ENDDO T_LOOP
      ENDDO BAND_LOOP
C
      CALL RCDEALLOC
C
      ENDIF CO2_EXISTS
C
      ENDIF MIXTURE_FRACTION_IF
!
!-----------------------------------------------------
!
!     Tables for droplet absorption coefficients
!     ******************************************
!
      SPRINKLERS: IF (WATER_EVAPORATION .OR. FUEL_EVAPORATION) THEN
C
      NDG= 33
      ALLOCATE(WQABS(0:NDG,1:NSB))
      CALL ChkMemErr('INIT','WQABS',IZERO)
      ALLOCATE(WQSCA(0:NDG,1:NSB))
      CALL ChkMemErr('INIT','WQSCA',IZERO)
      ALLOCATE(KWR(0:NDG))
      CALL ChkMemErr('INIT','KWR',IZERO)
C
      WQABS = 0.
      WQSCA = 0.
      KWR = 0.
C
      IF (FUEL_EVAPORATION) THEN
C
C     Fuel absorption derived from:
C     A. Tuntomo, C.W. Tien, and S.H. Park, "Optical Constants of Liquid
C     Hydrocarbon Fuels," Comb. Sci. and Tech., 84, 133-140, 1992.
      WQABS(:,1) = (/0,10,16,28,52,98,191,386,792,1629,3272,6163,
     .        10389,15588,20807,23011,22123,22342,22200,22241,21856,
     .        22795,23633,24427,25285,26207,27006,27728,28364,28866,
     .        29260,29552,29748,30000/)
C
      KWR(0) = 0.
      DO I=1,NDG
         KWR(I)=EXP(I/2.5-4.)
      ENDDO
      WQABS=WQABS/10000.
      KWR=KWR/1000000.
C
      ENDIF
C
      IF (WATER_EVAPORATION) CALL MEAN_CROSS_SECTIONS
C
      ENDIF SPRINKLERS
C
C-----------------------------------------------------
C
C     Constants
C
      R4PI       = 1./(4.*PI)
      FOUR_SIGMA = 4.*SIGMA
      RPI_SIGMA  = RPI*SIGMA
C
C     In axially symmetric case, each angle represent two symmetric
C     angles. So, weight the intensities by two.
      W_AXI = 1.
      IF (CYLINDRICAL) W_AXI = 2.
C
      END SUBROUTINE INIT_RADIATION
C
C
C
      SUBROUTINE RADIATION_FVM(NM)
C
      USE PACKER
C
      REAL(EB) TNOW_RAD, ZZ,
     .     RAP, AX, AXU, AXD, AY, AYU, AYD, AZ, VC, RU, RD, RP,
     .     ILXU, ILYU, ILZU, KAPW, MFW,QVAL,
     .     BBF, BBFA, NCSDROP, RSA_RAT, WAXIDLN
      INTEGER N,IIG,JJG,KKG,I,J,K,IW,II,JJ,KK,IOR,IC,
     .     IWUP, IWDOWN,
     .     ISTART, IEND, ISTEP,
     .     JSTART, JEND, JSTEP,
     .     KSTART, KEND, KSTEP,
     .     NSTART, NEND, NSTEP,
     .     I_UIID, N_UPDATES,
     .     IBND, TYY, NOM
C
C Counters
C
      LOGICAL UPDATE_INTENSITY
C
      REAL(EB), POINTER, DIMENSION(:,:,:) :: KFST4, IL, UIIOLD, KAPPAW,
     .                                       KFST4W, EXTCOE, SCAEFF
      REAL(EB), POINTER, DIMENSION(:)     :: ST4_W, INRAD_W
      INTEGER, INTENT(IN) :: NM
      TYPE (OMESH_TYPE), POINTER :: M2
C
      IF (EVACUATION_ONLY(NM)) RETURN
C
      CALL UNPACK_VAR(NM)
C
      KFST4   => WORK1
      IL      => WORK2
      UIIOLD  => WORK3
      EXTCOE  => WORK4
      KAPPAW  => WORK5
      SCAEFF  => WORK6
      KFST4W  => WORK7
      ST4_W   => WALL_WORK1
      INRAD_W => WALL_WORK2
!
      TNOW_RAD=SECOND()
!
!     Ratio of solid angle, used in scattering
!
      RSA_RAT = 1./(1.-1./NRA)
!
!     Check if it time to update radiation intensity field
!
      RAD_CALL_COUNTER  = RAD_CALL_COUNTER+1
      IF ( (MOD(RAD_CALL_COUNTER,TIME_STEP_INCREMENT).EQ.1) .OR.
     .     (TIME_STEP_INCREMENT.EQ.1)) THEN
         UPDATE_INTENSITY = .TRUE.
         ELSE
         UPDATE_INTENSITY = .FALSE.
         ENDIF

      IF (WIDE_BAND_MODEL) QR = 0.
!
!     Loop over spectral bands
!
      BANDLOOP: DO IBND = 1,NUMBER_SPECTRAL_BANDS
!
      KAPPAW = 0.
      KFST4W = 0.
      SCAEFF = 0.
!
!     Calculate fraction on ambient black body radiation
!
      IF (NUMBER_SPECTRAL_BANDS.GT.1) THEN
         BBFA = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMPA)
         ELSE
         BBFA = 1.0
         ENDIF
!
!     Generate water absorption and scattering coefficients
!
      IF_DROPLETS_INCLUDED: IF (NLP.GT.0 .AND. 
     .            (WATER_EVAPORATION .OR. FUEL_EVAPORATION) ) THEN
!
      DO K=1,KBAR
      DO J=1,JBAR
      ZLOOPM: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K)) .OR. AVG_DROP_DEN(I,J,K).EQ.0.) CYCLE ZLOOPM
!
!     Calculate number density * cross section product
!
      NCSDROP = THFO*AVG_DROP_DEN(I,J,K)/
     .          (LAGRANGIAN(NPC)%DENSITY*AVG_DROP_RAD(I,J,K))
!
!     Interpolate absorption and scattering efficiency
!
      CALL INTERPOLATE1D(NDG,KWR,WQABS(:,IBND),
     .                      0.95*AVG_DROP_RAD(I,J,K),QVAL)
      KAPPAW(I,J,K) = NCSDROP*QVAL
      IF (WATER_EVAPORATION) THEN
         CALL INTERPOLATE1D(NDG,KWR,WQSCA(:,IBND),
     .                      0.95*AVG_DROP_RAD(I,J,K),QVAL)
         SCAEFF(I,J,K) = NCSDROP*QVAL
         ENDIF
!
      ENDDO ZLOOPM
      ENDDO
      ENDDO
C
      IF (NUMBER_SPECTRAL_BANDS.EQ.1) THEN
         BBF = 1.
         ELSE
         BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),RADTMP)
         ENDIF
      KFST4W = BBF*KAPPAW*FOUR_SIGMA*AVG_DROP_TMP**4
      QR_W = 0.
!
      ENDIF IF_DROPLETS_INCLUDED
!
!     Compute absorption coefficient KAPPA and source kappa*4*sigma*T**4
!
      BBF = 1.
C
      DO K=1,KBAR
      DO J=1,JBAR
      ZLOOP: DO I=1,IBAR
      IF (SOLID(ICA(I,J,K))) CYCLE ZLOOP
C
      TYY = NINT(NKAPPAT*(TMP(I,J,K)-300.)/2100.)
      TYY = MAX(0,MIN(NKAPPAT,TYY))
C
      MIXTURE_FRACTION_IF: IF (MIXTURE_FRACTION) THEN
C
         ZZ  = MAX(0._EB,YY(I,J,K,IFUEL))
         ZZ  = MIN(ZZ,1._EB)
         KAPPA(I,J,K) = ZZ2KAPPA(ZZ,Z_F,NKAPPAZ,TYY,IBND)
         IF (WATER_EVAPORATION) THEN
            MFW  = YY(I,J,K,IWATER)
            KAPW = KAPWT(TYY,IBND)
            KAPPA(I,J,K) = KAPPA(I,J,K)*(1.-MFW)+MFW*KAPW
            ENDIF
C
      ELSE MIXTURE_FRACTION_IF  ! no mixture fraction
C
      CO2_EXISTS: IF (ICO2.NE.0 .AND. NUN(ICO2).NE.0.) THEN 
         ZZ  = YY(I,J,K,ICO2)
         ZZ  = MIN(ZZ,Y_CO2_MAX)
         KAPPA(I,J,K) = ZZ2KAPPA(ZZ,Y_CO2_MAX,2*NKAPPAZ,TYY,IBND)
         ELSE CO2_EXISTS
         KAPPA(I,J,K) = KAPPA0
         ENDIF CO2_EXISTS
      ENDIF MIXTURE_FRACTION_IF
C
      IF (WIDE_BAND_MODEL)  
     .   BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMP(I,J,K))
C
      KFST4(I,J,K) = BBF*KAPPA(I,J,K)*FOUR_SIGMA*TMP(I,J,K)**4
C
      IF (BBF*CHI_R*Q(I,J,K).GT.KFST4(I,J,K)) THEN
         KFST4(I,J,K) = BBF*CHI_R*Q(I,J,K)
         KAPPA(I,J,K) = 0.
         ENDIF
C
      ENDDO ZLOOP
      ENDDO
      ENDDO
C
C     Calculate extinction coefficient
C
      EXTCOE = KAPPA + KAPPAW + SCAEFF*RSA_RAT
C
C     Update intensity field
C
      INTENSITY_UPDATE: IF (UPDATE_INTENSITY) THEN
C
      IF (WIDE_BAND_MODEL) THEN
         UIIOLD = UIID(:,:,:,IBND)
         ELSE
         UIIOLD = UII
         ENDIF
      UII = 0.
C
C     Compute boundary condition term 4*sigma*Tw**4
C
      IF (WIDE_BAND_MODEL) THEN
         DO IW = 1,NWC
         BBF = BLACKBODY_FRACTION(WL_LOW(IBND),WL_HIGH(IBND),TMP_F(IW))
         ST4_W(IW) = BBF*RPI_SIGMA*TMP_F(IW)**4
         ENDDO
         ELSE
         ST4_W(1:NWC) = RPI_SIGMA*TMP_F(1:NWC)**4
         ENDIF
C
C     Compute boundary condition term incoming radiation integral
C
      DO IW = 1,NWC
         IF (IV(IW).NE.1) CYCLE
         IOR = IJKW(4,IW)
         INRAD_W(IW) = SUM(-W_AXI*DLN(IOR,:)* WALL(IW)%ILW(:,IBND),1,
     .                            DLN(IOR,:).LT.0.)
      ENDDO
C
C     If updating intensities first time, sweep ALL angles
C
      N_UPDATES = 1
      IF (RAD_CALL_COUNTER.EQ.1) N_UPDATES = ANGLE_INCREMENT
      UIIDIMLOOP: DO I_UIID = 1,N_UPDATES
C
C     Update counters inside the radiation routine
C
      ANGLE_INC_COUNTER = MOD(ANGLE_INC_COUNTER,ANGLE_INCREMENT) + 1
C
      IF (WIDE_BAND_MODEL) THEN
         UIID(:,:,:,IBND) = 0.
      ELSE
         UIID(:,:,:,ANGLE_INC_COUNTER) = 0.
      ENDIF
C
C     Set the bounds and increment for the angleloop
C     Stepping downdard, because in cylindrical case
C     the Nth angle boundary condition comes from (N+1)th angle.
C
      NSTART = NRA - ANGLE_INC_COUNTER + 1
      NEND   = 1
      NSTEP  = -ANGLE_INCREMENT
C
C     Sweep through control angles
C
      IL(:,:,:) = BBFA*RPI_SIGMA*TMPA4
      ANGLELOOP: DO N = NSTART,NEND,NSTEP
C
C     Boundary conditions: Intensities leaving the boundaries.
C
      WLOOP1: DO IW=1,NWC
         IF (IV(IW).EQ.0)       CYCLE WLOOP1
C
         IOR = IJKW(4,IW)
         IF (DLN(IOR,N) .LT. 0.) CYCLE WLOOP1
C
         II  = IJKW(1,IW)
         JJ  = IJKW(2,IW)
         KK  = IJKW(3,IW)
C
         IF (JBAR.GT.1 .OR. ABS(IOR).NE.2) THEN
C
C        Select boundary type
C        Open, Mirror, interpolated or Solid wall
C        Set boundary value to the dummy cell
C
         SELECT CASE(IV(IW))
         CASE(2) ! open
            IL(II,JJ,KK) = BBFA*RPI_SIGMA*TMPA4
         CASE(3) ! mirror
            WALL(IW)%ILW(N,IBND) = WALL(IW)%ILW(DLM(N,ABS(IOR)),IBND)
            IL(II,JJ,KK) = WALL(IW)%ILW(N,IBND)
         CASE(4) ! interpolated
            IL(II,JJ,KK) = WALL(IW)%ILW(N,IBND)
         CASE DEFAULT ! solid wall
            WALL(IW)%ILW(N,IBND) = E_WALL(IW) * ST4_W(IW) +
     .          RPI*(1.-E_WALL(IW))* INRAD_W(IW)
         END SELECT
         ELSEIF (CYLINDRICAL) THEN
            IF (IV(IW).EQ.2) CYCLE WLOOP1
            IL(II,JJ,KK) = WALL(IW)%ILW(N,IBND)
         ENDIF
C
      ENDDO WLOOP1
C
C     Determine sweep direction in physical space
C
      ISTART = 1 ; IEND   = IBAR ; ISTEP  = 1
      JSTART = 1 ; JEND   = JBAR ; JSTEP  = 1
      KSTART = 1 ; KEND   = KBAR ; KSTEP  = 1
      IF (DLX(N) .LT. 0.) THEN
         ISTART = IBAR ; IEND   = 1 ; ISTEP  = -1
      ENDIF
      IF (DLY(N) .LT. 0.) THEN
         JSTART = JBAR ; JEND   = 1 ; JSTEP  = -1
         ENDIF
      IF (DLZ(N) .LT. 0.) THEN
         KSTART = KBAR ; KEND   = 1 ; KSTEP  = -1
         ENDIF
C
      IF (CYLINDRICAL) THEN
C
C     Sweep in axisymmetric geometry
C
      J = 1
      CKLOOP: DO K=KSTART,KEND,KSTEP
      CILOOP: DO I=ISTART,IEND,ISTEP
C
         IC = ICA(I,J,K)
         IF (SOLID(IC)) CYCLE CILOOP
C
         ILXU = IL(I-ISTEP,J,K)
         ILYU = IL(I,J-JSTEP,K)
         ILZU = IL(I,J,K-KSTEP)
C
         IF (IC.NE.0) THEN
            IW = IWA(IC,-ISTEP)
            IF (IV(IW).EQ.1) ILXU = WALL(IW)%ILW(N,IBND)
            IW = IWA(IC,-JSTEP*2)
            IF (IV(IW).EQ.1) ILYU = WALL(IW)%ILW(N,IBND)
            IW = IWA(IC,-KSTEP*3)
            IF (IV(IW).EQ.1) ILZU = WALL(IW)%ILW(N,IBND)
         ENDIF
C
         IF (DLX(N).GE.0.) THEN
            RU  = R(I-1)
            RD  = R(I)
         ELSE
            RU  = R(I)
            RD  = R(I-1)
         ENDIF
         RP  = SQRT(0.5*(RU**2+RD**2))
         VC  = DX(I)  * RP*DPHI0 * DZ(K)
         AXU =          RU       * DZ(K) * ABS(DLX(N))
         AXD =          RD       * DZ(K) * ABS(DLX(N))
         AYU = DX(I)             * DZ(K) * ABS(DLB(N))
         AYD = DX(I)             * DZ(K) * ABS(DLY(N))
         AZ  = DX(I)  * RP*DPHI0         * ABS(DLZ(N))
C
C     Zero out the terms involving symmetric overhang
C
         IF (MODULO(N,NRP(1)).EQ.1) AYD = 0.
         IF (MODULO(N,NRP(1)).EQ.0) AYU = 0.
C
         RAP = 1./(AXD+AYD+AZ+EXTCOE(I,J,K)*VC*RSA(N))
         IL(I,J,K) = MAX(0._EB, RAP * (
     .      AXU*ILXU + AYU*ILYU + AZ*ILZU + 
     .      VC*RSA(N)*R4PI*( KFST4(I,J,K)+KFST4W(I,J,K) +
     .      RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
      ENDDO CILOOP
      ENDDO CKLOOP
C
      ELSEIF (JBAR.EQ.1) THEN
C
C     Sweep in 2D cartesian geometry
C
      J = 1
      K2LOOP: DO K=KSTART,KEND,KSTEP
      I2LOOP: DO I=ISTART,IEND,ISTEP
C
         IC = ICA(I,J,K)
         IF (SOLID(IC)) CYCLE I2LOOP
C
         ILXU  = IL(I-ISTEP,J,K)
         ILZU  = IL(I,J,K-KSTEP)
C
         IF (IC.NE.0) THEN
            IW = IWA(IC,-ISTEP)
            IF (IV(IW).EQ.1) ILXU = WALL(IW)%ILW(N,IBND)
            IW = IWA(IC,-KSTEP*3)
            IF (IV(IW).EQ.1) ILZU = WALL(IW)%ILW(N,IBND)
         ENDIF
C
         VC  = DX(I) * DZ(K)
         AX  =         DZ(K) * ABS(DLX(N))
         AZ  = DX(I)         * ABS(DLZ(N))
         RAP = 1./(AX+AZ+EXTCOE(I,J,K)*VC*RSA(N))
         IL(I,J,K) = MAX(0._EB, RAP * (
     .      AX*ILXU + AZ*ILZU + 
     .      VC*RSA(N)*R4PI*(KFST4(I,J,K)+KFST4W(I,J,K) + 
     .      RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) ) 
      ENDDO I2LOOP
      ENDDO K2LOOP
C
      ELSE
C
C     Sweep in 3D cartesian geometry
C
      KLOOP: DO K=KSTART,KEND,KSTEP
      JLOOP: DO J=JSTART,JEND,JSTEP
      ILOOP: DO I=ISTART,IEND,ISTEP
C
         IC = ICA(I,J,K)
         IF (SOLID(IC)) CYCLE ILOOP
C
         ILXU  = IL(I-ISTEP,J,K)
         ILYU  = IL(I,J-JSTEP,K)
         ILZU  = IL(I,J,K-KSTEP)
C
         IF (IC.NE.0) THEN
            IW = IWA(IC,-ISTEP)
            IF (IV(IW).EQ.1) ILXU = WALL(IW)%ILW(N,IBND)
            IW = IWA(IC,-JSTEP*2)
            IF (IV(IW).EQ.1) ILYU = WALL(IW)%ILW(N,IBND)
            IW = IWA(IC,-KSTEP*3)
            IF (IV(IW).EQ.1) ILZU = WALL(IW)%ILW(N,IBND)
            ENDIF
C
         VC  = DX(I) * DY(J) * DZ(K)
         AX  =         DY(J) * DZ(K) * ABS(DLX(N))
         AY  = DX(I)         * DZ(K) * ABS(DLY(N))
         AZ  = DX(I) * DY(J)         * ABS(DLZ(N))
         RAP = 1./(AX+AY+AZ+EXTCOE(I,J,K)*VC*RSA(N))
         IL(I,J,K) = MAX(0._EB, RAP * (
     .      AX*ILXU + AY*ILYU + AZ*ILZU + 
     .      VC*RSA(N)*R4PI*( KFST4(I,J,K)+KFST4W(I,J,K) +
     .      RSA_RAT*SCAEFF(I,J,K)*UIIOLD(I,J,K) ) ) )
C
      ENDDO ILOOP
      ENDDO JLOOP
      ENDDO KLOOP
C
      ENDIF
C
C Boundary values: Incoming radiation
C
      WLOOP2: DO IW=1,NWC
C
      IF (IV(IW).EQ.0) CYCLE WLOOP2     ! not a wall
      IF (IV(IW).EQ.2) CYCLE WLOOP2     ! open boundary
      IOR = IJKW(4,IW)
      IF (TWO_D .AND. .NOT.CYLINDRICAL
     .     .AND. ABS(IOR).EQ.2) CYCLE WLOOP2  ! 2-D non cylindrical
      IF (DLN(IOR,N).GE.0.) CYCLE WLOOP2     ! outgoing
C
      WAXIDLN = - W_AXI*DLN(IOR,N)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      INRAD_W(IW) = INRAD_W(IW) - WAXIDLN * WALL(IW)%ILW(N,IBND) ! update incoming radiation,step 1
      WALL(IW)%ILW(N,IBND) = IL(IIG,JJG,KKG)
      INRAD_W(IW) = INRAD_W(IW) + WAXIDLN * WALL(IW)%ILW(N,IBND) ! update incoming radiation,step 2
C
      ENDDO WLOOP2
C
C     Copy the Y-downwind intensities to Y-upwind in cylindrical case
C
      IF (CYLINDRICAL) THEN
C
         J=1
         DO K=1,KBAR
         DO I=1,IBAR
            IWUP   = IWA(ICA(I,0,K),2)
            IWDOWN = IWA(ICA(I,2,K),-2)
            WALL(IWUP)%ILW(MAX(1,N-1),IBND) = WALL(IWDOWN)%ILW(N,IBND)
         ENDDO
         ENDDO
C
      ENDIF
C
C Calculate integrated intensity UIID
C
      IF (WIDE_BAND_MODEL) THEN
         UIID(:,:,:,IBND) =
     .   UIID(:,:,:,IBND) + W_AXI*RSA(N)*IL
      ELSE
         UIID(:,:,:,ANGLE_INC_COUNTER) =
     .   UIID(:,:,:,ANGLE_INC_COUNTER) + W_AXI*RSA(N)*IL
      ENDIF
C
C Interpolate boundary intensities onto other meshes
C
      INTERPOLATION_LOOP: DO NOM=1,NMESHES
      IF (NM.EQ.NOM) CYCLE INTERPOLATION_LOOP
      IF (NIC(NOM,NM).EQ.0) CYCLE INTERPOLATION_LOOP
      M2=>OMESH(NOM)
      OTHER_WALL_LOOP: DO IW=1,MESH(NOM)%NEWC
cc    IF (M2%IJKW(9,IW).NE.NM) CYCLE OTHER_WALL_LOOP
      IF (M2%IJKW(9,IW).NE.NM .OR. M2%IV(IW).NE.4) CYCLE OTHER_WALL_LOOP
      IOR = M2%IJKW(4,IW)
      IF (DLN(IOR,N).LE.0.) CYCLE OTHER_WALL_LOOP
      M2%WALL(IW)%ILW(N,IBND)=
     .   IL(M2%IJKW(10,IW),M2%IJKW(11,IW),M2%IJKW(12,IW))
      ENDDO OTHER_WALL_LOOP
      ENDDO INTERPOLATION_LOOP
C
      ENDDO ANGLELOOP
      ENDDO UIIDIMLOOP
C
      ENDIF INTENSITY_UPDATE
C
C     Save source term for the energy equation (QR = -DIV Q)
C
      IF (WIDE_BAND_MODEL) THEN
         QR = QR + KAPPA*UIID(:,:,:,IBND)-KFST4
         IF (NLP.GT.0 .AND. 
     .      (WATER_EVAPORATION .OR. FUEL_EVAPORATION) ) 
     .      QR_W = QR_W + KAPPAW*UIID(:,:,:,IBND) - KFST4W
      ENDIF
C
      ENDDO BANDLOOP
C
C     Sum up the parts of the intensity
C
      IF (UPDATE_INTENSITY) UII = SUM(UIID, DIM = 4)
C
C     Save source term for the energy equation (QR = -DIV Q)
C     Done only in one-band (gray gas) case.
C
      IF (.NOT. WIDE_BAND_MODEL) THEN
         QR  = KAPPA*UII - KFST4
         IF (NLP.GT.0 .AND. 
     .          (WATER_EVAPORATION .OR. FUEL_EVAPORATION) ) 
     .      QR_W = QR_W + KAPPAW*UII - KFST4W
         ENDIF
C
C     Radiation flux boundary conditions for temperatures
C
      QRAD = 0.
      CALL RADIATION_BC
C
      TUSED(9,NM)=TUSED(9,NM)+SECOND()-TNOW_RAD
C
      CONTAINS
C
C
      SUBROUTINE RADIATION_BC
C
C     Get radiative flux on each boundary cell
C
      INTEGER :: IBC
C
      WCLOOP: DO IW=1,NWC
C
      IF (IV(IW).NE.1) CYCLE WCLOOP
C
      IOR = IJKW(4,IW)
      IBC = IJKW(5,IW)
      QRAD(IW) = -W_AXI*SUM(DLN(IOR,:)*SUM(WALL(IW)%ILW(:,:),DIM=2))
      QRAD(IW) = QRAD(IW) + E_WALL(IW)*QEXT(IBC)
C
      ENDDO WCLOOP
C
      END SUBROUTINE RADIATION_BC
C
      END SUBROUTINE RADIATION_FVM
C
      REAL (EB) FUNCTION ZZ2KAPPA(ZZ,ZZF,NIYY,TTD,BAND)
C
C     Calculate index IYY as a function of ZZ.
C
      REAL(EB) ZZ,ZZF,IYY
      INTEGER TTD,BAND,IY1,IY2,ZINT,NIYY
C
      IF (ZZ.LE.ZZF) THEN
        IYY = ZZ*REAL(NIYY/2,EB)/ZZF
      ELSE
        IYY = NIYY/2 + REAL(NIYY-NIYY/2,EB)*(ZZ-ZZF)/(1.-ZZF) 
      ENDIF
      ZINT   = NINT(IYY)
C
      IF(IYY.EQ.ZINT) THEN
      ZZ2KAPPA=KAPZT(ZINT,TTD,BAND)
      ELSE
      ZINT = INT(IYY)
      IY2  = ZINT + 1
      IY1  = ZINT
      ZZ2KAPPA=KAPZT(IY1,TTD,BAND)+(IYY-REAL(IY1,EB))*
     .      (KAPZT(IY2,TTD,BAND)-KAPZT(IY1,TTD,BAND))
      ENDIF
      END FUNCTION ZZ2KAPPA
C
C
      REAL(EB) FUNCTION IYY2ZZ(IYY,ZZF,NIYY)
C
C     Function gives the mixture fraction as a function of index
C     IYY 
C
      REAL(EB) ZZF
      INTEGER IYY, NIYY
C
      IF (IYY.LE.NIYY/2) THEN
        IYY2ZZ = REAL(IYY,EB)*ZZF/REAL(NIYY/2,EB)
      ELSE
        IYY2ZZ = ZZF+(1.-ZZF)*REAL(IYY-NIYY/2,EB)/REAL(NIYY-NIYY/2,EB)
      ENDIF
C
      END FUNCTION IYY2ZZ
c
      REAL(EB) FUNCTION BLACKBODY_FRACTION(L1,L2,TEMP)
C
C     Calculates the fraction of black body radiation
C     between wavelengths L1 and L2 (micron) in Temperature TEMP
C
      REAL(EB) L1, L2, TEMP
      REAL(EB) LT1, LT2, BBFLOW, BBFHIGH
      INTEGER IYY

      LT1    =   L1 * TEMP/LTSTEP
      LT2    =   L2 * TEMP/LTSTEP
      IYY = MIN(NLAMBDAT,MAX(0,NINT(LT1)))
      BBFLOW = BBFRAC(IYY)
      IYY = MIN(NLAMBDAT,MAX(0,NINT(LT2)))
      BBFHIGH = BBFRAC(IYY)
      BLACKBODY_FRACTION = BBFHIGH - BBFLOW

      END FUNCTION BLACKBODY_FRACTION
C
      END MODULE RAD
