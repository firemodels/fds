!> \brief Global constants, parameters, variables 
!>
!> \details Each MPI process stores a copy of the GLOBAL_CONSTANTS. It cannot be
!> assumed that these values are same for each MPI process.

MODULE GLOBAL_CONSTANTS

USE PRECISION_PARAMETERS
IMPLICIT NONE

INTEGER, PARAMETER :: DNS_MODE=1                 !< SIM_MODE is Direct Numerical Simulation
INTEGER, PARAMETER :: LES_MODE=2                 !< SIM_MODE is Large Eddy Simulation
INTEGER, PARAMETER :: VLES_MODE=3                !< SIM_MODE is Very Large Eddy Simulation
INTEGER, PARAMETER :: SVLES_MODE=4               !< SIM_MODE is Simple Very Large Eddy Simulation
INTEGER, PARAMETER :: GAS_SPECIES=2              !< Flag for SPECIES\%MODE indicating a gaseous species
INTEGER, PARAMETER :: AEROSOL_SPECIES=3          !< Flag for SPECIES\%MODE indicating an aerosol species
INTEGER, PARAMETER :: EXPLICIT_EULER=1           !< COMBUSTION_ODE_SOLVER is explicit first-order Euler
INTEGER, PARAMETER :: RK2=2                      !< COMBUSTION_ODE_SOLVER is second-order Runge-Kutta
INTEGER, PARAMETER :: RK3=3                      !< COMBUSTION_ODE_SOLVER is third-order Runge-Kutta
INTEGER, PARAMETER :: RK2_RICHARDSON=4           !< COMBUSTION_ODE_SOLVER is second-order Runge-Kutta with Richardson extrapolation
INTEGER, PARAMETER :: EXTINCTION_1=1             !< Flag for EXTINCT_MOD (EXTINCTION MODEL 1)
INTEGER, PARAMETER :: EXTINCTION_2=2             !< Flag for EXTINCT_MOD (EXTINCTION MODEL 2)
INTEGER, PARAMETER :: NO_TURB_MODEL=0            !< No turbulence model (DNS)
INTEGER, PARAMETER :: CONSMAG=1                  !< Constant Smagorinsky turbulence model
INTEGER, PARAMETER :: DYNSMAG=2                  !< Dynamic Smagorinsky turbulence model
INTEGER, PARAMETER :: DEARDORFF=3                !< Deardorff turbulence model
INTEGER, PARAMETER :: VREMAN=4                   !< Vreman turbulence model
INTEGER, PARAMETER :: RNG=5                      !< ReNormalization Group turbulence model
INTEGER, PARAMETER :: WALE=6                     !< Wall-Adapting Local Eddy viscosity turbulence model
INTEGER, PARAMETER :: MU_TURB_INTERP=7           !< Flag for NEAR_WALL_TURB_MODEL that avoids jump in viscosity, \f$ \mu \f$
INTEGER, PARAMETER :: CONVECTIVE_FLUX_BC=-1,NET_FLUX_BC=0,SPECIFIED_TEMPERATURE=1,&
                      NO_CONVECTION=2,THERMALLY_THICK=3,INFLOW_OUTFLOW=4,&
                      INTERPOLATED_BC=6,THERMALLY_THICK_HT3D=7
INTEGER, PARAMETER :: H_CUSTOM=-1,H_DEFAULT=0,H_LOGLAW=1,H_ABL=2,H_RAYLEIGH=3,H_YUAN=4,&
                      H_FREE_HORIZONTAL_CYLINDER=5,H_BLOWING_SPHERE=6,H_BLOWING_PLATE=7     ! Heat transfer coefficient model
INTEGER, PARAMETER :: WALL_MODEL_BC=2,FREE_SLIP_BC=3,NO_SLIP_BC=4,BOUNDARY_FUEL_MODEL_BC=5,&
                      INTERPOLATED_VELOCITY_BC=6                                            ! Velocity BC
INTEGER, PARAMETER :: EXPOSED=0,VOID=1,INSULATED=2                                          ! Surface backing
INTEGER, PARAMETER :: SURF_CARTESIAN=0,SURF_CYLINDRICAL=1,SURF_SPHERICAL=2,SURF_BLOWING_PLATE=3 ! Surface geometry
INTEGER, PARAMETER :: NO_MASS_FLUX=1,SPECIFIED_MASS_FRACTION=2,SPECIFIED_MASS_FLUX=3,&
                      INFLOW_OUTFLOW_MASS_FLUX=4                                            ! Mass transfer BC
INTEGER, PARAMETER :: NULL_BOUNDARY=0,SOLID_BOUNDARY=1,OPEN_BOUNDARY=2,MIRROR_BOUNDARY=3, & ! Boundary Types
                      INTERPOLATED_BOUNDARY=6,PERIODIC_BOUNDARY=7,SLICE_BOUNDARY=8
INTEGER, PARAMETER :: FISHPAK_BC_PERIODIC=0, &                                              ! Fishpak (FFT solver) BC
                      FISHPAK_BC_DIRICHLET_DIRICHLET=1, &
                      FISHPAK_BC_DIRICHLET_NEUMANN=2, &
                      FISHPAK_BC_NEUMANN_NEUMANN=3, &
                      FISHPAK_BC_NEUMANN_DIRICHLET=4
INTEGER, PARAMETER :: DIRICHLET=1,NEUMANN=2,INTERNAL=3                                      ! Pressure Boundary Conditions
INTEGER, PARAMETER :: PYROLYSIS_NONE=0,PYROLYSIS_SOLID=1,PYROLYSIS_LIQUID=2,&
                      PYROLYSIS_PREDICTED=3,PYROLYSIS_SPECIFIED=4,PYROLYSIS_VEGETATION=5     ! Pyrolysis model
INTEGER, PARAMETER :: ATMOSPHERIC_PROFILE=1, PARABOLIC_PROFILE=2,&
                      BOUNDARY_LAYER_PROFILE=3, RAMP_PROFILE=4                              ! Surface velocity profile
INTEGER, PARAMETER :: CELL_CENTER=1, CELL_FACE=2, CELL_EDGE=3                               ! Cell position for output quantities
INTEGER, PARAMETER :: NO_STOP=0, INSTABILITY_STOP=1, USER_STOP=2, SETUP_STOP=3, &
                      SETUP_ONLY_STOP=4, CTRL_STOP=5, TGA_ANALYSIS_STOP=6, &
                      LEVELSET_STOP=7, REALIZABILITY_STOP=8, EVACUATION_STOP=9, &
                      VERSION_STOP=10                                                       ! Program Stop Status
INTEGER, PARAMETER :: SPHERE_DRAG=1,CYLINDER_DRAG=2,USER_DRAG=3,SCREEN_DRAG=4,POROUS_DRAG=5 ! Options for drag computation
INTEGER, PARAMETER :: OLD=1,NEW=2,GUESS=3,PREVIOUS=4                                        ! Network solver indices
INTEGER, PARAMETER :: NODE1=1,NODE2=2,LINEAR_INTERPOLATION=-1                               ! Duct mass transport interpolation
INTEGER, PARAMETER :: HVAC_BOUNDARY=42
INTEGER, PARAMETER :: OBST_SPHERE_TYPE=1,OBST_CYLINDER_TYPE=2,OBST_CONE_TYPE=3,OBST_BOX_TYPE=4
INTEGER :: N_SIMPLE_CHEMISTRY_REACTIONS=1

INTEGER :: FUEL_INDEX=0,O2_INDEX=0,N2_INDEX=0,H2O_INDEX=0,CO2_INDEX=0,CO_INDEX=0,H2_INDEX=0,SOOT_INDEX=0
INTEGER :: H2O_SMIX_INDEX = -1
INTEGER :: HCN_INDEX=0,NO_INDEX=0,NO2_INDEX=0

INTEGER :: ZETA_INDEX=0                    !< Index of the unmixed fuel fraction, ZETA

INTEGER :: STOP_STATUS=NO_STOP             !< Indicator of whether and why to stop the job
INTEGER :: INPUT_FILE_LINE_NUMBER=0        !< Indicator of what line in the input file is being read

REAL(EB) :: FUEL_C_TO_CO_FRACTION=0.6667_EB !< Fraction of carbon atoms in the fuel that are converted to CO
REAL(EB) :: FUEL_H_TO_H2_FRACTION=0._EB     !< Fraction of hydrogen atoms in the fuel that are converted to H2

! Miscellaneous logical constants

LOGICAL :: RADIATION=.TRUE.                 !< Perform radiation transport
LOGICAL :: RADIATION_COMPLETED=.TRUE.       !< A complete radiation update has completed
LOGICAL :: EXCHANGE_RADIATION=.FALSE.       !< Do an MPI radiation exchange at this time step
LOGICAL :: CYLINDRICAL=.FALSE.              !< Cylindrical domain option
LOGICAL :: NOISE=.TRUE.                     !< Initialize velocity field with a small amount of divergence-free motion
LOGICAL :: PREDICTOR                        !< The first half of the second-order accurate time-step
LOGICAL :: CORRECTOR                        !< The second half of the second-order accurate time-step
LOGICAL :: INITIALIZATION_PHASE=.TRUE.      !< The set-up phase before the time-stepping loop
LOGICAL :: APPEND=.FALSE.                   !< For a RESTARTed calculation, APPEND the exising output files
LOGICAL :: PARTICLE_FILE=.FALSE.            !< Indicates the existence of Lagrangian particles
LOGICAL :: RESTART=.FALSE.                  !< Indicates if a former calculation is to be RESTARTed
LOGICAL :: SUPPRESSION=.TRUE.               !< Indicates if gas-phase combustion extinction is modeled
LOGICAL :: ACCUMULATE_WATER=.FALSE.         !< Indicates that integrated liquid outputs are specified
LOGICAL :: WRITE_XYZ=.FALSE.                !< Indicates that a Plot3D geometry file is specified by user
LOGICAL :: CHECK_POISSON=.FALSE.            !< Check the accuracy of the Poisson solver
LOGICAL :: TWO_D=.FALSE.                    !< Perform a 2-D simulation
LOGICAL :: SETUP_ONLY=.FALSE.               !< Indicates that the calculation should be stopped before time-stepping
LOGICAL :: CHECK_MESH_ALIGNMENT=.FALSE.     !< Indicates that the user wants to check the mesh alignment and then stop
LOGICAL :: SMOKE3D=.TRUE.                   !< Indicates that the 3D smoke and fire output is desired
LOGICAL :: STATUS_FILES=.FALSE.             !< Produce an output file CHID.notready which is deleted if the simulation completes
LOGICAL :: LOCK_TIME_STEP=.FALSE.           !< Do not allow time step to change for diagnostic purposes
LOGICAL :: RESTRICT_TIME_STEP=.TRUE.        !< Do not let the time step increase above its intial value
LOGICAL :: FLUSH_FILE_BUFFERS=.TRUE.        !< Periodically flush the output buffers during simulation for better Smokeviewing
LOGICAL :: CLIP_RESTART_FILES=.TRUE.        !< Append RESTARTed output files at the time the former calculation terminated
LOGICAL :: COLUMN_DUMP_LIMIT=.FALSE.        !< Limit the number of columns in output files
LOGICAL :: MASS_FILE=.FALSE.                !< Output a comma-delimited file of gas species masses
LOGICAL :: STRATIFICATION=.TRUE.            !< Assume that the atmosphere decreases in pressure with height
LOGICAL :: SOLID_PHASE_ONLY=.FALSE.         !< Only perform a solid phase heat transfer and pyrolysis simulation
LOGICAL :: AEROSOL_AL2O3=.FALSE.            !< Assume that the SOOT is Al_2 O_3
LOGICAL :: SHARED_FILE_SYSTEM=.TRUE.        !< Assume that FDS is being run on computers with a shared file system
LOGICAL :: FREEZE_VELOCITY=.FALSE.          !< Hold velocity fixed, do not perform a velocity update
LOGICAL :: BNDF_DEFAULT=.TRUE.              !< Output boundary output files
LOGICAL :: SPATIAL_GRAVITY_VARIATION=.FALSE.!< Assume gravity varies as a function of the \f$ x \f$ coordinate
LOGICAL :: PROJECTION=.FALSE.               !< Apply the projection method for the divergence
LOGICAL :: CHECK_VN=.TRUE.                  !< Check the Von Neumann number
LOGICAL :: SOLID_PARTICLES=.FALSE.          !< Indicates the existence of solid particles
LOGICAL :: HVAC=.FALSE.                     !< Perform an HVAC calculation
LOGICAL :: BAROCLINIC=.TRUE.                !< Include the baroclinic terms in the momentum equation
LOGICAL :: GRAVITATIONAL_DEPOSITION=.TRUE.  !< Allow aerosol gravitational deposition
LOGICAL :: GRAVITATIONAL_SETTLING=.TRUE.    !< Allow aerosol gravitational settling
LOGICAL :: THERMOPHORETIC_DEPOSITION=.TRUE. !< Allow aerosol thermophoretic deposition
LOGICAL :: THERMOPHORETIC_SETTLING=.TRUE.   !< Allow aerosol thermophoretic settling
LOGICAL :: TURBULENT_DEPOSITION=.TRUE.      !< Allow aerosol turbulent deposition
LOGICAL :: DEPOSITION=.TRUE.                !< Allow aerosol deposition
LOGICAL :: AEROSOL_SCRUBBING=.FALSE.        !< Allow aerosol scrubbing
LOGICAL :: VELOCITY_ERROR_FILE=.FALSE.      !< Generate a diagnostic output file listing velocity and pressure errors
LOGICAL :: CFL_FILE=.FALSE.                 !< Generate a diagnostic output file listing quantities related to CFL and VN
LOGICAL :: CONSTANT_SPECIFIC_HEAT_RATIO=.FALSE. !< Assume that the ratio of specific heats is constant, \f$ \gamma=1.4 \f$
LOGICAL :: MEAN_FORCING(3)=.FALSE.          !< Apply mean forcing to wind in each coordinate direction
LOGICAL :: CHECK_HT=.FALSE.                 !< Apply heat transfer stability condition
LOGICAL :: PATCH_VELOCITY=.FALSE.           !< Assume user-defined velocity patches
LOGICAL :: OVERWRITE=.TRUE.                 !< Overwrite old output files
LOGICAL :: INIT_HRRPUV=.FALSE.              !< Assume an initial spatial distribution of HRR per unit volume
LOGICAL :: INIT_HRRPUV_RAMP=.FALSE.         !< Assume that initial spatial distribution of HRR per unit volume is ramped up or down
LOGICAL :: TENSOR_DIFFUSIVITY=.FALSE.
LOGICAL :: SYNTHETIC_EDDY_METHOD=.FALSE.
LOGICAL :: UVW_RESTART=.FALSE.              !< Initialize velocity field with values from a file
LOGICAL :: PARTICLE_CFL=.FALSE.             !< Include particle velocity as a constraint on time step
LOGICAL :: IBM_FEM_COUPLING=.FALSE.
LOGICAL :: ENTHALPY_TRANSPORT=.TRUE.
LOGICAL :: CONSTANT_H_SOLID=.TRUE.
LOGICAL :: POTENTIAL_TEMPERATURE_CORRECTION=.FALSE.
LOGICAL :: RTE_SOURCE_CORRECTION=.TRUE.     !< Apply a correction to the radiation source term to achieve desired rad fraction
LOGICAL :: LAPLACE_PRESSURE_CORRECTION=.FALSE.
LOGICAL :: CHECK_REALIZABILITY=.FALSE.
LOGICAL :: MIN_DEVICES_EXIST=.FALSE.
LOGICAL :: MAX_DEVICES_EXIST=.FALSE.
LOGICAL :: SUPPRESS_DIAGNOSTICS=.FALSE.     !< Do not print detailed mesh-specific output in the .out file
LOGICAL :: WRITE_GEOM_FIRST=.TRUE.
LOGICAL :: SIMPLE_CHEMISTRY=.FALSE.         !< Use simple chemistry combustion model
LOGICAL :: FIRST_PASS                       !< The point in the time step before the CFL constraint is applied
LOGICAL :: VERBOSE=.FALSE.                  !< Add extra output in the .err file
LOGICAL :: SOLID_HT3D=.FALSE.
LOGICAL :: SOLID_MT3D=.FALSE.
LOGICAL :: SOLID_PYRO3D=.FALSE.
LOGICAL :: HVAC_MASS_TRANSPORT=.FALSE.
LOGICAL :: EXTERNAL_BOUNDARY_CORRECTION=.FALSE.
LOGICAL :: DUCT_HT=.FALSE.
LOGICAL :: DUCT_HT_INSERTED=.FALSE.
LOGICAL :: STORE_Q_DOT_PPP_S=.FALSE.
LOGICAL :: STORE_OLD_VELOCITY=.FALSE.
LOGICAL :: QFAN_BETA_TEST=.FALSE.
LOGICAL :: USE_ATMOSPHERIC_INTERPOLATION=.FALSE.
LOGICAL :: POSITIVE_ERROR_TEST=.FALSE.
LOGICAL :: OBST_SHAPE_AREA_ADJUST=.FALSE.
LOGICAL :: TRI_MODEL=.FALSE.
LOGICAL :: FLAME_INDEX_MODEL=.FALSE.
LOGICAL :: STORE_SPECIES_FLUX=.FALSE.
LOGICAL :: CHAR_OXIDATION=.FALSE.
LOGICAL :: STORE_DIVERGENCE_CORRECTION=.FALSE.
LOGICAL :: PERIODIC_DOMAIN_X=.FALSE.                !< The domain is periodic \f$ x \f$
LOGICAL :: PERIODIC_DOMAIN_Y=.FALSE.                !< The domain is periodic \f$ y \f$
LOGICAL :: PERIODIC_DOMAIN_Z=.FALSE.                !< The domain is periodic \f$ z \f$

INTEGER :: BNDF_TIME_INTEGRALS=0

INTEGER, ALLOCATABLE, DIMENSION(:) :: CHANGE_TIME_STEP_INDEX,SETUP_PRESSURE_ZONES_INDEX

! Miscellaneous character strings

CHARACTER(255) :: TITLE,RENDER_FILE,UVW_FILE='null'
CHARACTER(CHID_LENGTH) :: CHID,RESTART_CHID
CHARACTER(255) :: TERRAIN_IMAGE

! Dates, version numbers, revision numbers

REAL(FB) :: VERSION_NUMBER=6.0
CHARACTER(255) :: REVISION='unknown'
CHARACTER(255) :: REVISION_DATE='unknown'
CHARACTER(255) :: COMPILE_DATE='unknown'

! Global EVACuation parameters

LOGICAL, ALLOCATABLE, DIMENSION(:) :: EVACUATION_ONLY, EVACUATION_SKIP
REAL(EB) :: EVAC_DT_FLOWFIELD,EVAC_DT_STEADY_STATE,T_EVAC,T_EVAC_SAVE
INTEGER :: EVAC_PRESSURE_ITERATIONS,EVAC_TIME_ITERATIONS,EVAC_N_QUANTITIES,I_EVAC
INTEGER :: EVAC_AVATAR_NCOLOR
LOGICAL :: EVACUATION_MC_MODE=.FALSE.,EVACUATION_DRILL=.FALSE.,NO_EVACUATION=.FALSE.
CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: EVAC_CLASS_NAME, EVAC_CLASS_PROP
INTEGER, ALLOCATABLE, DIMENSION(:) :: EVAC_QUANTITIES_INDEX
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: EVAC_CLASS_RGB,EVAC_AVATAR_RGB
REAL(EB), ALLOCATABLE, DIMENSION(:) :: EVACUATION_Z_OFFSET

! Miscellaneous real constants

REAL(EB) :: CPOPR,RSC,RPR,TMPA,TMPA4,RHOA,P_INF,RADIATIVE_FRACTION, &
            CP_GAMMA,GAMMA,GM1OG,U0,V0,W0,H0,GVEC(3),FVEC(3)=0._EB,OVEC(3)=0._EB, &
            C_SMAGORINSKY=0.2_EB,C_DEARDORFF=0.1_EB,C_VREMAN=0.07_EB,C_RNG=0.2_EB,C_RNG_CUTOFF=10._EB,C_WALE=0.60_EB, &
            LAPSE_RATE,TEX_ORI(3),KAPPA0, &
            ASSUMED_GAS_TEMPERATURE,PR_ONTH,MU_AIR_0,CP_AIR_0,PR_AIR,K_AIR_0, &
            CHARACTERISTIC_VELOCITY,CFL_MAX,CFL_MIN,VN_MAX,VN_MIN,PR,SC,H_V_W,GROUND_LEVEL=0._EB, &
            LIMITING_DT_RATIO=1.E-4_EB,NOISE_VELOCITY=0.005_EB, &
            SCALAR_ENERGY_TOLERANCE=1.E-12_EB,DT_MEAN_FORCING=1._EB, &
            TAU_DEFAULT=1._EB,TAU_CHEM=1.E-10_EB,TAU_FLAME=1.E10_EB,SMOKE_ALBEDO=0.3_EB, &
            Y_WERNER_WENGLE=11.81_EB,PARTICLE_CFL_MAX=1.0_EB,PARTICLE_CFL_MIN=0.8_EB,DELTA_IBM,GRAV=9.80665_EB,&
            H_V_H2O(0:5000),CHI_R_MIN=0._EB,CHI_R_MAX=1._EB,EVAP_FILM_FAC=1._EB/3._EB,&
            LATITUDE=10000._EB,GEOSTROPHIC_WIND(2)=0._EB

REAL(EB), PARAMETER :: TMPM=273.15_EB, P_STP=101325._EB,R0=8314.472_EB,R1=1.986257E-03_EB,&
                       SIGMA=5.670373E-8_EB,K_BOLTZMANN=1.3806488E-23_EB,&
                       EARTH_OMEGA=7.272205216643040e-05_EB ! [radians/s] = 2*pi/(24*3600)

! Parameters associated with parallel mode

INTEGER :: MYID=0,N_MPI_PROCESSES=1,EVAC_PROCESS=-1,LOWER_MESH_INDEX=1000000000,UPPER_MESH_INDEX=-1000000000
LOGICAL :: PROFILING=.FALSE.
INTEGER, ALLOCATABLE, DIMENSION(:) :: PROCESS,FILE_COUNTER

! Time parameters

REAL(EB) :: DT_INITIAL,T_BEGIN,T_END,T_END_GEOM,TWFIN,TIME_SHRINK_FACTOR,RELAXATION_FACTOR=1._EB,MPI_TIMEOUT=300._EB,&
            DT_END_MINIMUM=2._EB*EPSILON(1._EB),DT_END_FILL=1.E-6_EB
EQUIVALENCE(T_END,TWFIN)

! Combustion parameters

REAL(EB) :: Y_O2_INFTY=0.232378_EB,Y_CO2_INFTY=0.000595_EB,Y_H2O_INFTY=0._EB,&
            MW_AIR=28.84852_EB,MW_N2,MW_O2,MW_CO2,MW_H2O,MW_CO,MW_H2,VISIBILITY_FACTOR, &
            EC_LL,ZZ_MIN_GLOBAL=1.E-10_EB,&
            FIXED_MIX_TIME=-1._EB,INITIAL_UNMIXED_FRACTION=1._EB,RICHARDSON_ERROR_TOLERANCE=1.E-6_EB,&
            H_F_REFERENCE_TEMPERATURE=25._EB,FREE_BURN_TEMPERATURE=600._EB, &
            AUTO_IGNITION_TEMPERATURE=0._EB,AIT_EXCLUSION_ZONE(6,MAX_AIT_EXCLUSION_ZONES)=-1.E6_EB
REAL(FB) :: HRRPUV_MAX_SMV=1200._FB, TEMP_MAX_SMV=2000._FB
INTEGER :: N_SPECIES=0,N_REACTIONS,I_PRODUCTS=-1,I_WATER=-1,I_CO2=-1,N_TRACKED_SPECIES=0,N_SURFACE_DENSITY_SPECIES=0,&
           COMBUSTION_ODE_SOLVER=-1,EXTINCT_MOD=-1,MAX_CHEMISTRY_SUBSTEPS=20,MAX_PRIORITY=1,&
           N_PASSIVE_SCALARS=0,N_TOTAL_SCALARS=0,N_FIXED_CHEMISTRY_SUBSTEPS=-1
LOGICAL :: SIMPLE_CHEM=.TRUE.,OUTPUT_CHEM_IT=.FALSE.,REAC_SOURCE_CHECK=.FALSE.
REAL(EB) :: RSUM0
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: Z2Y,CP_Z,CPBAR_Z,K_RSQMW_Z,MU_RSQMW_Z,D_Z,CP_AVG_Z,G_F_Z,H_SENS_Z
REAL(EB), ALLOCATABLE, DIMENSION(:) :: MWR_Z,RSQ_MW_Z
CHARACTER(LABEL_LENGTH) :: EXTINCTION_MODEL='null'

! Radiation parameters

INTEGER :: NUMBER_SPECTRAL_BANDS=0,NUMBER_RADIATION_ANGLES=0,ANGLE_INCREMENT=0,RADIATION_ITERATIONS=1, &
           INITIAL_RADIATION_ITERATIONS,NUMBER_FSK_POINTS=1
REAL(EB) :: RTE_SOURCE_CORRECTION_FACTOR=1._EB,RAD_Q_SUM=0._EB,KFST4_SUM=0._EB,QR_CLIP,C_MAX=100._EB,C_MIN=1._EB

! Ramping parameters

CHARACTER(LABEL_LENGTH), POINTER, DIMENSION(:) :: RAMP_ID,RAMP_TYPE
INTEGER :: I_RAMP_AGT,I_RAMP_U0_T,I_RAMP_V0_T,I_RAMP_W0_T,I_RAMP_U0_Z,I_RAMP_V0_Z,I_RAMP_W0_Z,I_RAMP_GX,I_RAMP_GY,I_RAMP_GZ,&
           I_RAMP_FVX_T,I_RAMP_FVY_T,I_RAMP_FVZ_T,N_RAMP=0,I_RAMP_TMP0_Z=0,I_RAMP_P0_Z=0,I_RAMP_SPEED=0,I_RAMP_DIRECTION=0
INTEGER, PARAMETER :: TIME_HEAT=-1,TIME_VELO=-2,TIME_TEMP=-3,TIME_EFLUX=-4,TIME_PART=-5,TANH_RAMP=-2,TSQR_RAMP=-1,&
                      VELO_PROF_X=-6,VELO_PROF_Y=-7,VELO_PROF_Z=-8

! TABLe parameters

CHARACTER(LABEL_LENGTH) :: TABLE_ID(1000)
INTEGER :: N_TABLE=0,TABLE_TYPE(1000)
INTEGER, PARAMETER :: SPRAY_PATTERN=1,PART_RADIATIVE_PROPERTY=2,POINTWISE_INSERTION=3,TABLE_2D_TYPE=4

! Variables related to meshes

INTEGER :: NMESHES=1,IBAR_MAX=0,JBAR_MAX=0,KBAR_MAX=0,MESH_LIST_EMB(100)
REAL(EB) :: XS_MIN=1.E6_EB,XF_MAX=-1.E6_EB,YS_MIN=1.E6_EB,YF_MAX=-1.E6_EB,ZS_MIN=1.E6_EB,ZF_MAX=-1.E6_EB
CHARACTER(LABEL_LENGTH), DIMENSION(:), ALLOCATABLE :: MESH_NAME

! Variables related to pressure solver

LOGICAL :: ITERATE_PRESSURE=.FALSE.,ITERATE_BAROCLINIC_TERM,SUSPEND_PRESSURE_ITERATIONS=.TRUE.
REAL(EB) :: VELOCITY_TOLERANCE=0._EB,PRESSURE_TOLERANCE=0._EB,ITERATION_SUSPEND_FACTOR=0.95_EB
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VELOCITY_ERROR_MAX,PRESSURE_ERROR_MAX
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: VELOCITY_ERROR_MAX_LOC,PRESSURE_ERROR_MAX_LOC
INTEGER :: PRESSURE_ITERATIONS=0,MAX_PRESSURE_ITERATIONS=10,TOTAL_PRESSURE_ITERATIONS=0
CHARACTER(LABEL_LENGTH):: PRES_METHOD = 'FFT'

!! Parameters for the Laplace solver (under construction)
!
!INTEGER, PARAMETER :: SCI_KM1=1,SCI_JM1=2,SCI_IM1=3,SCI_IJK=4,SCI_IP1=5,SCI_JP1=6,SCI_KP1=7 ! (S)parse (C)olumn (I)ndex
!INTEGER, PARAMETER :: NDIAG=7

! Miscellaneous integer constants

INTEGER :: ICYC,ICYC_RESTART=0,NFRAMES,PERIODIC_TEST=0,SIM_MODE=3,TURB_MODEL=0,SLIP_CONDITION=2,FISHPAK_BC(3)=-1,&
           STOP_AT_ITER=0,HT3D_TEST=0,WALL_INCREMENT=2,WALL_INCREMENT_HT3D=1,NEAR_WALL_TURB_MODEL=1

! Clocks for output file dumps

REAL(EB), ALLOCATABLE, DIMENSION(:) :: PART_CLOCK,SLCF_CLOCK,SL3D_CLOCK,&
                                       PL3D_CLOCK,BNDF_CLOCK,ISOF_CLOCK,PROF_CLOCK,RADF_CLOCK
REAL(EB) :: MINT_CLOCK,DEVC_CLOCK,HRR_CLOCK,EVAC_CLOCK=1.E6_EB,CTRL_CLOCK,FLUSH_CLOCK,CPU_CLOCK,RESTART_CLOCK,&
            BNDC_CLOCK,GEOC_CLOCK,GEOM_CLOCK,UVW_CLOCK,TURB_INIT_CLOCK=-1.E10_EB,&
            UVW_CLOCK_CBC(1:4)=(/0._EB,0.28_EB,0.67_EB,1.E10_EB/)
REAL(EB) :: UVW_TIMER(10),MMS_TIMER=1.E10_EB
REAL(EB) :: DT_SLCF,DT_BNDF,DT_DEVC,DT_PL3D,DT_PART,DT_RESTART,DT_ISOF,DT_HRR,DT_MASS,DT_PROF,DT_CTRL,&
            DT_FLUSH,DT_SL3D,DT_GEOM,DT_CPU,DT_RADF,DT_MOM
REAL(EB) :: T_RADF_BEGIN,T_RADF_END
LOGICAL  :: UPDATE_DEVICES_AGAIN=.FALSE.

! Miscellaneous mesh dimensions

REAL(EB) :: CHARACTERISTIC_CELL_SIZE=1.E6_EB,MESH_SEPARATION_DISTANCE,NEIGHBOR_SEPARATION_DISTANCE

! Logical units and output file names

INTEGER                              :: LU_ERR=0,LU_END=2,LU_GIT=3,LU_SMV=4,LU_INPUT=5,LU_OUTPUT=6,LU_STOP=7,LU_CPU=8,&
                                        LU_CATF=9
INTEGER                              :: LU_MASS,LU_HRR,LU_STEPS,LU_NOTREADY,LU_VELOCITY_ERROR,LU_CFL,LU_LINE=-1,LU_CUTCELL
INTEGER                              :: LU_EVACCSV,LU_EVACEFF,LU_EVACFED,LU_EVACOUT,LU_HISTOGRAM,LU_EVAC_CB
INTEGER                              :: LU_BNDC=-1,LU_GEOC=-1,LU_TGA,LU_INFO
INTEGER, ALLOCATABLE, DIMENSION(:)   :: LU_PART,LU_PROF,LU_XYZ,LU_TERRAIN,LU_PL3D,LU_DEVC,LU_STATE,LU_CTRL,LU_CORE,LU_RESTART
INTEGER, ALLOCATABLE, DIMENSION(:)   :: LU_VEG_OUT,LU_GEOM
INTEGER                              :: LU_GEOM_TRAN
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LU_SLCF,LU_SLCF_GEOM,LU_BNDF,LU_BNDF_GEOM,LU_BNDG,LU_ISOF,LU_ISOF2, &
                                        LU_SMOKE3D,LU_BNDF_SLCF,LU_RADF
INTEGER                              :: DEVC_COLUMN_LIMIT=254,CTRL_COLUMN_LIMIT=254

CHARACTER(250)                             :: FN_INPUT='null'
CHARACTER(80)                              :: FN_STOP='null',FN_CPU,FN_CFL,FN_OUTPUT='null'
CHARACTER(80)                              :: FN_MASS,FN_HRR,FN_STEPS,FN_SMV,FN_END,FN_ERR,FN_NOTREADY,FN_VELOCITY_ERROR,FN_GIT
CHARACTER(80)                              :: FN_EVACCSV,FN_EVACEFF,FN_EVACFED,FN_EVACOUT,FN_LINE,FN_HISTOGRAM,FN_CUTCELL,FN_TGA
CHARACTER(80), ALLOCATABLE, DIMENSION(:)   :: FN_PART,FN_PROF,FN_XYZ,FN_TERRAIN,FN_PL3D,FN_DEVC,FN_STATE,FN_CTRL,FN_CORE,FN_RESTART
CHARACTER(80), ALLOCATABLE, DIMENSION(:)   :: FN_VEG_OUT,FN_GEOM
CHARACTER(80), ALLOCATABLE, DIMENSION(:,:) :: FN_SLCF,FN_SLCF_GEOM,FN_BNDF,FN_BNDF_GEOM,FN_BNDG, &
                                              FN_ISOF,FN_ISOF2,FN_SMOKE3D,FN_BNDF_SLCF,FN_RADF,FN_GEOM_TRNF

CHARACTER(9) :: FMT_R
LOGICAL :: OUT_FILE_OPENED=.FALSE.

! Boundary condition arrays

CHARACTER(LABEL_LENGTH) :: MATL_NAME(1:1000)
INTEGER :: N_SURF,N_SURF_RESERVED,N_MATL,MIRROR_SURF_INDEX,OPEN_SURF_INDEX,INTERPOLATED_SURF_INDEX,DEFAULT_SURF_INDEX=0, &
           INERT_SURF_INDEX=0,PERIODIC_SURF_INDEX,PERIODIC_WIND_SURF_INDEX,HVAC_SURF_INDEX=-1,EVACUATION_SURF_INDEX=-1,&
           MASSLESS_TRACER_SURF_INDEX, MASSLESS_TARGET_SURF_INDEX,DROPLET_SURF_INDEX,VEGETATION_SURF_INDEX,NWP_MAX,&
           SLICE_SURF_INDEX
REAL(EB), ALLOCATABLE, DIMENSION(:) :: AAS,BBS,CCS,DDS,DDT,DX_S,RDX_S,RDXN_S,DX_WGT_S, &
                                       K_S,RHOCBAR,C_S,RHO_S,Q_S,TWO_DX_KAPPA_S,X_S_NEW,R_S,MF_FRAC,REGRID_FACTOR,R_S_NEW
INTEGER,  ALLOCATABLE, DIMENSION(:) :: LAYER_INDEX,CELL_COUNT

! Divergence Arrays

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: DSUM,USUM,PSUM

! Level Set vegetation fire spread

LOGICAL  :: VEG_LEVEL_SET=.FALSE.,VEG_LEVEL_SET_COUPLED=.TRUE.
LOGICAL  :: LSET_ELLIPSE,LSET_TAN2

! Parameters for Terrain and Wind simulation needs

LOGICAL :: TERRAIN_CASE=.FALSE.
INTEGER :: N_VENT_TOTAL=0,SPONGE_CELLS,N_MEAN_FORCING_BINS
REAL(EB), ALLOCATABLE, DIMENSION(:) :: MEAN_FORCING_SUM_U_VOL,MEAN_FORCING_SUM_V_VOL,MEAN_FORCING_SUM_W_VOL, &
                                       MEAN_FORCING_SUM_VOL_X,MEAN_FORCING_SUM_VOL_Y,MEAN_FORCING_SUM_VOL_Z, &
                                       U_MEAN_FORCING,V_MEAN_FORCING,W_MEAN_FORCING

! Sprinkler Variables

REAL(EB) :: C_DIMARZO=6.E6_EB
INTEGER :: N_ACTUATED_SPRINKLERS=0,N_OPEN_NOZZLES=0,EVAP_MODEL=0
INTEGER, PARAMETER :: NDC=1000,NDC2=100
LOGICAL :: POROUS_FLOOR=.TRUE.,ALLOW_UNDERSIDE_PARTICLES=.FALSE.,ALLOW_SURFACE_PARTICLES=.TRUE.

! Particles and PARTICLEs

INTEGER :: MAXIMUM_PARTICLES,N_LAGRANGIAN_CLASSES,N_EVAC,N_LP_ARRAY_INDICES=0
REAL(EB) :: CNF_CUTOFF=0.005_EB
LOGICAL :: EB_PART_FILE=.FALSE.,PL3D_PARTICLE_FLUX=.FALSE.,SLCF_PARTICLE_FLUX=.FALSE.,DEVC_PARTICLE_FLUX=.FALSE.
LOGICAL :: OMESH_PARTICLES=.FALSE.

INTEGER :: MOMENTUM_INTERPOLATION_METHOD=0

! Soot oxidation

REAL(EB), ALLOCATABLE, DIMENSION(:) :: NU_SOOT_OX

! Agglomeration model

LOGICAL :: AGGLOMERATION = .TRUE.,SOOT_OXIDATION=.FALSE.
INTEGER :: N_PARTICLE_BINS(MAX_SPECIES)=0,AGGLOMERATION_SPEC_INDEX(MAX_SPECIES)=-1,AGGLOMERATION_SMIX_INDEX(MAX_SPECIES)=-1,&
           N_AGGLOMERATION_SPECIES=0
REAL(EB) :: MIN_PARTICLE_DIAMETER(MAX_SPECIES),MAX_PARTICLE_DIAMETER(MAX_SPECIES),&
            NUCLEATION_SITES=1.E7_EB !1E7 is 10 nucleation sites per cm^3

! Number of initial value, pressure zone, and multiplier derived types

INTEGER :: N_INIT,N_ZONE,N_MULT,N_MOVE
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: CONNECTED_ZONES
LOGICAL, ALLOCATABLE, DIMENSION(:) :: SEALED_MESH
REAL(EB) :: PRESSURE_RELAX_TIME=1._EB

! Clipping values

REAL(EB) :: TMPMIN,TMPMAX,RHOMIN,RHOMAX

! Flux limiter

INTEGER, PARAMETER :: CENTRAL_LIMITER=0,GODUNOV_LIMITER=1,SUPERBEE_LIMITER=2,MINMOD_LIMITER=3,CHARM_LIMITER=4,MP5_LIMITER=5
INTEGER :: I_FLUX_LIMITER=SUPERBEE_LIMITER,CFL_VELOCITY_NORM=0

! Numerical quadrature (used in TEST_FILTER)

INTEGER, PARAMETER :: TRAPAZOID_QUADRATURE=0, SIMPSON_QUADRATURE=1, MIDPOINT_QUADRATURE=2
INTEGER :: TEST_FILTER_QUADRATURE=TRAPAZOID_QUADRATURE

INTEGER, PARAMETER :: N_TIMERS=14                   !< Number of subroutine timers
REAL(EB), ALLOCATABLE, DIMENSION(:) :: T_USED       !< Array of subroutine timings
REAL(EB) :: WALL_CLOCK_START                        !< MPI_WTIME i.e. wall clock time when FDS starts
REAL(EB) :: WALL_CLOCK_START_ITERATIONS=0._EB       !< MPI_WTIME i.e. wall clock time when main iteration loop starts
REAL(EB) :: CPU_TIME_START                          !< CPU_TIME when FDS starts

INTEGER :: OPENMP_AVAILABLE_THREADS = 1         !< OpenMP parameter
INTEGER :: OPENMP_USED_THREADS      = 1         !< OpenMP parameter
LOGICAL :: OPENMP_USER_SET_THREADS  = .FALSE.   !< OpenMP parameter
LOGICAL :: USE_OPENMP               = .FALSE.   !< OpenMP parameter

INTEGER :: N_CSVF=0  !< Number of external velocity (.csv) files

INTEGER :: N_FACE=0,N_GEOM=0
REAL(EB):: DT_BNDC=1.E10_EB

LOGICAL :: CC_IBM=.FALSE.
LOGICAL :: CHECK_MASS_CONSERVE =.FALSE.
LOGICAL :: GLMAT_SOLVER =.FALSE.
LOGICAL :: GLMAT_VERBOSE=.FALSE.
LOGICAL :: PRES_ON_WHOLE_DOMAIN=.TRUE.
LOGICAL :: PRES_ON_CARTESIAN=.TRUE.
LOGICAL :: DO_IMPLICIT_CCREGION=.FALSE.
LOGICAL :: COMPUTE_CUTCELLS_ONLY=.FALSE.
LOGICAL :: CC_ZEROIBM_VELO=.FALSE.
LOGICAL :: CC_SLIPIBM_VELO=.FALSE.

! Threshold factor for volume of cut-cells respect to volume of Cartesian cells:
! Currently used in the thermo div definition of cut-cells.

REAL(EB) :: CCVOL_LINK=0.15_EB
LOGICAL  :: GET_CUTCELLS_VERBOSE=.FALSE.

INTEGER, PARAMETER :: LOW_IND   = 1
INTEGER, PARAMETER :: HIGH_IND  = 2

INTEGER, PARAMETER :: MAX_DIM   = 3 !< Maximum number of spatial dimensions for a problem.
INTEGER, PARAMETER :: IAXIS = 1
INTEGER, PARAMETER :: JAXIS = 2
INTEGER, PARAMETER :: KAXIS = 3

INTEGER, PARAMETER :: NOD1 = 1
INTEGER, PARAMETER :: NOD2 = 2
INTEGER, PARAMETER :: NOD3 = 3
INTEGER, PARAMETER :: NOD4 = 4
INTEGER, PARAMETER :: EDG1 = 1
INTEGER, PARAMETER :: EDG2 = 2
INTEGER, PARAMETER :: EDG3 = 3

INTEGER :: MAXIMUM_GEOMETRY_ZVALS=100000, MAXIMUM_GEOMETRY_VOLUS=350000, &
           MAXIMUM_GEOMETRY_FACES=100000, MAXIMUM_GEOMETRY_VERTS=100000, &
           MAXIMUM_GEOMETRY_IDS  =1000,   MAXIMUM_GEOMETRY_SURFIDS=100

! Allocation increment parameters:

INTEGER, PARAMETER :: IBM_ALLOC_DVERT = 10
INTEGER, PARAMETER :: IBM_ALLOC_DELEM = 10

INTEGER, PARAMETER :: IBM_MAX_WSTRIANG_SEG =  2  !< Up to two ws-triangles related to a segment
INTEGER, PARAMETER :: IBM_MAX_WSTRIANG_TRI =  1  !< Up to 1 ws-triangle per BODINT_PLANE triangle

INTEGER, ALLOCATABLE, DIMENSION(:)   :: CELL_COUNT_CC
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: N_EDGES_DIM_CC

! HVAC Parameters

INTEGER :: N_DUCTNODES = 0, N_DUCTS = 0, N_FANS = 0, N_FILTERS = 0, N_AIRCOILS = 0,N_NETWORKS=0, N_DUCTRUNS=0
INTEGER , ALLOCATABLE, DIMENSION(:) :: DUCT_NE,DUCTNODE_NE,DUCT_DR,DUCTNODE_DR
REAL(EB) :: HVAC_PRES_RELAX=0.5_EB
LOGICAL :: HVAC_SOLVE=.FALSE.,HVAC_LOCAL_PRESSURE=.TRUE.

REAL(EB), POINTER, DIMENSION(:,:) :: ORIENTATION_VECTOR !< Global array of orientation vectors
INTEGER :: N_ORIENTATION_VECTOR                         !< Number of ORIENTATION_VECTORs

INTEGER :: TGA_SURF_INDEX=-100             !< Surface properties to use for special TGA calculation
INTEGER :: TGA_WALL_INDEX=-100             !< Wall index to use for special TGA calculation
INTEGER :: TGA_PARTICLE_INDEX=-100         !< Particle index to use for special TGA calculation
REAL(EB) :: TGA_HEATING_RATE=5._EB         !< Heat rate (K/min) to use for special TGA calculation
REAL(EB) :: TGA_FINAL_TEMPERATURE=800._EB  !< Final Temperature (C) to use for special TGA calculation

LOGICAL :: IBLANK_SMV=.TRUE.  !< Parameter passed to smokeview (in .smv file) to control generation of blockages

END MODULE GLOBAL_CONSTANTS
