!> \brief Global constants, parameters, variables
!>
!> \details Each MPI process stores a copy of the parameters within this module. It cannot be
!> assumed that these values are the same for each MPI process.

MODULE GLOBAL_CONSTANTS


USE PRECISION_PARAMETERS
USE MPI_F08
USE ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
IMPLICIT NONE (TYPE,EXTERNAL)

! constants used by smokeview for drawing surfaces

INTEGER, PARAMETER :: SMV_REGULAR=0              !< Parameter for smokeview drawing: surface drawn as a solid
INTEGER, PARAMETER :: SMV_TEXTURE=1              !< Parameter for smokeview drawing: texture image drawn over the surface
INTEGER, PARAMETER :: SMV_OUTLINE=2              !< Parameter for smokeview drawing: surface drawn as an outline
INTEGER, PARAMETER :: SMV_HIDDEN=-2              !< Parameter for smokeview drawing: surface not drawn

INTEGER, PARAMETER :: DNS_MODE=1                 !< Flag for SIM_MODE: Direct Numerical Simulation
INTEGER, PARAMETER :: LES_MODE=2                 !< Flag for SIM_MODE: Large Eddy Simulation
INTEGER, PARAMETER :: VLES_MODE=3                !< Flag for SIM_MODE: Very Large Eddy Simulation
INTEGER, PARAMETER :: SVLES_MODE=4               !< Flag for SIM_MODE: Simple Very Large Eddy Simulation

INTEGER, PARAMETER :: GAS_SPECIES=2              !< Flag for SPECIES\%MODE indicating a gaseous species
INTEGER, PARAMETER :: AEROSOL_SPECIES=3          !< Flag for SPECIES\%MODE indicating an aerosol species

INTEGER, PARAMETER :: EXPLICIT_EULER=1           !< Flag for COMBUSTION_ODE_SOLVER: explicit first-order Euler
INTEGER, PARAMETER :: RK2_RICHARDSON=2           !< Flag for COMBUSTION_ODE_SOLVER: second-order Runge-Kutta, Richardson extrap.
INTEGER, PARAMETER :: CVODE_SOLVER=3             !< Flag for COMBUSTION_ODE_SOLVER: SUNDIALS CVODE

INTEGER, PARAMETER :: EXTINCTION_1=1             !< Flag for EXTINCT_MOD (EXTINCTION MODEL 1)
INTEGER, PARAMETER :: EXTINCTION_2=2             !< Flag for EXTINCT_MOD (EXTINCTION MODEL 2)

INTEGER, PARAMETER :: NO_TURB_MODEL=0            !< Flag for TURB_MODEL: No turbulence model (DNS)
INTEGER, PARAMETER :: CONSMAG=1                  !< Flag for TURB_MODEL: Constant Smagorinsky turbulence model
INTEGER, PARAMETER :: DYNSMAG=2                  !< Flag for TURB_MODEL: Dynamic Smagorinsky turbulence model
INTEGER, PARAMETER :: DEARDORFF=3                !< Flag for TURB_MODEL: Deardorff turbulence model
INTEGER, PARAMETER :: VREMAN=4                   !< Flag for TURB_MODEL: Vreman turbulence model
INTEGER, PARAMETER :: WALE=5                     !< Flag for TURB_MODEL: Wall-Adapting Local Eddy viscosity turbulence model
INTEGER, PARAMETER :: CONSTANT_EDDY_VISCOSITY=6  !< Flag for NEAR_WALL_TURB_MODEL: specified near wall eddy viscosity

INTEGER, PARAMETER :: MEAN_LES_FILTER=1          !< Flag for LES_FILTER_WIDTH_TYPE: geometric mean filter width
INTEGER, PARAMETER :: MAX_LES_FILTER=2           !< Flag for LES_FILTER_WIDTH_TYPE: take max of DX,DY,DZ
INTEGER, PARAMETER :: FIXED_LES_FILTER=3         !< Flag for LES_FILTER_WIDTH_TYPE: constant value
INTEGER            :: LES_FILTER_WIDTH_TYPE=MEAN_LES_FILTER !< Default filter width type is mean
REAL(EB)           :: FIXED_LES_FILTER_WIDTH=-1._EB !< Constant LES filter width (m)

INTEGER, PARAMETER :: CONVECTIVE_FLUX_BC=-1      !< Flag for SF\%THERMAL_BC_INDEX: Specified convective flux
INTEGER, PARAMETER :: NET_FLUX_BC=0              !< Flag for SF\%THERMAL_BC_INDEX: Specified net heat flux
INTEGER, PARAMETER :: SPECIFIED_TEMPERATURE=1    !< Flag for SF\%THERMAL_BC_INDEX: Specified surface temperature
INTEGER, PARAMETER :: NO_CONVECTION=2            !< Flag for SF\%THERMAL_BC_INDEX: No heat transfer at MIRROR boundary
INTEGER, PARAMETER :: THERMALLY_THICK=3          !< Flag for SF\%THERMAL_BC_INDEX: Thermally thick 1-D solid
INTEGER, PARAMETER :: INFLOW_OUTFLOW=4           !< Flag for SF\%THERMAL_BC_INDEX: OPEN boundary
INTEGER, PARAMETER :: INTERPOLATED_BC=6          !< Flag for SF\%THERMAL_BC_INDEX: Interface between two meshes

INTEGER, PARAMETER :: DEFAULT_HTC_MODEL=0          !< Flag for SF\%HEAT_TRANSFER_MODEL
INTEGER, PARAMETER :: LOGLAW_HTC_MODEL=1           !< Flag for SF\%HEAT_TRANSFER_MODEL
INTEGER, PARAMETER :: RAYLEIGH_HTC_MODEL=3         !< Flag for SF\%HEAT_TRANSFER_MODEL
INTEGER, PARAMETER :: IMPINGING_JET_HTC_MODEL=4    !< Flag for SF\%HEAT_TRANSFER_MODEL
INTEGER, PARAMETER :: FM_HTC_MODEL=5               !< Flag for SF\%HEAT_TRANSFER_MODEL
INTEGER, PARAMETER :: UGENT_HTC_MODEL=6            !< Flag for SF\%HEAT_TRANSFER_MODEL

INTEGER, PARAMETER :: WALL_MODEL_BC=2              !< Flag for SF\%VELOCITY_BC_INDEX
INTEGER, PARAMETER :: FREE_SLIP_BC=3               !< Flag for SF\%VELOCITY_BC_INDEX
INTEGER, PARAMETER :: NO_SLIP_BC=4                 !< Flag for SF\%VELOCITY_BC_INDEX
INTEGER, PARAMETER :: BOUNDARY_FUEL_MODEL_BC=5     !< Flag for SF\%VELOCITY_BC_INDEX
INTEGER, PARAMETER :: INTERPOLATED_VELOCITY_BC=6   !< Flag for SF\%VELOCITY_BC_INDEX

INTEGER, PARAMETER :: EXPOSED=0                    !< Flag for SF\%BACKING: Exposed to conditions on the other side
INTEGER, PARAMETER :: VOID=1                       !< Flag for SF\%BACKING: Exposed to ambient void
INTEGER, PARAMETER :: INSULATED=2                  !< Flag for SF\%BACKING: Insulated, no heat transfer out the back

INTEGER, PARAMETER :: SURF_CARTESIAN=0             !< Flag for SF\%GEOMETRY: Flat surface
INTEGER, PARAMETER :: SURF_CYLINDRICAL=1           !< Flag for SF\%GEOMETRY: Outer surface of cylinder
INTEGER, PARAMETER :: SURF_INNER_CYLINDRICAL=-1    !< Flag for SF\%GEOMETRY: Inner surface of cylinder
INTEGER, PARAMETER :: SURF_SPHERICAL=2             !< Flag for SF\%GEOMETRY: Outer surface of sphere

INTEGER, PARAMETER :: NO_MASS_FLUX=1               !< Flag for SF\%SPECIES_BC_INDEX
INTEGER, PARAMETER :: SPECIFIED_MASS_FRACTION=2    !< Flag for SF\%SPECIES_BC_INDEX
INTEGER, PARAMETER :: SPECIFIED_MASS_FLUX=3        !< Flag for SF\%SPECIES_BC_INDEX
INTEGER, PARAMETER :: INFLOW_OUTFLOW_MASS_FLUX=4   !< Flag for SF\%SPECIES_BC_INDEX

INTEGER, PARAMETER :: NULL_BOUNDARY=0              !< Flag for SF\%BOUNDARY_TYPE, VT\%BOUNDARY_TYPE, WC\%BOUNDARY_TYPE
INTEGER, PARAMETER :: SOLID_BOUNDARY=1             !< Flag for SF\%BOUNDARY_TYPE, VT\%BOUNDARY_TYPE, WC\%BOUNDARY_TYPE
INTEGER, PARAMETER :: OPEN_BOUNDARY=2              !< Flag for SF\%BOUNDARY_TYPE, VT\%BOUNDARY_TYPE, WC\%BOUNDARY_TYPE
INTEGER, PARAMETER :: MIRROR_BOUNDARY=3            !< Flag for SF\%BOUNDARY_TYPE, VT\%BOUNDARY_TYPE, WC\%BOUNDARY_TYPE
INTEGER, PARAMETER :: INTERPOLATED_BOUNDARY=6      !< Flag for SF\%BOUNDARY_TYPE, VT\%BOUNDARY_TYPE, WC\%BOUNDARY_TYPE
INTEGER, PARAMETER :: PERIODIC_BOUNDARY=7          !< Flag for SF\%BOUNDARY_TYPE, VT\%BOUNDARY_TYPE, WC\%BOUNDARY_TYPE
INTEGER, PARAMETER :: HVAC_BOUNDARY=42             !< Flag for SF\%THERMAL_BC_INDEX, SF\%SPECIES_BC_INDEX, VT\%BOUNDARY_TYPE

INTEGER, PARAMETER :: FISHPAK_BC_PERIODIC=0            !< Flag for FISHPAK_BC(I) I=1,3: Period BC for Poisson solver
INTEGER, PARAMETER :: FISHPAK_BC_DIRICHLET_DIRICHLET=1 !< Flag for FISHPAK_BC(I) I=1,3: Dirichlet at both sides
INTEGER, PARAMETER :: FISHPAK_BC_DIRICHLET_NEUMANN=2   !< Flag for FISHPAK_BC(I) I=1,3: Dirichlet at lower, Neumann at upper
INTEGER, PARAMETER :: FISHPAK_BC_NEUMANN_NEUMANN=3     !< Flag for FISHPAK_BC(I) I=1,3: Neumann at both sides
INTEGER, PARAMETER :: FISHPAK_BC_NEUMANN_DIRICHLET=4   !< Flag for FISHPAK_BC(I) I=1,3: Neumann at lower, Dirichlet at upper

INTEGER, PARAMETER :: DIRICHLET=1                      !< Flag for WC\%PRESSURE_BC_INDEX
INTEGER, PARAMETER :: NEUMANN=2                        !< Flag for WC\%PRESSURE_BC_INDEX
INTEGER, PARAMETER :: INTERNAL=3                       !< Flag for WC\%PRESSURE_BC_INDEX

INTEGER, PARAMETER :: PYROLYSIS_NONE=0                 !< Flag for SF\%PYROLYSIS_MODEL, ML\%PYROLYSIS_MODEL
INTEGER, PARAMETER :: PYROLYSIS_SOLID=1                !< Flag for ML\%PYROLYSIS_MODEL
INTEGER, PARAMETER :: PYROLYSIS_LIQUID=2               !< Flag for ML\%PYROLYSIS_MODEL
INTEGER, PARAMETER :: PYROLYSIS_PREDICTED=3            !< Flag for SF\%PYROLYSIS_MODEL
INTEGER, PARAMETER :: PYROLYSIS_SPECIFIED=4            !< Flag for SF\%PYROLYSIS_MODEL
INTEGER, PARAMETER :: PYROLYSIS_SURFACE_OXIDATION=5    !< Flag for ML\%PYROLYSIS_MODEL

INTEGER, PARAMETER :: ATMOSPHERIC_PROFILE=1            !< Flag for SF\%PROFILE
INTEGER, PARAMETER :: PARABOLIC_PROFILE=2              !< Flag for SF\%PROFILE
INTEGER, PARAMETER :: BOUNDARY_LAYER_PROFILE=3         !< Flag for SF\%PROFILE
INTEGER, PARAMETER :: RAMP_PROFILE=4                   !< Flag for SF\%PROFILE

INTEGER, PARAMETER :: CELL_CENTER=1                    !< Flag for OUTPUT_QUANTITY()\%CELL_POSITION
INTEGER, PARAMETER :: CELL_FACE=2                      !< Flag for OUTPUT_QUANTITY()\%CELL_POSITION
INTEGER, PARAMETER :: CELL_EDGE=3                      !< Flag for OUTPUT_QUANTITY()\%CELL_POSITION

INTEGER, PARAMETER :: NO_STOP=0                        !< Flag for STATUS_STOP
INTEGER, PARAMETER :: INSTABILITY_STOP=1               !< Flag for STATUS_STOP
INTEGER, PARAMETER :: USER_STOP=2                      !< Flag for STATUS_STOP
INTEGER, PARAMETER :: SETUP_STOP=3                     !< Flag for STATUS_STOP
INTEGER, PARAMETER :: SETUP_ONLY_STOP=4                !< Flag for STATUS_STOP
INTEGER, PARAMETER :: CTRL_STOP=5                      !< Flag for STATUS_STOP
INTEGER, PARAMETER :: TGA_ANALYSIS_STOP=6              !< Flag for STATUS_STOP
INTEGER, PARAMETER :: LEVELSET_STOP=7                  !< Flag for STATUS_STOP
INTEGER, PARAMETER :: REALIZABILITY_STOP=8             !< Flag for STATUS_STOP
INTEGER, PARAMETER :: VERSION_STOP=10                  !< Flag for STATUS_STOP
INTEGER, PARAMETER :: ODE_STOP=11                      !< Flag for STATUS_STOP
INTEGER, PARAMETER :: HEARTBEAT_STOP=12                !< Flag for STATUS_STOP
INTEGER, PARAMETER :: CVODE_SUBSTEP_STOP=13            !< Flag for STATUS_STOP

INTEGER, PARAMETER :: SPHERE_DRAG=1                    !< Flag for LPC\%DRAG_LAW (LPC means LAGRANGIAN_PARTICLE_CLASS)
INTEGER, PARAMETER :: CYLINDER_DRAG=2                  !< Flag for LPC\%DRAG_LAW
INTEGER, PARAMETER :: USER_DRAG=3                      !< Flag for LPC\%DRAG_LAW: User-specified constant drag coefficient
INTEGER, PARAMETER :: SCREEN_DRAG=4                    !< Flag for LPC\%DRAG_LAW: Special drag model for screens
INTEGER, PARAMETER :: POROUS_DRAG=5                    !< Flag for LPC\%DRAG_LAW: Special drag model for porous media
INTEGER, PARAMETER :: DISK_DRAG=6                      !< Flag for LPC\%DRAG_LAW: Disk drag for cartesian particles

INTEGER, PARAMETER :: OLD=1                            !< Argument for DUCT()\%VEL()
INTEGER, PARAMETER :: NEW=2                            !< Argument for DUCT()\%VEL()
INTEGER, PARAMETER :: GUESS=3                          !< Argument for DUCT()\%VEL()
INTEGER, PARAMETER :: PREVIOUS=4                       !< Argument for DUCT()\%VEL()

INTEGER, PARAMETER :: OBST_SPHERE_TYPE=1               !< Flag for OB\%SHAPE_TYPE
INTEGER, PARAMETER :: OBST_CYLINDER_TYPE=2             !< Flag for OB\%SHAPE_TYPE
INTEGER, PARAMETER :: OBST_CONE_TYPE=3                 !< Flag for OB\%SHAPE_TYPE
INTEGER, PARAMETER :: OBST_BOX_TYPE=4                  !< Flag for OB\%SHAPE_TYPE

INTEGER, PARAMETER :: ARRHENIUS_TYPE            = 1    !< Flag for RN\%REACTYPE
INTEGER, PARAMETER :: THREE_BODY_ARRHENIUS_TYPE = 2    !< Flag for RN\%REACTYPE
INTEGER, PARAMETER :: FALLOFF_TROE_TYPE         = 3    !< Flag for RN\%REACTYPE
INTEGER, PARAMETER :: FALLOFF_LINDEMANN_TYPE    = 4    !< Flag for RN\%REACTYPE

INTEGER :: FUEL_INDEX=0                    !< Index for FUEL in SIMPLE_CHEMISTRY model
INTEGER :: O2_INDEX=0                      !< Index for O2 in SIMPLE_CHEMISTRY model
INTEGER :: N2_INDEX=0                      !< Index for N2 in SIMPLE_CHEMISTRY model
INTEGER :: H2O_INDEX=0                     !< Index for H2O in SIMPLE_CHEMISTRY model
INTEGER :: CO2_INDEX=0                     !< Index for CO2 in SIMPLE_CHEMISTRY model
INTEGER :: CO_INDEX=0                      !< Index for CO in SIMPLE_CHEMISTRY model
INTEGER :: H2_INDEX=0                      !< Index for H2 in SIMPLE_CHEMISTRY model
INTEGER :: SOOT_INDEX=0                    !< Index for SOOT in SIMPLE_CHEMISTRY model
INTEGER :: H2O_SMIX_INDEX = -1             !< Index for H2O
INTEGER :: HCN_INDEX=0                     !< Index for HCN
INTEGER :: NO_INDEX=0                      !< Index for NO
INTEGER :: NO2_INDEX=0                     !< Index for NO2
INTEGER :: ZETA_INDEX=0                    !< Index for unmixed fuel fraction, ZETA
INTEGER :: MOISTURE_INDEX=0                !< Index for MATL MOISTURE
INTEGER :: CHAR_INDEX=0                    !< Index for MATL CHAR

INTEGER :: STOP_STATUS=NO_STOP             !< Indicator of whether and why to stop the job
INTEGER :: INPUT_FILE_LINE_NUMBER=0        !< Indicator of what line in the input file is being read

INTEGER :: RND_SEED=0                      !< User RANDOM_SEED

! Miscellaneous logical constants

LOGICAL :: HVAC_DEBUG=.FALSE.               !< Output known hvac values to smokeview
LOGICAL :: RADIATION=.TRUE.                 !< Perform radiation transport
LOGICAL :: UPDATE_ALL_ANGLES=.FALSE.        !< Update all radiation angles the next time the solver is called
LOGICAL :: INCLUDE_PYROLYSIS=.FALSE.        !< Solid phase pyrolysis is included in the simulation
LOGICAL :: EXCHANGE_RADIATION=.FALSE.       !< Do an MPI radiation exchange at this time step
LOGICAL :: EXCHANGE_OBST_MASS=.FALSE.       !< Exchange mass loss information for obstructions bordering interpolated meshes
LOGICAL :: CYLINDRICAL=.FALSE.              !< Cylindrical domain option
LOGICAL :: NOISE=.TRUE.                     !< Initialize velocity field with a small amount of divergence-free motion
LOGICAL :: PREDICTOR                        !< The first half of the second-order accurate time-step
LOGICAL :: CORRECTOR                        !< The second half of the second-order accurate time-step
LOGICAL :: INITIALIZATION_PHASE=.TRUE.      !< The set-up phase before the time-stepping loop
LOGICAL :: APPEND=.FALSE.                   !< For a RESTARTed calculation, APPEND the exising output files
LOGICAL :: PARTICLE_FILE=.FALSE.            !< Indicates the existence of Lagrangian particles
LOGICAL :: PARTICLE_DRAG=.FALSE.            !< Indicates there are particles that drag the gas
LOGICAL :: RESTART=.FALSE.                  !< Indicates if a former calculation is to be RESTARTed
LOGICAL :: SUPPRESSION=.TRUE.               !< Indicates if gas-phase combustion extinction is modeled
LOGICAL :: ACCUMULATE_WATER=.FALSE.         !< Indicates that integrated liquid outputs are specified
LOGICAL :: WRITE_XYZ=.FALSE.                !< Indicates that a Plot3D geometry file is specified by user
LOGICAL :: CHECK_POISSON=.FALSE.            !< Check the accuracy of the Poisson solver
LOGICAL :: WRITE_PARCSRPCG_MATRIX=.FALSE.   !< If true, write out matrix for UGLMAT HYPRE solver
LOGICAL :: TWO_D=.FALSE.                    !< Perform a 2-D simulation
LOGICAL :: SETUP_ONLY=.FALSE.               !< Indicates that the calculation should be stopped before time-stepping
LOGICAL :: CHECK_MESH_ALIGNMENT=.FALSE.     !< Indicates that the user wants to check the mesh alignment and then stop
LOGICAL :: SMOKE3D=.TRUE.                   !< Indicates that the 3D smoke and fire output is desired
LOGICAL :: SMV_PARALLEL_WRITE=.FALSE.       !< If true, the CHID.smv file is written in parallel using MPI-IO.
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
LOGICAL :: FREEZE_VELOCITY=.FALSE.          !< Hold velocity fixed, do not perform a velocity update
LOGICAL :: BNDF_DEFAULT=.TRUE.              !< Output boundary output files
LOGICAL :: SPATIAL_GRAVITY_VARIATION=.FALSE.!< Assume gravity varies as a function of the \f$ x \f$ coordinate
LOGICAL :: CHECK_VN=.TRUE.                  !< Check the Von Neumann number
LOGICAL :: CHECK_FO=.FALSE.                 !< Check the solid phase Fourier number
LOGICAL :: LIQUID_DROPLETS=.FALSE.          !< Indicates the existence of liquid droplets
LOGICAL :: SOLID_PARTICLES=.FALSE.          !< Indicates the existence of solid particles
LOGICAL :: ORIENTED_PARTICLES=.FALSE.       !< Indicates the existence of particles with a specified orientation
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
LOGICAL :: CHECK_HT=.FALSE.                 !< Apply heat transfer stability condition
LOGICAL :: PATCH_VELOCITY=.FALSE.           !< Assume user-defined velocity patches
LOGICAL :: OVERWRITE=.TRUE.                 !< Overwrite old output files
LOGICAL :: INIT_HRRPUV=.FALSE.              !< Assume an initial spatial distribution of HRR per unit volume
LOGICAL :: SYNTHETIC_EDDY_METHOD=.FALSE.
LOGICAL :: UVW_RESTART=.FALSE.              !< Initialize velocity field with values from a file
LOGICAL :: TMP_RESTART=.FALSE.              !< Initialize temperature field with values from a file
LOGICAL :: SPEC_RESTART=.FALSE.             !< Initialize tracked species field with values from a file
LOGICAL :: PARTICLE_CFL=.FALSE.             !< Include particle velocity as a constraint on time step
LOGICAL :: RTE_SOURCE_CORRECTION=.TRUE.     !< Apply a correction to the radiation source term to achieve desired rad fraction
LOGICAL :: OBST_CREATED_OR_REMOVED=.FALSE.  !< An obstruction has just been created or removed and wall cells must be reassigned
LOGICAL :: CHECK_REALIZABILITY=.FALSE.
LOGICAL :: MIN_DEVICES_EXIST=.FALSE.
LOGICAL :: MAX_DEVICES_EXIST=.FALSE.
LOGICAL :: SUPPRESS_DIAGNOSTICS=.FALSE.     !< Do not print detailed mesh-specific output in the .out file
LOGICAL :: WRITE_GEOM_FIRST=.TRUE.
LOGICAL :: SIMPLE_CHEMISTRY=.FALSE.         !< Use simple chemistry combustion model
LOGICAL :: FIRST_PASS                       !< The point in the time step before the CFL constraint is applied
LOGICAL :: VERBOSE=.FALSE.                  !< Add extra output in the .err file
LOGICAL :: SOLID_HEAT_TRANSFER_3D=.FALSE.
LOGICAL :: HVAC_MASS_TRANSPORT=.FALSE.
LOGICAL :: DUCT_HT=.FALSE.
LOGICAL :: DUCT_HT_INSERTED=.FALSE.
LOGICAL :: HVAC_QFAN=.FALSE.
LOGICAL :: USE_ATMOSPHERIC_INTERPOLATION=.FALSE.
LOGICAL :: NEAR_WALL_PARTICLE_INTERPOLATION=.FALSE.
LOGICAL :: POSITIVE_ERROR_TEST=.FALSE.
LOGICAL :: OBST_SHAPE_AREA_ADJUST=.FALSE.
LOGICAL :: STORE_SPECIES_FLUX=.FALSE.
LOGICAL :: STORE_RADIATION_TERMS=.FALSE.
LOGICAL :: STORE_PRESSURE_POISSON_RESIDUAL=.FALSE.
LOGICAL :: OXIDATION_REACTION=.FALSE.
LOGICAL :: PERIODIC_DOMAIN_X=.FALSE.                !< The domain is periodic \f$ x \f$
LOGICAL :: PERIODIC_DOMAIN_Y=.FALSE.                !< The domain is periodic \f$ y \f$
LOGICAL :: PERIODIC_DOMAIN_Z=.FALSE.                !< The domain is periodic \f$ z \f$
LOGICAL :: OPEN_WIND_BOUNDARY=.FALSE.               !< There is a prevailing wind
LOGICAL :: WRITE_DEVC_CTRL=.FALSE.                  !< Flag for writing DEVC and CTRL logfile
LOGICAL :: INIT_INVOKED_BY_SURF=.FALSE.             !< Flag indicating that a SURF line specifies an INIT line
LOGICAL :: NO_PRESSURE_ZONES=.FALSE.                !< Flag to suppress pressure zones
LOGICAL :: CTRL_DIRECT_FORCE=.FALSE.                !< Allow adjustable direct force via CTRL logic
LOGICAL :: REACTING_THIN_OBSTRUCTIONS=.FALSE.       !< Thin obstructions that off-gas are present
LOGICAL :: TENSOR_DIFFUSIVITY=.FALSE.               !< If true, use experimental tensor diffusivity model for spec and tmp
LOGICAL :: OXPYRO_MODEL=.FALSE.                     !< Flag to use oxidative pyrolysis mass transfer model
LOGICAL :: OUTPUT_WALL_QUANTITIES=.FALSE.           !< Flag to force call to WALL_MODEL
LOGICAL :: STORE_FIRE_ARRIVAL=.FALSE.               !< Flag for tracking arrival of spreading fire front over a surface
LOGICAL :: STORE_FIRE_RESIDENCE=.FALSE.             !< Flag for tracking residence time of spreading fire front over a surface
LOGICAL :: STORE_LS_SPREAD_RATE=.FALSE.             !< Flag for outputting local level set spread rate magnitude
LOGICAL :: TEST_NEW_CHAR_MODEL=.FALSE.              !< Flag to envoke new char model
LOGICAL :: FLUX_LIMITER_MW_CORRECTION=.TRUE.        !< Flag for MW correction ensure consistent equation of state at face

INTEGER, ALLOCATABLE, DIMENSION(:) :: CHANGE_TIME_STEP_INDEX      !< Flag to indicate if a mesh needs to change time step
INTEGER, ALLOCATABLE, DIMENSION(:) :: SETUP_PRESSURE_ZONES_INDEX  !< Flag to indicate if a mesh needs to keep searching for ZONEs
REAL(EB), ALLOCATABLE, DIMENSION(:) :: MAX_CELL_ASPECT_RATIO      !< Max cell aspect ratio for each mesh

! Miscellaneous character strings

CHARACTER(FORMULA_LENGTH) :: TITLE                                        !< Job title from HEAD line
CHARACTER(FORMULA_LENGTH) :: RENDER_FILE                                  !< Third party geometry file for Smokeview
CHARACTER(FORMULA_LENGTH) :: UVW_FILE='null'                              !< Velocity initialization field file
INTEGER :: N_TERRAIN_IMAGES=0                                             !< Number of terrain images
CHARACTER(FORMULA_LENGTH), DIMENSION(MAX_TERRAIN_IMAGES) :: TERRAIN_IMAGE !< Name of external terrain image
CHARACTER(CHID_LENGTH) :: CHID                                            !< Job ID
CHARACTER(CHID_LENGTH) :: RESTART_CHID                                    !< Job ID for a restarted case
CHARACTER(FILE_LENGTH) :: RESULTS_DIR                                     !< Custom directory for output
CHARACTER(FILE_LENGTH) :: BINGEOM_DIR                                     !< Custom directory for writing binary geometry files
CHARACTER(5) :: DECIMAL_SPECIFIER='POINT'                                 !< Use point or comma for real outputs
CHARACTER(1) :: SEPARATOR                                                 !< Decimal point or comma

! Dates, version numbers, revision numbers

REAL(FB) :: VERSION_NUMBER=6.0                                            !< Version number solely for a few Smokeview files
CHARACTER(FORMULA_LENGTH) :: REVISION='unknown'                           !< Git revision hash
CHARACTER(FORMULA_LENGTH) :: REVISION_DATE='unknown'                      !< Date of last code change
CHARACTER(FORMULA_LENGTH) :: COMPILE_DATE='unknown'                       !< Date of last code compilation

! Miscellaneous real constants

REAL(EB) :: CPOPR                              !< Specific heat divided by the Prandtl number (J/kg/K)
REAL(EB) :: RSC_T                              !< Reciprocal of the turbulent Schmidt number
REAL(EB) :: RPR_T                              !< Reciprocal of the turbulent Prandtl number
REAL(EB) :: TMPA                               !< Ambient temperature (K)
REAL(EB) :: TMPA4                              !< Ambient temperature to the fourth power (K^4)
REAL(EB) :: RHOA                               !< Ambient density (kg/m3)
REAL(EB) :: P_INF=101325._EB                   !< Ambient pressure at the ground (Pa)
REAL(EB) :: CP_GAMMA                           !< \f$ \frac{\gamma}{\gamma-1} \frac{R}{W_1} \f$
REAL(EB) :: GAMMA=1.4_EB                       !< Ratio of specific heats, typically 1.4
REAL(EB) :: GM1OG                              !< \f$ (\gamma-1)/\gamma \f$
REAL(EB) :: U0                                 !< Wind speed in the \f$ x \f$ direction (m/s)
REAL(EB) :: V0                                 !< Wind speed in the \f$ y \f$ direction (m/s)
REAL(EB) :: W0                                 !< Wind speed in the \f$ z \f$ direction (m/s)
REAL(EB) :: GVEC(3)                            !< Gravity vector (m/s2)
REAL(EB) :: FVEC(3)=0._EB                      !< Force vector (N/m3)
REAL(EB) :: OVEC(3)=0._EB                      !< Coriolis vector (1/s)
REAL(EB) :: C_SMAGORINSKY=0.2_EB               !< Coefficient in turbulence model
REAL(EB) :: C_DEARDORFF=0.1_EB                 !< Coefficient in turbulence model
REAL(EB) :: C_VREMAN=0.07_EB                   !< Coefficient in turbulence model
REAL(EB) :: C_WALE=0.60_EB                     !< Coefficient in turbulence model
REAL(EB) :: LAPSE_RATE                         !< Temperature change with height (K/m)
REAL(EB) :: TEX_ORI(3)                         !< Origin of the texture map for Smokeview (m)
REAL(EB) :: KAPPA0                             !< Background gas radiative absorption coefficient (1/m)
REAL(EB) :: PR_ONTH                            !< Prandtl number to the 1/3 power
REAL(EB) :: MU_AIR_0=1.8E-5_EB                 !< Dynamic Viscosity of Air at 20 C (kg/m/s)
REAL(EB) :: PR_AIR=0.7_EB                      !< Prandtl number for Air
REAL(EB) :: CFL_MAX=1.0_EB                     !< Upper bound of CFL constraint
REAL(EB) :: CFL_MIN=0.8_EB                     !< Lower bound of CFL constraint
REAL(EB) :: VN_MAX=1.0_EB                      !< Upper bound of von Neumann constraint
REAL(EB) :: VN_MIN=0.8_EB                      !< Lower bound of von Neumann constraint
REAL(EB) :: PR_T                               !< Turbulent Prandtl number
REAL(EB) :: SC_T                               !< Turbulent Schmidt number
REAL(EB) :: GROUND_LEVEL=0._EB                 !< Height of the ground, used for establishing atmospheric profiles (m)
REAL(EB) :: LIMITING_DT_RATIO=1.E-4_EB         !< Ratio of current to initial time step when code is stopped
REAL(EB) :: NOISE_VELOCITY=0.005_EB            !< Velocity of random noise vectors (m/s)
REAL(EB) :: TAU_DEFAULT=1._EB                  !< Default ramp-up time (s)
REAL(EB) :: TAU_CHEM=1.E-5_EB                  !< Smallest reaction mixing time scale (s)
REAL(EB) :: TAU_FLAME=1.E10_EB                 !< Largest reaction mixing time scale (s)
REAL(EB) :: SMOKE_ALBEDO=0.3_EB                !< Parmeter used by Smokeview
REAL(EB) :: Y_WERNER_WENGLE=11.81_EB           !< Limit of y+ in Werner-Wengle model
REAL(EB) :: PARTICLE_CFL_MAX=1.0_EB            !< Upper limit of CFL constraint based on particle velocity
REAL(EB) :: PARTICLE_CFL_MIN=0.8_EB            !< Lower limit of CFL constraint based on particle velocity
REAL(EB) :: GRAV=9.80665_EB                    !< Acceleration of gravity (m/s2)
REAL(EB), ALLOCATABLE, DIMENSION(:) :: H_V_H2O !< Heat of vaporization for water (J/kg)
REAL(EB) :: CHI_R_MIN=0._EB                    !< Lower bound for radiative fraction
REAL(EB) :: CHI_R_MAX=1._EB                    !< Upper bound for radiative fraction
REAL(EB) :: SPHERE_FILM_FACTOR=ONTH            !< Weighting factor used in droplet evaporation algorithm for droplets
REAL(EB) :: PLATE_FILM_FACTOR=ONTH             !< Weighting factor used in droplet evaporation algorithm for walls
REAL(EB) :: ORIGIN_LAT=-1.E6_EB                !< Latitude of terrain map
REAL(EB) :: ORIGIN_LON=-1.E6_EB                !< Longitude of terrain map
REAL(EB) :: NORTH_BEARING=0._EB                !< North bearing for terrain map
REAL(EB) :: LATITUDE=10000._EB                 !< Latitude for geostrophic calculation
REAL(EB) :: GEOSTROPHIC_WIND(2)=0._EB          !< Wind vector (m/s)
REAL(EB) :: DY_MIN_BLOWING=1.E-8_EB            !< Parameter in liquid evaporation algorithm (m)
REAL(EB) :: SOOT_DENSITY=1800._EB              !< Density of solid soot (kg/m3)
REAL(EB) :: HVAC_MASS_TRANSPORT_CELL_L=-1._EB  !< Global cell size for HVAC mass transport (m)

REAL(EB), PARAMETER :: TMPM=273.15_EB                       !< Melting temperature of water, conversion factor (K)
REAL(EB), PARAMETER :: P_STP=101325._EB                     !< Standard pressure (Pa)
REAL(EB), PARAMETER :: R0=8314.472_EB                       !< Gas constant (J/K/kmol)
REAL(EB), PARAMETER :: SIGMA=5.670373E-8_EB                 !< Stefan-Boltzmann constant (W/m2/K4)
REAL(EB), PARAMETER :: K_BOLTZMANN=1.3806488E-23_EB         !< Parameter in soot algorithm
REAL(EB), PARAMETER :: EARTH_OMEGA=7.272205216643040e-05_EB !< Earth rotation rate [radians/s] = 2*pi/(24*3600)
REAL(EB), PARAMETER :: VON_KARMAN_CONSTANT=0.41_EB          !< von Karman constant
REAL(EB), PARAMETER :: BTILDE_ROUGH=8.5_EB                  !< Fully rough B(s+) in wall model

! Parameters associated with parallel mode

INTEGER :: MY_RANK=0                                           !< The MPI process index, starting at 0
INTEGER :: N_MPI_PROCESSES=1                                !< Number of MPI processes
INTEGER :: LOWER_MESH_INDEX=1000000000                      !< Lower bound of meshes controlled by the current MPI process
INTEGER :: UPPER_MESH_INDEX=-1000000000                     !< Upper bound of meshes controlled by the current MPI process
LOGICAL :: PROFILING=.FALSE.
INTEGER, ALLOCATABLE, DIMENSION(:) :: PROCESS               !< The MPI process of the given mesh index
INTEGER, ALLOCATABLE, DIMENSION(:) :: FILE_COUNTER          !< Counter for the number of output files currently opened

TYPE (MPI_COMM), ALLOCATABLE, DIMENSION(:) :: MPI_COMM_NEIGHBORS       !< MPI communicator for the a given mesh and its neighbors
TYPE (MPI_COMM), ALLOCATABLE, DIMENSION(:) :: MPI_COMM_CLOSE_NEIGHBORS !< MPI communicator for the a given mesh and its neighbors
INTEGER, ALLOCATABLE, DIMENSION(:) :: MPI_COMM_NEIGHBORS_ROOT          !< The rank of the given mesh within the MPI communicator
INTEGER, ALLOCATABLE, DIMENSION(:) :: MPI_COMM_CLOSE_NEIGHBORS_ROOT    !< The rank of the given mesh within the MPI communicator
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS,DISPLS,COUNTS_10,DISPLS_10,COUNTS_20,DISPLS_20,COUNTS_VENT,DISPLS_VENT

! Time parameters

REAL(EB) :: DT_INITIAL                                      !< Initial time step size (s)
REAL(EB) :: T_BEGIN                                         !< Beginning time of simulation (s)
REAL(EB) :: T_END                                           !< Ending time of simulation (s)
REAL(EB) :: TIME_SHRINK_FACTOR                              !< Factor to reduce specific heat and total run time
REAL(EB) :: RELAXATION_FACTOR=1._EB                         !< Factor used to relax normal velocity nudging at immersed boundaries
REAL(EB) :: MPI_TIMEOUT=600._EB                             !< Time to wait for MPI messages to be received (s)
REAL(EB) :: DT_END_MINIMUM=TWO_EPSILON_EB                   !< Smallest possible final time step (s)
REAL(EB) :: DT_END_FILL=1.E-6_EB
INTEGER  :: DIAGNOSTICS_INTERVAL                            !< Number of time steps between diagnostic outputs
REAL(EB) :: UNFREEZE_TIME                                   !< Time to unfreeze a simulation

! Combustion parameters

REAL(EB) :: Y_O2_INFTY=0.232378_EB                                  !< Ambient mass fraction of oxygen
REAL(EB) :: Y_CO2_INFTY=0.000595_EB                                 !< Ambient mass fraction of carbon dioxide
REAL(EB) :: Y_H2O_INFTY=0._EB                                       !< Ambient mass fraction of water vapor
REAL(EB) :: MW_AIR=28.84852_EB                                      !< Molecular weight of air (g/mol)
REAL(EB) :: MW_N2                                                   !< Molecular weight of nitrogen (g/mol)
REAL(EB) :: MW_O2                                                   !< Molecular weight of oxygen (g/mol)
REAL(EB) :: MW_CO2                                                  !< Molecular weight of carbon dioxide (g/mol)
REAL(EB) :: MW_H2O                                                  !< Molecular weight of water vapor (g/mol)
REAL(EB) :: MW_CO                                                   !< Molecular weight of carbon monoxide (g/mol)
REAL(EB) :: MW_H2                                                   !< Molecular weight of hydrogen (g/mol)
REAL(EB) :: MW_HCN                                                  !< Molecular weight of hydrogen cyanide (g/mol)
REAL(EB) :: MW_SOOT                                                 !< Molecular weight of soot (g/mol)
REAL(EB) :: VISIBILITY_FACTOR=3._EB                                 !< Parameter in light extinction calculation
REAL(EB) :: EC_LL                                                   !< Extinction Coefficient, Lower Limit (1/m)
REAL(EB) :: ZZ_MIN_GLOBAL=1.E-10_EB                                 !< Minimum lumped species mass fraction
REAL(EB) :: FIXED_MIX_TIME=-1._EB                                   !< User-specified reaction mixing time (s)
REAL(EB) :: INITIAL_UNMIXED_FRACTION=1._EB                          !< Initial amount of mixed air-fuel in combustion chamber
REAL(EB) :: GLOBAL_ODE_REL_ERROR=1.E-4_EB                           !< Error tolerance in Richardson extrapolation
REAL(EB) :: H_F_REFERENCE_TEMPERATURE=25._EB                        !< Heat of formation reference temperature (C->K)
REAL(EB) :: FREE_BURN_TEMPERATURE=600._EB                           !< Temperature above which fuel and oxygen burn freely (C->K)
REAL(EB) :: FINITE_RATE_MIN_TEMP=-273.15                            !< When FR is present, min temp. to compute combustion (C->K)
REAL(FB) :: HRRPUV_MAX_SMV=1200._FB                                 !< Clipping value used by Smokeview (kW/m3)
REAL(FB) :: TEMP_MAX_SMV=2000._FB                                   !< Clipping value used by Smokeview (C)
REAL(FB) :: TEMP_MIN_SMV=20._FB                                     !< Clipping value used by Smokeview (C)

INTEGER :: N_SPECIES=0                                              !< Number of total gas phase primitive species
INTEGER :: N_REACTIONS                                              !< Number of gas phase reactions
INTEGER :: I_WATER=-1                                               !< Index of the 'WATER VAPOR' tracked species
INTEGER :: N_TRACKED_SPECIES=0                                      !< Number of lumped or tracked (computed) gas species
INTEGER :: N_SURFACE_DENSITY_SPECIES=0
INTEGER :: COMBUSTION_ODE_SOLVER=-1                                 !< Indicator of ODE solver
INTEGER :: EXTINCT_MOD=-1                                           !< Indicator of extinction model
INTEGER :: MAX_CHEMISTRY_SUBSTEPS=20                                !< Limit on combustion iterations
INTEGER :: MAX_PRIORITY=1                                           !< Maximum numbers of serial fast reactions
INTEGER :: N_PASSIVE_SCALARS=0                                      !< Number of passive scalars
INTEGER :: N_TOTAL_SCALARS=0                                        !< Number of total scalars, tracked and passive
INTEGER :: N_FIXED_CHEMISTRY_SUBSTEPS=-1                            !< Number of chemistry substeps in combustion routine
INTEGER :: ZETA_0_RAMP_INDEX=0                                      !< Ramp index for initial unmixed fraction

LOGICAL :: OUTPUT_CHEM_IT=.FALSE.
LOGICAL :: REAC_SOURCE_CHECK=.FALSE.
LOGICAL :: COMPUTE_ADIABATIC_FLAME_TEMPERATURE=.FALSE.              !< Report adiabatic flame temperature per REAC in LU_OUTPUT

REAL(EB) :: RSUM0                                     !< Initial specific gas constant, \f$ R \sum_i Z_{i,0}/W_i \f$

INTEGER :: I_MAX_TEMP=5000 !< Maximum dimension in K for temperature arrays
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: Z2Y          !< Matrix that converts lumped species vector to primitive, \f$ AZ=Y \f$
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: CP_Z         !< CP_Z(I,J) Specific heat (J/kg/K) of lumped species J at temperature I (K)
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: CPBAR_Z
!< CPBAR_Z(I,J) Average specific heat (J/kg/K) of lumped species J at temperature I (K). Includes reference enthalpy.
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: K_RSQMW_Z
!< K_RSQMW_Z(I,J) Conductivty (W/m/K) of lumped species J at temperature I (K) divided by SM%MW^0.5
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: MU_RSQMW_Z
!< MU_RSQMW_Z(I,J) Viscosity (m^2/s)  of lumped species J at temperature I (K) divided by SM%MW^0.5
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: D_Z          !< D_Z(I,J) Diffusivity (m^2/s) of lumped species J at temp I (K)
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: G_F_Z        !< CP_Z(I,J) Gibbs free energy (J/kg) of lumped species J at temp I (K)
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: S_Z          !< Entropy (J/K) of lumped species J at temp I (K)
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: H_SENS_Z     !< H_SENS(I,J) Sensible enthalpy (J/kg) of lumped species J at temp I (K)
REAL(EB) :: DZZ_CLIP                                  !< Value for processing DZZ in combustion

REAL(EB), ALLOCATABLE, DIMENSION(:) :: MWR_Z,RSQ_MW_Z
CHARACTER(LABEL_LENGTH) :: EXTINCTION_MODEL='null'

! Radiation parameters

LOGICAL, ALLOCATABLE, DIMENSION(:) :: RADIATION_COMPLETED  !< Indicates that the radiation field is completely updated

INTEGER :: NUMBER_SPECTRAL_BANDS=0                         !< Number of wavelength bands for rad solver (1 for gray gas)
INTEGER :: NUMBER_RADIATION_ANGLES=0                       !< Number of solid angles over which radiation is solved
INTEGER :: ANGLE_INCREMENT=0                               !< Indicates how many radiation angles are updated in one time step
INTEGER :: RADIATION_ITERATIONS=1                          !< Number of times to repeat radiation solve in a single time step
INTEGER :: INITIAL_RADIATION_ITERATIONS                    !< Number of radiation solves before time stepping starts

REAL(EB) :: RTE_SOURCE_CORRECTION_FACTOR=1._EB   !< Multiplicative factor used in correcting RTE source term
REAL(EB) :: RAD_Q_SUM=0._EB   !< \f$ \sum_{ijk} \left( \chi_{\rm r} \dot{q}_{ijk}''' + \kappa_{ijk} U_{ijk} \right) V_{ijk} \f$
REAL(EB) :: KFST4_SUM=0._EB   !< \f$ \sum_{ijk} 4 \kappa_{ijk} \sigma T_{ijk}^4 V_{ijk} \f$
REAL(EB) :: QR_CLIP=1._EB     !< Lower bound of \f$ \chi_{\rm r} \dot{q}_{ijk}''' \f$ below which no source correction is made
REAL(EB) :: C_MAX=100._EB     !< Maximum value of RAD_Q_SUM/KFST4_SUM
REAL(EB) :: C_MIN=0.1_EB      !< Minimum value of RAD_Q_SUM/KFST4_SUM

! Ramping parameters

CHARACTER(LABEL_LENGTH), POINTER, DIMENSION(:) :: RAMP_ID,RAMP_TYPE
INTEGER :: MAX_RAMPS=100,I_RAMP_GX,I_RAMP_GY,I_RAMP_GZ,&
           I_RAMP_PGF_T,I_RAMP_FVX_T,I_RAMP_FVY_T,I_RAMP_FVZ_T,N_RAMP=0,I_RAMP_TMP0_Z=0,I_RAMP_P0_Z=0,&
           I_RAMP_SPEED_T=0,I_RAMP_SPEED_Z=0,I_RAMP_DIRECTION_T=0,I_RAMP_DIRECTION_Z=0,&
           I_RAMP_UX,I_RAMP_UY,I_RAMP_UZ,I_RAMP_VX,I_RAMP_VY,I_RAMP_VZ,I_RAMP_WX,I_RAMP_WY,I_RAMP_WZ
INTEGER, PARAMETER :: MAX_QDOTPP_REF=10                    !< Maximum number of REFERENCE_HEAT_FLUX curves for Spyro
INTEGER, PARAMETER :: TIME_HEAT=-11,TIME_VELO=-2,TIME_TEMP=-3,TIME_EFLUX=-4,TIME_PART=-5,TANH_RAMP=-2,TSQR_RAMP=-1,&
                      VELO_PROF_X=-6,VELO_PROF_Y=-7,VELO_PROF_Z=-8,TIME_TGF=-9,TIME_TGB=-10,TIME_TB=-1,&
                      N_SURF_RAMPS=11+MAX_QDOTPP_REF-1

! TABLe parameters

CHARACTER(LABEL_LENGTH) :: TABLE_ID(1000)
INTEGER :: N_TABLE=0,TABLE_TYPE(1000)
INTEGER, PARAMETER :: SPRAY_PATTERN=1,PART_RADIATIVE_PROPERTY=2,POINTWISE_INSERTION=3

! Variables related to meshes

INTEGER :: NMESHES=1,IBAR_MAX=0,JBAR_MAX=0,KBAR_MAX=0,MESH_LIST_EMB(100)
REAL(EB) :: XS_MIN=1.E6_EB,XF_MAX=-1.E6_EB,YS_MIN=1.E6_EB,YF_MAX=-1.E6_EB,ZS_MIN=1.E6_EB,ZF_MAX=-1.E6_EB, &
            DZS_MAX=-1._EB,DZF_MAX=-1._EB
CHARACTER(LABEL_LENGTH), DIMENSION(:), ALLOCATABLE :: MESH_NAME

! Variables related to pressure solver

LOGICAL :: ITERATE_PRESSURE=.FALSE.                              !< Flag indicating if pressure solution is iterated
LOGICAL :: ITERATE_BAROCLINIC_TERM                               !< Flag indicating if baroclinic term is iterated
LOGICAL :: SUSPEND_PRESSURE_ITERATIONS=.FALSE.                   !< Flag for stopping pressure iterations if solution seems stuck
REAL(EB) :: VELOCITY_TOLERANCE=0._EB                             !< Error tolerance for normal velocity at solids or boundaries
REAL(EB) :: PRESSURE_TOLERANCE=0._EB                             !< Error tolerance for iteration of baroclinic pressure term
REAL(EB) :: ITERATION_SUSPEND_FACTOR=0.95_EB                     !< If new velocity error is not this value of old, stop iteration
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VELOCITY_ERROR_MAX        !< Max velocity error of entire domain
REAL(EB), ALLOCATABLE, DIMENSION(:) :: PRESSURE_ERROR_MAX        !< Max pressure error of entire domain
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: VELOCITY_ERROR_MAX_LOC   !< Indices of max velocity error
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: PRESSURE_ERROR_MAX_LOC   !< Indices of max pressure error
INTEGER :: PRESSURE_ITERATIONS=0                                 !< Counter for pressure iterations
INTEGER :: MAX_PREDICTOR_PRESSURE_ITERATIONS=-1                  !< Max pressure iterations per pressure solve in predictor
INTEGER :: MAX_PRESSURE_ITERATIONS=10                            !< Max pressure iterations per pressure solve
INTEGER :: TOTAL_PRESSURE_ITERATIONS=0                           !< Counter for total pressure iterations
CHARACTER(LABEL_LENGTH) :: PRES_METHOD='FFT'                     !< Pressure solver method
INTEGER, PARAMETER :: FFT_FLAG=0                                 !< Integer pressure solver parameter FFT
INTEGER, PARAMETER :: GLMAT_FLAG=1                               !< Integer pressure solver parameter GLMAT
INTEGER, PARAMETER :: UGLMAT_FLAG=2                              !< Integer pressure solver parameter UGLMAT
INTEGER, PARAMETER :: ULMAT_FLAG=3                               !< Integer pressure solver parameter ULMAT
INTEGER, PARAMETER :: MKL_PARDISO_FLAG=1                         !< Integer matrix solver library flag for MKL PARDISO
INTEGER, PARAMETER :: MKL_CPARDISO_FLAG=1                        !< Integer matrix solver library flag for MKL CLUSTER PARDISO
INTEGER, PARAMETER :: HYPRE_FLAG=2                               !< Integer matrix solver library flag for HYPRE
INTEGER :: ULMAT_SOLVER_LIBRARY=MKL_PARDISO_FLAG                 !< Integer ULMAT library flag (defaults to MKL PARDISO)
INTEGER :: UGLMAT_SOLVER_LIBRARY=MKL_CPARDISO_FLAG               !< Integer UGLMAT library flag (defaults to MKL CPARDISO)
INTEGER :: PRES_FLAG = FFT_FLAG                                  !< Pressure solver
LOGICAL :: TUNNEL_PRECONDITIONER=.FALSE.                         !< Use special pressure preconditioner for tunnels
INTEGER :: TUNNEL_NXP                                            !< Number of x points in the entire tunnel
REAL(EB), ALLOCATABLE, DIMENSION(:) :: TP_AA                     !< Upper off-diagonal of tri-diagonal matrix for tunnel pressure
REAL(EB), ALLOCATABLE, DIMENSION(:) :: TP_BB                     !< Lower off-diagonal of tri-diagonal matrix for tunnel pressure
REAL(EB), ALLOCATABLE, DIMENSION(:) :: TP_CC                     !< Right hand side of 1-D tunnel pressure linear system
REAL(EB), ALLOCATABLE, DIMENSION(:) :: TP_DD                     !< Diagonal of tri-diagonal matrix for tunnel pressure solver
REAL(EB), ALLOCATABLE, DIMENSION(:) :: TP_RDXN                   !< Reciprocal of the distance between tunnel precon points
REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:) :: H_BAR             !< Pressure solution of 1-D tunnel pressure solver, predictor step
INTEGER, ALLOCATABLE, DIMENSION(:) :: COUNTS_TP                  !< Counter for MPI calls used for 1-D tunnel pressure solver
INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS_TP                  !< Displacements for MPI calls used for 1-D tunnel pressure solver
INTEGER, ALLOCATABLE, DIMENSION(:) :: I_OFFSET                   !< Spatial index of tunnel

! Miscellaneous integer constants

INTEGER :: ICYC,NFRAMES,PERIODIC_TEST=0,SIM_MODE=3,TURB_MODEL=0,FISHPAK_BC(3)=-1,&
           STOP_AT_ITER=0,WALL_INCREMENT=2,WALL_COUNTER=0,&
           CLIP_DT_RESTRICTIONS_MAX=5,BNDF_TIME_INTEGRALS=0

LOGICAL  :: UPDATE_DEVICES_AGAIN=.FALSE.

! Miscellaneous mesh dimensions

REAL(EB) :: CHARACTERISTIC_CELL_SIZE=1.E6_EB  !< \f$ \min \left( \delta \xi \, \delta \eta \, \delta \zeta \right)^{1/3} \f$
REAL(EB) :: MESH_SEPARATION_DISTANCE          !< Meshes separated if gap greater than min(0.001,0.05*CHARACTERISTIC_CELL_SIZE) (m)
REAL(EB) :: NEIGHBOR_SEPARATION_DISTANCE=-1.  !< No message passing beyond 5*CHARACTERISTIC_CELL_SIZE (m)
REAL(EB) :: ALIGNMENT_TOLERANCE=0.001_EB      !< Maximum ratio of sizes of abutting grid cells

! Logical units and output file names

INTEGER                              :: LU_ERR=ERROR_UNIT,LU_END=2,LU_GIT=3,LU_SMV=4,LU_INPUT=5,LU_OUTPUT=6,LU_STOP=7,LU_CPU=8,&
                                        LU_CATF=9,LU_RDIR=10,LU_GDIR=11,LU_SETCC=12,LU_BINGEOM=13,LU_PARCSRPCG_MATRIX=14
INTEGER                              :: LU_MASS,LU_HRR,LU_STEPS,LU_NOTREADY,LU_VELOCITY_ERROR,LU_CFL,LU_LINE=-1,LU_CUTCELL, &
                                        LU_CVODE_SUBSTEPS
INTEGER                              :: LU_HISTOGRAM,LU_HVAC
INTEGER                              :: LU_GEOC=-1,LU_TGA,LU_INFO,LU_DEVC_CTRL=-1
INTEGER, ALLOCATABLE, DIMENSION(:)   :: LU_PART,LU_PROF,LU_XYZ,LU_TERRAIN,LU_PL3D,LU_DEVC,LU_STATE,LU_CTRL,LU_CORE,LU_RESTART
INTEGER, ALLOCATABLE, DIMENSION(:)   :: LU_VEG_OUT,LU_GEOM,LU_CFACE_GEOM
INTEGER                              :: LU_GEOM_TRAN
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LU_SLCF,LU_SLCF_GEOM,LU_BNDF,LU_BNDG,LU_ISOF,LU_ISOF2, &
                                        LU_SMOKE3D,LU_RADF
INTEGER                              :: DEVC_COLUMN_LIMIT=254,CTRL_COLUMN_LIMIT=254

CHARACTER(FN_LENGTH) :: FN_INPUT='null',FN_STOP='null',FN_CPU,FN_CFL,FN_OUTPUT='null',FN_MASS,FN_HRR,FN_STEPS,FN_SMV,FN_END, &
                        FN_ERR,FN_NOTREADY,FN_VELOCITY_ERROR,FN_GIT,FN_LINE,FN_HISTOGRAM,FN_CUTCELL,FN_TGA,FN_DEVC_CTRL,FN_HVAC, &
                        FN_CVODE_SUBSTEPS
CHARACTER(FN_LENGTH), ALLOCATABLE, DIMENSION(:) :: FN_PART,FN_PROF,FN_XYZ,FN_TERRAIN,FN_PL3D,FN_DEVC,FN_STATE,FN_CTRL,FN_CORE, &
                                                   FN_RESTART,FN_VEG_OUT,FN_GEOM,FN_CFACE_GEOM
CHARACTER(FN_LENGTH), ALLOCATABLE, DIMENSION(:,:) :: FN_SLCF,FN_SLCF_GEOM,FN_BNDF,FN_BNDG,FN_ISOF,FN_ISOF2,FN_SMOKE3D,FN_RADF

CHARACTER(9) :: FMT_R
CHARACTER(25) :: REAL_LIST
CHARACTER( 9) :: CHAR_LIST
CHARACTER(11) :: INTG_LIST
LOGICAL :: OUT_FILE_OPENED=.FALSE.

! Boundary condition arrays

CHARACTER(LABEL_LENGTH) :: MATL_NAME(1:1000)
INTEGER :: N_SURF,N_SURF_RESERVED,N_MATL,MIRROR_SURF_INDEX,OPEN_SURF_INDEX,INTERPOLATED_SURF_INDEX,DEFAULT_SURF_INDEX=0, &
           INERT_SURF_INDEX=0,PERIODIC_SURF_INDEX,PERIODIC_FLOW_ONLY_SURF_INDEX,HVAC_SURF_INDEX=-1,&
           MASSLESS_TRACER_SURF_INDEX, MASSLESS_TARGET_SURF_INDEX,DROPLET_SURF_INDEX,VEGETATION_SURF_INDEX,NWP_MAX
INTEGER,  ALLOCATABLE, DIMENSION(:) :: CELL_COUNT,CELL_COUNT_INTEGERS,CELL_COUNT_LOGICALS
INTEGER,  ALLOCATABLE, DIMENSION(:) :: EDGE_COUNT

! Divergence Arrays

REAL(EB), ALLOCATABLE, DIMENSION(:) :: DSUM,USUM,PSUM

! Level Set vegetation fire spread

INTEGER :: LEVEL_SET_MODE=0               !< Indicator of the type of level set calculation to be done
LOGICAL :: LEVEL_SET_COUPLED_FIRE=.TRUE.  !< Indicator for fire and wind level set coupling
LOGICAL :: LEVEL_SET_COUPLED_WIND=.TRUE.  !< Indicator for fire and wind level set coupling
LOGICAL :: LEVEL_SET_ELLIPSE=.TRUE.       !< Indicator of Richards elliptical level set formulation
LOGICAL :: LSET_TAN2

! Parameters for Terrain and Wind simulation needs

LOGICAL :: TERRAIN_CASE=.FALSE.
INTEGER :: N_VENT_TOTAL=0

! Sprinkler Variables

REAL(EB) :: C_DIMARZO=6.E6_EB
INTEGER :: N_ACTUATED_SPRINKLERS=0
INTEGER, PARAMETER :: NDC=1000,NDC2=100
INTEGER, PARAMETER :: RM_NO_B        = -1 !< Ranz-Marshall no B number
INTEGER, PARAMETER :: RM_B           =  0 !< Ranz-Marshall with B number
INTEGER, PARAMETER :: RM_LEWIS_B     =  1 !< Ranz-Marshall with Lewis number based B Number
INTEGER, PARAMETER :: RM_FL_LEWIS_B  =  2 !< Ranz-Marshall with flux limited, Lewis number based B Number
LOGICAL :: POROUS_FLOOR=.TRUE.

! Particles

INTEGER :: MAXIMUM_PARTICLES,N_LAGRANGIAN_CLASSES,N_LP_ARRAY_INDICES=0
REAL(EB) :: CNF_CUTOFF=0.005_EB
LOGICAL :: PL3D_PARTICLE_FLUX=.FALSE.,SLCF_PARTICLE_FLUX=.FALSE.,DEVC_PARTICLE_FLUX=.FALSE.
LOGICAL :: MPI_PARTICLE_EXCHANGE=.FALSE.,EXCHANGE_INSERTED_PARTICLES=.FALSE.

INTEGER :: MOMENTUM_INTERPOLATION_METHOD=0

! Soot oxidation

REAL(EB), ALLOCATABLE, DIMENSION(:) :: NU_SOOT_OX

! Agglomeration model

LOGICAL :: AGGLOMERATION = .TRUE.,SOOT_OXIDATION=.FALSE.
INTEGER :: N_AGGLOMERATION_SPECIES=0
INTEGER, ALLOCATABLE, DIMENSION(:) :: N_PARTICLE_BINS,AGGLOMERATION_SPEC_INDEX,AGGLOMERATION_SMIX_INDEX
REAL(EB) :: NUCLEATION_SITES=1.E7_EB !1E7 is 10 nucleation sites per cm^3
REAL(EB), ALLOCATABLE, DIMENSION(:) :: MIN_PARTICLE_DIAMETER,MAX_PARTICLE_DIAMETER

! Number of initial value, pressure zone, and multiplier derived types

INTEGER :: N_INIT,N_ZONE,N_MULT,N_MOVE
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: CONNECTED_ZONES
REAL(EB) :: MINIMUM_ZONE_VOLUME=0._EB
REAL(EB) :: PRESSURE_RELAX_TIME=1._EB

! Clipping values

REAL(EB) :: TMPMIN                              !< Minimum gas phase temperature (K)
REAL(EB) :: TMPMAX                              !< Maximum gas phase temperature (K)
REAL(EB) :: RHOMIN                              !< Minimum gas density (kg/m3)
REAL(EB) :: RHOMAX                              !< Maximum gas density (kg/m3)

! Flux limiter

INTEGER, PARAMETER :: CENTRAL_LIMITER=0,GODUNOV_LIMITER=1,SUPERBEE_LIMITER=2,MINMOD_LIMITER=3,CHARM_LIMITER=4,MP5_LIMITER=5
INTEGER :: I_FLUX_LIMITER=SUPERBEE_LIMITER,CFL_VELOCITY_NORM=-999
LOGICAL :: CFL_VELOCITY_NORM_USER_SPECIFIED=.FALSE.

! Numerical quadrature (used in TEST_FILTER)

INTEGER, PARAMETER :: TRAPEZOID_QUADRATURE=0, SIMPSON_QUADRATURE=1, MIDPOINT_QUADRATURE=2
INTEGER :: TEST_FILTER_QUADRATURE=TRAPEZOID_QUADRATURE

INTEGER, PARAMETER :: N_TIMERS=16                   !< Number of subroutine timers
REAL(EB), ALLOCATABLE, DIMENSION(:) :: T_USED       !< Array of subroutine timings
REAL(EB) :: WALL_CLOCK_START                        !< MPI_WTIME i.e. wall clock time when FDS starts
REAL(EB) :: WALL_CLOCK_START_ITERATIONS=0._EB       !< MPI_WTIME i.e. wall clock time when main iteration loop starts
REAL(EB) :: CPU_TIME_START                          !< CPU_TIME when FDS starts

INTEGER :: OPENMP_AVAILABLE_THREADS = 1         !< OpenMP parameter
LOGICAL :: USE_OPENMP               = .FALSE.   !< OpenMP parameter

INTEGER :: N_FACE=0,N_GEOM=0

LOGICAL :: STORE_CUTCELL_DIVERGENCE = .FALSE.
LOGICAL :: STORE_CARTESIAN_DIVERGENCE=.FALSE.

LOGICAL :: CC_IBM=.FALSE.
LOGICAL :: GLMAT_VERBOSE=.FALSE.
LOGICAL :: PRES_ON_WHOLE_DOMAIN=.TRUE.
LOGICAL :: CC_ONLY_IBEDGES_FLAG=.TRUE.
LOGICAL :: ONE_UNKH_PER_CUTCELL=.FALSE.
LOGICAL :: ONE_CC_PER_CARTESIAN_CELL=.TRUE.

! Threshold factor for volume of cut-cells respect to volume of Cartesian cells:
! Currently used in the thermo div definition of cut-cells.

REAL(EB) :: CCVOL_LINK=0.95_EB
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
INTEGER, PARAMETER :: EDG4 = 4

INTEGER :: MAXIMUM_GEOMETRY_ZVALS= 100, MAXIMUM_GEOMETRY_VOLUS=2400,  &
           MAXIMUM_GEOMETRY_FACES=1000, MAXIMUM_GEOMETRY_VERTS=1000,  &
           MAXIMUM_GEOMETRY_IDS  =1000, MAXIMUM_GEOMETRY_SURFIDS=100, &
           MAXIMUM_POLY_VERTS    =1000

INTEGER, PARAMETER :: EXTERNAL_CFACE=-HUGE(1)

! Allocation increment parameters:

INTEGER, PARAMETER :: CC_ALLOC_DVERT = 10
INTEGER, PARAMETER :: CC_ALLOC_DELEM = 10

INTEGER, PARAMETER :: CC_MAX_WSTRIANG_SEG =  2  !< Up to two ws-triangles related to a segment
INTEGER, PARAMETER :: CC_MAX_WSTRIANG_TRI =  1  !< Up to 1 ws-triangle per BODINT_PLANE triangle

INTEGER, ALLOCATABLE, DIMENSION(:,:) :: N_EDGES_DIM_CC

! HVAC Parameters

INTEGER :: N_DUCTNODES = 0, N_DUCTS = 0, N_FANS = 0, N_FILTERS = 0, N_AIRCOILS = 0,N_NETWORKS=0, N_DUCTRUNS=0,&
           N_CONNECTIVITY_INDICES, N_NODE_VARS, N_DUCT_VARS


INTEGER , ALLOCATABLE, DIMENSION(:) :: DUCT_NE,DUCTNODE_NE,DUCT_DR,DUCTNODE_DR
REAL(EB) :: HVAC_PRES_RELAX=1.0_EB,NODE_Z_MIN,NODE_Z_MAX
LOGICAL :: HVAC_SOLVE=.FALSE.,HVAC_LOCAL_PRESSURE=.TRUE.

REAL(EB), POINTER, DIMENSION(:,:) :: ORIENTATION_VECTOR       !< Global array of orientation vectors
INTEGER, ALLOCATABLE, DIMENSION(:) :: NEAREST_RADIATION_ANGLE !< Index of the rad angle most opposite the given ORIENTATION_VECTOR
REAL(EB), POINTER, DIMENSION(:) :: COS_HALF_VIEW_ANGLE     !< View angle of the given ORIENTATION_VECTOR
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VIEW_ANGLE_FACTOR        !< View angle area ORIENTATION_VECTOR
INTEGER :: N_ORIENTATION_VECTOR                               !< Number of ORIENTATION_VECTORs

INTEGER :: TGA_MESH_INDEX=HUGE(INTEGER_ONE) !< Mesh for the special TGA calculation
INTEGER :: TGA_SURF_INDEX=-100              !< Surface properties to use for special TGA calculation
INTEGER :: TGA_WALL_INDEX=-100              !< Wall index to use for special TGA calculation
INTEGER :: TGA_PARTICLE_INDEX=-100          !< Particle index to use for special TGA calculation
REAL(EB) :: TGA_DT=-1._EB                   !< Time step (s) to use for special TGA calculation
REAL(EB) :: TGA_DUMP=1._EB                  !< Temperature output interval (K), starting at TMPA, to use for special TGA calculation
REAL(EB) :: TGA_HEATING_RATE=5._EB          !< Heat rate (K/min) to use for special TGA calculation
REAL(EB) :: TGA_FINAL_TEMPERATURE=800._EB   !< Final Temperature (C) to use for special TGA calculation
REAL(EB) :: TGA_CONVERSION_FACTOR=1._EB     !< Conversion factor for TGA output
REAL(EB) :: MCC_CONVERSION_FACTOR=1._EB     !< Conversion factor for MCC output
REAL(EB) :: DSC_CONVERSION_FACTOR=1._EB     !< Conversion factor for DSC output

LOGICAL :: IBLANK_SMV=.TRUE.  !< Parameter passed to smokeview (in .smv file) to control generation of blockages

! External file control
CHARACTER(250) :: EXTERNAL_FILENAME='null',EXTERNAL_HEARTBEAT_FILENAME='null'
LOGICAL :: READ_EXTERNAL = .FALSE.,HEARTBEAT_FAIL=.TRUE.
INTEGER :: LU_EXTERNAL,LU_EXTERNAL_HEARTBEAT
REAL(EB) :: DT_EXTERNAL=0._EB, T_EXTERNAL,DT_EXTERNAL_HEARTBEAT=0._EB
REAL(EB), ALLOCATABLE, DIMENSION(:) :: EXTERNAL_RAMP
LOGICAL, ALLOCATABLE, DIMENSION(:) :: EXTERNAL_CTRL

! VENT array
REAL(EB), ALLOCATABLE, DIMENSION(:) :: VENT_TOTAL_AREA  !< Array holding grid-snapped areas for all vents

END MODULE GLOBAL_CONSTANTS


!> \brief Clocks for output file dumps

MODULE OUTPUT_CLOCKS

USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)

INTEGER :: RAMP_BNDF_INDEX=0  !< Ramp index for boundary file time series
INTEGER :: RAMP_CTRL_INDEX=0  !< Ramp index for control file time series
INTEGER :: RAMP_CPU_INDEX =0  !< Ramp index for CPU file time series
INTEGER :: RAMP_DEVC_INDEX=0  !< Ramp index for device file time series
INTEGER :: RAMP_FLSH_INDEX=0  !< Ramp index for flush time series
INTEGER :: RAMP_GEOM_INDEX=0  !< Ramp index for geometry output
INTEGER :: RAMP_HRR_INDEX =0  !< Ramp index for hrr file time series
INTEGER :: RAMP_HVAC_INDEX=0  !< Ramp index for hvac file time series
INTEGER :: RAMP_ISOF_INDEX=0  !< Ramp index for isosurface file time series
INTEGER :: RAMP_MASS_INDEX=0  !< Ramp index for mass file time series
INTEGER :: RAMP_PART_INDEX=0  !< Ramp index for particle file time series
INTEGER :: RAMP_PL3D_INDEX=0  !< Ramp index for Plot3D file time series
INTEGER :: RAMP_PROF_INDEX=0  !< Ramp index for profile file time series
INTEGER :: RAMP_RADF_INDEX=0  !< Ramp index for radiation file time series
INTEGER :: RAMP_RSRT_INDEX=0  !< Ramp index for restart file time series
INTEGER :: RAMP_SLCF_INDEX=0  !< Ramp index for slice file time series
INTEGER :: RAMP_SL3D_INDEX=0  !< Ramp index for 3D slice file time series
INTEGER :: RAMP_SM3D_INDEX=0  !< Ramp index for smoke3d file time series
INTEGER :: RAMP_SPEC_INDEX=0  !< Ramp index for species file time series
INTEGER :: RAMP_TIME_INDEX=0  !< Ramp index for specified simulation time steps
INTEGER :: RAMP_DT_INDEX=0    !< Ramp index for specified minimum simulation time step
INTEGER :: RAMP_TMP_INDEX =0  !< Ramp index for temperature file time series
INTEGER :: RAMP_UVW_INDEX =0  !< Ramp index for velocity file time series
REAL(EB), ALLOCATABLE, DIMENSION(:) :: BNDF_CLOCK, CPU_CLOCK,CTRL_CLOCK,DEVC_CLOCK,FLSH_CLOCK,GEOM_CLOCK, HRR_CLOCK,HVAC_CLOCK,&
                                       ISOF_CLOCK,MASS_CLOCK,PART_CLOCK,PL3D_CLOCK,PROF_CLOCK,RADF_CLOCK,RSRT_CLOCK,&
                                       SLCF_CLOCK,SL3D_CLOCK,SM3D_CLOCK,UVW_CLOCK ,TMP_CLOCK ,SPEC_CLOCK
INTEGER, ALLOCATABLE, DIMENSION(:) :: BNDF_COUNTER, CPU_COUNTER,CTRL_COUNTER,DEVC_COUNTER,FLSH_COUNTER,GEOM_COUNTER, HRR_COUNTER,&
                                      HVAC_COUNTER,ISOF_COUNTER,MASS_COUNTER,PART_COUNTER,PL3D_COUNTER,PROF_COUNTER,RADF_COUNTER,&
                                      RSRT_COUNTER,SLCF_COUNTER,SL3D_COUNTER,SM3D_COUNTER,UVW_COUNTER ,TMP_COUNTER ,SPEC_COUNTER
REAL(EB) :: TURB_INIT_CLOCK=-1.E10_EB
REAL(EB) :: MMS_TIMER=1.E10_EB
REAL(EB) :: DT_SLCF,DT_BNDF,DT_DEVC,DT_PL3D,DT_PART,DT_RESTART,DT_ISOF,DT_HRR,DT_HVAC,DT_MASS,DT_PROF,DT_CTRL,&
            DT_FLUSH,DT_SL3D,DT_GEOM,DT_CPU,DT_RADF,DT_SMOKE3D,DT_UVW ,DT_TMP,DT_SPEC
REAL(EB) :: DT_SLCF_SPECIFIED =-1._EB,DT_BNDF_SPECIFIED   =-1._EB,DT_DEVC_SPECIFIED=-1._EB,DT_PL3D_SPECIFIED=-1._EB,&
            DT_PART_SPECIFIED =-1._EB,DT_RESTART_SPECIFIED=-1._EB,DT_ISOF_SPECIFIED=-1._EB,DT_HRR_SPECIFIED =-1._EB,&
            DT_HVAC_SPECIFIED =-1._EB,DT_MASS_SPECIFIED   =-1._EB,DT_PROF_SPECIFIED=-1._EB,DT_CTRL_SPECIFIED=-1._EB,&
            DT_FLUSH_SPECIFIED=-1._EB,DT_SL3D_SPECIFIED   =-1._EB,DT_GEOM_SPECIFIED=-1._EB,DT_CPU_SPECIFIED =-1._EB,&
            DT_RADF_SPECIFIED =-1._EB,DT_SMOKE3D_SPECIFIED=-1._EB,DT_UVW_SPECIFIED =-1._EB,DT_TMP_SPECIFIED =-1._EB,&
            DT_SPEC_SPECIFIED =-1._EB

END MODULE OUTPUT_CLOCKS


!> \brief Radiation parameters

MODULE RADCONS

USE PRECISION_PARAMETERS
IMPLICIT NONE (TYPE,EXTERNAL)

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: DLN                !< Wall-normal matrix
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: DLANG              !< Angles
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: ORIENTATION_FACTOR !< Fraction of radiation angle corresponding to a particular direction
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: BBFRAC             !< Fraction of blackbody radiation
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: WL_LOW             !< Lower wavelength limit of the spectral band
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: WL_HIGH            !< Upper wavelength limit of the spectral band
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: DLX                !< Mean x-component of the control angle vector
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: DLY                !< Mean y-component of the control angle vector
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: DLZ                !< Mean z-component of the control angle vector
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: DLB                !< Mean bottom component of RAYN vector (cylindrical case)
REAL(EB), ALLOCATABLE, DIMENSION(:)   :: RSA                !< Array of solid angles

INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: DLM                !< Mirroring indices
INTEGER, ALLOCATABLE, DIMENSION(:)    :: NRP                !< Number of radiation phi angles at each theta band

REAL(EB) :: RADTMP                !< Radiation temperature (K) for absorption properties (Mie)
REAL(EB) :: PATH_LENGTH           !< Mean path length for the gray gas absorption coefficient (m)
REAL(EB) :: DGROUP_A
REAL(EB) :: DGROUP_B
REAL(EB) :: WEIGH_CYL             !< Weight factor for cylindrical coordinates where all intensities represent 2 angles
REAL(EB) :: DPHI0                 !< Opening angle of the cylindrical domain
REAL(EB) :: FOUR_SIGMA            !< \f$ 4\sigma \f$
REAL(EB) :: RPI_SIGMA             !< \f$ \sigma/\pi \f$
REAL(EB) :: LTSTEP                !< Maximum LAMBDA*T/NLAMBDAT
REAL(EB) :: RTMPMAX               !< Maximum temperature (K) for tabulation of radiative properties
REAL(EB) :: RTMPMIN               !< Minimum temperature (K) for tabulation of radiative properties
REAL(EB) :: MIE_MINIMUM_DIAMETER  !< Minimum droplet size (micron) considered in Mie initialization
REAL(EB) :: MIE_MAXIMUM_DIAMETER  !< Maximum droplet size (micron) considered in Mie initialization

INTEGER :: TIME_STEP_INCREMENT    !< Frequency of calls to radiation solver
INTEGER :: NMIEANG                !< Number of angle bins in forward scattering integration
INTEGER :: NRDMIE                 !< Number of particle radii in Mie calculations
INTEGER :: NLMBDMIE               !< Number of wavelengths in Mie calculations
INTEGER :: MIE_NDG                !< Number of particle radii in WQABS and WQSCA arrays
INTEGER :: NRT                    !< Number of radiation theta angles
INTEGER :: NCO
INTEGER :: UIIDIM
INTEGER :: NLAMBDAT               !< Number of wavelength subdivisions
INTEGER :: N_RADCAL_ARRAY_SIZE
INTEGER :: RADCAL_SPECIES_INDEX(16)
INTEGER :: N_KAPPA_T=44           !< Number of temperature points in absorption coefficient look-up table
INTEGER :: N_KAPPA_Y=50           !< Number of species points in absorption coefficient look-up table

LOGICAL :: WIDE_BAND_MODEL        !< Non-gray gas, wide band model
LOGICAL :: WSGG_MODEL             !< Weighted Sum of Gray Gas model

REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: Z2RADCAL_SPECIES
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: RADCAL_SPECIES2KAPPA
CHARACTER(LABEL_LENGTH) :: RADCAL_SPECIES_ID(16)='NULL'

END MODULE RADCONS

!> \brief Variables for DVODE solver usage
MODULE CHEMCONS
USE PRECISION_PARAMETERS

!> \brief Parameters associated with IGNITION_ZONES
TYPE IGNITION_ZONE_TYPE
   REAL(EB) :: X1            !< Lower x bound of Ignition Zone
   REAL(EB) :: X2            !< Upper x bound of Ignition Zone
   REAL(EB) :: Y1            !< Lower y bound of Ignition Zone
   REAL(EB) :: Y2            !< Upper y bound of Ignition Zone
   REAL(EB) :: Z1            !< Lower z bound of Ignition Zone
   REAL(EB) :: Z2            !< Upper z bound of Ignition Zone
   INTEGER :: DEVC_INDEX=0   !< Index of device controlling the status of the zone
   CHARACTER(LABEL_LENGTH) :: DEVC_ID='null'  !< Name of device controlling the status of the zone
END TYPE IGNITION_ZONE_TYPE

INTEGER, ALLOCATABLE, DIMENSION(:) :: YP2ZZ
REAL(EB) :: ODE_MIN_ATOL= -1._EB
LOGICAL  :: EQUIV_RATIO_CHECK = .FALSE.
REAL(EB) :: MIN_EQUIV_RATIO=0.1_EB
REAL(EB) :: MAX_EQUIV_RATIO=20.0_EB
LOGICAL  :: DO_CHEM_LOAD_BALANCE = .FALSE.
INTEGER  :: MAX_CVODE_SUBSTEPS=100000
INTEGER  :: CVODE_MAX_TRY=4
INTEGER  :: CVODE_ORDER=0
INTEGER  :: CVODE_ERR_CODE_MIN=-100
INTEGER  :: CVODE_ERR_CODE_MAX=100
INTEGER  :: CVODE_WARNING_CELLS(-100:100)! Index of the array is error code and value is Cell count
CHARACTER(LEN=100) :: CVODE_WARN_MESSAGES(-100:100)

! FOR WRITING CVODE SUBSTEPS
LOGICAL  :: WRITE_CVODE_SUBSTEPS = .FALSE.
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: CVODE_SUBSTEP_DATA
INTEGER :: TOTAL_SUBSTEPS_TAKEN

! Adiabatic flame temperature calculation
CHARACTER(LABEL_LENGTH) :: FUEL_ID_FOR_AFT='null'
INTEGER :: I_FUEL,I_CO2,I_H2O,I_O2 ! Store the index of the species in the ZZ array.
LOGICAL  :: USE_MIXED_ZN_AFT_TMP = .FALSE.

! Mixing
REAL(EB) :: ZETA_ARTIFICAL_MIN_LIMIT=0.99_EB
REAL(EB) :: ZETA_ARTIFICAL_MAX_LIMIT=0.9999_EB
REAL(EB) :: ZETA_FIRST_STEP_DIV=10._EB

! IGNITION ZONES (mainly for premixed flame)
INTEGER :: N_IGNITION_ZONES = 0
TYPE(IGNITION_ZONE_TYPE), DIMENSION(MAX_IGNITION_ZONES) :: IGNITION_ZONES !< Coordinates of ignition zones

CONTAINS
   SUBROUTINE INIT_CVODE_WARN_MESSAGES()
      CVODE_WARN_MESSAGES = 'CVODE didn''t finish ODE solution with this code.'
      CVODE_WARN_MESSAGES(-1) = 'CVODE took all internal substeps.'
      CVODE_WARN_MESSAGES(-3) = 'Minimum step size was reached.'
      CVODE_WARN_MESSAGES(-4) = 'Convergence test failure.'
   END SUBROUTINE INIT_CVODE_WARN_MESSAGES

END MODULE CHEMCONS
