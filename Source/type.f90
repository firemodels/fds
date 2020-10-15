!> \brief A collection of major derived types used in FDS.

MODULE TYPES

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS, ONLY : IAXIS,JAXIS,KAXIS,MAX_DIM,LOW_IND,HIGH_IND

IMPLICIT NONE

!> \brief Parameters associated with an entire class of Lagrangian particles

TYPE LAGRANGIAN_PARTICLE_CLASS_TYPE

   CHARACTER(LABEL_LENGTH) :: ID                                               !< Name of particle class
   CHARACTER(LABEL_LENGTH) :: SPEC_ID='null'                                   !< Name of evaporating gas species
   CHARACTER(LABEL_LENGTH) :: DEVC_ID='null'                                   !< Name of controlling device
   CHARACTER(LABEL_LENGTH) :: CTRL_ID='null'                                   !< Name of control function
   CHARACTER(LABEL_LENGTH) :: SURF_ID='null'                                   !< Name of SURFace type
   CHARACTER(LABEL_LENGTH) :: PROP_ID='null'                                   !< Name of PROPerty type
   CHARACTER(LABEL_LENGTH) :: RADIATIVE_PROPERTY_TABLE_ID='null'               !< Name of radiative property table
   CHARACTER(LABEL_LENGTH) :: CNF_RAMP_ID='null'                     !< Cumulative Number Fraction (CNF) function
   CHARACTER(LABEL_LENGTH) :: BREAKUP_CNF_RAMP_ID='null'             !< User-defined cumulative number fraction after break-up
   CHARACTER(LABEL_LENGTH) :: DISTRIBUTION='ROSIN-RAMMLER-LOGNORMAL'           !< Droplet size distribution
   CHARACTER(LABEL_LENGTH) :: BREAKUP_DISTRIBUTION='ROSIN-RAMMLER-LOGNORMAL'   !< Droplet size distribution after break-up
   CHARACTER(LABEL_LENGTH) :: QUANTITIES(10)                                   !< Names of output quantities.
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_LABEL(10)                              !< Smokeview file label for output quantities
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_BAR_LABEL(10)                          !< Short Smokeview label for output quantities
   CHARACTER(LABEL_LENGTH) :: QUANTITIES_SPEC_ID(10)                           !< SPECies IDs for output quantities

   REAL(EB) :: HEAT_OF_COMBUSTION         !< Heat of Combustion (J/kg) of the evaporated gas
   REAL(EB) :: ADJUST_EVAPORATION         !< LPC\%HEAT_OF_COMBUSTION/RN(1)\%HEAT_OF_COMBUSTION
   REAL(EB) :: LIFETIME                   !< Time (s) after insertion when particle is to be removed
   REAL(EB) :: DIAMETER                   !< Median volumetric diameter (m) of the particles
   REAL(EB) :: MINIMUM_DIAMETER           !< Minimum particle diameter (m) in distribution
   REAL(EB) :: MAXIMUM_DIAMETER           !< Maximum particle diameter (m) in distribution
   REAL(EB) :: GAMMA                      !< Parameter in Rosin-Rommler distribution
   REAL(EB) :: KILL_RADIUS                !< Radius (m) below which particle is killed (removed)
   REAL(EB) :: TMP_INITIAL                !< Initial temperature (K) of the particles
   REAL(EB) :: SIGMA                      !< Parameter in Rosin-Rammler distribution
   REAL(EB) :: VERTICAL_VELOCITY          !< Vertical velocity component (m/s) of particle stuck to a solid surface
   REAL(EB) :: HORIZONTAL_VELOCITY        !< Horizontal velocity component (m/s) of particle stuck to a solid surface
   REAL(EB) :: MASS_TRANSFER_COEFFICIENT        !< Mass transfer coefficient from gas to liquid droplet (kg/m2/s)
   REAL(EB) :: HEAT_TRANSFER_COEFFICIENT_GAS    !< Heat transfer coefficient from gas to liquid droplet (W/m2/K)
   REAL(EB) :: HEAT_TRANSFER_COEFFICIENT_SOLID  !< Heat transfer coefficient from solid surface to liquid droplet (W/m2/K)
   REAL(EB) :: DRAG_COEFFICIENT(3)        !< Drag coefficient in 3 coordinate directions
   REAL(EB) :: SURFACE_DIAMETER           !< Effective liquid droplet diameter (m) on a solid surface
   REAL(EB) :: SURFACE_TENSION            !< Surface tension (N/m) of liquid droplets
   REAL(EB) :: BREAKUP_RATIO              !< Ratio of child Sauter mean to parent size in Bag breakup regime
   REAL(EB) :: BREAKUP_GAMMA              !< Rosin-Rammler size distribution parameter for break-up distribution
   REAL(EB) :: BREAKUP_SIGMA              !< Rosin-Rammler size distribution parameter for break-up distribution
   REAL(EB) :: DENSE_VOLUME_FRACTION      !< Limiting volume fraction for drag reduction
   REAL(EB) :: PERMEABILITY(3)            !< Parameter in porous media drag model, \f$K\f$ (\f$ {\rm m}^2 \f$)
   REAL(EB) :: REAL_REFRACTIVE_INDEX      !< Radiative property of liquid droplet
   REAL(EB) :: COMPLEX_REFRACTIVE_INDEX   !< Radiative property of liquid droplet
   REAL(EB) :: DENSITY=-1._EB             !< Density of liquid droplet (kg/m\f$^3\f$)
   REAL(EB) :: FTPR                       !< 4/3 * PI * SPECIES(N)\%DENSITY_LIQUID (kg/m3)
   REAL(EB) :: FREE_AREA_FRACTION         !< Area fraction of cell open for flow in SCREEN_DRAG model
   REAL(EB) :: POROUS_VOLUME_FRACTION     !< Volume fraction of cell open to flow in porous media model
   REAL(EB) :: MEAN_DROPLET_VOLUME=0._EB  !< Mean droplet volume
   REAL(EB) :: RUNNING_AVERAGE_FACTOR     !< Fraction of older value to use for particle statistics summations
   REAL(EB) :: SHAPE_FACTOR               !< Ratio of particle cross sectional area to surface area
   REAL(EB) :: EMBER_DENSITY_THRESHOLD    !< Density at which vegetative particle becomes a flying ember
   REAL(EB) :: EMBER_VELOCITY_THRESHOLD   !< Velocity at which vegetative particle becomes a flying ember
   REAL(EB) :: PRIMARY_BREAKUP_TIME       !< Time (s) after insertion when droplet breaks up
   REAL(EB) :: PRIMARY_BREAKUP_DRAG_REDUCTION_FACTOR   !< Drag reduction factor
   REAL(EB) :: RUNNING_AVERAGE_FACTOR_WALL             !< Fraction of old value used in summations of droplets stuck to walls

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: R_CNF         !< Independent variable (radius) in particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CNF           !< Cumulative Number Fraction particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CVF           !< Cumulative Volume Fraction particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: BREAKUP_R_CNF !< R_CNF of new distribution after particle break-up
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: BREAKUP_CNF   !< CNF of new distribution after particle break-up
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: BREAKUP_CVF   !< CVF of new distribution after particle break-up
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: W_CNF         !< Weighting factor in particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: R50           !< Array of median particle diameters for Mie calculation
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: SOLID_ANGLE   !< Array of solid angles for particle with multiple orientations

   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: WQABS       !< Absorption efficiency factor array
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: WQSCA       !< Scattering efficiency factor array

   INTEGER :: SAMPLING_FACTOR             !< Reduce particle output by this factor
   INTEGER :: N_QUANTITIES                !< Number of output quantities for this class of particles
   INTEGER :: QUANTITIES_INDEX(10)        !< Array of indices of output quantities for this class of particles
   INTEGER :: QUANTITIES_Y_INDEX(10)=-1   !< Array of species indices associated with the output quantities
   INTEGER :: QUANTITIES_Z_INDEX(10)=-1   !< Array of species mixture indices associated with the output quantities
   INTEGER :: ARRAY_INDEX=0               !< Array of indices corresponding to solid or liquid particles
   INTEGER :: RGB(3)                      !< Color indices for default particle class color in Smokeview
   INTEGER :: RADIATIVE_PROPERTY_INDEX=0  !< Index for this class of particles in radiative property table
   INTEGER :: SURF_INDEX=0                !< Surface properties for solid particle
   INTEGER :: DRAG_LAW=1                  !< Code indicating type of drag law
   INTEGER :: DEVC_INDEX=0                !< Index of device that governs this class of particles
   INTEGER :: CTRL_INDEX=0                !< Index of controller that governs this class of particles
   INTEGER :: ORIENTATION_INDEX=0         !< Starting position of the particle class orientation vector within the master array
   INTEGER :: N_ORIENTATION               !< Number of orientations (directions) corresponding to this class of particles
   INTEGER :: Z_INDEX=-1                  !< Species mixture index for this class
   INTEGER :: Y_INDEX=-1                  !< Species index for this class
   INTEGER :: N_STRATA=6                  !< Number of bins in subdivision of size distribution
   INTEGER :: NEAREST_RAD_ANGLE_INDEX=0   !< Index of the radiation angle nearest the given orientation vector
   INTEGER :: CNF_RAMP_INDEX=-1           !< Ramp index for Cumulative Number Fraction function
   INTEGER :: BREAKUP_CNF_RAMP_INDEX=-1   !< Ramp index for break-up Cumulative Number Fraction function
   INTEGER :: N_STORAGE_REALS             !< Number of reals to store for this particle class
   INTEGER :: N_STORAGE_INTEGERS          !< Number of integers to store for this particle class
   INTEGER :: N_STORAGE_LOGICALS          !< Number of logicals to store for this particle class

   INTEGER, ALLOCATABLE, DIMENSION(:) :: STRATUM_INDEX_LOWER  !< Lower index of size distribution band
   INTEGER, ALLOCATABLE, DIMENSION(:) :: STRATUM_INDEX_UPPER  !< Upper index of size distribution band

   LOGICAL :: STATIC=.FALSE.                !< Flag indicating if particles move or not
   LOGICAL :: MASSLESS_TRACER=.FALSE.       !< Flag indicating if particles are just tracers for visualization
   LOGICAL :: MASSLESS_TARGET=.FALSE.       !< Flag indicating if particles are just targets for an output quantity
   LOGICAL :: LIQUID_DROPLET=.FALSE.        !< Flag indicating if particles are liquid droplets
   LOGICAL :: SOLID_PARTICLE=.FALSE.        !< Flag indicating if particles are solid, not liquid
   LOGICAL :: MONODISPERSE=.FALSE.          !< Flag indicating if particle size is monodisperse
   LOGICAL :: TURBULENT_DISPERSION=.FALSE.  !< Flag indicating if subgrid-scale turbulence is applied
   LOGICAL :: BREAKUP=.FALSE.               !< Flag indicating if paricles or droplets break-up
   LOGICAL :: CHECK_DISTRIBUTION=.FALSE.    !< Flag indicating if diagnostic output on size distribution is specified
   LOGICAL :: FUEL=.FALSE.                  !< Flag indicating if droplets evaporate into fuel gas
   LOGICAL :: DUCT_PARTICLE=.FALSE.         !< Flag indicating if particles can pass through a duct
   LOGICAL :: EMBER_PARTICLE=.FALSE.        !< Flag indicating if particles can become flying embers
   LOGICAL :: ADHERE_TO_SOLID=.FALSE.       !< Flag indicating if particles can stick to a solid

END TYPE LAGRANGIAN_PARTICLE_CLASS_TYPE

TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: LAGRANGIAN_PARTICLE_CLASS


!> \brief Solid material density for 1-D pyrolysis/conduction algorithm

TYPE MATL_COMP_TYPE
   REAL(EB), POINTER, DIMENSION(:) :: RHO !< (1:NWP) Solid density (kg/m3)
   REAL(EB), POINTER, DIMENSION(:) :: RHO_DOT !< (1:NWP) Change in solid density (kg/m3/s)
END TYPE MATL_COMP_TYPE

!> \brief Gas mass concentration in solid for 1-D mass transfer
TYPE SPEC_COMP_TYPE
   REAL(EB), POINTER, DIMENSION(:) :: RHO !< (1:NWP) Gas density (kg/m3)
!   REAL(EB), POINTER, DIMENSION(:) :: RHO_DOT !< (1:NWP) Change in gas density (kg/m3/s)
END TYPE SPEC_COMP_TYPE


!> \brief Radiation intensity at a boundary for a given wavelength band

TYPE BAND_TYPE
   REAL(EB), POINTER, DIMENSION(:) :: ILW !< (1:NRA) Radiation intensity (W/m2/sr)
END TYPE BAND_TYPE

! Note: If you change the number of scalar variables in ONE_D_M_AND_E_XFER_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_ONE_D_SCALAR_REALS=36,N_ONE_D_SCALAR_INTEGERS=12,N_ONE_D_SCALAR_LOGICALS=1

!> \brief Variables associated with a WALL, PARTICLE, or CFACE boundary cell

TYPE ONE_D_M_AND_E_XFER_TYPE

   REAL(EB), POINTER, DIMENSION(:) :: TMP                 !< Temperature in center of each solid cell, \f$ T_{{\rm s},i} \f$
   REAL(EB), POINTER, DIMENSION(:) :: LAYER_THICKNESS     !< (1:SF\%N_LAYERS) Thickness of layer (m)
   REAL(EB), POINTER, DIMENSION(:) :: X                   !< (0:NWP) Depth (m), \f$ x_{{\rm s},i} \f$
   REAL(EB), POINTER, DIMENSION(:) :: M_DOT_G_PP_ACTUAL   !< (1:N_TRACKED_SPECIES) Actual mass production rate per unit area
   REAL(EB), POINTER, DIMENSION(:) :: M_DOT_S_PP          !< (1:SF\%N_MATL) Mass production rate of solid species
   REAL(EB), POINTER, DIMENSION(:) :: M_DOT_G_PP_ADJUST   !< (1:N_TRACKED_SPECIES) Adjusted mass production rate per unit area
   REAL(EB), POINTER, DIMENSION(:) :: IL                  !< (1:NSB) Radiance (W/m2/sr); output only
   REAL(EB), POINTER, DIMENSION(:) :: ZZ_G                !< (1:N_TRACKED_SPECIES) Species mixture mass fraction in gas grid cell
   REAL(EB), POINTER, DIMENSION(:) :: ZZ_F                !< (1:N_TRACKED_SPECIES) Species mixture mass fraction at surface
   REAL(EB), POINTER, DIMENSION(:) :: RHO_D_F             !< (1:N_TRACKED_SPECIES) Diffusion at surface, \f$ \rho D_\alpha \f$
   REAL(EB), POINTER, DIMENSION(:) :: RHO_D_DZDN_F        !< \f$ \rho D_\alpha \partial Z_\alpha / \partial n \f$
   REAL(EB), POINTER, DIMENSION(:) :: A_LP_MPUA           !< Accumulated liquid droplet mass per unit area (kg/m2)
   REAL(EB), POINTER, DIMENSION(:) :: AWM_AEROSOL         !< Accumulated aerosol mass per unit area (kg/m2)
   REAL(EB), POINTER, DIMENSION(:) :: LP_CPUA             !< Liquid droplet cooling rate unit area (W/m2)
   REAL(EB), POINTER, DIMENSION(:) :: LP_MPUA             !< Liquid droplet mass per unit area (kg/m2)
   REAL(EB), POINTER, DIMENSION(:) :: LP_TEMP             !< Liquid droplet mean temperature (K)
   REAL(EB), POINTER, DIMENSION(:) :: RHO_C_S             !< Solid density times specific heat (J/m3/K)
   REAL(EB), POINTER, DIMENSION(:) :: K_S                 !< Solid conductivity (W/m/K)

   TYPE(MATL_COMP_TYPE), ALLOCATABLE, DIMENSION(:) :: MATL_COMP !< (1:SF\%N_MATL) Material component
   TYPE(SPEC_COMP_TYPE), ALLOCATABLE, DIMENSION(:) :: SPEC_COMP !< (1:SF\%N_SPEC) Gas component
   TYPE(BAND_TYPE), ALLOCATABLE, DIMENSION(:) :: BAND           !< 1:NSB) Radiation wavelength band
   INTEGER, POINTER, DIMENSION(:) :: N_LAYER_CELLS              !< (1:SF\%N_LAYERS) Number of cells in the layer

   INTEGER, POINTER :: ARRAY_INDEX    !< WALL, LAGRANGIAN_PARTICLE, or CFACE index
   INTEGER, POINTER :: STORAGE_INDEX  !< Index in the WALL, LP, or CFACE storate array
   INTEGER, POINTER :: II             !< Ghost cell \f$ x \f$ index
   INTEGER, POINTER :: JJ             !< Ghost cell \f$ y \f$ index
   INTEGER, POINTER :: KK             !< Ghost cell \f$ z \f$ index
   INTEGER, POINTER :: IIG            !< Gas cell \f$ x \f$ index
   INTEGER, POINTER :: JJG            !< Gas cell \f$ y \f$ index
   INTEGER, POINTER :: KKG            !< Gas cell \f$ z \f$ index
   INTEGER, POINTER :: IOR            !< Index of orientation of the WALL cell
   INTEGER, POINTER :: PRESSURE_ZONE  !< Pressure ZONE of the adjacent gas phase cell
   INTEGER, POINTER :: NODE_INDEX     !< HVAC node index associated with surface
   INTEGER, POINTER :: N_SUBSTEPS     !< Number of substeps in the 1-D conduction/reaction update

   REAL(EB), POINTER :: AREA            !< Face area (m2)
   REAL(EB), POINTER :: HEAT_TRANS_COEF !< Heat transfer coefficient (W/m2/K)
   REAL(EB), POINTER :: Q_CON_F         !< Convective heat flux at surface (W/m2)
   REAL(EB), POINTER :: Q_RAD_IN        !< Incoming radiative flux (W/m2)
   REAL(EB), POINTER :: Q_RAD_OUT       !< Outgoing radiative flux (W/m2)
   REAL(EB), POINTER :: EMISSIVITY      !< Surface emissivity
   REAL(EB), POINTER :: AREA_ADJUST     !< Ratio of actual surface area to grid cell face area
   REAL(EB), POINTER :: T_IGN           !< Ignition time (s)
   REAL(EB), POINTER :: TMP_F           !< Surface temperature (K)
   REAL(EB), POINTER :: TMP_F_OLD       !< Holding value for surface temperature (K)
   REAL(EB), POINTER :: TMP_B           !< Back surface temperature (K)
   REAL(EB), POINTER :: U_NORMAL        !< Normal component of velocity (m/s) at surface, start of time step
   REAL(EB), POINTER :: U_NORMAL_S      !< Estimated normal component of velocity (m/s) at next time step
   REAL(EB), POINTER :: U_NORMAL_0      !< Initial or specified normal component of velocity (m/s) at surface
   REAL(EB), POINTER :: RSUM_G          !< \f$ R_0 \sum_\alpha Z_\alpha/W_\alpha \f$ in first gas phase cell
   REAL(EB), POINTER :: TMP_G           !< Temperature (K) in adjacent gas phase cell
   REAL(EB), POINTER :: RHO_G           !< Gas density (kg/m3) in adjacent gas phase cell
   REAL(EB), POINTER :: U_TANG          !< Tangential velocity (m/s) near surface
   REAL(EB), POINTER :: RHO_F           !< Gas density at the wall (kg/m3)
   REAL(EB), POINTER :: RDN             !< \f$ 1/ \delta n \f$ at the surface (1/m)
   REAL(EB), POINTER :: MU_G            !< Viscosity, \f$ \mu \f$, in adjacent gas phase cell
   REAL(EB), POINTER :: K_G             !< Thermal conductivity, \f$ k \f$, in adjacent gas phase cell
   REAL(EB), POINTER :: U_TAU           !< Friction velocity (m/s)
   REAL(EB), POINTER :: Y_PLUS          !< Dimensionless boundary layer thickness unit
   REAL(EB), POINTER :: Z_STAR          !< Dimensionless boundary layer unit
   REAL(EB), POINTER :: PHI_LS          !< Level Set value for output only
   REAL(EB), POINTER :: WORK1           !< Work array
   REAL(EB), POINTER :: WORK2           !< Work array
   REAL(EB), POINTER :: Q_DOT_G_PP      !< Heat release rate per unit area (W/m2)
   REAL(EB), POINTER :: Q_DOT_O2_PP     !< Heat release rate per unit area (W/m2) due to oxygen consumption
   REAL(EB), POINTER :: Q_CONDENSE      !< Heat release rate per unit area (W/m2) due to gas condensation
   REAL(EB), POINTER :: K_SUPPRESSION   !< Suppression coefficent (m2/kg/s)
   REAL(EB), POINTER :: BURN_DURATION   !< Duration of a specified fire (s)
   REAL(EB), POINTER :: T_SCALE         !< Scaled time for a surface with CONE_HEAT_FLUX (s)
   REAL(EB), POINTER :: Q_SCALE         !< Scaled integrated heat release for a surface with CONE_HEAT_FLUX
   REAL(EB), POINTER :: L_OBUKHOV       !< Obukhov length (m)

   LOGICAL, POINTER :: BURNAWAY         !< Indicater if cell can burn away when fuel is exhausted

END TYPE ONE_D_M_AND_E_XFER_TYPE

! Note: If you change the number of scalar variables in LAGRANGIAN_PARTICLE_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_PARTICLE_SCALAR_REALS=17,N_PARTICLE_SCALAR_INTEGERS=10,N_PARTICLE_SCALAR_LOGICALS=4

!> \brief Variables assoicated with a single Lagrangian particle

TYPE LAGRANGIAN_PARTICLE_TYPE

   TYPE (ONE_D_M_AND_E_XFER_TYPE) :: ONE_D  !< Most of the particle properties are contained within this derived type

   LOGICAL, POINTER :: SHOW                 !< Show the particle in Smokeview
   LOGICAL, POINTER :: SPLAT                !< The liquid droplet has hit a solid
   LOGICAL, POINTER :: EMBER                !< The particle can break away and become a burning ember
   LOGICAL, POINTER :: PATH_PARTICLE

   REAL(EB), POINTER :: X                   !< \f$ x \f$ coordinate of particle (m)
   REAL(EB), POINTER :: Y                   !< \f$ y \f$ coordinate of particle (m)
   REAL(EB), POINTER :: Z                   !< \f$ z \f$ coordinate of particle (m)
   REAL(EB), POINTER :: U                   !< \f$ x \f$ velocity component of particle (m/s)
   REAL(EB), POINTER :: V                   !< \f$ y \f$ velocity component of particle (m/s)
   REAL(EB), POINTER :: W                   !< \f$ z \f$ velocity component of particle (m/s)
   REAL(EB), POINTER :: PWT                 !< Weight factor of particle; i.e. the number of real particles it represents
   REAL(EB), POINTER :: ACCEL_X             !< Contribution to acceleration of gas in \f$ x \f$ direction (m/s2)
   REAL(EB), POINTER :: ACCEL_Y             !< Contribution to acceleration of gas in \f$ y \f$ direction (m/s2)
   REAL(EB), POINTER :: ACCEL_Z             !< Contribution to acceleration of gas in \f$ z \f$ direction (m/s2)
   REAL(EB), POINTER :: RE                  !< Reynolds number based on particle diameter
   REAL(EB), POINTER :: MASS                !< Particle mass (kg)
   REAL(EB), POINTER :: T_INSERT            !< Time when particle was inserted (s)
   REAL(EB), POINTER :: DX                  !< Length scale used in POROUS_DRAG calculation (m)
   REAL(EB), POINTER :: DY                  !< Length scale used in POROUS_DRAG calculation (m)
   REAL(EB), POINTER :: DZ                  !< Length scale used in POROUS_DRAG calculation (m)
   REAL(EB), POINTER :: M_DOT               !< Particle mass evaporation rate (kg/s)

   INTEGER, POINTER :: TAG                  !< Unique integer identifier for the particle
   INTEGER, POINTER :: ARRAY_INDEX          !< Index in the array of evaporating particles
   INTEGER, POINTER :: STORAGE_INDEX        !< Index in the large storage array of all particles
   INTEGER, POINTER :: CLASS_INDEX          !< LAGRANGIAN_PARTICLE_CLASS of particle
   INTEGER, POINTER :: ORIENTATION_INDEX    !< Index in the array of all ORIENTATIONs
   INTEGER, POINTER :: WALL_INDEX           !< If liquid droplet has stuck to a wall, this is the WALL cell index
   INTEGER, POINTER :: DUCT_INDEX           !< Index of duct
   INTEGER, POINTER :: INIT_INDEX           !< Index of INIT line
   INTEGER, POINTER :: DUCT_CELL_INDEX      !< Index of duct cell
   INTEGER, POINTER :: CFACE_INDEX          !< Index of immersed boundary CFACE that the droplet has attached to

END TYPE LAGRANGIAN_PARTICLE_TYPE


TYPE STORAGE_TYPE
   INTEGER :: N_STORAGE_SLOTS=0,NEXT_AVAILABLE_SLOT=1
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: REALS
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INTEGERS
   LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LOGICALS
END TYPE STORAGE_TYPE


! Note: If you change the number of scalar variables in WALL_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_WALL_SCALAR_REALS=7,N_WALL_SCALAR_INTEGERS=14,N_WALL_SCALAR_LOGICALS=0

!> \brief Variables associated with a WALL cell

TYPE WALL_TYPE

   TYPE (ONE_D_M_AND_E_XFER_TYPE) :: ONE_D     !< Derived type carrying most of the solid boundary conditions

   REAL(EB), POINTER :: DUNDT                  !< \f$ \partial u_n / \partial t \f$
   REAL(EB), POINTER :: Q_LEAK                 !< Heat production of leaking gas (W/m3)
   REAL(EB), POINTER :: V_DEP                  !< Deposition velocity (m/s)
   REAL(EB), POINTER :: VEL_ERR_NEW            !< Velocity mismatch at mesh or solid boundary (m/s)
   REAL(EB), POINTER :: X                      !< \f$ x \f$ coordinate of boundary cell center
   REAL(EB), POINTER :: Y                      !< \f$ y \f$ coordinate of boundary cell center
   REAL(EB), POINTER :: Z                      !< \f$ z \f$ coordinate of boundary cell center

   INTEGER, POINTER :: BACK_INDEX              !< WALL index of back side of obstruction or exterior wall cell
   INTEGER, POINTER :: BACK_MESH               !< Mesh number on back side of obstruction or exterior wall cell
   INTEGER, POINTER :: BOUNDARY_TYPE           !< Descriptor: SOLID, MIRROR, OPEN, INTERPOLATED, etc
   INTEGER, POINTER :: OBST_INDEX              !< Index of the OBSTruction
   INTEGER, POINTER :: PRESSURE_BC_INDEX       !< Poisson boundary condition, NEUMANN or DIRICHLET
   INTEGER, POINTER :: SURF_INDEX              !< Index of the SURFace conditions
   INTEGER, POINTER :: SURF_INDEX_ORIG         !< Original SURFace index for this cell
   INTEGER, POINTER :: VENT_INDEX              !< Index of the VENT containing this cell
   INTEGER, POINTER :: WALL_INDEX              !< Self-identifier
   INTEGER, POINTER :: LAPLACE_BC_INDEX
   INTEGER, POINTER :: JD11_INDEX
   INTEGER, POINTER :: JD12_INDEX
   INTEGER, POINTER :: JD21_INDEX
   INTEGER, POINTER :: JD22_INDEX

END TYPE WALL_TYPE


!> \brief Variables associated with the external boundary of a mesh

TYPE EXTERNAL_WALL_TYPE
   INTEGER :: NOM                                     !< Number of the adjacent (Other) Mesh
   INTEGER :: NIC_MIN                                 !< Start of indices for the cell in the other mesh
   INTEGER :: NIC_MAX                                 !< End of indices for the cell in the other mesh
   INTEGER :: IIO_MIN                                 !< Minimum I index of adjacent cell in other mesh
   INTEGER :: IIO_MAX                                 !< Maximum I index of adjacent cell in other mesh
   INTEGER :: JJO_MIN                                 !< Minimum J index of adjacent cell in other mesh
   INTEGER :: JJO_MAX                                 !< Maximum J index of adjacent cell in other mesh
   INTEGER :: KKO_MIN                                 !< Minimum K index of adjacent cell in other mesh
   INTEGER :: KKO_MAX                                 !< Maximum K index of adjacent cell in other mesh
   REAL(EB) :: AREA_RATIO                             !< Ratio of face areas of adjoining cells
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: FVN         !< Flux-limited \f$ \int \rho Y_\alpha u_n \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: FVNS        !< Estimated value of FVN at next time step
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_D_DZDN  !< Species diffusive flux as computed in divg.f90
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_D_DZDNS !< RHO_D_DZDN estimated at next time step
END TYPE EXTERNAL_WALL_TYPE


!> \brief Derived type used to hold back side wall properties

TYPE EXPOSED_WALL_TYPE
   REAL(EB) :: Q_RAD_IN  !< Incoming radiation heat flux (W/m2)
   REAL(EB) :: TMP_GAS   !< Gas temperature (K)
END TYPE EXPOSED_WALL_TYPE


!> \brief Variables associated with a single primitive gas species

TYPE SPECIES_TYPE

   REAL(EB) :: MW=0._EB                           !< Molecular weight (g/mol)
   REAL(EB) :: YY0=0._EB                          !< Inital mass fraction
   REAL(EB) :: RCON                               !< Gas constant divided by molecular weight, \f$ R_0/W \f$ (J/kg/K)
   REAL(EB) :: MASS_EXTINCTION_COEFFICIENT=0._EB  !< Light extinction coefficient (m2/kg)
   REAL(EB) :: SPECIFIC_HEAT=-1._EB               !< Specific heat input by user (J/kg/K)
   REAL(EB) :: REFERENCE_ENTHALPY=-1._EB          !< Enthalpy at reference temperature (J/kg)
   REAL(EB) :: REFERENCE_TEMPERATURE              !< Basis temperature (K) for sensible enthalpy calculation
   REAL(EB) :: MU_USER=-1._EB                     !< User-specified viscosity (kg/m/s)
   REAL(EB) :: K_USER=-1._EB                      !< User-specified thermal conductivity (W/m/K)
   REAL(EB) :: D_USER=-1._EB                      !< User-specified diffusivity (m2/s)
   REAL(EB) :: EPSK=-1._EB                        !< Lennard-Jones \f$ \epsilon/k \f$ (K)
   REAL(EB) :: SIG=-1._EB                         !< Lennard_Jones hard-sphere diameter (Angstroms)
   REAL(EB) :: PR_USER=-1._EB                     !< User-specified Prandtl number
   REAL(EB) :: FLD_LETHAL_DOSE=0._EB
   REAL(EB) :: FIC_CONCENTRATION=0._EB
   REAL(EB) :: SPECIFIC_HEAT_LIQUID=-1            !< Liquid specific heat (J/kg/K)
   REAL(EB) :: DENSITY_LIQUID                     !< Liquid density (kg/m3)
   REAL(EB) :: HEAT_OF_VAPORIZATION=-1._EB        !< Heat of vaporization (J/kg)
   REAL(EB) :: H_F                                !< Heat of fusion (J/kg)
   REAL(EB) :: H_V_REFERENCE_TEMPERATURE=-1._EB   !< Heat of vaporization reference temperature (K)
   REAL(EB) :: TMP_V=-1._EB                       !< Vaporization temperature (K)
   REAL(EB) :: TMP_MELT=-1._EB                    !< Melting temperature (K)
   REAL(EB) :: ATOMS(118)=0._EB                   !< Atom count for molecular formula
   REAL(EB) :: MEAN_DIAMETER=1.E-6_EB             !< Diameter for aerosol (m)
   REAL(EB) :: CONDUCTIVITY_SOLID                 !< Thermal conductivity of solid (W/m/K)
   REAL(EB) :: DENSITY_SOLID                      !< Densith of solid (kg/m3)
   REAL(EB) :: BETA_LIQUID                        !< Coefficient of thermal expansion of the liquid (1/K)
   REAL(EB) :: MU_LIQUID                          !< Viscosity of the liquid (kg/m/s)
   REAL(EB) :: K_LIQUID                           !< Conductivity of the liquid (W/m/K)
   REAL(EB) :: PR_LIQUID                          !< Prandtl number of the liquid
   REAL(EB) :: THERMOPHORETIC_DIAMETER=0.03E-6_EB !< For use in aerosol deposition (m)

   LOGICAL ::  ISFUEL=.FALSE.                     !< Fuel species
   LOGICAL ::  LISTED=.FALSE.                     !< Properties are known to FDS
   LOGICAL ::  AGGLOMERATING=.FALSE.              !< Can form a particle or aerosol
   LOGICAL ::  EXPLICIT_H_F=.FALSE.               !< Heat of Formation is explicitly specified
   LOGICAL ::  CONDENSABLE=.FALSE.                !< Species can condense to liquid form

   CHARACTER(LABEL_LENGTH) :: ID                  !< Species name
   CHARACTER(LABEL_LENGTH) :: RAMP_CP             !< Name of specific heat ramp
   CHARACTER(LABEL_LENGTH) :: RAMP_CP_L           !< Name of liquid specific heat rame
   CHARACTER(LABEL_LENGTH) :: RAMP_K              !< Name of conductivity ramp
   CHARACTER(LABEL_LENGTH) :: RAMP_MU             !< Name of viscosity rame
   CHARACTER(LABEL_LENGTH) :: RAMP_D              !< Name of diffusivity ramp
   CHARACTER(LABEL_LENGTH) :: RADCAL_ID           !< Name of closest species with RADCAL properties
   CHARACTER(LABEL_LENGTH) :: RAMP_G_F
   CHARACTER(LABEL_LENGTH) :: PROP_ID             !< Name of PROPerty parameters
   CHARACTER(FORMULA_LENGTH) :: FORMULA           !< Chemical formula
   INTEGER :: MODE=2
   INTEGER :: RAMP_CP_INDEX=-1                    !< Index of specific heat ramp
   INTEGER :: RAMP_CP_L_INDEX=-1                  !< Index of liquid specific heat ramp
   INTEGER :: RAMP_K_INDEX=-1                     !< Index of conductivity ramp
   INTEGER :: RAMP_MU_INDEX=-1                    !< Index of viscosity ramp
   INTEGER :: RAMP_D_INDEX=-1                     !< Index of diffusivity heat ramp
   INTEGER :: RADCAL_INDEX=-1                     !< Index of nearest species with RADCAL properties
   INTEGER :: RAMP_G_F_INDEX=-1
   INTEGER :: AWM_INDEX=-1
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: H_V       !< Heat of vaporization as a function of temperature
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: C_P_L     !< Liquid specific heat as a function of temperature
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: C_P_L_BAR !< Average liquid specific heat as a function of temperture
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: H_L       !< Heat of liquification as a function of temperature

END TYPE SPECIES_TYPE

TYPE (SPECIES_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: SPECIES


!> \brief Properties of lumped species

TYPE SPECIES_MIXTURE_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MASS_FRACTION    !< Mass fractions of components
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: VOLUME_FRACTION  !< Volume fractions of components

   REAL(EB) :: MW                                  !< Molecular weight (g/mol)
   REAL(EB) :: RCON                                !< Specific gas constant, \f$ R_0 \sum_\alpha Y_\alpha/W_\alpha \f$ (J/kg/K)
   REAL(EB) :: ZZ0=0._EB                           !< Initial mass fraction of lumped species
   REAL(EB) :: MASS_EXTINCTION_COEFFICIENT=0._EB   !< Absorption coefficient of visible light (m2/kg)
   REAL(EB) :: ADJUST_NU=1._EB                     !< Adjustment factor if stoichiometric coefficients given for non-normalized VF
   REAL(EB) :: ATOMS(118)=0._EB                    !< Count of each atom in the mixture
   REAL(EB) :: MEAN_DIAMETER
   REAL(EB) :: SPECIFIC_HEAT=-1._EB                !< Specific heat (J/kg/K)
   REAL(EB) :: REFERENCE_ENTHALPY=-2.E20_EB        !< Enthalpy at REFERENCE_TEMPERATURE (J/kg)
   REAL(EB) :: THERMOPHORETIC_DIAMETER
   REAL(EB) :: REFERENCE_TEMPERATURE               !< Reference temperature of mixture (K)
   REAL(EB) :: MU_USER=-1._EB                      !< User-specified viscosity (kg/m/s)
   REAL(EB) :: K_USER=-1._EB                       !< User-specified thermal conductivity (W/m/K)
   REAL(EB) :: D_USER=-1._EB                       !< User-specified diffusion coefficient (m2/s)
   REAL(EB) :: PR_USER=-1._EB                      !< User-specified Prandhl number
   REAL(EB) :: EPSK=-1._EB
   REAL(EB) :: SIG=-1._EB
   REAL(EB) :: FLD_LETHAL_DOSE=0._EB
   REAL(EB) :: FIC_CONCENTRATION=0._EB
   REAL(EB) :: DENSITY_SOLID
   REAL(EB) :: CONDUCTIVITY_SOLID
   REAL(EB) :: H_F=-1.E30_EB

   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: SPEC_ID  !< Array of component species names
   CHARACTER(LABEL_LENGTH) :: ID='null'                           !< Name of lumped species
   CHARACTER(LABEL_LENGTH) :: RAMP_CP                             !< Name of specific heat ramp
   CHARACTER(LABEL_LENGTH) :: RAMP_CP_L                           !< Name of liquid specific heat ramp
   CHARACTER(LABEL_LENGTH) :: RAMP_K                              !< Name of conductivity ramp
   CHARACTER(LABEL_LENGTH) :: RAMP_MU                             !< Name of viscosity ramp
   CHARACTER(LABEL_LENGTH) :: RAMP_D                              !< Name of diffusion coefficient ramp
   CHARACTER(LABEL_LENGTH) :: RAMP_G_F
   CHARACTER(FORMULA_LENGTH) :: FORMULA='null'                    !< Chemical formula of lumped species

   INTEGER :: AWM_INDEX = -1,RAMP_CP_INDEX=-1,SINGLE_SPEC_INDEX=-1,RAMP_K_INDEX=-1,RAMP_MU_INDEX=-1,RAMP_D_INDEX=-1,&
              RAMP_G_F_INDEX=-1,CONDENSATION_SMIX_INDEX=-1,EVAPORATION_SMIX_INDEX=-1,AGGLOMERATION_INDEX=-1
   LOGICAL :: DEPOSITING=.FALSE.,VALID_ATOMS=.TRUE.,EVAPORATING=.FALSE.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: WQABS,WQSCA
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: R50

END TYPE SPECIES_MIXTURE_TYPE

TYPE (SPECIES_MIXTURE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: SPECIES_MIXTURE

TYPE REACTION_TYPE
   CHARACTER(LABEL_LENGTH) :: FUEL,OXIDIZER,PRODUCTS,ID,RAMP_CHI_R
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: SPEC_ID_NU,SPEC_ID_NU_READ,SPEC_ID_N_S,SPEC_ID_N_S_READ
   REAL(EB) :: C,H,N,O,EPUMO2,HEAT_OF_COMBUSTION,HOC_COMPLETE,A,A_PRIME,A_PRIME_FAST,A_IN,E,E_IN,K,&
               Y_O2_INFTY,Y_N2_INFTY=0._EB,MW_FUEL,MW_SOOT,Y_O2_MIN,&
               CO_YIELD,SOOT_YIELD,H2_YIELD,HCN_YIELD,SOOT_H_FRACTION,RHO_EXPONENT,RHO_EXPONENT_FAST,CRIT_FLAME_TMP,&
               NU_O2=0._EB,NU_N2=0._EB,NU_H2O=0._EB,NU_H2=0._EB,NU_HCN=0._EB,NU_CO2=0._EB,NU_CO=0._EB,NU_SOOT=0._EB,S=0._EB,&
               COMPLETE_HEAT_OF_COMBUSTION
   REAL(EB) :: N_T=0._EB,CHI_R
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: NU,NU_READ,NU_SPECIES,N_S,N_S_READ,NU_MW_O_MW_F
   INTEGER :: FUEL_SMIX_INDEX=-1,AIR_SMIX_INDEX=-1,N_SMIX,N_SPEC,RAMP_CHI_R_INDEX=0,PRIORITY=1
   LOGICAL :: IDEAL,CHECK_ATOM_BALANCE,FAST_CHEMISTRY=.FALSE.,REVERSE=.FALSE.,THIRD_BODY=.FALSE.
   CHARACTER(MESSAGE_LENGTH) :: FYI='null'
   CHARACTER(255) :: EQUATION
   CHARACTER(100) :: FWD_ID
END TYPE REACTION_TYPE

TYPE (REACTION_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: REACTION

TYPE MATERIAL_TYPE
   REAL(EB) :: K_S,C_S,RHO_S,EMISSIVITY,THERMAL_DIFFUSIVITY,KAPPA_S,TMP_BOIL,REFRACTIVE_INDEX,H(0:5000)=0._EB
   INTEGER :: PYROLYSIS_MODEL
   CHARACTER(LABEL_LENGTH) :: ID
   CHARACTER(LABEL_LENGTH) :: RAMP_H_R(MAX_REACTIONS),RAMP_K_S,RAMP_C_S
   INTEGER :: N_REACTIONS,PROP_INDEX=-1,H_R_I(MAX_REACTIONS)=0
   INTEGER, DIMENSION(MAX_REACTIONS) :: N_RESIDUE
   INTEGER, DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: RESIDUE_MATL_INDEX
   INTEGER, DIMENSION(3) :: RGB
   REAL(EB), DIMENSION(MAX_REACTIONS) :: TMP_REF,TMP_THR,RATE_REF,THR_SIGN
   REAL(EB), DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: NU_RESIDUE=0._EB
   REAL(EB), DIMENSION(MAX_REACTIONS) :: A,E,H_R,N_S,N_T,N_O2,GAS_DIFFUSION_DEPTH,NU_O2_CHAR,BETA_CHAR
   REAL(EB), DIMENSION(MAX_REACTIONS) :: HEATING_RATE,PYROLYSIS_RANGE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: NU_GAS,ADJUST_BURN_RATE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DIFFUSIVITY_GAS
   REAL(EB), DIMENSION(MAX_SPECIES,MAX_REACTIONS) :: NU_SPEC,HEAT_OF_COMBUSTION
   REAL(EB), DIMENSION(MAX_SPECIES) :: DIFFUSIVITY_SPEC
   LOGICAL, DIMENSION(MAX_REACTIONS) :: PCR = .FALSE.
   LOGICAL :: ALLOW_SHRINKING, ALLOW_SWELLING
   CHARACTER(LABEL_LENGTH), DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: RESIDUE_MATL_NAME
   CHARACTER(LABEL_LENGTH), DIMENSION(MAX_SPECIES,MAX_REACTIONS) :: SPEC_ID
   CHARACTER(MESSAGE_LENGTH) :: FYI='null'
END TYPE MATERIAL_TYPE

TYPE (MATERIAL_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: MATERIAL

!> \brief Variables associated with a surface type

TYPE SURFACE_TYPE
   REAL(EB) :: AREA_MULTIPLIER=1._EB                     !< Factor for manual surface area adjustment
   REAL(EB) :: TMP_FRONT=-1._EB,TMP_BACK=-1._EB,VEL,VEL_GRAD,PLE, &
               Z0,Z_0,CONVECTIVE_HEAT_FLUX,NET_HEAT_FLUX, &
               VOLUME_FLOW,HRRPUA,MLRPUA,T_IGN,SURFACE_DENSITY,CELL_SIZE_FACTOR, &
               E_COEFFICIENT,TEXTURE_WIDTH,TEXTURE_HEIGHT,THICKNESS,EXTERNAL_FLUX, &
               DXF,DXB,MASS_FLUX_TOTAL,PARTICLE_MASS_FLUX,EMISSIVITY,MAX_PRESSURE, &
               TMP_IGN,TMP_EXT,H_V,LAYER_DIVIDE,ROUGHNESS,LENGTH=-1._EB,WIDTH=-1._EB, &
               DT_INSERT,H_FIXED=-1._EB,H_FIXED_B=-1._EB,HM_FIXED=-1._EB,EMISSIVITY_BACK, &
               CONV_LENGTH,XYZ(3),FIRE_SPREAD_RATE, &
               MINIMUM_LAYER_THICKNESS,INNER_RADIUS=0._EB,MASS_FLUX_VAR=-1._EB,VEL_BULK, &
               PARTICLE_SURFACE_DENSITY=-1._EB,DRAG_COEFFICIENT=2.8_EB,SHAPE_FACTOR=0.25_EB,&
               MINIMUM_BURNOUT_TIME=1.E6_EB,DELTA_TMP_MAX=10._EB,BURN_DURATION=1.E6_EB,CONE_HEAT_FLUX=-1._EB
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DX,RDX,RDXN,X_S,DX_WGT,MF_FRAC,PARTICLE_INSERT_CLOCK
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: RHO_0
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MASS_FRACTION,MASS_FLUX,TAU,ADJUST_BURN_RATE
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: RAMP_INDEX
   INTEGER, DIMENSION(3) :: RGB
   REAL(EB) :: TRANSPARENCY
   REAL(EB), DIMENSION(2) :: VEL_T
   INTEGER, DIMENSION(2) :: LEAK_PATH,DUCT_PATH
   INTEGER :: THERMAL_BC_INDEX,NPPC,SPECIES_BC_INDEX,VELOCITY_BC_INDEX,SURF_TYPE,N_CELLS_INI,N_CELLS_MAX=0, &
              PART_INDEX,PROP_INDEX=-1,RAMP_T_I_INDEX=-1, RAMP_T_B_INDEX=0
   INTEGER :: PYROLYSIS_MODEL,NRA,NSB
   INTEGER :: N_LAYERS,N_MATL,SUBSTEP_POWER=2,N_SPEC=0
   INTEGER :: N_ONE_D_STORAGE_REALS,N_ONE_D_STORAGE_INTEGERS,N_ONE_D_STORAGE_LOGICALS
   INTEGER :: N_WALL_STORAGE_REALS,N_WALL_STORAGE_INTEGERS,N_WALL_STORAGE_LOGICALS
   INTEGER :: N_CFACE_STORAGE_REALS,N_CFACE_STORAGE_INTEGERS,N_CFACE_STORAGE_LOGICALS
   INTEGER, DIMENSION(30) :: ONE_D_REALS_ARRAY_SIZE=0,ONE_D_INTEGERS_ARRAY_SIZE=0,ONE_D_LOGICALS_ARRAY_SIZE=0
   INTEGER, ALLOCATABLE, DIMENSION(:) :: N_LAYER_CELLS,LAYER_INDEX,MATL_INDEX
   INTEGER, DIMENSION(MAX_LAYERS,MAX_MATERIALS) :: LAYER_MATL_INDEX
   INTEGER, DIMENSION(MAX_LAYERS) :: N_LAYER_MATL,N_LAYER_CELLS_MAX
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MIN_DIFFUSIVITY
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LAYER_THICKNESS,INTERNAL_HEAT_SOURCE
   REAL(EB), DIMENSION(MAX_LAYERS) :: LAYER_DENSITY,TMP_INNER,STRETCH_FACTOR,&
                                      MOISTURE_FRACTION,SURFACE_VOLUME_RATIO,PACKING_RATIO,KAPPA_S=-1._EB
   REAL(EB), DIMENSION(MAX_NUMBER_FSK_POINTS) :: FSK_K, FSK_W, FSK_A
   INTEGER :: NUMBER_FSK_POINTS
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: MATL_NAME
   CHARACTER(LABEL_LENGTH), DIMENSION(MAX_LAYERS,MAX_MATERIALS) :: LAYER_MATL_NAME
   REAL(EB), DIMENSION(MAX_LAYERS,MAX_MATERIALS) :: LAYER_MATL_FRAC
   LOGICAL :: BURN_AWAY,ADIABATIC,INTERNAL_RADIATION,USER_DEFINED=.TRUE., &
              FREE_SLIP=.FALSE.,NO_SLIP=.FALSE.,SPECIFIED_NORMAL_VELOCITY=.FALSE.,SPECIFIED_TANGENTIAL_VELOCITY=.FALSE., &
              SPECIFIED_NORMAL_GRADIENT=.FALSE.,CONVERT_VOLUME_TO_MASS=.FALSE.,SPECIFIED_HEAT_SOURCE=.FALSE.,&
              IMPERMEABLE=.FALSE.,BOUNDARY_FUEL_MODEL=.FALSE.,BLOWING=.FALSE.,BLOWING_2=.FALSE.,ABL_MODEL=.FALSE., &
              HT3D=.FALSE., MT1D=.FALSE.
   INTEGER :: GEOMETRY,BACKING,PROFILE,HEAT_TRANSFER_MODEL=0
   CHARACTER(LABEL_LENGTH) :: PART_ID,RAMP_Q,RAMP_V,RAMP_T,RAMP_EF,RAMP_PART,RAMP_V_X,RAMP_V_Y,RAMP_V_Z,RAMP_T_B,RAMP_T_I
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: RAMP_MF
   CHARACTER(LABEL_LENGTH) :: ID,TEXTURE_MAP,LEAK_PATH_ID(2)
   CHARACTER(MESSAGE_LENGTH) :: FYI='null'

   ! Level Set Firespread

   LOGICAL :: VEG_LSET_SPREAD,VEG_LSET_TAN2
   REAL(EB) :: VEG_LSET_IGNITE_T,VEG_LSET_ROS_HEAD,VEG_LSET_ROS_00,VEG_LSET_QCON,VEG_LSET_ROS_FLANK,VEG_LSET_ROS_BACK, &
               VEG_LSET_WIND_EXP,VEG_LSET_SIGMA,VEG_LSET_HT,VEG_LSET_BETA,&
               VEG_LSET_M1,VEG_LSET_M10,VEG_LSET_M100,VEG_LSET_MLW,VEG_LSET_MLH,VEG_LSET_SURF_LOAD,VEG_LSET_FIREBASE_TIME, &
               VEG_LSET_CHAR_FRACTION
   INTEGER :: VEG_LSET_FUEL_INDEX

   !HTC Custom
   REAL(EB) :: C_FORCED_CONSTANT=0._EB,C_FORCED_PR_EXP=0._EB,C_FORCED_RE=0._EB,C_FORCED_RE_EXP=0._EB,C_HORIZONTAL,C_VERTICAL

END TYPE SURFACE_TYPE

TYPE (SURFACE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: SURFACE

TYPE OMESH_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU,RHO,RHOS,TMP,U,V,W,US,VS,WS,H,HS,FVX,FVY,FVZ,D,DS,KRES,IL_S,IL_R,Q
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: ZZ,ZZS
   INTEGER, ALLOCATABLE, DIMENSION(:) :: IIO_R,JJO_R,KKO_R,IOR_R,IIO_S,JJO_S,KKO_S,IOR_S
   INTEGER, ALLOCATABLE, DIMENSION(:) :: N_PART_ORPHANS,N_PART_ADOPT,EXPOSED_WALL_CELL_BACK_INDICES,WALL_CELL_INDICES_SEND
   INTEGER :: I_MIN_R=-10,I_MAX_R=-10,J_MIN_R=-10,J_MAX_R=-10,K_MIN_R=-10,K_MAX_R=-10,NIC_R=0,N_WALL_CELLS_SEND=0, &
              I_MIN_S=-10,I_MAX_S=-10,J_MIN_S=-10,J_MAX_S=-10,K_MIN_S=-10,K_MAX_S=-10,NIC_S=0,N_EXPOSED_WALL_CELLS=0
   INTEGER, DIMENSION(7) :: INTEGER_SEND_BUFFER,INTEGER_RECV_BUFFER
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: &
         REAL_SEND_PKG1,REAL_SEND_PKG2,REAL_SEND_PKG3,REAL_SEND_PKG4,REAL_SEND_PKG5,REAL_SEND_PKG6,REAL_SEND_PKG7,REAL_SEND_PKG8,&
         REAL_RECV_PKG1,REAL_RECV_PKG2,REAL_RECV_PKG3,REAL_RECV_PKG4,REAL_RECV_PKG5,REAL_RECV_PKG6,REAL_RECV_PKG7,REAL_RECV_PKG8
   INTEGER :: N_EXTERNAL_OBST=0,N_INTERNAL_OBST=0
   TYPE (STORAGE_TYPE), ALLOCATABLE, DIMENSION(:) :: ORPHAN_PARTICLE_STORAGE,ADOPT_PARTICLE_STORAGE
   TYPE (EXPOSED_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: EXPOSED_WALL

   ! CC_IBM data exchange arrays:
   INTEGER :: NICC_S(2)=0, NICC_R(2)=0, NFCC_S(2)=0, NFCC_R(2)=0, NCC_INT_R=0, NFEP_R(5)=0, NFEP_R_G=0
   REAL(EB), ALLOCATABLE, DIMENSION(:) ::                &
         REAL_SEND_PKG11,REAL_SEND_PKG12,REAL_SEND_PKG13,&
         REAL_RECV_PKG11,REAL_RECV_PKG12,REAL_RECV_PKG13
   INTEGER, ALLOCATABLE, DIMENSION(:) :: ICC_UNKZ_CT_S, ICC_UNKZ_CC_S, ICC_UNKZ_CT_R, ICC_UNKZ_CC_R
   INTEGER, ALLOCATABLE, DIMENSION(:) :: UNKZ_CT_S, UNKZ_CC_S, UNKZ_CT_R, UNKZ_CC_R

   ! Face variables data (velocities):
   INTEGER, ALLOCATABLE, DIMENSION(:) :: IIO_FC_R,JJO_FC_R,KKO_FC_R,AXS_FC_R,IIO_FC_S,JJO_FC_S,KKO_FC_S,AXS_FC_S
   INTEGER, ALLOCATABLE, DIMENSION(:) :: IIO_CC_R,JJO_CC_R,KKO_CC_R,IIO_CC_S,JJO_CC_S,KKO_CC_S
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IFEP_R_1, IFEP_R_2, IFEP_R_3, IFEP_R_4, IFEP_R_5

   ! Level Set
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: PHI_LS,PHI1_LS,U_LS,V_LS,Z_LS
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: REAL_SEND_PKG14,REAL_RECV_PKG14

END TYPE OMESH_TYPE


!> \brief Variables associated with a rectangular OBSTruction

TYPE OBSTRUCTION_TYPE

   CHARACTER(LABEL_LENGTH) :: DEVC_ID='null'  !< Name of controlling device
   CHARACTER(LABEL_LENGTH) :: CTRL_ID='null'  !< Name of controller
   CHARACTER(LABEL_LENGTH) :: PROP_ID='null'  !< Name of PROPerty type
   CHARACTER(LABEL_LENGTH) :: MATL_ID='null'  !< Name of material type
   CHARACTER(LABEL_LENGTH) :: ID='null'       !< Name of obstruction

   INTEGER, DIMENSION(-3:3) :: SURF_INDEX=0   !< SURFace properties for each face
   INTEGER, DIMENSION(3) :: RGB=(/0,0,0/)     !< Color indices for Smokeview

   REAL(EB) :: TRANSPARENCY=1._EB             !< Transparency index for Smokeview, 0=invisible, 1=solid
   REAL(EB) :: VOLUME_ADJUST=1._EB            !< Effective volume divided by user specified volume
   REAL(EB) :: BULK_DENSITY=-1._EB            !< Mass per unit volume (kg/m3) of specified OBST
   REAL(EB) :: INTERNAL_HEAT_SOURCE=0._EB     !< Energy generation rate per unit volume (W/m3)
   REAL(EB) :: X1=0._EB                       !< Lower specified \f$ x \f$ boundary (m)
   REAL(EB) :: X2=1._EB                       !< Upper specified \f$ x \f$ boundary (m)
   REAL(EB) :: Y1=0._EB                       !< Lower specified \f$ y \f$ boundary (m)
   REAL(EB) :: Y2=1._EB                       !< Upper specified \f$ y \f$ boundary (m)
   REAL(EB) :: Z1=0._EB                       !< Lower specified \f$ z \f$ boundary (m)
   REAL(EB) :: Z2=1._EB                       !< Upper specified \f$ z \f$ boundary (m)
   REAL(EB) :: MASS=1.E6_EB                   !< Actual mass of the obstruction (kg)

   REAL(EB), DIMENSION(3) :: INPUT_AREA=-1._EB           !< Specified area of x, y, and z faces (m2)
   REAL(EB), DIMENSION(3) :: UNDIVIDED_INPUT_AREA=-1._EB !< Area of x, y, z faces (m2) unbroken by mesh boundaries
   REAL(EB), DIMENSION(3) :: SHAPE_AREA=0._EB            !< Area of idealized top, sides, bottom (m2)
   REAL(EB), DIMENSION(3) :: TEXTURE=0._EB               !< Origin of texture map (m)
   REAL(EB), DIMENSION(3) :: FDS_AREA=-1._EB             !< Effective areas of x, y, and z faces (m2)

   INTEGER :: I1=-1               !< Lower I node
   INTEGER :: I2=-1               !< Upper I node
   INTEGER :: J1=-1               !< Lower J node
   INTEGER :: J2=-1               !< Upper J node
   INTEGER :: K1=-1               !< Lower K node
   INTEGER :: K2=-1               !< Upper K node
   INTEGER :: COLOR_INDICATOR=-1  !< Coloring code: -3=use specified color, -2=invisible, -1=no color specified
   INTEGER :: TYPE_INDICATOR=-1   !< Smokeview code: 2=outline, -1=solid
   INTEGER :: ORDINAL=0           !< Order of OBST in input file
   INTEGER :: SHAPE_TYPE=-1       !< Indicator of shape carved out of larger obstruction
   INTEGER :: DEVC_INDEX=-1       !< Index of controlling device
   INTEGER :: CTRL_INDEX=-1       !< Index of controlling controller
   INTEGER :: PROP_INDEX=-1       !< Index of PROPerty type
   INTEGER :: DEVC_INDEX_O=-1     !< Original DEVC_INDEX
   INTEGER :: CTRL_INDEX_O=-1     !< Original CTRL_INDEX
   INTEGER :: MATL_INDEX=-1       !< Index of material
   INTEGER :: MULT_INDEX=-1       !< Index of multiplier function
   INTEGER :: RAMP_Q_INDEX=0      !< Index of HRR ramp

   LOGICAL, DIMENSION(-3:3) :: SHOW_BNDF=.TRUE. !< Show boundary quantities in Smokeview
   LOGICAL :: HIDDEN=.FALSE.                    !< Hide obstruction in Smokeview and ignore in simulation
   LOGICAL :: PERMIT_HOLE=.TRUE.                !< Allow the obstruction to have a hole cutout
   LOGICAL :: ALLOW_VENT=.TRUE.                 !< Allow a VENT to sit on the OBST
   LOGICAL :: CONSUMABLE=.FALSE.                !< The obstruction can burn away
   LOGICAL :: REMOVABLE=.FALSE.                 !< The obstruction can be removed from the simulation
   LOGICAL :: HOLE_FILLER=.FALSE.               !< The obstruction fills a HOLE
   LOGICAL :: OVERLAY=.TRUE.                    !< The obstruction can have another obstruction overlap a surface
   LOGICAL :: SCHEDULED_FOR_REMOVAL=.FALSE.     !< The obstruction is scheduled for removal during the current time step
   LOGICAL :: SCHEDULED_FOR_CREATION=.FALSE.    !< The obstruction is scheduled for creation during the current time step

   ! 3D pyrolysis:
   LOGICAL :: PYRO3D=.FALSE.
   LOGICAL :: MT3D=.FALSE.
   LOGICAL :: HT3D=.FALSE.
   LOGICAL :: PYRO3D_LIQUID=.FALSE.
   INTEGER :: MATL_SURF_INDEX=-1
   INTEGER :: PYRO3D_IOR=0
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: RHO

END TYPE OBSTRUCTION_TYPE


TYPE TRIBIN_TYPE
   REAL(EB):: X1_LOW, X1_HIGH
   INTEGER :: NTL
   INTEGER, ALLOCATABLE, DIMENSION(:) :: TRI_LIST
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: XYZV
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FACECUBE
END TYPE TRIBIN_TYPE

TYPE TBAXIS_TYPE
   INTEGER :: N_BINS
   REAL(EB):: DELBIN
   TYPE(TRIBIN_TYPE), ALLOCATABLE, DIMENSION(:) :: TRIBIN
END TYPE TBAXIS_TYPE

TYPE TRANSFORM_TYPE
   CHARACTER(LABEL_LENGTH) :: GEOM_ID, ID
   CHARACTER(80) :: FILE
   REAL(EB) :: Z_OFFSET
   INTEGER :: GEOM_INDEX=-1
END TYPE TRANSFORM_TYPE

TYPE GEOMETRY_TYPE
   CHARACTER(LABEL_LENGTH) :: ID,MATL_ID,DEVC_ID,PROP_ID,MOVE_ID
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: SURF_ID
   CHARACTER(LABEL_LENGTH) :: BNDC_FILENAME='null',GEOC_FILENAME='null',TEXTURE_MAPPING
   LOGICAL :: COMPONENT_ONLY,IS_DYNAMIC=.TRUE.,HAVE_SURF,HAVE_MATL,AUTO_TEXTURE,HIDDEN,REMOVEABLE,SHOW_BNDF=.TRUE., &
              READ_BINARY=.FALSE.,SNAP_TO_GRID=.FALSE.,IS_TERRAIN=.FALSE.
   INTEGER :: N_VERTS_BASE,N_FACES_BASE,N_VOLUS_BASE,N_VERTS,N_EDGES,N_FACES,N_VOLUS,NSUB_GEOMS,GEOM_TYPE,IJK(3),N_LEVELS,&
              DEVC_INDEX=-1,CTRL_INDEX=-1,PROP_INDEX=-1,DEVC_INDEX_O=-1,CTRL_INDEX_O=-1,MATL_INDEX=-1,&
              CYLINDER_NSEG_THETA,CYLINDER_NSEG_AXIS
   INTEGER, DIMENSION(3) :: RGB=-1
   REAL(EB) :: GEOM_VOLUME,GEOM_AREA,GEOM_XYZCEN(3),OMEGA=0._EB,XYZ0(3),AZIM_BASE,ELEV_BASE,SCALE_BASE(3),XYZ_BASE(3),&
               AZIM,ELEV,SCALE(3),XYZ(3),AZIM_DOT,ELEV_DOT,SCALE_DOT(3),XYZ_DOT(3),GAXIS(3),GROTATE,GROTATE_DOT,GROTATE_BASE,&
               XB(6),XB_ORIG(6),SPHERE_ORIGIN(3),SPHERE_RADIUS,TEXTURE_ORIGIN(3),TEXTURE_SCALE(2),MIN_LEDGE,MAX_LEDGE,MEAN_LEDGE,&
               GEOM_BOX(LOW_IND:HIGH_IND,IAXIS:KAXIS),TRANSPARENCY,CYLINDER_ORIGIN(3),CYLINDER_AXIS(3),&
               CYLINDER_RADIUS,CYLINDER_LENGTH
   INTEGER,ALLOCATABLE,DIMENSION(:) :: FACES,VOLUS,SUB_GEOMS,SURFS,MATLS
   INTEGER,ALLOCATABLE,DIMENSION(:,:) :: EDGES,FACE_EDGES,EDGE_FACES
   REAL(EB),ALLOCATABLE,DIMENSION(:) :: FACES_AREA,VERTS_BASE,VERTS,TFACES,DAZIM,DELEV,ZVALS
   REAL(EB),ALLOCATABLE,DIMENSION(:,:) :: FACES_NORMAL,DSCALE,DXYZ0,DXYZ
   REAL(EB),ALLOCATABLE, DIMENSION(:,:,:) :: FACECUBE
   TYPE(TBAXIS_TYPE) :: TBAXIS(IAXIS:KAXIS)
END TYPE GEOMETRY_TYPE

INTEGER :: N_GEOMETRY=0, GEOMETRY_CHANGE_STATE=0, N_TRNF=0
TYPE(GEOMETRY_TYPE),  ALLOCATABLE, TARGET, DIMENSION(:)   :: GEOMETRY
TYPE(GEOMETRY_TYPE),  ALLOCATABLE, TARGET, DIMENSION(:,:) :: GEOMETRY_TRANSFORM
TYPE(TRANSFORM_TYPE), ALLOCATABLE, TARGET, DIMENSION(:,:) :: TRANSFORM
LOGICAL :: IS_GEOMETRY_DYNAMIC

TYPE VERTEX_TYPE
   REAL(EB) :: X,Y,Z
   REAL(EB) :: XYZ(3) ! this form may be more efficient
END TYPE VERTEX_TYPE

TYPE(VERTEX_TYPE), TARGET, ALLOCATABLE, DIMENSION(:) :: VERTEX

! --- CC_IBM types:
! Edge crossings data structure:
INTEGER, PARAMETER :: IBM_MAXCROSS_EDGE = 10 ! Size definition parameter. Max number of crossings per Cartesian Edge.
TYPE IBM_EDGECROSS_TYPE
   INTEGER :: NCROSS   ! Number of BODINT_PLANE segments - Cartesian edge crossings.
   REAL(EB), DIMENSION(1:IBM_MAXCROSS_EDGE)   ::  SVAR ! Locations along x2 axis of NCROSS intersections.
   INTEGER,  DIMENSION(1:IBM_MAXCROSS_EDGE)   :: ISVAR ! Type of intersection (i.e. SG, GS or GG).
   INTEGER,  DIMENSION(MAX_DIM+1)             ::   IJK ! [ i j k X2AXIS]
END TYPE IBM_EDGECROSS_TYPE

! Cartesian Edge Cut-Edges data structure:
TYPE IBM_CUTEDGE_TYPE
   INTEGER :: NVERT, NEDGE, NEDGE1, STATUS         ! Local Vertices, cut-edges and status of this Cartesian edge.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           :: XYZVERT  ! Locations of vertices.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  CEELEM  ! Cut-Edge connectivities.
   INTEGER,  DIMENSION(MAX_DIM+2)                  ::     IJK  ! [ i j k X2AXIS cetype]
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  INDSEG  ! [ntr tr1 tr2 ibod]
   INTEGER,  ALLOCATABLE, DIMENSION(:)             ::NOD_PERM  ! Permutation array for INSERT_FACE_VERT.
END TYPE IBM_CUTEDGE_TYPE

! IBM_EDGE type, used for computation of wall model turbulent viscosity, shear stress, vorticity.
TYPE IBM_EDGE_TYPE
   INTEGER,  DIMENSION(MAX_DIM+1)                  ::     IJK  ! [ i j k X1AXIS]
   INTEGER :: IE=0
   ! Here: VIND=IAXIS:KAXIS, EP=1:INT_N_EXT_PTS,
   ! INT_VEL_IND = 1; INT_VELS_IND = 2; INT_FV_IND = 3; INT_DHDX_IND = 4; N_INT_FVARS = 4;
   ! INT_NPE_LO = INT_NPE(LOW,VIND,EP,IFACE); INT_NPE_LO = INT_NPE(HIGH,VIND,EP,IEDGE).
   ! IEDGE = 0, Cartesian GASPHASE EDGE.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_IJK        ! (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:)        :: INT_COEF       ! (INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_DCOEF      ! (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_XYZBF ,INT_NOUT     ! (IAXIS:KAXIS,IEDGE)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_INBFC      ! (1:3,IEDGE)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:,:,:)  :: INT_NPE        ! (LOW:HIGH,VIND,EP,IEDGE)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_XN,INT_CN  ! (0:INT_N_EXT_PTS,IEDGE)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_FVARS      ! (1:N_INT_FVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_CVARS      ! (1:N_INT_CVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_NOMIND     ! (LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:)        :: XB_IB
   INTEGER,  ALLOCATABLE, DIMENSION(:)        :: SURF_INDEX
   LOGICAL,  ALLOCATABLE, DIMENSION(:)        :: PROCESS_EDGE_ORIENTATION,EDGE_IN_MESH
END TYPE IBM_EDGE_TYPE

! Cartesian Faces Cut-Faces data structure:
INTEGER, PARAMETER :: IBM_MAXVERTS_FACE  =3072 ! Size definition parameter. Max number of vertices per Cartesian Face.
INTEGER, PARAMETER :: IBM_MAXCEELEM_FACE =3072 ! Size definition parameter. Max segments per face.
INTEGER, PARAMETER :: IBM_MAXCFELEM_FACE =3072 ! Size definition parameter. Max number of cut faces per Cartesian Face.
INTEGER, PARAMETER :: IBM_MAXVERT_CUTFACE=  24 ! Size definition parameter.
INTEGER, PARAMETER :: MAX_INTERP_POINTS_PLANE = 4
TYPE IBM_CUTFACE_TYPE
   INTEGER :: IWC=0
   INTEGER :: NVERT=0, NSVERT=0, NFACE=0, NSFACE=0, STATUS !Local Vertices, cut-faces and status of this Cartesian face.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           :: XYZVERT  ! Locations of vertices.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  CFELEM  ! Cut-faces connectivities.
   INTEGER,  DIMENSION(MAX_DIM+1)                  ::     IJK  ! [ i j k X1AXIS]
   REAL(EB), ALLOCATABLE, DIMENSION(:)             ::    AREA  ! Cut-faces areas.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  XYZCEN  ! Cut-faces centroid locations.
   !Integrals to be used in cut-cell volume and centroid computations.
   REAL(EB), ALLOCATABLE, DIMENSION(:)             ::  INXAREA, INXSQAREA, JNYSQAREA, KNZSQAREA
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  BODTRI
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  UNKH, UNKZ
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  XCENLOW, XCENHIGH
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  RHO
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  ZZ_FACE, DIFF_FACE, RHO_D, VELD
   REAL(EB), ALLOCATABLE, DIMENSION(:)             :: TMP_FACE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:)         :: RHO_D_DZDN
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           :: H_RHO_D_DZDN
   REAL(EB), ALLOCATABLE, DIMENSION(:)             ::  VEL, VELS, DHDX, FN, VELNP1, VELINT
   INTEGER,  ALLOCATABLE, DIMENSION(:,:,:)         ::  JDZ, JDH
   REAL(EB) :: VELN_CRF, VELD_CRF, DHDX_CRF, FN_CRF, VELNP1_CRF, VELINT_CRF
   INTEGER,  ALLOCATABLE, DIMENSION(:,:,:)                         ::      CELL_LIST ! [RC_TYPE I J K ]

   ! Here: VIND=IAXIS:KAXIS, EP=1:INT_N_EXT_PTS,
   ! INT_VEL_IND = 1; INT_VELS_IND = 2; INT_FV_IND = 3; INT_DHDX_IND = 4; N_INT_FVARS = 4;
   ! INT_NPE_LO = INT_NPE(LOW,VIND,EP,IFACE); INT_NPE_HI = INT_NPE(HIGH,VIND,EP,IFACE).
   ! IFACE = 0, undelying Cartesian face, IFACE=1:NFACE GASPHASE cut-faces.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_IJK        ! (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:)        :: INT_COEF       ! (INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_DCOEF      ! (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_XYZBF,INT_NOUT ! (IAXIS:KAXIS,0:NFACE)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_INBFC      ! (1:3,0:NFACE)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:,:,:)  :: INT_NPE        ! (LOW:HIGH,VIND,EP,0:NFACE)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_XN,INT_CN  ! (0:INT_N_EXT_PTS,0:NFACE)  ! 0 is interpolation point.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_FVARS      ! (1:N_INT_FVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_NOMIND     ! (LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

   ! Fields used in INBOUNDARY faces:
   ! Here: VIND=0, EP=1:INT_N_EXT_PTS,
   ! INT_H_IND=1, etc. N_INT_CVARS=INT_P_IND+N_TRACKED_SPECIES
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_CVARS      ! (1:N_INT_CVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

   REAL(EB), ALLOCATABLE, DIMENSION(:,:)            :: RHOPVN
   INTEGER :: NOMICF(2)=0  ! [NOM icf]
   INTEGER,  ALLOCATABLE, DIMENSION(:)              :: CFACE_INDEX, SURF_INDEX
END TYPE IBM_CUTFACE_TYPE


TYPE RAD_CFACE_TYPE
   INTEGER :: N_ASSIGNED_CFACES_RADI=0
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ASSIGNED_CFACES_RADI
END TYPE RAD_CFACE_TYPE

! Note: If you change the number of scalar variables in CFACE_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_CFACE_SCALAR_REALS=10,N_CFACE_SCALAR_INTEGERS=6,N_CFACE_SCALAR_LOGICALS=0

TYPE CFACE_TYPE
   TYPE (ONE_D_M_AND_E_XFER_TYPE) :: ONE_D
   REAL(EB), POINTER :: AREA,X,Y,Z,NVEC(:),VEL_ERR_NEW,V_DEP,Q_LEAK
   INTEGER, POINTER :: CFACE_INDEX,SURF_INDEX,VENT_INDEX,BOUNDARY_TYPE,CUT_FACE_IND1,CUT_FACE_IND2
END TYPE CFACE_TYPE

! Cartesian Cells Cut-Cells data structure:
INTEGER, PARAMETER :: IBM_MAXVERTS_CELL   =3072
INTEGER, PARAMETER :: IBM_NPARAM_CCFACE   =   5 ! [face_type side iaxis cei icf]

TYPE IBM_CUTCELL_TYPE
   INTEGER :: NCELL, NFACE_CELL
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)                     ::    CCELEM ! Cut-cells faces connectivities in FACE_LIST.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)                     :: FACE_LIST ! List of faces, cut-faces.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::    VOLUME ! Cut-cell volumes.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     ::    XYZCEN ! Cut-cell centroid locaitons.
   INTEGER,  DIMENSION(MAX_DIM)                              ::       IJK ! [ i j k ]
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: RHO, RHOS ! Cut cells densities.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::  RSUM,TMP ! Cut cells temperatures.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::    D,  DS, DVOL ! Cut cell thermodynamic divg.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::Q,QR,D_SOURCE ! Q,Thermo divg reaction component.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: CHI_R,MIX_TIME ! Cut-cell combustion
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     ::    Q_REAC          ! variables.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     :: REAC_SOURCE_TERM   !
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     ::   ZZ, ZZS, M_DOT_PPP ! Cut cells species mass
                                                                                ! fractions and rho*D_z,reaction source.
   INTEGER,  ALLOCATABLE, DIMENSION(:)                       :: UNKH,UNKZ ! Unknown number for pressure H,
                                                                          ! and scalars.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::      H,HS ! Pressure H containers.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::  RTRM,R_H_G,RHO_0,WVEL

   ! Here: VIND=0, EP=1:INT_N_EXT_PTS
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_IJK        ! (IAXIS:KAXIS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:)        :: INT_COEF       ! (INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_XYZBF,INT_NOUT ! (IAXIS:KAXIS,0:NCELL)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_INBFC      ! (1:3,0:NCELL)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:,:,:)  :: INT_NPE        ! (LOW:HIGH,0,EP,0:NCELL)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_XN,INT_CN  ! (0:INT_N_EXT_PTS,0:NCELL)  ! 0 is interpolation point.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: INT_CCVARS      ! (1:N_INT_CCVARS,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)      :: INT_NOMIND     ! (LOW_IND:HIGH_IND,INT_NPE_LO+1:INT_NPE_LO+INT_NPE_HI)

   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     :: DEL_RHO_D_DEL_Z_VOL, U_DOT_DEL_RHO_Z_VOL
   LOGICAL,  ALLOCATABLE, DIMENSION(:)                       :: USE_CC_VOL
   INTEGER :: NOMICC(2)=0
   REAL(EB):: DIVVOL_BC=0._EB
END TYPE IBM_CUTCELL_TYPE


TYPE IBM_REGFACE_TYPE
   INTEGER,  DIMENSION(MAX_DIM)                                    ::       IJK
   INTEGER,  DIMENSION(1:2,1:2)                                    ::        JD
END TYPE IBM_REGFACE_TYPE

TYPE IBM_REGFACEZ_TYPE
   INTEGER :: IWC=0
   INTEGER,  DIMENSION(MAX_DIM)                                    ::       IJK
   INTEGER,  DIMENSION(1:2,1:2)                                    ::        JD
   REAL(EB), DIMENSION(MAX_SPECIES)                                ::   DIFF_FACE=0._EB, RHO_D=0._EB, VELD=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES,LOW_IND:HIGH_IND)               ::   RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES)                                :: H_RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(-1:0)                                       ::    RHOPVN=0._EB
END TYPE IBM_REGFACEZ_TYPE

TYPE IBM_RCFACE_TYPE
   INTEGER,  DIMENSION(MAX_DIM+1)                                  ::       IJK ! [ I J K x1axis]
   INTEGER,  DIMENSION(LOW_IND:HIGH_IND)                           ::       UNK
   REAL(EB), DIMENSION(MAX_DIM,LOW_IND:HIGH_IND)                   ::      XCEN
   INTEGER,  DIMENSION(1:2,1:2)                                    ::        JD
END TYPE IBM_RCFACE_TYPE

TYPE IBM_RCFACE_LST_TYPE
   INTEGER :: IWC=0
   REAL(EB):: TMP_FACE=0._EB
   INTEGER,  DIMENSION(MAX_DIM+1)                                  ::       IJK ! [ I J K x1axis]
   INTEGER,  DIMENSION(LOW_IND:HIGH_IND)                           ::       UNK
   REAL(EB), DIMENSION(MAX_DIM,LOW_IND:HIGH_IND)                   ::      XCEN
   INTEGER,  DIMENSION(1:2,1:2)                                    ::        JD
   INTEGER,  DIMENSION(MAX_DIM+1,LOW_IND:HIGH_IND)                 :: CELL_LIST ! [RC_TYPE I J K ]
   REAL(EB), DIMENSION(MAX_SPECIES)                              :: ZZ_FACE=0._EB,DIFF_FACE=0._EB,RHO_D=0._EB,VELD=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES,LOW_IND:HIGH_IND)               :: RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES)                                :: H_RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(-1:0)                                       ::    RHOPVN=0._EB
END TYPE IBM_RCFACE_LST_TYPE

TYPE IBM_EXIMFACE_TYPE
   INTEGER :: LHFACE, UNKZ, IWC=0
   INTEGER,  DIMENSION(MAX_DIM+1)                                  ::       IJK ! [ I J K x1axis]
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                           ::       FLX
   REAL(EB) :: AREA,FN_H_S
   REAL(EB), DIMENSION(MAX_SPECIES)                                ::H_RHO_D_DZDN=0._EB,FN_ZZ=0._EB
END TYPE IBM_EXIMFACE_TYPE

TYPE CSVF_TYPE
    CHARACTER(255) :: CSVFILE,UVWFILE
END TYPE CSVF_TYPE

TYPE(CSVF_TYPE), ALLOCATABLE, DIMENSION(:) :: CSVFINFO

TYPE VENTS_TYPE
   INTEGER :: I1=-1,I2=-1,J1=-1,J2=-1,K1=-1,K2=-1,BOUNDARY_TYPE=0,IOR=0,SURF_INDEX=0,DEVC_INDEX=-1,CTRL_INDEX=-1, &
              COLOR_INDICATOR=99,TYPE_INDICATOR=0,ORDINAL=0,PRESSURE_RAMP_INDEX=0,TMP_EXTERIOR_RAMP_INDEX=0,NODE_INDEX=-1, &
              OBST_INDEX=0
   INTEGER, DIMENSION(3) :: RGB=-1
   REAL(EB) :: TRANSPARENCY = 1._EB
   REAL(EB), DIMENSION(3) :: TEXTURE=0._EB
   REAL(EB) :: X1=0._EB,X2=0._EB,Y1=0._EB,Y2=0._EB,Z1=0._EB,Z2=0._EB,FDS_AREA=0._EB, &
               X1_ORIG=0._EB,X2_ORIG=0._EB,Y1_ORIG=0._EB,Y2_ORIG=0._EB,Z1_ORIG=0._EB,Z2_ORIG=0._EB, &
               X0=-9.E6_EB,Y0=-9.E6_EB,Z0=-9.E6_EB,FIRE_SPREAD_RATE,UNDIVIDED_INPUT_AREA=0._EB,INPUT_AREA=0._EB,&
               TMP_EXTERIOR=-1000._EB,DYNAMIC_PRESSURE=0._EB,UVW(3)=-1.E12_EB,RADIUS=-1._EB
   LOGICAL :: ACTIVATED=.TRUE.,GHOST_CELLS_ONLY=.FALSE.,GEOM=.FALSE.
   CHARACTER(LABEL_LENGTH) :: DEVC_ID='null',CTRL_ID='null',ID='null'
   ! turbulent inflow (experimental)
   INTEGER :: N_EDDY=0
   REAL(EB) :: R_IJ(3,3)=0._EB,A_IJ(3,3)=0._EB,SIGMA_IJ(3,3)=0._EB,EDDY_BOX_VOLUME=0._EB, &
               X_EDDY_MIN=0._EB,X_EDDY_MAX=0._EB, &
               Y_EDDY_MIN=0._EB,Y_EDDY_MAX=0._EB, &
               Z_EDDY_MIN=0._EB,Z_EDDY_MAX=0._EB
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: U_EDDY,V_EDDY,W_EDDY
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: X_EDDY,Y_EDDY,Z_EDDY,CU_EDDY,CV_EDDY,CW_EDDY
END TYPE VENTS_TYPE

TYPE TABLES_TYPE
   INTEGER :: NUMBER_ROWS,NUMBER_COLUMNS
   REAL(EB) :: LX,LY,UX,UY
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: X,Y
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TABLE_DATA,Z
END TYPE TABLES_TYPE

TYPE (TABLES_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: TABLES

TYPE RAMPS_TYPE
   REAL(EB) :: SPAN,RDT,T_MIN,T_MAX
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: INDEPENDENT_DATA,DEPENDENT_DATA,INTERPOLATED_DATA
   INTEGER :: NUMBER_DATA_POINTS,NUMBER_INTERPOLATION_POINTS,DEVC_INDEX=-1,CTRL_INDEX=-1
   CHARACTER(LABEL_LENGTH) :: DEVC_ID='null',CTRL_ID='null'
   LOGICAL :: DEP_VAR_UNITS_CONVERTED=.FALSE.,RESERVED=.FALSE.
END TYPE RAMPS_TYPE

TYPE (RAMPS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: RAMPS

TYPE RESERVED_RAMPS_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: INDEPENDENT_DATA,DEPENDENT_DATA
   INTEGER :: NUMBER_DATA_POINTS
END TYPE RESERVED_RAMPS_TYPE

INTEGER :: N_RESERVED_RAMPS=0
TYPE (RESERVED_RAMPS_TYPE), DIMENSION(10), TARGET :: RESERVED_RAMPS


TYPE HUMAN_TYPE
   CHARACTER(LABEL_LENGTH) :: NODE_NAME='null'
   CHARACTER(LABEL_LENGTH) :: FFIELD_NAME='null'
   REAL(EB) :: X=0._EB,Y=0._EB,Z=0._EB,U=0._EB,V=0._EB,W=0._EB,F_X=0._EB,F_Y=0._EB,&
               X_old=0._EB,Y_old=0._EB,X_group=0._EB,Y_group=0._EB, U_CB=0.0_EB, V_CB=0.0_EB
   REAL(EB) :: UBAR=0._EB, VBAR=0._EB, UBAR_Center=0._EB, VBAR_Center=0._EB, T_LastRead_CB=0.0_EB
   REAL(EB) :: Speed=1.25_EB, Radius=0.255_EB, Mass=80.0_EB, Tpre=1._EB, Tau=1._EB, &
               Eta=0._EB, Ksi=0._EB, Tdet=0._EB, Speed_ave=0._EB, T_Speed_ave=0._EB, &
               X_CB=0.0_EB, Y_CB=0.0_EB, X2_CB=0.0_EB, Y2_CB=0.0_EB, UAVE_CB=0.0_EB, VAVE_CB=0.0_EB
   REAL(EB) :: r_torso=0.15_EB, r_shoulder=0.095_EB, d_shoulder=0.055_EB, angle=0._EB, &
               torque=0._EB, m_iner=4._EB
   REAL(EB) :: tau_iner=0.2_EB, angle_old=0._EB, omega=0._EB
   REAL(EB) :: A=2000._EB, B=0.08_EB, C_Young=120000._EB, Gamma=16000._EB, Kappa=40000._EB, &
               Lambda=0.5_EB, Commitment=0._EB
   REAL(EB) :: SumForces=0._EB, IntDose=0._EB, DoseCrit1=0._EB, DoseCrit2=0._EB, SumForces2=0._EB
   REAL(EB) :: TempMax1=0._EB, FluxMax1=0._EB, TempMax2=0._EB, FluxMax2=0._EB, Density=0._EB, DensityR=0._EB, DensityL=0._EB
   REAL(EB) :: P_detect_tot=0._EB, v0_fac=1._EB, D_Walls=0._EB
   REAL(EB) :: T_FallenDown=0._EB, F_FallDown=0._EB, Angle_FallenDown=0._EB, SizeFac_FallenDown=0._EB, T_CheckFallDown=0._EB
   INTEGER  :: IOR=-1, ILABEL=0, COLOR_INDEX=0, INODE=0, IMESH=-1, IPC=0, IEL=0, I_FFIELD=0, I_Target2=0, ID_CB=-1
   INTEGER  :: GROUP_ID=0, DETECT1=0, GROUP_SIZE=0, I_Target=0, I_DoorAlgo=0, I_Door_Mode=0, STRS_Direction = 1
   INTEGER  :: STR_SUB_INDX, SKIP_WALL_FORCE_IOR
   LOGICAL  :: SHOW=.TRUE., NewRnd=.TRUE., CROWBAR_READ_IN=.FALSE., CROWBAR_UPDATE_V0=.FALSE.
   LOGICAL  :: SeeDoorXB1=.FALSE., SeeDoorXB2=.FALSE., SeeDoorXYZ1=.FALSE., SeeDoorXYZ2=.FALSE.
END TYPE HUMAN_TYPE

TYPE HUMAN_GRID_TYPE
! (x,y,z) Centers of the grid cells in the main evacuation meshes
! SOOT_DENS: Smoke density at the center of the cell (mg/m3)
! FED_CO_CO2_O2: Purser's FED for co, co2, and o2
   REAL(EB) :: X,Y,Z,SOOT_DENS,FED_CO_CO2_O2,TMP_G,RADFLUX
   INTEGER :: N, N_old, IGRID, IHUMAN, ILABEL
! IMESH: (x,y,z) which fire mesh, if any
! II,JJ,KK: Fire mesh cell reference
   INTEGER  :: IMESH,II,JJ,KK
END TYPE HUMAN_GRID_TYPE

TYPE SLICE_TYPE
   INTEGER :: I1,I2,J1,J2,K1,K2,GEOM_INDEX=-1,TRNF_INDEX=-1,INDEX,INDEX2=0,Z_INDEX=-999,Y_INDEX=-999,MATL_INDEX=-999,&
              PART_INDEX=0,VELO_INDEX=0,PROP_INDEX=0,REAC_INDEX=0,SLCF_INDEX
   INTEGER, ALLOCATABLE, DIMENSION(:) :: REORDER_TO_KJI
   REAL(FB), DIMENSION(2) :: MINMAX
   REAL(FB) :: RLE_MIN, RLE_MAX
   REAL(EB):: AGL_SLICE
   LOGICAL :: TERRAIN_SLICE=.FALSE.,CELL_CENTERED=.FALSE.,FACE_CENTERED=.FALSE.,RLE=.FALSE.,MULTI_RES=.FALSE.
   CHARACTER(LABEL_LENGTH) :: SLICETYPE='STRUCTURED',SMOKEVIEW_LABEL
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_BAR_LABEL,ID='null',MATL_ID='null',TRNF_ID='null'
END TYPE SLICE_TYPE

TYPE RAD_FILE_TYPE
   INTEGER :: I1,I2,J1,J2,K1,K2,I_STEP=1,J_STEP=1,K_STEP=1,N_POINTS
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: IL_SAVE
END TYPE RAD_FILE_TYPE

TYPE PATCH_TYPE
   INTEGER :: I1,I2,J1,J2,K1,K2,IG1,IG2,JG1,JG2,KG1,KG2,IOR,OBST_INDEX
END TYPE PATCH_TYPE

TYPE BOUNDARY_FILE_TYPE
   INTEGER :: DEBUG=0,INDEX,PROP_INDEX,Z_INDEX=-999,Y_INDEX=-999,PART_INDEX=0,TIME_INTEGRAL_INDEX=0
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_LABEL
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_BAR_LABEL,UNITS,MATL_ID='null'
   LOGICAL :: CELL_CENTERED=.FALSE.
END TYPE BOUNDARY_FILE_TYPE

TYPE (BOUNDARY_FILE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: BOUNDARY_FILE

TYPE BOUNDARY_ELEMENT_FILE_TYPE
   INTEGER :: INDEX,PROP_INDEX,Z_INDEX=-999,Y_INDEX=-999,PART_INDEX=0
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_LABEL
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_BAR_LABEL
   LOGICAL :: CELL_CENTERED=.TRUE.
END TYPE BOUNDARY_ELEMENT_FILE_TYPE

TYPE (BOUNDARY_ELEMENT_FILE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: BOUNDARY_ELEMENT_FILE

TYPE ISOSURFACE_FILE_TYPE
   INTEGER :: INDEX=1,N_VALUES=1,Y_INDEX=-999,Z_INDEX=-999,VELO_INDEX=0,SKIP=1,&
              INDEX2=-1,Y_INDEX2=-999,Z_INDEX2=-999,VELO_INDEX2=0
   REAL(FB) :: VALUE(10), DELTA
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_LABEL, SMOKEVIEW_LABEL2
   CHARACTER(LABEL_LENGTH) :: SMOKEVIEW_BAR_LABEL, SMOKEVIEW_BAR_LABEL2
END TYPE ISOSURFACE_FILE_TYPE

TYPE (ISOSURFACE_FILE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: ISOSURFACE_FILE

TYPE PROFILE_TYPE
   REAL(EB) :: X,Y,Z
   INTEGER  :: IOR=0,WALL_INDEX=0,LP_TAG=0,ORDINAL,MESH,FORMAT_INDEX=1,PART_CLASS_INDEX=-1
   CHARACTER(LABEL_LENGTH) :: ID='null',QUANTITY='TEMPERATURE',INIT_ID='null'
END TYPE PROFILE_TYPE

TYPE (PROFILE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: PROFILE


!> \brief Parameters associated with initialization of particles or regions of the domain

TYPE INITIALIZATION_TYPE
   REAL(EB) :: TEMPERATURE      !< Temperature (K) of the initialized region
   REAL(EB) :: DENSITY          !< Density (kg/m3) of the initialized region
   REAL(EB) :: X1               !< Lower x boundary of the initialized region (m)
   REAL(EB) :: X2               !< Upper x boundary of the initialized region (m)
   REAL(EB) :: Y1               !< Lower y boundary of the initialized region (m)
   REAL(EB) :: Y2               !< Upper y boundary of the initialized region (m)
   REAL(EB) :: Z1               !< Lower z boundary of the initialized region (m)
   REAL(EB) :: Z2               !< Upper z boundary of the initialized region (m)
   REAL(EB) :: MASS_PER_VOLUME  !< Mass per unit volume of particles (kg/m3)
   REAL(EB) :: MASS_PER_TIME    !< Mass (kg) per time (s) of particles
   REAL(EB) :: DT_INSERT        !< Time increment between particle inserts (s)
   REAL(EB) :: T_INSERT         !< Time to start inserting particles (s)
   REAL(EB) :: X0               !< x origin of initialization region (m)
   REAL(EB) :: Y0               !< y origin of initialization region (m)
   REAL(EB) :: Z0               !< z origin of initialization region (m)
   REAL(EB) :: U0               !< Initial u component of velocity of particles (m/s)
   REAL(EB) :: V0               !< Initial v component of velocity of particles (m/s)
   REAL(EB) :: W0               !< Initial w component of velocity of particles (m/s)
   REAL(EB) :: VOLUME           !< Volume of region (m3)
   REAL(EB) :: HRRPUV=0._EB     !< Heat Release Rate Per Unit Volume (W/m3)
   REAL(EB) :: DX=0._EB         !< Spacing (m) of an array of particles
   REAL(EB) :: DY=0._EB         !< Spacing (m) of an array of particles
   REAL(EB) :: DZ=0._EB         !< Spacing (m) of an array of particles
   REAL(EB) :: HEIGHT           !< Height of initialization region (m)
   REAL(EB) :: RADIUS           !< Radius of initialization region, like a cone (m)
   REAL(EB) :: DIAMETER=-1._EB  !< Diameter of liquid droplets specified on an INIT line (m)
   REAL(EB) :: PARTICLE_WEIGHT_FACTOR !< Multiplicative factor for particles specified on the INIT line
   REAL(EB) :: PACKING_RATIO    !< Volume of particles divided by the volume of gas
   REAL(EB) :: CHI_R            !< Radiative fraction of HRRPUV
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: PARTICLE_INSERT_CLOCK  !< Time of last particle insertion (s)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MASS_FRACTION          !< Mass fraction of gas components
   INTEGER  :: PART_INDEX=0     !< Particle class index of inserted particles
   INTEGER  :: N_PARTICLES      !< Number of particles to insert
   INTEGER  :: DEVC_INDEX=0     !< Index of the device that uses this INITIALIZATION variable
   INTEGER  :: CTRL_INDEX=0     !< Index of the controller that uses this INITIALIZATION variable
   INTEGER  :: N_PARTICLES_PER_CELL=0 !< Number of particles to insert in each cell
   INTEGER  :: PATH_RAMP_INDEX(3)=0   !< Ramp index of a particle path
   INTEGER  :: RAMP_Q_INDEX=0         !< Ramp of HRRPUV
   LOGICAL :: ADJUST_DENSITY=.FALSE.
   LOGICAL :: ADJUST_TEMPERATURE=.FALSE.
   LOGICAL :: SINGLE_INSERTION=.TRUE.
   LOGICAL :: CELL_CENTERED=.FALSE.
   LOGICAL :: UNIFORM=.FALSE.
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: ALREADY_INSERTED
   CHARACTER(LABEL_LENGTH) :: SHAPE
   CHARACTER(LABEL_LENGTH) :: DEVC_ID
   CHARACTER(LABEL_LENGTH) :: CTRL_ID
   CHARACTER(LABEL_LENGTH) :: ID
END TYPE INITIALIZATION_TYPE

TYPE (INITIALIZATION_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: INITIALIZATION


!> \brief Parameters used to initalize particles used as gas phase devices

TYPE INIT_RESERVED_TYPE
   INTEGER :: N_PARTICLES  !< Number of particles for the device, usually 1
   INTEGER :: DEVC_INDEX   !< The index of the device
END TYPE INIT_RESERVED_TYPE

TYPE (INIT_RESERVED_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: INIT_RESERVED


!> \brief Parameters associated with a pressure ZONE

TYPE P_ZONE_TYPE
   REAL(EB) :: X                                                   !< x coordinate of ZONE specifier (m)
   REAL(EB) :: Y                                                   !< y coordinate of ZONE specifier (m)
   REAL(EB) :: Z                                                   !< z coordinate of ZONE specifier (m)
   REAL(EB) :: DPSTAR=0._EB
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LEAK_AREA                !< Array of leak areas to other ZONEs
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LEAK_REFERENCE_PRESSURE  !< Array of leak reference pressures
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LEAK_PRESSURE_EXPONENT   !< Array of leak reference exponents
   INTEGER :: N_DUCTNODES                                          !< Number of duct nodes in the ZONE
   INTEGER :: MESH_INDEX=0                                         !< Index of the MESH where the ZONE is located
   LOGICAL :: EVACUATION=.FALSE.
   LOGICAL :: PERIODIC=.FALSE.                                     !< Indicator if the ZONE boundary is periodic
   INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_INDEX                !< Array of NODE indices connected to the ZONE
   CHARACTER(LABEL_LENGTH) :: ID                                   !< Identifier
END TYPE P_ZONE_TYPE

TYPE (P_ZONE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: P_ZONE


!> \brief Parameters associated with a MOVE line, used to rotate or translate OBSTructions

TYPE MOVEMENT_TYPE
   REAL(EB) :: SCALE              !< Overall scale factor
   REAL(EB) :: SCALEX             !< Scale factor for x dimension
   REAL(EB) :: SCALEY             !< Scale factor for y dimension
   REAL(EB) :: SCALEZ             !< Scale factor for z dimension
   REAL(EB) :: DX                 !< Translation in x direction (m)
   REAL(EB) :: DY                 !< Translation in y direction (m)
   REAL(EB) :: DZ                 !< Translation in z direction (m)
   REAL(EB) :: X0                 !< x origin of rotation axis (m)
   REAL(EB) :: Y0                 !< y origin of rotation axis (m)
   REAL(EB) :: Z0                 !< z origin of rotation axis (m)
   REAL(EB) :: ROTATION_ANGLE     !< Angle of rotation around AXIS (deg)
   REAL(EB) :: AXIS(3)            !< Vector originating at (X0,Y0,Z0) that defines axis of rotation
   REAL(EB) :: T34(3,4)           !< Transformation matrix
   INTEGER :: INDEX               !< Integer identifier
   CHARACTER(LABEL_LENGTH) :: ID  !< Character identifier
END TYPE MOVEMENT_TYPE

TYPE (MOVEMENT_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: MOVEMENT


!> \brief Parameters associated with a MULT line for multiplying OBSTs, VENTs, etc

TYPE MULTIPLIER_TYPE
   REAL(EB) :: DXB(6)                               !< Array of deltas for six coordinates (m)
   REAL(EB) :: DX0                                  !< Translation in x direction (m)
   REAL(EB) :: DY0                                  !< Translation in x direction (m)
   REAL(EB) :: DZ0                                  !< Translation in x direction (m)
   REAL(EB) :: FDS_AREA(6)=0._EB
   INTEGER  :: I_LOWER                              !< Lower bound of i index
   INTEGER  :: I_UPPER                              !< Upper bound of i index
   INTEGER  :: J_LOWER                              !< Lower bound of j index
   INTEGER  :: J_UPPER                              !< Upper bound of j index
   INTEGER  :: K_LOWER                              !< Lower bound of k index
   INTEGER  :: K_UPPER                              !< Upper bound of k index
   INTEGER  :: N_LOWER                              !< Lower bound of n index
   INTEGER  :: N_UPPER                              !< Upper bound of n index
   INTEGER  :: N_COPIES                             !< Total number of copies of the original item
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: SKIP   !< Flag to skip copies in one of three directions
   CHARACTER(LABEL_LENGTH) :: ID                    !< Identifying string
   LOGICAL :: SEQUENTIAL=.FALSE.                    !< Flag indicating of the MULT creates a series of copies or 3D array
END TYPE MULTIPLIER_TYPE

TYPE (MULTIPLIER_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: MULTIPLIER

TYPE FILTER_TYPE
   INTEGER :: RAMP_INDEX = -1
   CHARACTER(LABEL_LENGTH) :: ID,TABLE_ID
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: EFFICIENCY,MULTIPLIER
   REAL(EB) :: CLEAN_LOSS,LOADING_LOSS
END TYPE FILTER_TYPE

TYPE (FILTER_TYPE), DIMENSION(:), ALLOCATABLE,  TARGET :: FILTER

TYPE DUCTNODE_TYPE
   INTEGER :: FILTER_INDEX=-1, N_DUCTS,VENT_INDEX = -1, MESH_INDEX = -1,ZONE_INDEX=-1
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DUCT_INDEX
   CHARACTER(LABEL_LENGTH) :: ID,TABLE_ID,VENT_ID='null'
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: LOSS_ARRAY,FILTER_LOADING
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ,DIR,ZZ_V
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_OLD
   REAL(EB) :: LOSS, P,TMP=273.15_EB,RHO,RSUM,CP,XYZ(3),FILTER_LOSS,TMP_V,RHO_V,RSUM_V,CP_V
   REAL(EB) :: P_OLD,TMP_OLD,RHO_OLD,CP_OLD,RSUM_OLD
   LOGICAL :: UPDATED, READ_IN, FIXED, AMBIENT = .FALSE.,LEAKAGE=.FALSE.,VENT=.FALSE.
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: IN_MESH
END TYPE DUCTNODE_TYPE

TYPE (DUCTNODE_TYPE), DIMENSION(:), ALLOCATABLE,  TARGET :: DUCTNODE

TYPE NODE_BC_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_V
   REAL(EB) :: TMP_V,RHO_V,RSUM_V,CP_V, P
END TYPE NODE_BC_TYPE

TYPE DUCT_TYPE
   INTEGER :: NODE_INDEX(2)=-1,DEVC_INDEX=-1,CTRL_INDEX=-1,FAN_INDEX=-1,AIRCOIL_INDEX=-1,RAMP_INDEX=0,RAMP_LOSS_INDEX=0,&
              N_CELLS=1,DUCT_INTERP_TYPE_INDEX,SURF_INDEX=-1,LPC_INDEX=-1,QFAN_N=-1
   REAL(EB) :: AREA,AREA_INITIAL,COIL_Q=0._EB,CP_D,DIAMETER,DP_FAN(2)=0._EB,DX=-1._EB,LENGTH,MASS_FLOW_INITIAL=1.E6_EB,&
               RHO_D,ROUGHNESS,RSUM_D=0._EB,TAU=-1._EB,TMP_D=273.15_EB,TOTAL_LOSS=0._EB,VEL(4)=0._EB,VOLUME_FLOW=1.E6_EB,&
               VOLUME_FLOW_INITIAL=1.E6_EB,LOSS(2)=0._EB
   REAL(EB) :: CP_D_OLD,RHO_D_OLD,RSUM_D_OLD
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CP_C,RHO_C,TMP_C,ZZ,PART_INDEX,PART_CELL_INDEX
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CP_C_OLD,RHO_C_OLD,TMP_C_OLD,ZZ_OLD
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: ZZ_C
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: ZZ_C_OLD
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: VEL_SYSTEM
   LOGICAL :: ROUND = .TRUE.,SQUARE = .FALSE.,DAMPER = .FALSE.,DAMPER_OPEN = .TRUE.,FAN_OPERATING=.TRUE.,COIL_OPERATING=.TRUE.,&
              FIXED=.FALSE.,REVERSE=.FALSE.,UPDATED,LEAKAGE=.FALSE.,LEAK_ENTHALPY=.FALSE.,LOCALIZED_LEAKAGE=.FALSE.
   CHARACTER(LABEL_LENGTH) :: ID
   REAL(EB) :: FAN_ON_TIME = 1.E10_EB,COIL_ON_TIME = 1.E10_EB
END TYPE DUCT_TYPE

TYPE (DUCT_TYPE), DIMENSION(:), ALLOCATABLE,  TARGET :: DUCT

TYPE FAN_TYPE
   INTEGER :: FAN_TYPE,RAMP_INDEX,SPIN_INDEX=0
   REAL(EB) :: VOL_FLOW,MAX_FLOW,MAX_PRES,OFF_LOSS=0._EB,TAU=0._EB
   CHARACTER(LABEL_LENGTH) :: ID,FAN_RAMP
END TYPE FAN_TYPE

TYPE(FAN_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: FAN

TYPE AIRCOIL_TYPE
   REAL(EB) :: COOLANT_TEMPERATURE=293.15_EB,COOLANT_SPECIFIC_HEAT=4186._EB, EFFICIENCY, COOLANT_MASS_FLOW=-9999._EB, &
               FIXED_Q=-1.E10_EB,TAU=0._EB
   INTEGER :: RAMP_INDEX=0
   CHARACTER(LABEL_LENGTH) :: ID
END TYPE AIRCOIL_TYPE

TYPE(AIRCOIL_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: AIRCOIL

TYPE NETWORK_TYPE
   INTEGER :: N_DUCTS,N_DUCTNODES,N_MATRIX
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DUCT_INDEX,NODE_INDEX,MATRIX_INDEX
END TYPE NETWORK_TYPE

TYPE(NETWORK_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: NETWORK

TYPE DUCTRUN_TYPE
   INTEGER :: N_DUCTS,N_DUCTNODES,DUCT_INTERP_TYPE_INDEX,N_MATRIX_SYSTEM=-1,N_QFANS=-1,N_FANS
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DUCT_INDEX,NODE_INDEX,MATRIX_SYSTEM_INDEX,QFAN_INDEX
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: RHS_SYSTEM_2,LHS_SYSTEM_1,LHS_SYSTEM_2
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHS_SYSTEM_1
END TYPE DUCTRUN_TYPE

TYPE(DUCTRUN_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: DUCTRUN


END MODULE TYPES
