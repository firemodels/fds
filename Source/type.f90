!> \brief A collection of major derived types used in FDS.
!> \details There are several TYPEs that require special attention. WALL, CFACE, and PARTICLE TYPES reference other
!> derived types that start with BOUNDARY, like BOUNDARY_COORD or BOUNDARY_ONE_D. The number of real, integer, and logical scalar
!> and array components of these derived types are denoted like this, for example: N_BOUNDARY_COORD_SCALAR_INTEGERS. You must
!> adjust this value if you add or subtract components from the derived type. You should then use an existing component
!> as a guide and trace it through func.f90 to see how to initialize this component.

MODULE TYPES

USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS, ONLY : IAXIS,JAXIS,KAXIS,MAX_DIM,LOW_IND,HIGH_IND

#ifdef WITH_MKL
USE MKL_PARDISO
#endif /* WITH_MKL */

IMPLICIT NONE (TYPE,EXTERNAL)

!> \brief Arrays to hold derived type (WALL, CFACE, or PARTICLE) components for I/O or MPI exchanges

TYPE STORAGE_TYPE
   INTEGER :: N_STORAGE_SLOTS=0                      !< The second dimension of the REALS, INTEGERS and LOGICALS arrays
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: REALS    !< Array of reals
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INTEGERS  !< Array of integers
   LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LOGICALS  !< Array of logicals
END TYPE STORAGE_TYPE

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
   REAL(EB) :: LENGTH                     !< Cylinder or plate length used for POROUS_DRAG or SCREEN_DRAG (m)

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: R_CNF         !< Independent variable (radius) in particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CNF           !< Cumulative Number Fraction particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CVF           !< Cumulative Volume Fraction particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: BREAKUP_R_CNF !< R_CNF of new distribution after particle break-up
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: BREAKUP_CNF   !< CNF of new distribution after particle break-up
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: BREAKUP_CVF   !< CVF of new distribution after particle break-up
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: W_CNF         !< Weighting factor in particle size distribution
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: R50           !< Array of median particle diameters for Mie calculation

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
   INTEGER :: N_STORAGE_REALS=0           !< Number of reals to store for this particle class
   INTEGER :: N_STORAGE_INTEGERS=0        !< Number of integers to store for this particle class
   INTEGER :: N_STORAGE_LOGICALS=0        !< Number of logicals to store for this particle class
   INTEGER :: NEW_PARTICLE_INCREMENT=50   !< Number of new storage slots to allocate when NPLDIM is exceeded

   INTEGER, ALLOCATABLE, DIMENSION(:) :: STRATUM_INDEX_LOWER  !< Lower index of size distribution band
   INTEGER, ALLOCATABLE, DIMENSION(:) :: STRATUM_INDEX_UPPER  !< Upper index of size distribution band

   LOGICAL :: STATIC=.FALSE.                         !< Flag indicating if particles move or not
   LOGICAL :: MASSLESS_TRACER=.FALSE.                !< Flag indicating if particles are just tracers for visualization
   LOGICAL :: MASSLESS_TARGET=.FALSE.                !< Flag indicating if particles are just targets for an output quantity
   LOGICAL :: LIQUID_DROPLET=.FALSE.                 !< Flag indicating if particles are liquid droplets
   LOGICAL :: SOLID_PARTICLE=.FALSE.                 !< Flag indicating if particles are solid, not liquid
   LOGICAL :: MONODISPERSE=.FALSE.                   !< Flag indicating if particle size is monodisperse
   LOGICAL :: TURBULENT_DISPERSION=.FALSE.           !< Flag indicating if subgrid-scale turbulence is applied
   LOGICAL :: BREAKUP=.FALSE.                        !< Flag indicating if paricles or droplets break-up
   LOGICAL :: CHECK_DISTRIBUTION=.FALSE.             !< Flag indicating if diagnostic output on size distribution is specified
   LOGICAL :: FUEL=.FALSE.                           !< Flag indicating if droplets evaporate into fuel gas
   LOGICAL :: DUCT_PARTICLE=.FALSE.                  !< Flag indicating if particles can pass through a duct
   LOGICAL :: EMBER_PARTICLE=.FALSE.                 !< Flag indicating if particles can become flying embers
   LOGICAL :: TRACK_EMBERS=.TRUE.                    !< Flag indicating if flying embers are tracked or removed immediately
   LOGICAL :: ADHERE_TO_SOLID=.FALSE.                !< Flag indicating if particles can stick to a solid
   LOGICAL :: INCLUDE_BOUNDARY_COORD_TYPE=.TRUE.     !< This particle requires basic coordinate information
   LOGICAL :: INCLUDE_BOUNDARY_PROPS_TYPE=.FALSE.    !< This particle requires surface variables for heat and mass transfer
   LOGICAL :: INCLUDE_BOUNDARY_ONE_D_TYPE=.TRUE.     !< This particle requires in-depth 1-D conduction/reaction arrays
   LOGICAL :: INCLUDE_BOUNDARY_RADIA_TYPE=.FALSE.    !< This particle requires angular-specific radiation intensities

   TYPE(STORAGE_TYPE) :: PARTICLE_STORAGE

END TYPE LAGRANGIAN_PARTICLE_CLASS_TYPE

TYPE (LAGRANGIAN_PARTICLE_CLASS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: LAGRANGIAN_PARTICLE_CLASS


!> \brief Solid material density for 1-D pyrolysis/conduction algorithm

TYPE MATL_COMP_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO !< (1:NWP) Solid density (kg/m3)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_DOT !< (1:NWP) Change in solid density (kg/m3/s)
END TYPE MATL_COMP_TYPE

!> \brief Gas mass concentration in solid for 1-D mass transfer

TYPE SPEC_COMP_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_ZZ !< (0:NWP+1) Gas concentratoin (kg/m3)
END TYPE SPEC_COMP_TYPE

!> \brief Radiation intensity at a boundary for a given wavelength band

TYPE BAND_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ILW !< (1:NRA) Radiation intensity (W/m2/sr)
END TYPE BAND_TYPE

!> \brief Coordinate variables associated with a WALL or CFACE boundary cell
!> \details If you change the number of scalar variables in BOUNDARY_COORD_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_BOUNDARY_COORD_SCALAR_REALS=3     !< Number of scalar reals in BOUNDARY_COORD
INTEGER, PARAMETER :: N_BOUNDARY_COORD_SCALAR_INTEGERS=7  !< Number of scalar integers in BOUNDARY_COORD
INTEGER, PARAMETER :: N_BOUNDARY_COORD_SCALAR_LOGICALS=0  !< Number of scalar logicals in BOUNDARY_COORD
INTEGER :: N_BOUNDARY_COORD_STORAGE_REALS                 !< Total number of reals in BOUNDARY_COORD
INTEGER :: N_BOUNDARY_COORD_STORAGE_INTEGERS              !< Total number of integers in BOUNDARY_COORD
INTEGER :: N_BOUNDARY_COORD_STORAGE_LOGICALS              !< Total number of logicals in BOUNDARY_COORD

TYPE BOUNDARY_COORD_TYPE

   INTEGER :: II             !< Ghost cell \f$ x \f$ index
   INTEGER :: JJ             !< Ghost cell \f$ y \f$ index
   INTEGER :: KK             !< Ghost cell \f$ z \f$ index
   INTEGER :: IIG            !< Gas cell \f$ x \f$ index
   INTEGER :: JJG            !< Gas cell \f$ y \f$ index
   INTEGER :: KKG            !< Gas cell \f$ z \f$ index
   INTEGER :: IOR=0          !< Index of orientation of the WALL cell

   REAL(EB) :: X             !< \f$ x \f$ coordinate of boundary cell center
   REAL(EB) :: Y             !< \f$ y \f$ coordinate of boundary cell center
   REAL(EB) :: Z             !< \f$ z \f$ coordinate of boundary cell center

END TYPE BOUNDARY_COORD_TYPE


!> \brief Variables associated with a WALL, PARTICLE, or CFACE boundary cell
!> \details If you change the number of scalar variables in BOUNDARY_ONE_D_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_ONE_D_SCALAR_REALS=27
INTEGER, PARAMETER :: N_ONE_D_SCALAR_INTEGERS=4
INTEGER, PARAMETER :: N_ONE_D_SCALAR_LOGICALS=1

TYPE BOUNDARY_ONE_D_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: M_DOT_G_PP_ACTUAL   !< (1:N_TRACKED_SPECIES) Actual mass production rate per unit area
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: M_DOT_S_PP          !< (1:SF\%N_MATL) Mass production rate of solid species
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: M_DOT_G_PP_ADJUST   !< (1:N_TRACKED_SPECIES) Adjusted mass production rate per unit area
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: X                   !< (0:NWP) Depth (m), \f$ x_{{\rm s},i} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: TMP                 !< Temperature in center of each solid cell, \f$ T_{{\rm s},i} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LAYER_THICKNESS     !< (1:SF\%N_LAYERS) Thickness of layer (m)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_F                !< (1:N_TRACKED_SPECIES) Species mixture mass fraction at surface
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_D_F             !< (1:N_TRACKED_SPECIES) Diffusion at surface, \f$ \rho D_\alpha \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_D_DZDN_F        !< \f$ \rho D_\alpha \partial Z_\alpha / \partial n \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_C_S             !< Solid density times specific heat (J/m3/K)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: K_S                 !< Solid conductivity (W/m/K)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: AWM_AEROSOL         !< Accumulated aerosol mass per unit area (kg/m2)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DDSUM               !< Scaling factor to get minimum cell size
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: SMALLEST_CELL_SIZE  !< Minimum cell size (m)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: PART_MASS           !< Accumulated mass of particles waiting to be injected (kg/m2)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: PART_ENTHALPY       !< Accumulated enthalpy of particles waiting to be injected (kJ/m2)

   TYPE(MATL_COMP_TYPE), ALLOCATABLE, DIMENSION(:) :: MATL_COMP !< (1:SF\%N_MATL) Material component
   TYPE(SPEC_COMP_TYPE), ALLOCATABLE, DIMENSION(:) :: SPEC_COMP !< (1:SF\%N_SPEC) Gas component

   INTEGER, ALLOCATABLE, DIMENSION(:) :: N_LAYER_CELLS              !< (1:SF\%N_LAYERS) Number of cells in the layer

   INTEGER :: SURF_INDEX=-1    !< SURFACE index
   INTEGER :: PRESSURE_ZONE=0  !< Pressure ZONE of the adjacent gas phase cell
   INTEGER :: NODE_INDEX=0     !< HVAC node index associated with surface
   INTEGER :: N_SUBSTEPS=1     !< Number of substeps in the 1-D conduction/reaction update

   REAL(EB) :: AREA=0._EB            !< Face area (m2)
   REAL(EB) :: HEAT_TRANS_COEF=0._EB !< Heat transfer coefficient (W/m2/K)
   REAL(EB) :: Q_CON_F=0._EB         !< Convective heat flux at surface (W/m2)
   REAL(EB) :: Q_RAD_IN=0._EB        !< Incoming radiative flux (W/m2)
   REAL(EB) :: Q_RAD_OUT=0._EB       !< Outgoing radiative flux (W/m2)
   REAL(EB) :: EMISSIVITY=1._EB      !< Surface emissivity
   REAL(EB) :: AREA_ADJUST=1._EB     !< Ratio of actual surface area to grid cell face area
   REAL(EB) :: T_IGN=0._EB           !< Ignition time (s)
   REAL(EB) :: TMP_F                 !< Surface temperature (K)
   REAL(EB) :: TMP_F_OLD             !< Holding value for surface temperature (K)
   REAL(EB) :: TMP_B                 !< Back surface temperature (K)
   REAL(EB) :: U_NORMAL=0._EB        !< Normal component of velocity (m/s) at surface, start of time step
   REAL(EB) :: U_NORMAL_S=0._EB      !< Estimated normal component of velocity (m/s) at next time step
   REAL(EB) :: U_NORMAL_0=0._EB      !< Initial or specified normal component of velocity (m/s) at surface
   REAL(EB) :: U_TANG=0._EB          !< Tangential velocity (m/s) near surface
   REAL(EB) :: RHO_F                 !< Gas density at the wall (kg/m3)
   REAL(EB) :: RDN=1._EB             !< \f$ 1/ \delta n \f$ at the surface (1/m)
   REAL(EB) :: K_G=0.1_EB            !< Thermal conductivity, \f$ k \f$, in adjacent gas phase cell
   REAL(EB) :: Q_DOT_G_PP=0._EB      !< Heat release rate per unit area (W/m2)
   REAL(EB) :: Q_DOT_O2_PP=0._EB     !< Heat release rate per unit area (W/m2) due to oxygen consumption
   REAL(EB) :: Q_CONDENSE=0._EB      !< Heat release rate per unit area (W/m2) due to gas condensation
   REAL(EB) :: BURN_DURATION=0._EB   !< Duration of a specified fire (s)
   REAL(EB) :: T_SCALE=0._EB         !< Scaled time for a surface with CONE_HEAT_FLUX (s)
   REAL(EB) :: Q_SCALE=0._EB         !< Scaled integrated heat release for a surface with CONE_HEAT_FLUX
   REAL(EB) :: T_MATL_PART=0._EB     !< Time interval for current value in PART_MASS and PART_ENTHALPY arrays (s)
   REAL(EB) :: B_NUMBER=0._EB        !< B number for droplet or wall
   REAL(EB) :: M_DOT_PART_ACTUAL     !< Mass flux of all particles (kg/m2/s)

   LOGICAL :: BURNAWAY=.FALSE.       !< Indicater if cell can burn away when fuel is exhausted

END TYPE BOUNDARY_ONE_D_TYPE

!> \brief Property variables associated with a WALL or CFACE boundary cell
!> \details If you change the number of scalar variables in BOUNDARY_PROPS_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_BOUNDARY_PROPS_SCALAR_REALS=7
INTEGER, PARAMETER :: N_BOUNDARY_PROPS_SCALAR_INTEGERS=1
INTEGER, PARAMETER :: N_BOUNDARY_PROPS_SCALAR_LOGICALS=0
INTEGER, DIMENSION(10) :: BOUNDARY_PROPS_REALS_ARRAY_SIZE=0, &
                          BOUNDARY_PROPS_INTEGERS_ARRAY_SIZE=0, &
                          BOUNDARY_PROPS_LOGICALS_ARRAY_SIZE=0
INTEGER :: N_BOUNDARY_PROPS_STORAGE_REALS,N_BOUNDARY_PROPS_STORAGE_INTEGERS,N_BOUNDARY_PROPS_STORAGE_LOGICALS

TYPE BOUNDARY_PROPS_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: A_LP_MPUA           !< Accumulated liquid droplet mass per unit area (kg/m2)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LP_CPUA             !< Liquid droplet cooling rate unit area (W/m2)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LP_MPUA             !< Liquid droplet mass per unit area (kg/m2)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LP_TEMP             !< Liquid droplet mean temperature (K)

   REAL(EB) :: U_TAU=0._EB           !< Friction velocity (m/s)
   REAL(EB) :: Y_PLUS=1._EB          !< Dimensionless boundary layer thickness unit
   REAL(EB) :: Z_STAR=1._EB          !< Dimensionless boundary layer unit
   REAL(EB) :: PHI_LS=-1._EB         !< Level Set value for output only
   REAL(EB) :: WORK1=0._EB           !< Work array
   REAL(EB) :: WORK2=0._EB           !< Work array
   REAL(EB) :: K_SUPPRESSION=0._EB   !< Suppression coefficent (m2/kg/s)

   INTEGER  :: SURF_INDEX=-1         !< Surface index

END TYPE BOUNDARY_PROPS_TYPE


!> \brief Angular radiation intensities associated with a WALL, CFACE, or LAGRANGIAN_PARTICLE
!> \details If you change the number of scalar variables in BOUNDARY_RADIA_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_BOUNDARY_RADIA_SCALAR_REALS=0
INTEGER, PARAMETER :: N_BOUNDARY_RADIA_SCALAR_INTEGERS=0
INTEGER, PARAMETER :: N_BOUNDARY_RADIA_SCALAR_LOGICALS=0
INTEGER, DIMENSION(10) :: BOUNDARY_RADIA_REALS_ARRAY_SIZE=0, &
                          BOUNDARY_RADIA_INTEGERS_ARRAY_SIZE=0, &
                          BOUNDARY_RADIA_LOGICALS_ARRAY_SIZE=0
INTEGER :: N_BOUNDARY_RADIA_STORAGE_REALS,N_BOUNDARY_RADIA_STORAGE_INTEGERS,N_BOUNDARY_RADIA_STORAGE_LOGICALS

TYPE BOUNDARY_RADIA_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: IL              !< (1:NSB) Radiance (W/m2/sr); output only
   TYPE(BAND_TYPE), ALLOCATABLE, DIMENSION(:) :: BAND     !< (1:NSB) Radiation wavelength band
END TYPE BOUNDARY_RADIA_TYPE


!> \brief Variables associated with a single Lagrangian particle
!> \details If you change the number of scalar variables in LAGRANGIAN_PARTICLE_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_PARTICLE_SCALAR_REALS=14
INTEGER, PARAMETER :: N_PARTICLE_SCALAR_INTEGERS=14
INTEGER, PARAMETER :: N_PARTICLE_SCALAR_LOGICALS=4

TYPE LAGRANGIAN_PARTICLE_TYPE

   INTEGER :: LP_INDEX=0             !< Self-identifier
   INTEGER :: BC_INDEX=0             !< Coordinate variables
   INTEGER :: OD_INDEX=0             !< Variables devoted to 1-D heat conduction in depth
   INTEGER :: BP_INDEX=0             !< Variables devoted to surface properties
   INTEGER :: BR_INDEX=0             !< Variables devoted to radiation intensities
   INTEGER :: TAG                    !< Unique integer identifier for the particle
   INTEGER :: CLASS_INDEX=0          !< LAGRANGIAN_PARTICLE_CLASS of particle
   INTEGER :: ORIENTATION_INDEX=0    !< Index in the array of all ORIENTATIONs
   INTEGER :: WALL_INDEX=0           !< If liquid droplet has stuck to a wall, this is the WALL cell index
   INTEGER :: DUCT_INDEX=0           !< Index of duct
   INTEGER :: INIT_INDEX=0           !< Index of INIT line
   INTEGER :: DUCT_CELL_INDEX=0      !< Index of duct cell
   INTEGER :: CFACE_INDEX=0          !< Index of immersed boundary CFACE that the droplet has attached to
   INTEGER :: PROP_INDEX=0           !< Index of the PROPERTY_TYPE assigned to the particle

   LOGICAL :: SHOW=.FALSE.         !< Show the particle in Smokeview
   LOGICAL :: SPLAT=.FALSE.        !< The liquid droplet has hit a solid
   LOGICAL :: EMBER=.FALSE.        !< The particle can break away and become a burning ember
   LOGICAL :: PATH_PARTICLE =.FALSE. !< Indicator of a particle with a specified path

   REAL(EB) :: U=0._EB             !< \f$ x \f$ velocity component of particle (m/s)
   REAL(EB) :: V=0._EB             !< \f$ y \f$ velocity component of particle (m/s)
   REAL(EB) :: W=0._EB             !< \f$ z \f$ velocity component of particle (m/s)
   REAL(EB) :: PWT=1._EB           !< Weight factor of particle; i.e. the number of real particles it represents
   REAL(EB) :: ACCEL_X=0._EB       !< Contribution to acceleration of gas in \f$ x \f$ direction (m/s2)
   REAL(EB) :: ACCEL_Y=0._EB       !< Contribution to acceleration of gas in \f$ y \f$ direction (m/s2)
   REAL(EB) :: ACCEL_Z=0._EB       !< Contribution to acceleration of gas in \f$ z \f$ direction (m/s2)
   REAL(EB) :: RE=0._EB            !< Reynolds number based on particle diameter
   REAL(EB) :: MASS=0._EB          !< Particle mass (kg)
   REAL(EB) :: T_INSERT=0._EB      !< Time when particle was inserted (s)
   REAL(EB) :: DX=1.               !< Length factor used in POROUS_DRAG calculation (m)
   REAL(EB) :: DY=1.               !< Length factor used in POROUS_DRAG calculation (m)
   REAL(EB) :: DZ=1.               !< Length factor used in POROUS_DRAG calculation (m)
   REAL(EB) :: C_DRAG=0._EB        !< Drag coefficient

END TYPE LAGRANGIAN_PARTICLE_TYPE


!> \brief Variables associated with a WALL cell
!> \details If you change the number of scalar variables in WALL_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_WALL_SCALAR_REALS=4
INTEGER, PARAMETER :: N_WALL_SCALAR_INTEGERS=18
INTEGER, PARAMETER :: N_WALL_SCALAR_LOGICALS=0

TYPE WALL_TYPE

   REAL(EB) :: DUNDT=0._EB            !< \f$ \partial u_n / \partial t \f$
   REAL(EB) :: Q_LEAK=0._EB           !< Heat production of leaking gas (W/m3)
   REAL(EB) :: V_DEP=0._EB            !< Deposition velocity (m/s)
   REAL(EB) :: VEL_ERR_NEW=0._EB      !< Velocity mismatch at mesh or solid boundary (m/s)

   INTEGER :: WALL_INDEX=0            !< Index of itself -- used to determine if the WALL cell has been assigned
   INTEGER :: BC_INDEX=0              !< Index within the array BOUNDARY_COORD
   INTEGER :: OD_INDEX=0              !< Index within the array BOUNDARY_ONE_D
   INTEGER :: BP_INDEX=0              !< Index within the array BOUNDARY_PROPS
   INTEGER :: BR_INDEX=0              !< Index within the array BOUNDARY_RADIA
   INTEGER :: SURF_INDEX=0            !< Index of the SURFace conditions
   INTEGER :: BACK_INDEX=0            !< WALL index of back side of obstruction or exterior wall cell
   INTEGER :: BACK_MESH               !< Mesh number on back side of obstruction or exterior wall cell
   INTEGER :: BOUNDARY_TYPE=0         !< Descriptor: SOLID, MIRROR, OPEN, INTERPOLATED, etc
   INTEGER :: SURF_INDEX_ORIG=0       !< Original SURFace index for this cell
   INTEGER :: OBST_INDEX=0            !< Index of the OBSTruction
   INTEGER :: PRESSURE_BC_INDEX       !< Poisson boundary condition, NEUMANN or DIRICHLET
   INTEGER :: VENT_INDEX=0            !< Index of the VENT containing this cell
   INTEGER :: JD11_INDEX=0
   INTEGER :: JD12_INDEX=0
   INTEGER :: JD21_INDEX=0
   INTEGER :: JD22_INDEX=0
   INTEGER :: CUT_FACE_INDEX=0

END TYPE WALL_TYPE


!> \brief Variables associated with the external boundary of a mesh

TYPE EXTERNAL_WALL_TYPE
   INTEGER :: NOM                                     !< Number of the adjacent (Other) Mesh
   INTEGER :: NIC_MIN                                 !< Start of indices for the cell in the other mesh
   INTEGER :: NIC_MAX                                 !< End of indices for the cell in the other mesh
   INTEGER :: NIC                                     !< NIC_MAX-NIC_MIN
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
   CHARACTER(LABEL_LENGTH) :: FUEL        !< Name of reaction fuel species
   CHARACTER(LABEL_LENGTH) :: OXIDIZER    !< Name of reaction oxidizer (lumped) species
   CHARACTER(LABEL_LENGTH) :: PRODUCTS    !< Name of reaction product (lumped) species
   CHARACTER(LABEL_LENGTH) :: ID          !< Identifer of reaction
   CHARACTER(LABEL_LENGTH) :: RAMP_CHI_R  !< Name of ramp for radiative fraction
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: SPEC_ID_NU       !< Array of species names corresponding to stoich coefs
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: SPEC_ID_NU_READ  !< Holding array for SPEC_ID_NU
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: SPEC_ID_N_S      !< Array of finite rate species exponents
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: SPEC_ID_N_S_READ !< Holding array of finite rate species exponents
   CHARACTER(MESSAGE_LENGTH) :: FYI='null'  !< User comment
   CHARACTER(FORMULA_LENGTH) :: EQUATION    !< Reaction equation
   CHARACTER(LABEL_LENGTH) :: FWD_ID        !< ID of forward reaction
   REAL(EB) :: C                            !< Number of carbon atoms in the fuel molecule (SIMPLE_CHEMISTRY)
   REAL(EB) :: H                            !< Number of hydrogen atoms in the fuel molecule (SIMPLE_CHEMISTRY)
   REAL(EB) :: N                            !< Number of nitrogen atoms in the fuel molecule (SIMPLE_CHEMISTRY)
   REAL(EB) :: O                            !< Number of oxygen atoms in the fuel molecule (SIMPLE_CHEMISTRY)
   REAL(EB) :: EPUMO2                       !< Energy Per Unit Mass Oxygen consumed (J/kg)
   REAL(EB) :: HEAT_OF_COMBUSTION           !< Energy per unit mass fuel consumed (J/kg)
   REAL(EB) :: HOC_COMPLETE                 !< Complete heat of combustion for two step SIMPLE_CHEMISTRY (J/kg)
   REAL(EB) :: A_PRIME                      !< Adjusted pre-exponential reaction kinetic parameter
   REAL(EB) :: A_IN                         !< Unajusted pre-exponential reaction kinetic parameter
   REAL(EB) :: E                            !< Activation energy (J/kmol)
   REAL(EB) :: E_IN                         !< User-specified activation energy (J/mol)
   REAL(EB) :: MW_FUEL                      !< Molecular weight of fuel (g/mol)
   REAL(EB) :: MW_SOOT                      !< Molecular weight of soot surrogate gas (g/mol)
   REAL(EB) :: Y_O2_MIN                     !< Lower oxygen limit in terms of mass fraction
   REAL(EB) :: CO_YIELD                     !< CO yield in SIMPLE_CHEMISTRY model
   REAL(EB) :: SOOT_YIELD                   !< Soot yield in SIMPLE_CHEMISTRY model
   REAL(EB) :: H2_YIELD                     !< H2 yield in SIMPLE_CHEMISTRY model
   REAL(EB) :: HCN_YIELD                    !< HCN yield in SIMPLE_CHEMISTRY model
   REAL(EB) :: SOOT_H_FRACTION              !< Mass fraction of hydrogen within soot
   REAL(EB) :: RHO_EXPONENT                 !< Exponent of density in reaction expression
   REAL(EB) :: CRIT_FLAME_TMP               !< Critical Flame Temperature (K)
   REAL(EB) :: AUTO_IGNIT_TMP               !< Reaction specific Auto Ignition Temperature (K)
   REAL(EB) :: NU_O2=0._EB                  !< Oxygen coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: NU_N2=0._EB                  !< Nitrogen coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: NU_H2O=0._EB                 !< Water coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: NU_H2=0._EB                  !< Hydrogen coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: NU_HCN=0._EB                 !< Hydrogen cyanide coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: NU_CO2=0._EB                 !< Carbon dioxide coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: NU_CO=0._EB                  !< Carbon monoxide coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: NU_SOOT=0._EB                !< Soot coefficient in SIMPLE_CHEMISTRY model
   REAL(EB) :: S=0._EB                      !< Stoichiometric coefficient for MIXTURE FRACTION output
   REAL(EB) :: N_T=0._EB                    !< Temperature exponent in reaction expression
   REAL(EB) :: CHI_R                        !< Radiative fraction
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: NU              !< Array of stoichiometric coefficients for lumped species equation
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: NU_READ         !< Holding array of stoichiometric coefficients
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: NU_SPECIES      !< Array of stoichiometric coefficients for primitive species equation
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: N_S             !< Array of species exponents
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: N_S_READ        !< Holding array of species exponents
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: NU_MW_O_MW_F    !< Species mol. weight times stoich. coef. over fuel MW
   INTEGER :: FUEL_SMIX_INDEX=-1            !< Lumped species index for fuel
   INTEGER :: AIR_SMIX_INDEX=-1             !< Lumped species index for air
   INTEGER :: N_SMIX                        !< Number of lumped species in reaction equation
   INTEGER :: N_SPEC                        !< Number of primitive species in reaction equation
   INTEGER :: RAMP_CHI_R_INDEX=0            !< Index of radiative fraction ramp
   INTEGER :: PRIORITY=1                    !< Index used in fast-fast SIMPLE_CHEMISTRY two step reaction
   LOGICAL :: IDEAL                         !< Indicator that the given HEAT_OF_COMBUSTION is the ideal value
   LOGICAL :: CHECK_ATOM_BALANCE            !< Indicator for diagnostic output
   LOGICAL :: FAST_CHEMISTRY=.FALSE.        !< Indicator of fast reaction
   LOGICAL :: REVERSE=.FALSE.               !< Indicator of a reverse reaction
   LOGICAL :: THIRD_BODY=.FALSE.            !< Indicator of catalyst
END TYPE REACTION_TYPE

TYPE (REACTION_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: REACTION

TYPE MATERIAL_TYPE
   REAL(EB) :: RHO_S                                    !< Density (kg/m3) of the pure material
   REAL(EB) :: EMISSIVITY                               !< Emissivity of surface
   REAL(EB) :: THERMAL_DIFFUSIVITY                      !< Thermal diffusivity (m2/s)
   REAL(EB) :: KAPPA_S                                  !< Absorption coefficient (1/m)
   REAL(EB) :: TMP_BOIL                                 !< Boiling temperature (K) of a liquid
   REAL(EB) :: REFRACTIVE_INDEX
   REAL(EB) :: POROSITY=0._EB                           !< Porosity
   REAL(EB) :: MW=-1._EB                                !< Molecular weight (g/mol)
   REAL(EB) :: HEAT_OF_GASIFICATION                     !< Heat of gasification (J/kg)
   INTEGER :: PYROLYSIS_MODEL                           !< Type of pyrolysis model (SOLID, LIQUID, VEGETATION)
   CHARACTER(LABEL_LENGTH) :: ID                        !< Identifier
   CHARACTER(LABEL_LENGTH) :: RAMP_H_R(MAX_REACTIONS)   !< Name of RAMP for Heat of Reaction
   CHARACTER(LABEL_LENGTH) :: RAMP_K_S                  !< Name of RAMP for thermal conductivity of solid
   CHARACTER(LABEL_LENGTH) :: RAMP_C_S                  !< Name of RAMP for specific heat of solid
   INTEGER :: N_REACTIONS                               !< Number of solid phase reactions
   INTEGER :: PROP_INDEX=-1
   INTEGER :: I_RAMP_K_S=-1                             !< Index of conductivity RAMP
   INTEGER :: I_RAMP_C_S=-1                             !< Index of specific heat RAMP
   INTEGER, DIMENSION(MAX_REACTIONS) :: N_RESIDUE       !< Number of residue materials
   INTEGER, DIMENSION(MAX_REACTIONS) :: I_RAMP_H_R=-1   !< Index of the heat of reaction RAMP
   INTEGER, DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: RESIDUE_MATL_INDEX !< Index of the residue
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LPC_INDEX    !< Lagrangian Particle Class for material conversion
   INTEGER, DIMENSION(3) :: RGB                         !< Color indices for material
   REAL(EB), DIMENSION(MAX_REACTIONS) :: TMP_REF        !< Reference temperature used for calculating kinetic constants (K)
   REAL(EB), DIMENSION(MAX_REACTIONS) :: RATE_REF       !< Reference rate used for calculating kinetic constants (1/s)
   REAL(EB), DIMENSION(MAX_REACTIONS) :: MAX_REACTION_RATE !< Maximum reaction rate (kg/m3/s)
   REAL(EB), DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: NU_RESIDUE=0._EB !< Mass stoichiometric coefficient of residue
   REAL(EB), DIMENSION(MAX_REACTIONS) :: A              !< Pre-exponential constant (1/s)
   REAL(EB), DIMENSION(MAX_REACTIONS) :: E              !< Activation energy (J/kmol)
   REAL(EB), DIMENSION(MAX_REACTIONS) :: N_S            !< Reaction order
   REAL(EB), DIMENSION(MAX_REACTIONS) :: N_T            !< Optional exponent for temperature in reaction expression
   REAL(EB), DIMENSION(MAX_REACTIONS) :: N_O2           !< Optional exponent for oxygen term in reaction expression
   REAL(EB), DIMENSION(MAX_REACTIONS) :: GAS_DIFFUSION_DEPTH !< Length scale used in char oxidation calculation (m)
   REAL(EB), DIMENSION(MAX_REACTIONS) :: NU_O2_CHAR     !< Mass stoichiometric coefficient for oxygen in char reaaction
   REAL(EB), DIMENSION(MAX_REACTIONS) :: BETA_CHAR      !< Constant used in char oxidation model
   REAL(EB), DIMENSION(MAX_REACTIONS) :: HEATING_RATE   !< Heating rate (K/s) used in calculation of kinetic constants
   REAL(EB), DIMENSION(MAX_REACTIONS) :: PYROLYSIS_RANGE !< Temperature range (K) over which pyrolysis occurs
   REAL(EB), DIMENSION(MAX_LPC,MAX_REACTIONS) :: NU_PART !< Mass stoichiometric coefficient that dictates fraction of mass to parts
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: NU_GAS       !< Mass stoichiometric coefficient for solid to gas conversion
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: ADJUST_BURN_RATE !< Adjustment to pyrolysis rate to account for different HoC
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: NU_LPC       !< Mass stoichiometric coefficient for particles
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: H_R          !< Heat of Reaction (J/kg)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DIFFUSIVITY_GAS
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: H
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: K_S
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: C_S
   REAL(EB), DIMENSION(MAX_SPECIES,MAX_REACTIONS) :: NU_SPEC
   REAL(EB), DIMENSION(MAX_SPECIES,MAX_REACTIONS) :: HEAT_OF_COMBUSTION
   REAL(EB), DIMENSION(MAX_SPECIES) :: DIFFUSIVITY_SPEC
   LOGICAL :: ALLOW_SHRINKING
   LOGICAL :: ALLOW_SWELLING
   LOGICAL :: CONST_C=.TRUE.
   CHARACTER(LABEL_LENGTH), DIMENSION(MAX_MATERIALS,MAX_REACTIONS) :: RESIDUE_MATL_NAME
   CHARACTER(LABEL_LENGTH), DIMENSION(MAX_SPECIES,MAX_REACTIONS) :: SPEC_ID
   CHARACTER(LABEL_LENGTH), DIMENSION(MAX_LPC,MAX_REACTIONS) :: PART_ID
   CHARACTER(MESSAGE_LENGTH) :: FYI='null'
END TYPE MATERIAL_TYPE

TYPE (MATERIAL_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: MATERIAL

!> \brief Variables associated with a surface type

TYPE SURFACE_TYPE

   REAL(EB) :: AREA_MULTIPLIER=1._EB                     !< Factor for manual surface area adjustment
   REAL(EB) :: TMP_FRONT=-1._EB                          !< Specified front surface temperture (K)
   REAL(EB) :: TMP_BACK=-1._EB                           !< Specified back surface gas temperature (K)
   REAL(EB) :: TMP_INNER_HT3D=-1._EB                     !< Specified inner temperature for 3D heating (K)
   REAL(EB) :: VEL                                       !< Specified normal velocity (m/s)
   REAL(EB) :: VEL_GRAD
   REAL(EB) :: PLE                                       !< Exponent for boundary layer velocity profile
   REAL(EB) :: Z0                                        !< Reference height for boundary layer profile (m)
   REAL(EB) :: Z_0                                       !< Surface roughness (m)
   REAL(EB) :: NEAR_WALL_EDDY_VISCOSITY                  !< Used with CONSTANT KINEMATIC EDDY VISCOSITY (m2/s)
   REAL(EB) :: CONVECTIVE_HEAT_FLUX                      !< Specified convective heat flux at surface (W/m2)
   REAL(EB) :: NET_HEAT_FLUX                             !< Specified net heat flux at surface (W/m2)
   REAL(EB) :: VOLUME_FLOW                               !< Specified volume flow (m3/s)
   REAL(EB) :: HRRPUA                                    !< Specified Heat Release Rate Per Unit Volume (W/m3)
   REAL(EB) :: MLRPUA                                    !< Specified Mass Loss Rate Per Unit Area (kg/m2/s)
   REAL(EB) :: T_IGN                                     !< Specified ignition time (s)
   REAL(EB) :: SURFACE_DENSITY                           !< Mass per unit area (kg/m2)
   REAL(EB) :: CELL_SIZE_FACTOR
   REAL(EB) :: E_COEFFICIENT
   REAL(EB) :: TEXTURE_WIDTH
   REAL(EB) :: TEXTURE_HEIGHT
   REAL(EB) :: THICKNESS
   REAL(EB) :: EXTERNAL_FLUX
   REAL(EB) :: DXF
   REAL(EB) :: DXB
   REAL(EB) :: MASS_FLUX_TOTAL
   REAL(EB) :: PARTICLE_MASS_FLUX
   REAL(EB) :: EMISSIVITY
   REAL(EB) :: MAX_PRESSURE
   REAL(EB) :: TMP_IGN
   REAL(EB) :: TMP_EXT
   REAL(EB) :: H_V
   REAL(EB) :: LAYER_DIVIDE
   REAL(EB) :: ROUGHNESS
   REAL(EB) :: LENGTH=-1._EB
   REAL(EB) :: WIDTH=-1._EB
   REAL(EB) :: DT_INSERT
   REAL(EB) :: H_FIXED=-1._EB
   REAL(EB) :: H_FIXED_B=-1._EB
   REAL(EB) :: HM_FIXED=-1._EB
   REAL(EB) :: EMISSIVITY_BACK
   REAL(EB) :: CONV_LENGTH
   REAL(EB) :: XYZ(3)
   REAL(EB) :: FIRE_SPREAD_RATE
   REAL(EB) :: MINIMUM_LAYER_THICKNESS
   REAL(EB) :: INNER_RADIUS=0._EB
   REAL(EB) :: MASS_FLUX_VAR=-1._EB
   REAL(EB) :: VEL_BULK
   REAL(EB) :: VEL_PART
   REAL(EB) :: PARTICLE_SURFACE_DENSITY=-1._EB
   REAL(EB) :: DRAG_COEFFICIENT=2.8_EB
   REAL(EB) :: SHAPE_FACTOR=0.25_EB
   REAL(EB) :: MINIMUM_BURNOUT_TIME=1.E6_EB
   REAL(EB) :: DELTA_TMP_MAX=10._EB
   REAL(EB) :: BURN_DURATION=1.E6_EB
   REAL(EB) :: CONE_HEAT_FLUX=-1._EB
   REAL(EB) :: PARTICLE_EXTRACTION_VELOCITY=1.E6_EB
   REAL(EB) :: INIT_PER_AREA=0._EB

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DX,RDX,RDXN,X_S,DX_WGT,MF_FRAC,PARTICLE_INSERT_CLOCK
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: RHO_0
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MASS_FRACTION,MASS_FLUX,TAU,ADJUST_BURN_RATE,DDSUM,SMALLEST_CELL_SIZE
   INTEGER,  ALLOCATABLE, DIMENSION(:) :: RAMP_INDEX
   INTEGER, DIMENSION(3) :: RGB
   REAL(EB) :: TRANSPARENCY
   REAL(EB), DIMENSION(2) :: VEL_T,EMBER_GENERATION_HEIGHT=-1._EB
   INTEGER, DIMENSION(2) :: LEAK_PATH,DUCT_PATH
   INTEGER :: THERMAL_BC_INDEX,NPPC,SPECIES_BC_INDEX,VELOCITY_BC_INDEX,SURF_TYPE,N_CELLS_INI,N_CELLS_MAX=0, &
              PART_INDEX,PROP_INDEX=-1,RAMP_T_I_INDEX=-1, RAMP_T_B_INDEX=0
   INTEGER, DIMENSION(10) :: INIT_INDICES=0
   INTEGER :: PYROLYSIS_MODEL
   INTEGER :: N_LAYERS,N_MATL,SUBSTEP_POWER=2,N_SPEC=0,N_LPC=0
   INTEGER :: N_ONE_D_STORAGE_REALS,N_ONE_D_STORAGE_INTEGERS,N_ONE_D_STORAGE_LOGICALS
   INTEGER, DIMENSION(30) :: ONE_D_REALS_ARRAY_SIZE=0,ONE_D_INTEGERS_ARRAY_SIZE=0,ONE_D_LOGICALS_ARRAY_SIZE=0
   INTEGER, ALLOCATABLE, DIMENSION(:) :: N_LAYER_CELLS,LAYER_INDEX,MATL_INDEX,MATL_PART_INDEX
   INTEGER, DIMENSION(MAX_LAYERS,MAX_MATERIALS) :: LAYER_MATL_INDEX
   INTEGER, DIMENSION(MAX_LAYERS) :: N_LAYER_MATL,N_LAYER_CELLS_MAX
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MIN_DIFFUSIVITY
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: LAYER_THICKNESS,INTERNAL_HEAT_SOURCE
   REAL(EB), DIMENSION(MAX_LAYERS) :: LAYER_DENSITY,TMP_INNER,STRETCH_FACTOR,&
                                      MOISTURE_FRACTION,SURFACE_VOLUME_RATIO,PACKING_RATIO,KAPPA_S=-1._EB,RENODE_DELTA_T
   REAL(EB), DIMENSION(MAX_NUMBER_FSK_POINTS) :: FSK_K, FSK_W, FSK_A
   REAL(EB), DIMENSION(MAX_LAYERS,MAX_MATERIALS) :: DENSITY_ADJUST_FACTOR=1._EB,RHO_S
   INTEGER :: NUMBER_FSK_POINTS
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: MATL_NAME
   CHARACTER(LABEL_LENGTH), DIMENSION(MAX_LAYERS,MAX_MATERIALS) :: LAYER_MATL_NAME
   REAL(EB), DIMENSION(MAX_LAYERS,MAX_MATERIALS) :: LAYER_MATL_FRAC
   LOGICAL :: BURN_AWAY,ADIABATIC,INTERNAL_RADIATION,USER_DEFINED=.TRUE., &
              FREE_SLIP=.FALSE.,NO_SLIP=.FALSE.,SPECIFIED_NORMAL_VELOCITY=.FALSE.,SPECIFIED_TANGENTIAL_VELOCITY=.FALSE., &
              SPECIFIED_NORMAL_GRADIENT=.FALSE.,CONVERT_VOLUME_TO_MASS=.FALSE.,SPECIFIED_HEAT_SOURCE=.FALSE.,&
              IMPERMEABLE=.FALSE.,BOUNDARY_FUEL_MODEL=.FALSE., &
              HT3D=.FALSE., MT1D=.FALSE.,SET_H=.FALSE.
   LOGICAL :: INCLUDE_BOUNDARY_COORD_TYPE=.TRUE.     !< This surface requires basic coordinate information
   LOGICAL :: INCLUDE_BOUNDARY_PROPS_TYPE=.TRUE.  !< This surface requires surface variables for heat and mass transfer
   LOGICAL :: INCLUDE_BOUNDARY_ONE_D_TYPE=.TRUE.     !< This surface requires in-depth 1-D conduction/reaction arrays
   LOGICAL :: INCLUDE_BOUNDARY_RADIA_TYPE=.TRUE.     !< This surface requires angular-specific radiation intensities
   LOGICAL :: HORIZONTAL=.FALSE.                     !< Indicates if a cylinder is horizontally oriented
   INTEGER :: N_WALL_STORAGE_REALS=0,N_WALL_STORAGE_INTEGERS=0,N_WALL_STORAGE_LOGICALS=0
   INTEGER :: N_CFACE_STORAGE_REALS=0,N_CFACE_STORAGE_INTEGERS=0,N_CFACE_STORAGE_LOGICALS=0
   INTEGER :: GEOMETRY,BACKING,PROFILE,HEAT_TRANSFER_MODEL=0,NEAR_WALL_TURB_MODEL=5
   CHARACTER(LABEL_LENGTH) :: PART_ID,RAMP_Q,RAMP_V,RAMP_T,RAMP_EF,RAMP_PART,RAMP_V_X,RAMP_V_Y,RAMP_V_Z,RAMP_T_B,RAMP_T_I
   CHARACTER(LABEL_LENGTH), ALLOCATABLE, DIMENSION(:) :: RAMP_MF
   CHARACTER(LABEL_LENGTH) :: ID,TEXTURE_MAP,LEAK_PATH_ID(2)
   CHARACTER(MESSAGE_LENGTH) :: FYI='null'
   CHARACTER(LABEL_LENGTH), DIMENSION(10) :: INIT_IDS='null'

   ! 1D mass transfer

   REAL(EB), DIMENSION(MAX_LAYERS,MAX_SPECIES) :: LAYER_SPEC_FRAC
   REAL(EB), DIMENSION(MAX_LAYERS) :: LAYER_POROSITY
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: PHIRHOZ_0

   ! Level Set Firespread

   LOGICAL :: VEG_LSET_SPREAD,VEG_LSET_TAN2
   REAL(EB) :: VEG_LSET_IGNITE_T,VEG_LSET_ROS_HEAD,VEG_LSET_ROS_00,VEG_LSET_QCON,VEG_LSET_ROS_FLANK,VEG_LSET_ROS_BACK, &
               VEG_LSET_WIND_EXP,VEG_LSET_SIGMA,VEG_LSET_HT,VEG_LSET_BETA,&
               VEG_LSET_M1,VEG_LSET_M10,VEG_LSET_M100,VEG_LSET_MLW,VEG_LSET_MLH,VEG_LSET_SURF_LOAD,VEG_LSET_FIREBASE_TIME, &
               VEG_LSET_CHAR_FRACTION
   INTEGER :: VEG_LSET_FUEL_INDEX

   TYPE(STORAGE_TYPE) :: WALL_STORAGE,CFACE_STORAGE

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
   LOGICAL :: HT3D=.FALSE.,HT3D_RESTART=.FALSE.
   LOGICAL :: PYRO3D_LIQUID=.FALSE.
   INTEGER :: MATL_SURF_INDEX=-1
   INTEGER :: PYRO3D_IOR=0
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: RHO

END TYPE OBSTRUCTION_TYPE


TYPE TRIBIN_TYPE
   REAL(EB):: X1_LOW, X1_HIGH
   INTEGER :: NTL
   INTEGER, ALLOCATABLE, DIMENSION(:) :: TRI_LIST
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
   CHARACTER(LABEL_LENGTH) :: GEOC_FILENAME='null',TEXTURE_MAPPING
   LOGICAL :: COMPONENT_ONLY,IS_DYNAMIC=.TRUE.,HAVE_SURF,HAVE_MATL,AUTO_TEXTURE,HIDDEN,REMOVEABLE,SHOW_BNDF=.TRUE., &
              READ_BINARY=.FALSE.,SNAP_TO_GRID=.FALSE.,IS_TERRAIN=.FALSE.,EXPAND_CUTCELLS=.FALSE.
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
   INTEGER :: IE=0
   INTEGER :: NVERT, NEDGE, NEDGE1, STATUS         ! Local Vertices, cut-edges and status of this Cartesian edge.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  XYZVERT  ! Locations of vertices.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::   CEELEM  ! Cut-Edge connectivities.
   INTEGER,  DIMENSION(MAX_DIM+2)                  ::      IJK  ! [ i j k X2AXIS cetype]
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::   INDSEG  ! [ntr tr1 tr2 ibod]
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::      DXX  ! [DXX(1,JEC) DXX(2,JEC)]
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:)         ::FACE_LIST  ! [1:2, -2:2, JEC] Cut-face connected to edge.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::   DUIDXJ, MU_DUIDXJ ! Unstructured VelGrad components.
   INTEGER,  ALLOCATABLE, DIMENSION(:)             :: NOD_PERM  ! Permutation array for INSERT_FACE_VERT.
END TYPE IBM_CUTEDGE_TYPE

! IBM_EDGE type, used for computation of wall model turbulent viscosity, shear stress, vorticity.
TYPE IBM_EDGE_TYPE
   INTEGER,  DIMENSION(MAX_DIM+1)                  ::     IJK  ! [ i j k X1AXIS]
   INTEGER :: IE=0
   INTEGER, ALLOCATABLE, DIMENSION(:,:)            ::FACE_LIST  ! [1:3, -2:2] Reg/Cut face connected to edge.
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
   REAL(EB), ALLOCATABLE, DIMENSION(:)        :: XB_IB,DUIDXJ,MU_DUIDXJ
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
   INTEGER :: IWC=0,PRES_ZONE=-1
   INTEGER :: NVERT=0, NSVERT=0, NFACE=0, NSFACE=0, STATUS !Local Vertices, cut-faces and status of this Cartesian face.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           :: XYZVERT  ! Locations of vertices.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  CFELEM  ! Cut-faces connectivities.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  CEDGES  ! Cut-Edges. Points to EDGE_LIST.
   INTEGER,  DIMENSION(MAX_DIM+1)                  ::     IJK  ! [ i j k X1AXIS]
   REAL(EB), ALLOCATABLE, DIMENSION(:)             ::    AREA  ! Cut-faces areas.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  XYZCEN  ! Cut-faces centroid locations.
   LOGICAL,  ALLOCATABLE, DIMENSION(:)             ::  SHARED
   INTEGER,  ALLOCATABLE, DIMENSION(:)             ::  LINK_LEV ! Level in local Face Linking Hierarchy.
   !Integrals to be used in cut-cell volume and centroid computations.
   REAL(EB), ALLOCATABLE, DIMENSION(:)             ::  INXAREA, INXSQAREA, JNYSQAREA, KNZSQAREA
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  BODTRI
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  UNKH, UNKZ
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  XCENLOW, XCENHIGH
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  RHO
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  ZZ_FACE, RHO_D
   REAL(EB), ALLOCATABLE, DIMENSION(:)             ::  TMP_FACE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)           ::  RHO_D_DZDN, H_RHO_D_DZDN
   REAL(EB), ALLOCATABLE, DIMENSION(:)             ::  VEL, VELS, FN, FN_B, VELINT
   INTEGER,  ALLOCATABLE, DIMENSION(:,:,:)         ::  JDH
   REAL(EB) :: VELINT_CRF=0._EB,FV=0._EB,FV_B=0._EB,ALPHA_CF=1._EB,VEL_CF=0._EB,VEL_CRT=0._EB
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)           ::  EDGE_LIST ! [CE_TYPE IEC JEC] or [RG_TYPE SIDE LOHI_AXIS]
   INTEGER,  ALLOCATABLE, DIMENSION(:,:,:)         ::  CELL_LIST ! [RC_TYPE I J K  ]

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

   REAL(EB), ALLOCATABLE, DIMENSION(:,:)      :: RHOPVN
   INTEGER :: NOMICF(2)=0  ! [NOM icf]
   INTEGER,  ALLOCATABLE, DIMENSION(:)        :: CFACE_INDEX, SURF_INDEX, UNKF
END TYPE IBM_CUTFACE_TYPE


TYPE RAD_CFACE_TYPE
   INTEGER :: N_ASSIGNED_CFACES_RADI=0
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ASSIGNED_CFACES_RADI
END TYPE RAD_CFACE_TYPE


! Note: If you change the number of scalar variables in CFACE_TYPE, adjust the numbers below

INTEGER, PARAMETER :: N_CFACE_SCALAR_REALS=13
INTEGER, PARAMETER :: N_CFACE_SCALAR_INTEGERS=12
INTEGER, PARAMETER :: N_CFACE_SCALAR_LOGICALS=0

TYPE CFACE_TYPE
   INTEGER :: CFACE_INDEX=0              !< Index of itself -- used to determine if the CFACE cell has been assigned
   INTEGER :: BC_INDEX=0                 !< Derived type carrying coordinate variables
   INTEGER :: OD_INDEX=0                 !< Derived type carrying 1-D solid info
   INTEGER :: BP_INDEX=0                 !< Derived type carrying most of the surface boundary conditions
   INTEGER :: BR_INDEX=0                 !< Derived type carrying angular-specific radiation intensities
   INTEGER :: SURF_INDEX=0
   INTEGER :: VENT_INDEX=0
   INTEGER :: BACK_MESH=0
   INTEGER :: BACK_INDEX=0
   INTEGER :: BOUNDARY_TYPE=0
   INTEGER :: CUT_FACE_IND1=-11
   INTEGER :: CUT_FACE_IND2=-11
   REAL(EB) :: AREA=0._EB
   REAL(EB) :: NVEC(3)=0._EB
   REAL(EB) :: VEL_ERR_NEW=0._EB
   REAL(EB) :: V_DEP=0._EB
   REAL(EB) :: Q_LEAK=0._EB
   REAL(EB) :: DUNDT=0._EB
   REAL(EB) :: RSUM_G=0._EB                     !< \f$ R_0 \sum_\alpha Z_\alpha/W_\alpha \f$ in first gas phase cell
   REAL(EB) :: TMP_G                            !< Temperature (K) in adjacent gas phase cell
   REAL(EB) :: RHO_G                            !< Gas density (kg/m3) in adjacent gas phase cell
   REAL(EB) :: MU_G=0.1_EB                      !< Viscosity, \f$ \mu \f$, in adjacent gas phase cell
   REAL(EB) :: PRES_BXN
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_G  !< (1:N_TRACKED_SPECIES) Species mixture mass fraction in gas
END TYPE CFACE_TYPE

! Cartesian Cells Cut-Cells data structure:

INTEGER, PARAMETER :: IBM_MAXVERTS_CELL   =3072
INTEGER, PARAMETER :: IBM_NPARAM_CCFACE   =   6 ! [face_type side iaxis cei icf to_master]

TYPE IBM_CUTCELL_TYPE
   INTEGER :: NCELL, NFACE_CELL
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)                     ::    CCELEM ! Cut-cells faces connectivities in FACE_LIST.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)                     :: FACE_LIST ! List of faces, cut-faces.
   INTEGER,  ALLOCATABLE, DIMENSION(:,:)                     ::  IJK_LINK ! Cell/cut-cell each cut-cell is linked to.
   INTEGER,  ALLOCATABLE, DIMENSION(:)                       ::  LINK_LEV ! Level in local Linking Hierarchy tree.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::    VOLUME ! Cut-cell volumes.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     ::    XYZCEN ! Cut-cell centroid locaitons.
   INTEGER,  DIMENSION(MAX_DIM)                              ::       IJK ! [ i j k ]
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: RHO, RHOS ! Cut cells densities.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       ::  RSUM,TMP ! Cut cells temperatures.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: D,DS,DVOL ! Cut cell thermodynamic divg.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: Q,QR,D_SOURCE ! Q,Thermo divg reaction component.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: CHI_R,MIX_TIME ! Cut-cell combustion
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     ::    Q_REAC      ! variables.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     :: REAC_SOURCE_TERM   !
   REAL(EB), ALLOCATABLE, DIMENSION(:,:)                     ::   ZZ, ZZS, M_DOT_PPP ! Cut cells species mass
                                                                                ! fractions and rho*D_z,reaction source.
   INTEGER,  ALLOCATABLE, DIMENSION(:)                       :: UNKH,UNKZ ! Unknown number for pressure H,
                                                                          ! and scalars.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: KRES,H,HS ! Kinetic Energy, Pressure H containers.
   REAL(EB), ALLOCATABLE, DIMENSION(:)                       :: RTRM,R_H_G,RHO_0,WVEL

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
   INTEGER :: N_NOMICC=0
   INTEGER,  ALLOCATABLE, DIMENSION(:,:) :: NOMICC
   REAL(EB):: DIVVOL_PR=0._EB, DIVVOL_BC=0._EB, DDDTVOL=0._EB
END TYPE IBM_CUTCELL_TYPE


TYPE IBM_REGFACE_TYPE
   INTEGER:: PRES_ZONE=-1
   INTEGER,  DIMENSION(MAX_DIM)                                    ::       IJK
   INTEGER,  DIMENSION(1:2,1:2)                                    ::        JD
END TYPE IBM_REGFACE_TYPE

TYPE IBM_REGFACEZ_TYPE
   LOGICAL :: DO_LO_IND=.FALSE.,DO_HI_IND=.FALSE.
   INTEGER :: IWC=0
   REAL(EB):: FN_H_S=0._EB
   INTEGER,  DIMENSION(MAX_DIM)                                    ::       IJK
   INTEGER,  DIMENSION(1:2,1:2)                                    ::        JD
   REAL(EB), DIMENSION(MAX_SPECIES)                                ::   RHO_D=0._EB,RHOZZ_U=0._EB,FN_ZZ=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES)                                ::   RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES)                                :: H_RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(-1:0)                                       ::    RHOPVN=0._EB
END TYPE IBM_REGFACEZ_TYPE

TYPE IBM_RCFACE_TYPE
   LOGICAL:: SHAREDH=.FALSE.
   INTEGER:: PRES_ZONE=-1
   INTEGER,  DIMENSION(MAX_DIM+1)                                  ::       IJK ! [ I J K x1axis]
   INTEGER,  DIMENSION(LOW_IND:HIGH_IND)                           ::      UNKH
   REAL(EB), DIMENSION(MAX_DIM,LOW_IND:HIGH_IND)                   ::      XCEN
   INTEGER,  DIMENSION(1:2,1:2)                                    ::       JDH
END TYPE IBM_RCFACE_TYPE

TYPE IBM_RCFACE_LST_TYPE
   LOGICAL :: SHAREDZ=.FALSE.
   INTEGER :: IWC=0, UNKF=0
   REAL(EB):: TMP_FACE=0._EB
   INTEGER,  DIMENSION(MAX_DIM+1)                                  ::       IJK ! [ I J K x1axis]
   INTEGER,  DIMENSION(LOW_IND:HIGH_IND)                           ::      UNKZ
   REAL(EB), DIMENSION(MAX_DIM,LOW_IND:HIGH_IND)                   ::      XCEN
   INTEGER,  DIMENSION(MAX_DIM+1,LOW_IND:HIGH_IND)                 :: CELL_LIST ! [RC_TYPE I J K ]
   REAL(EB), DIMENSION(MAX_SPECIES)                                :: ZZ_FACE=0._EB,RHO_D=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES)                                :: RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(MAX_SPECIES)                                :: H_RHO_D_DZDN=0._EB
   REAL(EB), DIMENSION(-1:0)                                       ::    RHOPVN=0._EB
END TYPE IBM_RCFACE_LST_TYPE

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
   INTEGER :: NUMBER_DATA_POINTS,NUMBER_INTERPOLATION_POINTS,DEVC_INDEX=-1,CTRL_INDEX=-1,RESERVED_RAMP_INDEX=0
   CHARACTER(LABEL_LENGTH) :: DEVC_ID='null',CTRL_ID='null'
   LOGICAL :: DEP_VAR_UNITS_CONVERTED=.FALSE.
END TYPE RAMPS_TYPE

TYPE (RAMPS_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: RAMPS

TYPE RESERVED_RAMPS_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: INDEPENDENT_DATA,DEPENDENT_DATA
   INTEGER :: NUMBER_DATA_POINTS
   CHARACTER(LABEL_LENGTH) :: ID
END TYPE RESERVED_RAMPS_TYPE

INTEGER :: N_RESERVED_RAMPS=0
TYPE (RESERVED_RAMPS_TYPE), DIMENSION(10), TARGET :: RESERVED_RAMPS

TYPE SLICE_TYPE
   INTEGER :: I1,I2,J1,J2,K1,K2,GEOM_INDEX=-1,TRNF_INDEX=-1,INDEX,INDEX2=0,Z_INDEX=-999,Y_INDEX=-999,MATL_INDEX=-999,&
              PART_INDEX=0,VELO_INDEX=0,PROP_INDEX=0,REAC_INDEX=0,SLCF_INDEX
   REAL(FB), DIMENSION(2) :: MINMAX
   REAL(FB) :: RLE_MIN, RLE_MAX
   REAL(EB):: AGL_SLICE
   LOGICAL :: TERRAIN_SLICE=.FALSE.,CELL_CENTERED=.FALSE.,FACE_CENTERED=.FALSE.,RLE=.FALSE.,DEBUG=.FALSE.
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
   LOGICAL :: CELL_CENTERED=.FALSE.
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
   REAL(EB) :: CHI_R            !< Radiative fraction of HRRPUV
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: PARTICLE_INSERT_CLOCK  !< Time of last particle insertion (s)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MASS_FRACTION          !< Mass fraction of gas components
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: VOLUME_ADJUST          !< Multiplicative factor to account for FDS grid snap
   INTEGER  :: PART_INDEX=0     !< Particle class index of inserted particles
   INTEGER  :: N_PARTICLES      !< Number of particles to insert
   INTEGER  :: DEVC_INDEX=0     !< Index of the device that uses this INITIALIZATION variable
   INTEGER  :: CTRL_INDEX=0     !< Index of the controller that uses this INITIALIZATION variable
   INTEGER  :: N_PARTICLES_PER_CELL=0 !< Number of particles to insert in each cell
   INTEGER  :: PATH_RAMP_INDEX(3)=0   !< Ramp index of a particle path
   INTEGER  :: RAMP_Q_INDEX=0         !< Ramp index for HRRPUV
   INTEGER  :: RAMP_PART_INDEX=0         !< Ramp index for MASS_PER_TIME or MASS_PER_VOLUME
   LOGICAL :: ADJUST_DENSITY=.FALSE.
   LOGICAL :: ADJUST_TEMPERATURE=.FALSE.
   LOGICAL :: SINGLE_INSERTION=.TRUE.
   LOGICAL :: CELL_CENTERED=.FALSE.
   LOGICAL :: UNIFORM=.FALSE.
   LOGICAL :: INVOKED_BY_SURF=.FALSE.  ! Invoked by a SURF line for repeated insertion
   LOGICAL :: DRY=.FALSE.              ! Indicates if MASS_PER_VOLUME refers to dry material only
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: ALREADY_INSERTED
   CHARACTER(LABEL_LENGTH) :: SHAPE
   CHARACTER(LABEL_LENGTH) :: DEVC_ID
   CHARACTER(LABEL_LENGTH) :: CTRL_ID
   CHARACTER(LABEL_LENGTH) :: NODE_ID
   CHARACTER(LABEL_LENGTH) :: ID
   CHARACTER(MESSAGE_LENGTH) :: BULK_DENSITY_FILE
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
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DISCHARGE_COEFFICIENT    !< Array of leak disharge coefficients
   INTEGER :: N_DUCTNODES                                          !< Number of duct nodes in the ZONE
   INTEGER :: MESH_INDEX=0                                         !< Index of the MESH where the ZONE is located
   INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_INDEX                !< Array of NODE indices connected to the ZONE
   CHARACTER(LABEL_LENGTH) :: ID='null'                            !< Identifier
END TYPE P_ZONE_TYPE

TYPE (P_ZONE_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: P_ZONE


!> \brief Parameters associated with a ZONE within a MESH, used in LOCMAT_SOLVER unstructured pressure solver

TYPE ZONE_MESH_TYPE
#ifdef WITH_MKL
   TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: PT_H(:)  !< Internal solver memory pointer
#else
   INTEGER, ALLOCATABLE :: PT_H(:)
#endif /* WITH_MKL */
   INTEGER :: NUNKH=0                                 !< Number of unknowns in pressure solution for a given ZONE_MESH
   INTEGER :: NCVLH=0                                 !< Number of pressure control volumes for a given ZONE_MESH
   INTEGER :: ICVL=0                                  !< Control volume counter for parent ZONE
   INTEGER :: IROW=0                                  !< Parent ZONE matrix row index
   INTEGER :: NUNKH_CART=0                            !< Number of unknowns in Cartesian cells of ZONE_MESH
   INTEGER :: NCVLH_CART=0                            !< Number of pressure CVs in Cartesian cells of ZONE_MESH
   INTEGER :: MTYPE=0                                 !< Matrix type (symmetric indefinite, or symm positive definite)
   INTEGER :: CONNECTED_ZONE_PARENT=0                 !< Index of first zone in a connected zone list
   LOGICAL :: USE_FFT=.TRUE.                          !< Flag for use of FFT solver
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: MESH_IJK   !< I,J,K positions of cell with unknown row IROW (1:3,1:NUNKH)
   REAL(EB),ALLOCATABLE, DIMENSION(:)   :: A_H        !< Matrix coefficients for ZONE_MESH, up triang part, CSR format.
   INTEGER ,ALLOCATABLE, DIMENSION(:)   :: IA_H,JA_H  !< Matrix indexes for ZONE_MESH, up triang part, CSR format.
   REAL(EB),ALLOCATABLE, DIMENSION(:)   :: F_H,X_H    !< RHS and Solution containers for the ZONE_MESH.
END TYPE ZONE_MESH_TYPE


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

!> \brief Parameters associated with HVAC TYPE_ID FILTER

TYPE FILTER_TYPE
   INTEGER :: RAMP_INDEX = -1                             !< Table of loading vs. loss coefficient
   CHARACTER(LABEL_LENGTH) :: ID                          !< Name of the filter
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: EFFICIENCY      !< Array giving species filtering efficiency
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: MULTIPLIER      !< Array containing species multiplier used to compute loading
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: INITIAL_LOADING !< Array of containing initial filter loading (kg) for each species
   REAL(EB) :: CLEAN_LOSS                                 !< Loss coefficient of a filter when the loading is 0 kg.
   REAL(EB) :: LOADING_LOSS                               !< Multiplier applied to loading to determine filter loss
   REAL(EB) :: AREA                                       !< Cross-sectional area of filter
END TYPE FILTER_TYPE

TYPE (FILTER_TYPE), DIMENSION(:), ALLOCATABLE,  TARGET :: FILTER

!> \brief Parameters associated with HVAC TYPE_ID NODE

TYPE DUCTNODE_TYPE
   INTEGER :: FILTER_INDEX=-1                              !< Index of a filter contained in the node
   INTEGER :: N_DUCTS                                      !< Number of ducts attached to the node
   INTEGER :: VENT_INDEX = -1                              !< Index of a VENT the node is attached to
   INTEGER :: ZONE_INDEX=-1                                !< Pressure zone containing the node
   INTEGER :: DUCTRUN                                      !< Ductrun node belongs to
   INTEGER :: DUCTRUN_INDEX=-1                             !< Index in ductrun node belongs to
   INTEGER :: DUCTRUN_M_INDEX=-1                           !< Index of node in ductrun solution matrix
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DUCT_INDEX        !< List of ducts attached to the node
   CHARACTER(LABEL_LENGTH) :: ID                           !< Name of the node
   CHARACTER(LABEL_LENGTH) :: VENT_ID='null'               !< Name of a VENT the node is attached to
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: LOSS_ARRAY     !< (i,j) Flow loss for flow from duct i to duct j
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: FILTER_LOADING !< (loading for each species, old/new/guess)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ               !< Current species mass fractions at the node
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DIR              !< (i) Duct i lists the node first (-1) or second (+1)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_V             !< Current species mass fractions at an attached CENT
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_OLD           !< Prior iteration species at the node
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ0              !< INIT node mass fraction
   REAL(EB) :: LOSS                                        !< Current node flow loss
   REAL(EB) :: TMP=273.15_EB                               !< Current node temperature (K)
   REAL(EB) :: TMP_OLD                                     !< Prior timestep node temperature (K)
   REAL(EB) :: TMP_V                                       !< Temperature of VENT connected to node (K)
   REAL(EB) :: TMP0                                        !< INIT derived node temperature (K)
   REAL(EB) :: RHO                                         !< Current node density (kg/m3)
   REAL(EB) :: RHO_OLD                                     !< Prior timestep node density (kg/m3)
   REAL(EB) :: RHO_V                                       !< Node density of VENT connected to node (kg/m3)
   REAL(EB) :: RSUM                                        !< Specfiic gas constant of node (J/kg/K)
   REAL(EB) :: RSUM_OLD                                    !< Prior timestep specfic gas constant (J/kg/K)
   REAL(EB) :: RSUM_V                                      !< Specific gas constant of VENT connected to node (J/kg/K)
   REAL(EB) :: CP                                          !< Specific heat of node (J/kg/K)
   REAL(EB) :: CP_OLD                                      !< Prior timestep specific heat of node (J/kg/K)
   REAL(EB) :: CP_V                                        !< Specific heat of VENT connected to node (J/kg/K)
   REAL(EB) :: XYZ(3)=(/0._EB,0._EB,0._EB/)                !< Node (x,y,z) coordinate (m)
   REAL(EB) :: FILTER_LOSS                                 !< Current loss for a filter in the node
   REAL(EB) :: P                                           !< Pressure at the node (Pa)
   REAL(EB) :: P_OLD                                       !< Prior timestep pressure at the node (Pa)
   LOGICAL :: UPDATED                                      !< Node has been updated in the current timestep
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_SUM           !< Holds vent total species flow for HVAC_MASS_TRANSPORT
   REAL(EB) :: E_SUM                                       !< Holds vent total energy flow for HVAC_MASS_TRANSPORT
   LOGICAL :: READ_IN                                      !< Node defintion is explicit in the input file
   LOGICAL :: AMBIENT = .FALSE.                            !< Node is connected to the ambient
   LOGICAL :: LEAKAGE=.FALSE.                              !< Node is being used for leakage
   LOGICAL :: VENT=.FALSE.                                 !< Node has an attached vent
   LOGICAL :: HMT_FILTER=.FALSE.                           !< Filter is in mass transport ductrun
   LOGICAL :: TRANSPORT_PARTICLES=.FALSE.                  !< Particles will be transported through the vent attached to the node
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: IN_MESH           !< (i) Flag indicating node is present in mesh i
END TYPE DUCTNODE_TYPE

TYPE (DUCTNODE_TYPE), DIMENSION(:), ALLOCATABLE,  TARGET :: DUCTNODE

!> \brief Parameters associated with aggregating node boundary conditions over the area of a VENT

TYPE NODE_BC_TYPE
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_V !< Stores aggregated species mass fraction
   REAL(EB) :: TMP_V                           !< Stores aggreated temperature
   REAL(EB) :: RHO_V                           !< Stores aggreated density
   REAL(EB) :: RSUM_V                          !< Stores aggreated gas constant
   REAL(EB) :: CP_V                            !< Stores aggreated specific heat
   REAL(EB) :: P                               !< Stores aggreated pressure
END TYPE NODE_BC_TYPE

!> \brief Parameters associated with HVAC TYPE_ID DUCT

TYPE DUCT_TYPE
   INTEGER :: NODE_INDEX(2)=-1                            !< Nodes connected to ends of duct
   INTEGER :: DEVC_INDEX=-1                               !< Index for DEVC controlling fan, damper, or aircoil
   INTEGER :: CTRL_INDEX=-1                               !< Index for CTRL controlling fan, damper, or aircoil
   INTEGER :: FAN_INDEX=-1                                !< Index of fan present in duct
   INTEGER :: AIRCOIL_INDEX=-1                            !< Type of aircoil present in fuct
   INTEGER :: SURF_INDEX=-1                               !< Surface type for duct wall
   INTEGER :: RAMP_INDEX=0                                !< Index of RAMP used for MASS_FLOW or VOLUME_FLOW
   INTEGER :: RAMP_LOSS_INDEX=0                           !< Ramp for damper flow loss
   INTEGER :: N_CELLS=0                                   !< Number of duct segments for HVAC mass transport
   INTEGER :: N_WAYPOINTS=0                               !< Number of waypoints for a duct
   INTEGER :: LPC_INDEX=-1                                !< Particle class for visualization
   INTEGER :: QFAN_INDEX=-1                               !< Index of duct fan in ductrun fan array
   INTEGER :: DUCTRUN=-1                                  !< Ductrun duct belongs to
   INTEGER :: DUCTRUN_INDEX=-1                            !< Index in ductrun duct belongs to
   INTEGER :: DUCTRUN_M_INDEX=-1                          !< Index of duct in ductrun solution matrix
   INTEGER :: N_HT_SEGMENTS=0                             !< Number of heat trasnfer segments
   REAL(EB) :: AREA                                       !< Current duct cross sectional area (m2)
   REAL(EB) :: AREA_OLD                                   !< Prior timestep duct cross sectional area (m2)
   REAL(EB) :: AREA_INITIAL                               !< Input duct cross sectional area (m2)
   REAL(EB) :: COIL_Q=0._EB                               !< Current heat rate from an aircoil (W)
   REAL(EB) :: HT_Q=0._EB                                 !< Current duct wall heat transfer rate (W)
   REAL(EB) :: CP_D                                       !< Upstream specific heat (J/kg/K)
   REAL(EB) :: CP_D_OLD                                   !< Prior timestep upstream specific heat (J/kg/K)
   REAL(EB) :: DIAMETER                                   !< Duct diamater (m)
   REAL(EB) :: PERIMETER                                  !< Duct perimeter (m)
   REAL(EB) :: DP_FAN(2)=0._EB                            !< Prior and current fan pressure (Pa)
   REAL(EB) :: DX=-1._EB                                  !< Duct segment length (m)
   REAL(EB) :: LENGTH                                     !< Duct length (m)
   REAL(EB) :: MASS_FLOW_INITIAL=1.E6_EB                  !< Input duct mass flow (kg/m3)
   REAL(EB) :: RHO_D                                      !< Upstream density (kg/m3)
   REAL(EB) :: RHO_D_OLD                                  !< Prior timestep upstream density (kg/m3)
   REAL(EB) :: ROUGHNESS                                  !< Wall roughness (m)
   REAL(EB) :: RSUM_D=0._EB                               !< Upstream specific gas constant (J/kg/K)
   REAL(EB) :: RSUM_D_OLD                                 !< Prior timestep upstream specific gas constant (J/kg/K)
   REAL(EB) :: TAU=-1._EB                                 !< TANH or t2 ramp for flow
   REAL(EB) :: TMP_D=273.15_EB                            !< Upstream duct temperature (K)
   REAL(EB) :: TOTAL_LOSS=0._EB                           !< Current flow loss in duct
   REAL(EB) :: VEL(4)=0._EB                               !< Velocity in duct (old,new,guess,previous) (m/s)
   REAL(EB) :: VOLUME_FLOW=1.E6_EB                        !< Current duct volume flow (m3/s)
   REAL(EB) :: VOLUME_FLOW_INITIAL=1.E6_EB                !< Input duct volume flow (m3/s)
   REAL(EB) :: LOSS(2)=0._EB                              !< Upstream specific heat (J/kg/K)
   REAL(EB) :: LEAK_PRESSURE_EXPONENT=0.5_EB              !< Pressure exponent in leakage expression
   REAL(EB) :: LEAK_REFERENCE_PRESSURE=4._EB              !< Reference pressure in leakage expression (Pa)
   REAL(EB) :: DISCHARGE_COEFFICIENT=1._EB                !< Discharge coefficient, C, in leakage expression
   INTEGER, ALLOCATABLE, DIMENSION(:) :: HT_INDEX         !< Index to a duct heat transfer segment
   INTEGER, ALLOCATABLE, DIMENSION(:) :: SEGMENT_INDEX    !< Mapping of heat transfer segements to mass flow segments
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CP_C            !< Current specific heat in each duct segment (J/kg/K)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CP_C_OLD        !< Prior timestep specific heat in each duct segment (J/kg/K)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_C           !< Current density in each duct segment (kg/m3)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_C_OLD       !< Prior timestep density in each duct segment (kg/m3)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: TMP_C           !< Current temperature in each duct segment (K)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: TMP_C_OLD       !< Prior timestep temperature in each duct segment (K)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ              !< Current species mass fractions in duct
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZZ_OLD          !< Prior timestep species mass fractions in duct
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: ZZ_C          !< Current species mass fractions in each duct segment
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: ZZ_C_OLD      !< Prior timestep species mass fractions in each duct segment
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: WAYPOINTS_XYZ !< Waypoints for a duct (m)
   REAL(EB), ALLOCATABLE, DIMENSION(:)   :: WAYPOINTS_L   !< Length of duct between waypoints (m)
   LOGICAL :: DAMPER = .FALSE.                            !< Duct contains a damper
   LOGICAL :: DAMPER_OPEN = .TRUE.                        !< Duct damper is open
   LOGICAL :: FAN_OPERATING=.TRUE.                        !< Duct fan is operating
   LOGICAL :: COIL_OPERATING=.TRUE.                       !< Duct aircoil is operating
   LOGICAL :: FIXED=.FALSE.                               !< Duct has flow defined
   LOGICAL :: UPDATED                                     !< Duct solution has been updated in the current timestep
   LOGICAL :: REVERSE=.FALSE.                             !< Fan direction is opposite that implied by the node order
   LOGICAL :: LEAK_ENTHALPY=.FALSE.                       !< Conserve enthalpy for leakage flow
   LOGICAL :: LEAKAGE=.FALSE.                             !< Duct is being used for bulk leakage
   LOGICAL :: LOCALIZED_LEAKAGE=.FALSE.                   !< Duct is being used for localized leakage
   LOGICAL :: DUCT_HT=.FALSE.                             !< Duct heat transfer
   LOGICAL :: ROUND=.TRUE.                                !< Duct is a round duct
   CHARACTER(LABEL_LENGTH) :: ID                          !< Name of duct
   REAL(EB) :: FAN_ON_TIME = 1.E10_EB                     !< Time fan in duct is turned on
   REAL(EB) :: COIL_ON_TIME = 1.E10_EB                    !< Time aircoil in duct is turned on
END TYPE DUCT_TYPE

TYPE (DUCT_TYPE), DIMENSION(:), ALLOCATABLE,  TARGET :: DUCT

!>\brief Parameters associated with a duct heat transfer segement
TYPE DUCT_HT_SEGMENT_TYPE
   INTEGER:: DUCT_INDEX                       !< Index of duct HT segment is in
   INTEGER:: SEGMENT_INDEX                    !< Index of duct mass transfer segment HT segement is in
   INTEGER:: MESH_INDEX                       !< Mesh location of segment
   INTEGER :: IIG                             !< Gas cell \f$ x \f$ index
   INTEGER :: JJG                             !< Gas cell \f$ y \f$ index
   INTEGER :: KKG                             !< Gas cell \f$ z \f$ index
   INTEGER :: N_SUBSTEPS                      !< Number of sub timesteps for heat transfer solution
   REAL(EB) :: X                              !< \f$ x \f$ coordinate of segment center in gas cell (m)
   REAL(EB) :: Y                              !< \f$ y \f$ coordinate of segment center in gas cell (m)
   REAL(EB) :: Z                              !< \f$ z \f$ coordinate of segment center in gas cell (m)
   REAL(EB), ALLOCATABLE, DIMENSION(:):: TMP  !< Duct wall node temperatures (K)
   REAL(EB):: TMP_I                           !< Inside wall temperature (K)
   REAL(EB):: TMP_O                           !< Outside wall temperature (K)
   REAL(EB):: TMP_G                           !< Gas temperature outside the duct (K)
   REAL(EB):: RHO_G                           !< Density outside the duct (K)
   REAL(EB):: VEL_G                           !< Gas velocity outside of the duct (m/s)
   REAL(EB):: HTC_I                           !< Inside heat trasnfer coefficient (W/m/K)
   REAL(EB):: HTC_O                           !< Outside heat trasnfer coefficient (W/m/K)
   REAL(EB):: Q_CON_I                         !< Inside convectve heat flux (W/m2)
   REAL(EB):: Q_CON_O                         !< Outside convectve heat flux (W/m2)
   REAL(EB):: Q_RAD_I                         !< Net radiation heat flux insde (W/m2)
   REAL(EB):: Q_RAD_O                         !< Incident radiation heat flux outside(W/m2)
   REAL(EB):: KAPPA                           !< Duct absorptivity
END TYPE DUCT_HT_SEGMENT_TYPE

TYPE(DUCT_HT_SEGMENT_TYPE), POINTER, DIMENSION(:) :: DUCT_HT_SEGMENT

!> \brief array for passing HT segment data between meshes

INTEGER, PARAMETER :: N_HT_SEGMENT_REALS = 10
INTEGER, PARAMETER :: N_HT_SEGMENT_INTEGERS = 6

TYPE HT_SEGMENT_TRANSFER_TYPE
   INTEGER:: DUCT_INDEX                       !< Index of duct HT segment is in
   INTEGER:: SEGMENT_INDEX                    !< Index of duct mass transfer segment HT segement is in
   INTEGER:: MESH_INDEX                       !< Mesh location of segment
   INTEGER :: IIG                             !< Gas cell \f$ x \f$ index
   INTEGER :: JJG                             !< Gas cell \f$ y \f$ index
   INTEGER :: KKG                             !< Gas cell \f$ z \f$ index
   REAL(EB) :: X                              !< \f$ x \f$ coordinate of segment center in gas cell (m)
   REAL(EB) :: Y                              !< \f$ y \f$ coordinate of segment center in gas cell (m)
   REAL(EB) :: Z                              !< \f$ z \f$ coordinate of segment center in gas cell (m)
   REAL(EB):: TMP_G                           !< Gas temperature outside the duct (K)
   REAL(EB):: RHO_G                           !< Density outside the duct (K)
   REAL(EB):: VEL_G                           !< Gas velocity outside of the duct (m/s)
   REAL(EB),DIMENSION(:),ALLOCATABLE:: ZZ_G   !< Species mass fractions outside of the duct (kg/kg)
   REAL(EB):: Q_RAD_O                         !< Incident radiation heat flux outside(W/m2)
   REAL(EB):: TMP_O                           !< Outside wall temperature (K)
   REAL(EB):: HTC_O                           !< Outside heat trasnfer coefficient (W/m/K)
   REAL(EB):: Q_CON_O                         !< Outside convectve heat flux (W/m2)
END TYPE HT_SEGMENT_TRANSFER_TYPE

!> \brief Parameters associated with HVAC TYPE_ID FAN

TYPE FAN_TYPE
   INTEGER :: FAN_TYPE                 !< Constant flow = 1, quadratic flow = 2, fan curve = 3
   INTEGER :: RAMP_INDEX               !< Index of RAMP containing fan curve used for FAN_TYPE=3
   INTEGER :: SPIN_INDEX=0             !< Indicates use of t2 or tanh ramp
   REAL(EB) :: VOL_FLOW                !< Constant volume flow (m3/s) used for FAN_TYPE=1
   REAL(EB) :: MAX_FLOW                !< Maximum fan flow (m3/s) used for FAN_TYPE=2
   REAL(EB) :: MAX_PRES                !< Maximum fan pressure (Pa) used for FAN_TYPE=2
   REAL(EB) :: OFF_LOSS=1._EB          !< Flow loss through fan when it is not running
   REAL(EB) :: TAU=0._EB               !< t2 or tanh time constant for spinning up fan speed
   CHARACTER(LABEL_LENGTH) :: ID       !< Name of fan
   CHARACTER(LABEL_LENGTH) :: FAN_RAMP !< Name of RAMP containing fan curve
END TYPE FAN_TYPE

TYPE(FAN_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: FAN

!> \brief Parameters associated with HVAC TYPE_ID AIRCOIL

TYPE AIRCOIL_TYPE
   REAL(EB) :: COOLANT_TEMPERATURE=293.15_EB  !< Temperature of working fluid in aircoil (K)
   REAL(EB) :: COOLANT_SPECIFIC_HEAT=4186._EB !< Specific heat of working fluid in aircoil (J/kg/K)
   REAL(EB) :: COOLANT_MASS_FLOW=-9999._EB    !< Mass flow rate of working fluid in aircoil (kg/s)
   REAL(EB) :: FIXED_Q=-1.E10_EB              !< Fixed heating/cooling rate for aircoil (W)
   REAL(EB) :: TAU=0._EB                      !< TANH or t2 ramp time for aircoil when turned on
   REAL(EB) :: EFFICIENCY                     !< Aircoil heat exchange efficiency
   INTEGER :: RAMP_INDEX=0                    !< Index of RAMP for FIXED_Q
   CHARACTER(LABEL_LENGTH) :: ID              !< Name of aircoil
END TYPE AIRCOIL_TYPE

TYPE(AIRCOIL_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: AIRCOIL

!> \brief Parameters defining independent HVAC networks

TYPE NETWORK_TYPE
   INTEGER :: N_DUCTS                                 !< Number of ducts in network
   INTEGER :: N_DUCTNODES                             !< Number of nodes in netowrk
   INTEGER :: N_MATRIX                                !< Number of element ins solution matrix for network
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DUCT_INDEX   !< List of ducts in network
   INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_INDEX   !< List of nodes in netwwork
   INTEGER, ALLOCATABLE, DIMENSION(:) :: MATRIX_INDEX !< Position of ducts and nodes in solution matrix
END TYPE NETWORK_TYPE

TYPE(NETWORK_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: NETWORK

TYPE DUCTRUN_TYPE
   INTEGER :: N_DUCTS = 0                              !< Number of ducts in ductrun
   INTEGER :: N_DUCTNODES = 0                          !< Number of ductnodes in ductrun
   INTEGER :: N_M_DUCTS = 0                            !< Number of ducts in ductrun solution matrix
   INTEGER :: N_M_DUCTNODES = 0                        !< Number of ductnodes in ductrun solution matrix
   INTEGER :: N_QFANS = 0                              !< Number of fans (FAN%FAN_TYPE/=1) in ductrun
   LOGICAL :: MASS_TRANSPORT = .FALSE.                 !< A mass transport duct is in the ductrun
   REAL(EB) :: DT_CFL = 0._EB                          !< CFL for a ductrun with HVAC mass transport
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: FAN_OPERATING !< If a QFAN is operating
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DUCT_INDEX    !< List of ducts in ductrun
   INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_INDEX    !< List of node in ductrun
   INTEGER, ALLOCATABLE, DIMENSION(:) :: FAN_INDEX     !< List of fans in ductrun
   INTEGER, ALLOCATABLE, DIMENSION(:) :: DUCT_M_INDEX  !< Position of ducts being solved for in ductrun solution matrix
   INTEGER, ALLOCATABLE, DIMENSION(:) :: NODE_M_INDEX  !< Position of ductnodes being solved for in ductrun solution matrix
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: RHO_D      !< Ductrun upstream duct density (kg/m3) (duct,fan)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: RHO_N      !< Ductrun node density (kg/m3) (node,fan)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: ZZ_D     !< Ductrun upstream duct mass fraction (kg/kg) (duct,fan,species)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: ZZ_N     !< Ductrun node mass fraction (kg/kg) (node,fan,species)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TMP_D      !< Ductrun upstream duct temperature (K) (duct,fan)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TMP_N      !< Ductrun node temperature (K) (node,fan)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: CP_D       !< Ductrun upstream duct specific heat (J/kg/K) (duct,fan)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: CP_N       !< Ductrun node specific heat (J/kg/K) (node,fan)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: LOSS       !< Ductrun duct loss (duct)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: VEL      !< Ductrun duct elocity (m/s) (duct,fan,old/new/guess/previous)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: P        !< Ductrun node pressure (Pa) (node,fan,old/new)
END TYPE DUCTRUN_TYPE

TYPE(DUCTRUN_TYPE), DIMENSION(:), ALLOCATABLE, TARGET :: DUCTRUN

END MODULE TYPES
