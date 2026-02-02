!> \brief Set important parameters having to do with variable precision and array allocations

MODULE PRECISION_PARAMETERS

IMPLICIT NONE (TYPE,EXTERNAL)

INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)     !< Precision of "Four Byte" reals
INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)    !< Precision of "Eight Byte" reals
INTEGER, PARAMETER :: QB = SELECTED_REAL_KIND(33,4931) !< Precision of "Sixteen Byte" reals
INTEGER, PARAMETER :: MAX_LPC=20                     !< Maximum number of declared particle classes
INTEGER, PARAMETER :: MAX_SPECIES=20                 !< Maximum number of declared species
INTEGER, PARAMETER :: MAX_LAYERS=20                  !< Maximum number of solid material layers
INTEGER, PARAMETER :: MAX_LAYERS_HT3D=500            !< Maximum number of solid material layers for an HT3D solid
INTEGER, PARAMETER :: MAX_MATERIALS=20               !< Maximum number of solid material components
INTEGER, PARAMETER :: MAX_MATERIALS_TOTAL=400        !< Dimension of material work array
INTEGER, PARAMETER :: MAX_CONE_CURVES=10             !< Maximum number of cone calorimeter curves
INTEGER, PARAMETER :: MAX_REACTIONS=10               !< Maximum number of chemical reactions
INTEGER, PARAMETER :: MAX_STEPS=20                   !< Maximum steps in processing of material residues
INTEGER, PARAMETER :: MAX_NUMBER_SPECTRAL_BANDS=9    !< Maximum number of radiation spectral bands
INTEGER, PARAMETER :: MAX_TERRAIN_IMAGES=10          !< Maximum number of images to paste onto complex terrain
INTEGER, PARAMETER :: MAX_INPUT_ID=40                !< Maximum number of CTRL INPUT_IDs
INTEGER, PARAMETER :: N_ZONE_POINTS=100              !< Maximum number of declared ZONE points (deprecated)
INTEGER, PARAMETER :: MAX_AIT_EXCLUSION_ZONES=10     !< Maximum number of AUTO_IGNITION_TEMPERATURE exclusion zones
INTEGER, PARAMETER :: MAX_IGNITION_ZONES=10          !< Maximum number of Ignition zones
INTEGER, PARAMETER :: SMOKEVIEW_OBJECTS_DIMENSION=20 !< Number of parameters that can be passed to Smokeview to describe objects
INTEGER, PARAMETER :: LABEL_LENGTH=60                !< Maximum length of most labels
INTEGER, PARAMETER :: MESSAGE_LENGTH=200             !< Maximum length of error and warning labels
INTEGER, PARAMETER :: FORMULA_LENGTH=255             !< Maximum length of chemical formulae
INTEGER, PARAMETER :: CHID_LENGTH=50                 !< Maximum length of job ID
INTEGER, PARAMETER :: FILE_LENGTH=200                !< Maximum length of filenames ID
INTEGER, PARAMETER :: NAMELIST_LENGTH=300            !< Maximum length of NAMELIST line in input file (see SEARCH_INPUT_FILE)
INTEGER, PARAMETER :: MESH_STRING_LENGTH=LABEL_LENGTH + 100 !< Length for storage of strings
INTEGER, PARAMETER :: FN_LENGTH=FILE_LENGTH+CHID_LENGTH     !< Length for output filename strings (includes output directory)
INTEGER, PARAMETER :: N_OUTPUT_QUANTITIES=600        !< Dimension of array that holds names of output quantities
INTEGER, PARAMETER :: POINTS_ARRAY_DIM=100           !< Dimension of arrays of linear device coordinates

REAL(EB), PARAMETER :: ONE_M_EPS=1._EB-100._EB*EPSILON(1._EB)  !< Number that is slightly less than 1
REAL(EB), PARAMETER :: ONE_P_EPS=1._EB+100._EB*EPSILON(1._EB)  !< Number that is slightly greater than 1
REAL(EB), PARAMETER :: MICRON=1.E-6_EB                         !< A relatively small length (m)
REAL(EB), PARAMETER :: TWENTY_EPSILON_EB=20._EB*EPSILON(1._EB) !< A very small number 8 byte accuracy
REAL(EB), PARAMETER :: TEN_EPSILON_EB   =10._EB*EPSILON(1._EB) !< A very small number 8 byte accuracy
REAL(EB), PARAMETER :: TWO_EPSILON_EB   = 2._EB*EPSILON(1._EB) !< A very small number 8 byte accuracy
REAL(EB), PARAMETER :: TINY_EB=TINY(1._EB)                     !< The smallest resolvable 8 byte real number
REAL(EB), PARAMETER :: HUGE_EB=HUGE(1._EB)                     !< The largest resolvable 8 btye real number
REAL(EB), PARAMETER :: HUGE_FB=HUGE(1._FB)                     !< A large number but not too large for various operations

! Often used numbers

REAL(EB), PARAMETER :: ONTH=1._EB/3._EB,FOTH=4._EB/3._EB,TWTH=2._EB/3._EB,ONSI=1._EB/6._EB,&
                       SR2=SQRT(2._EB),SR3=SQRT(3._EB),FTTOT=4._EB*(2._EB/3._EB)**(1._EB/3._EB),EIONTH=18._EB**(1._EB/3._EB)
REAL(EB), PARAMETER :: PI=4._EB*ATAN(1.0_EB), SQRTPI=SQRT(PI), RPI=1._EB/PI, TWOPI=2._EB*PI, PIO2=PI/2._EB, &
                       RFPI=1._EB/(4._EB*PI), FOTHPI = FOTH*PI, CR2=2._EB**(1._EB/3._EB)
INTEGER, PARAMETER  :: INTEGER_ZERO=0,INTEGER_ONE=1,INTEGER_TWO=2,INTEGER_THREE=3
REAL(EB), PARAMETER :: DEG2RAD=PI/180.0_EB

END MODULE PRECISION_PARAMETERS
