!> \brief Set important parameters having to do with variable precision and array allocations

MODULE PRECISION_PARAMETERS

IMPLICIT NONE

INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)     !< Precision of "Four Byte" reals
INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)    !< Precision of "Eight Byte" reals
INTEGER, PARAMETER :: MAX_SPECIES=20                 !< Maximum number of declared species
INTEGER, PARAMETER :: MAX_LAYERS=20                  !< Maximum number of solid material layers
INTEGER, PARAMETER :: MAX_MATERIALS=20               !< Maximum number of solid material components
INTEGER, PARAMETER :: MAX_MATERIALS_TOTAL=400        !< Dimension of material work array
INTEGER, PARAMETER :: MAX_REACTIONS=10               !< Maximum number of chemical reactions
INTEGER, PARAMETER :: MAX_STEPS=20                   !< Maximum steps in processing of material residues
INTEGER, PARAMETER :: MAX_NUMBER_SPECTRAL_BANDS=9    !< Maximum number of radiation spectral bands
INTEGER, PARAMETER :: MAX_NUMBER_FSK_POINTS=32
INTEGER, PARAMETER :: MAX_TERRAIN_IMAGES=10          !< Maximum number of images to paste onto complex terrain
INTEGER, PARAMETER :: MAX_INPUT_ID=40                !< Maximum number of CTRL INPUT_IDs
INTEGER, PARAMETER :: N_ZONE_POINTS=100              !< Maximum number of declared ZONE points (deprecated)
INTEGER, PARAMETER :: MAX_AIT_EXCLUSION_ZONES=10     !< Maximum number of AUTO_IGNITION_TEMPERATURE exclusion zones
INTEGER, PARAMETER :: SMOKEVIEW_OBJECTS_DIMENSION=20 !< Number of parameters that can be passed to Smokeview to describe objects
INTEGER, PARAMETER :: LABEL_LENGTH=60                !< Maximum length of most labels
INTEGER, PARAMETER :: MESSAGE_LENGTH=200             !< Maximum length of error and warning labels
INTEGER, PARAMETER :: FORMULA_LENGTH=255             !< Maximum length of chemical formulae
INTEGER, PARAMETER :: CHID_LENGTH=50                 !< Maximum length of job ID
INTEGER, PARAMETER :: MESH_STRING_LENGTH=LABEL_LENGTH + 80 !< Length for storage of strings
INTEGER, PARAMETER :: N_OUTPUT_QUANTITIES=550        !< Dimension of array that holds names of output quantities

REAL(EB), PARAMETER :: ALMOST_ONE=1._EB-EPSILON(1._EB)     !< Largest number that is < 1 in 8 byte accuracy
REAL(EB), PARAMETER :: MICRON=1.E-6_EB                     !< A relatively small length (m)
REAL(EB), PARAMETER :: TWO_EPSILON_EB=2._EB*EPSILON(1._EB) !< A very small number
REAL(EB), PARAMETER :: TINY_EB=TINY(1._EB)                 !< The smallest resolvable 8 byte real number
REAL(EB), PARAMETER :: HUGE_EB=HUGE(1._EB)                 !< The largest resolvable 8 btye real number

! Often used numbers

REAL(EB), PARAMETER :: ONTH=1._EB/3._EB,FOTH=4._EB/3._EB,TWTH=2._EB/3._EB,ONSI=1._EB/6._EB,&
                       SR3=SQRT(3._EB),FTTOT=4._EB*(2._EB/3._EB)**(1._EB/3._EB),EIONTH=18._EB**(1._EB/3._EB)
REAL(EB), PARAMETER :: PI=4._EB*ATAN(1.0_EB), SQRTPI=SQRT(PI), RPI=1._EB/PI, TWOPI=2._EB*PI, PIO2=PI/2._EB, &
                       RFPI=1._EB/(4._EB*PI), FOTHPI = FOTH*PI
INTEGER, PARAMETER  :: INTEGER_ONE=1,INTEGER_ZERO=0
REAL(EB), PARAMETER :: DEG2RAD=4.0_EB*ATAN(1.0_EB)/180.0_EB

END MODULE PRECISION_PARAMETERS
