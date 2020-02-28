MODULE PRECISION_PARAMETERS

! Set important parameters having to do with variable precision and array allocations

! A partial list of the parameters defined here:
! ONTH ... Self evident ratios
! SR3      square root of 3
! EIONTH   cubic root of 18
! PI       the number pi = 3.14159265359
! TWOPI    the number 2pi = 6.28318530718
! PIO2     pi over 2, approx. 1.57079632679
! SQRTPI   square root of pi, approx. 1.77245385091
! RPI      reciprocal pi, approx. 0.31830988618


IMPLICIT NONE

! Precision of "Four Byte" and "Eight Byte" reals

INTEGER, PARAMETER :: FB = SELECTED_REAL_KIND(6)
INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)

! Single- and double-precision complex

INTEGER, PARAMETER :: SPC = KIND((1._FB,1._FB))
INTEGER, PARAMETER :: DPC = KIND((1._EB,1._EB))

! Hardwired bounds for certain species arrays

INTEGER, PARAMETER :: MAX_SPECIES=20

! Hardwired bounds for surface and material arrays

INTEGER, PARAMETER :: MAX_LAYERS=20, MAX_MATERIALS=20, MAX_MATERIALS_TOTAL=400, MAX_REACTIONS=10, MAX_STEPS=20
INTEGER, PARAMETER :: MAX_NUMBER_SPECTRAL_BANDS=9, MAX_NUMBER_FSK_POINTS=32, N_ZONE_POINTS=100, MAX_AIT_EXCLUSION_ZONES=10

! Hardwired number of parameters that can be passed to Smokeview to describe a drawn object or device

INTEGER, PARAMETER :: SMOKEVIEW_OBJECTS_DIMENSION=20

! Hardwired length of most labels

INTEGER, PARAMETER :: LABEL_LENGTH=60, MESSAGE_LENGTH=200, FORMULA_LENGTH=255, CHID_LENGTH=50
! Size of M%STRING to accomodate SMV SLCF lines data + ~60 character SLCF ID.
INTEGER, PARAMETER :: MESH_STRING_LENGTH=LABEL_LENGTH + 80

! Hardwired length of output quantities

INTEGER, PARAMETER :: N_OUTPUT_QUANTITIES=550

! Special numbers
! Numbers such as the largest number that is < 1 in 8-byte accuracy (ALMOST_ONE) are defined here

REAL(EB), PARAMETER :: ALMOST_ONE=1._EB-EPSILON(1._EB),MICRON=1.E-6_EB,NANOMETER=1.E-9_EB,&
                       TWO_EPSILON_EB=2._EB*EPSILON(1._EB),TINY_EB=TINY(1._EB),HUGE_EB=HUGE(1._EB)

! Often used numbers

REAL(EB), PARAMETER :: ONTH=1._EB/3._EB,FOTH=4._EB/3._EB,TWTH=2._EB/3._EB,ONSI=1._EB/6._EB,&
                       SR3=SQRT(3._EB),FTTOT=4._EB*(2._EB/3._EB)**(1._EB/3._EB),EIONTH=18._EB**(1._EB/3._EB)
REAL(EB), PARAMETER :: PI=4._EB*ATAN(1.0_EB), SQRTPI=SQRT(PI), RPI=1._EB/PI, TWOPI=2._EB*PI, PIO2=PI/2._EB, &
                       RFPI=1._EB/(4._EB*PI), FOTHPI = FOTH*PI
INTEGER, PARAMETER  :: INTEGER_ONE=1
REAL(EB), PARAMETER :: DEG2RAD=4.0_EB*ATAN(1.0_EB)/180.0_EB

END MODULE PRECISION_PARAMETERS
