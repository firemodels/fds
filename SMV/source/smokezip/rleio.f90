#ifdef pp_RLETEST
MODULE RLEIO

USE RLE

  ! RLE FILE FORMAT

  !*** header
  ! endian
  ! fileversion, slice version
  ! global min max (used to perform conversion)
  ! i1,i2,j1,j2,k1,k2


  ! *** frame
   !time, compressed frame size                        for each frame
   !compressed buffer

CONTAINS

SUBROUTINE RLEHEADER(RLEFILE, MINMAX, IJKMINMAX)

!DEC$ ATTRIBUTES ALIAS:'_rleheader@16' :: rleheader

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN) :: RLEFILE
  REAL, INTENT(IN),  DIMENSION(2) :: MINMAX
  INTEGER, INTENT(IN), DIMENSION(6) :: IJKMINMAX

  INTEGER, PARAMETER :: ONE=1
  INTEGER :: I

  OPEN(UNIT=10,FILE=RLEFILE,FORM='UNFORMATTED',STATUS='REPLACE')

 ! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvV
 ! put these four lines into FDS using your unit numbers and variable names

  WRITE(10)ONE
  WRITE(10)ONE,ONE     ! FOR FUTURE COMPATIBILITY
  WRITE(10)MINMAX(1), MINMAX(2)
  WRITE(10)(IJKMINMAX(I),I=1,6)

  ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  CLOSE(10)
END SUBROUTINE RLEHEADER

! STDCALL FORTrleframe(char *rlefile, float *t, float *q, int *nq, float *minmax, unsigned char *cval, int *ncval, int lenrlefile, int lencval);
SUBROUTINE RLEFRAME(RLEFILE,T,Q, NQ, MINMAX, CVAL,NRLE)

!DEC$ ATTRIBUTES ALIAS:'_rleframe@36' :: rleframe

IMPLICIT NONE

CHARACTER(LEN=*), INTENT(IN) :: RLEFILE
REAL, INTENT(IN) :: T
INTEGER, INTENT(IN) :: NQ
REAL, INTENT(IN), DIMENSION(NQ) :: Q
REAL, INTENT(IN), DIMENSION(2) :: MINMAX
INTEGER, INTENT(OUT) :: NRLE
CHARACTER(LEN=1), DIMENSION(NQ), INTENT(OUT) :: CVAL

CHARACTER(LEN=1024) :: SIZEFILE

INTEGER :: I
INTEGER, PARAMETER :: ZERO=0.0

OPEN(UNIT=10,FILE=RLEFILE,FORM='UNFORMATTED',STATUS='OLD',POSITION='APPEND')



! NOTE:  CVAL MUST BE DIMENSIONED AT LEAST AS LARGE AS Q (IE NQ)
!        NRLE IS THE AMOUNT OF SPACE ACTUALLY USED
!        DUE TO COMPRESSION)

! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
!  put these FOUR lines into FDS using your unit numbers and variable names

NRLE = RLECOMPRESS(Q, NQ, MINMAX(1), MINMAX(2), CVAL)

WRITE(10)T                    ! keeping data of the same type on one line makes things easier
WRITE(10)NRLE
WRITE(10)(CVAL(I),I=1,NRLE)
!  WRITE(11,*)T,ZERO,NRLE     ! also write out info to a sizing file (the zero is a placeholder for the zlib compression info)

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
CLOSE(10)


END SUBROUTINE RLEFRAME
END MODULE RLEIO
#else
MODULE RLEIO
END MODULE RLEIO
#endif
