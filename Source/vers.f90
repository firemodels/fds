
#ifndef GITHASH_PP
#define GITHASH_PP "unknown"
#endif
#ifndef GITDATE_PP
#define GITDATE_PP "unknown"
#endif
#ifndef BUILDDATE_PP
#define BUILDDATE_PP "unknown"
#endif

!> \brief Defines the various parameters that require the local git hash

MODULE GIT_VERSION

IMPLICIT NONE (TYPE,EXTERNAL)

  ! GITHASH_PP, etc, are set by the makefile with -DGITHASH_PP="..."

  CHARACTER(LEN=*), PARAMETER :: GIT_HASH = GITHASH_PP
  CHARACTER(LEN=*), PARAMETER :: GIT_DATE = GITDATE_PP
  CHARACTER(LEN=*), PARAMETER :: BUILD_DATE = BUILDDATE_PP

END MODULE GIT_VERSION

