PROGRAM LOOP3D

USE OMP_LIB
IMPLICIT NONE (TYPE,EXTERNAL)

! Miscellaneous declarations
INTEGER, PARAMETER :: EB = SELECTED_REAL_KIND(12)
INTEGER, PARAMETER :: IBAR = 256, JBAR = 256, KBAR = 256
INTEGER, PARAMETER :: NEDGE = 12
INTEGER, PARAMETER :: IBP1 = IBAR+1, JBP1 = JBAR+1, KBP1 = KBAR+1
REAL(EB), PARAMETER :: FOTH = 4.0_EB/3.0_EB

REAL(EB), ALLOCATABLE, DIMENSION(:) :: GX, GY, GZ, RHO_0, X, Y, Z, RDXN, RDYN, RDZN, RDX, RDY, RDZ
REAL(EB), ALLOCATABLE, TARGET, DIMENSION(:,:,:) :: U, V, W, D, RHO, WORK1, WORK2, WORK3, WORK4, WORK5, WORK6
REAL(EB), POINTER, DIMENSION(:,:,:) :: DP, RHOP, UU, VV, WW, OMY, OMX, OMZ, TXZ, TXY, TYZ
REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU, FVX, FVY, FVZ
INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: CELL_INDEX
REAL(EB) :: MUX,MUY,MUZ,UP,UM,VP,VM,WP,WM,VTRM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXYP,TXYM,TXZP,TXZM,TYZP,TYZM, &
            DTXYDY,DTXZDZ,DTYZDZ,DTXYDX,DTXZDX,DTYZDY, &
            DUDX,DVDY,DWDZ,DUDY,DUDZ,DVDX,DVDZ,DWDX,DWDY, &
            VOMZ,WOMY,UOMY,VOMX,UOMZ,WOMX, &
            RRHO,TXXP,TXXM,TYYP,TYYM,TZZP,TZZM,DTXXDX,DTYYDY,DTZZDZ,T_NOW,T_END
INTEGER :: I,J,K,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,IC,IC1,IC2,IE,MAX_EDGE,NT
CHARACTER(LEN=50) :: FILENAME

TYPE CELL_TYPE
   INTEGER :: EDGE_INDEX(NEDGE)=0
END TYPE CELL_TYPE
TYPE(CELL_TYPE), ALLOCATABLE, DIMENSION(:) :: CELL

TYPE EDGE_TYPE
   REAL(EB) :: OMEGA(-2:2)=-1.E6_EB
   REAL(EB) :: TAU(-2:2)=-1.E6_EB
END TYPE EDGE_TYPE
TYPE(EDGE_TYPE), ALLOCATABLE, DIMENSION(:) :: EDGE

! Write out Starting:
!$OMP PARALLEL
!$OMP MASTER
!$ NT = OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP BARRIER
!$OMP END PARALLEL

WRITE(FILENAME,'(A,I4.4,A,I2.2,A)') 'loop3d_',IBAR,'_',NT,'THR.txt'
WRITE(*,*) 'Starting Loop3D, out file: ',TRIM(FILENAME)
OPEN(UNIT=10,FILE=TRIM(FILENAME),STATUS='UNKNOWN')
WRITE(10,*) 'Starting Loop3D'
WRITE(10,*) 'IBAR=',IBAR,' JBAR=',JBAR,' KBAR=',KBAR,' OMP_NUM_THREADS=',OMP_GET_NUM_THREADS()

! Allocate vars in CPU:
ALLOCATE(X(0:IBAR), Y(0:JBAR), Z(0:KBAR), RDXN(0:IBAR), RDYN(0:JBAR), RDZN(0:KBAR), RDX(0:IBAR), RDY(0:JBAR), RDZ(0:KBAR))
ALLOCATE(U(0:IBP1,0:JBP1,0:KBP1), V(0:IBP1,0:JBP1,0:KBP1), W(0:IBP1,0:JBP1,0:KBP1), MU(0:IBP1,0:JBP1,0:KBP1) )
ALLOCATE(D(0:IBP1,0:JBP1,0:KBP1), RHO(0:IBP1,0:JBP1,0:KBP1), WORK1(0:IBP1,0:JBP1,0:KBP1))
ALLOCATE(WORK2(0:IBP1,0:JBP1,0:KBP1), WORK3(0:IBP1,0:JBP1,0:KBP1), WORK4(0:IBP1,0:JBP1,0:KBP1))
ALLOCATE(WORK5(0:IBP1,0:JBP1,0:KBP1), WORK6(0:IBP1,0:JBP1,0:KBP1), DP(0:IBP1,0:JBP1,0:KBP1), RHOP(0:IBP1,0:JBP1,0:KBP1))
ALLOCATE(RHO_0(0:KBAR), GX(0:IBAR), GY(0:IBAR), GZ(0:IBAR))
ALLOCATE(FVX(0:IBP1,0:JBP1,0:KBP1), FVY(0:IBP1,0:JBP1,0:KBP1), FVZ(0:IBP1,0:JBP1,0:KBP1))
ALLOCATE(CELL_INDEX(0:IBP1,0:JBP1,0:KBP1))

! Initialize:
DO I=0,IBAR
   X(I)    = REAL(I,EB)
   RDXN(I) = 1.0_EB
   RDX(I)  = 1.0_EB
ENDDO
DO J=0,JBAR
   Y(J)    = REAL(J,EB)
   RDYN(J) = 1.0_EB
   RDY(J)  = 1.0_EB
ENDDO
DO K=0,KBAR
   Z(K)    = REAL(K,EB)
   RDZN(K) = 1.0_EB
   RDZ(K)  = 1.0_EB
ENDDO

! Cell Index, CELL and EDGE:
CELL_INDEX = 0
IC = 0
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
         IC = IC + 1
         CELL_INDEX(I,J,K) = IC
      ENDDO
   ENDDO
ENDDO
ALLOCATE(CELL(0:IC))

IC = 0
MAX_EDGE=-1
DO K=1,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         IF( .NOT. (ANY( K==(/1,KBAR/) ) .OR. ANY( J==(/1,JBAR/) ) .OR. ANY( I==(/1,IBAR/) )) ) CYCLE
         IC = IC + 1
         CELL(IC)%EDGE_INDEX(1)  = CELL_INDEX(I  ,J  ,K  ) + 1
         CELL(IC)%EDGE_INDEX(2)  = CELL_INDEX(I+1,J  ,K  ) + 1
         CELL(IC)%EDGE_INDEX(3)  = CELL_INDEX(I  ,J  ,K  ) + 2
         CELL(IC)%EDGE_INDEX(4)  = CELL_INDEX(I  ,J+1,K  ) + 2
         CELL(IC)%EDGE_INDEX(5)  = CELL_INDEX(I  ,J  ,K  ) + 3
         CELL(IC)%EDGE_INDEX(6)  = CELL_INDEX(I  ,J  ,K-1) + 3
         CELL(IC)%EDGE_INDEX(7)  = CELL_INDEX(I  ,J  ,K  ) + 4
         CELL(IC)%EDGE_INDEX(8)  = CELL_INDEX(I  ,J+1,K  ) + 4
         CELL(IC)%EDGE_INDEX(9)  = CELL_INDEX(I  ,J  ,K  ) + 5
         CELL(IC)%EDGE_INDEX(10) = CELL_INDEX(I  ,J  ,K-1) + 5
         CELL(IC)%EDGE_INDEX(11) = CELL_INDEX(I  ,J  ,K  ) + 6
         CELL(IC)%EDGE_INDEX(12) = CELL_INDEX(I+1,J  ,K  ) + 6
         DO IE=1,NEDGE
            MAX_EDGE = MAX(MAX_EDGE,CELL(IC)%EDGE_INDEX(IE))
         ENDDO
      ENDDO
   ENDDO
ENDDO
ALLOCATE(EDGE(0:MAX_EDGE))
DO IE=1,MAX_EDGE,2
   EDGE(IE)%OMEGA = 1.5E-4_EB
   EDGE(IE)%TAU   = 2.5E-4_EB
ENDDO

RHO   = 1.19_EB
RHO_0 = 1.19_EB
D     = 0.0015_EB
MU    = 0.0019_EB
WORK1 = 0.0_EB; WORK2 = 0.0_EB; WORK3 = 0.0_EB; WORK4 = 0.0_EB; WORK5 = 0.0_EB; WORK6 = 0.0_EB
GX(:) = 0.0_EB; GY(:) = 0.0_EB; GZ(:) = 1.0_EB

! U, V, W:
DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         ! Some Trig functions:
         U(I,J,K) =  SIN(X(I))*COS(Y(J))*COS(Z(K))
         V(I,J,K) = -COS(X(I))*SIN(Y(J))*COS(Z(K))
         W(I,J,K) =  COS(X(I))*COS(Y(J))*SIN(Z(K))
      ENDDO
   ENDDO
ENDDO

! Compute Tau OMG:
UU => U
VV => V
WW => W
DP => D
RHOP => RHO
TXY => WORK1
TXZ => WORK2
TYZ => WORK3
OMX => WORK4
OMY => WORK5
OMZ => WORK6

DO K=0,KBAR
   DO J=0,JBAR
      DO I=0,IBAR
         DUDY = RDYN(J)*(UU(I,J+1,K)-UU(I,J,K))
         DVDX = RDXN(I)*(VV(I+1,J,K)-VV(I,J,K))
         DUDZ = RDZN(K)*(UU(I,J,K+1)-UU(I,J,K))
         DWDX = RDXN(I)*(WW(I+1,J,K)-WW(I,J,K))
         DVDZ = RDZN(K)*(VV(I,J,K+1)-VV(I,J,K))
         DWDY = RDYN(J)*(WW(I,J+1,K)-WW(I,J,K))
         OMX(I,J,K) = DWDY - DVDZ
         OMY(I,J,K) = DUDZ - DWDX
         OMZ(I,J,K) = DVDX - DUDY
         MUX = 0.25_EB*(MU(I,J+1,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I,J+1,K+1))
         MUY = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J,K+1)+MU(I+1,J,K+1))
         MUZ = 0.25_EB*(MU(I+1,J,K)+MU(I,J,K)+MU(I,J+1,K)+MU(I+1,J+1,K))
         TXY(I,J,K) = MUZ*(DVDX + DUDY)
         TXZ(I,J,K) = MUY*(DUDZ + DWDX)
         TYZ(I,J,K) = MUX*(DVDZ + DWDY)
      ENDDO
   ENDDO
ENDDO

CALL CPU_TIME(T_NOW)
CALL LOOP3D_OMP_CPU()
CALL CPU_TIME(T_END)

WRITE(10,*) 'Time=',T_END-T_NOW
WRITE(10,*) 'mean FVX =',SUM(FVX(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
WRITE(10,*) 'mean FVY =',SUM(FVY(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
WRITE(10,*) 'mean FVZ =',SUM(FVZ(1:IBAR,1:JBAR,1:KBAR))/(IBAR*JBAR*KBAR)
WRITE(10,*) 'Ending Loop3D'
CLOSE(10)
WRITE(*,*) 'Loop3D done.'

CONTAINS

SUBROUTINE LOOP3D_OMP_CPU()

! Compute x-direction flux term FVX

!$OMP PARALLEL PRIVATE(WP,WM,VP,VM,UP,UM,OMXP,OMXM,OMYP,OMYM,OMZP,OMZM,TXZP,TXZM,TXYP,TXYM,TYZP,TYZM, &
!$OMP& IC,IEXP,IEXM,IEYP,IEYM,IEZP,IEZM,RRHO,DUDX,DVDY,DWDZ,VTRM)

!$OMP DO SCHEDULE(STATIC) PRIVATE(WOMY, VOMZ, TXXP, TXXM, DTXXDX, DTXYDY, DTXZDZ)
DO K=1,KBAR
   DO J=1,JBAR
      DO I=0,IBAR
         WP    = WW(I,J,K)   + WW(I+1,J,K)
         WM    = WW(I,J,K-1) + WW(I+1,J,K-1)
         VP    = VV(I,J,K)   + VV(I+1,J,K)
         VM    = VV(I,J-1,K) + VV(I+1,J-1,K)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I,J-1,K)
         IC    = CELL_INDEX(I,J,K)
         IEYP  = CELL(IC)%EDGE_INDEX(8)
         IEYM  = CELL(IC)%EDGE_INDEX(6)
         IEZP  = CELL(IC)%EDGE_INDEX(12)
         IEZM  = CELL(IC)%EDGE_INDEX(10)
         IF (EDGE(IEYP)%OMEGA(-1)>-1.E5_EB) THEN
            OMYP = EDGE(IEYP)%OMEGA(-1)
            TXZP = EDGE(IEYP)%TAU(-1)
         ENDIF
         IF (EDGE(IEYM)%OMEGA( 1)>-1.E5_EB) THEN
            OMYM = EDGE(IEYM)%OMEGA( 1)
            TXZM = EDGE(IEYM)%TAU( 1)
         ENDIF
         IF (EDGE(IEZP)%OMEGA(-2)>-1.E5_EB) THEN
            OMZP = EDGE(IEZP)%OMEGA(-2)
            TXYP = EDGE(IEZP)%TAU(-2)
         ENDIF
         IF (EDGE(IEZM)%OMEGA( 2)>-1.E5_EB) THEN
            OMZM = EDGE(IEZM)%OMEGA( 2)
            TXYM = EDGE(IEZM)%TAU( 2)
         ENDIF
         WOMY  = WP*OMYP + WM*OMYM
         VOMZ  = VP*OMZP + VM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I+1,J,K))
         DVDY  = (VV(I+1,J,K)-VV(I+1,J-1,K))*RDY(J)
         DWDZ  = (WW(I+1,J,K)-WW(I+1,J,K-1))*RDZ(K)
         TXXP  = MU(I+1,J,K)*( FOTH*DP(I+1,J,K) - 2._EB*(DVDY+DWDZ) )
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         TXXM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DVDY+DWDZ) )
         DTXXDX= RDXN(I)*(TXXP-TXXM)
         DTXYDY= RDY(J) *(TXYP-TXYM)
         DTXZDZ= RDZ(K) *(TXZP-TXZM)
         VTRM  = DTXXDX + DTXYDY + DTXZDZ
         FVX(I,J,K) = 0.25_EB*(WOMY - VOMZ) - GX(I) + RRHO*(GX(I)*RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

! Compute y-direction flux term FVY

!$OMP DO SCHEDULE(STATIC) PRIVATE(WOMX, UOMZ, TYYP, TYYM, DTXYDX, DTYYDY, DTYZDZ)
DO K=1,KBAR
   DO J=0,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J+1,K)
         UM    = UU(I-1,J,K) + UU(I-1,J+1,K)
         WP    = WW(I,J,K)   + WW(I,J+1,K)
         WM    = WW(I,J,K-1) + WW(I,J+1,K-1)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J,K-1)
         OMZP  = OMZ(I,J,K)
         OMZM  = OMZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J,K-1)
         TXYP  = TXY(I,J,K)
         TXYM  = TXY(I-1,J,K)
         IC    = CELL_INDEX(I,J,K)
         IEXP  = CELL(IC)%EDGE_INDEX(4)
         IEXM  = CELL(IC)%EDGE_INDEX(2)
         IEZP  = CELL(IC)%EDGE_INDEX(12)
         IEZM  = CELL(IC)%EDGE_INDEX(11)
         IF (EDGE(IEXP)%OMEGA(-2)>-1.E5_EB) THEN
            OMXP = EDGE(IEXP)%OMEGA(-2)
            TYZP = EDGE(IEXP)%TAU(-2)
         ENDIF
         IF (EDGE(IEXM)%OMEGA( 2)>-1.E5_EB) THEN
            OMXM = EDGE(IEXM)%OMEGA( 2)
            TYZM = EDGE(IEXM)%TAU( 2)
         ENDIF
         IF (EDGE(IEZP)%OMEGA(-1)>-1.E5_EB) THEN
            OMZP = EDGE(IEZP)%OMEGA(-1)
            TXYP = EDGE(IEZP)%TAU(-1)
         ENDIF
         IF (EDGE(IEZM)%OMEGA( 1)>-1.E5_EB) THEN
            OMZM = EDGE(IEZM)%OMEGA( 1)
            TXYM = EDGE(IEZM)%TAU( 1)
         ENDIF
         WOMX  = WP*OMXP + WM*OMXM
         UOMZ  = UP*OMZP + UM*OMZM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J+1,K))
         DUDX  = (UU(I,J+1,K)-UU(I-1,J+1,K))*RDX(I)
         DWDZ  = (WW(I,J+1,K)-WW(I,J+1,K-1))*RDZ(K)
         TYYP  = MU(I,J+1,K)*( FOTH*DP(I,J+1,K) - 2._EB*(DUDX+DWDZ) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
         DWDZ  = (WW(I,J,K)-WW(I,J,K-1))*RDZ(K)
         TYYM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DWDZ) )
         DTXYDX= RDX(I) *(TXYP-TXYM)
         DTYYDY= RDYN(J)*(TYYP-TYYM)
         DTYZDZ= RDZ(K) *(TYZP-TYZM)
         VTRM  = DTXYDX + DTYYDY + DTYZDZ
         FVY(I,J,K) = 0.25_EB*(UOMZ - WOMX) - GY(I) + RRHO*(GY(I)*RHO_0(K) - VTRM)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

! Compute z-direction flux term FVZ

!$OMP DO SCHEDULE(STATIC) PRIVATE(UOMY, VOMX, TZZP, TZZM, DTXZDX, DTYZDY, DTZZDZ)
DO K=0,KBAR
   DO J=1,JBAR
      DO I=1,IBAR
         UP    = UU(I,J,K)   + UU(I,J,K+1)
         UM    = UU(I-1,J,K) + UU(I-1,J,K+1)
         VP    = VV(I,J,K)   + VV(I,J,K+1)
         VM    = VV(I,J-1,K) + VV(I,J-1,K+1)
         OMYP  = OMY(I,J,K)
         OMYM  = OMY(I-1,J,K)
         OMXP  = OMX(I,J,K)
         OMXM  = OMX(I,J-1,K)
         TXZP  = TXZ(I,J,K)
         TXZM  = TXZ(I-1,J,K)
         TYZP  = TYZ(I,J,K)
         TYZM  = TYZ(I,J-1,K)
         IC    = CELL_INDEX(I,J,K)
         IEXP  = CELL(IC)%EDGE_INDEX(4)
         IEXM  = CELL(IC)%EDGE_INDEX(3)
         IEYP  = CELL(IC)%EDGE_INDEX(8)
         IEYM  = CELL(IC)%EDGE_INDEX(7)
         IF (EDGE(IEXP)%OMEGA(-1)>-1.E5_EB) THEN
            OMXP = EDGE(IEXP)%OMEGA(-1)
            TYZP = EDGE(IEXP)%TAU(-1)
         ENDIF
         IF (EDGE(IEXM)%OMEGA( 1)>-1.E5_EB) THEN
            OMXM = EDGE(IEXM)%OMEGA( 1)
            TYZM = EDGE(IEXM)%TAU( 1)
         ENDIF
         IF (EDGE(IEYP)%OMEGA(-2)>-1.E5_EB) THEN
            OMYP = EDGE(IEYP)%OMEGA(-2)
            TXZP = EDGE(IEYP)%TAU(-2)
         ENDIF
         IF (EDGE(IEYM)%OMEGA( 2)>-1.E5_EB) THEN
            OMYM = EDGE(IEYM)%OMEGA( 2)
            TXZM = EDGE(IEYM)%TAU( 2)
         ENDIF
         UOMY  = UP*OMYP + UM*OMYM
         VOMX  = VP*OMXP + VM*OMXM
         RRHO  = 2._EB/(RHOP(I,J,K)+RHOP(I,J,K+1))
         DUDX  = (UU(I,J,K+1)-UU(I-1,J,K+1))*RDX(I)
         DVDY  = (VV(I,J,K+1)-VV(I,J-1,K+1))*RDY(J)
         TZZP  = MU(I,J,K+1)*( FOTH*DP(I,J,K+1) - 2._EB*(DUDX+DVDY) )
         DUDX  = (UU(I,J,K)-UU(I-1,J,K))*RDX(I)
         DVDY  = (VV(I,J,K)-VV(I,J-1,K))*RDY(J)
         TZZM  = MU(I,J,K)  *( FOTH*DP(I,J,K)   - 2._EB*(DUDX+DVDY) )
         DTXZDX= RDX(I) *(TXZP-TXZM)
         DTYZDY= RDY(J) *(TYZP-TYZM)
         DTZZDZ= RDZN(K)*(TZZP-TZZM)
         VTRM  = DTXZDX + DTYZDY + DTZZDZ
         FVZ(I,J,K) = 0.25_EB*(VOMX - UOMY) - GZ(I) + RRHO*(GZ(I)*0.5_EB*(RHO_0(K)+RHO_0(K+1)) - VTRM)
      ENDDO
   ENDDO
ENDDO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

END SUBROUTINE LOOP3D_OMP_CPU

END PROGRAM LOOP3D

