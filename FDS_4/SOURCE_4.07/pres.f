      MODULE PRES
C
C Find the perturbation pressure by solving Poisson's Equation
C
      USE PREC
      USE CONS
      USE PACKER
      USE POIS
      IMPLICIT NONE
      PRIVATE
      PUBLIC PRESSURE_SOLVER
C
C
      CONTAINS
C
      SUBROUTINE PRESSURE_SOLVER(NM)
C
      INTEGER, INTENT(IN) :: NM
      REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW
      REAL(EB), POINTER, DIMENSION(:) :: UWP
      INTEGER :: I,J,K,IW,IOR
      REAL(EB) :: TRM1,TRM2,TRM3,TRM4,RES,LHSS,RHSS,HH,
     .            DWDT,DVDT,DUDT,HQ2,RFODT,U2,V2,W2,HFAC,H0RR(6)
      LOGICAL :: GET_H
      CHARACTER(10) :: BC_TYPE
C
      CALL UNPACK_VAR(NM)
C
      TNOW=SECOND()
C
      IF (PREDICTOR) THEN
         UU => U
         VV => V
         WW => W
         UWP=> UWS
         ELSE
         UU => US
         VV => VS
         WW => WS
         UWP=> UW
         ENDIF
C
      RFODT = RF/DT
      HFAC  = 1.-RHOA/RHO_AVG
C
      H0RR = 0.
      IF (U0.GE.0.) H0RR(1) = H0*RHOA/RHO_AVG
      IF (U0.LE.0.) H0RR(2) = H0*RHOA/RHO_AVG
      IF (V0.GE.0.) H0RR(3) = H0*RHOA/RHO_AVG
      IF (V0.LE.0.) H0RR(4) = H0*RHOA/RHO_AVG
      IF (W0.GE.0.) H0RR(5) = H0*RHOA/RHO_AVG
      IF (W0.LE.0.) H0RR(6) = H0*RHOA/RHO_AVG
C
C Apply pressure boundary conditions at external cells.
C If Neumann, BXS, BXF, etc., contain dH/dx(x=XS), dH/dx(x=XF), etc.
C If Dirichlet, BXS, BXF, etc., contain H(x=XS), H(x=XF), etc.
C LBC, MBC and NBC are codes used be Poisson solver to denote type
C of boundary condition at x, y and z boundaries. See Crayfishpak
C manual for details.
C
      WALL_CELL_LOOP: DO IW=1,NEWC
C
      I   = IJKW(1,IW)
      J   = IJKW(2,IW)
      K   = IJKW(3,IW)
      IOR = IJKW(4,IW)
C
      SELECT CASE(IOR)
      CASE( 1)
      IF (LBC.EQ.3 .OR. LBC.EQ.4 .OR. LBC.EQ.6) BC_TYPE = 'NEUMANN'
      IF (LBC.EQ.1 .OR. LBC.EQ.2 .OR. LBC.EQ.5) BC_TYPE = 'DIRICHLET'
      CASE(-1)
      IF (LBC.EQ.2 .OR. LBC.EQ.3 .OR. LBC.EQ.6) BC_TYPE = 'NEUMANN'
      IF (LBC.EQ.1 .OR. LBC.EQ.4 .OR. LBC.EQ.5) BC_TYPE = 'DIRICHLET'
      CASE( 2)
      IF (MBC.EQ.3 .OR. MBC.EQ.4) BC_TYPE = 'NEUMANN'
      IF (MBC.EQ.1 .OR. MBC.EQ.2) BC_TYPE = 'DIRICHLET'
      CASE(-2)
      IF (MBC.EQ.3 .OR. MBC.EQ.2) BC_TYPE = 'NEUMANN'
      IF (MBC.EQ.1 .OR. MBC.EQ.4) BC_TYPE = 'DIRICHLET'
      CASE( 3)
      IF (NBC.EQ.3 .OR. NBC.EQ.4) BC_TYPE = 'NEUMANN'
      IF (NBC.EQ.1 .OR. NBC.EQ.2) BC_TYPE = 'DIRICHLET'
      CASE(-3)
      IF (NBC.EQ.3 .OR. NBC.EQ.2) BC_TYPE = 'NEUMANN'
      IF (NBC.EQ.1 .OR. NBC.EQ.4) BC_TYPE = 'DIRICHLET'
      END SELECT
C
      NEUMANN: IF (BC_TYPE.EQ.'NEUMANN') THEN
C
      SELECT CASE(IOR)
      CASE( 1) ; BXS(J,K) = HX(0)   *(-FVX(0,J,K)    + DUWDT(IW))
      CASE(-1) ; BXF(J,K) = HX(IBP1)*(-FVX(IBAR,J,K) - DUWDT(IW))
      CASE( 2) ; BYS(I,K) = HY(0)   *(-FVY(I,0,K)    + DUWDT(IW))
      CASE(-2) ; BYF(I,K) = HY(JBP1)*(-FVY(I,JBAR,K) - DUWDT(IW))
      CASE( 3) ; BZS(I,J) = HZ(0)   *(-FVZ(I,J,0)    + DUWDT(IW))
      CASE(-3) ; BZF(I,J) = HZ(KBP1)*(-FVZ(I,J,KBAR) - DUWDT(IW))
      END SELECT
C
      ENDIF NEUMANN
C
      DIRICHLET: IF (BC_TYPE.EQ.'DIRICHLET') THEN
C
      NOT_OPEN: IF (IV(IW).NE.2) THEN
C
      SELECT CASE(IOR)
      CASE( 1) ; DUDT = -RFODT*(UU(0,J,K)   +UWP(IW)) ; HH = H(1,J,K)
      CASE(-1) ; DUDT = -RFODT*(UU(IBAR,J,K)-UWP(IW)) ; HH = H(IBAR,J,K)
      CASE( 2) ; DVDT = -RFODT*(VV(I,0,K)   +UWP(IW)) ; HH = H(I,1,K)
      CASE(-2) ; DVDT = -RFODT*(VV(I,JBAR,K)-UWP(IW)) ; HH = H(I,JBAR,K)
      CASE( 3) ; DWDT = -RFODT*(WW(I,J,0)   +UWP(IW)) ; HH = H(I,J,1)
      CASE(-3) ; DWDT = -RFODT*(WW(I,J,KBAR)-UWP(IW)) ; HH = H(I,J,KBAR)
      END SELECT
C
      GET_H = .FALSE.
C
      IF (IV(IW).EQ.4) THEN
      SELECT CASE(IOR)
      CASE( 1); IF (UU(0,J,K)   .GT.0. .AND. UWP(IW).LT.0.) GET_H=.TRUE.
      CASE(-1); IF (UU(IBAR,J,K).LT.0. .AND. UWP(IW).LT.0.) GET_H=.TRUE.
      CASE( 2); IF (VV(I,0,K)   .GT.0. .AND. UWP(IW).LT.0.) GET_H=.TRUE.
      CASE(-2); IF (VV(I,JBAR,K).LT.0. .AND. UWP(IW).LT.0.) GET_H=.TRUE.
      CASE( 3); IF (WW(I,J,0)   .GT.0. .AND. UWP(IW).LT.0.) GET_H=.TRUE.
      CASE(-3); IF (WW(I,J,KBAR).LT.0. .AND. UWP(IW).LT.0.) GET_H=.TRUE.
      END SELECT
      IF (GET_H)
     .HH = OMESH(IJKW(9,IW))%H(IJKW(10,IW),IJKW(11,IW),IJKW(12,IW))
      ENDIF
C
      SELECT CASE(IOR)
      CASE( 1) ; BXS(J,K) = HH + 0.5*DX(0)   *(DUDT+FVX(0,J,K))
      IF (GET_H.AND.IV(IW).EQ.4) BXS(J,K) = HH 
      CASE(-1) ; BXF(J,K) = HH - 0.5*DX(IBP1)*(DUDT+FVX(IBAR,J,K))
      IF (GET_H.AND.IV(IW).EQ.4) BXF(J,K) = HH
      CASE( 2) ; BYS(I,K) = HH + 0.5*DY(0)   *(DVDT+FVY(I,0,K))
      IF (GET_H.AND.IV(IW).EQ.4) BYS(I,K) = HH
      CASE(-2) ; BYF(I,K) = HH - 0.5*DY(JBP1)*(DVDT+FVY(I,JBAR,K))
      IF (GET_H.AND.IV(IW).EQ.4) BYF(I,K) = HH
      CASE( 3) ; BZS(I,J) = HH + 0.5*DZ(0)   *(DWDT+FVZ(I,J,0))
      IF (GET_H.AND.IV(IW).EQ.4) BZS(I,J) = HH
      CASE(-3) ; BZF(I,J) = HH - 0.5*DZ(KBP1)*(DWDT+FVZ(I,J,KBAR))
      IF (GET_H.AND.IV(IW).EQ.4) BZF(I,J) = HH
      END SELECT
C
      ENDIF NOT_OPEN
C
      OPEN: IF (IV(IW).EQ.2) THEN
C
      SELECT CASE(IOR)
      CASE( 1)
         U2  = UU(0,J,K)**2
         V2  = .25*(VV(1,J,K)+VV(1,J-1,K))**2
         W2  = .25*(WW(1,J,K)+WW(1,J,K-1))**2
         HQ2 = MIN(5000._EB,0.5*(U2+V2+W2))
         IF (UU(0,J,K).LT.0.) THEN
            BXS(J,K) = HQ2
            ELSE
            BXS(J,K) = H0RR(1) + HQ2*HFAC
            ENDIF
      CASE(-1)
         U2  = UU(IBAR,J,K)**2
         V2  = .25*(VV(IBAR,J,K)+VV(IBAR,J-1,K))**2
         W2  = .25*(WW(IBAR,J,K)+WW(IBAR,J,K-1))**2
         HQ2 = MIN(5000._EB,0.5*(U2+V2+W2))
         IF (UU(IBAR,J,K).GT.0.) THEN
            BXF(J,K) = HQ2
            ELSE
            BXF(J,K) = H0RR(2) + HQ2*HFAC
            ENDIF
      CASE( 2)
         U2  = .25*(UU(I,1,K)+UU(I-1,1,K))**2
         V2  = VV(I,0,K)**2
         W2  = .25*(WW(I,1,K)+WW(I,1,K-1))**2
         HQ2 = MIN(5000._EB,0.5*(U2+V2+W2))
         IF (VV(I,0,K).LT.0.) THEN
            BYS(I,K) = HQ2
            ELSE
            BYS(I,K) = H0RR(3) + HQ2*HFAC
            ENDIF
      CASE(-2)
         U2  = .25*(UU(I,JBAR,K)+UU(I-1,JBAR,K))**2
         V2  = VV(I,JBAR,K)**2
         W2  = .25*(WW(I,JBAR,K)+WW(I,JBAR,K-1))**2
         HQ2 = MIN(5000._EB,0.5*(U2+V2+W2))
         IF (VV(I,JBAR,K).GT.0.) THEN
            BYF(I,K) = HQ2
            ELSE
            BYF(I,K) = H0RR(4) + HQ2*HFAC
            ENDIF
      CASE( 3)
         U2  = .25*(UU(I,J,1)+UU(I-1,J,1))**2
         V2  = .25*(VV(I,J,1)+VV(I,J-1,1))**2
         W2  = WW(I,J,0)**2
         HQ2 = MIN(5000._EB,0.5*(U2+V2+W2))
         IF (WW(I,J,0).LT.0.) THEN
            BZS(I,J) = HQ2
            ELSE
            BZS(I,J) = H0RR(5) + HQ2*HFAC
            ENDIF
      CASE(-3)
         U2  = .25*(UU(I,J,KBAR)+UU(I-1,J,KBAR))**2
         V2  = .25*(VV(I,J,KBAR)+VV(I,J-1,KBAR))**2
         W2  = WW(I,J,KBAR)**2
         HQ2 = MIN(5000._EB,0.5*(U2+V2+W2))
         IF (WW(I,J,KBAR).GT.0.) THEN
            BZF(I,J) = HQ2
            ELSE
            BZF(I,J) = H0RR(6) + HQ2*HFAC
            ENDIF
      END SELECT
C
      ENDIF OPEN
C
      ENDIF DIRICHLET
C
      ENDDO WALL_CELL_LOOP
C
C Compute the RHS of the Poisson equation
C
      IF (CYLINDRICAL) THEN
      DO K=1,KBAR
      DO I=1,IBAR
      TRM1 = (R(I-1)*FVX(I-1,1,K)-R(I)*FVX(I,1,K))*RDX(I)*RRN(I)
      TRM3 = (FVZ(I,1,K-1)-FVZ(I,1,K))*RDZ(K)
      TRM4 = -DDDT(I,1,K)
      PRHS(I,1,K) = TRM1 + TRM3 + TRM4
      ENDDO
      ENDDO
      END IF
C
      IF ( (IPS.LE.1 .OR. IPS.EQ.4 .OR. IPS.EQ.7) .AND. 
     .      .NOT.CYLINDRICAL) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
      TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J) 
      TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
      TRM4 = -DDDT(I,J,K)
      PRHS(I,J,K) = TRM1 + TRM2 + TRM3 + TRM4
      ENDDO
      ENDDO
      ENDDO
      END IF
C
      IF (IPS.EQ.2) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
      TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J)
      TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
      TRM4 = -DDDT(I,J,K)
      PRHS(J,I,K) = TRM1 + TRM2 + TRM3 + TRM4
      ENDDO
      ENDDO
      ENDDO
      BZST = TRANSPOSE(BZS)
      BZFT = TRANSPOSE(BZF)
      END IF
C
      IF (IPS.EQ.3 .OR. IPS.EQ.6) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
      TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J)
      TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
      TRM4 = -DDDT(I,J,K)
      PRHS(K,J,I) = TRM1 + TRM2 + TRM3 + TRM4
      ENDDO
      ENDDO
      ENDDO
      BXST = TRANSPOSE(BXS)
      BXFT = TRANSPOSE(BXF)
      BYST = TRANSPOSE(BYS)
      BYFT = TRANSPOSE(BYF)
      BZST = TRANSPOSE(BZS)
      BZFT = TRANSPOSE(BZF)
      END IF
C
      IF (IPS.EQ.5) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      TRM1 = (FVX(I-1,J,K)-FVX(I,J,K))*RDX(I)
      TRM2 = (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J)
      TRM3 = (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
      TRM4 = -DDDT(I,J,K)
      PRHS(I,K,J) = TRM1 + TRM2 + TRM3 + TRM4
      ENDDO
      ENDDO
      ENDDO
      BXST = TRANSPOSE(BXS)
      BXFT = TRANSPOSE(BXF)
      END IF
C
C Poisson solve with stretching the 1st coordinate
C
      IF (IPS.LE.1) THEN
      IF (.NOT.TWO_D)
     .CALL H3CZSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,
     .            PRHS,POIS_PTB,SAVE,WORK,HX)
      IF (TWO_D .AND. .NOT. CYLINDRICAL)
     .CALL H2CZSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HX)
      IF (TWO_D .AND. CYLINDRICAL)
     .CALL H2CYSS(BXS,BXF,BZS,BZF,ITRN,PRHS,POIS_PTB,SAVE,WORK)
      ENDIF
C
      IF (IPS.EQ.2)    
     .CALL H3CZSS(BYS,BYF,BXS,BXF,BZST,BZFT,ITRN,JTRN,
     .            PRHS,POIS_PTB,SAVE,WORK,HY)
C
      IF (IPS.EQ.3) THEN
      IF (.NOT.TWO_D) 
     .CALL H3CZSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,
     .            PRHS,POIS_PTB,SAVE,WORK,HZ)
      IF (TWO_D)
     .CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HZ)
      ENDIF
C
C Poisson solve with stretching the 1st and 2nd coordinate
C
      IF (IPS.EQ.4)    
     .CALL H3CSSS(BXS,BXF,BYS,BYF,BZS,BZF,ITRN,JTRN,
     .            PRHS,POIS_PTB,SAVE,WORK,HX,HY)
C
      IF (IPS.EQ.5) THEN
      IF (.NOT.TWO_D)
     .CALL H3CSSS(BXST,BXFT,BZS,BZF,BYS,BYF,ITRN,JTRN,
     .            PRHS,POIS_PTB,SAVE,WORK,HX,HZ)
      IF (TWO_D)
     .CALL H2CZSS(BZS,BZF,BXS,BXF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HZ)
      ENDIF
C
      IF (IPS.EQ.6)    
     .CALL H3CSSS(BZST,BZFT,BYST,BYFT,BXST,BXFT,ITRN,JTRN,
     .            PRHS,POIS_PTB,SAVE,WORK,HZ,HY)
C
      IF (IPS.EQ.7)
     .CALL H2CZSS(BXS,BXF,BYS,BYF,ITRN,PRHS,POIS_PTB,SAVE,WORK,HX)
C
      IF (IPS.LE.1 .OR. IPS.EQ.4 .OR. IPS.EQ.7) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      H(I,J,K) = PRHS(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      END IF
C
      IF (IPS.EQ.2) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      H(I,J,K) = PRHS(J,I,K)
      ENDDO
      ENDDO
      ENDDO
      END IF
C
      IF (IPS.EQ.3 .OR. IPS.EQ.6) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      H(I,J,K) = PRHS(K,J,I)
      ENDDO
      ENDDO
      ENDDO
      END IF
C
      IF (IPS.EQ.5) THEN
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      H(I,J,K) = PRHS(I,K,J)
      ENDDO
      ENDDO
      ENDDO
      END IF
C
C Apply boundary conditions to H
C
      DO K=1,KBAR
      DO J=1,JBAR
      IF (LBC.EQ.3 .OR. LBC.EQ.4) 
     .    H(0,J,K)    = H(1,J,K)    - DXI*BXS(J,K)
      IF (LBC.EQ.3 .OR. LBC.EQ.2 .OR. LBC.EQ.6) 
     .    H(IBP1,J,K) = H(IBAR,J,K) + DXI*BXF(J,K)
      IF (LBC.EQ.1 .OR. LBC.EQ.2) 
     .    H(0,J,K)    =-H(1,J,K)    + 2.*BXS(J,K)
      IF (LBC.EQ.1 .OR. LBC.EQ.4 .OR. LBC.EQ.5) 
     .    H(IBP1,J,K) =-H(IBAR,J,K) + 2.*BXF(J,K)
      IF (LBC.EQ.5 .OR. LBC.EQ.6)
     .   H(0,J,K) = H(1,J,K)
      ENDDO
      ENDDO
C
      DO K=1,KBAR
      DO I=1,IBAR
      IF (MBC.EQ.3 .OR. MBC.EQ.4)
     .    H(I,0,K)    = H(I,1,K)    - DETA*BYS(I,K)
      IF (MBC.EQ.3 .OR. MBC.EQ.2)
     .    H(I,JBP1,K) = H(I,JBAR,K) + DETA*BYF(I,K)
      IF (MBC.EQ.1 .OR. MBC.EQ.2)
     .    H(I,0,K)    =-H(I,1,K)    + 2.*BYS(I,K)
      IF (MBC.EQ.1 .OR. MBC.EQ.4)
     .    H(I,JBP1,K) =-H(I,JBAR,K) + 2.*BYF(I,K)
      ENDDO
      ENDDO
C
      DO J=1,JBAR
      DO I=1,IBAR
      IF (NBC.EQ.3 .OR. NBC.EQ.4)
     .    H(I,J,0)    = H(I,J,1)    - DZETA*BZS(I,J)
      IF (NBC.EQ.3 .OR. NBC.EQ.2)
     .    H(I,J,KBP1) = H(I,J,KBAR) + DZETA*BZF(I,J)
      IF (NBC.EQ.1 .OR. NBC.EQ.2)
     .    H(I,J,0)    =-H(I,J,1)    + 2.*BZS(I,J)
      IF (NBC.EQ.1 .OR. NBC.EQ.4)
     .    H(I,J,KBP1) =-H(I,J,KBAR) + 2.*BZF(I,J)
      ENDDO
      ENDDO
C
C ************************* Check the Solution *************************
C
      IF (CHECK_POISSON) THEN     
      POIS_ERR = 0.
      DO K=1,KBAR
      DO J=1,JBAR
      DO I=1,IBAR
      RHSS = (R(I-1)*FVX(I-1,J,K)-R(I)*FVX(I,J,K))*RDX(I)*RRN(I) +
     .       (FVY(I,J-1,K)-FVY(I,J,K))*RDY(J) +
     .       (FVZ(I,J,K-1)-FVZ(I,J,K))*RDZ(K)
     .     - DDDT(I,J,K)
      LHSS = ((H(I+1,J,K)-H(I,J,K))*RDXN(I)*R(I) -
     .        (H(I,J,K)-H(I-1,J,K))*RDXN(I-1)*R(I-1) )*RDX(I)*RRN(I)
     .      +((H(I,J+1,K)-H(I,J,K))*RDYN(J) -
     .        (H(I,J,K)-H(I,J-1,K))*RDYN(J-1) )*RDY(J)
     .      +((H(I,J,K+1)-H(I,J,K))*RDZN(K) -
     .        (H(I,J,K)-H(I,J,K-1))*RDZN(K-1) )*RDZ(K)
      RES = ABS(RHSS-LHSS)
      POIS_ERR = MAX(RES,POIS_ERR)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
C
C **********************************************************************
C
      TUSED(5,NM)=TUSED(5,NM)+SECOND()-TNOW
      END SUBROUTINE PRESSURE_SOLVER
C
C
      END MODULE PRES
