      MODULE PYRO
C
      USE PREC
      USE CONS, ONLY : NMESHES,IOXYGEN,MHC,IFUEL,R1,SHUTDOWN,
     .                 LU5,RADIATION,SIGMA,TMPA,CHECKREAD
      USE VARS, ONLY : MESH
C
      IMPLICIT NONE
      PRIVATE
      PUBLIC INIT_SOLID_PHASE_REACTIONS,SOLID_PHASE_REACTIONS
C
      REAL(EB) A_P,A_OX,A_CHAR,E_P,E_OX,E_CHAR,N_P,
     .         N_OX,N_F_OX,N_O2_CHAR,N_CHAR,NU_CHAR_P,NU_CHAR_OX,
     .         NU_ASH_CHAR,QH_P,QH_OX,QH_CHAR,C_S_0,C_S,DELTA,RHO_S_0,
     .         NU_SF_P,NU_SF_OX,NU_SF_CHAR,NU_O2_OX,NU_O2_CHAR,
     .         QMAX,QEXT_ON,QEXT_OFF,XRAD,YRAD,ZRAD,RMAX
      INTEGER IBC,IW
C
      TYPE PYRO_TYPE
      REAL(EB), POINTER, DIMENSION(:) :: RHO_S,Y_SF,Y_CHAR,
     .         RHO_S_ORIG,Y_SF_ORIG,Y_CHAR_ORIG
      END TYPE PYRO_TYPE
C
      TYPE (PYRO_TYPE), ALLOCATABLE, TARGET, DIMENSION(:) :: PMESH
C
C
      CONTAINS
C
C
      SUBROUTINE INIT_SOLID_PHASE_REACTIONS
C
      INTEGER IOS,NM
      NAMELIST /PYRO/ A_P,A_OX,A_CHAR,E_P,E_OX,E_CHAR,N_P,
     .         N_OX,N_F_OX,N_O2_CHAR,N_CHAR,NU_CHAR_P,NU_CHAR_OX,
     .         NU_ASH_CHAR,QH_P,QH_OX,QH_CHAR,C_S_0,DELTA,RHO_S_0,
     .         NU_SF_P,NU_SF_OX,NU_SF_CHAR,NU_O2_OX,NU_O2_CHAR,
     .         QMAX,QEXT_ON,QEXT_OFF,XRAD,YRAD,ZRAD,RMAX
C
      ALLOCATE(PMESH(NMESHES))
C
      DO NM=1,NMESHES
      ALLOCATE(PMESH(NM)%RHO_S(MESH(NM)%NDWC))
      ALLOCATE(PMESH(NM)%Y_SF(MESH(NM)%NDWC))
      ALLOCATE(PMESH(NM)%Y_CHAR(MESH(NM)%NDWC))
      ALLOCATE(PMESH(NM)%RHO_S_ORIG(MESH(NM)%NDWC))
      ALLOCATE(PMESH(NM)%Y_SF_ORIG(MESH(NM)%NDWC))
      ALLOCATE(PMESH(NM)%Y_CHAR_ORIG(MESH(NM)%NDWC))
      ENDDO
C
C Whatman ashless filter 44, see Kashiwagi, 26th Symposium, pp 1345
C
      A_P          = 3.0E18          ! 1/min
      A_OX         = 5.0E18          ! 1/min
      A_CHAR       = 1.0E11          ! 1/min
      E_P          = 56.6            ! kcal/mol
      E_OX         = 53.6            ! kcal/mol
      E_CHAR       = 39.7            ! kcal/mol
      N_P          = 1.20
      N_OX         = 0.40
      N_F_OX       = 1.00
      N_O2_CHAR    = 1.00
      N_CHAR       = 1.00
      NU_CHAR_P    = 0.14
      NU_CHAR_OX   = 0.20
      NU_ASH_CHAR  = 0.02
      NU_SF_P      = 0.23
      NU_SF_OX     = 0.16
      NU_SF_CHAR   = 0.52
      NU_O2_OX     = 0.41
      NU_O2_CHAR   = 1.65
      QH_P         =  -64.           ! kJ/kg
      QH_OX        =-1000.           ! kJ/kg
      QH_CHAR      =-2400.           ! kJ/kg
      C_S_0        = 0.96            ! kJ/kg/K
      DELTA        = 0.000065        ! m
      RHO_S_0      = 438.5           ! kg/m^3
C
C External radiative flux
C
      QMAX         = 50.             ! kW/m^2
      QEXT_ON      = 0.              ! s
      QEXT_OFF     = 1000000.        ! s
      RMAX         = 0.0025          ! m
      XRAD         = 0.              ! m
      YRAD         = 0.              ! m
      ZRAD         = 0.              ! m
C
C Read in parameters from CHID.data file, NAMELIST: PYRO
C
      REWIND(LU5)
      READ_PYRO_LOOP: DO
      CALL CHECKREAD('PYRO',LU5,IOS) ; IF (IOS.EQ.1) EXIT READ_PYRO_LOOP
      READ(LU5,NML=PYRO,END=11,ERR=12,IOSTAT=IOS)
   12 IF (IOS.GT.0) THEN
         CALL SHUTDOWN(' ERROR: Problem with PYRO NAMELIST line')
         ENDIF
      ENDDO READ_PYRO_LOOP
   11 REWIND(LU5)
C
C Convert input parameters into appropriate units
C
      QMAX    = QMAX*1000.         !  kW/m^2  --> W/m^2
      QH_P    = QH_P*1000.         !  kJ/kg   --> J/kg
      QH_OX   = QH_OX*1000.        !  kJ/kg   --> J/kg
      QH_CHAR = QH_CHAR*1000.      !  kJ/kg   --> J/kg
      C_S_0   = C_S_0*1000.        !  kJ/kg/K --> J/kg/K
      A_P     = A_P/60.            !  1/min   --> 1/s
      A_OX    = A_OX/60.           !  1/min   --> 1/s
      A_CHAR  = A_CHAR/60.         !  1/min   --> 1/s
C
C Set initial values for solid phase parameters
C
      MESH_LOOP: DO NM=1,NMESHES
C
      WALL_CELL_LOOP: DO IW=1,MESH(NM)%NWC
      IBC = MESH(NM)%IJKW(5,IW)
      IF (MHC(IBC).NE.4 .OR. MESH(NM)%IV(IW).EQ.0) CYCLE WALL_CELL_LOOP
      PMESH(NM)%RHO_S(IW)  = RHO_S_0
      PMESH(NM)%Y_SF(IW)   = 1.0 
      PMESH(NM)%Y_CHAR(IW) = 0.0 
      ENDDO WALL_CELL_LOOP
C
      ENDDO MESH_LOOP
C
      END SUBROUTINE INIT_SOLID_PHASE_REACTIONS
C
C
      SUBROUTINE SOLID_PHASE_REACTIONS(T,NM)
C
      USE PACKER, ONLY : UNPACK_VAR,TMP_F,
     .                 YY,IJKW,NWC,DT,MASSFLUX,
     .                 QRAD,QCONF,XW,YW,ZW,KW,
     .                 TMP,RDN,RHO,
     .                 RHODW,TMP_W,YY_W,RHO_W,IV,NEW_TIME_STEP
      USE CONS, ONLY : TMPMIN,TMPMAX
C
      REAL(EB) DV(4),RHS(4),T
      REAL(EB) K_P,K_OX,K_CHAR,Y_O2,TMPS,DT_S,
     .         QEXT,DTMP,QREAC,R2,MDOT,MDOT_SF,MDOT_O2,
     .         RHO_GAS,RHO_GHOST,Y_O2_GHOST,Y_O2_GAS,UN,FLUX,
     .         FLUX_O2,FLUX_SF,RHOD,DENOM,R_RHO
      INTEGER II,JJ,KK,IIG,JJG,KKG,NITER,ITER,IOR
      INTEGER, INTENT(IN) :: NM
C
      CALL UNPACK_VAR(NM)
C
      IF (.NOT.NEW_TIME_STEP) THEN
         PMESH(NM)%RHO_S_ORIG  = PMESH(NM)%RHO_S
         PMESH(NM)%Y_SF_ORIG   = PMESH(NM)%Y_SF
         PMESH(NM)%Y_CHAR_ORIG = PMESH(NM)%Y_CHAR
         ENDIF
C
      WALL_CELL_LOOP: DO IW=1,NWC
C
      IBC = IJKW(5,IW)
      IF (MHC(IBC).NE.4 .OR. IV(IW).NE.1) CYCLE WALL_CELL_LOOP
C
      II  = IJKW(1,IW)
      JJ  = IJKW(2,IW)
      KK  = IJKW(3,IW)
      IIG = IJKW(6,IW)
      JJG = IJKW(7,IW)
      KKG = IJKW(8,IW)
      IOR = IJKW(4,IW)
C
      DV(1) = PMESH(NM)%RHO_S_ORIG(IW)/RHO_S_0
      DV(2) = DV(1)*PMESH(NM)%Y_SF_ORIG(IW)
      DV(3) = DV(1)*PMESH(NM)%Y_CHAR_ORIG(IW)
      DV(4) = TMP_F(IW)
C
      IF (T.GE.QEXT_ON .AND. T.LT.QEXT_OFF) THEN
         SELECT CASE(ABS(IOR))
         CASE(1) ; R2 = (YW(IW)-YRAD)**2 + (ZW(IW)-ZRAD)**2
         CASE(2) ; R2 = (XW(IW)-XRAD)**2 + (ZW(IW)-ZRAD)**2
         CASE(3) ; R2 = (XW(IW)-XRAD)**2 + (YW(IW)-YRAD)**2
         END SELECT
         QEXT = QMAX*EXP(-R2/RMAX**2)
         ELSE
         QEXT = 0.
         ENDIF
C
      RHO_GHOST  = RHO_W(IW)
      RHO_GAS    = RHO(IIG,JJG,KKG)
      R_RHO      = 2./(RHO_GHOST+RHO_GAS)
      Y_O2_GHOST = YY_W(IW,IOXYGEN)
      Y_O2_GAS   = YY(IIG,JJG,KKG,IOXYGEN)
      Y_O2  = 0.5*(Y_O2_GHOST+Y_O2_GAS)
      Y_O2  = MAX(0.0_EB,Y_O2)
      RHOD  = RHODW(IW,IOXYGEN)*RDN(IW)
C
      NITER = 50
      DT_S  = DT/REAL(NITER)
      MDOT    = 0.
      MDOT_SF = 0.
      MDOT_O2 = 0.
C
      ODE_SOLVER: DO ITER=1,NITER
C
      TMPS   = DV(4)
      C_S    = C_S_0 + 4.19*(DV(4)-TMPA)
C     C_S    = C_S_0
      K_P    = A_P   *EXP(-E_P   /(R1*TMPS))*DV(2)**N_P
      K_OX   = A_OX  *EXP(-E_OX  /(R1*TMPS))*DV(2)**N_F_OX*Y_O2**N_OX
      K_CHAR = A_CHAR*EXP(-E_CHAR/(R1*TMPS))*DV(3)**N_CHAR*
     .         Y_O2**N_O2_CHAR
      QREAC = (-QH_P*K_P - QH_OX*K_OX - QH_CHAR*K_CHAR)*DELTA*RHO_S_0
      DTMP  = TMP(IIG,JJG,KKG) - TMPS
      QCONF(IW) = 2.*KW(IW)*DTMP*RDN(IW)
C
      RHS(1) = (NU_CHAR_P-1.)  *K_P + 
     .         (NU_CHAR_OX-1.) *K_OX + 
     .         (NU_ASH_CHAR-1.)*K_CHAR 
      RHS(2) = -K_P - K_OX
      RHS(3) = NU_CHAR_P*K_P + NU_CHAR_OX*K_OX - K_CHAR
      IF (.NOT.RADIATION) QRAD(IW) = -SIGMA*(TMPS**4-TMPA**4)
      RHS(4) = (QREAC + DV(1)*QEXT + QRAD(IW) + QCONF(IW))/
     .         (DELTA*C_S*RHO_S_0*DV(1))
C
      DV(1:4) = DV(1:4) + DT_S*RHS(1:4)
C
      DV(1)   = MAX(0.01_EB,DV(1))
      DV(2)   = MAX(0.00_EB,DV(2))
      DV(3)   = MAX(0.00_EB,DV(3))
      DV(4)   = MIN(TMPMAX,MAX(TMPMIN,DV(4)))
C
      FLUX    = -RHS(1)*DELTA*RHO_S_0
      FLUX_SF = (NU_SF_P*K_P+NU_SF_OX*K_OX+NU_SF_CHAR*K_CHAR)*
     .          DELTA*RHO_S_0
      FLUX_O2 = (-NU_O2_OX*K_OX-NU_O2_CHAR*K_CHAR)*DELTA*RHO_S_0
C
      MDOT    = MDOT    + FLUX*DT_S
      MDOT_SF = MDOT_SF + FLUX_SF*DT_S
      MDOT_O2 = MDOT_O2 + FLUX_O2*DT_S
C
      UN    = FLUX*R_RHO
      DENOM = RHOD + 0.5*UN*RHO_GHOST
      Y_O2_GHOST = ( FLUX_O2 + Y_O2_GAS*(RHOD-0.5*UN*RHO_GAS) )/DENOM
      Y_O2  = 0.5*(Y_O2_GHOST+Y_O2_GAS)
      Y_O2  = MAX(0.0_EB,Y_O2)
C
      ENDDO ODE_SOLVER
C
C Set surface temperature (TMP_F) and ghost cell temperature (TMP_W)
C
      TMP_F(IW) = TMPS
      TMP_W(IW) = TMP(IIG,JJG,KKG) - QCONF(IW)/(RDN(IW)*KW(IW))
C
C Set mass fluxes at the surface
C
      MASSFLUX(IW,IFUEL)   = MDOT_SF/DT
      MASSFLUX(IW,IOXYGEN) = MDOT_O2/DT
      MASSFLUX(IW,0)       = (MDOT-MDOT_SF-MDOT_O2)/DT
C
C Update values of RHO_S, Y_SF, and Y_CHAR
C
      PMESH(NM)%RHO_S(IW) = DV(1)*RHO_S_0
C
      PMESH(NM)%Y_SF(IW)   = DV(2)/DV(1)
      PMESH(NM)%Y_CHAR(IW) = DV(3)/DV(1)
C
      ENDDO WALL_CELL_LOOP
C
      END SUBROUTINE SOLID_PHASE_REACTIONS
C 
      END MODULE PYRO
