!> \brief Variables that are defined on each mesh.

MODULE MESH_VARIABLES

USE PRECISION_PARAMETERS
USE TYPES
IMPLICIT NONE

!> \brief Derived type containing the bulk of the variables defined on each mesh.
!>
!> \details The variables listed in this module are defined on each mesh that is
!> owned by a given MPI process. For example, MESHES(NM)%U(I,J,K) is the \f$u\f$
!> component of velocity at the forward \f$x\f$ face of grid cell (I,J,K) and
!> mesh NM.

TYPE MESH_TYPE

   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: U  !< Velocity component at current time step, \f$u_{ijk}^n\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: V  !< Velocity component at current time step, \f$v_{ijk}^n\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: W  !< Velocity component at current time step, \f$w_{ijk}^n\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: US !< Velocity component estimated at next time step, \f$u_{ijk}^*\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: VS !< Velocity component estimated at next time step, \f$v_{ijk}^*\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WS !< Velocity component estimated at next time step, \f$w_{ijk}^*\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: DDDT    !< \f$(\partial D/\partial t)_{ijk}\f$ where \f$D_{ijk}\f$ is the divergence
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: D       !< Divergence at current time step, \f$D_{ijk}^n\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: DS      !< Divergence estimate next time step, \f$D_{ijk}^*\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: H       !< \f$ \tilde{p}_{ijk}/\rho_{ijk} + |\mathbf{u}|^2_{ijk}/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: HS      !< H estimated at next time step
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: H_PRIME !< Experimental pressure correction
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: KRES    !< Resolved kinetic energy, \f$ |\mathbf{u}|^2_{ijk}/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVX     !< Momentum equation flux terms, \f$ F_{{\rm A},x,ijk}+F_{{\rm B},x,ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVY     !< Momentum equation flux terms, \f$ F_{{\rm A},y,ijk}+F_{{\rm B},y,ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVZ     !< Momentum equation flux terms, \f$ F_{{\rm A},z,ijk}+F_{{\rm B},z,ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVX_B   !< Momentum equation flux terms, \f$ F_{{\rm B},x,ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVY_B   !< Momentum equation flux terms, \f$ F_{{\rm B},y,ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FVZ_B   !< Momentum equation flux terms, \f$ F_{{\rm B},z,ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: RHO     !< Density (kg/m3) at current time step, \f$ \rho_{ijk}^n \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: RHOS    !< Density (kg/m3) at next time step, \f$ \rho_{ijk}^* \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU      !< Turbulent viscosity (kg/m/s), \f$ \mu_{{\rm t},ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MU_DNS  !< Laminar viscosity (kg/m/s)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: TMP     !< Gas temperature, \f$ T_{ijk} \f$ (K)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: Q       !< Heat release rate per unit volume, \f$ \dot{q}_{ijk}''' \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: KAPPA_GAS !< Radiation absorption coefficient by gas, \f$ \kappa_{ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: CHI_R   !< Radiative fraction, \f$ \chi_{{\rm r},ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: QR      !< Radiation source term, \f$ -\nabla \cdot \dot{\mathbf{q}}_{\rm r}'' \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: QR_W    !< Radiation source term, particles and droplets
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: UII     !< Integrated intensity, \f$ U_{ijk}=\sum_{l=1}^N I_{ijk}^l\delta\Omega^l\f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: RSUM    !< \f$ R_0 \sum_\alpha Z_{\alpha,ijk}/W_\alpha \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: D_SOURCE!< Source terms in the expression for the divergence
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: U_OLD   !< Value of \f$ u_{ijk} \f$ at the previous time step, used for output only
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: V_OLD   !< Value of \f$ v_{ijk} \f$ at the previous time step, used for output only
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: W_OLD   !< Value of \f$ w_{ijk} \f$ at the previous time step, used for output only
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: CSD2    !< \f$ C_s \Delta^2 \f$ in Smagorinsky turbulence expression
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: CHEM_SUBIT  !< Number of chemistry sub-iterations
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: MIX_TIME    !< Mixing-controlled combustion reaction time (s)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: STRAIN_RATE !< Strain rate \f$ |S|_{ijk} \f$ (1/s)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: D_Z_MAX     !< \f$ \max D_\alpha \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: Q_DOT_PPP_S !< Heat release rate per unit volume in 3D pyrolysis model
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: PR_T        !< Turbulent Prandtl number (experimental)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: TMP_FLAME   !< Flame temperature (K) (experimental)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: FLAME_INDEX

   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: ZZ               !< Lumped species, current time step, \f$ Z_{\alpha,ijk}^n \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: ZZS              !< Lumped species, next time step, \f$ Z_{\alpha,ijk}^* \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: REAC_SOURCE_TERM !< \f$ \dot{m}_{\alpha,ijk}''' \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: DEL_RHO_D_DEL_Z  !< \f$ (\nabla \cdot \rho D_\alpha \nabla Z_\alpha)_{ijk} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: FX               !< \f$ \rho Z_{\alpha,ijk} \f$ at \f$ x \f$ face of cell
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: FY               !< \f$ \rho Z_{\alpha,ijk} \f$ at \f$ y \f$ face of cell
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: FZ               !< \f$ \rho Z_{\alpha,ijk} \f$ at \f$ z \f$ face of cell
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: Q_REAC           !< \f$ \dot{q}_{ijk}''' \f$ for a specified reaction
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: AVG_DROP_DEN     !< Droplet mass per unit volume for a certain droplet type
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: AVG_DROP_TMP     !< Average temperature for a certain droplet type
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: AVG_DROP_RAD     !< Average radius for a certain droplet type
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: AVG_DROP_AREA    !< Average area for a certain droplet type
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: M_DOT_PPP        !< Mass source term, \f$ \dot{m}_{\alpha,ijk}''' \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: M_DOT_G_PPP_S    !< Mass source term, \f$ \dot{m}_{\alpha,ijk}''' \f$, 3D solid
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: RHO_ZZ_G_S
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: TRI_COR
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: ADV_FX
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: ADV_FY
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: ADV_FZ
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: DIF_FX
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: DIF_FY
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: DIF_FZ
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: DIF_FXS
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: DIF_FYS
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: DIF_FZS

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: U_EDGE_Y,U_EDGE_Z,V_EDGE_X,V_EDGE_Z,W_EDGE_X,W_EDGE_Y

   REAL(EB) :: POIS_PTB,POIS_ERR,LAPLACE_PTB,LAPLACE_ERR
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: SAVE1,SAVE2,WORK
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: PRHS !< Right hand side of Poisson (pressure) equation
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: BXS,BXF,BYS,BYF,BZS,BZF, BXST,BXFT,BYST,BYFT,BZST,BZFT
   INTEGER :: LSAVE,LWORK,LBC,MBC,NBC,LBC2,MBC2,NBC2,ITRN,JTRN,KTRN,IPS

   REAL(EB), ALLOCATABLE, DIMENSION(:) :: P_0         !< Ambient pressure profile, \f$ \overline{p}_0(z) \f$ (Pa)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RHO_0       !< Ambient density profile, \f$ \overline{\rho}_0(z) \f$ (kg/m\f$^3\f$)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: TMP_0       !< Ambient temperature profile, \f$ \overline{T}_0(z) \f$ (K)
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: D_PBAR_DT   !< \f$ (\partial \overline{p}_m/\partial t)^n \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: D_PBAR_DT_S !< \f$ (\partial \overline{p}_m/\partial t)^* \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: U_LEAK
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: U_DUCT

   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: PBAR      !< Background pressure, current, \f$ \overline{p}_m^n(z,t) \f$ (Pa)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: PBAR_S    !< Background pressure, estimated, \f$ \overline{p}_m^*(z,t) \f$ (Pa)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: R_PBAR    !< \f$ 1/\overline{p}_m(z,t) \f$
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: PRESSURE_ZONE !< Index of the pressure zone for cell (I,J,K)
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: DCOR !< Divergence correction term, \f$ \nabla \cdot (\bar{\mathbf{u}}-\mathbf{u}) \f$

   !! Laplace solve, sparse LU

   !INTEGER :: A_N_ROWS,A_N_COLS,A_N_ELEMENTS
   !INTEGER, ALLOCATABLE, DIMENSION(:) :: A_ROW_INDEX,A_COLUMNS
   !REAL(EB), ALLOCATABLE, DIMENSION(:) :: A_VALUES

   ! Work arrays

   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: SCALAR_WORK1
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: SCALAR_WORK2
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: SCALAR_WORK3
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: SCALAR_WORK4
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: WORK1,WORK2,WORK3,WORK4,WORK5,WORK6,WORK7,WORK8,WORK9
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IWORK1
   REAL(EB),     ALLOCATABLE, DIMENSION(:,:,:) :: PWORK1,PWORK2,PWORK3,PWORK4
   COMPLEX(EB),  ALLOCATABLE, DIMENSION(:,:,:) :: PWORK5,PWORK6,PWORK7,PWORK8
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: TURB_WORK1,TURB_WORK2,TURB_WORK3,TURB_WORK4
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: TURB_WORK5,TURB_WORK6,TURB_WORK7,TURB_WORK8
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: TURB_WORK9,TURB_WORK10
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:) :: IBM_SAVE1,IBM_SAVE2,IBM_SAVE3,IBM_SAVE4,IBM_SAVE5,IBM_SAVE6
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: WALL_WORK1,WALL_WORK2,FACE_WORK1,FACE_WORK2,FACE_WORK3
   REAL(FB), ALLOCATABLE, DIMENSION(:,:,:,:) :: QQ, QQ2
   REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: PP,PPN,BNDF_TIME_INTEGRAL
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IBK
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IBLK

   REAL(EB) :: CFL,DIVMX,DIVMN,VN,RESMAX,PART_UVWMAX=0._EB
   INTEGER  :: ICFL,JCFL,KCFL,IMX,JMX,KMX,IMN,JMN,KMN, I_VN,J_VN,K_VN,IRM,JRM,KRM
   INTEGER  :: IUVW=1
   LOGICAL, DIMENSION(-3:3) :: APPLY_SPONGE_LAYER=.FALSE.
   LOGICAL :: BAROCLINIC_TERMS_ATTACHED=.FALSE.
   CHARACTER(LABEL_LENGTH) :: TRNX_ID,TRNY_ID,TRNZ_ID

   INTEGER :: N_EDGES=0
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IJKE,EDGE_INDEX
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: TAU_E,OME_E

   INTEGER :: MESH_LEVEL,LBC_EMB,MBC_EMB,NBC_EMB
   INTEGER :: IBAR !< Number of cells in the \f$ x \f$ direction, \f$ I \f$
   INTEGER :: JBAR !< Number of cells in the \f$ y \f$ direction, \f$ J \f$
   INTEGER :: KBAR !< Number of cells in the \f$ z \f$ direction, \f$ K \f$
   INTEGER :: IBM1 !< IBAR minus 1
   INTEGER :: JBM1 !< JBAR minus 1
   INTEGER :: KBM1 !< KBAR minus 1
   INTEGER :: IBP1 !< IBAR plus 1
   INTEGER :: JBP1 !< JBAR plus 1
   INTEGER :: KBP1 !< KBAR plus 1
   INTEGER :: N_NEIGHBORING_MESHES !< Number of meshing abutting the current one
   INTEGER, ALLOCATABLE, DIMENSION(:) :: NEIGHBORING_MESH  !< Array listing the indices of neighboring meshes
   INTEGER, ALLOCATABLE, DIMENSION(:) :: RGB               !< Color indices of the mesh for Smokeview

   ! Mesh coordinate variables

   REAL(EB) :: DXI                               !< \f$ \delta \xi = (x_I-x_0)/I \f$
   REAL(EB) :: DETA                              !< \f$ \delta \eta = (y_J-y_0)/J \f$
   REAL(EB) :: DZETA                             !< \f$ \delta \zeta = (z_K-z_0)/K \f$
   REAL(EB) :: RDXI                              !< \f$ 1/ \delta \xi \f$
   REAL(EB) :: RDETA                             !< \f$ 1/ \delta \eta \f$
   REAL(EB) :: RDZETA                            !< \f$ 1/ \delta \zeta \f$
   REAL(EB) :: DXMIN                             !< \f$ \min_i \delta x_i \f$
   REAL(EB) :: DXMAX                             !< \f$ \max_i \delta x_i \f$
   REAL(EB) :: DYMIN                             !< \f$ \min_j \delta y_j \f$
   REAL(EB) :: DYMAX                             !< \f$ \max_j \delta y_j \f$
   REAL(EB) :: DZMIN                             !< \f$ \min_k \delta z_k \f$
   REAL(EB) :: DZMAX                             !< \f$ \max_k \delta z_k \f$
   REAL(EB) :: XS                                !< Lower extent of mesh x coordinate, \f$ x_0 \f$
   REAL(EB) :: XF                                !< Upper extent of mesh x coordinate, \f$ x_I \f$
   REAL(EB) :: YS                                !< Lower extent of mesh y coordinate, \f$ y_0 \f$
   REAL(EB) :: YF                                !< Upper extent of mesh y coordinate, \f$ y_J \f$
   REAL(EB) :: ZS                                !< Lower extent of mesh z coordinate, \f$ z_0 \f$
   REAL(EB) :: ZF                                !< Upper extent of mesh z coordinate, \f$ z_K \f$
   REAL(EB) :: RDXINT                            !< \f$ 500/\delta \xi \f$
   REAL(EB) :: RDYINT                            !< \f$ 500/\delta \eta \f$
   REAL(EB) :: RDZINT                            !< \f$ 500/\delta \zeta \f$
   REAL(EB) :: CELL_SIZE                         !< Approximate cell size, \f$ (\delta\xi\,\delta\eta\,\delta\zeta)^{1/3} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: R      !< Radial coordinate, \f$ r_i \f$, for CYLINDRICAL geometry
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RC     !< Radial coordinate, cell center, \f$ (r_i+r_{i-1})/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RRN    !< \f$ 2/(r_i+r_{i-1}) \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: X      !< Position of forward x face of cell (I,J,K), \f$ x_i \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: Y      !< Position of forward y face of cell (I,J,K), \f$ y_j \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: Z      !< Position of forward z face of cell (I,J,K), \f$ z_k \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: XC     !< x coordinate of cell center, \f$ (x_i+x_{i-1})/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: YC     !< y coordinate of cell center, \f$ (y_j+y_{j-1})/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: ZC     !< z coordinate of cell center, \f$ (z_k+z_{k-1})/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: HX     !< Grid stretch factor, \f$ (x_i-x_{i-1})/\delta \xi \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: HY     !< Grid stretch factor, \f$ (y_j-y_{j-1})/\delta \eta \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: HZ     !< Grid stretch factor, \f$ (z_k-z_{k-1})/\delta \zeta \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DX     !< \f$ \delta x_i = x_i-x_{i-1} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DY     !< \f$ \delta y_j = y_j-y_{j-1} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DZ     !< \f$ \delta z_k = z_k-z_{k-1} \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDX    !< \f$ 1/\delta x_i \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDY    !< \f$ 1/\delta y_j \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDZ    !< \f$ 1/\delta z_k \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DXN    !< \f$ (x_i+x_{i+1})/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DYN    !< \f$ (y_j+y_{j+1})/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: DZN    !< \f$ (z_k+z_{k+1})/2 \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDXN   !< \f$ 2/(x_i+x_{i+1}) \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDYN   !< \f$ 2/(y_j+y_{j+1}) \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: RDZN   !< \f$ 2/(z_k+z_{k+1}) \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CELLSI !< Array used to locate the cell index of \f$ x \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CELLSJ !< Array used to locate the cell index of \f$ y \f$
   REAL(EB), ALLOCATABLE, DIMENSION(:) :: CELLSK !< Array used to locate the cell index of \f$ z \f$
   INTEGER :: CELLSI_LO,CELLSI_HI                !< hold CELLSI array bounds
   INTEGER :: CELLSJ_LO,CELLSJ_HI                !< hold CELLSJ array bounds
   INTEGER :: CELLSK_LO,CELLSK_HI                !< hold CELLSK array bounds
   REAL(FB), ALLOCATABLE, DIMENSION(:) :: XPLT   !< 4 byte real array holding \f$ x \f$ mesh coordinates
   REAL(FB), ALLOCATABLE, DIMENSION(:) :: YPLT   !< 4 byte real array holding \f$ y \f$ mesh coordinates
   REAL(FB), ALLOCATABLE, DIMENSION(:) :: ZPLT   !< 4 byte real array holding \f$ z \f$ mesh coordinates

   INTEGER :: N_OBST=0                                              !< Number of obstructions in the mesh
   TYPE(OBSTRUCTION_TYPE), ALLOCATABLE, DIMENSION(:) :: OBSTRUCTION !< Derived type variable holding obstruction information

   INTEGER :: N_VENT=0                                              !< Number of vents in the mesh
   TYPE(VENTS_TYPE), ALLOCATABLE, DIMENSION(:) :: VENTS             !< Derived type variable holding vent information

   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: CELL_INDEX             !< Unique integer identifier for grid cell (I,J,K)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: I_CELL                     !< I index of cell with identifier CELL_INDEX(I,J,K)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: J_CELL                     !< J index of cell with identifier CELL_INDEX(I,J,K)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: K_CELL                     !< K index of cell with identifier CELL_INDEX(I,J,K)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: OBST_INDEX_C               !< Index of obstruction occupying cell with CELL_INDEX(I,J,K)

   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: WALL_INDEX               !< Wall index of 6 faces of cell with CELL_INDEX(I,J,K)
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: WALL_INDEX_HT3D
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: SOLID                      !< T or F if cell with CELL_INDEX(I,J,K) is solid
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: EXTERIOR                   !< T or F if cell with CELL_INDEX(I,J,K) is outside mesh
   LOGICAL, ALLOCATABLE, DIMENSION(:) :: CONNECTED_MESH             !< T or F if cell is within another mesh
   LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: MEAN_FORCING_CELL      !< T or F if cell undergoes mean forcing wind
   INTEGER, ALLOCATABLE, DIMENSION(:) :: K_MEAN_FORCING             !< Vertical index of master mean forcing wind profile

   INTEGER :: NREGFACE_H(MAX_DIM)
   TYPE(IBM_REGFACE_TYPE), ALLOCATABLE, DIMENSION(:) :: REGFACE_IAXIS_H, &
                                                        REGFACE_JAXIS_H, &
                                                        REGFACE_KAXIS_H

   !---------------------- STR: CC_IBM mesh Arrays ------------------------------------------
   INTEGER :: N_EDGE_CROSS=0,  N_CUTEDGE_MESH=0, N_CUTFACE_MESH=0, N_CUTCELL_MESH=0
   INTEGER :: N_BBCUTFACE_MESH=0, N_GCCUTFACE_MESH=0, N_GCCUTCELL_MESH=0
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: VERTVAR, CCVAR
   INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: ECVAR, FCVAR
   TYPE(IBM_EDGECROSS_TYPE), ALLOCATABLE, DIMENSION(:) :: EDGE_CROSS
   TYPE(IBM_CUTEDGE_TYPE),   ALLOCATABLE, DIMENSION(:) :: CUT_EDGE
   TYPE(IBM_CUTFACE_TYPE),   ALLOCATABLE, DIMENSION(:) :: CUT_FACE
   TYPE(IBM_CUTCELL_TYPE),   ALLOCATABLE, DIMENSION(:) :: CUT_CELL

   INTEGER :: IBM_NREGFACE_Z(MAX_DIM)=0, IBM_NBBREGFACE_Z(MAX_DIM)=0
   TYPE(IBM_REGFACEZ_TYPE), ALLOCATABLE, DIMENSION(:) :: IBM_REGFACE_IAXIS_Z, &
                                                         IBM_REGFACE_JAXIS_Z, &
                                                         IBM_REGFACE_KAXIS_Z
   INTEGER :: IBM_NRCFACE_Z=0, IBM_NBBRCFACE_Z=0, IBM_NRCFACE_H=0
   TYPE(IBM_RCFACE_TYPE), ALLOCATABLE, DIMENSION(:) :: IBM_RCFACE_H
   TYPE(IBM_RCFACE_LST_TYPE), ALLOCATABLE, DIMENSION(:) :: IBM_RCFACE_Z
   INTEGER :: IBM_NEXIMFACE_MESH=0, IBM_NBBEXIMFACE_MESH=0
   TYPE(IBM_EXIMFACE_TYPE), ALLOCATABLE, DIMENSION(:) :: IBM_EXIM_FACE
   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: RHO_ZZN

   ! CFACE to be used in conjunction with solid side solvers:
   TYPE(CFACE_TYPE), ALLOCATABLE, DIMENSION(:) :: CFACE

   ! Edges connected to regular Gas faces to receive wall model shear stress and vorticity.
   INTEGER :: IBM_NRCEDGE=0, IBM_NIBEDGE=0
   TYPE(IBM_EDGE_TYPE), ALLOCATABLE, DIMENSION(:) :: IBM_RCEDGE, IBM_IBEDGE

   ! RAD_CFACE container that has CFACES assigned to Cartesian faces:
   INTEGER :: N_RAD_CFACE_CELLS_DIM
   TYPE(RAD_CFACE_TYPE), ALLOCATABLE, DIMENSION(:) :: RAD_CFACE

   ! Array with maximum height (Z) of geometry intersections with vertical grid lines in the mesh.
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: GEOM_ZMAX

   ! Arrays for special cut-cells:
   INTEGER :: N_SPCELL=0, N_SPCELL_CF=0
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: SPCELL_LIST

   !---------------------- END: CC_IBM mesh Arrays ------------------------------------------

   INTEGER :: N_WALL_CELLS,N_WALL_CELLS_DIM,N_INTERNAL_WALL_CELLS,N_EXTERNAL_WALL_CELLS,WALL_COUNTER,WALL_COUNTER_HT3D
   INTEGER :: N_CFACE_CELLS=0,N_CFACE_CELLS_DIM
   REAL(EB) :: BC_CLOCK,BC_CLOCK_HT3D
   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: EDGE_INTERPOLATION_FACTOR
   REAL(EB), ALLOCATABLE, DIMENSION(:)   :: D_CORR,DS_CORR,UVW_SAVE,U_GHOST,V_GHOST,W_GHOST
   TYPE(WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: WALL
   TYPE(EXTERNAL_WALL_TYPE), ALLOCATABLE, DIMENSION(:) :: EXTERNAL_WALL
   TYPE(OMESH_TYPE), ALLOCATABLE, DIMENSION(:) :: OMESH
   TYPE(LAGRANGIAN_PARTICLE_TYPE), ALLOCATABLE, DIMENSION(:) :: LAGRANGIAN_PARTICLE
   TYPE(STORAGE_TYPE), ALLOCATABLE, DIMENSION(:) :: PARTICLE_STORAGE,WALL_STORAGE,CFACE_STORAGE
   INTEGER :: NLP,NLPDIM,PARTICLE_TAG
   TYPE(HUMAN_TYPE), ALLOCATABLE, DIMENSION(:) :: HUMAN
   INTEGER :: N_HUMANS,N_HUMANS_DIM
   TYPE(HUMAN_GRID_TYPE), ALLOCATABLE, DIMENSION(:,:) :: HUMAN_GRID

   INTEGER :: N_SLCF=0
   TYPE(SLICE_TYPE), ALLOCATABLE, DIMENSION(:) :: SLICE

   INTEGER :: N_RADF=0
   TYPE(RAD_FILE_TYPE), ALLOCATABLE, DIMENSION(:) :: RAD_FILE

   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INC
   INTEGER :: N_PATCH
   TYPE(PATCH_TYPE), ALLOCATABLE, DIMENSION(:) :: PATCH

   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: UIID
   INTEGER :: RAD_CALL_COUNTER,ANGLE_INC_COUNTER

   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: INTERPOLATED_MESH

   CHARACTER(MESH_STRING_LENGTH), ALLOCATABLE, DIMENSION(:) :: STRING
   INTEGER :: N_STRINGS,N_STRINGS_MAX

   INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: K_AGL_SLICE
   INTEGER, ALLOCATABLE, DIMENSION(:,:) :: LS_KLO_TERRAIN,K_LS,LS_SURF_INDEX
   INTEGER :: N_TERRAIN_SLCF=0

   REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: FLUX0_LS,FLUX1_LS,PHI_LS,PHI1_LS,ROS_BACKU, &
                                            ROS_HEAD,ROS_FLANK,WIND_EXP, &
                                            SR_X_LS,SR_Y_LS,U_LS,V_LS,Z_LS,DZTDX,DZTDY,MAG_ZT, &
                                            PHI_WS,UMF,THETA_ELPS,PHI_S,PHI_S_X,PHI_S_Y,PHI_W,LS_WORK1,LS_WORK2

   ! Embedded Mesh

   REAL(EB), ALLOCATABLE, DIMENSION(:,:,:,:) :: SCALAR_SAVE1,SCALAR_SAVE2,SCALAR_SAVE3

END TYPE MESH_TYPE

TYPE (MESH_TYPE), SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: MESHES

END MODULE MESH_VARIABLES


MODULE MESH_POINTERS

USE PRECISION_PARAMETERS
USE MESH_VARIABLES
IMPLICIT NONE

REAL(EB), POINTER, DIMENSION(:,:,:) :: &
   U,V,W,US,VS,WS,DDDT,D,DS,H,HS,H_PRIME,KRES,FVX,FVY,FVZ,FVX_B,FVY_B,FVZ_B,RHO,RHOS, &
   MU,MU_DNS,TMP,Q,KAPPA_GAS,CHI_R,QR,QR_W,UII,RSUM,D_SOURCE,U_OLD,V_OLD,W_OLD, &
   CSD2,MTR,MSR,WEM,MIX_TIME,CHEM_SUBIT,STRAIN_RATE,D_Z_MAX,Q_DOT_PPP_S,PR_T,TMP_FLAME,FLAME_INDEX
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: ZZ,ZZS,REAC_SOURCE_TERM,DEL_RHO_D_DEL_Z,FX,FY,FZ, &
                                         SCALAR_WORK1,SCALAR_WORK2,SCALAR_WORK3,SCALAR_WORK4, &
                                         Q_REAC,AVG_DROP_DEN,AVG_DROP_TMP,AVG_DROP_RAD,AVG_DROP_AREA, &
                                         M_DOT_PPP,M_DOT_G_PPP_S,RHO_ZZ_G_S,TRI_COR, &
                                         ADV_FX,ADV_FY,ADV_FZ,DIF_FX,DIF_FY,DIF_FZ,DIF_FXS,DIF_FYS,DIF_FZS
REAL(EB), POINTER, DIMENSION(:) :: U_EDGE_Y,U_EDGE_Z,V_EDGE_X,V_EDGE_Z,W_EDGE_X,W_EDGE_Y
REAL(EB), POINTER :: POIS_PTB,POIS_ERR,LAPLACE_PTB,LAPLACE_ERR
REAL(EB), POINTER, DIMENSION(:) :: SAVE1,SAVE2,WORK
REAL(EB), POINTER, DIMENSION(:,:,:) :: PRHS
REAL(EB), POINTER, DIMENSION(:,:) :: BXS,BXF,BYS,BYF,BZS,BZF, BXST,BXFT,BYST,BYFT,BZST,BZFT
INTEGER, POINTER :: LSAVE,LWORK,LBC,MBC,NBC,LBC2,MBC2,NBC2,ITRN,JTRN,KTRN,IPS
REAL(EB), POINTER, DIMENSION(:) :: P_0,RHO_0,TMP_0,D_PBAR_DT,D_PBAR_DT_S,U_LEAK,U_DUCT
REAL(EB), POINTER, DIMENSION(:,:) :: PBAR,PBAR_S,R_PBAR
INTEGER, POINTER, DIMENSION(:,:,:) :: PRESSURE_ZONE
REAL(EB), POINTER, DIMENSION(:,:,:) :: WORK1,WORK2,WORK3,WORK4,WORK5,WORK6,WORK7,WORK8,WORK9
REAL(EB), POINTER, DIMENSION(:,:,:) :: DCOR
INTEGER, POINTER, DIMENSION(:,:,:) :: IWORK1

!! Laplace solve, sparse LU
!
!INTEGER, POINTER :: A_N_ROWS,A_N_COLS,A_N_ELEMENTS
!INTEGER, POINTER, DIMENSION(:) :: A_ROW_INDEX,A_COLUMNS
!REAL(EB), POINTER, DIMENSION(:) :: A_VALUES

REAL(EB),     POINTER, DIMENSION(:,:,:) :: PWORK1,PWORK2,PWORK3,PWORK4
COMPLEX(EB),  POINTER, DIMENSION(:,:,:) :: PWORK5,PWORK6,PWORK7,PWORK8

REAL(EB), POINTER, DIMENSION(:,:,:) :: TURB_WORK1,TURB_WORK2,TURB_WORK3,TURB_WORK4
REAL(EB), POINTER, DIMENSION(:,:,:) :: TURB_WORK5,TURB_WORK6,TURB_WORK7,TURB_WORK8
REAL(EB), POINTER, DIMENSION(:,:,:) :: TURB_WORK9,TURB_WORK10

REAL(EB), POINTER, DIMENSION(:,:,:) :: IBM_SAVE1,IBM_SAVE2,IBM_SAVE3,IBM_SAVE4,IBM_SAVE5,IBM_SAVE6

REAL(EB), POINTER, DIMENSION(:) :: WALL_WORK1,WALL_WORK2,FACE_WORK1,FACE_WORK2,FACE_WORK3
REAL(FB), POINTER, DIMENSION(:,:,:,:) :: QQ, QQ2
REAL(FB), POINTER, DIMENSION(:,:) :: PP,PPN,BNDF_TIME_INTEGRAL
INTEGER, POINTER, DIMENSION(:,:) :: IBK
INTEGER, POINTER, DIMENSION(:,:,:) :: IBLK
REAL(EB), POINTER :: CFL,DIVMX,DIVMN,VN,RESMAX,PART_UVWMAX
INTEGER, POINTER :: ICFL,JCFL,KCFL,IMX,JMX,KMX,IMN,JMN,KMN,I_VN,J_VN,K_VN,IRM,JRM,KRM,IUVW
LOGICAL, POINTER, DIMENSION(:) :: APPLY_SPONGE_LAYER
LOGICAL, POINTER :: BAROCLINIC_TERMS_ATTACHED
INTEGER, POINTER :: N_EDGES
INTEGER, POINTER, DIMENSION(:,:) :: IJKE,EDGE_INDEX
REAL(EB), POINTER, DIMENSION(:,:) :: TAU_E,OME_E

INTEGER, POINTER :: MESH_LEVEL,LBC_EMB,MBC_EMB,NBC_EMB
INTEGER, POINTER :: IBAR,JBAR,KBAR,IBM1,JBM1,KBM1,IBP1,JBP1,KBP1
INTEGER, POINTER :: N_NEIGHBORING_MESHES
INTEGER, POINTER, DIMENSION(:) :: NEIGHBORING_MESH

INTEGER, POINTER, DIMENSION(:) :: RGB
REAL(EB), POINTER :: DXI,DETA,DZETA,RDXI,RDETA,RDZETA, &
   DXMIN,DXMAX,DYMIN,DYMAX,DZMIN,DZMAX, &
   XS,XF,YS,YF,ZS,ZF,RDXINT,RDYINT,RDZINT,CELL_SIZE
REAL(EB), POINTER, DIMENSION(:) :: R,RC,X,Y,Z,XC,YC,ZC,HX,HY,HZ, &
   DX,RDX,DXN,RDXN,DY,RDY,DYN,RDYN,DZ,RDZ,DZN,RDZN, &
   CELLSI,CELLSJ,CELLSK,RRN
INTEGER, POINTER :: CELLSI_LO,CELLSI_HI,CELLSJ_LO,CELLSJ_HI,CELLSK_LO,CELLSK_HI
REAL(FB), POINTER, DIMENSION(:) :: XPLT,YPLT,ZPLT
INTEGER, POINTER :: N_OBST
TYPE(OBSTRUCTION_TYPE), POINTER, DIMENSION(:) :: OBSTRUCTION
INTEGER, POINTER :: N_VENT
TYPE(VENTS_TYPE), POINTER, DIMENSION(:) :: VENTS
INTEGER, POINTER, DIMENSION(:,:,:) :: CELL_INDEX
INTEGER, POINTER, DIMENSION(:) :: I_CELL,J_CELL,K_CELL,OBST_INDEX_C
TYPE(IBM_REGFACE_TYPE), POINTER, DIMENSION(:) :: REGFACE_IAXIS_H, REGFACE_JAXIS_H, REGFACE_KAXIS_H
INTEGER, POINTER :: IBM_NEXIMFACE_MESH
INTEGER, POINTER, DIMENSION(:,:,:,:,:) :: FCVAR
INTEGER, POINTER, DIMENSION(:,:,:,:)   :: CCVAR
TYPE(IBM_CUTFACE_TYPE),   POINTER, DIMENSION(:) :: CUT_FACE
TYPE(IBM_CUTCELL_TYPE),   POINTER, DIMENSION(:) :: CUT_CELL
TYPE(IBM_REGFACEZ_TYPE),  POINTER, DIMENSION(:) :: IBM_REGFACE_IAXIS_Z, &
                                                   IBM_REGFACE_JAXIS_Z, &
                                                   IBM_REGFACE_KAXIS_Z
TYPE(IBM_RCFACE_LST_TYPE), POINTER, DIMENSION(:):: IBM_RCFACE_Z
TYPE(IBM_EXIMFACE_TYPE), POINTER, DIMENSION(:)  :: IBM_EXIM_FACE
TYPE(IBM_EDGE_TYPE), POINTER, DIMENSION(:):: IBM_RCEDGE, IBM_IBEDGE
TYPE(CFACE_TYPE), POINTER, DIMENSION(:) :: CFACE
TYPE(RAD_CFACE_TYPE), POINTER, DIMENSION(:) :: RAD_CFACE
REAL(EB), POINTER, DIMENSION(:,:) :: GEOM_ZMAX

INTEGER, POINTER, DIMENSION(:,:) :: WALL_INDEX,WALL_INDEX_HT3D
LOGICAL, POINTER, DIMENSION(:) :: SOLID,EXTERIOR,CONNECTED_MESH
LOGICAL, POINTER, DIMENSION(:,:,:) :: MEAN_FORCING_CELL
INTEGER, POINTER, DIMENSION(:) :: K_MEAN_FORCING
INTEGER, POINTER :: N_WALL_CELLS,N_WALL_CELLS_DIM,N_INTERNAL_WALL_CELLS,N_EXTERNAL_WALL_CELLS,WALL_COUNTER,WALL_COUNTER_HT3D
INTEGER, POINTER :: N_CFACE_CELLS,N_CFACE_CELLS_DIM
REAL(EB),POINTER :: BC_CLOCK,BC_CLOCK_HT3D
REAL(EB), POINTER, DIMENSION(:,:) :: EDGE_INTERPOLATION_FACTOR
REAL(EB), POINTER, DIMENSION(:)   :: D_CORR,DS_CORR,UVW_SAVE,U_GHOST,V_GHOST,W_GHOST
TYPE(WALL_TYPE), POINTER, DIMENSION(:) :: WALL
TYPE(EXTERNAL_WALL_TYPE), POINTER, DIMENSION(:) :: EXTERNAL_WALL
TYPE(OMESH_TYPE), POINTER, DIMENSION(:) :: OMESH
TYPE(LAGRANGIAN_PARTICLE_TYPE), POINTER, DIMENSION(:) :: LAGRANGIAN_PARTICLE
TYPE(STORAGE_TYPE), POINTER, DIMENSION(:) :: WALL_STORAGE,PARTICLE_STORAGE,CFACE_STORAGE
INTEGER, POINTER :: NLP,NLPDIM,PARTICLE_TAG
TYPE(HUMAN_TYPE), POINTER, DIMENSION(:) :: HUMAN
INTEGER, POINTER :: N_HUMANS,N_HUMANS_DIM
TYPE(HUMAN_GRID_TYPE), POINTER, DIMENSION(:,:) :: HUMAN_GRID
INTEGER, POINTER :: N_SLCF
TYPE(SLICE_TYPE), POINTER, DIMENSION(:) :: SLICE
INTEGER, POINTER :: N_RADF
TYPE(RAD_FILE_TYPE), POINTER, DIMENSION(:) :: RAD_FILE
INTEGER, POINTER, DIMENSION(:,:) :: INC
INTEGER, POINTER :: N_PATCH
TYPE(PATCH_TYPE), POINTER, DIMENSION(:) :: PATCH
REAL(EB), POINTER, DIMENSION(:,:,:,:) :: UIID
INTEGER,  POINTER :: RAD_CALL_COUNTER,ANGLE_INC_COUNTER
INTEGER, POINTER, DIMENSION(:,:,:) :: INTERPOLATED_MESH
CHARACTER(MESH_STRING_LENGTH), POINTER, DIMENSION(:) :: STRING
INTEGER, POINTER :: N_STRINGS,N_STRINGS_MAX
INTEGER, POINTER, DIMENSION(:,:,:) :: K_AGL_SLICE
INTEGER, POINTER, DIMENSION(:,:) :: LS_KLO_TERRAIN,K_LS,LS_SURF_INDEX
INTEGER, POINTER :: N_TERRAIN_SLCF
REAL(EB), POINTER, DIMENSION(:,:) :: FLUX0_LS,FLUX1_LS,PHI_LS,PHI1_LS,ROS_BACKU, &
                                     ROS_HEAD,ROS_FLANK,WIND_EXP, &
                                     SR_X_LS,SR_Y_LS,U_LS,V_LS,Z_LS,DZTDX,DZTDY,MAG_ZT, &
                                     PHI_WS,UMF,THETA_ELPS,PHI_S,PHI_S_X,PHI_S_Y,PHI_W,LS_WORK1,LS_WORK2

! Embedded Mesh

REAL(EB), POINTER, DIMENSION(:,:,:,:) :: SCALAR_SAVE1,SCALAR_SAVE2,SCALAR_SAVE3


CONTAINS


SUBROUTINE POINT_TO_MESH(NM)

! Local names for MESH variables point to Global names

INTEGER, INTENT(IN) ::  NM
TYPE (MESH_TYPE), POINTER :: M=>NULL()

M=>MESHES(NM)
U=>M%U
V=>M%V
W=>M%W
US=>M%US
VS=>M%VS
WS=>M%WS
DDDT=>M%DDDT
D=>M%D
DS=>M%DS
H=>M%H
HS=>M%HS
H_PRIME=>M%H_PRIME
KRES=>M%KRES
FVX=>M%FVX
FVY=>M%FVY
FVZ=>M%FVZ
FVX_B=>M%FVX_B
FVY_B=>M%FVY_B
FVZ_B=>M%FVZ_B
RHO=>M%RHO
RHOS=>M%RHOS
TMP=>M%TMP
CHEM_SUBIT=>M%CHEM_SUBIT
MU=>M%MU
PR_T=>M%PR_T
MU_DNS=>M%MU_DNS
D_Z_MAX=>M%D_Z_MAX
CSD2=>M%CSD2
STRAIN_RATE=>M%STRAIN_RATE
MIX_TIME=>M%MIX_TIME
Q=>M%Q
Q_REAC=>M%Q_REAC
Q_DOT_PPP_S=>M%Q_DOT_PPP_S
M_DOT_G_PPP_S=>M%M_DOT_G_PPP_S
RHO_ZZ_G_S=>M%RHO_ZZ_G_S
CHI_R => M%CHI_R
QR=>M%QR
QR_W=>M%QR_W
KAPPA_GAS=>M%KAPPA_GAS
UII=>M%UII
TRI_COR=>M%TRI_COR
FLAME_INDEX=>M%FLAME_INDEX
TMP_FLAME=>M%TMP_FLAME
M_DOT_PPP=>M%M_DOT_PPP
AVG_DROP_DEN=>M%AVG_DROP_DEN
AVG_DROP_AREA=>M%AVG_DROP_AREA
AVG_DROP_TMP=>M%AVG_DROP_TMP
AVG_DROP_RAD=>M%AVG_DROP_RAD
D_SOURCE=>M%D_SOURCE
RSUM=>M%RSUM
ZZ=>M%ZZ
ZZS=>M%ZZS
REAC_SOURCE_TERM=>M%REAC_SOURCE_TERM
DEL_RHO_D_DEL_Z=>M%DEL_RHO_D_DEL_Z
FX=>M%FX
FY=>M%FY
FZ=>M%FZ
U_EDGE_Y=>M%U_EDGE_Y
U_EDGE_Z=>M%U_EDGE_Z
V_EDGE_X=>M%V_EDGE_X
V_EDGE_Z=>M%V_EDGE_Z
W_EDGE_X=>M%W_EDGE_X
W_EDGE_Y=>M%W_EDGE_Y
POIS_PTB=>M%POIS_PTB
POIS_ERR=>M%POIS_ERR
LAPLACE_PTB=>M%LAPLACE_PTB
LAPLACE_ERR=>M%LAPLACE_ERR
SAVE1=>M%SAVE1
SAVE2=>M%SAVE2
!! Laplace solve
!A_N_ROWS=>M%A_N_ROWS
!A_N_COLS=>M%A_N_COLS
!A_N_ELEMENTS=>M%A_N_ELEMENTS
!A_ROW_INDEX=>M%A_ROW_INDEX
!A_COLUMNS=>M%A_COLUMNS
!A_VALUES=>M%A_VALUES
!! -------------
SCALAR_SAVE1=>M%SCALAR_SAVE1
SCALAR_SAVE2=>M%SCALAR_SAVE2
SCALAR_SAVE3=>M%SCALAR_SAVE3
ADV_FX=>M%ADV_FX
ADV_FY=>M%ADV_FY
ADV_FZ=>M%ADV_FZ
DIF_FX=>M%DIF_FX
DIF_FY=>M%DIF_FY
DIF_FZ=>M%DIF_FZ
DIF_FXS=>M%DIF_FXS
DIF_FYS=>M%DIF_FYS
DIF_FZS=>M%DIF_FZS
SCALAR_WORK1=>M%SCALAR_WORK1
SCALAR_WORK2=>M%SCALAR_WORK2
SCALAR_WORK3=>M%SCALAR_WORK3
SCALAR_WORK4=>M%SCALAR_WORK4
WORK=>M%WORK
LSAVE=>M%LSAVE
LWORK=>M%LWORK
PRHS=>M%PRHS
BXS=>M%BXS
BXF=>M%BXF
BYS=>M%BYS
BYF=>M%BYF
BZS=>M%BZS
BZF=>M%BZF
BXST=>M%BXST
BXFT=>M%BXFT
BYST=>M%BYST
BYFT=>M%BYFT
BZST=>M%BZST
BZFT=>M%BZFT
LBC=>M%LBC
MBC=>M%MBC
NBC=>M%NBC
LBC2=>M%LBC2
MBC2=>M%MBC2
NBC2=>M%NBC2
ITRN=>M%ITRN
JTRN=>M%JTRN
KTRN=>M%KTRN
IPS=>M%IPS
LBC_EMB=>M%LBC_EMB
MBC_EMB=>M%MBC_EMB
NBC_EMB=>M%NBC_EMB
U_LEAK=>M%U_LEAK
U_DUCT=>M%U_DUCT
D_PBAR_DT=>M%D_PBAR_DT
D_PBAR_DT_S=>M%D_PBAR_DT_S
PBAR=>M%PBAR
PBAR_S=>M%PBAR_S
R_PBAR=>M%R_PBAR
P_0=>M%P_0
RHO_0=>M%RHO_0
TMP_0=>M%TMP_0
PRESSURE_ZONE=>M%PRESSURE_ZONE
WORK1=>M%WORK1
WORK2=>M%WORK2
WORK3=>M%WORK3
WORK4=>M%WORK4
WORK5=>M%WORK5
WORK6=>M%WORK6
WORK7=>M%WORK7
WORK8=>M%WORK8
WORK9=>M%WORK9
DCOR=>M%DCOR
IWORK1=>M%IWORK1
PWORK1=>M%PWORK1
PWORK2=>M%PWORK2
PWORK3=>M%PWORK3
PWORK4=>M%PWORK4
PWORK5=>M%PWORK5
PWORK6=>M%PWORK6
PWORK7=>M%PWORK7
PWORK8=>M%PWORK8
TURB_WORK1=>M%TURB_WORK1
TURB_WORK2=>M%TURB_WORK2
TURB_WORK3=>M%TURB_WORK3
TURB_WORK4=>M%TURB_WORK4
TURB_WORK5=>M%TURB_WORK5
TURB_WORK6=>M%TURB_WORK6
TURB_WORK7=>M%TURB_WORK7
TURB_WORK8=>M%TURB_WORK8
TURB_WORK9=>M%TURB_WORK9
TURB_WORK10=>M%TURB_WORK10
IBM_SAVE1=>M%IBM_SAVE1
IBM_SAVE2=>M%IBM_SAVE2
IBM_SAVE3=>M%IBM_SAVE3
IBM_SAVE4=>M%IBM_SAVE4
IBM_SAVE5=>M%IBM_SAVE5
IBM_SAVE6=>M%IBM_SAVE6
U_OLD=>M%U_OLD
V_OLD=>M%V_OLD
W_OLD=>M%W_OLD
MESH_LEVEL=>M%MESH_LEVEL
WALL_WORK1=>M%WALL_WORK1
WALL_WORK2=>M%WALL_WORK2
FACE_WORK1=>M%FACE_WORK1
FACE_WORK2=>M%FACE_WORK2
FACE_WORK3=>M%FACE_WORK3
QQ=>M%QQ
QQ2=>M%QQ2
PP=>M%PP
PPN=>M%PPN
BNDF_TIME_INTEGRAL=>M%BNDF_TIME_INTEGRAL
IBK=>M%IBK
IBLK=>M%IBLK
CFL=>M%CFL
DIVMX=>M%DIVMX
VN=>M%VN
RESMAX=>M%RESMAX
DIVMN=>M%DIVMN
PART_UVWMAX=>M%PART_UVWMAX
ICFL=>M%ICFL
JCFL=>M%JCFL
KCFL=>M%KCFL
IMX=>M%IMX
JMX=>M%JMX
KMX=>M%KMX
IMN=>M%IMN
JMN=>M%JMN
KMN=>M%KMN
IRM=>M%IRM
JRM=>M%JRM
KRM=>M%KRM
IUVW=>M%IUVW
APPLY_SPONGE_LAYER=>M%APPLY_SPONGE_LAYER
BAROCLINIC_TERMS_ATTACHED=>M%BAROCLINIC_TERMS_ATTACHED
I_VN=>M%I_VN
J_VN=>M%J_VN
K_VN=>M%K_VN
N_EDGES=>M%N_EDGES
IJKE=>M%IJKE
EDGE_INDEX=>M%EDGE_INDEX
TAU_E=>M%TAU_E
OME_E=>M%OME_E
IBAR=>M%IBAR
JBAR=>M%JBAR
KBAR=>M%KBAR
IBM1=>M%IBM1
JBM1=>M%JBM1
KBM1=>M%KBM1
IBP1=>M%IBP1
JBP1=>M%JBP1
KBP1=>M%KBP1
N_NEIGHBORING_MESHES=>M%N_NEIGHBORING_MESHES
NEIGHBORING_MESH=>M%NEIGHBORING_MESH
RGB=>M%RGB
DXI=>M%DXI
DETA=>M%DETA
DZETA=>M%DZETA
RDXI=>M%RDXI
RDETA=>M%RDETA
RDZETA=>M%RDZETA
DXMIN=>M%DXMIN
DXMAX=>M%DXMAX
DYMIN=>M%DYMIN
DYMAX=>M%DYMAX
DZMIN=>M%DZMIN
DZMAX=>M%DZMAX
CELL_SIZE=>M%CELL_SIZE
XS=>M%XS
XF=>M%XF
YS=>M%YS
YF=>M%YF
ZS=>M%ZS
ZF=>M%ZF
RDXINT=>M%RDXINT
RDYINT=>M%RDYINT
RDZINT=>M%RDZINT
R=>M%R
RC=>M%RC
X=>M%X
Y=>M%Y
Z=>M%Z
XC=>M%XC
YC=>M%YC
ZC=>M%ZC
HX=>M%HX
HY=>M%HY
HZ=>M%HZ
DX=>M%DX
DY=>M%DY
DZ=>M%DZ
DXN=>M%DXN
DYN=>M%DYN
DZN=>M%DZN
RDX=>M%RDX
RDY=>M%RDY
RDZ=>M%RDZ
RDXN=>M%RDXN
RDYN=>M%RDYN
RDZN=>M%RDZN
CELLSI=>M%CELLSI
CELLSJ=>M%CELLSJ
CELLSK=>M%CELLSK
CELLSI_LO=>M%CELLSI_LO
CELLSJ_LO=>M%CELLSJ_LO
CELLSK_LO=>M%CELLSK_LO
CELLSI_HI=>M%CELLSI_HI
CELLSJ_HI=>M%CELLSJ_HI
CELLSK_HI=>M%CELLSK_HI
RRN=>M%RRN
XPLT=>M%XPLT
YPLT=>M%YPLT
ZPLT=>M%ZPLT
N_OBST=>M%N_OBST
OBSTRUCTION=>M%OBSTRUCTION
N_VENT=>M%N_VENT
VENTS=>M%VENTS
CELL_INDEX=>M%CELL_INDEX
I_CELL=>M%I_CELL
J_CELL=>M%J_CELL
K_CELL=>M%K_CELL
REGFACE_IAXIS_H=>M%REGFACE_IAXIS_H
REGFACE_JAXIS_H=>M%REGFACE_JAXIS_H
REGFACE_KAXIS_H=>M%REGFACE_KAXIS_H
FCVAR=>M%FCVAR
CCVAR=>M%CCVAR
CUT_FACE=>M%CUT_FACE
CUT_CELL=>M%CUT_CELL
IBM_REGFACE_IAXIS_Z=>M%IBM_REGFACE_IAXIS_Z
IBM_REGFACE_JAXIS_Z=>M%IBM_REGFACE_JAXIS_Z
IBM_REGFACE_KAXIS_Z=>M%IBM_REGFACE_KAXIS_Z
IBM_RCFACE_Z=>M%IBM_RCFACE_Z
IBM_NEXIMFACE_MESH=>M%IBM_NEXIMFACE_MESH
IBM_EXIM_FACE=>M%IBM_EXIM_FACE
IBM_RCEDGE=>M%IBM_RCEDGE
IBM_IBEDGE=>M%IBM_IBEDGE
CFACE=>M%CFACE
RAD_CFACE=>M%RAD_CFACE
GEOM_ZMAX=>M%GEOM_ZMAX
OBST_INDEX_C=>M%OBST_INDEX_C
WALL_INDEX=>M%WALL_INDEX
WALL_INDEX_HT3D=>M%WALL_INDEX_HT3D
SOLID=>M%SOLID
EXTERIOR=>M%EXTERIOR
CONNECTED_MESH=>M%CONNECTED_MESH
MEAN_FORCING_CELL=>M%MEAN_FORCING_CELL
K_MEAN_FORCING=>M%K_MEAN_FORCING
N_CFACE_CELLS=>M%N_CFACE_CELLS
N_CFACE_CELLS_DIM=>M%N_CFACE_CELLS_DIM
N_WALL_CELLS=>M%N_WALL_CELLS
N_WALL_CELLS_DIM=>M%N_WALL_CELLS_DIM
N_INTERNAL_WALL_CELLS=>M%N_INTERNAL_WALL_CELLS
N_EXTERNAL_WALL_CELLS=>M%N_EXTERNAL_WALL_CELLS
WALL_COUNTER=>M%WALL_COUNTER
WALL_COUNTER_HT3D=>M%WALL_COUNTER_HT3D
BC_CLOCK=>M%BC_CLOCK
BC_CLOCK_HT3D=>M%BC_CLOCK_HT3D
UVW_SAVE=>M%UVW_SAVE
D_CORR=>M%D_CORR
DS_CORR=>M%DS_CORR
U_GHOST=>M%U_GHOST
V_GHOST=>M%V_GHOST
W_GHOST=>M%W_GHOST
EDGE_INTERPOLATION_FACTOR=>M%EDGE_INTERPOLATION_FACTOR
WALL=>M%WALL
EXTERNAL_WALL=>M%EXTERNAL_WALL
OMESH=>M%OMESH
LAGRANGIAN_PARTICLE =>M%LAGRANGIAN_PARTICLE
PARTICLE_STORAGE =>M%PARTICLE_STORAGE
WALL_STORAGE =>M%WALL_STORAGE
CFACE_STORAGE =>M%CFACE_STORAGE
NLP=>M%NLP
NLPDIM=>M%NLPDIM
PARTICLE_TAG=>M%PARTICLE_TAG
HUMAN =>M%HUMAN
N_HUMANS=>M%N_HUMANS
N_HUMANS_DIM=>M%N_HUMANS_DIM
HUMAN_GRID =>M%HUMAN_GRID
N_SLCF=>M%N_SLCF
SLICE=>M%SLICE
N_RADF=>M%N_RADF
RAD_FILE=>M%RAD_FILE
INC=>M%INC
N_PATCH=>M%N_PATCH
PATCH=>M%PATCH
UIID=>M%UIID
RAD_CALL_COUNTER=>M%RAD_CALL_COUNTER
ANGLE_INC_COUNTER=>M%ANGLE_INC_COUNTER
INTERPOLATED_MESH => M%INTERPOLATED_MESH
STRING=>M%STRING
N_STRINGS=>M%N_STRINGS
N_STRINGS_MAX=>M%N_STRINGS_MAX
K_AGL_SLICE   =>M%K_AGL_SLICE
LS_KLO_TERRAIN => M%LS_KLO_TERRAIN
N_TERRAIN_SLCF=>M%N_TERRAIN_SLCF
K_LS  =>M%K_LS
LS_SURF_INDEX =>M%LS_SURF_INDEX
FLUX0_LS => M%FLUX0_LS
FLUX1_LS => M%FLUX1_LS
PHI_LS => M%PHI_LS
PHI1_LS => M%PHI1_LS
ROS_BACKU => M%ROS_BACKU
ROS_HEAD => M%ROS_HEAD
ROS_FLANK => M%ROS_FLANK
WIND_EXP => M%WIND_EXP
SR_X_LS => M%SR_X_LS
SR_Y_LS => M%SR_Y_LS
U_LS => M%U_LS
V_LS => M%V_LS
Z_LS => M%Z_LS
DZTDX => M%DZTDX
DZTDY => M%DZTDY
MAG_ZT => M%MAG_ZT
PHI_WS => M%PHI_WS
UMF => M%UMF
THETA_ELPS => M%THETA_ELPS
PHI_S => M%PHI_S
PHI_S_X => M%PHI_S_X
PHI_S_Y => M%PHI_S_Y
PHI_W => M%PHI_W
LS_WORK1 => M%LS_WORK1
LS_WORK2 => M%LS_WORK2

END SUBROUTINE POINT_TO_MESH

END MODULE MESH_POINTERS
