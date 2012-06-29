PROGRAM FDTs

IMPLICIT NONE

CHARACTER(60) :: INPUT_FILE,OUTPUT_FILE
REAL :: ACTIVATION_TEMPERATURE,A_C,ALPHA,AREA,A_T,A_V,C_I,CHI_R,CONDUIT_DIAMETER,CONDUIT_THICKNESS,CUTOFF_TIME,C_PL,C_CJ, &
        C_S,C_STEEL,D,DELTA,DELTA_T_C,DT,EPSILON,F_V,H,H_C,H_K,H_V,JACKET_THICKNESS,K_I,K_S,L,L_F,LEAK_AREA,MASS_PER_LENGTH,M_DOT,P,Q,Q_STAR,Q_STEP, &
        R,RADIATIVE_FRACTION,RHO_AIR,RHO_I,RHO_S,RTI,SCALING_FACTOR,T,t_activation,T_END,T_P,TMP_A,TMP_G,T_CLOCK,U_JET, &
        V,V_DOT,W,W_D,W_V
REAL, DIMENSION(20) :: X,Z,T_PLUME
REAL, DIMENSION(30) :: TIME_RAMP,Q_RAMP
REAL, DIMENSION(0:5) :: TMP_RAMP,T_RAMP
CHARACTER(20), DIMENSION(20) :: Z_LABEL
LOGICAL :: PROFILE=.FALSE.,ITER=.TRUE.,TIME_OUTPUT=.FALSE.,FLIP_AXIS=.FALSE.,STEEL_UNPROTECTED=.FALSE.,STEEL_PROTECTED=.FALSE.
INTEGER :: I,IOS
REAL, PARAMETER :: C_P=1.,G=9.81,GAMMA=1.4,RHO_A=1.2,PI=3.141592654,P_0=101325.
NAMELIST /FPA/ C_S,DELTA,H,K_S,L,M_DOT,OUTPUT_FILE,Q,RHO_S,T_END,TMP_A,W
NAMELIST /DB/  C_S,DELTA,H,K_S,L,LEAK_AREA,M_DOT,OUTPUT_FILE,Q,RHO_S,T_END,TMP_A,W
NAMELIST /MQH/ C_S,DELTA,H,H_V,K_S,L,M_DOT,OUTPUT_FILE,PROFILE,Q,RHO_S,T_END,TMP_A,W,W_V, &
               STEEL_UNPROTECTED,STEEL_PROTECTED,F_V,RHO_S,C_STEEL,H_C,EPSILON,K_I,H,W_D,C_I,RHO_I
NAMELIST /BEYLER/ C_S,DELTA,H,K_S,L,OUTPUT_FILE,Q,RHO_S,T_END,TMP_A,W
NAMELIST /RAD/ AREA,CHI_R,OUTPUT_FILE,Q,X,Z,Z_LABEL,TIME_OUTPUT,FLIP_AXIS
NAMELIST /PRESSURE/ H,L,LEAK_AREA,M_DOT,OUTPUT_FILE,Q,T_END,W
NAMELIST /THIEF/ CONDUIT_DIAMETER,CONDUIT_THICKNESS,D,JACKET_THICKNESS,MASS_PER_LENGTH,OUTPUT_FILE,T_END,TMP_A,TMP_RAMP,T_RAMP
NAMELIST /ALPERT/ ALPHA,CUTOFF_TIME,RTI,ACTIVATION_TEMPERATURE,H,R,TMP_A,RADIATIVE_FRACTION,SCALING_FACTOR,OUTPUT_FILE
NAMELIST /HESKESTAD/ TIME_RAMP,Q_RAMP,Z,A_C,TMP_A,RHO_AIR,RADIATIVE_FRACTION,OUTPUT_FILE,Z_LABEL
NAMELIST /MCCAFFREY/ TIME_RAMP,Q_RAMP,Z,TMP_A,OUTPUT_FILE,Z_LABEL,PROFILE
NAMELIST /MOWRER/ ALPHA,CUTOFF_TIME,C_PL,C_CJ,H,R,OUTPUT_FILE
NAMELIST /MILKE/ ALPHA,CUTOFF_TIME,DELTA_T_C,H,R,RADIATIVE_FRACTION,OUTPUT_FILE

CALL GET_COMMAND_ARGUMENT(1,INPUT_FILE)

WRITE(0,FMT='(/A,A)') 'Processing ', INPUT_FILE

OPEN(10,FILE=TRIM(INPUT_FILE),FORM='FORMATTED',STATUS='OLD')

! Process FPA lines

DO
   READ(10,NML=FPA,END=100,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_FPA
ENDDO
100 WRITE(0,*) 'Done FPA'
REWIND(10)

! Process DB (Deal and Beyler)
DO
   READ(10,NML=DB,END=101,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_DB
ENDDO
101 WRITE(0,*) 'Done Deal and Beyler'
REWIND(10)

! Process MQH lines

DO
   READ(10,NML=MQH,END=102,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_MQH
ENDDO
102 WRITE(0,*) 'Done MQH'
REWIND(10)

! Process Beyler lines

DO
   READ(10,NML=BEYLER,END=103,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_BEYLER
ENDDO
103 WRITE(0,*) 'Done Beyler'
REWIND(10)

! Process Point Source Radiation lines

DO
   X=-1. ; Z=-1.
   READ(10,NML=RAD,END=104,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_POINT_SOURCE_RADIATION
ENDDO
104 WRITE(0,*) 'Done Point Source Radiation'
REWIND(10)

! Process THIEF lines

DO
   CONDUIT_DIAMETER=0. ; CONDUIT_THICKNESS=0. ; T_RAMP=0.
   READ(10,NML=THIEF,END=105,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_THIEF
ENDDO
105 WRITE(0,*) 'Done THIEF'
REWIND(10)

! Process ALPERT (Sprinkler Activation) lines

DO
   READ(10,NML=ALPERT,END=106,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_ALPERT
ENDDO
106 WRITE(0,*) 'Done Sprinkler Activation (Alpert)'
REWIND(10)

! Process HESKESTAD (Plume Centerline Temperature) lines

DO
   TIME_RAMP=-1 ; Q_RAMP=-1 ; Z=-1.
   READ(10,NML=HESKESTAD,END=107,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_HESKESTAD
ENDDO
107 WRITE(0,*) 'Done Heskestad Plume'
REWIND(10)

! Process MCCAFFREY (Plume Centerline Temperature) lines

DO
   TIME_RAMP=-1 ; Q_RAMP=-1 ; Z=-1.
   READ(10,NML=MCCAFFREY,END=108,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_MCCAFFREY
ENDDO
108 WRITE(0,*) 'Done McCaffrey Plume'
REWIND(10)

! Process MOWRER (Smoke Detector Activation) lines

DO
   READ(10,NML=MOWRER,END=109,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_MOWRER
ENDDO
109 WRITE(0,*) 'Done Smoke Detector Activation (Mowrer)'
REWIND(10)

! Process MILKE (Smoke Detector Activation) lines

DO
   READ(10,NML=MILKE,END=110,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_MILKE
ENDDO
110 WRITE(0,*) 'Done Smoke Detector Activation (Milke)'
REWIND(10)

! Error message

99 IF (IOS>0) THEN
      WRITE(0,*) 'ERROR: Problem input line.'
      STOP
   ENDIF

CONTAINS

SUBROUTINE COMPUTE_FPA

REAL :: SIGMA,Q_RAD

SIGMA = 5.67e-11

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

T_P = (RHO_S*C_S/K_S) * (DELTA/2.)**2
A_T = 2.*L*W + 2.*L*H + 2.*W*H
TMP_A = TMP_A + 273.

WRITE(11,'(A)') 'Time, Temp, Immersed HGL Radiation Heat Flux (kW/m2)'

DO I=0,50

   T = I*T_END/50.

   IF (T<T_P) THEN
      H_K = SQRT(K_S*RHO_S*C_S/(T+0.000000001))
   ELSE
      H_K = K_S/DELTA
   ENDIF

   TMP_G = TMP_A*(1. + 0.63*(Q/(M_DOT*C_P*TMP_A))**0.72 * (H_K*A_T/(M_DOT*C_P))**-0.36)

   Q_RAD = SIGMA*(TMP_G**4-TMP_A**4)

   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') T,',',TMP_G-273.,',',Q_RAD

ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_FPA


SUBROUTINE COMPUTE_DB

REAL :: HRR,SIGMA,Q_RAD

SIGMA = 5.67e-11

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

A_T = 2.*L*W + 2.*L*H + 2.*W*H
V = L*W*H
TMP_A = TMP_A + 273.

WRITE(11,'(A)') 'Time, Temp, Pres, Immersed HGL Radiation Heat Flux (kW/m2)'

DT = 0.05
P  = P_0
T  = -60.
T_CLOCK = T 

DO 

   T = T + DT

   IF (T>0.) THEN
      H_K = 0.4*MAX( SQRT(K_S*RHO_S*C_S/(T+0.000000001)) , K_S/DELTA )
      TMP_G = TMP_A + (Q/(M_DOT*C_P + H_K*A_T))
   ELSE
      H_K = 0.
      TMP_G = TMP_A
   ENDIF

   V_DOT = SIGN(1.,P-P_0)*LEAK_AREA*SQRT(2.*ABS(P-P_0)/RHO_A) + M_DOT/RHO_A
   IF (T<=0.) THEN
      HRR = 0.
   ELSE
      HRR = Q
   ENDIF
   P = P + DT*( (GAMMA-1)/V*1000.*(HRR-H_K*A_T*(TMP_G-TMP_A)) - GAMMA*P*V_DOT/V )

   Q_RAD = SIGMA*(TMP_G**4-TMP_A**4)

   IF (T>T_CLOCK) THEN
      WRITE(11,'(F6.1,A1,F6.1,A1,F8.1,A1,F6.2)') T,',',TMP_G-273.,',',P-P_0,',',Q_RAD
      T_CLOCK = T_CLOCK + T_END/500.
   ENDIF

   IF (T>T_END) EXIT

ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_DB


SUBROUTINE COMPUTE_MQH

REAL :: RHO_G,Z,SIGMA,Q_RAD
REAL :: T_STEEL(9999),DELTA_T
INTEGER :: J

SIGMA = 5.67e-11

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

T_P = (RHO_S*C_S/K_S) * (DELTA/2.)**2
A_V = H_V*W_V
A_T = 2.*L*W + 2.*L*H + 2.*W*H - A_V
TMP_A = TMP_A + 273.
T_STEEL = 273

IF (PROFILE) THEN
   WRITE(11,'(A)') 'Height, Temp, Immersed HGL Radiation Heat Flux (kW/m2)'
ELSE
   WRITE(11,'(A)') 'Time, Temp, Height, Immersed HGL Radiation Heat Flux (kW/m2),Steel Temperature (C)'
ENDIF

DO I=0,50

   T = I*T_END/50.

   IF (T<T_P) THEN
      H_K = SQRT(K_S*RHO_S*C_S/(T+0.000000001))
   ELSE
      H_K = K_S/DELTA
   ENDIF

   TMP_G = TMP_A + 6.85*( Q**2/(A_V*SQRT(H_V)*A_T*H_K) )**(1./3.)

   RHO_G = 353./TMP_G
   Z = MAX( H_V , H*(1. + 2.*(0.05/RHO_G)*Q**0.333*T*H**0.667/(3.*L*W) )**-1.5 )

   Q_RAD = SIGMA*(TMP_G**4-TMP_A**4)

   ! Compute steel temperatures

   IF (STEEL_UNPROTECTED) THEN
      DELTA_T = 1
      T_STEEL = 293
      DO J=2,T_END
         T_STEEL(J) = (F_V * (1/(RHO_S*C_STEEL)) * ((H_C * ((TMP_G)-(T_STEEL(J-1)))) + ((SIGMA*1000) * EPSILON * ((TMP_G)**4 - (T_STEEL(J-1))**4))) * DELTA_T) + (T_STEEL(J-1))
      ENDDO
   ENDIF

   IF (STEEL_PROTECTED) THEN
      DELTA_T = 1
      T_STEEL = 293
      DO J=2,T_END
         T_STEEL(J) = (((K_I/(C_STEEL*H*W_D + 0.5*C_I*RHO_I*H**2)) * ((TMP_G)-(T_STEEL(J-1)))) * DELTA_T) + (T_STEEL(J-1))
      ENDDO
   ENDIF

   IF ((.NOT.PROFILE) .AND. (I==0)) WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2,A1,F7.2)') T,',',TMP_G-273.,',',Z,',',Q_RAD,',',T_STEEL(1)-273
   IF ((.NOT.PROFILE) .AND. (I/=0)) WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2,A1,F7.2)') T,',',TMP_G-273.,',',Z,',',Q_RAD,',',T_STEEL(T)-273

ENDDO

IF (PROFILE) THEN
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2)') 0.,',',TMP_A-273.,',',Q_RAD
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2)') Z ,',',TMP_A-273.,',',Q_RAD
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2)') Z ,',',TMP_G-273.,',',Q_RAD
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2)') H ,',',TMP_G-273.,',',Q_RAD
ENDIF

CLOSE(11)

STEEL_UNPROTECTED = .FALSE.
STEEL_PROTECTED = .FALSE.

END SUBROUTINE COMPUTE_MQH


SUBROUTINE COMPUTE_BEYLER

REAL :: K1,K2,M,SIGMA,Q_RAD

SIGMA = 5.67e-11

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

A_T = 2.*L*W + 2.*L*H + 2.*W*H - A_V
TMP_A = TMP_A + 273.
M = L*W*H*RHO_A

WRITE(11,'(A)') 'Time, Temp, Immersed HGL Radiation Heat Flux (kW/m2)'

DO I=0,50

   T = I*T_END/50.

   K1 = 2*(0.4*SQRT(K_S*RHO_S*C_S))*A_T/(M*C_P)
   K2 = Q/(M*C_P)

   TMP_G = TMP_A + (2.*K2/K1**2)*(K1*SQRT(T) - 1. + EXP(-K1*SQRT(T)))

   Q_RAD = SIGMA*(TMP_G**4-TMP_A**4)

   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') T,',',TMP_G-273.,',',Q_RAD

ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_BEYLER


SUBROUTINE COMPUTE_POINT_SOURCE_RADIATION

REAL :: R,Q_RAD(20)
INTEGER :: I,K,N_X,N_Z
CHARACTER(30) :: FMT

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

IF (FLIP_AXIS) THEN
   WRITE(11,'(A)') 'Height (cm),Vertical Flux'
ENDIF

D = SQRT(4.*AREA/PI)
Q_STAR = Q/(RHO_A*C_P*293.*SQRT(G)*D**2.5)
L_F = D*(3.7*Q_STAR**0.4 - 1.02)

N_X = 0
DO I=1,20
   IF (X(I)<0.) EXIT
   N_X = N_X + 1
   N_Z = 0
   DO K=1,20
      IF (Z(K)<0) EXIT
      N_Z = N_Z + 1
      R = SQRT(X(I)**2+(Z(K)-L_F/3.)**2)
      Q_RAD(K) = (X(I)/R)*CHI_R*Q/(4.*PI*R**2)
      IF (FLIP_AXIS) THEN
         WRITE(11,'(F6.2,A1,F6.2)') Z(K)*100,',',Q_RAD(K)
      ENDIF
   ENDDO

   IF ((I==1) .AND. (.NOT. FLIP_AXIS)) THEN
      WRITE(FMT,'(A,I2.1,5A)') "(",N_Z,"(","A",",','),","A",")"
      WRITE(11,FMT) 'x',(TRIM(Z_LABEL(K)),K=1,N_Z)
      IF (TIME_OUTPUT) THEN
         WRITE(FMT,'(A,I2.1,5A)') "(",N_Z,"(","F7.2",",','),","F7.2",")"
         WRITE(11,FMT) 0.,(Q_RAD(K),K=1,N_Z)
         WRITE(11,FMT) 9999.,(Q_RAD(K),K=1,N_Z)
      ENDIF
   ENDIF
   
   IF ((.NOT. TIME_OUTPUT) .AND. (.NOT. FLIP_AXIS)) THEN
      WRITE(FMT,'(A,I2.1,5A)') "(",N_Z,"(","F7.2",",','),","F7.2",")"
      WRITE(11,FMT) X(I),(Q_RAD(K),K=1,N_Z)
   ENDIF
ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_POINT_SOURCE_RADIATION


SUBROUTINE COMPUTE_THIEF

real, allocatable, dimension(:) :: tmp,tmp_next,r
real :: alpha,radius,dt,t,flux,dr,ravg,h,eps,sigma,tmp_s, &
        area,r_critical,c_rho_delta,tmp_conduit,tmp_exp, &
        flux_in,flux_out,tmp_gas,eps_steel,conduit_radius,view_factor,T_WRITE
integer :: i,n_cells,n_steps,i_critical,n
integer, parameter :: n_tmp_g=5
logical :: conduit

! Convert inputs to m,K,kg,J,W, etc.

radius = D/2000.                
jacket_thickness = jacket_thickness/1000.
conduit_thickness = conduit_thickness/1000.
conduit_radius = conduit_diameter/2000. - conduit_thickness
TMP_A = TMP_A + 273.

! Exposing temperature profile

T_RAMP(0) = 0.
TMP_RAMP(0) = TMP_A
do i=1,n_tmp_g
   IF (T_RAMP(I)==0.) EXIT
   TMP_RAMP(i) = TMP_RAMP(i) + 273.                 ! Exposing gas temperature, K
enddo

! Set various constants

area  = pi*radius**2                 ! Cross-sectional area, m2
rho_s = mass_per_length/area         ! Density, kg/m3
c_s   = 1500.                        ! Specific Heat, J/kg/K, fixed
k_s   = 0.2                          ! Conductivity, W/m/K, fixed
eps   = 0.95                         ! Emissivity, fixed
sigma = 5.67e-8                      ! Stefan-Boltzmann, W/m2/K4, fixed
h     = 10.                          ! Convective heat transfer coefficient, W/m2/K
alpha = k_s/(rho_s*c_s)              ! Thermal diffusivity, m2/s
conduit = .false.
if (conduit_thickness>0.) conduit = .true.
c_rho_delta = 460.*7850.*conduit_thickness  ! Steel conduit specific heat*density*thickness
eps_steel = 0.85

! Determine the radial increment and time step

dr  = 0.0002                           ! Desired radial increment, m
n_cells = nint(radius/dr)              ! Number of radial increments
dr  = radius/n_cells                   ! Actual radial increment, m
r_critical = radius - jacket_thickness ! Where failure is assumed to occur, just under the jacket
i_critical = nint(r_critical/dr)       ! Radial cell where failure is to occur
dt  = dr**2/(2.*alpha)                 ! Time step depends on the radial increment

! Open output file

open(11,file=TRIM(OUTPUT_FILE),form='formatted',status='replace')
write(11,'(a)') 'Time,Exposing Temp,Cable Temp,Conduit Temp'

! Allocate arrays

allocate(tmp(0:n_cells+1))
allocate(tmp_next(0:n_cells+1))
allocate(r(0:n_cells))

! Initial conditions

t   = 0.
tmp = TMP_A
tmp_conduit = TMP_A
do i=0,n_cells
   r(i) = i*dr  ! m
enddo

T_WRITE = 0.

! March forward in time

time_loop: do 

   t = t + dt  ! Advance the time

   ! Update the temperature of the cable at the next time step

   do i=1,n_cells
      ravg = 0.5*(r(i)+r(i-1))
      tmp_next(i) = tmp(i) + dt*alpha/(ravg*dr)*( r(i)*(tmp(i+1)-tmp(i))/dr - r(i-1)*(tmp(i)-tmp(i-1))/dr )
   enddo
   tmp_next(0) = tmp_next(1) ! Fill in ghost cell at the cable center, just to avoid numerical issues
   tmp_s = tmp_next(n_cells) ! Assume the surface temperature is the center of the outermost radial increment

   ! Get the exposing gas temperature at time t by linearly interpolating the input temperature profile

   exposure_loop: do n=0,n_tmp_g+1
      if (t>=T_RAMP(n) .and. t<=T_RAMP(n+1)) then
         tmp_gas = TMP_RAMP(n) + (TMP_RAMP(n+1)-TMP_RAMP(n))*(t-T_RAMP(n))/(T_RAMP(n+1)-T_RAMP(n))
         exit exposure_loop
      endif
   enddo exposure_loop

   ! Determine the exposing temperature that the cable sees

   if (conduit) then
      view_factor = 1./(1./eps + (radius/conduit_radius)*(1.-eps_steel)/eps_steel)
      flux_in  = eps_steel*sigma*(tmp_gas**4-tmp_conduit**4) + h*(tmp_gas-tmp_conduit) 
      flux_out = view_factor*(radius/conduit_radius)*sigma*(tmp_conduit**4-tmp_s**4)  + h*(tmp_conduit-tmp_s)
      tmp_conduit = tmp_conduit + dt*(flux_in-flux_out)/c_rho_delta
      tmp_exp = tmp_conduit
   else
      view_factor = eps
      tmp_exp = tmp_gas
   endif

   ! Get the heat flux to the cable surface

   flux = view_factor*sigma*(tmp_exp**4-tmp_s**4) + h*(tmp_exp-tmp_s)  

   ! Set the ghost cell value of the cable temperature just off the surface

   tmp_next(n_cells+1) = tmp_next(n_cells) + dr*flux/k_s   ! Boundary condition at exterior surface

   ! Move the temperatures at the next time step back into the tmp array

   tmp = tmp_next  ! Move update temperatures back to the tmp array

   ! Write output

   IF (T>T_WRITE) THEN
      write(11,'(f6.1,a,f6.1,a,f6.1,a,f6.1)') t,',',tmp_gas-273.,',',tmp(i_critical)-273.,',',tmp_conduit-273.
      T_WRITE = T_WRITE + T_END/50.
   ENDIF

   IF (T>T_END) EXIT time_loop

enddo time_loop

CLOSE(11)

END SUBROUTINE COMPUTE_THIEF


SUBROUTINE COMPUTE_ALPERT

REAL :: Q,Q_C,T_JET

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

WRITE(11,'(A)') 'Time,Ceiling jet temperature,Activation time,Total HRR'

ITER = .TRUE.
T = 0

DO WHILE (ITER)
   IF (T<=CUTOFF_TIME) THEN
      Q = ALPHA * T**2 * SCALING_FACTOR
      Q_C = ALPHA * T**2 * SCALING_FACTOR * (1 - RADIATIVE_FRACTION)
   ELSEIF (T>CUTOFF_TIME) THEN
      Q = ALPHA * CUTOFF_TIME**2 * SCALING_FACTOR
      Q_C = ALPHA * CUTOFF_TIME**2 * SCALING_FACTOR * (1 - RADIATIVE_FRACTION)
   ENDIF

   ! Compute ceiling jet temperature

   IF (R/H<=0.18) THEN
      T_JET = (16.9 * Q_C**(2./3.) / H**(5./3.)) + (TMP_A+273)
   ELSEIF (R/H>0.18) THEN
      T_JET = (5.38 * (Q_C/R)**(2./3.) / H) + (TMP_A+273)
   ENDIF

   ! Compute ceiling jet velocity

   IF (R/H<=0.15) THEN
           U_JET = 0.96 * (Q/H)**(1./3.)
   ELSEIF (R/H>0.15) THEN
           U_JET = 0.195 * Q**(1./3.) * H**(1./2.) / R**(5./6.)
   ENDIF

   ! Compute sprinkler activation time

   t_activation = (RTI / SQRT(U_JET)) * LOG((T_JET-273 - TMP_A)/(T_JET-273 - ACTIVATION_TEMPERATURE))

   IF ((((T_JET-273 - TMP_A)/(T_JET-273 - ACTIVATION_TEMPERATURE))<=0) .OR. (t_activation<=0)) THEN
      WRITE(11,'(F6.1,A1,F6.1,A5,F6.1)') T,',',T_JET-273.,',NaN,',Q
   ELSE
      WRITE(11,'(F6.1,A1,F6.1,A1,F6.1,A1,F6.1)') T,',',T_JET-273.,',',t_activation,',',Q
   ENDIF

   IF ((t_activation>0) .AND. (t_activation<=T)) THEN
      ITER = .FALSE.
   ELSE
      T = T + 1
   ENDIF
ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_ALPERT

SUBROUTINE COMPUTE_HESKESTAD

INTEGER :: I,J,K,N_T,N_Z
REAL :: Q_C,D,Z_0,T_P
CHARACTER(30) :: FMT

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

DO I=1,30
   IF (TIME_RAMP(I)<0.) EXIT
   N_T = N_T + 1
   
   N_Z = 0
   DO J=1,20
   IF (Z(J)<0) EXIT
      N_Z = N_Z + 1
      
      ! Compute convective HRR

      Q_C = Q_RAMP(I) * (1 - RADIATIVE_FRACTION)

      ! Compute fire diameter

      D = SQRT(4.*A_C/PI)

      ! Compute hypothetical virtual origin

      Z_0 = -1.02*D + 0.083*(Q_RAMP(I))**(2./5.)

      ! Compute plume centerline temperature

      T_PLUME(J) = 9.1 * ((TMP_A+273)/(G*(C_P**2.)*(RHO_AIR)**2.))**(1./3.) * (Q_C)**(2./3.) * (Z(J)-Z_0)**(-5./3.) + (TMP_A+273)
   ENDDO
   
   IF (I==1) THEN
      WRITE(FMT,'(A,I1.1,5A)') "(",N_Z,"(","A",",','),","A",")"
      WRITE(11,FMT) 'Time',(TRIM(Z_LABEL(K)),K=1,N_Z)
   ENDIF
   WRITE(FMT,'(A,I1.1,5A)') "(",N_Z,"(","F7.2",",','),","F7.2",")"
   WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z)
ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_HESKESTAD


SUBROUTINE COMPUTE_MCCAFFREY

INTEGER :: I,J,K,N_T,N_Z,REGION
REAL :: Q,Z_Q_2_5,KAPPA,ETA,T_P,SIGMA,Q_RAD(20)
CHARACTER(30) :: FMT

SIGMA = 5.67e-11

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

DO I=1,30
   IF ((TIME_RAMP(I)<0.) .AND. (.NOT. PROFILE)) EXIT
   IF ((Z(I)<0.) .AND. (PROFILE)) EXIT
   N_T = N_T + 1
   
   N_Z = 0
   DO J=1,20
      IF (Z(J)<0) EXIT
      N_Z = N_Z + 1

      Q = Q_RAMP(I)

      ! Compute correlation constants depending on continuous, intermittent, or plume regions

      Z_Q_2_5 = Z(J)/(Q)**(2./5.)

      IF (Z_Q_2_5 < 0.08) THEN
         KAPPA = 6.8
         ETA = 0.5
         REGION = 1
      ELSEIF ((Z_Q_2_5 >= 0.08) .AND. (Z_Q_2_5 <= 0.20)) THEN
         KAPPA = 1.9
         ETA = 0
         REGION = 2
      ELSE
         KAPPA = 1.1
         ETA = -(1./3.)
         REGION = 3
      ENDIF

      ! Compute plume centerline temperature

      T_PLUME(J) = (((KAPPA)/(0.9*SQRT(2*G)))**(2.) * (Z_Q_2_5)**(2*ETA-1) * (TMP_A+273)) + (TMP_A+273)

      ! Compute immersed plume radiation heat flux

      Q_RAD(J) = SIGMA*(T_PLUME(J)**4)

      ! Compute steel temperatures

      ! T_STEEL(J+1) = (F/V * (1/(RHO_S*C_S)) * (H_C * ((T_F+273)-(T_STEEL(J)+273)) + SIGMA * EPSILON * ((T_F+273)**4 - (T_STEEL(J)+273)**4)) * DELTA_T) + (TMP_A+273)

      IF ((I==1) .AND. (J==1) .AND. (PROFILE)) THEN
         WRITE(11,'(A)') 'Height,Flux,Region (1=cont.;2=inter.;3=plume)'
         WRITE(11,'(A5,A1,F6.2,A1,I1)') Z_LABEL(J), ',', Q_RAD(J), ',', REGION
      ELSEIF ((I==1) .AND. (PROFILE)) THEN
         WRITE(11,'(A5,A1,F6.2,A1,I1)') Z_LABEL(J), ',', Q_RAD(J), ',', REGION
      ENDIF

   ENDDO

   IF ((I==1) .AND. (.NOT. PROFILE)) THEN
      WRITE(FMT,'(A,I2.1,5A)') "(",N_Z*2,"(","A",",','),","A",")"
      WRITE(11,FMT) 'Time', (TRIM(Z_LABEL(K)),K=1,N_Z), (TRIM('Radiation '//Z_LABEL(K)),K=1,N_Z)
      WRITE(FMT,'(A,I2.1,5A)') "(",N_Z*2,"(","F7.2",",','),","F7.2",")"
      WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z), (Q_RAD(K),K=1,N_Z)
   ELSEIF (.NOT. PROFILE) THEN
      WRITE(FMT,'(A,I2.1,5A)') "(",N_Z*2,"(","F7.2",",','),","F7.2",")"
      WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z), (Q_RAD(K),K=1,N_Z)
   ENDIF
ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_MCCAFFREY


SUBROUTINE COMPUTE_MOWRER

REAL :: Q,T_PL,T_CJ

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

WRITE(11,'(A)') 'Time,Activation time,Total HRR'

ITER = .TRUE.
T = 1

DO WHILE (ITER)
   IF (T<=CUTOFF_TIME) THEN
      Q = ALPHA * T**2
   ELSEIF (T>CUTOFF_TIME) THEN
      Q = ALPHA * CUTOFF_TIME**2
   ENDIF

   T_PL = C_PL*H**(4./3.) / (Q)**(1./3.)

   T_CJ = R**(11./6.) / (C_CJ*(Q)**(1./3.)*H**(1./2.))

   t_activation = T_PL + T_CJ

   WRITE(11,'(F6.1,A1,F6.1,A1,F6.1)') T,',',t_activation,',',Q

   IF ((t_activation>0) .AND. (t_activation<=T)) THEN
      ITER = .FALSE.
   ELSE
      T = T + 1
   ENDIF
ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_MOWRER


SUBROUTINE COMPUTE_MILKE

REAL :: Q,Y,X

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

WRITE(11,'(A)') 'Time,Activation time,Total HRR'

ITER = .TRUE.
T = 1

DO WHILE (ITER)
   IF (T<=CUTOFF_TIME) THEN
      Q = ALPHA * T**2
   ELSEIF (T>CUTOFF_TIME) THEN
      Q = ALPHA * CUTOFF_TIME**2
   ENDIF

   Y = (DELTA_T_C*1.8)*(H*(1./0.3048))**(5./3.)/(Q*(1./1.05505585))**(2./3.)

   X = 4.6*1E-4*Y**2 + 2.7*1E-15*Y**6

   t_activation = X*(H*(1./0.3048))**(4./3.)/(Q*(1./1.05505585))**(1./3.)

   IF (t_activation>9999) THEN
      WRITE(11,'(F6.1,A5,F6.1)') T,',NaN,',Q
   ELSE
      WRITE(11,'(F6.1,A1,F6.1,A1,F6.1)') T,',',t_activation,',',Q
   ENDIF

   IF ((t_activation>0) .AND. (t_activation<=T)) THEN
      ITER = .FALSE.
   ELSE
      T = T + 1
   ENDIF
ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_MILKE

END PROGRAM
