PROGRAM correlations

IMPLICIT NONE

CHARACTER(60) :: INPUT_FILE,OUTPUT_FILE
REAL :: ACTIVATION_TEMPERATURE,A_C,ALPHA,AREA,A_T,A_V,C_I,CONDUIT_DIAMETER,CONDUIT_THICKNESS,CUTOFF_TIME,C_PL,C_CJ, &
        C_S,C_STEEL,D,DELTA,DELTA_T_C,DT,EPSILON,FUEL_HEIGHT,F_V,H,HEAT_LOSS_FRACTION,H_C,H_I,H_K,H_V,JACKET_THICKNESS, &
        K_I,K_S,L,LOCATION_FACTOR,L_F,MASS_PER_LENGTH,M_DOT,P,Q,Q_STAR,Q_STEP,R,RADIATIVE_FRACTION,RHO_A, &
        RHO_I,RHO_S,RHO_STEEL,RTI,T,t_activation,T_END,T_P,TMP_A,TMP_G,T_CLOCK,U_JET,V,V_ENT, &
        V_EXP,V_UL,V_DOT,W,W_D,W_V,Z_YT,Z_ASET
REAL, DIMENSION(20) :: X,Z,T_PLUME,R_VALUES,H_VALUES
REAL, DIMENSION(9999) :: TIME_RAMP,Q_RAMP
REAL, DIMENSION(0:5) :: TMP_RAMP,T_RAMP
CHARACTER(20), DIMENSION(20) :: Z_LABEL,LABEL
LOGICAL :: ITER=.TRUE.,PROFILE=.FALSE.,STEEL_PROTECTED=.FALSE., &
           STEEL_UNPROTECTED=.FALSE.,TIME_OUTPUT=.FALSE.,T_SQUARED=.FALSE.
INTEGER :: I,IOS
INTEGER, DIMENSION(20) :: IOR
REAL, PARAMETER :: C_P=1.,G=9.80665,PI=3.141592654

NAMELIST /FPA/ C_S,DELTA,H,K_S,L,M_DOT,OUTPUT_FILE,Q,RHO_S,T_END,TMP_A,W
NAMELIST /DB/  C_S,DELTA,H,K_S,L,M_DOT,OUTPUT_FILE,Q,RHO_S,T_END,TMP_A,W
NAMELIST /MQH/ C_S,DELTA,H,H_V,K_S,L,M_DOT,OUTPUT_FILE,PROFILE,Q,RHO_S,T_END,TMP_A,W,W_V, &
               STEEL_UNPROTECTED,STEEL_PROTECTED,F_V,RHO_STEEL,C_STEEL,H_C,EPSILON,W_D,K_I,RHO_I,C_I,H_I
NAMELIST /BEYLER/ C_S,DELTA,FUEL_HEIGHT,H,HEAT_LOSS_FRACTION,K_S,L,LOCATION_FACTOR,OUTPUT_FILE,Q,RHO_S,T_END,TMP_A,W
NAMELIST /RAD/ AREA,RADIATIVE_FRACTION,OUTPUT_FILE,Q,X,Z,IOR,Z_LABEL,TIME_OUTPUT
NAMELIST /THIEF/ CONDUIT_DIAMETER,CONDUIT_THICKNESS,D,JACKET_THICKNESS,MASS_PER_LENGTH,OUTPUT_FILE,T_END,TMP_A,TMP_RAMP,T_RAMP
NAMELIST /ALPERT/ Q,LOCATION_FACTOR,T_END,R_VALUES,H_VALUES,LABEL,TMP_A,OUTPUT_FILE,T_SQUARED,ALPHA
NAMELIST /SPRINKLER/ T_SQUARED,ALPHA,CUTOFF_TIME,RTI,ACTIVATION_TEMPERATURE,H,R,TMP_A,LOCATION_FACTOR,OUTPUT_FILE
NAMELIST /HESKESTAD/ TIME_RAMP,Q_RAMP,Z,A_C,TMP_A,RADIATIVE_FRACTION,OUTPUT_FILE,Z_LABEL
NAMELIST /MCCAFFREY/ TIME_RAMP,Q_RAMP,Z,TMP_A,OUTPUT_FILE,Z_LABEL,PROFILE,STEEL_UNPROTECTED, &
                     STEEL_PROTECTED,F_V,RHO_STEEL,C_STEEL,H_C,EPSILON,W_D,K_I,RHO_I,C_I,H_I

CALL GET_COMMAND_ARGUMENT(1,INPUT_FILE)

WRITE(0,FMT='(/A,A)') 'Processing ', INPUT_FILE

OPEN(10,FILE=TRIM(INPUT_FILE),FORM='FORMATTED',STATUS='OLD')

! Process FPA lines

DO
   READ(10,NML=FPA,END=100,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_FPA
ENDDO
100 WRITE(0,*) 'Completed FPA.'
REWIND(10)

! Process DB (Deal and Beyler)

DO
   READ(10,NML=DB,END=101,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_DB
ENDDO
101 WRITE(0,*) 'Completed Deal and Beyler.'
REWIND(10)

! Process MQH lines

DO
   READ(10,NML=MQH,END=102,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_MQH
ENDDO
102 WRITE(0,*) 'Completed MQH.'
REWIND(10)

! Process Beyler lines

DO
   READ(10,NML=BEYLER,END=103,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_BEYLER
ENDDO
103 WRITE(0,*) 'Completed Beyler.'
REWIND(10)

! Process Radiation lines

DO
   X=-1. ; Z=-1.
   READ(10,NML=RAD,END=104,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_RADIATION
ENDDO
104 WRITE(0,*) 'Completed Radiation.'
REWIND(10)

! Process THIEF lines

DO
   CONDUIT_DIAMETER=0. ; CONDUIT_THICKNESS=0. ; T_RAMP=0.
   READ(10,NML=THIEF,END=105,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_THIEF
ENDDO
105 WRITE(0,*) 'Completed THIEF.'
REWIND(10)

! Process ALPERT (Ceiling Jet Temperature) lines

DO
   H_VALUES=-1. ; R_VALUES=-1.
   READ(10,NML=ALPERT,END=106,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_ALPERT
ENDDO
106 WRITE(0,*) 'Completed Ceiling Jet Temperatures (Alpert).'
REWIND(10)

! Process SPRINKLER (Sprinkler Activation) lines

DO
   READ(10,NML=SPRINKLER,END=107,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_SPRINKLER
ENDDO
107 WRITE(0,*) 'Completed Sprinkler Activation (Alpert).'
REWIND(10)

! Process HESKESTAD (Plume Centerline Temperature) lines

DO
   TIME_RAMP=-1 ; Q_RAMP=-1 ; Z=-1.
   READ(10,NML=HESKESTAD,END=108,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_HESKESTAD
ENDDO
108 WRITE(0,*) 'Completed Heskestad Plume.'
REWIND(10)

! Process MCCAFFREY (Plume Centerline Temperature) lines

DO
   TIME_RAMP=-1 ; Q_RAMP=-1 ; Z=-1.
   READ(10,NML=MCCAFFREY,END=109,ERR=99,IOSTAT=IOS)
   CALL COMPUTE_MCCAFFREY
ENDDO
109 WRITE(0,*) 'Completed McCaffrey Plume.'
REWIND(10)

! Error message

99 IF (IOS>0) THEN
      WRITE(0,*) 'ERROR: Problem with input line.'
      STOP
   ENDIF

CONTAINS

SUBROUTINE COMPUTE_FPA

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

T_P = (RHO_S*C_S/K_S) * (DELTA/2.)**2
A_T = 2.*L*W + 2.*L*H + 2.*W*H
TMP_A = TMP_A + 273.

WRITE(11,'(A)') 'Time,Temp'

DO I=0,50

   T = I*T_END/50.

   IF (T<=T_P) THEN
      H_K = SQRT(K_S*RHO_S*C_S/(T+0.000000001))
   ELSE
      H_K = K_S/DELTA
   ENDIF

   TMP_G = TMP_A*(1. + 0.63*(Q/(M_DOT*C_P*TMP_A))**0.72 * (H_K*A_T/(M_DOT*C_P))**-0.36)

   WRITE(11,'(F6.1,A1,F6.1)') T,',',TMP_G-273.

ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_FPA


SUBROUTINE COMPUTE_DB

REAL :: HRR

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

A_T = 2.*L*W + 2.*L*H + 2.*W*H
V = L*W*H
TMP_A = TMP_A + 273.
RHO_A = 353./(TMP_A)

WRITE(11,'(A)') 'Time,Temp'

DT = 0.05
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

   IF (T>T_CLOCK) THEN
      WRITE(11,'(F6.1,A1,F6.1)') T,',',TMP_G-273.
      T_CLOCK = T_CLOCK + T_END/500.
   ENDIF

   IF (T>T_END) EXIT

ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_DB


SUBROUTINE COMPUTE_MQH

REAL :: RHO_G,Z,SIGMA
REAL :: T_STEEL(9999),DELTA_T
INTEGER :: J

SIGMA = 5.67e-11

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

T_P = (RHO_S*C_S/K_S) * (DELTA/2.)**2
A_V = H_V*W_V
A_T = 2.*L*W + 2.*L*H + 2.*W*H - A_V
TMP_A = TMP_A + 273.

IF (PROFILE) THEN
   WRITE(11,'(A)') 'Height,Temp'
ELSEIF ((STEEL_UNPROTECTED == .TRUE.) .OR. (STEEL_PROTECTED == .TRUE.)) THEN
   WRITE(11,'(A)') 'Time,Temp,Height,Steel Temperature (C)'
ELSE
   WRITE(11,'(A)') 'Time,Temp,Height'
ENDIF

DO I=0,50

   T = I*T_END/50.

   IF (T<=T_P) THEN
      H_K = SQRT(K_S*RHO_S*C_S/(T+0.000000001))
   ELSE
      H_K = K_S/DELTA
   ENDIF

   TMP_G = TMP_A + 6.85*( Q**2/(A_V*SQRT(H_V)*A_T*H_K) )**(1./3.)

   RHO_G = 353./TMP_G
   Z = MAX( H_V , H*(1. + 2.*(0.05/RHO_G)*Q**0.333*T*H**0.667/(3.*L*W) )**-1.5 )

   ! Compute steel temperatures

   IF (STEEL_UNPROTECTED) THEN
      DELTA_T = 1
      T_STEEL = TMP_A
      DO J=1,T_END
         T_STEEL(J+1) = (F_V * (1/(RHO_STEEL*C_STEEL)) * ((H_C * ((TMP_G)-(T_STEEL(J)))) + ((SIGMA*1000) * EPSILON * ((TMP_G)**4 - (T_STEEL(J))**4))) * DELTA_T) + (T_STEEL(J))
      ENDDO
   ENDIF

   IF (STEEL_PROTECTED) THEN
      DELTA_T = 1
      T_STEEL = TMP_A
      IF (C_STEEL * W_D > 2 * C_I * RHO_I * H_I) THEN
         DO J=1,T_END
            T_STEEL(J+1) = (((K_I/(C_STEEL*H_I*W_D)) * ((TMP_G)-(T_STEEL(J)))) * DELTA_T) + (T_STEEL(J))
         ENDDO
      ELSE
         DO J=1,T_END
            T_STEEL(J+1) = ((((K_I/H_I)/(C_STEEL*W_D + 0.5*C_I*RHO_I*H_I)) * ((TMP_G)-(T_STEEL(J)))) * DELTA_T) + (T_STEEL(J))
         ENDDO
      ENDIF
   ENDIF

   IF ((STEEL_UNPROTECTED == .TRUE.) .OR. (STEEL_PROTECTED == .TRUE.)) THEN
      IF ((.NOT.PROFILE) .AND. (I==0)) WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F7.2)') T,',',TMP_G-273.,',',Z,',',T_STEEL(1)-273
      IF ((.NOT.PROFILE) .AND. (I/=0)) WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F7.2)') T,',',TMP_G-273.,',',Z,',',T_STEEL(T+1)-273
   ELSE
      IF ((.NOT.PROFILE) .AND. (I==0)) WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') T,',',TMP_G-273.,',',Z
      IF ((.NOT.PROFILE) .AND. (I/=0)) WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') T,',',TMP_G-273.,',',Z
   ENDIF

ENDDO

IF (PROFILE) THEN
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') 0.,',',TMP_A-273.
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') Z ,',',TMP_A-273.
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') Z ,',',TMP_G-273.
   WRITE(11,'(F6.1,A1,F6.1,A1,F6.2)') H ,',',TMP_G-273.
ENDIF

CLOSE(11)

STEEL_UNPROTECTED = .FALSE.
STEEL_PROTECTED = .FALSE.
PROFILE=.FALSE.

END SUBROUTINE COMPUTE_MQH


SUBROUTINE COMPUTE_BEYLER

REAL :: K1,K2,M,K,Z,DELTA_T,Z_ASET(9999)
INTEGER :: J

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

A_T = 2.*L*W + 2.*L*H + 2.*W*H
TMP_A = TMP_A + 273.
RHO_A = 353./(TMP_A)
M = L*W*H*RHO_A

WRITE(11,'(A)') 'Time,Temp,HGL Height Yamana Tanaka (m),HGL Height ASET (m)'

DO I=0,50

   T = I*T_END/50.

   K1 = 2*(0.4*SQRT(K_S*RHO_S*C_S))*A_T/(M*C_P)
   K2 = Q/(M*C_P)

   TMP_G = TMP_A + (2.*K2/K1**2)*(K1*SQRT(T) - 1. + EXP(-K1*SQRT(T)))

   ! Calculate HGL height using ASET correlation
   Z_ASET = H
   DELTA_T = 1.

   DO J=1,T_END
      V_EXP = (1 - HEAT_LOSS_FRACTION) * Q / 353.
      V_ENT = ((1 / LOCATION_FACTOR) * 0.21 * (G / (RHO_A * TMP_A))**(1./3.)) * (LOCATION_FACTOR * Q)**(1./3.) * (Z_ASET(J) - FUEL_HEIGHT)**(5./3.)
      V_UL = V_EXP + V_ENT
      Z_ASET(J+1) = Z_ASET(J) - ( (V_UL) / (L * W) ) * DELTA_T
      IF (Z_ASET(J+1) < FUEL_HEIGHT) THEN
         Z_ASET(J+1) = FUEL_HEIGHT
      ENDIF
   ENDDO

   ! Calculate HGL height using Yamana and Tanaka correlation (1985)
   K = 0.076/(353./TMP_G)
   Z_YT = (2*K*Q**(1./3.)*T/(3*L*W) + (1/H**(2./3.)))**(-3./2.)

   IF (I==0) THEN
      WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2)') T,',',TMP_G-273.,',',Z_YT,',',Z_ASET(1)
   ELSE
      WRITE(11,'(F6.1,A1,F6.1,A1,F6.2,A1,F6.2)') T,',',TMP_G-273.,',',Z_YT,',',Z_ASET(T+1)
   ENDIF

ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_BEYLER


! Combined point source radiation model and solid flame model
SUBROUTINE COMPUTE_RADIATION

REAL :: R,Q_RAD(20),Q_RAD_SOLID(20),E,S,H1,H2,A1,A2,F_12_V1,F_12_V2,F_12_H,F_12_V,F_12,A,B
INTEGER :: I,K,N_X,N_Z
CHARACTER(30) :: FMT

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

! Flame height calculation
RHO_A = 353./(273.+20.)
D = SQRT(4.*AREA/PI)
Q_STAR = Q/(RHO_A*C_P*293.*SQRT(G)*D**2.5)
L_F = D*(3.7*Q_STAR**0.4 - 1.02)

! Emissive power calculation
E = 58*(10**(-0.00823*D))

N_X = 0
DO I=1,20
   IF (X(I)<0.) EXIT
   N_X = N_X + 1
   N_Z = 0
   DO K=1,20
      IF (Z(K)<0) EXIT
      N_Z = N_Z + 1

      ! Point source radiation model
      R = SQRT(X(I)**2+(Z(K)-L_F/3.)**2)
      SELECT CASE(ABS(IOR(K)))
      ! Select cos term based on orientation of heat flux gauge
         CASE(1)
            Q_RAD(K) = (X(I)/R)*RADIATIVE_FRACTION*Q/(4.*PI*R**2)
         CASE(2)
            Q_RAD(K) = (X(I)/R)*RADIATIVE_FRACTION*Q/(4.*PI*R**2)
         CASE(3)
            Q_RAD(K) = (Z(K)/R)*RADIATIVE_FRACTION*Q/(4.*PI*R**2)
      END SELECT

      ! Solid flame radiation model
      S = 2*X(I)/D
      IF (Z(I)>0) THEN
         H1 = 2*Z(K)/D
         H2 = 2*(L_F-Z(K))/D
         A1 = (H1**2 + S**2 + 1 ) / (2 * S)
         A2 = (H2**2 + S**2 + 1 ) / (2 * S)
         F_12_V1 = (1/(PI*S)) * ATAN(H1/SQRT(S**2-1)) - (H1/(PI*S))*ATAN(SQRT((S-1)/(S+1))) + (A1*H1/(PI*S*SQRT(A1**2-1))) * ATAN(SQRT((A1+1)*(S-1)/((A1-1)*(S+1))))
         F_12_V2 = (1/(PI*S)) * ATAN(H2/SQRT(S**2-1)) - (H2/(PI*S))*ATAN(SQRT((S-1)/(S+1))) + (A2*H2/(PI*S*SQRT(A2**2-1))) * ATAN(SQRT((A2+1)*(S-1)/((A2-1)*(S+1))))
         F_12 = F_12_V1 + F_12_V2
      ELSE
         H = 2*L_F/D
         A = (H**2 + S**2 + 1 ) / (2 * S)
         B = (1+S**2) / (2*S)
         F_12_H = (1/(PI*S)) * ATAN(H/SQRT(S**2-1)) - (H/(PI*S))*ATAN(SQRT((S-1)/(S+1))) + (A*H/(PI*S*SQRT(A**2-1))) * ATAN(SQRT((A+1)*(S-1)/((A-1)*(S+1))))
         F_12_V = ((B-(1/S))/(PI*SQRT(B**2-1))) * ATAN(SQRT((B+1)*(S-1)/((B-1)*(S+1)))) - ((A-(1/S))/(PI*SQRT(A**2-1))) * ATAN(SQRT((A+1)*(S-1)/((A-1)*(S+1))))
         F_12 = SQRT(F_12_H**2 + F_12_V**2)
      ENDIF
      Q_RAD_SOLID(K) = E * F_12
   ENDDO

   IF (I==1) THEN
      IF (TIME_OUTPUT) THEN
         WRITE(FMT,'(A,I2.1,5A)') "(",2*N_Z,"(","A",",','),","A",")"
         WRITE(11,FMT) 'Time',(TRIM(Z_LABEL(K))//'_PS_RAD',K=1,N_Z),(TRIM(Z_LABEL(K))//'_SF_RAD',K=1,N_Z)
         WRITE(FMT,'(A,I2.1,5A)') "(",2*N_Z,"(","F7.2",",','),","F7.2",")"
         WRITE(11,FMT) 0.,(Q_RAD(K),K=1,N_Z),(Q_RAD_SOLID(K),K=1,N_Z)
         WRITE(11,FMT) 9999.,(Q_RAD(K),K=1,N_Z),(Q_RAD_SOLID(K),K=1,N_Z)
      ELSE
         WRITE(FMT,'(A,I2.1,5A)') "(",2*N_Z,"(","A",",','),","A",")"
         WRITE(11,FMT) 'x',(TRIM(Z_LABEL(K))//'_PS_RAD',K=1,N_Z),(TRIM(Z_LABEL(K))//'_SF_RAD',K=1,N_Z)
      ENDIF
   ENDIF
   
   IF (.NOT. TIME_OUTPUT) THEN
      WRITE(FMT,'(A,I2.1,5A)') "(",2*N_Z,"(","F7.2",",','),","F7.2",")"
      WRITE(11,FMT) X(I),(Q_RAD(K),K=1,N_Z),(Q_RAD_SOLID(K),K=1,N_Z)
   ENDIF
ENDDO

CLOSE(11)

TIME_OUTPUT=.FALSE.

END SUBROUTINE COMPUTE_RADIATION


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

INTEGER :: I,J,K,N_PTS
REAL, DIMENSION(20) :: T_JET, U_JET
CHARACTER(30) :: FMT

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

TMP_A = TMP_A + 273.

! Scaling factor for Q, 2Q, or 4Q (open, wall, or corner fire placement)
Q = Q * LOCATION_FACTOR

DO I=0,50

   T = I*T_END/50.

   IF (T_SQUARED) THEN
      ! Scaling factor for Q, 2Q, or 4Q (open, wall, or corner fire placement)
      Q = ALPHA * T**2 * LOCATION_FACTOR
   ENDIF

   N_PTS = 0
   DO J=1,30
      IF (R_VALUES(J)<0.) EXIT
      N_PTS = N_PTS + 1

      R = R_VALUES(J)
      H = H_VALUES(J)

      ! Compute ceiling jet temperature

      IF (R/H<=0.18) THEN
         T_JET(J) = (16.9 * Q**(2./3.) / H**(5./3.)) + (TMP_A)
      ELSEIF (R/H>0.18) THEN
         T_JET(J) = (5.38 * (Q/R)**(2./3.) / H) + (TMP_A)
      ENDIF

      ! Compute ceiling jet velocity

      IF (R/H<=0.15) THEN
          U_JET = 0.947 * (Q/H)**(1./3.)
      ELSEIF (R/H>0.15) THEN
          U_JET = 0.197 * Q**(1./3.) * H**(1./2.) / R**(5./6.)
      ENDIF
   ENDDO

   IF (I==0) THEN
         WRITE(FMT,'(A,I2.1,5A)') "(",N_PTS*2,"(","A",",','),","A",")"
         WRITE(11,FMT) 'Time',(TRIM(LABEL(K)),K=1,N_PTS),('Velocity '//(TRIM(LABEL(K))),K=1,N_PTS)
         WRITE(FMT,'(A,I2.1,5A)') "(",N_PTS*2,"(","F7.2",",','),","F7.2",")"
         WRITE(11,FMT) T, (T_JET(K)-273,K=1,N_PTS),(U_JET(K),K=1,N_PTS)
      ELSE
         WRITE(FMT,'(A,I2.1,5A)') "(",N_PTS*2,"(","F7.2",",','),","F7.2",")"
         WRITE(11,FMT) T, (T_JET(K)-273,K=1,N_PTS),(U_JET(K),K=1,N_PTS)
   ENDIF

ENDDO

CLOSE(11)

T_SQUARED=.FALSE.

END SUBROUTINE COMPUTE_ALPERT


SUBROUTINE COMPUTE_SPRINKLER

REAL :: Q,T_JET

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

WRITE(11,'(A)') 'Time,Activation,Ceiling jet temperature,Activation time,Total HRR'

TMP_A = TMP_A + 273
ACTIVATION_TEMPERATURE = ACTIVATION_TEMPERATURE + 273
ITER = .TRUE.
T = 0

DO WHILE (ITER)
   IF ((T_SQUARED) .AND. (T<=CUTOFF_TIME)) THEN
      ! Scaling factor for Q, 2Q, or 4Q (open, wall, or corner fire placement)
      Q = ALPHA * T**2 * LOCATION_FACTOR
   ELSEIF ((T_SQUARED) .AND. (T>CUTOFF_TIME)) THEN
      ! Scaling factor for Q, 2Q, or 4Q (open, wall, or corner fire placement)
      Q = ALPHA * CUTOFF_TIME**2 * LOCATION_FACTOR
   ENDIF

   ! Compute ceiling jet temperature

   IF (R/H<=0.18) THEN
      T_JET = (16.9 * Q**(2./3.) / H**(5./3.)) + (TMP_A)
   ELSEIF (R/H>0.18) THEN
      T_JET = (5.38 * (Q/R)**(2./3.) / H) + (TMP_A)
   ENDIF

   ! Compute ceiling jet velocity

   IF (R/H<=0.15) THEN
       U_JET = 0.947 * (Q/H)**(1./3.)
   ELSEIF (R/H>0.15) THEN
       U_JET = 0.197 * Q**(1./3.) * H**(1./2.) / R**(5./6.)
   ENDIF

   ! Compute sprinkler or detector activation time

   t_activation = (RTI / SQRT(U_JET)) * LOG((T_JET - TMP_A)/(T_JET - ACTIVATION_TEMPERATURE))

   IF ((((T_JET - TMP_A)/(T_JET - ACTIVATION_TEMPERATURE))<=0) .OR. (t_activation<=0)) THEN
      WRITE(11,'(F6.1,A1,I2,A1,F6.1,A5,F6.1)') T,',',-1,',',T_JET-273.,',NaN,',Q
   ELSEIF (t_activation>T) THEN
      WRITE(11,'(F6.1,A1,I2,A1,F6.1,A1,F6.1,A1,F6.1)') T,',',-1,',',T_JET-273.,',',t_activation,',',Q
   ELSE
      WRITE(11,'(F6.1,A1,I2,A1,F6.1,A1,F6.1,A1,F6.1)') T,',',1,',',T_JET-273.,',',t_activation,',',Q
   ENDIF

   IF ((t_activation>0) .AND. (t_activation<=T)) THEN
      ITER = .FALSE.
   ELSE
      T = T + 1
   ENDIF
ENDDO

CLOSE(11)

T_SQUARED=.FALSE.

END SUBROUTINE COMPUTE_SPRINKLER


SUBROUTINE COMPUTE_HESKESTAD

INTEGER :: I,J,K,N_T,N_Z
REAL :: Q_C,D,Z_0,T_P
CHARACTER(30) :: FMT

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

TMP_A = TMP_A + 273
RHO_A = 353./(TMP_A)

DO I=1,30
   IF (TIME_RAMP(I)<0.) EXIT
   N_T = N_T + 1

   T = TIME_RAMP(I)

   N_Z = 0
   DO J=1,9999
   IF (Z(J)<0) EXIT
      N_Z = N_Z + 1

      ! Compute convective HRR

      Q_C = Q_RAMP(I) * (1 - RADIATIVE_FRACTION)

      ! Compute fire diameter

      D = SQRT(4.*A_C/PI)

      ! Compute hypothetical virtual origin

      Z_0 = -1.02*D + 0.083*(Q_RAMP(I))**(2./5.)

      ! Compute plume centerline temperature

      T_PLUME(J) = 9.1 * ((TMP_A)/(G*(C_P**2.)*(RHO_A)**2.))**(1./3.) * (Q_C)**(2./3.) * (Z(J)-Z_0)**(-5./3.) + (TMP_A)
   ENDDO

   IF (I==1) THEN
      WRITE(FMT,'(A,I1.1,5A)') "(",N_Z,"(","A",",','),","A",")"
      WRITE(11,FMT) 'Time',(TRIM(Z_LABEL(K)),K=1,N_Z)
   ENDIF
   WRITE(FMT,'(A,I1.1,5A)') "(",N_Z,"(","F8.2",",','),","F8.2",")"
   WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z)
ENDDO

CLOSE(11)

END SUBROUTINE COMPUTE_HESKESTAD


SUBROUTINE COMPUTE_MCCAFFREY

INTEGER :: I,J,K,M,N_T,N_Z,REGION
REAL :: Q,Z_Q_2_5,KAPPA,ETA,T_P,SIGMA,MAX_STEEL_TEMP(20)
REAL :: T_STEEL(9999),DELTA_T
CHARACTER(30) :: FMT

SIGMA = 5.67e-11
TMP_A = TMP_A + 273.

OPEN(11,FILE=TRIM(OUTPUT_FILE),FORM='FORMATTED',STATUS='REPLACE')

DO I=1,9999
   IF ((TIME_RAMP(I)<0.) .AND. (.NOT. PROFILE)) EXIT
   IF ((Z(I)<0.) .AND. (PROFILE)) EXIT
   N_T = N_T + 1

   T = TIME_RAMP(I)
   
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

      T_PLUME(J) = (((KAPPA)/(0.9*SQRT(2*G)))**(2.) * (Z_Q_2_5)**(2*ETA-1) * (TMP_A)) + (TMP_A)

      ! Compute steel temperatures

      IF (STEEL_UNPROTECTED) THEN
         DELTA_T = 1
         T_STEEL = TMP_A
         DO M=1,TIME_RAMP(I)-1
            T_STEEL(M+1) = (F_V * (1/(RHO_STEEL*C_STEEL)) * ((H_C * ((T_PLUME(J))-(T_STEEL(M)))) + ((SIGMA*1000) * EPSILON * ((T_PLUME(J))**4 - (T_STEEL(M))**4))) * DELTA_T) + (T_STEEL(M))
         ENDDO
         MAX_STEEL_TEMP(J) = MAXVAL(T_STEEL)
      ENDIF

      IF (STEEL_PROTECTED) THEN
         DELTA_T = 1
         T_STEEL = TMP_A
         IF (C_STEEL * W_D > 2 * C_I * RHO_I * H_I) THEN
            DO M=1,TIME_RAMP(I)-1
               T_STEEL(M+1) = (((K_I/(C_STEEL*H_I*W_D)) * ((T_PLUME(J))-(T_STEEL(M)))) * DELTA_T) + (T_STEEL(M))
            ENDDO
         ELSE
            DO M=1,TIME_RAMP(I)-1
               T_STEEL(M+1) = ((((K_I/H_I)/(C_STEEL*W_D + 0.5*C_I*RHO_I*H_I)) * ((T_PLUME(J))-(T_STEEL(M)))) * DELTA_T) + (T_STEEL(M))
            ENDDO
         ENDIF
         MAX_STEEL_TEMP(J) = MAXVAL(T_STEEL)
      ENDIF

   ENDDO

   IF ((STEEL_UNPROTECTED == .TRUE.) .OR. (STEEL_PROTECTED == .TRUE.)) THEN
      IF ((I==1) .AND. (.NOT. PROFILE)) THEN
         WRITE(FMT,'(A,I2.1,5A)') "(",N_Z*2,"(","A",",','),","A",")"
         WRITE(11,FMT) 'Time', (TRIM(Z_LABEL(K)),K=1,N_Z), (TRIM('Steel Temperature '//Z_LABEL(K)),K=1,N_Z)
         WRITE(FMT,'(A,I2.1,5A)') "(",N_Z*2,"(","F7.2",",','),","F7.2",")"
         WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z), (MAX_STEEL_TEMP(K)-273,K=1,N_Z)
      ELSEIF (.NOT. PROFILE) THEN
         WRITE(FMT,'(A,I2.1,5A)') "(",N_Z*2,"(","F8.2",",','),","F8.2",")"
         WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z), (MAX_STEEL_TEMP(K)-273,K=1,N_Z)
      ENDIF
   ELSE
      IF ((I==1) .AND. (.NOT. PROFILE)) THEN
         WRITE(FMT,'(A,I2.1,5A)') "(",N_Z,"(","A",",','),","A",")"
         WRITE(11,FMT) 'Time', (TRIM(Z_LABEL(K)),K=1,N_Z)
         WRITE(FMT,'(A,I2.1,5A)') "(",N_Z,"(","F7.2",",','),","F7.2",")"
         WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z)
      ELSEIF (.NOT. PROFILE) THEN
         WRITE(FMT,'(A,I2.1,5A)') "(",N_Z,"(","F7.2",",','),","F7.2",")"
         WRITE(11,FMT) TIME_RAMP(I), (T_PLUME(K)-273,K=1,N_Z)
      ENDIF
   ENDIF
ENDDO

CLOSE(11)

STEEL_UNPROTECTED = .FALSE.
STEEL_PROTECTED = .FALSE.
PROFILE=.FALSE.

END SUBROUTINE COMPUTE_MCCAFFREY


END PROGRAM
