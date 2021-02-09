PROGRAM TUNNEL_E

! Estimate area of intersection of circle with origin (X0,Y0), radius, RAD, and rectangle (X1,Y1,X2,Y2)

Y0 = 0.
Z0 = 0.123
RAD = 0.121
JBAR = 27
KBAR = 24
YS = -0.137
YF = 0.137
ZS = 0.0
ZF = 0.244
PI = 3.141592654

DY = (YF-YS)/REAL(JBAR)
DZ = (ZF-ZS)/REAL(KBAR)
AREA = 0.

DO K=1,KBAR
   Z = ZS + DZ*(K-0.5)
   DO J=1,JBAR
      Y = YS + DY*(J-0.5)
      IF ( Z>0.138 .AND. ((Y-Y0)**2+(Z-Z0)**2)>RAD**2 ) THEN
         WRITE(6,'(A,5(F6.3,","),F6.3,A)') "&OBST XB=",0.00,15.0,Y-0.5*DY,Y+0.5*DY,Z-0.5*DZ,Z+0.5*DZ, " /"
         AREA = AREA + DY*DZ
      ELSEIF ( Z<=0.138 .AND. Z>(ABS(YS)-ABS(Y))*TAN(83*PI/180.)) THEN
         WRITE(6,'(A,5(F6.3,","),F6.3,A)') "&OBST XB=",0.00,15.0,Y-0.5*DY,Y+0.5*DY,Z-0.5*DZ,Z+0.5*DZ, " /"
         AREA = AREA + DY*DZ
      ENDIF
   ENDDO
ENDDO

WRITE(0,'(A,F8.5,A)') 'Area of solid is ',AREA,' m2'

END PROGRAM TUNNEL_E

