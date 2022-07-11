PROGRAM CROSS_SECTION

! Estimate area of intersection of circle with origin (X0,Y0), radius, RAD, and rectangle (X1,Y1,X2,Y2)

Y0 = 0.
Z0 = 3.475
RAD = 4.389
JBAR = 32
KBAR = 28
YS = -4.4
YF = 4.4
ZS = 0.0
ZF = 7.9
PI = 3.141592654
XS = 0.0
XF = 856.8

DY = (YF-YS)/REAL(JBAR)
DZ = (ZF-ZS)/REAL(KBAR)
AREA = 0.

DO K=1,KBAR
   Z = ZS + DZ*(K-0.5)
   DO J=1,JBAR
      Y = YS + DY*(J-0.5)
      IF ( Z>Z0 .AND. ((Y-Y0)**2+(Z-Z0)**2)>RAD**2 ) THEN
         WRITE(6,'(A,5(F6.2,","),F6.3,A)') "&OBST XB=",XS,XF,Y-0.5*DY,Y+0.5*DY,Z-0.5*DZ,Z+0.5*DZ, " /"
         AREA = AREA + DY*DZ
      ENDIF
   ENDDO
ENDDO

WRITE(0,'(A,F8.5,A)') 'Area of solid is ',AREA,' m2'

END PROGRAM CROSS_SECTION

