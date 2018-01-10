PROGRAM PAN_WITH_LIP

! Estimate area of intersection of circle with origin (X0,Y0), radius, RAD, and rectangle (X1,Y1,X2,Y2)

X0 = 0.
Y0 = 0.
Z1 = -0.04
Z2 = 0.
RAD = 0.1524
PAN_THICKNESS = 0.015
IBAR = 30
JBAR = 30
XS = -0.3
XF = 0.3
YS = -.3
YF = 0.3

RAD2 = RAD + PAN_THICKNESS
DX = (XF-XS)/REAL(IBAR)
DY = (YF-YS)/REAL(JBAR)
AREA = 0.

DO J=1,JBAR
   Y = YS + DY*(J-0.5)
   DO I=1,IBAR
      X = XS + DX*(I-0.5)
      IF ( ((X-X0)**2+(Y-Y0)**2)<RAD**2 ) THEN
         WRITE(6,'(A,6(F6.3,","),A)') "&OBST XB=",X-0.5*DX,X+0.5*DX,Y-0.5*DY,Y+0.5*DY,Z1,Z2, " SURF_IDS='POOL','PAN','PAN' /"
         AREA = AREA + DX*DY
      ELSEIF ( ((X-X0)**2+(Y-Y0)**2)<RAD2**2 ) THEN
         WRITE(6,'(A,6(F6.3,","),A)') "&OBST XB=",X-0.5*DX,X+0.5*DX,Y-0.5*DY,Y+0.5*DY,Z1,Z2+0.01, " SURF_ID='PAN' /"
      ENDIF
   ENDDO
ENDDO

WRITE(0,'(A,F8.5,A)') 'Area of pan is ',AREA,' m2'

END PROGRAM PAN_WITH_LIP

