 
# Generate OBST lines to model a circular pan with a lip.
# Estimate the area the pan surface.
# All dimensions in meters

X0 = 0.                 # x origin
Y0 = 0.                 # y origin
Z1 = -0.06              # height of pan bottom
Z2 = 0.                 # height of pan surface
RAD = 0.1525            # diameter of pan surface
PAN_THICKNESS = 0.015   # thickness of pan lip
LIP_HEIGHT = 0.01       # height of lip above pan surface
IBAR = 120 
JBAR = 120
XS = -0.3
XF = 0.3
YS = -0.3
YF = 0.3

RAD2 = RAD + PAN_THICKNESS
DX = (XF - XS) / float(IBAR)
DY = (YF - YS) / float(JBAR)
AREA = 0.

fid = open('pan_with_lip.txt', 'w')

for J in range(1, JBAR + 1):
    Y = YS + DY * (J - 0.5)
    for I in range(1, IBAR + 1):
        X = XS + DX * (I - 0.5)
        if ((X - X0)**2 + (Y - Y0)**2) < RAD**2:
            fid.write(f"&OBST XB= {X - 0.5 * DX:6.3f}, {X + 0.5 * DX:6.3f}, {Y - 0.5 * DY:6.3f}, {Y + 0.5 * DY:6.3f}, {Z1:6.3f}, {Z2:6.3f}, SURF_ID='POOL' /\n")
            AREA = AREA + DX * DY
        elif ((X - X0)**2 + (Y - Y0)**2) < RAD2**2:
            fid.write(f"&OBST XB= {X - 0.5 * DX:6.3f}, {X + 0.5 * DX:6.3f}, {Y - 0.5 * DY:6.3f}, {Y + 0.5 * DY:6.3f}, {Z1:6.3f}, {Z2 + LIP_HEIGHT:6.3f}, SURF_ID='PAN' /\n")

print(f'Area of pan is {AREA} m2')

fid.close()


