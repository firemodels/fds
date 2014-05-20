@echo off
set chid=%1
set casename=%chid%.fds



> %casename% (
@echo.^&HEAD CHID='%chid%', TITLE='plume fire'  / 
@echo.
@echo.^&MESH IJK=24,24,32, XB=0.0,1.6,0.0,1.6,0.0,3.2 /
@echo.
@echo.^&TIME T_END=10.0 /
@echo.^&DUMP NFRAMES=40 /
@echo.
@echo.^&REAC FUEL         = 'POLYURETHANE'
@echo.      FYI        = 'C_6.3 H_7.1 N O_2.1, NFPA Handbook, Babrauskas'
@echo.      SOOT_YIELD = 0.10
@echo.      N          = 1.0
@echo.      C          = 6.3
@echo.      H          = 7.1
@echo.      O          = 2.1  /
@echo.
@echo.^&SURF ID='BURNER', HRRPUA=600., PART_ID='smoke', COLOR='RASPBERRY' /
@echo.
@echo.^&MATL ID            = 'GYPSUM PLASTER'
@echo.      CONDUCTIVITY  = 0.48
@echo.      SPECIFIC_HEAT = 0.84
@echo.      DENSITY       = 1440. /
@echo.
@echo.^&SURF ID             = 'WALL'
@echo.      DEFAULT        = .TRUE.
@echo.      RGB            = 200,200,200
@echo.      MATL_ID        = 'GYPSUM PLASTER'
@echo.      THICKNESS      = 0.012 /
@echo.
@echo.^&PART ID='smoke', MASSLESS=.TRUE., SAMPLING_FACTOR=1 /
@echo.
@echo.^&VENT XB=0.5,1.1,0.5,1.1,0.1,0.1,SURF_ID='BURNER' /  
@echo.^&OBST XB=0.5,1.1,0.5,1.1,0.0,0.1,SURF_ID='WALL' /
@echo.
@echo.^&VENT MB='XMIN',SURF_ID='OPEN' / 
@echo.^&VENT MB='XMAX',SURF_ID='OPEN' / 
@echo.^&VENT MB='YMIN',SURF_ID='OPEN' / 
@echo.^&VENT MB='YMAX',SURF_ID='OPEN' / 
@echo.^&VENT MB='ZMAX',SURF_ID='OPEN' / 
@echo.
@echo.^&ISOF QUANTITY='TEMPERATURE',VALUE^(1:2^)=30.0,100.0 /  Show 3D contours of temperature at 30 C and 100 C
@echo.
@echo.^&SLCF PBX=0.8,QUANTITY='TEMPERATURE',VECTOR=.TRUE. /  Add vector slices colored by temperature
@echo.^&SLCF PBY=0.8,QUANTITY='TEMPERATURE',VECTOR=.TRUE. /
@echo.^&SLCF PBZ=1.6,QUANTITY='TEMPERATURE',VECTOR=.TRUE. /
@echo.^&SLCF PBZ=3.0,QUANTITY='TEMPERATURE',VECTOR=.TRUE. /
@echo.
@echo.^&BNDF QUANTITY='GAUGE_HEAT_FLUX' /   Common surface quantities. Good for monitoring fire spread.
@echo.^&BNDF QUANTITY='BURNING_RATE' /
@echo.^&BNDF QUANTITY='WALL_TEMPERATURE' /
@echo.
@echo.^&TAIL /
)
