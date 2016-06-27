@echo off
set chid=%1
set casename=%chid%.fds



> %casename% (
@echo.^&HEAD CHID='%chid%', TITLE='General purpose input file to test FDS timings, SVN $Revision: 19308 $' /
@echo.
@echo.^&MESH IJK=64,64,64, XB=0.0,1.0,0.0,1.0,0.0,1.0 /
@echo.
@echo.^&TIME T_END=1.0 /
@echo.
@echo.^&DUMP NFRAMES=2,DT_DEVC=0.1 /
@echo.
@echo.^&SPEC ID='METHANE' /
@echo.^&SPEC ID='WATER VAPOR' /
@echo.
@echo.^&SURF ID='HOT', VEL=-0.1, TMP_FRONT=100., COLOR='RED' /
@echo.
@echo.^&VENT MB='XMIN', SURF_ID='OPEN' /
@echo.^&VENT MB='XMAX', SURF_ID='OPEN' /
@echo.^&VENT MB='YMIN', SURF_ID='OPEN' /
@echo.^&VENT MB='YMAX', SURF_ID='OPEN' /
@echo.^&VENT PBZ=0.0,   SURF_ID='HOT' /
@echo.^&VENT MB='ZMAX', SURF_ID='OPEN' /
@echo.
@echo.^&DEVC XYZ=0.5,0.5,0.5, QUANTITY='CPU TIME', ID='cpu time' /
@echo.
@echo.^&TAIL /
)
