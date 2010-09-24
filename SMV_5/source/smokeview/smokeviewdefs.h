// $Date$ 
// $Revision$
// $Author$

#ifndef SET_SMOKEVIEWDEFS
#define SET_SMOKEVIEWDEFS
#ifdef _DEBUG
void _sniffErrors(char *whereat);
#define sniffErrors(f) _sniffErrors(f)
#else
#define sniffErrors(f)
#endif

#define TEPS 0.00

#define PART_POINTS 1
#define PART_SPHERES 2
#define PART_LINES 3
#define PART_SMV_DEVICE 4

#define DOUBLE_BUFFER 2
#define SINGLE_BUFFER 1

#define SCRIPT_RENDERONCE 101
#define SCRIPT_RENDERDOUBLEONCE 102
#define SCRIPT_RENDERALL 103
#define SCRIPT_RENDERDIR 104

#define SCRIPT_LOADFILE 201
#define SCRIPT_LOADVFILE 202
#define SCRIPT_LOADBOUNDARY 203
#define SCRIPT_LOAD3DSMOKE 204
#define SCRIPT_LOADISO 205
#define SCRIPT_LOADPARTICLES 206
#define SCRIPT_LOADSLICE 207
#define SCRIPT_LOADVSLICE 208
#define SCRIPT_LOADTOUR 209
#define SCRIPT_UNLOADTOUR 210
#define SCRIPT_PARTCLASSCOLOR 211
#define SCRIPT_PARTCLASSTYPE 212
#define SCRIPT_LOADINIFILE 213
#define SCRIPT_LOADPLOT3D 214
#define SCRIPT_SHOWPLOT3DDATA 215
#define SCRIPT_PLOT3DPROPS 216

#define SCRIPT_SETTIMEVAL 301
#define SCRIPT_SETVIEWPOINT 302
#define SCRIPT_UNLOADALL 303
#define SCRIPT_EXIT 304

#define SCRIPT_SLICE_FILE 0
#define SCRIPT_BOUNDARY_FILE 1
#define SCRIPT_SMOKE3D_FILE 2
#define SCRIPT_PART_FILE 3
#define SCRIPT_ISO_FILE 4

#define DIRX 1
#define DIRY 2
#define DIRZ 3
#define ISO 4

#define PROJECTION 24

#define PARTICLES 0
#define HUMANS 1

#define GLOBAL_INI 0
#define STDOUT_INI 1
#define LOCAL_INI  2
#define SCRIPT_INI 3

#define RESTORE_SAVED_VIEW 2
#define RESTORE_EXTERIOR_VIEW 0
#define RESTORE_INTERIOR_VIEW 1
#define SAVE_VIEW 3
#define RESTORE_LAST_VIEW 4
#define TOGGLE_TITLE_SAFE 5
#define RESTORE_EXTERIOR_VIEW_ZOOM 6

#define WORLD_CENTERED 0
#define EYE_CENTERED 1
#define WORLD_CENTERED_LEVEL 2

#define TO_BW 0
#define TO_COLOR 1

#define DOWN_Y 0
#define UP_X   1 
#define UP_Y   2
#define DOWN_X 3 
#define DOWN_Z 4
#define UP_Z   5
#define NOT_USED 6

#define NO_PLOTS 0
#define STATIC_PLOTS 1
#define STATIC_PLOTS_NORECURSE 3
#define DYNAMIC_PLOTS 2
#define DYNAMIC_PLOTS_NORECURSE 4

#define SHOWALL_PLOT3D 998
#define HIDEALL_PLOT3D 999
#define SHOWALL_BOUNDARY 998
#define HIDEALL_BOUNDARY 999
#define SHOW_CHAR 997
#define HIDEALL_PARTICLE 4
#define SHOWALL_PARTICLE
#define HIDEALL_ISO 10002
#define SHOWALL_ISO 10001


#define TEMP_IGNITION_MAX 100000.
#define SURFACE_TEMPMIN  -100000.
#define SURFACE_TEMPMAX   100000.

#define PERCENTILE_MIN 0
#define SET_MIN 1
#define GLOBAL_MIN 2
#define CHOP_MIN 3

#define PERCENTILE_MAX 0
#define SET_MAX 1
#define GLOBAL_MAX 2
#define CHOP_MAX 3

#define SHADED_CONTOURS 0
#define STEPPED_CONTOURS 1
#define LINE_CONTOURS 2

#define BLOCK_regular 0
#define BLOCK_texture 1
#define BLOCK_outline 2
#define BLOCK_smooth  3
#define BLOCK_hidden -2

#define BLOCK_face 0
#define VENT_face 1
#define OUTLINE_FRAME_face 2
#define SHADED_FRAME_face 3

#define visBLOCKAsInput 1
#define visBLOCKNormal 8
#define visBLOCKSolidOutline 12
//#define visBLOCKFacet 3
//#define visBLOCKSmooth 4
#define visBLOCKOutline 2
#define visBLOCKHide 0
#define visBLOCKSmoothAsNormal 9
#define visBLOCKTransparent 10
#define visSmoothBLOCKSolid 11

#define BLOCKlocation_grid 5
#define BLOCKlocation_exact 6
#define BLOCKlocation_cad 7
#define BLOCKtexture_cad 31
#define SMOOTH_BLOCKAGES 101
#define SMOOTH_ATSTART 102

#define WALL_1 0
#define WALL_3 1
#define WALL_6 2

#define IMIN 0
#define IMAX 1
#define JMIN 2
#define JMAX 3
#define KMIN 4
#define KMAX 5

#define CLOSE_WINDOW -2
#define UPDATE_WINDOW -3
#define CANCEL_WINDOW -4
#define UNDO_BLOCKAGE -5
#define UNDO_ALL_BLOCKAGES -6

#define VISIBLE   1
#define INVISIBLE 0

#define NBUCKETS 1000000 

#define MOVE 0
#define STRETCH_BLACK 1
#define STRETCH_WHITE 2

#define XDIR 0
#define YDIR 1
#define ZDIR 2

#define NTARGTIMES 100

#define RELOADALL 4
#define UNLOADALL 1
#define SHOWFILES 5

#define START_RECORDING_SCRIPT -2
#define STOP_RECORDING_SCRIPT -3

#define RENDER 1
#define SELECT 2

#define CORRECT 1
#define VIEW_LEFT 0
#define VIEW_RIGHT 1
#define VIEW_CENTER 2
#define IINT int
#define UIINT unsigned int
#define FFLOAT float
#define ALL_TRANSPARENT 1
#define ALL_SOLID 4
#define MIN_SOLID 2
#define MAX_SOLID 3
#define HIDE_ALL -1
#define SHOW_ALL -2
#define SHOWALL_SLICE SHOW_ALL
#define SHOWALL_SMOKE3D SHOW_ALL
#define HIDEALL_SLICE HIDE_ALL
#define HIDEALL_SMOKE3D HIDE_ALL
#define HIDEALL_VSLICE HIDE_ALL
#define SHOWALL_VSLICE SHOW_ALL

#define MAXPOINTS 50000000
#define INCFRAMES 20
#define MAXFRAMES 5001
#define PI 3.14159265359f
#define MAXRGB 256
#define StepOn 10000
#define RenderCancel 999
#define RenderOnce 998
#ifdef pp_JPEG
#define RenderJPEG 997
#endif
#ifdef pp_GDGIF
#define RenderGIF 988
#endif
#define RenderPNG 996
#define Render320 995
#define Render640 994
#define RenderWindow 993
#define Render2Window 992
#define LABELLEN 30
#define dwinW 150

#define EXTERIORwallmenu -1
#define INTERIORwallmenu -2
#define FRONTwallmenu -3
#define BACKwallmenu -4
#define LEFTwallmenu -5
#define RIGHTwallmenu -6
#define UPwallmenu -7
#define DOWNwallmenu -8
#define DUMMYwallmenu -9

#define INTERIORwall 0
#define FRONTwall 1
#define BACKwall 2
#define LEFTwall 3
#define RIGHTwall 4
#define UPwall 5
#define DOWNwall 6

#define offsetscale 100

#define NVECLENGTHS 7

#define LOAD 0
#define UNLOAD 1
#define RESETBOUNDS 2

#define MAXPLOT3DVARS 6
#define NRGB 12

#define SMALL_FONT 0
#define LARGE_FONT 1
#define LARGE_FONT_SAFE 2

#define FFALSE 0
#define TTRUE 1

#define BLOCKAGE_AS_INPUT 35
#define BLOCKAGE_AS_INPUT2 36


#define BLUE 0
#define GREEN 1
#define YELLOW 2
#define PINK 3
#define RED 4

#define DRAW_SOLID 0
#define DRAW_TRANSPARENT 1

#ifdef X64
#define STRUCTSTAT struct __stat64
#define STAT _stat64
#else
#define STRUCTSTAT struct stat
#define STAT stat
#endif

#endif

