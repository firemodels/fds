#ifndef SMOKEVIEWDEFS_H_DEFINED
#define SMOKEVIEWDEFS_H_DEFINED
#ifdef _DEBUG
void _Sniff_Errors(char *whereat);
#define SNIFF_ERRORS(f) _Sniff_Errors(f)
#else
#define SNIFF_ERRORS(f)
#endif

#define PARTFILE_MAP  0
#define PARTFILE_REMAP 1

#define PARTFILE_LOADALL -11
#define PARTFILE_RELOADALL -12

#define FIRST_TIME 1
#define NOT_FIRST_TIME 2

#define SET_SLICECOLOR 0
#define DEFER_SLICECOLOR 1

#define PARTDATA 0
#define HISTDATA 1

#ifdef pp_SLICEDUP
#define SLICEDUP_KEEPALL 0
#define SLICEDUP_KEEPFINE 1
#define SLICEDUP_KEEPCOARSE 2
#endif

#define SMOKESENSORS_HIDDEN 0
#define SMOKESENSORS_0255 1
#define SMOKESENSORS_01 2
#define SMOKESENSORS_SCALED 3
#define SMOKESENSORS_0INF 4

#define TOURINDEX_ALL  -3
#define TOURINDEX_MANUAL -1
#define TOURINDEX_DEFAULT -4

#define SURFACE_HIDDEN 0
#define SURFACE_SOLID 1
#define SURFACE_OUTLINE 2
#define SURFACE_POINTS 3

#define SHOWALL_FILES 0
#define SHOWONLY_FILE 1
#define HIDEALL_FILES 2

#define UNCOMPRESSED_ALLFRAMES 0
#define UNCOMPRESSED_BYFRAME 1
#define COMPRESSED_ALLFRAMES 2

#define UNCOMPRESSED 0
#define COMPRESSED_ZLIB 1

#define DISABLE 0
#define ENABLE 1

#define XWALLMIN -1
#define XWALLMAX 1
#define YWALLMIN -2
#define YWALLMAX 2
#define ZWALLMIN -3
#define ZWALLMAX 3

#define NOT_FDSBLOCK 0
#define FDSBLOCK 1

#define GEOM_GEOM 0
#define GEOM_ISO 1
#define GEOM_SLICE 2

#define PATCH_NODE_CENTER 0
#define PATCH_CELL_CENTER 1
#define PATCH_GEOMETRY 2

#define NODATA 0
#define HASDATA 1

#define NO_TEST 0
#define TRIANGLE_TEST 1
#define POLYGON_TEST 2
#define TETRAHEDRON_TEST 3

#ifndef UPDATE_SMOKEFIRE_COLORS
#define UPDATE_SMOKEFIRE_COLORS 54
#endif
#define SOOT 1
#define FIRE 2
#define WATER 3

#define NELEV_ZONE 100

#define UPDATE_ISO_OFF 0
#define UPDATE_ISO_ONE_NOW 1
#define UPDATE_ISO_ALL_NOW 2
#define UPDATE_ISO_START_ALL -1

#define MAX_ISO_COLORS 10

#define ZONETEMP_COLOR 0
#define ZONEHAZARD_COLOR 1
#define ZONESMOKE_COLOR 2

#define MAX_HSLABS 10
#define MAX_VSLABS 2
#define MAX_MSLABS 2

#define MAKE_MOVIE 28

#define PNG 0
#define JPEG 1

#define AVI 0
#define MP4 1
#define WMV 2

#define EXTERNAL_LIST_ID 1

#define TEXTURE_SPHERICAL 0
#define TEXTURE_RECTANGULAR 1

#define ADD_KEYFRAME 1
#define DELETE_KEYFRAME -1

#define REL_VIEW 0
#define ABS_VIEW 1

#define IS_AVATAR 1
#define IS_NOT_AVATAR 0

#define C_GENERATED 0
#define FORTRAN_GENERATED 1

#define SHOW_ALL_VENTS 10
#define HIDE_ALL_VENTS 22

#define VENT_SOLID 0
#define VENT_OUTLINE 2
#define VENT_HIDDEN -2

#define HFLOW_VENT 0
#define VFLOW_VENT 1
#define MFLOW_VENT 2

#define CLIP_ON_DENORMAL 2
#define CLIP_ON 1

#define RENDER_ON 1
#define RENDER_OFF 0

#ifndef TYPE_SMV
#define TYPE_SMV 0
#endif
#ifndef TYPE_INI
#define TYPE_INI 1
#endif

#define CLIP_OFF 0
#define CLIP_BLOCKAGES_DATA 1
#define CLIP_BLOCKAGES 2
#define CLIP_DATA 3

#define UNCLIP setClipPlanes(NULL,CLIP_OFF)
#define CLIP setClipPlanes(&clipinfo,CLIP_ON)

#define CLIP_GEOMETRY   \
  {int clip_geom=0;\
    if(clip_mode==CLIP_BLOCKAGES||clip_mode==CLIP_BLOCKAGES_DATA)clip_geom=1;\
    if( clipon==0&&clip_geom==1){CLIP;}\
    else if( clipon==1&&clip_geom==0){UNCLIP;}\
  }

#define CLIP_VALS   \
  {int clip_data=0;\
  if(clip_mode==CLIP_DATA||clip_mode==CLIP_BLOCKAGES_DATA)clip_data=1;\
  if( clipon==0&&clip_data==1){CLIP;}\
    else if( clipon==1&&clip_data==0){UNCLIP;}\
  }

#define GAS 1
#define SOLID 0
#define GASGAS 2
#define SOLIDSOLID 0
#define SOLIDGAS 1
#define GASSOLID 1

#define IN_GAS 0
#define IN_SOLID 1
#define IN_CUTCELL 2

#define EMBED_YES 0
#define EMBED_NO  1

#define XXX 0
#define YYY 1
#define ZZZ 2

#define KEY_ALT 0
#define KEY_CTRL 1
#define KEY_SHIFT 3
#define KEY_NONE 2

#define noGridnoProbe 0
#define GridnoProbe 1
#define GridProbe 2
#define noGridProbe 3

#define FROM_SMOKEVIEW 0
#define FROM_CALLBACK 1
#define FROM_SCRIPT 2

#define STEPS_PER_DEG 10.0

#define FED_SLICE 0
#define FED_ISO 1

#define UNKNOWN -1
#define RLE 0
#define ZLIB 1

#define SLICE_NODE_CENTER 1
#define SLICE_CELL_CENTER 2
#define SLICE_FIRELINE 3
#define SLICE_TERRAIN 4
#define SLICE_FACE_CENTER 5

#define TERRAIN_3D 0
#define TERRAIN_2D_STEPPED 1
#define TERRAIN_2D_LINE 2
#define TERRAIN_3D_MAP 3
#define TERRAIN_HIDDEN 4

#define CSV_FDS 0
#define CSV_CFAST 1
#define CSV_EXP 2

#define CSVTYPE_HRR 1
#define CSVTYPE_DEVC 2
#define CSVTYPE_EXT 3
#define CSVTYPE_NULL 0

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
#define SCRIPT_VOLSMOKERENDERALL 104
#define SCRIPT_RENDERDIR 105
#define SCRIPT_RENDERCLIP 106
#define SCRIPT_SCENECLIP 107
#define SCRIPT_XSCENECLIP 108
#define SCRIPT_YSCENECLIP 109
#define SCRIPT_ZSCENECLIP 110
#define SCRIPT_CBARFLIP 111
#define SCRIPT_CBARNORMAL 112
#define SCRIPT_RENDERSTART 113
#define SCRIPT_MAKEMOVIE 114
#define SCRIPT_RENDERTYPE 115
#define SCRIPT_RENDERSIZE 116
#define SCRIPT_MOVIETYPE 117
#define SCRIPT_ISORENDERALL 118

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
#define SCRIPT_LOADVOLSMOKE 217
#define SCRIPT_LOADVOLSMOKEFRAME 218
#define SCRIPT_LOADISOM 219
#define SCRIPT_LOADBOUNDARYM 220
#define SCRIPT_LOADSLICEM 221
#define SCRIPT_LOADVSLICEM 222

#define SCRIPT_SETTIMEVAL 301
#define SCRIPT_SETVIEWPOINT 302
#define SCRIPT_UNLOADALL 303
#define SCRIPT_KEYBOARD 304
#define SCRIPT_GSLICEVIEW 305
#define SCRIPT_GSLICEPOS 306
#define SCRIPT_GSLICEORIEN 307
#define SCRIPT_SETTOURVIEW 308
#define SCRIPT_SETTOURKEYFRAME 39
#define SCRIPT_EXIT 310
#define SCRIPT_LABEL 311

#define SCRIPT_SLICE_FILE 0
#define SCRIPT_BOUNDARY_FILE 1
#define SCRIPT_SMOKE3D_FILE 2
#define SCRIPT_PART_FILE 3
#define SCRIPT_ISO_FILE 4

#define SCRIPT_UNKNOWN -1

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

#define ROTATION_2AXIS 0
#define EYE_CENTERED 1
#define ROTATION_1AXIS 2
#define ROTATION_3AXIS 3

#define FIRSTCALL 1
#define NOT_FIRSTCALL 0

#define SELECT_BLOCKS 1
#define NOT_SELECT_BLOCKS 0

#define TO_BW 0
#define TO_COLOR 1

#define VENT_CIRCLE 0
#define VENT_RECTANGLE 1
#define VENT_HIDE 2

#define TETRA_CLIPPLANES 1
#define BOX_CLIPPLANES 0

#define DOWN_Y 0
#define UP_X   1
#define UP_Y   2
#define DOWN_X 3
#define DOWN_Z 4
#define UP_Z   5

#define GEOM_STATIC 0
#define GEOM_DYNAMIC 1

#define GEOM_UPDATE_ALL 0
#define GEOM_UPDATE_NORMALS 1

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
#define SHOWALL_PARTICLE 3
#define HIDEALL_ISO 10002
#define SHOWALL_ISO 10001
#define HIDEALL_EVAC 4
#define SHOWALL_EVAC 3


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

#define SLICE_LINE_CONTOUR 0
#define SLICE_SOLID_CONTOUR 1

#define BLOCK_regular 0
#define BLOCK_texture 1
#define BLOCK_outline 2
#define BLOCK_hidden -2

#define BLOCK_face 0
#define VENT_face 1
#define OUTLINE_FRAME_face 2
#define SHADED_FRAME_face 3

#define visBLOCKAsInput 1
#define visBLOCKAsInputOutline 13
#define visBLOCKNormal 8
#define visBLOCKSolidOutline 12
//#define visBLOCKFacet 3
#define visBLOCKOutline 2
#define visBLOCKHide 0
#define visBLOCKTransparent 10
#define visBLOCKAddOutline 14
#define visBLOCKOnlyOutline 15
#define visBLOCKOutlineColor 16
#define visCADOpaque 17

#define OUTLINE_NONE 0
#define OUTLINE_ONLY 1
#define OUTLINE_ADDED 2
#define BLOCKAGE_ASINPUT 0
#define BLOCKAGE_SOLID 1
#define BLOCKAGE_HIDDEN 2

#define BLOCKlocation_grid 5
#define BLOCKlocation_exact 6
#define BLOCKlocation_cad 7
#define BLOCKtexture_cad 31

#define WALL_1 0
#define WALL_3 1
#define WALL_6 2

// (front wall = 1, right wall = 2, back wall = 3, left wall = 4)

#define FRONT_WALL 1
#define RIGHT_WALL 2
#define BACK_WALL 3
#define LEFT_WALL 4
#define BOTTOM_WALL 5
#define TOP_WALL 6

#define XLEFT -1
#define XRIGHT 1
#define YFRONT -2
#define YBACK 2
#define ZBOTTOM -3
#define ZTOP 3

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

#define XDIR 1
#define YDIR 2
#define ZDIR 3
#define XDIRNEG -1
#define YDIRNEG -2
#define ZDIRNEG -3
#define ISO 4

#define NTARGTIMES 100

#define RELOAD_NOW 0
#define STOP_RENDERING -1

#define RELOADALL 4
#define UNLOADALL 1
#define SHOWFILES 5
#define REDIRECT 6

#define SCRIPT_START_RECORDING2 -6
#define SCRIPT_START_RECORDING -2
#define SCRIPT_STOP_RECORDING -3
#define SCRIPT_FILE_LOADING -4
#define SCRIPT_STEP -5
#define SCRIPT_CONTINUE -7
#define SCRIPT_CANCEL -8

#define DRAWSCENE 1
#define SELECTOBJECT 2

#define CORRECT 1

#define VIEW_LEFT 0
#define VIEW_RIGHT 1
#define VIEW_CENTER 2

#define STEREO_NONE 0
#define STEREO_TIME 1
#define STEREO_LR 2
#define STEREO_RB 3
#define STEREO_RC 4
#define STEREO_CUSTOM 5

#define LEFT_EYE 0
#define RIGHT_EYE 1
#define BOTH_EYES 2

#define IINT int
#define UIINT unsigned int
#define FFLOAT float
#define ALL_TRANSPARENT 1
#define ALL_SOLID 4
#define MIN_SOLID 2
#define MAX_SOLID 3
#define HIDE_ALL -1
#define SHOW_ALL -2
#define UNLOAD_ALL -1
#define LOAD_ALL -2
#define HAVE_LIGHT -3
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
#define MAXSMOKERGB 256
#define StepOn 10000
#define RenderCancel 999
#define RENDER_CURRENT_SINGLE 998
#define RENDER_CURRENT_MULTIPLE 978
#define RenderJPEG 997
#define RenderPNG 996
#define Render320 995
#define Render640 994
#define RenderWindow 993
#define RenderCustom 992
#define LABELLEN 30
#define RenderLABELframenumber 980
#define RenderLABELtime 979

#define EXTERIORwallmenu -1
#define INTERIORwallmenu -2
#define FRONTwallmenu -3
#define BACKwallmenu -4
#define LEFTwallmenu -5
#define RIGHTwallmenu -6
#define UPwallmenu -7
#define DOWNwallmenu -8
#define DUMMYwallmenu -9
#define SOLIDpatchmenu -10
#define OUTLINEpatchmenu -11
#define POINTSpatchmenu -12
#define INSOLIDpatchmenu -13
#define INGASpatchmenu -14
#define INCUTCELLpatchmenu -15


#define INTERIORwall 0
#define FRONTwall 1
#define BACKwall 2
#define LEFTwall 3
#define RIGHTwall 4
#define UPwall 5
#define DOWNwall 6

#define offsetscale 100

#define FIRECOLORMAP_DIRECT 0
#define FIRECOLORMAP_CONSTRAINT 1
#define FIRECOLORMAP_NOCONSTRAINT 2

#define RENDER_SLICE 0
#define RENDER_VOLUME 1

#define COLORBAR_FLIP -2
#define COLORBAR_TOGGLE_BW -12
#define COLORBAR_CONTINUOUS -17
#define COLORBAR_STEPPED -18
#define COLORBAR_LINES -19
#define COLORBAR_HIGHLIGHT_BELOW -7
#define COLORBAR_HIGHLIGHT_ABOVE -20
#define COLORBAR_TRANSPARENT -13
#define COLORBAR_RESET -4
#define COLORBAR_TOGGLE_BW_DATA -21

#define LOAD 0
#define UNLOAD 1
#define RESETBOUNDS 2

#define MAKE_SIZEFILE 0
#define GET_DATA 1

#define MAXPLOT3DVARS 6
#define NRGB 12

#define SMALL_FONT 0
#define LARGE_FONT 1
#define SCALED_FONT 2

#define FFALSE 0
#define TTRUE 1

#define BLOCKAGE_AS_INPUT 35
#define BLOCKAGE_AS_INPUT2 36

#define COLORBAR_INDEX_NONE -1

#define BLUE 0
#define GREEN 1
#define YELLOW 2
#define PINK 3
#define RED 4

#define DRAW_OPAQUE 0
#define DRAW_TRANSPARENT 1

#define VOL_READALL -1
#define VOL_UNLOAD -2
#define VOL_READNONE -3

#define MENU_LABEL_colorbar 0
#define MENU_LABEL_timebar 1
#define MENU_LABEL_title 2
#define MENU_LABEL_framerate 3
#define MENU_LABEL_axis 6
#define MENU_LABEL_textlabels 7
#define MENU_LABEL_timelabel 8
#define MENU_LABEL_meshlabel 10
#define MENU_LABEL_memload 11
#define MENU_LABEL_memusage 19
#define MENU_LABEL_fdsticks 12
#define MENU_LABEL_hmslabel 13
#define MENU_LABEL_grid 14
#define MENU_LABEL_sliceaverage 15
#define MENU_LABEL_hrrcutoff 17
#define MENU_LABEL_userticks 18
#define MENU_LABEL_gversion 20
#define MENU_LABEL_ShowAll 4
#define MENU_LABEL_HideAll 5
#define MENU_LABEL_framelabel 9
#define MENU_LABEL_hrr 16
#define MENU_LABEL_northangle 21

#define MENU_TRAINER_smoke 1
#define MENU_TRAINER_temp 2
#define MENU_TRAINER_oxy 3

#define ON 1
#define OFF 0

#define DIALOG_3DSMOKE 20
#define DIALOG_BOUNDS 14
#define DIALOG_CLIP 18
#define DIALOG_COLORBAR 23
#define DIALOG_DEVICE 28
#define DIALOG_DISPLAY 22
#define DIALOG_HIDEALL -2
#define DIALOG_MOTION 29
#define DIALOG_VIEW 30
#define DIALOG_RENDER 31
#define DIALOG_GEOMETRY 16
#define DIALOG_SHOOTER 27
#define DIALOG_SMOKEZIP 24
#define DIALOG_STEREO 19
#define DIALOG_TOUR 21
#define DIALOG_TRAINER 25
#define DIALOG_WUI 26
#define DIALOG_SHOWFILES 33
#define DIALOG_SCRIPT 32
#define DIALOG_CONFIG 34
#define DIALOG_FONTS 35
#define DIALOG_TICKS 36
#define DIALOG_LABELS 37
#define DIALOG_AUTOLOAD 38
#define DIALOG_TIME 39
#define DIALOG_COLORING 40
#define DIALOG_SCALING 41
#define DIALOG_WINDOW 42
#define DIALOG_MOVIE 43

#define UNLOAD_LAST -2

#define UPDATE_PROJECTION -2

#define MENU_TOUR_DEFAULT -1
#define MENU_TOUR_MANUAL -2
#define MENU_TOUR_SHOWALL -3
#define MENU_TOUR_SHOWDIALOG -4
#define MENU_TOUR_VIEWFROMROUTE -5
#define MENU_TOUR_NEW -12
#define MENU_TOUR_CLEARALL -13
#define MENU_TOUR_EDIT -14

#define MENU_TEXTURE_SHOWALL -1
#define MENU_TEXTURE_HIDEALL -2

#define MENU_SHOWHIDE_FLIP 15

#endif

