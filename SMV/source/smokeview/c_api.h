int set_slice_bound_min(const char *slice_type, int set, float value);
int set_slice_bound_max(const char *slice_type, int set, float value);
float get_slice_bound_min(const char *slice_type);
float get_slice_bound_max(const char *slice_type);

void printsliceinfo();

int loadsmvall(const char *input_filepath);
int loadsmv(char *input_filename, char *input_filename_ext);
void renderclip(int flag, int left, int right, int bottom, int top);
int render(const char *filename);
void gsliceview(int data, int show_triangles, int show_triangulation,
                int show_normal);
void gslicepos(float x, float y, float z);
void gsliceorien(float az, float elev);
void settourkeyframe(float keyframe_time);
void settourview(int edittourArg, int mode, int show_tourlocusArg,
                 float tour_global_tensionArg);
int getframe();
void setframe(int framenumber);
float gettime();
int settime(float timeval);
int loadfile(const char *filename);
void loadinifile(const char *filepath);
int loadvfile(const char *filepath);
void loadboundaryfile(const char *filepath);
void loadboundary(const char *filepath);
void label(const char *label);
void load3dsmoke(const char *smoke_type);
void set_slice_in_obst(int setting);
int get_slice_in_obst();
void loadvolsmoke(int meshnumber);
void loadvolsmokeframe(int meshnumber, int framenumber, int flag);
void rendertype(const char *type);
int get_rendertype(void);
void set_movietype(const char *type);
int get_movietype(void);
void makemovie(const char *name, const char *base, float framerate);
int loadtour(const char *tourname);
void loadparticles(const char *name);
void partclasscolor(const char *color);
void partclasstype(const char *part_type);
void plot3dprops(int variable_index, int showvector, int vector_length_index,
                 int display_type, float vector_length);
void loadplot3d(int meshnumber, float time_local);
void loadiso(const char *type);
void loadslice(const char *type, int axis, float distance);
void loadvslice(const char *type, int axis, float distance);
void unloadall();
void unloadtour();
void exit_smokeview();
void setcolorbarflip(int flip);
int getcolorbarflip();
int setviewpoint(const char *viewpoint);
int setrenderdir(const char *dir);
void setwindowsize(int width, int height);
void setgridvisibility(int selection);
void setgridparms(int x_vis, int y_vis, int z_vis, int x_plot, int y_plot, int z_plot);
void setcolorbarindex(int chosen_index);
int getcolorbarindex();
void camera_set_eyeview(int eyeview);
void camera_set_rotation_index(int rotation_index);
int camera_get_rotation_index();

void camera_set_rotation_type(int rotation_type);
int camera_get_rotation_type();

void camera_set_view_id(int view_id);
//void viewpoint_set_viewdir(float xcen, float ycen, float zcen);

void camera_mod_eyex(float delta);
void camera_mod_eyey(float delta);
void camera_mod_eyez(float delta);

void camera_set_eyex(float eyex);
void camera_set_eyey(float eyey);
void camera_set_eyez(float eyez);

float camera_get_eyex();
float camera_get_eyey();
float camera_get_eyez();

void camera_mod_az(float delta);
void camera_set_az(float az);
float camera_get_az();
void camera_mod_elev(float delta);
void camera_set_elev(float elev);
float camera_get_elev();

void camera_set_viewdir(float xcen, float ycen, float zcen);
float camera_get_xcen();
float camera_get_ycen();
float camera_get_zcen();
void camera_set_xcen(float xcen);
void camera_set_ycen(float ycen);
void camera_set_zcen(float zcen);

void camera_toggle_projection_type();
int camera_set_projection_type(int projection_type);
int camera_get_projection_type();

void set_sceneclip_x(int clipMin, float min, int clipMax, float max);
void set_sceneclip_x_min(int flag, float value);
void set_sceneclip_x_max(int flag, float value);
void set_sceneclip_y(int clipMin, float min, int clipMax, float max);
void set_sceneclip_y_min(int flag, float value);
void set_sceneclip_y_max(int flag, float value);
void set_sceneclip_z(int clipMin, float min, int clipMax, float max);
void set_sceneclip_z_min(int flag, float value);
void set_sceneclip_z_max(int flag, float value);
void set_clipping_mode(int mode);

int RenderFrameLua(int view_mode, const char *basename);
char* form_filename(int view_mode, char *renderfile_name, char *renderfile_dir,
                   char *renderfile_path, int woffset, int hoffset, int screenH,
                   const char *basename);

int parse_smv_filepath(const char *smv_filepath, char *fdsprefix,
                       char *input_filename_ext);

// --------- show/hide label options--------

// colorbar
void set_colorbar_visibility(int setting);
int get_colorbar_visibility();
void toggle_colorbar_visibility();

// timebar
void set_timebar_visibility(int setting);
int get_timebar_visibility();
void toggle_timebar_visibility();

// title
void set_title_visibility(int setting);
int get_title_visibility();
void toggle_title_visibility();

// axis
void set_axis_visibility(int setting);
int get_axis_visibility();
void toggle_axis_visibility();

// frame
void set_framelabel_visibility(int setting);
int get_framelabel_visibility();
void toggle_framelabel_visibility();

// framerate
void set_framerate_visibility(int setting);
int get_framerate_visibility();
void toggle_framerate_visibility();

// grid locations
void set_gridloc_visibility(int setting);
int get_gridloc_visibility();
void toggle_gridloc_visibility();

// hrrpuv cutoff
void set_hrrcutoff_visibility(int setting);
int get_hrrcutoff_visibility();
void toggle_hrrcutoff_visibility();

// hrr label
void set_hrrlabel_visibility(int setting);
int get_hrrlabel_visibility();
void toggle_hrrlabel_visibility();

// memory load
#ifdef pp_memstatus
void set_memload_visibility(int setting);
int get_memload_visibility();
void toggle_memload_visibility();
#endif

// mesh
void set_meshlabel_visibility(int setting);
int get_meshlabel_visibility();
void toggle_meshlabel_visibility();

// slice average
void set_slice_average_visibility(int setting);
int get_slice_average_visibility();
void toggle_slice_average_visibility();

// time
void set_time_visibility(int setting);
int get_time_visibility();
void toggle_time_visibility();

// user settable ticks
void set_user_ticks_visibility(int setting);
int get_user_ticks_visibility();
void toggle_user_ticks_visibility();

// version info
void set_version_info_visibility(int setting);
int get_version_info_visibility();
void toggle_version_info_visibility();

// set all
void set_all_label_visibility(int setting);

// --------- options--------

// Display Units
// time
void set_timehms(int setting);
int get_timehms();
void toggle_timehms();

// arbitrary
void set_units(int unitclass, int unit_index);
void set_units_default();
void set_unitclass_default(int unitclass);

int blockage_view_method(int setting);
int blockage_outline_color(int setting);
int blockage_locations(int setting);


// .ini config options
int set_ambientlight(float r, float g, float b); // AMBIENTLIGHT
int set_backgroundcolor(float r, float g, float b); // BACKGROUNDCOLOR
int set_blockcolor(float r, float g, float b); // BLOCKCOLOR
int set_blockshininess(float v); // BLOCKSHININESS
int set_blockspecular(float r, float g, float b); // BLOCKSPECULAR
int set_boundcolor(float r, float g, float b); // BOUNDCOLOR
int set_colorbar_textureflag(int v);
int get_colorbar_textureflag();
int set_colorbar_contourvalue(int v);
int set_colorbar_colors(int ncolors, float colors[][3]);
int set_color2bar_colors(int ncolors, float colors[][3]);
int set_diffuselight(float r, float g, float b); // DIFFUSELIGHT
int set_directioncolor(float r, float g, float b); // DIRECTIONCOLOR
int set_flip(int setting); // FLIP
int set_foregroundcolor(float r, float g, float b); // FOREGROUNDCOLOR
int set_heatoffcolor(float r, float g, float b); // HEATOFFCOLOR
int set_heatoncolor(float r, float g, float b); // HEATONCOLOR
// -- ISOCOLORS
// --  10.000000 0.800000 : shininess, default opaqueness
// --  0.700000 0.700000 0.700000 : specular
// --  10 : number of levels
// --  0.957255 0.000392 0.957255 0.800392 : red, green, blue, alpha (opaqueness)
// --  0.750000 0.800000 0.800000 0.800000
// --  0.000000 0.960000 0.280000 0.800000
// --  0.000000 0.000000 1.000000 0.800000
// --  0.000000 0.718750 1.000000 0.800000
// --  0.000000 1.000000 0.562500 0.800000
// --  0.171875 1.000000 0.000000 0.800000
// --  0.890625 1.000000 0.000000 0.800000
// --  1.000000 0.380952 0.000000 0.800000
// --  1.000000 0.000000 0.000000 0.800000
// -- COLORTABLE
// --  2
// --  210 180 140 255 % tan
// --  178 34 34 255 % firebrick
int set_light0(int setting); // LIGHT0
int set_light1(int setting); // LIGHT1
int set_lightmodellocalviewer(int setting); // LIGHTMODELLOCALVIEWER
int set_lightmodelseparatespecularcolor(int setting); // LIGHTMODELSEPARATESPECULARCOLOR
// -- LIGHTPOS0
// --  1.000000 1.000000 4.000000 0.000000
// -- LIGHTPOS1
// --  -1.000000 1.000000 4.000000 0.000000
int set_sensorcolor(float r, float g, float b); // SENSORCOLOR
int set_sensornormcolor(float r, float g, float b); // SENSORNORMCOLOR
int set_bw(int a, int b); // SETBW
int set_sprinkleroffcolor(float r, float g, float b); // SPRINKOFFCOLOR
int set_sprinkleroncolor(float r, float g, float b); // SPRINKONCOLOR
int set_staticpartcolor(float r, float g, float b); // STATICPARTCOLOR
int set_timebarcolor(float r, float g, float b); // TIMEBARCOLOR
int set_ventcolor(float r, float g, float b); // VENTCOLOR

// --    *** SIZES/OFFSETS ***
int set_gridlinewidth(float v); // GRIDLINEWIDTH
int set_isolinewidth(float v); // ISOLINEWIDTH
int set_isopointsize(float v); // ISOPOINTSIZE
int set_linewidth(float v); // LINEWIDTH
int set_partpointsize(float v); // PARTPOINTSIZE
int set_plot3dlinewidth(float v); // PLOT3DLINEWIDTH
int set_plot3dpointsize(float v); // PLOT3DPOINTSIZE
int set_sensorabssize(float v); // SENSORABSSIZE
int set_sensorrelsize(float v); // SENSORRELSIZE
int set_sliceoffset(float v); // SLICEOFFSET
int set_smoothlines(int v); // SMOOTHLINES
int set_spheresegs(int v); // SPHERESEGS
int set_sprinklerabssize(float v); // SPRINKLERABSSIZE
int set_streaklinewidth(float v); // STREAKLINEWIDTH
int set_ticklinewidth(float v); // TICKLINEWIDTH
int set_usenewdrawface(int v); // USENEWDRAWFACE
int set_vecontours(int v); // VECCONTOURS
int set_veclength(int a, float b, float c); // VECLENGTH
int set_vectorlinewidth(float a, float b); // VECTORLINEWIDTH
int set_vectorpointsize(float v); // VECTORPOINTSIZE
int set_ventlinewidth(float v); // VENTLINEWIDTH
int set_ventoffset(float v); // VENTOFFSET
int set_windowoffset(int v); // WINDOWOFFSET
int set_windowwidth(int v); // WINDOWWIDTH
int set_windowheight(int v); // WINDOWHEIGHT

// --  *** DATA LOADING ***
int set_boundzipstep(int v); // BOUNDZIPSTEP
int set_fed(int v); // FED
int set_fedcolorbar(char *name); // FEDCOLORBAR
int set_isozipstep(int v); // ISOZIPSTEP
int set_nopart(int v); // NOPART
int set_partpointstep(int v); // PARTPOINTSTEP
int set_showfedarea(int v); // SHOWFEDAREA
// -- SLICEAVERAGE
// --  0 10.000000 0
int set_slicedataout(int v); // SLICEDATAOUT
int set_slicezipstep(int v); // SLICEZIPSTEP
int set_smoke3dzipstep(int v); // SMOKE3DZIPSTEP
// -- USER_ROTATE
// -- 0 0 0.250000 0.250000 0.500000

// --  *** VIEW PARAMETERS ***

int set_aperature(int v); // APERTURE
int set_axissmooth(int v); // AXISSMOOTH
int set_blocklocation(int v); // BLOCKLOCATION
int set_boundarytwoside(int v); // BOUNDARYTWOSIDE
// -- CLIP
// --  0.001000 3.000000
int set_contourtype(int v); // CONTOURTYPE
int set_cullfaces(int v); // CULLFACES
int set_texturelighting(int v); // ENABLETEXTURELIGHTING
int set_eyeview(int v); // EYEVIEW
int set_eyex(float v); // EYEX
int set_eyey(float v); // EYEY
int set_eyez(float v); // EYEZ
int set_fontsize(int v); // FONTSIZE
int set_frameratevalue(int v); // FRAMERATEVALUE
int set_geomdiags(int a, int b, int c); // GEOMDIAGS
int set_gversion(int v); // GVERSION
int set_isotran2(int v); // ISOTRAN2
int set_meshvis(int n, int vals[]); // MESHVIS
int set_offsetslice(int v); // OFFSETSLICE
int set_outlinemode(int a, int b); // OUTLINEMODE
int set_p3dsurfacetype(int v); // P3DSURFACETYPE
int set_p3dsurfacesmooth(int v); // P3DSURFACESMOOTH
// -- P3VIEW
// --  0 4 1 4 0 8
// --  0 122 1 43 0 15
// --  0 122 1 55 0 15
// --  0 68 1 15 0 15
// --  0 64 1 30 0 30
// --  0 22 1 15 0 15
// --  0 1 1 46 0 15
// --  0 1 1 44 0 15
// --  0 61 1 102 0 15
// --  0 61 1 102 0 7
int set_projection(int v); // PROJECTION
int set_sbatstart(int v); // SBATSTART
// -- SCALEDFONT
// --  12 1.000000 1
// --  32 1.000000 1
int set_showalltextures(int v); // SHOWALLTEXTURES
int set_showsxislabels(int v); // SHOWAXISLABELS
int set_showblocklabel(int v); // SHOWBLOCKLABEL
int set_showblocks(int v); // SHOWBLOCKS
int set_showcadandgrid(int v); // SHOWCADANDGRID
int set_showcadopaque(int v); // SHOWCADOPAQUE
int set_showceiling(int v); // SHOWCEILING
int set_showcolorbars(int v); // SHOWCOLORBARS
int set_showcvents(int a, int b); // SHOWCVENTS
int set_showdummyvents(int v); // SHOWDUMMYVENTS
int set_showevacslices(int a, int b, int c); // SHOWEVACSLICES
int set_showfloor(int v); // SHOWFLOOR
int set_showframe(int v); // SHOWFRAME
int set_showframelabel(int v); // SHOWFRAMELABEL
int set_showframerate(int v); // SHOWFRAMERATE
int set_showgrid(int v); // SHOWGRID
int set_showgridloc(int v); // SHOWGRIDLOC
int set_showhmstimelabel(int v); // SHOWHMSTIMELABEL
int set_showhrrcutoff(int v); // SHOWHRRCUTOFF
int set_showiso(int v); // SHOWISO
int set_showisonormals(int v); // SHOWISONORMALS
int set_showlabels(int v); // SHOWLABELS
int set_showmemload(int v); // SHOWMEMLOAD
int set_shownormalwhensmooth(int v); // SHOWNORMALWHENSMOOTH
int set_showopenvents(int a, int b); // SHOWOPENVENTS
int set_showothervents(int v); // SHOWOTHERVENTS
int set_showsensors(int a, int b); // SHOWSENSORS
int set_showsliceinobst(int v); // SHOWSLICEINOBST
int set_showsmokepart(int v); // SHOWSMOKEPART
int set_showsprinkpart(int v); // SHOWSPRINKPART
int set_showstreak(int a, int b, int c, int d); // SHOWSTREAK
int set_showterrain(int v); // SHOWTERRAIN
int set_showtetras(int a, int b); // SHOWTETRAS
int set_showthreshold(int a, int b, float c); // SHOWTHRESHOLD
int set_showticks(int v); // SHOWTICKS
int set_showtimebar(int v); // SHOWTIMEBAR
int set_showtimelabel(int v); // SHOWTIMELABEL
int set_showtitle(int v); // SHOWTITLE
int set_showtracersalways(int v); // SHOWTRACERSALWAYS
int set_showtriangles(int a, int b, int c, int d, int e, int f); // SHOWTRIANGLES
int set_showtransparent(int v); // SHOWTRANSPARENT
int set_showtransparentvents(int v); // SHOWTRANSPARENTVENTS
int set_showtrianglecount(int v); // SHOWTRIANGLECOUNT
int set_showventflow(int a, int b, int c); // SHOWVENTFLOW
int set_showvents(int v); // SHOWVENTS
int set_showwalls(int v); // SHOWWALLS
int set_skipembedslice(int v); // SKIPEMBEDSLICE
int set_smokesensors(int a, int b); // SMOKESENSORS
int set_smoothblocksolid(int v); // SMOOTHBLOCKSOLID
int set_startuplang(char *lang); // STARTUPLANG
int set_stereo(int v); // STEREO
int set_surfinc(int v); // SURFINC
// -- TERRAINPARMS
// --  90 50 50
// --  200 200 200
// --  1.000000
int set_titlesafe(int v); // TITLESAFE
int set_trainerview(int v); // TRAINERVIEW
int set_transparent(int a, float b); // TRANSPARENT
int set_twosidedvents(int a, int b); // TWOSIDEDVENTS
int set_vectorskip(int v); // VECTORSKIP
// -- VOLSMOKE
// --  0 1 1 0 0
// --  20.000000 700.000000 1200.000000 3.000000 8700.000000 1.000000 1.000000
// -- ZOOM
// --  -2 1.000000

// --  *** MISC ***

int set_cellcentertext(int v); // CELLCENTERTEXT
int set_inputfile(char *filename); // INPUT_FILE
int set_labelstartupview(char *filename); // LABELSTARTUPVIEW
int set_pixelskip(int v); // PIXELSKIP
int set_renderclip(int a, int b, int c, int d, int e); // RENDERCLIP
int set_renderfilelabel(int v); // RENDERFILELABEL
int set_renderfiletype(int a, int b); // RENDERFILETYPE
int set_renderoption(int a, int b); // RENDEROPTION
int set_unticlasses(int n, int values[]); // UNITCLASSES

// --  *** 3D SMOKE INFO ***

int set_adjustalpha(int v); // ADJUSTALPHA
int set_colorbartype(int v, char *label); // COLORBARTYPE
int set_extremecolors(int a, int b, int c, int d, int e, int f); // EXTREMECOLORS
int set_firecolor(int r, int g, int b); // FIRECOLOR
int set_firecolormap(int a, int b); // FIRECOLORMAP
int set_firedepth(float v); // FIREDEPTH
int set_showextremedata(int a, int b, int c); // SHOWEXTREMEDATA
int set_smokecolor(int r, int g, int b); // SMOKECOLOR
int set_smokecull(int v); // SMOKECULL
int set_smokeskip(int v); // SMOKESKIP
int set_smokealbedo(float v); // SMOKEALBEDO
int set_smokerthick(float v); // SMOKERTHICK
int set_usegpu(int v); // USEGPU
// -- VOLSMOKE
// --  0 1 1 0 0
// --  20.000000 700.000000 1200.000000 3.000000 8700.000000

// --  *** ZONE FIRE PARAMETRES ***

int set_showhazardcolors(int v); // SHOWHAZARDCOLORS
int set_showhzone(int v); // SHOWHZONE
int set_showszone(int v); // SHOWSZONE
int set_showvzone(int v); // SHOWVZONE
int set_showzonefire(int v); // SHOWZONEFIRE

// --  *** TOUR INFO ***

int set_showpathnodes(int v); // SHOWPATHNODES
int set_showtourroute(int v); // SHOWTOURROUTE
// -- TOURCOLORS
// --  1.000000 0.000000 0.000000   :selected path line
// --  1.000000 0.000000 0.000000   :selected path line knots
// --  0.000000 1.000000 0.000000   :selected knot
// --  -1.000000 -1.000000 -1.000000   :path line
// --  -1.000000 -1.000000 -1.000000   :path knots
// --  -1.000000 -1.000000 -1.000000   :text
// --  1.000000 0.000000 0.000000   :avatar
int set_tourconstantvel(int v); // TOURCONSTANTVEL
int set_viewalltours(int v); // VIEWALLTOURS
int set_viewtimes(float a, float b, int c); // VIEWTIMES
int set_viewtourfrompath(int v); // VIEWTOURFROMPATH

// --  ------------ local ini settings ------------

int set_avatarevac(int v); // AVATAREVAC
// -- GEOMETRYTEST
// --  0 1 2.000000 10.000000
// --  1 1 1 1 1
// --  1 1 1 1 1
// --  106.000000 318.000000 50.000000 152.000000 4.000000 12.000000
// --  84.800003 39.799999 3.200000 339.200012 39.799999 3.200000
// --  212.000000 162.199997 3.200000 212.000000 101.000000 12.800000
// --  0.000000 0.000000 0.000000
int set_devicevectordimensions(float a, float b, float c, float d); // DEVICEVECTORDIMENSIONS
int set_devicebounds(float a, float b); // DEVICEBOUNDS
int set_deviceorientation(int a, float b); // DEVICEORIENTATION
// -- GRIDPARMS
// --  0 1 0
// --  398 225 32
// -- GSLICEPARMS
// --  0 0 0 0
// --  212.000000 101.000000 8.000000
// --  0.000000 90.000000
int set_loadfilesatstartup(int v); // LOADFILESATSTARTUP
int set_mscale(float a, float b, float c); // MSCALE
// -- SLICEAUTO
// --  9
// --  3
// --  8
// --  17
// --  26
// --  30
// --  42
// --  51
// --  81
// --  93
// -- MSLICEAUTO
// --  1
// --  9
int set_compressauto(int v); // COMPRESSAUTO
// -- PART5PROPDISP
// --   1  1  1  1  1
int set_part5color(int v); // PART5COLOR
// -- PART5CLASSVIS
// --  5
// --  1
// --  1
// --  1
// --  1
// --  1
// -- PROPINDEX
// --  3
// --  0 0
// --  1 0
// --  2 0
// -- SHOOTER
// --  0.500000 0.000000 0.018868
// --  0.250000 0.000000 0.000000
// --  0.000000 0.000000 0.000000
// --  1.000000 0.000000 4.000000
// --  10 1 100 0 0
// --  1.000000 1.000000
int set_showdevices(int n, char **names); // SHOWDEVICES
int set_showdevicevals(int showdeviceval, int showvdeviceval,
    int devicetypes_index, int colordeviceval, int vectortype, int vispilot,
    int showdevicetype, int showdeviceunit); // SHOWDEVICEVALS
int set_showmissingobjects(int v); // SHOWMISSINGOBJECTS
int set_tourindex(int v); // TOURINDEX

// -- USERTICKS
// --  0 1 5 1 1 1
// --  0.000000 -1.000000 0.000000
// --  0.000000 -1.000000 0.000000
// --  424.000000 203.000000 16.000000
// --  1.000000 1.000000 1.000000
// --  1 1 1
// -- XYZCLIP
// --  2
// --  0 -0.424000 0 424.424011
// --  0 94.497025 0 203.203995
// --  0 -0.016000 1 5.193398

// --  *** TIME/DATA BOUNDS ***
// --   (0/1 min max skip (1=set, 0=unset)

// -- C_PARTICLES
// --  0 1.000000 0 0.000000
// -- C_PARTICLES
// --  0 1.000000 0 0.000000 Uniform
// -- C_PLOT3D
// --  5
// --  1 0 1.000000 0 -0.000000
// --  2 0 1.000000 0 -0.000000
// --  3 0 1.000000 0 -0.000000
// --  4 0 1.000000 0 -0.000000
// --  5 0 1.000000 0 -0.000000
// -- C_SLICE
// --  0 1.000000 0 0.000000 X_CO
// -- C_SLICE
// --  0 1.000000 0 0.000000 temp
// -- C_SLICE
// --  0 1.000000 0 0.000000 VIS_C0.9H0.1
// -- CACHE_BOUNDARYDATA
// --  0
// -- CACHE_QDATA
// --  1
// -- PATCHDATAOUT
// --  0 1.000000 -1.000000 1.000000 -1.000000 1.000000 -1.000000 1.000000 -1.000000
// -- PERCENTILELEVEL
// --  0.010000
// -- TIMEOFFSET
// --  0.000000
// -- TLOAD
// --  0 0.000000 0 1.000000 0 0
// -- V_PARTICLES
// --  0 1.000000 0 0.000000
// -- V5_PARTICLES
// --  0 1.000000 0 0.000000 Uniform
// -- V_PLOT3D
// --  5
// --  1 0 1.000000 0 1.000000
// --  2 0 1.000000 0 1.000000
// --  3 0 1.000000 0 1.000000
// --  4 0 1.000000 0 1.000000
// --  5 0 1.000000 0 1.000000
// -- V_SLICE
// --  0 0.000000 0 0.000000 X_CO : 0.000000 1.000000 1
// -- V_SLICE
// --  0 1.000000 0 0.000000 temp : 0.000000 1.000000 1
// -- V_SLICE
// --  1 0.000000 1 20.000000 VIS_C0.9H0.1 : 0.000000 1.000000 1
// -- V_TARGET
// --  0 1.000000 0 0.000000
// -- VIEWPOINT5
// --  0 10 2
// --  0.490669 -2.257067 0.018868 1.000000 -2
// --  0.000000 0.000000 0.000000 1
// --  0.500000 0.240566 0.018868
// --  0.000000 90.000000
// --  1.000000 0.000000 0.000000 0.000000
// --  0.000000 1.000000 0.000000 0.000000
// --  0.000000 0.000000 1.000000 0.000000
// --  0.000000 0.000000 0.000000 1.000000
// --  2 0 0 0 0 0 1
// --  -0.424000 -1.204000 -0.016000 424.424011 203.203995 5.193398
// --  topDown