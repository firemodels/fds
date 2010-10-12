// $Date$ 
// $Revision$
// $Author$

#ifndef DEF_smokeheaders
#define DEF_smokeheaders

#include "isodefs.h"
#include "flowfiles.h"
#ifndef CPP
#include "egz_stdio.h"
#endif

#ifdef CPP
#define EXTERNCPP extern "C"
#else
#define EXTERNCPP
#endif

int SUB_portortho(int quad, 
                   GLint i_left, GLint i_down, GLsizei i_width, GLsizei i_height,
                   GLdouble x_left, GLdouble x_right, GLdouble x_bottom, GLdouble x_top,
                   GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height
                   );
int SUB_portfrustum(int quad, 
                   GLint i_left, GLint i_down, GLsizei i_width, GLsizei i_height,
                   GLdouble fleft, GLdouble fright, GLdouble fdown, GLdouble fup, GLdouble fnear, GLdouble ffar,
                   GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height
                   );

EXTERNCPP int last_slice_loadstack(void);
EXTERNCPP void push_slice_loadstack(int sliceindex);
EXTERNCPP void remove_slice_loadstack(int sliceindex);
EXTERNCPP int last_vslice_loadstack(void);
EXTERNCPP void push_vslice_loadstack(int sliceindex);
EXTERNCPP void remove_vslice_loadstack(int sliceindex);
EXTERNCPP int filecat(char *file_in1, char *file_in2, char *file_out);
EXTERNCPP void update_colorbar_smooth(void);
EXTERNCPP void update_tbounds(void);
EXTERNCPP void updateGluiTimeBounds(float time_min, float time_max);
EXTERNCPP void settimeval(float timeval);
EXTERNCPP int get_token_loc(char *var,sv_object_frame *frame);
EXTERNCPP void get_indep_var_indices(sv_object *smv_object,char **var_indep_strings, int nvars_indep,int *index);
EXTERNCPP void get_evac_indices(sv_object *smv_object, int *evac_index,int *nevac_index);
EXTERNCPP void update_glui_set_view_xyz(float *xyz);
EXTERNCPP void update_colorbar_list(void);
EXTERNCPP void update_colorbar_list2(void);

EXTERNCPP void update_glui_plot3dtype(void);
EXTERNCPP void update_glui_isotype(void);
EXTERNCPP void init_device(device *devicei, float *xyz, float *xyzn, int state0,int nparams, float *params, char *labelptr);
EXTERNCPP void init_device_plane(device *devicei);
EXTERNCPP void draw_devices_val(void);
EXTERNCPP void getsmokesensors(void);
EXTERNCPP float get_vecfactor(int *iveclengths);
EXTERNCPP void add_new_tour(void);
EXTERNCPP void cleanbuffer(char *buffer, char *buffer2);
EXTERNCPP void start_script(void);
EXTERNCPP int run_script(void);
EXTERNCPP int compile_script(char *scriptfile);
EXTERNCPP scriptfiledata *insert_scriptfile(char *file);
EXTERNCPP int file_exist(char *file);
EXTERNCPP char *get_inifilename(int id);
EXTERNCPP char *get_scriptfilename(int id);
EXTERNCPP inifiledata *insert_inifile(char *file);
EXTERNCPP void keyboard(unsigned char key, int x, int y);
EXTERNCPP void get_newscriptfilename(char *newscriptfilename);
EXTERNCPP void init_avatar(void);
EXTERNCPP void drawselect_avatars(void);
#ifdef pp_GPU
EXTERNCPP int log_base2(float xx);
#endif
EXTERNCPP void readterrain(char *file, int ifile, int flag, int *errorcode);
EXTERNCPP void initterrain_znode(mesh *meshi, terraindata *terri, float xmin, float xmax, int nx, float ymin, float ymax, int ny, 
                                 int allocate_memory);
EXTERNCPP void initterrain_all(void);
EXTERNCPP void update_terrain_colors(void);
EXTERNCPP void drawterrain(terraindata *terri, int only_geom);
EXTERNCPP void drawterrain_texture(terraindata *terri, int only_geom);
EXTERNCPP void drawtrees(void);
EXTERNCPP int createnulllabel(flowlabels *flowlabel);
#ifdef pp_CULL
EXTERNCPP void initcull(int cullflag);
EXTERNCPP void initcullplane(int cullflag);
EXTERNCPP void setPixelCount(void);
EXTERNCPP void setPixelCountOrthog(mesh *meshi);
EXTERNCPP void getPixelCount(void);
EXTERNCPP int init_cull_exts(void);
#endif
#ifdef pp_GPU
void getDepthTexture( void );
void createDepthTexture( void );
EXTERNCPP int init_shaders(void);
EXTERNCPP void LoadSmokeShaders(void);
EXTERNCPP void UnloadSmokeShaders(void);
#endif
EXTERNCPP void initspheresegs(int nlat, int nlong);
EXTERNCPP mesh *getmesh(float *xyz);
EXTERNCPP void update_glui_wui(void);
EXTERNCPP int have_terrain_slice(void);
EXTERNCPP void getfile_size(const char *filename, FILE_SIZE *filesize);
EXTERNCPP void getfile_modtime(char *filename, time_t *modtime);
EXTERNCPP float get_zcell_val_offset(mesh *meshi,float xval, float yval, int *loc);
EXTERNCPP void update_camera_ypos(camera *camera_data);
EXTERNCPP camera *get_camera(char *name);
EXTERNCPP char *get_camera_label(int index);
EXTERNCPP void clip2cam(camera *cam);
EXTERNCPP void cam2clip(camera *cam);
EXTERNCPP void to_lower(char *string);
EXTERNCPP void init_object_defs(void);
EXTERNCPP void update_device_textures(void);
EXTERNCPP char *get_device_label(char *buffer);
EXTERNCPP void get_elevaz(float *xyznorm,float *dtheta, float *rotate_axis, float *dpsi);
EXTERNCPP void drawTargetNorm(void);
EXTERNCPP void draw_devices(void);
EXTERNCPP void free_SVOBJECT(sv_object *object);
EXTERNCPP sv_object *parse_SVOBJECT(FILE *stream);
EXTERNCPP sv_object *init_SVOBJECT1(char *label, char *commands,int visible);
EXTERNCPP sv_object *init_SVOBJECT2(char *label, char *commandson, char *commandsoff,int visible);
EXTERNCPP sv_object *get_SVOBJECT_type(char *label, sv_object *default_object);
EXTERNCPP int read_object_defs(char *file);
EXTERNCPP void freeall_objects(void);
EXTERNCPP void parse_string(char *string, char **tokens, int *ntokens);
EXTERNCPP void update_partclass_depend(part5class *partclassi);

EXTERNCPP int get_plot3d_index(mesh *meshi, int dir, float val);
EXTERNCPP int plot3dlistcompare( const void *arg1, const void *arg2 );
EXTERNCPP int plot3dcompare( const void *arg1, const void *arg2 );
EXTERNCPP void update_plot_xyz(mesh *current_mesh);
EXTERNCPP void updateplotslice_mesh(mesh *mesh_in, int slicedir);

EXTERNCPP void printhrr(void);
EXTERNCPP void readhrr(int flag, int *errorcode);
EXTERNCPP char *get_chid(char *file);
EXTERNCPP void setColorbarClipPlanes(int flag);
EXTERNCPP void addcolorbar(int icolorbar);
EXTERNCPP void ReloadMenu(int value);
EXTERNCPP void ColorBarMenu(int val);
EXTERNCPP void initdefaultcolorbars(void);
EXTERNCPP void drawcolorbarpath(void);
EXTERNCPP void update_colorbar_splits(colorbardata *cbi);
EXTERNCPP void remapcolorbar(colorbardata *cbi);
EXTERNCPP void adjust_colorbar_splits(colorbardata *cbi);
EXTERNCPP void interpcolor(float *col1, float *col2,float *rgb,int npoints_seg, int jumpflag);
EXTERNCPP void freecolorbars(void);
EXTERNCPP void update_glui_stereo(void);
EXTERNCPP void escape_blanks(char *dirfrom, int maxlen);
EXTERNCPP void InitOpenGL(void);
EXTERNCPP void TextureShowMenu(int value);
EXTERNCPP void initcolors(void);
EXTERNCPP void copy_args(int *argc, char **aargv, char ***argv_sv);
EXTERNCPP void init_user_ticks(void);
EXTERNCPP void draw_user_ticks(void);
EXTERNCPP int get_tick_dir(float *mm);
EXTERNCPP void init_multi_threading(void);
#ifdef WIN32
EXTERNCPP void OpenSMVFile(char *filename,int filenamelength,int *openfile);
#endif
EXTERNCPP void get_smokezippath(char *progdir, char **zippath);
EXTERNCPP int AnySmoke(char *type);
EXTERNCPP int AnySlices(char *type);
EXTERNCPP void TrainerViewMenu(int var);

EXTERNCPP void delete_camera(camera *cam1);
EXTERNCPP void ShowAllSmoke(void);
EXTERNCPP void HideAllSmoke(void);
EXTERNCPP void HideAllSlices(void);
EXTERNCPP void ShowAllSlices(char *type1, char *type2);
EXTERNCPP void UnloadSliceMenu(int value);
EXTERNCPP void ViewpointMenu(int value);
EXTERNCPP void FrameRateMenu(int var);
EXTERNCPP void LoadUnloadMenu(int value);
EXTERNCPP void TourMenu(int var);
EXTERNCPP void ResetMenu(int var);
EXTERNCPP void LabelMenu(int value);
EXTERNCPP void ShadeMenu(int value);
EXTERNCPP void FontMenu(int value);
EXTERNCPP void ShowHideSliceMenu(int var);
EXTERNCPP void EvacShowMenu(int value);
EXTERNCPP void ParticleShowMenu(int value);
EXTERNCPP void Plot3DShowMenu(int value);
EXTERNCPP void IsoShowMenu(int value);
EXTERNCPP void ShowPatchMenu(int value);
EXTERNCPP void Smoke3DShowMenu(int value);
EXTERNCPP void ShowVSliceMenu(int value);
EXTERNCPP part5prop *get_part5prop_s(char *label);
EXTERNCPP int get_part5prop_index_s(char *shortlabel);
EXTERNCPP int get_part5prop_index(char *label);
EXTERNCPP void print_part5prop(void);
EXTERNCPP part5prop *get_part5prop(char *label);
EXTERNCPP void init_part5prop(void);
EXTERNCPP void update_streakvalue(float value);
EXTERNCPP void update_all_partvis2(void);
EXTERNCPP void ParticleMenu(int value);
EXTERNCPP void LoadPatchMenu(int value);
EXTERNCPP void LoadSliceMenu(int value);
EXTERNCPP void LoadVSliceMenu(int value);

EXTERNCPP void initvars1(void);
EXTERNCPP void initvars0(void);
EXTERNCPP void RenderState(int onoff);
EXTERNCPP void update_windowsizelist(void);
EXTERNCPP void ResizeWindow(int width, int height);
EXTERNCPP void update_glui_streakvalue(float rvalue);
EXTERNCPP void uncompress_slicedataframe(slice *sd,int iframe);
EXTERNCPP void glui_alert_setup(int main_window);
EXTERNCPP void show_load_alert(void);
EXTERNCPP void hide_load_alert(void);
EXTERNCPP void update_trainer_outline(void);
EXTERNCPP void update_trainer_moves(void);

#ifdef pp_MESSAGE
EXTERNCPP void message_message(char *message);
EXTERNCPP void warning_message(char *message);
EXTERNCPP void error_message(char *message);
EXTERNCPP void abort_message(char *message);
#endif
EXTERNCPP sv_object *get_object(char *label);
EXTERNCPP void get_labels(char *buffer, char **label1, char **label2);
EXTERNCPP void snap_view_angles(void);
#ifdef pp_SHOOTER
EXTERNCPP mesh *inmesh(float xyz[3]);
EXTERNCPP void get_plot3d_uvw(float xyz[3], float uvw[3]);
EXTERNCPP void solve_shooter_data(void);
EXTERNCPP void increment_shooter_data(shootpointdata *pold, shootpointdata *pnew, float dt);
EXTERNCPP void show_shooter(void);
EXTERNCPP void hide_shooter(void);
EXTERNCPP void draw_shooter(void);
#endif
EXTERNCPP void show_trainer(void);
EXTERNCPP void hide_trainer(void);
EXTERNCPP int get_trainee_location(void);
EXTERNCPP void set_trainer_controls(void);
EXTERNCPP void training_move(int  mode);
EXTERNCPP void load_startup_smoke(void);
EXTERNCPP void get_startup_vslice(int seq_id);
EXTERNCPP void get_startup_slice(int seq_id);
EXTERNCPP void get_startup_part(int seq_id);
EXTERNCPP void get_startup_plot3d(int seq_id);
EXTERNCPP void get_startup_smoke(int seq_id);
EXTERNCPP void get_startup_iso(int seq_id);
EXTERNCPP void get_startup_patch(int seq_id);
EXTERNCPP void set_3dsmoke_startup(void);
EXTERNCPP void clear_3dsmoke_startup(void);
EXTERNCPP void put_startup_smoke3d(FILE *fileout);
EXTERNCPP void setspeed(float tospeed);
EXTERNCPP void drawonlythreshold(const mesh *meshi);
EXTERNCPP void draw_transparent_faces(void);
EXTERNCPP smoothblockage *getsmoothblockage(mesh *meshi,float tt);
EXTERNCPP void freesmoothblocks(smoothblockage *sb);
EXTERNCPP int isblockagevisible(blockagedata *bc, float time);
EXTERNCPP float zoom2aperture(float zoom0);
EXTERNCPP float aperture2zoom(float ap);
EXTERNCPP int getZoneColor(float t, float tmin, float tmax, int nlevel);
EXTERNCPP void fill_roomdata(int izone);
EXTERNCPP void update_overwrite(void);
EXTERNCPP void compress_svzip(void);
EXTERNCPP char *getdir(char *argi);
EXTERNCPP char *which(char *prog);
EXTERNCPP char *lastname(char *argi);
EXTERNCPP void drawTargets(void);
EXTERNCPP void drawBlockages(int mode, int flag);
EXTERNCPP void drawLabels(void);
EXTERNCPP void Update_Tourlist(void);
EXTERNCPP void update_glui_speed(void);
EXTERNCPP void drawMovedir(void);
EXTERNCPP void getnewpos(float *oldpos, float dx, float dy, float dz, float speed_factor);
EXTERNCPP void free_skybox(void);
EXTERNCPP void draw_skybox(void);
EXTERNCPP void loadskytexture(char *filebase, texture *texti);
#ifdef USE_ZLIB
EXTERNCPP void uncompress_isodataframe(isosurface *asurface, int n);
EXTERNCPP void uncompress_patchdataframe(mesh *meshi,int iframe);
#endif
EXTERNCPP void getpatchdata_zlib(patch *patchi,unsigned char *data,int ndata, 
                       float *times, unsigned int *zipoffset, unsigned int *zipsize, int ntimes);
EXTERNCPP void getpatchsizeinfo(patch *patchi, int *nframes, int *buffersize);
EXTERNCPP void getpatchheader2(char *file, int *version, int *i1, int *i2, int *j1, int *j2, int *k1, int *k2, int *patchdir);
EXTERNCPP void getpatchheader(char *file,int *npatches,float *valmin, float *valmax);
EXTERNCPP int getsmoke3dversion(smoke3d *smoke3di);
EXTERNCPP void update_cadtextcoords(cadquad *quadi);
EXTERNCPP void open_smokepanel(void);
EXTERNCPP void open_smokezippanel(void);
EXTERNCPP void UpdateIndexColors(void);
EXTERNCPP void adjusttourtimes(tourdata *touri);
EXTERNCPP void update_tourindex(void);
EXTERNCPP void update_glui_zoom(void);
EXTERNCPP void SetTour(tourdata *thetour);
EXTERNCPP void freetrainers(void);
EXTERNCPP void update_plot3d_display(void);
EXTERNCPP void pauseSV(void);
EXTERNCPP void abortSV(char *message);
EXTERNCPP void ShellMenu(int var);
EXTERNCPP int getshellmenu_index(int menuid);
EXTERNCPP void makeshellmenus(int menuid,int flag);
EXTERNCPP void destroyshellmenus(void);
EXTERNCPP float cputime(void);
EXTERNCPP void update_smoke3dflags(void);
EXTERNCPP void mergesmoke3dcolors(void);
EXTERNCPP void setsmokecolorflags(void);
EXTERNCPP void sort_transparent_faces(float *mm);
EXTERNCPP void getsmokedir(float *mm);
EXTERNCPP void get_world_eyepos(float *mm, float eyepos[3]);
EXTERNCPP void ExtractFrustum(void);
EXTERNCPP int RectangleInFrustum( float *x11, float *x12, float *x22, float *x21);
EXTERNCPP unsigned char adjustalpha(unsigned char alpha, float factor);
EXTERNCPP unsigned int irle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
EXTERNCPP void updatesmoke3d(smoke3d *smoke3di);
EXTERNCPP void drawsmoke3d(smoke3d *smoke3di);
EXTERNCPP void TRANSLATE_CB(int var);
#ifdef pp_GPU
EXTERNCPP void drawsmoke3dGPU(smoke3d *smoke3di);
#endif
#ifdef pp_CULL
EXTERNCPP void drawsmoke3dCULL(void);
#endif
EXTERNCPP void get_drawing_parms(int *drawing_smooth, int *drawing_transparent, int *drawing_blockage_transparent, int *drawing_vent_transparent);
EXTERNCPP void updatesmoke3dmenulabels(void);
EXTERNCPP void Labels_CB(int value);
EXTERNCPP void bench_out(float frame_rate);
EXTERNCPP void output_Slicedata(void);
EXTERNCPP void init_Slicedata(void);
EXTERNCPP void update_camera_label(void);
EXTERNCPP void update_extreme(void);
EXTERNCPP void update_extreme2(void);
EXTERNCPP void update_colorbar_type(void);
EXTERNCPP void update_colorbar_label(void);
EXTERNCPP void init_camera_list(void);
EXTERNCPP camera *insert_camera(camera *cb,camera *source, char *name);
EXTERNCPP void add_default_views(void);
EXTERNCPP void update_view_gluilist(void);
EXTERNCPP void gluiIdle(void);
EXTERNCPP void gluiIdleNULL(void);
EXTERNCPP int interval_search(float *list, int nlist, float key, int guess);
EXTERNCPP void reset_gltime(void);
EXTERNCPP void updateshowtitles(void);
EXTERNCPP void enable_reset_saved_view(void);
EXTERNCPP void reset_glui_view(int ival);
EXTERNCPP void init_camera(camera *camera_data,char *name);
EXTERNCPP void copy_camera(camera *to, camera *from);
EXTERNCPP void set_camera_current(float angles[2], float eye[3], float zoom);
EXTERNCPP void update_camera(camera *ca);
EXTERNCPP void update_projection_type(void);
EXTERNCPP void update_eyerotate(void);
EXTERNCPP void update_cursor_checkbox(void);
EXTERNCPP void update_clip_all(void);
EXTERNCPP void backup_blockage(blockagedata *bc);
EXTERNCPP void update_blockpath(void);
EXTERNCPP void getinverse(float *m, float *mi);
EXTERNCPP void matmatmult(float *m1, float *m2, float *m3);
EXTERNCPP void matvecmult(double *m1, float *v1, float *v2);
EXTERNCPP void update_rotation_index(int val);
EXTERNCPP void update_meshlist1(int val);
EXTERNCPP void update_translate(void);
EXTERNCPP void BlockageMenu(int value);
EXTERNCPP void parsedatabase(char *databasefile);
EXTERNCPP char *STRSTR(char *c, const char *key);
EXTERNCPP void update_sorted_surfidlist(void);
EXTERNCPP void handle_plot3d_keys(int  key);
EXTERNCPP void handle_move_keys(int  key);
EXTERNCPP int get_interval(float val, float *array, int n);

EXTERNCPP void set_unit_vis(void);
EXTERNCPP int getrevision(char *svn);
EXTERNCPP void memorystatus(void);
EXTERNCPP int getnewfilename(void);
EXTERNCPP void showhide_translate(int var);
EXTERNCPP void updateallplotslices(void);
EXTERNCPP int makeiblank(void);
EXTERNCPP int makeiblank_carve(void);
EXTERNCPP void makeiblank_smoke3d(void);
EXTERNCPP void getunitinfo(const char *unitlabel, int *unitclass, int *unittype);
EXTERNCPP void update_unit_defs(void);

EXTERNCPP void SmoothIsoSurface(isosurface *surfacedata);
EXTERNCPP void updateslicefilenum(void);
EXTERNCPP void drawstaticiso(const isosurface *asurface,int surfacetype, 
                             int smoothnorms, int trans_flag, int data_type, 
                             float line_width);
EXTERNCPP void normalize(float *xyz, int n);
#ifndef CPP
EXTERNCPP void getisosizes(const char *isofile, int dataflag, EGZ_FILE **isostreamptr,
				 int *nvertices, int *ntriangles, float **levels, int *nisolevels,
				 int *nisosteps, int isoframestep, float *tmin, float *tmax, int endian);
#endif
EXTERNCPP void array2string(float *array, int narray, char *string);
EXTERNCPP void getisolevels(const char *isofile, int dataflag, float **levelsptr, int *nisolevels);
EXTERNCPP void getcisolevels(const char *isofile, float **levelsptr, int *nisolevels);
EXTERNCPP void getsmoothblockparms(mesh *gb, smoothblockage *sb);
EXTERNCPP void MakeIsoBlockages(mesh *gb, smoothblockage *sb);

EXTERNCPP int ifsmoothblock(void);
EXTERNCPP void updatevslices(void);
EXTERNCPP void updatepartmenulabels(void);
EXTERNCPP void updateisomenulabels(void);
EXTERNCPP void updatepatchmenulabels(void);
EXTERNCPP void updateslicemenulabels(void);
EXTERNCPP void updatevslicemenulabels(void);
EXTERNCPP void updateplot3dmenulabels(void);
EXTERNCPP void trimzeros(char *line);
EXTERNCPP void trimmzeros(char *line);
EXTERNCPP void handle_eyeview(int flag);

EXTERNCPP void checktimebound(void);
EXTERNCPP void getrgb(unsigned int val, unsigned char *rr, unsigned char *gg, unsigned char *bb);
EXTERNCPP unsigned char *readpicture(char *filename, int *width, int *height);
#ifdef pp_JPEG
EXTERNCPP unsigned char *readjpeg(const char *filename,int *width, int *height, int skip);
#endif
EXTERNCPP unsigned char *readrgb(const char *name, int *width, int *height);
EXTERNCPP unsigned char *readpng(const char *filename,int *width, int *height);

EXTERNCPP void update_whichface(int which_face);
EXTERNCPP void update_blockvals(int flag);
EXTERNCPP void show_glui_colorbar(void);
EXTERNCPP void hide_glui_colorbar(void);
EXTERNCPP void show_glui_motion(void);
EXTERNCPP void hide_glui_motion(void);
EXTERNCPP void show_glui_clip(void);
EXTERNCPP void update_glui_clip(void);

EXTERNCPP void hide_glui_clip(void);
EXTERNCPP void show_glui_wui(void);
EXTERNCPP void hide_glui_wui(void);
EXTERNCPP void show_glui_labels(void);
EXTERNCPP void set_labels_controls(void);
EXTERNCPP void hide_glui_labels(void);

EXTERNCPP void create_tourlist(void);
EXTERNCPP void delete_tourlist(void);
EXTERNCPP void updateviewtour(void);
EXTERNCPP void update_tourcontrols(void);
EXTERNCPP void adjustviewangle(keyframe *kf,float *azimuth, float *elevation);
EXTERNCPP void setup_tour(void);
EXTERNCPP void createtourpaths(void);
EXTERNCPP void drawtours(void);
EXTERNCPP void set_glui_keyframe(void);
EXTERNCPP void drawselect_tours(void);
EXTERNCPP void freetour(tourdata *touri);
EXTERNCPP void inittour(tourdata *touri);
EXTERNCPP void updatetourmenulabels(void);
EXTERNCPP void update_globaltension(void);
EXTERNCPP void defaulttour(void);
EXTERNCPP void show_glui_tour(void);
EXTERNCPP void hide_glui_tour(void);
EXTERNCPP void new_select(keyframe *newselect);
EXTERNCPP void delete_tour(int tour_index);
EXTERNCPP tourdata *add_tour(char *label);
EXTERNCPP void init_circulartour(void);
EXTERNCPP keyframe *delete_frame(keyframe *step);
EXTERNCPP void ReallocTourMemory(void);
EXTERNCPP keyframe *add_frame(keyframe *framei, float time, float *xyz, float key_azimuth, float elevation, float bank,
                    float params[3],int viewtype,float zoom,float view[3]);
EXTERNCPP void hide_glui_trainer(void);
EXTERNCPP void show_glui_stereo(void);
EXTERNCPP void hide_glui_stereo(void);
EXTERNCPP void show_glui_3dsmoke(void);
EXTERNCPP void hide_glui_3dsmoke(void);

EXTERNCPP void enable_boundary_glui(void);
EXTERNCPP void disable_boundary_glui(void);
EXTERNCPP void update_clipplanes(void);
EXTERNCPP void show_glui_bounds(void);
EXTERNCPP void hide_glui_bounds(void);

EXTERNCPP void get_blockvals(float *xmin, float *xmax,
                   float *ymin, float *ymax,
                   float *zmin, float *zmax,
                   int *imin, int *jmin, int *kmin);
EXTERNCPP void show_glui_edit(void);
EXTERNCPP void hide_glui_edit(void);

EXTERNCPP void updateslicelistindex(int sfn);
EXTERNCPP void updatepatchlistindex(int patchfilenum);
EXTERNCPP void updateplot3dlistindex(void);
EXTERNCPP int getindex(float key, const float *list, int nlist);
EXTERNCPP void transparentoff(void);
EXTERNCPP void transparenton(void);
EXTERNCPP void getlabels(const char *filein);
EXTERNCPP void initobst(blockagedata *bc,surface *surf,int index,int meshindex);
EXTERNCPP void setsurfaceindex(blockagedata *bc);
EXTERNCPP void updateusetextures(void);
EXTERNCPP int updatergbhist(int width, int height,int maketable, int *colortable, int nrgb);
EXTERNCPP void antialias(int flag);
EXTERNCPP void saveview(void);
EXTERNCPP void savelastview(void);
EXTERNCPP void setslicebounds(int islicetype);
EXTERNCPP void setisobounds(int islicetype);
EXTERNCPP void local2globalpatchbounds(const char *key);
EXTERNCPP void global2localpatchbounds(const char *key);
EXTERNCPP int getplotstate(int choice);
EXTERNCPP void update_loaded_lists(void);
EXTERNCPP void setsmokeviewvars(void);
EXTERNCPP void updateLights(int pos);
EXTERNCPP int mergescreenbuffers(GLubyte *screenbuffers[4]);
EXTERNCPP GLubyte *getscreenbuffer(void);
EXTERNCPP void ShowScene(int mode, int view_mode, int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
EXTERNCPP void updateShow(void);
EXTERNCPP int  inblockage(const mesh *gb,float x, float y, float z);
EXTERNCPP int inmesh_smoke(float x, float y, float z, int n, int flag);
EXTERNCPP mesh *get_mesh(float xyz[3]);
EXTERNCPP void initmesh(mesh *gb);
EXTERNCPP void updateglui(void);
EXTERNCPP void updateslicelist(int index);
EXTERNCPP void drawiso(const mesh *gb,int tranflag);
EXTERNCPP void drawtiso(const mesh *gb,int tranflag);
EXTERNCPP int getendian(void);
EXTERNCPP void drawplot3d(mesh *gb);
EXTERNCPP void drawplot3d_texture(mesh *gb);
EXTERNCPP void updateshowstep(int val, int slicedir);
EXTERNCPP void updateclip(int slicedir);
EXTERNCPP void updateclipbounds(int set_i0, int *i0, int set_i1, int *i1, int maxi);
EXTERNCPP void ClearBuffers(int mode);
EXTERNCPP void updatetimes(void);
EXTERNCPP void synctimes(void);
EXTERNCPP int compare( const void *arg1, const void *arg2 );
EXTERNCPP void updateplotslice(int slicedir);
EXTERNCPP void drawpatch(const mesh *gb);
EXTERNCPP void drawpatch_cellcenter(const mesh *gb);
EXTERNCPP void drawpatch_texture(const mesh *gb);
EXTERNCPP void drawpatch_texture_cellcenter(const mesh *gb);
EXTERNCPP void drawpatch_texture_threshold(const mesh *gb);
EXTERNCPP void drawpatch_threshold_cellcenter(const mesh *meshi);

EXTERNCPP void Render(int view_mode);
EXTERNCPP void updateslicebounds(void);
EXTERNCPP void updateisobounds(void);
EXTERNCPP void updateallslicecolors(int islicetype,int *errorcode);
EXTERNCPP void updateallisocolors(int iisotype,int *errorcode);
EXTERNCPP void updatevslicetypes(void);
EXTERNCPP int getvsliceindex(const vslice *vd);
EXTERNCPP int getvslicetype(const vslice *vd);
EXTERNCPP int getslicetype(const slice *sd);
EXTERNCPP int getslicetype_fromlabel(char *label);
EXTERNCPP void updateslicetypes(void);
EXTERNCPP int getsliceindex(const slice *sd);
EXTERNCPP void updatesliceboundlabels(void);
EXTERNCPP int getisotype(const iso *isoi);
EXTERNCPP int getisottype(const iso *isoi);
EXTERNCPP int getisoindex(const iso *isoi);
EXTERNCPP void update_isotype(void);
EXTERNCPP void updateisotypes(void);
EXTERNCPP int getpatchtype(const patch *patchi);
EXTERNCPP void update_patchtype(void);
EXTERNCPP void updatepatchtypes(void);

EXTERNCPP void update_mesh_terrain(void);
EXTERNCPP void update_terrain_options(void);
EXTERNCPP void update_glui_cellcenter(void);
EXTERNCPP void drawslice_cellcenter_interp(const slice *sd);
EXTERNCPP void update_plot3dtitle(void);
EXTERNCPP void LoadPlot3dMenu(int value);
EXTERNCPP void init_plot3dtimelist(void);
EXTERNCPP void update_framenumber(int changetime);
EXTERNCPP void update_iso_showlevels(void);
EXTERNCPP void update_current_mesh(mesh *meshi);
EXTERNCPP void DialogMenu(int value);
EXTERNCPP void ApertureMenu(int value);
EXTERNCPP void ZoomMenu(int value);
EXTERNCPP void setClipPlanes(int mode);
EXTERNCPP void unsetClipPlanes(void);
EXTERNCPP void setslicecolors(float slicemin, float slicemax, slice *sd, int *errorcode);
EXTERNCPP void setisocolors(float isomin, float isomax, iso *sd, int *errorcode);
EXTERNCPP void drawslice(const slice *sd);
EXTERNCPP void drawslice_cellcenter(const slice *sd);
EXTERNCPP void drawslice_texture(const slice *sd);
EXTERNCPP float MIN(float x, float y);
EXTERNCPP float MAX(float x, float y);
EXTERNCPP void drawslice_terrain(const slice *sd);
EXTERNCPP void drawvolslice_terrain(const slice *sd);
EXTERNCPP void drawvolslice_texture(const slice *sd);
EXTERNCPP int new_multi(slice *sdold,slice *sd);
EXTERNCPP int hide_slice2(slice *sdi,slice *sdj);
EXTERNCPP void drawvolslice(const slice *sd);
EXTERNCPP void drawvvolslice(const vslice *vd);
EXTERNCPP void drawvvolslice_terrain(const vslice *vd);
EXTERNCPP void drawvslice(const vslice *vd);
EXTERNCPP void drawvslice_terrain(const vslice *vd);
EXTERNCPP void drawTimeBar(void);
EXTERNCPP void drawColorBars(void);
EXTERNCPP void drawPart(const particle *parti);
EXTERNCPP void drawEvac(const particle *parti);
EXTERNCPP void drawStaticPart(const particle *parti);
EXTERNCPP void drawgrid(const mesh *gb);
EXTERNCPP void drawroomgeom(void);
EXTERNCPP void drawroomdata(void);
EXTERNCPP void drawventdata(void);
EXTERNCPP void Init(void);
EXTERNCPP void ResetView(int option);
EXTERNCPP void Reshape(int width, int height);
EXTERNCPP void UpdateTimeLabels(void);
EXTERNCPP void Idle(void);
#define IDLE() Idle();
EXTERNCPP void RenderFrame(int view_mode);
EXTERNCPP int readlabels(flowlabels *label, FILE *stream);
EXTERNCPP int readlabels_terrain(flowlabels *label, FILE *stream);
EXTERNCPP int readlabels_cellcenter(flowlabels *label, FILE *stream);
EXTERNCPP void update_terrain(int allocate_memory, float vertical_factor);
EXTERNCPP void PART_CB_INIT(void);
EXTERNCPP void Slice_CB(int var);
EXTERNCPP void RenderMenu(int value);
EXTERNCPP void LoadSmoke3DMenu(int value);
EXTERNCPP void Display(void);
EXTERNCPP void Visible(int state);
EXTERNCPP void Args(int argc, char **argv);
EXTERNCPP void usage(char **argv);
EXTERNCPP void version(void);
EXTERNCPP int getmaxrevision(void);
EXTERNCPP void draw_demo(int nlat, int nlong);
EXTERNCPP void draw_demo2(int option);
EXTERNCPP void init_demo(float rad, int nlat, int nlong);
EXTERNCPP void drawoutlines(void);
EXTERNCPP void drawcbox(float x, float y, float z, float size);
EXTERNCPP void specialkeyboard(int key, int x, int y);
EXTERNCPP void specialkeyboard_up(int key, int x, int y);
EXTERNCPP void handleiso(void);
EXTERNCPP void keyboard(unsigned char key, int x, int y);
EXTERNCPP void keyboard_up(unsigned char key, int x, int y);
EXTERNCPP void togglegridstate(int visg);
EXTERNCPP void updatesurface(void);
EXTERNCPP void WindowStatus(int state);
EXTERNCPP void mouse(int button, int state, int x, int y);
EXTERNCPP void motion(int xm, int ym);
EXTERNCPP void nodein_extvent(
                    int ipatch, 
                    int *patchblankcopy,const mesh *meshi,int i1,int i2, int j1, int j2, int k1, int k2);
EXTERNCPP void setventdirs(void);
EXTERNCPP int nodeinblockage(const mesh *meshi, int i,int j,int k, int *imesh, int *iblockage);
EXTERNCPP int nodeinvent(const mesh *gb, int i,int j,int k, int dir);
EXTERNCPP void MenuStatus(int status, int x, int y);
EXTERNCPP int readini2(char *inifile, int localfile);
EXTERNCPP int match(char *buffer, const char *key, unsigned int lenkey);
EXTERNCPP int matchonly(char *buffer, const char *key, unsigned int lenkey);
EXTERNCPP int match_upper(char *buffer, const char *key, unsigned int lenkey);
EXTERNCPP void obst_or_vent2faces(const mesh *gb,blockagedata *bc, ventdata *vi, facedata *faceptr,int facetype);
EXTERNCPP void initsurface(surface *surf);
EXTERNCPP void initventsurface(surface *surf);
EXTERNCPP void update_hidden_faces(void);
EXTERNCPP void update_selectfaces(void);
EXTERNCPP void update_selectblocks(void);
EXTERNCPP void draw_faces(void);
EXTERNCPP void drawselect_faces(void);
EXTERNCPP void allocate_faces(void);
EXTERNCPP void update_facelists(void);
EXTERNCPP void update_faces(void);
EXTERNCPP char *trim_front(char *line);
EXTERNCPP char *trim_string(char *buffer);
EXTERNCPP char *trim_front_back(char *line);
EXTERNCPP void trim(char *line);
EXTERNCPP void drawticks(void);
EXTERNCPP void set_startup_view(void);
EXTERNCPP void add_list_view(char *label_in);
EXTERNCPP float color2bw(const float *color);
EXTERNCPP float *getcolorptr(const float *color);
EXTERNCPP void colorconvert(int flag);
EXTERNCPP void initcadcolors(void);
EXTERNCPP void updatecolors(int colorindex);
EXTERNCPP void initrgb(void);
EXTERNCPP void updatechopcolors(void);
EXTERNCPP void freelabels(flowlabels *label);
EXTERNCPP int readini(int scriptconfigfile);
EXTERNCPP void writeini(int flag);
EXTERNCPP void DrawCone(float radius, float height);
EXTERNCPP int ispatchtype(int type);
EXTERNCPP void adjustdatabounds(const float *pdata, int skip, int ndata, int setpmin, float *pmin, int setpmax, float *pmax);
EXTERNCPP void adjustpartbounds(const float *pdata, int particle_type, int droplet_type, const unsigned char *isprink, 
                      int skip, int ndata, int setpmin, float *pmin, int setpmax, float *pmax);
EXTERNCPP void adjustpart5bounds(particle *parti);
EXTERNCPP void adjustPlot3Dbounds(int iplot3d, int setpmin, float *pmin, int setpmax, float *pmax);
EXTERNCPP void adjustslicebounds(const slice *sd, float *pmin, float *pmax);
EXTERNCPP void getslicedatabounds(const slice *sd, float *pmin, float *pmax);
EXTERNCPP void scalefloat2string(float floatfrom, char *stringto, const float *scale, float range);
EXTERNCPP void scalestring(const char *stringfrom, char *stringto, const float *scale, float range);
EXTERNCPP void num2string(char *string, float tval, float range);
EXTERNCPP float frexp10(float x, int *exp10);
EXTERNCPP int initcase_c(int argc, char **argv);


EXTERNCPP void freecadinfo(void);

EXTERNCPP void init_unit_defs(void);
EXTERNCPP void InitUnits(void);
EXTERNCPP void readcad2geom(cadgeom *cd);
EXTERNCPP void readcadgeom(cadgeom *cd);
EXTERNCPP void drawcadgeom(const cadgeom *cd);
EXTERNCPP void drawcad2geom_opaque(const cadgeom *cd,int trans_flag);

EXTERNCPP char *newtextptr(char ***texture_list,int *n_texture_list,char *texturebuffer,char *lastbuffer);

EXTERNCPP void readplot3d(char *file, int ifile, int flag,int *errorcode);
EXTERNCPP void readpatch(int ifile, int flag, int *errorcode);
EXTERNCPP void readpart(char *file, int ifile, int flag, int *errorcode);
EXTERNCPP void readzone(char *file, int ifile, int flag, int *errorcode);
EXTERNCPP void readvslice(int ivslice, int flag, int *errorcode);

EXTERNCPP void smooth_blockages(void);
EXTERNCPP void freesmoke3d(smoke3d *smoke3di);
EXTERNCPP void readsmoke(int ifile,int flag, int *errorcode);
EXTERNCPP void readsmoke3d(int ifile,int flag, int *errorcode);
EXTERNCPP int getsmoke3d_sizes(char *smokefile, int version, 
                      float **timelist, int **use_smokeframe,
                      int *nchars_uncompressed, 
                      int **nchars_compressed,
                      int **nchars_compressed_full,
                      int *nframes, int *nframes_full);
EXTERNCPP void readslice(char *file, int ifile, int flag, int *errorcode);
EXTERNCPP void readtarget(const char *file, int ifile, int flag, int *errorcode);
EXTERNCPP void update_smooth_blockages(void);
EXTERNCPP void readtarget2(const char *file, int ifile, int flag, int *errorcode);
EXTERNCPP void readiso(const char *file, int ifile, int flag, int *errorcode);

EXTERNCPP void InitMenus(int unload);
EXTERNCPP void smoothlabel(float *min, float *max, int n);
EXTERNCPP int readsmv(char *file, char *file2);
EXTERNCPP int STRCMP(const char *s1, const char *s2);
EXTERNCPP int gettargetposition(int itarget, float time, float *x, float *y, float *z);
EXTERNCPP void outputAxisLabels(void);
EXTERNCPP void outputLargeText(float x, float y, const char *string);
EXTERNCPP void outputText(float x, float y, const char *string);
EXTERNCPP void output3Text(float *color, float x, float y, float z, const char *string);
EXTERNCPP void outputBarText(float x, float y, const GLfloat *color, const char *string);
EXTERNCPP void getzonebounds(const float *pdata, int ndata, int setpmin, float *pmin, int setpmax, float *pmax);
EXTERNCPP void updatechar(void);
EXTERNCPP void updatetracers(void);
EXTERNCPP void getPart5Colors(particle *parti, int nlevels);
EXTERNCPP void getPartColors(const float *t, int skip, int nt, unsigned char *it,
                   const unsigned char *isprink, int particle_type, int droplet_type,
              const float *tmin, const float *tmax, int nlevel,
              char **labels, char *scale, float *partlevels256);
EXTERNCPP void getBoundaryColors(float *t, int nt, unsigned char *it, 
              int settmin, float *tmin, int settmax, float *tmax, 
              float *tmin_global, float *tmax_global,
              int ndatalevel, int nlevel,
              char **labels, char *scale, float *tvals256,
              int *extreme_min, int *extreme_max);
EXTERNCPP void getBoundaryColors2(float *t, int nt, unsigned char *it, 
              int settmin, float *ttmin, int settmax, float *ttmax,
              float *tmin_global, float *tmax_global,
              int ndatalevel,
              int *extreme_min, int *extreme_max
              );
EXTERNCPP void getBoundaryColors3(patch *patchi, float *t, int nt, unsigned char *it, 
              int settmin, float *tmin, int settmax, float *tmax, 
              float *tmin_global, float *tmax_global,
              int ndatalevel, int nlevel,
              char **labels, char *scale, float *tvals256,
              int *extreme_min, int *extreme_max);
EXTERNCPP void getBoundaryLabels(
              float tmin, float tmax,
              char **labels, char *scale, float *tvals256, int nlevel);
EXTERNCPP void getZoneColors(const float *t, int nt, unsigned char *it,
               const float *tmin, const float *tmax, int nlevel, int nlevel_full,
               char **labels, char *scale, float *tvals256
               );
EXTERNCPP void getPlot3DColors(int iplot, int settmin, float *ttmin, int settmax, float *ttmax, 
              int ndatalevel, int nlevel,
              char **labels,char **labelsiso, char **scale, float *tlevels, float *tlevels256,
              int *extreme_min, int *extreme_max
              );
EXTERNCPP float getsliceval(slice *sd, unsigned char ival);
EXTERNCPP void updateallslicelabels(int slicetype, int *errorcode);
EXTERNCPP void updateallisolabels(int slicetype, int *errorcode);
EXTERNCPP void setslicelabels(float smin, float smax, 
                    slice *sd, int *errorcode);
EXTERNCPP void getSliceLabels(float tmin, float tmax, int nlevel,
              char labels[12][11],char **scale, float *tlevels256);
EXTERNCPP void updatePart5extremes(void);
EXTERNCPP void getSliceColors(const float *t, int nt, unsigned char *it,
              float tmin, float tmax, 
              int ndatalevel, int nlevel,
              char labels[12][11],char **scale, float *tlevels2,
              int *extreme_min, int *extreme_max
              );
EXTERNCPP void setisolabels(float smin, float smax, 
                    iso *sd, int *errorcode);
EXTERNCPP void getIsoLabels(float tmin, float tmax, int nlevel,
              char labels[12][11],char **scale, float *tlevels256);
EXTERNCPP int SVimage2file(char *GIFfilename, int rendertype, int width, int height);
EXTERNCPP void update_showhidebuttons(void);
EXTERNCPP void update_fileload(void);

#ifndef CPP
#include "smokefortheaders.h"
#endif

#endif
