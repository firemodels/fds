#ifndef FLOWFILES_H_DEFINED
#define FLOWFILES_H_DEFINED

/* --------------------------  circdata ------------------------------------ */

typedef struct {
  float *xcirc, *ycirc;
  int ncirc;
} circdata;

/* --------------------------  langlistdata ------------------------------------ */

typedef struct {
  char *file;
  char lang_code[3];
  char lang_name[32];
} langlistdata;

/* --------------------------  procdata ------------------------------------ */
#ifdef CPP
typedef struct {
  GLUI_Rollout *rollout;
  int rollout_id;
} procdata;
#endif
/* --------------------------  csvdata ------------------------------------ */

typedef struct {
  char *file;
  int loaded, display;
  int type;
} csvdata;

/* --------------------------  point ------------------------------------ */

typedef struct {
  float xyz[3],point_norm[3],texture_xy[3];
  int itriangle,ntriangles,nused;
  unsigned char on_mesh_boundary;
  struct _triangle **triangles;
} point;

/* --------------------------  triangle ------------------------------------ */

typedef struct _triangle {
  unsigned char skinny;
  float distance, *color, tpoints[6], tri_norm[3];
  struct _texturedata *textureinfo;
  struct _surfdata *surf;
  int vert_index[3], interior;
  point *points[3];
} triangle;

/* --------------------------  tetrahedron ------------------------------------ */

typedef struct _tetrahedron {
  float distance, *color, face_norm[4];
  int vert_index[4],exterior[4],faces[12];
  struct _matldata *matl;
  point *points[4];
} tetrahedron;

/* --------------------------  geomlistdata ------------------------------------ */

typedef struct {
  int npoints,ntriangles,nvolus;
  point *points;
  triangle *triangles;
  tetrahedron *volumes;
} geomlistdata;

/* --------------------------  geomobjdata ------------------------------------ */

typedef struct {
  struct _surfdata *surf;
  struct _texturedata *texture;
  char *texture_name;
  float texture_width, texture_height, texture_center[3];
  int texture_mapping;
} geomobjdata;

/* --------------------------  geomdata ------------------------------------ */

typedef struct {
  char *file;
  int memory_id;
  int loaded, display;
  struct _surfdata *surf;
  geomlistdata *geomlistinfo,*geomlistinfo_0, *currentframe;
  float *float_vals;
  int *int_vals, nfloat_vals, nint_vals;
  float *times;
  int ntimes,itime,*timeslist;
  int ngeomobjinfo;
  geomobjdata *geomobjinfo;
} geomdata;

/* --------------------------  bounddata ------------------------------------ */

typedef struct {
  float percentile_min, percentile_max;
  float global_min, global_max;
  int defined;
} bounddata;

/* --------------------------  propdata ------------------------------------ */
#define PROPVARMAX 100
typedef struct {
  char *label;
  int menu_id;
  float *rotate_axis, rotate_angle;
  int inblockage,blockvis;
  int nsmokeview_ids,ismokeview_ids;
  char *smokeview_id,**smokeview_ids;
  struct _sv_object *smv_object, **smv_objects;
  int ntextures;
  char **texturefiles, **vars_indep, **svals;
  int *vars_indep_index, vars_dep_index[PROPVARMAX], fvars_evac_index[PROPVARMAX];
  int nvars_indep,      nvars_dep,                   nvars_evac;
  float *fvals, fvars_evac[PROPVARMAX], fvars_dep[PROPVARMAX];
  int draw_evac;
  int tag_number;
} propdata;

/* --------------------------  shootpointdata ------------------------------------ */

typedef struct _shootpointdata {
  struct _shootpointdata *prev;
  int visible;
  float xyz[3], uvw[3], uvw_air[3], val;
} shootpointdata;

/* --------------------------  shoottimedata ------------------------------------ */

typedef struct {
  float time;
  int frame;
  shootpointdata *beg, *end;
} shoottimedata;

/* --------------------------  infiledata ------------------------------------ */

typedef struct _inifiledata {
  struct _inifiledata *prev, *next;
  int id;
  char *file;
} inifiledata;

/* --------------------------  scriptfiledata ------------------------------------ */

typedef struct _scriptfiledata {
  struct _scriptfiledata *prev, *next;
  int id;
  int recording;
  char *file;
} scriptfiledata;

/* --------------------------  scriptdata ------------------------------------ */

typedef struct {
  int command;
  char command_label[32];
  int ival,ival2,ival3,ival4,ival5;
  char *cval,*cval2;
  float fval,fval2,fval3;
  int exit,first,remove_frame;
} scriptdata;

/* --------------------------  treedata ------------------------------------ */

typedef struct {
  float xyz[3];
  float time_char, time_complete;
  float trunk_diam;
  float tree_height;
  float base_diam;
  float base_height;
  int state;
} treedata;

/* --------------------------  colorbardata ------------------------------------ */

typedef struct {
  char label[1024], *label_ptr ;        // menu label
  int nnodes,nodehilight,nsplits;
  unsigned char rgb_node[3*256];
  unsigned char alpha[256];
  unsigned char index_node[256];  // colorbar index
  unsigned char splits[256];
  float colorbar[3*256];
} colorbardata;

/* --------------------------  colortabledata ------------------------------------ */

typedef struct {
  char label[1024];
  int color[4];
} colortabledata;

/* --------------------------  surfid ------------------------------------ */

typedef struct {
  char *label;
  int show;
  int location;
} surfid;

/* --------------------------  colordata ------------------------------------ */

typedef struct _colordata {
  float color[4], full_color[4], bw_color[4];
  struct _colordata *nextcolor;
} colordata;

/* --------------------------  outline ------------------------------------ */

typedef struct {
  int nlines;
  float *x1, *y1, *z1;
  float *x2, *y2, *z2;
} outline;

/* --------------------------  labeldata ------------------------------------ */

typedef struct _labeldata {
  struct _labeldata *prev, *next;
  char name[300];
  float xyz[3],frgb[4],tstart_stop[2];
  int rgb[4], glui_id, labeltype; // smv or ini
  int useforegroundcolor,show_always;
} labeldata;


/* --------------------------  texture ------------------------------------ */

typedef struct _texturedata {
  char *file;
  int loaded, display, used;
  GLuint name;
} texturedata;

/* --------------------------  terraincell ------------------------------------ */

typedef struct {
  int nallocated, nstates;
  float *time;
  int interval;
  unsigned char *state;
} terraincell;

/* --------------------------  terraindata ------------------------------------ */

typedef struct {
  char *file;
  unsigned char *state;
  int *timeslist;
  int loaded, display;
  int autoload;
  texturedata *ter_texture;
  int nx, ny;
  float xmin, xmax, ymin, ymax;
  float *x, *y;
  float levels[13];
  float *zcell, *znode, *znode_scaled, *znode_offset;
  unsigned char *uc_znormal;
  float *times;
  terraincell *tcell;
  struct _mesh *terrain_mesh;
  int ntimes;
} terraindata;

/* --------------------------  matldata ------------------------------------ */

typedef struct _matldata {
  char *matllabel;
  float *color;
} matldata;

/* --------------------------  surfdata ------------------------------------ */

typedef struct _surfdata {
  char *surfacelabel,*texturefile;
  int type; /* 
               0 - regular block non-textured 
               1 - regular block textured
               2 - outline
               3 - smoothed block
               4 - invisible
             */
  float *color, emis, temp_ignition;
  float transparent_level;
  int iso_level;
  float t_width, t_height;
  texturedata *textureinfo;
  int obst_surface;
  int invisible;
  int location;
  int transparent;
  int used_by_obst,used_by_vent;
} surfdata;

/* --------------------------  facedata ------------------------------------ */

typedef struct {
  int type,type2;
  float approx_vertex_coords[12];
  float exact_vertex_coords[12];
  float approx_texture_coords[8];
  float exact_texture_coords[8];
  float *texture_origin;
  float normal[3];
  float *color;
  float *linecolor,*linewidth;
  float approx_center_coord[3];
  float dist2eye;
  int meshindex, blockageindex;
  int imin, imax, jmin, jmax, kmin, kmax;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int dir,hidden,dup;
  int del;
  int invisible;
  int transparent;
  int patchpresent;
  struct _culldata *cullport;
  int **showtimelist_handle;
  int thinface;
  int show_bothsides, is_interior;
  struct _blockagedata *bc;
  surfdata *surfinfo;
  texturedata *textureinfo;
} facedata;

/* --------------------------  selectdata ------------------------------------ */

typedef struct {
  int mesh, blockage, side,dir;
  facedata *facei;
  int type;
} selectdata;

/* -------------------------- blockagedata ------------------------------------ */

typedef struct _blockagedata {
  int ijk[6],ijkORIG[6];
  float xmin, xmax, ymin, ymax, zmin, zmax, xyzORIG[6];
  float xyzEXACT[6];
  surfdata *surf[6],*surfORIG[6];
  propdata *prop;
  int walltype,walltypeORIG;
  int surf_index[6],surf_indexORIG[6];
  int patchvis[7];
  int usecolorindex;
  int blockage_id,dup;
  int is_wuiblock;
  int hole;
  int nnodes;
  int hidden,invisible;
  int transparent;
  int meshindex;
  int del;
  int changed,changed_surface;
  int type;
  float *showtime;
  int *showtimelist;
  unsigned char *showhide;
  int nshowtime, show;
  char *label;
  float *color;
  int colorindex;
  int useblockcolor;
  facedata *faceinfo[6];
  float texture_origin[3];
} blockagedata;

/* --------------------------  cadlook ------------------------------------ */

typedef struct {
  int index;
  float texture_width, texture_height, texture_origin[3];
  float rgb[4], shininess;
  texturedata textureinfo;
  int onesided;
} cadlook;

/* --------------------------  cadquad ------------------------------------ */

typedef struct {
  float xyzpoints[12];
  float txypoints[8];
  float normals[3];
  int colorindex;
  float colors[4];
  float time_show;
  cadlook *cadlookq;
} cadquad;

/* --------------------------  clipdata ------------------------------------ */

typedef struct {
  int option;
  GLdouble clipvals[24];
  int clip_xmin, clip_xmax;
  int clip_ymin, clip_ymax;
  int clip_zmin, clip_zmax;
  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;
} clipdata;

/* --------------------------  cadgeom ------------------------------------ */

typedef struct {
  char *file;
  int *order;
  int version;
  int ncadlookinfo;
  cadlook *cadlookinfo;
  int nquads;
  cadquad *quad;
} cadgeom;

/* --------------------------  cventdata ------------------------------------ */

typedef struct _cventdata {
  int dir,type,colorindex,cvent_id,isOpenvent;
  float boxmin[3], boxmax[3], texture_origin[3];
  float xmin, xmax, ymin, ymax, zmin, zmax;
  unsigned char *blank;
  int   imin, imax, jmin, jmax, kmin, kmax;
  int useventcolor,hideboundary;
  float origin[3], radius;
  float *color;
  surfdata *surf[1];
  texturedata *textureinfo[1];
} cventdata;

/* --------------------------  ventdata ------------------------------------ */

typedef struct _ventdata {
  int type,dummy;
  struct _ventdata *dummyptr;
  int hideboundary;
  int dir,dir2,vent_id;
  int useventcolor;
  int isOpenvent;
  float xvent1, xvent2;
  float yvent1, yvent2;
  float zvent1, zvent2;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int imin, imax, jmin, jmax, kmin, kmax;
  float xvent1plot, xvent2plot;
  float yvent1plot, yvent2plot;
  float zvent1plot, zvent2plot;
  float *showtime;
  int *showtimelist;
  unsigned char *showhide;
  int nshowtime;
  float *color,*color_bak;
  int transparent;
  int colorindex;
  int usecolorindex;
  surfdata *surf[1];
  texturedata *textureinfo[1];
  float texture_origin[3];
  float *linewidth;
} ventdata;

/* --------------------------  tickdata ------------------------------------ */

typedef struct {
  float begin[3],end[3],length;
  float dxyz[3],dlength;
  int dir,nbars,useforegroundcolor;
  float width, rgb[3];
} tickdata;


/* --------------------------  feddata ------------------------------------ */

typedef struct {
  struct _slicedata *co,*co2,*o2,*fed_slice;
  struct _isodata *fed_iso;
  int co_index, co2_index, o2_index, fed_index;
  int loaded,display;
} feddata;

/* --------------------------  iso ------------------------------------ */

typedef struct _isodata {
  int seq_id, autoload;
  int isof_index;
  char *reg_file, *size_file;
  short *normaltable;
  int memory_id;
  int nnormaltable; 
  char *file,*tfile;
  int dataflag,geomflag;
  int is_fed;
  feddata *fedptr;
  int type;
  int num_memblocks;
  int setvalmin, setvalmax;
  float valmin, valmax;
  int firstshort;
  flowlabels surface_label, color_label;
  geomdata *geominfo;
  int blocknumber,display,loaded;
  float tmin,tmax;
  float valmin_data, valmax_data;
  int extreme_min, extreme_max;
  int isoupdate_timestep;
  float *levels, **colorlevels;
  int nlevels;
  char menulabel[128];
} isodata;

/* --------------------------  smoothblockage ------------------------------------ */

typedef struct {
  int nsmoothblockagecolors;
  float *smoothblockagecolors;
  isosurface **smoothblockagesurfaces;
  float time;
} smoothblockage;

/* --------------------------  volrenderdata ------------------------------------ */

typedef struct _volrenderdata {
  char *rendermeshlabel;
  struct _slicedata *smokeslice, *fireslice;
  int is_compressed;
  unsigned char *c_smokedata_view, *c_firedata_view;
  float *smokedata_full, *firedata_full;
  float *smokedata_view, *firedata_view;
  LINT *firepos, *smokepos;
  void *smokedataptr, *firedataptr;
  void **smokedataptrs, **firedataptrs;
  int *nsmokedata_compressed, *nfiredata_compressed;
  float *times;
  int *dataready;
  int itime, ntimes, times_defined;
  int *timeslist;
  float *smokecolor_yz0, *smokecolor_xz0, *smokecolor_xy0;
  float *smokecolor_yz1, *smokecolor_xz1, *smokecolor_xy1;
  int loaded, display;
} volrenderdata;

/* --------------------------  mesh ------------------------------------ */

typedef struct _mesh {
  int ibar, jbar, kbar;
  float cellsize;
  int ncvents,nvents,ndummyvents;
  int nbptrs;
  int is_bottom;

  int *cutcells, ncutcells;
  int update_firehalfdepth;
  terraindata *terrain;
  int mesh_type;
#ifdef pp_GPU
  GLuint smoke_texture_id,fire_texture_id,blockage_texture_id;
  float *smoke_texture_buffer,*fire_texture_buffer;
  GLuint slice3d_texture_id;
  float *slice3d_texture_buffer,*slice3d_c_buffer;
#endif
  float meshrgb[3], *meshrgb_ptr;
  float mesh_offset[3], *mesh_offset_ptr;
  int blockvis;
  float *xplt, *yplt, *zplt;
  int ivolbar, jvolbar, kvolbar;
  float *xvolplt, *yvolplt, *zvolplt;
  float *xplt_cen, *yplt_cen, *zplt_cen;
  float *xplt_orig, *yplt_orig, *zplt_orig;
  float x0, x1, y0, y1, z0, z1;
  int drawsides[7];
  int extsides[7]; // 1 if on exterior side of a supermesh, 0 otherwise
  int inside;
  float boxmin[3], boxmax[3], dbox[3], boxeps[3], dcell, dcell3[3];
  float slice_min[3], slice_max[3];
  float boxmin_scaled[3], boxmax_scaled[3];
  float *zcell;
  float xyz_bar0[3], xyz_bar[3];
  float xcen, ycen, zcen;
  float face_centers[18];
  float offset[3];
  float xyzmaxdiff;
  float boxoffset;
  int plot3dfilenum,isofilenum,patchfilenum;
  int *iplotx_all, *iploty_all, *iplotz_all;
  int plotx, ploty, plotz;
  int slicedir;
  int plotn;
  char *c_iblank_node,*c_iblank_cell,*c_iblank_x,*c_iblank_y,*c_iblank_z;
  float *f_iblank_cell;
  char *c_iblank_embed;
  float *block_zdist;
  int zdist_flag;
  unsigned char *iblank_smoke3d;
  int iblank_smoke3d_defined;
  blockagedata **blockageinfoptrs;
  int *obst_bysize;
  ventdata *ventinfo;
  cventdata *cventinfo;
  unsigned char *is_block_terrain;
  unsigned char *iqdata;
  float *qdata, *udata, *vdata, *wdata;
  unsigned char *yzcolorbase, *xzcolorbase, *xycolorbase; 
  float *yzcolorfbase, *xzcolorfbase, *xycolorfbase;
  float *yzcolortbase, *xzcolortbase, *xycolortbase;
  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  char *c_iblank_xy, *c_iblank_xz, *c_iblank_yz;
  float plot3d_speedmax;
  contour plot3dcontour1,plot3dcontour2,plot3dcontour3;
  contour terrain_contour;
  isosurface currentsurf,currentsurf2;
  isosurface *blockagesurface;
  isosurface **blockagesurfaces;
  float *smoothblockagecolors;
  int nsmoothblockagecolors;
  int ntc;
  int nspr;
  float *xsprplot, *ysprplot, *zsprplot, *tspr;
  int nheat;
  float *xheatplot, *yheatplot, *zheatplot, *theat;
  float *xspr, *yspr, *zspr;
  float *xheat, *yheat, *zheat;

  isosurface *animatedsurfaces;
  int nisolevels, *showlevels;
  float *isolevels;
  int isomin_index, isomax_index;
  int niso_times;
  float *iso_times;
  int *iso_timeslist;
  int iso_itime;
  int smokedir,smokedir_old;
  float dx, dy, dz, dxy,dxz,dyz;
  float norm[3];

  int *patchtype;
  int *patchdir,*patch_surfindex;
  int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2;
  contour **patch_contours;
  int *blockonpatch;
  struct _mesh **meshonpatch;
  struct _mesh *nabors[6];
  struct _supermesh *super;
  int *ptype;
  int *patchrow, *patchcol, *blockstart;
  unsigned int *zipoffset, *zipsize;
  int *visPatches;
  float *xyzpatch, *xyzpatch_threshold;
  unsigned char *cpatchval_zlib, *cpatchval_iframe_zlib;
  unsigned char *cpatchval, *cpatchval_iframe;
  float *patch_times, *patch_timesi, *patchval, *patchval_iframe;
  float **patchventcolors;
  float *thresholdtime;
  int *patchblank;
  int npatch_times,npatches;
  int patch_itime;
  int *patch_timeslist;
  int npatchsize;
  int visInteriorPatches;
  float surface_tempmin, surface_tempmax;

  smoothblockage *smoothblockages_list;
  int nsmoothblockages_list;

  smoothblockage **showsmoothtimelist;

  int nface_textures, nface_outlines, nfaces;
  int nface_normals_single, nface_normals_double, nface_transparent_double;
  facedata *faceinfo, **face_normals_single, **face_normals_double, **face_transparent_double, **face_textures, **face_outlines;
  facedata **face_normals_single_DOWN_X,**face_normals_single_UP_X;
  facedata **face_normals_single_DOWN_Y,**face_normals_single_UP_Y;
  facedata **face_normals_single_DOWN_Z,**face_normals_single_UP_Z;
  int nface_normals_single_DOWN_X,nface_normals_single_UP_X;
  int nface_normals_single_DOWN_Y,nface_normals_single_UP_Y;
  int nface_normals_single_DOWN_Z,nface_normals_single_UP_Z;

  int itextureoffset;

  int mxpatch_frames;
  float vent_offset[3];
  int select_min, select_max;

  clipdata box_clipinfo;

  unsigned char *merge_color,*merge_alpha;

  char *label;

  int ncullgeominfo,nxyzgeomcull[3],nxyzskipgeomcull[3];
  struct _culldata *cullgeominfo;

#ifdef pp_CULL
  int ncullinfo;
  struct _culldata *cullinfo;
  GLuint *cullQueryId;
  int culldefined;
  struct _smoke3ddata *cull_smoke3d;
#endif

  volrenderdata volrenderinfo;

  float gslice_verts[6*3];
  int gslice_nverts,gslice_triangles[4*3],gslice_ntriangles;
  int s_offset[3];
} mesh;

/* --------------------------  supermesh ------------------------------------ */

typedef struct _supermesh {
#ifdef pp_GPU
  GLuint smoke_texture_id,fire_texture_id,blockage_texture_id;
  float *smoke_texture_buffer,*fire_texture_buffer;
#endif
  float *f_iblank_cell;
  float boxmin_scaled[3], boxmax_scaled[3];
  int drawsides[7];
  int nmeshes;
  mesh **meshes;
  int ibar, jbar, kbar;
} supermesh;

/* --------------------------  volfacelistdata ------------------------------------ */

typedef struct {
  float *xyz,dist2;
  int iwall;
  mesh *facemesh;
} volfacelistdata;

/* --------------------------  culldata ------------------------------------ */

typedef struct _culldata {
  float xbeg, xend, ybeg, yend, zbeg, zend;
  int   ibeg, iend, jbeg, jend, kbeg, kend;
  int iskip, jskip, kskip;
  int vis;
  mesh *cull_mesh;
  int npixels,npixels_old;
} culldata;

#ifdef pp_CULL
/* --------------------------  cullplanedata ------------------------------------ */

typedef struct {
  int   ibeg, iend, jbeg, jend, kbeg, kend;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  float norm[3];
  int dir;
  culldata *cull;
  mesh *cull_mesh;
} cullplanedata;

#endif


/* --------------------------  pathdata ------------------------------------ */

typedef struct _pathdata {
  float time, eye[3], xyz_view_abs[3], xyz_view_rel[3], tour_view[3];
  float az_path,zoom,elev_path;
  struct _pathdata *keysnap;
} pathdata;

/* --------------------------  keyframe ------------------------------------ */

typedef struct _keyframe {
  int selected,viewtype,npoints;
  float noncon_time, con_time, disp_time;
  pathdata nodeval;
  float keyview_xyz[3],keyview_xyz2[3];
  float total_distance, distance;
  float s_eye[3], s_az, s_elev, s_zoom, s_xyz_view[3];
  float d_eye[3], d_az, d_elev, d_zoom, d_xyz_view[3];
  float bias, continuity, tension;
  float bank;
  float az_path;
  float s1, s2, d1, d2;
  struct _keyframe *next,*prev;
} keyframe;

/* --------------------------  tourdata ------------------------------------ */

typedef struct _tourdata {
  char label[300],menulabel[128];
  keyframe first_frame,last_frame;
  keyframe **keyframe_list;
  pathdata *pathnodes;
  int glui_avatar_index, display2;
  float *path_times,*keyframe_times;
  float global_dist, local_dist;
  int *timeslist;
  int ntimes,nkeyframes;
  int display,periodic,global_tension_flag;
  int startup;
  int isDefault;
  float global_tension;
} tourdata;

/* --------------------------  tokendata ------------------------------------ */

typedef struct _tokendata {
  float var,*varptr,default_val,evac_var;
  int command,loc,type,reads,nvars,noutvars,is_label,is_string,is_texturefile;
  struct _sv_object *included_object;
  int included_frame;
  int texture_index;
  struct _tokendata *next,*elsenext;
  char token[64],tokenlabel[64],tokenfulllabel[64];
  char string[256],default_string[256],*stringptr;
} tokendata;

/* --------------------------  sv_object_frame ------------------------------------ */
#define NEVAC_TOKENS 12
typedef struct _sv_object_frame {
  int use_bw;
  int error;
  int display_list_ID;
  int *symbols, nsymbols;
  tokendata *tokens, **command_list;
  tokendata *evac_tokens[NEVAC_TOKENS];
  int ntokens,ncommands,ntextures,nevac_tokens;
  struct _sv_object *device;
  struct _sv_object_frame *prev, *next;
} sv_object_frame;

/* --------------------------  sv_object ------------------------------------ */

typedef struct _sv_object {
  char label[256];
  int type;
  int visible;
  int used, used_by_device;
  int use_displaylist;
  int select_mode;
  int nframes;
  sv_object_frame **obj_frames, first_frame, last_frame;
  struct _sv_object *prev, *next;
} sv_object;

/* --------------------------  device ------------------------------------ */

typedef struct _device{
  int active;
  int screenijk[3], visval;
  char label[30], *labelptr;
  char quantity[30], unit[30];
  float *times, *vals;
  int *valids;
  int ival,nvals,type2,type2vis;
  int in_devc_csv;
  mesh *device_mesh;
  texturedata *textureinfo;
  char *texturefile;
  int ntextures;
  float xyz[3], eyedist;
  float val;
  float xyzplot[3];
  float xyznorm[3];
  float dtheta, rotate_axis[3];
  float act_time;
  float *color, line_width;
  int filetype;
  float *act_times;
  int *state_values;
  int nparams;
  float *params;
  int istate_changes, nstate_changes, state0;
  int *showstatelist;
  int in_zone_csv;
  isosurface **plane_surface;
  propdata *prop;
  sv_object *object;
  struct _vdevicedata *vdevice;
  int type;
} devicedata;

#ifdef pp_PILOT
/* --------------------------  pilot ------------------------------------ */

typedef struct {
  float total;
  float fraction[8],vel[8];
} pilotdata;
#endif

/* --------------------------  vdevicedata ------------------------------------ */

typedef struct _vdevicedata {
  int unique;
  int filetype;
#ifdef pp_PILOT
  pilotdata pilotinfo;
#endif
  devicedata *udev,*vdev,*wdev,*valdev,*colordev,*veldev,*angledev,*sd_veldev,*sd_angledev;
} vdevicedata;

/* --------------------------  treedevice ------------------------------------ */

typedef struct {
  int nvdevices;
  vdevicedata **vdevices;
} treedevicedata;

/* --------------------------  camviewdata ------------------------------------ */

typedef struct {
  float time, eye0[3], view0[3], aperture, up[3];
  float eye[3], view[3];
} camviewdata;

/* --------------------------  camdata ------------------------------------ */

typedef struct {
  char *file, *label;
  int ncamviews;
  int available;
  float *time;
  camviewdata *camviews;
  char menulabel[128];
} camdata;

/* --------------------------  portdata ------------------------------------ */

typedef struct {
GLint left, right, down, top, width, height;
int text_height, text_width;
int doit;
} portdata;


/* --------------------------  mousedata ------------------------------------ */

typedef struct {
  int current[2], last[2], direction[2];
  float xcurrent[2], xdirection[2];
  int region;
  float angle, lastangle;
  float lasttime;
} mousedata;

/* --------------------------  camera ------------------------------------ */

typedef struct _camera {
  int defined,dirty;
  int projection_type;
  int rotation_type, rotation_index;
  float eye[3], view[3], up[3], eye_save[3];
  float isometric_y;
  float az_elev[2];
  float view_angle, azimuth, elevation;
  float xcen, ycen, zcen;
  float zoom;
  int quat_defined;
  float quaternion[4];

  int clip_mode;
  int clip_xmin, clip_xmax;
  int clip_ymin, clip_ymax;
  int clip_zmin, clip_zmax;
  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;

  int view_id;
  struct _camera *next,*prev;
  char name[301];
} camera;

/* --------------------------  part5class ------------------------------------ */

typedef struct {
  char *name;
  int kind;
  int col_diameter, col_length, col_azimuth, col_elevation;
  int col_u_vel, col_v_vel, col_w_vel;
  float dx, dy, dz;
  float diameter, length, azimuth, elevation;
  char *device_name;
  propdata *prop;
  sv_object *sphere, *smv_device;
  int vis_type;
  int maxpoints, ntypes;
  float *xyz, *rgb;
  int nvars_dep;
  int vars_dep_index[PROPVARMAX];
  float fvars_dep[PROPVARMAX];
  char *vars_dep[PROPVARMAX];
  flowlabels *labels;
} part5class;


/* --------------------------  part5prop ------------------------------------ */

typedef struct {
  flowlabels *label;
  char **partlabels, *scale;
  unsigned char *class_present, *class_vis;
  unsigned int *class_types;
  int human_property, particle_property;
  int display;
  float ppartlevels256[256];
  float valmin, valmax;
  int imin, imax;
  float global_min, global_max;
  int set_global_bounds;
  float percentile_min, percentile_max;
  float user_min, user_max;
  int setvalmin, setvalmax;
  float chopmin, chopmax;
  int setchopmin, setchopmax;
  int extreme_min, extreme_max;
  int *buckets;
} part5prop;

/* --------------------------  part5data ------------------------------------ */

typedef struct {
  part5class *partclassbase;
  float time;
  int npoints,n_rtypes, n_itypes;
  short *sx, *sy, *sz;
  float *dsx, *dsy, *dsz;
  float *avatar_angle, *avatar_width, *avatar_depth, *avatar_height;
  int humancolor_varindex;
  int *tags,*sort_tags;
  unsigned char *vis_part;
  float *rvals;
  unsigned char *irvals;
  unsigned char **cvals;
} part5data;

/* --------------------------  partdata ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file;
  char *comp_file, *size_file, *reg_file;
  int sort_tags_loaded;
  int compression_type;
  int loaded, display, reload;
  int evac;
  float zoffset;
  int blocknumber;
  int num_memblocks;
  float *times;
  int *timeslist;
  float *xpart, *ypart, *zpart, *tpart;
  short *xparts, *yparts, *zparts;
  unsigned char *xpartb, *ypartb, *zpartb;

  unsigned char *itpart,*isprink;
  int *sframe, *sprframe, *bframe, ntimes, itime;
  int particle_type, droplet_type;

  flowlabels label;
  char menulabel[128];
  int version;
  int nclasses;
  part5class **partclassptr;
  part5data *data5;
} partdata;

/* --------------------------  targ ------------------------------------ */

 typedef struct {
  int loaded,display,type;
  char *file;
} targ;

/* --------------------------  targpos ------------------------------------ */

typedef struct {
  int nsteps;
  float *x, *y, *z, *t;
  float *x2, *y2, *z2;
  float rgb[3];
  float *vals;
  float valmin,valmax;
  unsigned char *color;
} targpos;

/* --------------------------  compinfo ------------------------------------ */

typedef struct {
  int offset, size;
} compinfo;

/* --------------------------  menudata ------------------------------------ */

typedef struct {
  int menuvar;
  char label[256];
} menudata;

/* --------------------------  hrrdata ------------------------------------ */

typedef struct {
  char *file, hrrlabel[256];
  int loaded, display, *timeslist, itime;
  float *times_csv, *times, *hrrval_csv, *hrrval;
  int ntimes, ntimes_csv;
} hrrdata;

/* --------------------------  slice ------------------------------------ */

typedef struct _slicedata {
  int mesh_type;
  int seq_id, autoload;
  char *file;
  char *size_file;
  char *comp_file, *reg_file, *vol_file;
  char *slicelabel;
  int compression_type;
  int ncompressed;
  int slicetype;
  struct _multislicedata *mslice;
  int is_fed;
  feddata *fedptr;
  int menu_show;
  float *constant_color;
  float qval256[256];
  char slicedir[50];
  int loaded, display;
  int loaded_save, display_save;
  int num_memblocks;
  float position_orig;
  int blocknumber;
  int firstshort;
  int vec_comp;
  int setvalmin, setvalmax;
  float valmin, valmax;
  float globalmin, globalmax;
  float valmin_data, valmax_data;
  float diff_valmin,  diff_valmax;
  flowlabels label;
  float *qslicedata, *qsliceframe, *times, *qslice;
  unsigned char *qslicedata_compressed;
  unsigned char *slicecomplevel;
  contour *line_contours;
  int nline_contours;
  float *contour_areas;
  int *contour_areas_percen;
  int ncontour_areas;
  compinfo *compindex;
  unsigned char *slicelevel;
  char menulabel[128];
  char menulabel2[128];
  float *rgb_slice_ptr[256];
  int ntimes,itime;
  unsigned char *iqsliceframe;
  float above_ground_level;
  int volslice;
  int is1, is2, js1, js2, ks1, ks2;
  int ijk_min[3], ijk_max[3];
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float xyz_min[3], xyz_max[3];
  int nsliceii;
  int *timeslist;
  int idir;
  float sliceoffset;
  int nslicei, nslicej, nslicek;
  int nslicex, nslicey;
  int ndirxyz[4];
  int nslicetotal;
  int type;
  int vloaded;
  int reload;
  float delta_orig;
  int extreme_min, extreme_max;
} slicedata;

/* --------------------------  multislice ------------------------------------ */

typedef struct _multislicedata {
  int mesh_type;
  int seq_id, autoload;
  int loaded,display,type;
  int ndirxyz[4];
  int nslices;
  int *islices;
  float *contour_areas;
  int *contour_areas_percen;
  int ncontour_areas;
  char menulabel[128];
  char menulabel2[128];
} multislicedata;

/* --------------------------  multivslicedata ------------------------------------ */

typedef struct {
  int mesh_type;
  int seq_id, autoload;
  int loaded,display,type;
  int nvslices;
  int ndirxyz[4];
  int *ivslices;
  char menulabel[128];
  char menulabel2[128];
} multivslicedata;

/* --------------------------  databounds ------------------------------------ */

typedef struct {
  char *datalabel;
  int setvalmin, setvalmax;
  int setchopmin, setchopmax;
  float line_contour_min;
  float line_contour_max;
  int line_contour_num;
  float valmin, valmax;
  float chopmin, chopmax;
  float valmin_data,valmax_data;
  char colorlabels[12][11];
  float levels256[256];
  float fscale;
  char scale[31];
  flowlabels *label;
} databounds;

/* --------------------------  vslice ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  slicedata *u,*v,*w,*val;
  int volslice;
  int iu, iv, iw, ival;
  int loaded,display;
  float valmin, valmax;
  int type,vec_type;
  int slicetype;
  char menulabel[128];
  char menulabel2[128];
} vslicedata;

/* --------------------------  smokedata ------------------------------------ */

typedef struct {
  int ncomp_total;
  int *nchars_compressed, *nchars_compressed_full;
  unsigned char *frame_in, *frame_out, *view_tmp, *comp_all, **frame_comp_list;
} smokedata;

/* --------------------------  smoke3d ------------------------------------ */

typedef struct _smoke3ddata {
  int seq_id,autoload;
  char *file;
  char *comp_file, *reg_file;
  int filetype;
  int loaded, display, d_display;
  int is_zlib;
  int soot_loaded,water_loaded,hrrpuv_loaded;
  int blocknumber;
  int type;
  int is1, is2, js1, js2, ks1, ks2;
  int compression_type, have_light;
  flowlabels label;
  char menulabel[128];
  float *times;
  int *use_smokeframe;
  int fire_alpha;
  int *timeslist;
  int ntimes,ismoke3d_time,lastiframe,ntimes_full;
  int nchars_uncompressed;

  int ncomp_smoke_total;
  int *nchars_compressed_smoke, *nchars_compressed_smoke_full;
  unsigned char *smokeframe_in, *lightframe_in, *smokeframe_out, **smokeframe_comp_list;
  unsigned char *smokeview_tmp;
  unsigned char *smoke_comp_all;
  unsigned char *frame_all_zeros;
  smokedata smoke, light;
  unsigned char *hrrpuv_color, *water_color, *soot_color;
  int hrrpuv_index, water_index, soot_index;
  int dir;
} smoke3ddata;

/* --------------------------  patch ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file,*size_file;
  char *comp_file, *reg_file;
  char *geomfile;
  geomdata *geominfo;
  //int *patchsize;
  int version;
  int filetype;
  int type;
  int inuse,inuse_getbounds;
  int unit_start;
  int firstshort;
  int compression_type;
  int num_memblocks;
  int setvalmin, setvalmax;
  float valmin, valmax;
  int setchopmin, setchopmax;
  float chopmin, chopmax;
  float local_valmin, local_valmax;
  float diff_valmin, diff_valmax;
  int blocknumber,loaded,display;
  float *geom_times, *geom_vals;
  int *geom_timeslist,geom_itime;
  unsigned char *geom_ivals, **geom_ivals_static, **geom_ivals_dynamic;
  unsigned char *geom_ival_static, *geom_ival_dynamic;
  int geom_nval_static, geom_nval_dynamic;
  int *geom_nstatics, *geom_ndynamics;
  int geom_nvals, ngeom_times;
  flowlabels label;
  char scale[31];
  char menulabel[128];
  int extreme_min, extreme_max;
  time_t modtime;
  histogramdata *histogram;
  bounddata bounds;
} patchdata;

/* --------------------------  plot3ddata ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file,*reg_file,*comp_file;
  int compression_type;
  float time;
  int num_memblocks;
  int u, v, w, nvars;
  float diff_valmin[5], diff_valmax[5];
  int extreme_min[6], extreme_max[6];
  int blocknumber,loaded,display;
  flowlabels label[6];
  char menulabel[256],longlabel[256];
} plot3ddata;

/* --------------------------  zonedata ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file,*basefile;
  int loaded,display;
  int csv;
  flowlabels label[4];
  int setvalmin, setvalmax;
  float valmin, valmax;
  char scale[31];
} zonedata;

/* --------------------------  roomdata ------------------------------------ */

typedef struct {
  int valid;
  float dx,dy,dz;
  float x0,y0,z0;
  float x1,y1,z1;
  int drawsides[7];
  float pfloor, ylay, tl, tu, rho_L, rho_U, od_L, od_U;
  int itl, itu;
  int zoneinside;
} roomdata;

/* --------------------------  zvent ------------------------------------ */

typedef struct {
  float x1,x2,y1,y2,z1,z2,yy,zz;
  int dir;
  float dpmin, dpmax;
  float g_dpmin, g_dpmax;
  roomdata *room1, *room2;
  float area, area_fraction;
#ifdef pp_ZONEVENT
  float slab_bot[MAXSLABS], slab_top[MAXSLABS], slab_vel[MAXSLABS], slab_temp[MAXSLABS];
  int nslab;
#endif
  float *color;
  float vdata[20];
  int itempdata[20];
  int vent_orien, vent_type, face;
} zvent;

/* --------------------------  firedata ------------------------------------ */

typedef struct {
  float x, y, z, dz;
  float absx,absy,absz;
  int valid,roomnumber;
} firedata;

/* --------------------------  f_unit ------------------------------------ */

typedef struct {
  char unit[10];   /* m/s, mph etc - appears in color bar */
  float scale[2];  /* newval=scale[0]*oldval+scale[1] */
  char rel_val[20];
  int rel_defined;
} f_unit;

/* --------------------------  f_units ------------------------------------ */

typedef struct {
  int nunits;
  int unit_index,submenuid,visible;
  char unitclass[30]; /* ie: velocity, temperature */
  int diff_index;
  f_unit *units;
} f_units;

/* --------------------------  skyboxdata ------------------------------------ */

typedef struct {
  texturedata face[6];
} skyboxdata;


#endif
