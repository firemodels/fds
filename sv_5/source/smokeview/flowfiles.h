#ifndef DEF_FLOWFILES
#define DEF_FLOWFILES
#include "contourdefs.h"
#include "isodefs.h"

#ifdef pp_WUI
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
#endif
/* --------------------------  colorbardata ------------------------------------ */

typedef struct {
  char label[1024];    // menu label
  char *label_ptr;
  int npoints;      // number of legs
  int pointindex;   // leg being edited
  unsigned char *c_index;
  float *rgbnodes;
  float *c_vals;
  int *jumpflag;
  float valmin, valmax;
  float *rgb;     // color bar
} colorbardata;

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

/* --------------------------  thread_args ------------------------------------ */

typedef struct {
  char *file;
  int ifile;
  int flag;
  int errorcode;
} thread_args;

/* --------------------------  flowlabels ------------------------------------ */

typedef struct {
  char *longlabel, *shortlabel, *unit;
} flowlabels;

/* --------------------------  labeldata ------------------------------------ */

typedef struct {
	float xyz[3], rgb[4];
	char label[256];
  float tstart_stop[2];
  int useforegroundcolor;
} labeldata;


/* --------------------------  texture ------------------------------------ */

typedef struct {
  char *file;
  int loaded, display, used;
  GLuint name;
} texture;

#ifdef pp_WUI
/* --------------------------  terraindata ------------------------------------ */

typedef struct {
  int nallocated, nstates;
  float *time;
  int interval;
  unsigned char *state;
} terraincell;

typedef struct {
  char *file;
  unsigned char *state;
  int *timeslist;
  int loaded, display;
  int autoload;
  texture *ter_texture;
  int nx, ny;
  float xmin, xmax, ymin, ymax;
  float *x, *y;
  float *zcell, *znode, *znormal;
  float *times;
  terraincell *tcell;
  int ntimes;
} terraindata;
#endif

/* --------------------------  surface ------------------------------------ */

typedef struct {
  char *surfacelabel,*texturefile;
  int type; /* 
               0 - regular block non-textured 
               1 - regular block textured
               2 - outline
               3 - smoothed block
               4 - invisible
             */
  float *color, emis, temp_ignition;
  float t_width, t_height;
  texture  *textureinfo;
  int obst_surface;
  int invisible;
  int location;
  int transparent;
} surface;

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
  float *linewidth;
  float approx_center_coord[3];
  float dist2eye;
  int meshindex, blockageindex;
  int imin, imax, jmin, jmax, kmin, kmax;
  int dir,hidden;
  int del;
  int invisible;
  int transparent;
  int patchpresent;
  int **showtimelist_handle;
  int thinface;
  int show_bothsides, is_interior;
  surface *surfaceinfo;
  texture *textureinfo;
} facedata;

/* --------------------------  selectdata ------------------------------------ */

typedef struct {
  int mesh, blockage, side,dir;
  facedata *facei;
  int type;
} selectdata;

/* -------------------------- blockagedata ------------------------------------ */

typedef struct {
  int ijk[6],ijkORIG[6];
  float xmin, xmax, ymin, ymax, zmin, zmax, xyzORIG[6];
  surface *surf[6],*surfORIG[6];
  int walltype,walltypeORIG;
  int surf_index[6],surf_indexORIG[6];
  int patchvis[7];
  int usecolorindex;
  int id;
#ifdef pp_WUI
  int is_wuiblock;
#endif
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
  float texture_width, texture_height;
  float rgb[4];
  texture textureinfo;
} cadlook;

/* --------------------------  cadquad ------------------------------------ */

typedef struct {
  float xyzpoints[12];
  float txypoints[8];
  float normals[3];
  int colorindex;
  float colors[4];
  cadlook *cadlookq;
} cadquad;



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

/* --------------------------  ventdata ------------------------------------ */

typedef struct {
  int type,dummy;
  int hideboundary;
  int dir,dir2,id;
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
  float *color;
  int transparent;
  int colorindex;
  int usecolorindex;
  surface *surf[1];
  texture *textureinfo[1];
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

/* --------------------------  iso ------------------------------------ */

typedef struct {
  int seq_id, autoload;
#ifdef pp_ISO
  int compression_type;
  char *comp_file, *reg_file, *size_file;
  short *normaltable;
  unsigned char *comp_buffer, *full_bufferframe, *comp_bufferframe;
  int nfull_bufferframe, ncomp_bufferframe, maxfull_buffer;
  int nnormaltable; 
#endif
  char *file;
  int dataflag;
  int type;
  int num_memblocks;
  flowlabels label;
  int blocknumber,display,loaded;
  float tmin,tmax;
  char menulabel[128];

} iso;

/* --------------------------  smoothblockage ------------------------------------ */

typedef struct {
  int nsmoothblockagecolors;
  float *smoothblockagecolors;
  isosurface **smoothblockagesurfaces;
  float time;
} smoothblockage;

/* --------------------------  mesh ------------------------------------ */

typedef struct mesh_ {
#ifdef pp_PART5
  int mesh_type;
#endif
  float meshrgb[3], *meshrgb_ptr;
  float cellsize;
  float *xplt, *yplt, *zplt;
  float *xplt_orig, *yplt_orig, *zplt_orig;
  float xbar0, xbar, ybar0, ybar, zbar0, zbar;
  float xcen, ycen, zcen;
  float offset[3];
  float xyzmaxdiff;
  float boxoffset;
  float hrrpuv_cutoff;
  int plot3dfilenum,isofilenum,patchfilenum;
  int ibar, jbar, kbar;
  int plotx, ploty, plotz;
  int visx, visy, visz;
  int visx2, visy2, visz2;
  int slicedir;
  int plotn;
  int *iblank,*iblank_cell,*iblank_x,*iblank_y,*iblank_z;
  float *block_zdist;
  int zdist_flag;
  unsigned char *iblank_smoke3d;
  blockagedata **blockageinfoptrs,**deletelist,**carveblockptrs;
  int *obst_bysize;
  int ndeletelist,ncarveblocks;
  ventdata *ventinfo;
  int nvents,ndummyvents;
  int nbptrs;
  int *iqdata;
  float *qdata;
  int *yzcolorbase, *xzcolorbase, *xycolorbase; 
  float *yzcolorfbase, *xzcolorfbase, *xycolorfbase;
#ifdef pp_PLOTTEXTURE
  float *yzcolortbase, *xzcolortbase, *xycolortbase;
#endif
  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  int *iblank_xy, *iblank_xz, *iblank_yz;
  float plot3d_speedmax;
  contour plot3dcontour1,plot3dcontour2,plot3dcontour3;
  isosurface currentsurf,currentsurf2;
  isosurface *blockagesurface;
  isosurface **blockagesurfaces;
  float *smoothblockagecolors;
  int nsmoothblockagecolors;
  int ntc;
#ifndef pp_DEVICE
  float *xtcplot, *ytcplot, *ztcplot;
  float *xtc, *ytc, *ztc;
  float *xtcnorm, *ytcnorm, *ztcnorm;
  int *has_tcnorm;
#endif
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
  int nisosteps;
  float *isotimes;
  int *isotimeslist;
  int iiso;
  int smokedir;
  float dx, dy, dz, dxy,dxz,dyz;
  float norm[3];

  int *patchtype;
//  int *patchfacevis;
  int *patchdir,*patch_surfindex;
  int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2;
  contour **patch_contours;
  int *blockonpatch;
  struct mesh_ **meshonpatch;
  int *ptype;
  int *patchrow, *patchcol, *blockstart;
  unsigned int *zipoffset, *zipsize;
  int *visPatches;
  float *xyzpatch, *xyzpatch_ignited;
  unsigned char *ipqq_zlib, *ipqqi_zlib;
  unsigned char *ipqq, *ipqqi;
  float *patchtimes, *patchtimesi, *pqq, *pqqi;
  float *chartime;
  int *patchblank;
  int npatch_frames,npatches;
  int ipatch;
  int *patchtimeslist;
  int npatchsize;
  int patchfacevis2;
  int visInteriorPatches;
  float surface_tempmin, surface_tempmax;

  smoothblockage *smoothblockages_list;
  int nsmoothblockages_list;

  smoothblockage **showsmoothtimelist;

  int nface_textures, nface_outlines, nfaces;
  int nface_normals_single, nface_normals_double;
  facedata *faceinfo, **face_normals_single, **face_normals_double, **face_textures, **face_outlines;


  int itextureoffset;

  int mxpatch_frames;
  float vent_offset[3];
  int select_min, select_max;

  int *outline_blank;
  float *outline_list;
  unsigned char *merge_color,*merge_alpha;

  char *label;
  int smokeloaded;

} mesh;



/* --------------------------  pathdata ------------------------------------ */

typedef struct _pathdata {
  float time, eye[4], aview[3], oview[3];
  float zoom,elevation;
  struct _pathdata *keysnap;
} pathdata;

/* --------------------------  keyframe ------------------------------------ */

typedef struct _keyframe {
  int selected,viewtype,npoints;
  float noncon_time, con_time, disp_time;
  pathdata nodeval;
  float keyview_x, keyview_y;
  float cumdist, dist;
  float s_eye[6], d_eye[6];
  float s_aview[3], d_aview[3];
  float bias, continuity, tension;
  float bank;
  float azimuth;
  float s1, s2, d1, d2;
  struct _keyframe *next,*prev;
} keyframe;

/* --------------------------  tourdata ------------------------------------ */

typedef struct _tourdata {
  char label[300],menulabel[128];
  keyframe first_frame,last_frame;
  keyframe **keyframe_list;
  pathdata *pathnodes;
  float *path_times,*keyframe_times;
  float global_dist, local_dist;
  int *path_timeslist;
  int npath,nkeyframes;
  int display,periodic,global_tension_flag;
  int startup;
  int isDefault;
  float global_tension;
} tourdata;

#ifdef pp_DEVICE

/* --------------------------  sv_object_parts ------------------------------------ */

  // part/op types
  // 0 translatexyz/3, 1 rotatex/1, 2 rotatey/1, 3 rotatez/1, 4 scalexyz/3, 5 scale/1 6 drawcube/1, 7 drawdisk/2, 8 drawsphere/1, 9 drawline  10 drawelliposid/3

//typedef struct {
//  float *rgb4;
//  int nops,nargs;
//  float *args;
//  int *ops;
//} sv_object_part;

/* --------------------------  sv_object_frame ------------------------------------ */

typedef struct _sv_object_frame {
//  int nparts;
//  sv_object_part *parts;
  int *ops, nops, nargs;
  int error;
  int display_list_ID;
  float *args;
  struct _sv_object_frame *prev, *next;
} sv_object_frame;

/* --------------------------  sv_object ------------------------------------ */

typedef struct _sv_object {
  char label[256];
  int visible;
  int used;
  int nframes;
  sv_object_frame **obj_frames, first_frame, last_frame;
 struct _sv_object *prev, *next;
} sv_object;

/* --------------------------  device ------------------------------------ */

typedef struct {
  int active;
  float xyz[3];
  float val;
  float xyzplot[3];
  float xyznorm[3];
  float angle_elev, angle_az;
  float act_time;
  sv_object *object;
  int type;
} device;
#endif

/* --------------------------  casedata ------------------------------------ */

typedef struct {
  mesh *meshinfo;
  char label[128];
  int nobst,nmeshes,meshoffset;
} casedata;

/* --------------------------  camera ------------------------------------ */

typedef struct {
  float time, eye0[3], view0[3], aperture, up[3];
  float eye[3], view[3];
} camviewdata;


/* --------------------------  camera ------------------------------------ */

typedef struct {
  char *file, *label;
  int ncamviews;
  int available;
  float *time;
  camviewdata *camviews;
  char menulabel[128];
} camdata;

/* --------------------------  camera ------------------------------------ */

typedef struct _camera {
  int defined,dirty;
  int eyeview, rotation_index;
  float eye[3], view[3], up[3];
  float isometric_y;
  float angle_zx[2];
  float view_angle, direction_angle, elevation_angle;
  float sin_view_angle, cos_view_angle;
  float sin_direction_angle, cos_direction_angle;
  float sin_elevation_angle, cos_elevation_angle;
  float xcen, ycen, zcen;
  float zoom;
  float modelview[16];
  int view_id;
  struct _camera *next,*prev;
  char name[32];
} camera;

#ifdef pp_PART5
/* --------------------------  partclass ------------------------------------ */

typedef struct {
  char *name;
  int kind;
  int maxpoints, ntypes;
  float *xyz, *rgb;
  flowlabels *labels;
} part5class;


/* --------------------------  partvarprop ------------------------------------ */

typedef struct {
  flowlabels *label;
  char **partlabels, *scale;
  unsigned char *class_present, *class_vis;
  unsigned int *class_types;
  int display;
  float ppartlevels256[256];
  float valmin, valmax;
  float global_min, global_max;
  float percentile_min, percentile_max;
  float user_min, user_max;
  int setvalmin, setvalmax;
  float chopmin, chopmax;
  int setchopmin, setchopmax;
  int *buckets;
} part5prop;

/* --------------------------  part5data ------------------------------------ */

typedef struct {
  part5class *partclassbase;
  float time;
  int npoints,n_rtypes, n_itypes;
  short *sx, *sy, *sz;
  int *tags,*sort_tags;
  unsigned char *vis_part;
  float *rvals;
  unsigned char *irvals;
  int **ivals;
  unsigned char **cvals;
} part5data;
#endif

/* --------------------------  particle ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file;
#ifdef pp_PART5
  char *comp_file, *size_file, *reg_file;
  int sort_tags_loaded;
  int compression_type;
#endif
  int loaded, display;
  int evac;
  int blocknumber;
  int num_memblocks;
  float *ptimes;
  int *ptimeslist;
  float *xpart, *ypart, *zpart, *tpart;
  short *xparts, *yparts, *zparts;
  unsigned char *xpartb, *ypartb, *zpartb;

  unsigned char *itpart,*isprink;
  int *sframe, *sprframe, *bframe, nframes, iframe;
  int particle_type, droplet_type;

  flowlabels label;
  char menulabel[128];
  int version;
#ifdef pp_PART5
  int nclasses;
  part5class **partclassptr;
  part5data *data5;
#endif

} particle;

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

typedef struct {
  int offset, size;
} compinfo;


/* --------------------------  menudata ------------------------------------ */

typedef struct {
  int menuvar;
  char label[256];
} menudata;

/* --------------------------  slice ------------------------------------ */

typedef struct {
#ifdef pp_PART5
  int mesh_type;
#endif
  int seq_id, autoload;
  char *file;
#ifdef pp_SLICE
  char *size_file;
  char *rle_file, *comp_file, *reg_file;
  int compression_type;
  int ncompressed;
  float qval256[256];
#endif
  char slicedir[50];
  int loaded, display, benchvis;
  int num_memblocks;
  float position;
  int *iblank;
  int blocknumber;
  int firstshort;
  int vec_comp;
  int setvalmin, setvalmax;
  float valmin, valmax;
  float globalmin, globalmax;
  float valmin_data, valmax_data;
  flowlabels label;
  float *qslicedata, *slicetimes, *qslice;
#ifdef pp_SLICE
  unsigned char *qslicedata_compressed;
  unsigned char *slicecomplevel;
  compinfo *compindex;
  unsigned char *slicelevel;
#else
  int *slicelevel;
#endif
  char menulabel[128];
  char menulabel2[128];
  int nsteps,islice;
  float *slicedata;
#ifdef pp_SLICE
  unsigned char *slicepoint;
  int volslice;
#else
  int *slicepoint;
#endif
  int is1, is2, js1, js2, ks1, ks2;
  int nsliceii;
  int *slicetimeslist;
  int idir;
  float sliceoffset;
  int nslicei, nslicej, nslicek;
  int nslicex, nslicey;
  int nslicetotal;
  int type;
  int vloaded;
  int reload;
  float delta;

} slice;

/* --------------------------  multislice ------------------------------------ */

typedef struct {
#ifdef pp_PART5
  int mesh_type;
#endif
  int seq_id, autoload;
  int loaded,display,type;
  int nslices;
  int *islices;
  char menulabel[128];
  char menulabel2[128];
} multislice;

/* --------------------------  multivslice ------------------------------------ */

typedef struct {
#ifdef pp_PART5
  int mesh_type;
#endif
  int seq_id, autoload;
  int loaded,display,type;
  int nvslices;
  int *ivslices;
  char menulabel[128];
  char menulabel2[128];
} multivslice;

/* --------------------------  databounds ------------------------------------ */

typedef struct {
  char *datalabel;
  int setvalmin, setvalmax;
  int setchopmin, setchopmax;
  float valmin, valmax;
  float chopmin, chopmax;
  float valmin_data,valmax_data;
  char colorlabels[12][11];
  float levels256[256];
  char scale[31];
  flowlabels *label;
} databounds;

/* --------------------------  vslice ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  slice *u,*v,*w,*val;
#ifdef pp_SLICE
  int volslice;
#endif
  int iu, iv, iw, ival;
  int loaded,display;
  float valmin, valmax;
  int type,vec_type;
  char menulabel[128];
  char menulabel2[128];
} vslice;

typedef struct {
  int ncomp_total;
  int *nchars_compressed, *nchars_compressed_full;
  unsigned char *frame_in, *frame_out, *view_tmp, *comp_all, **frame_comp_list;
} smokedata;

/* --------------------------  smoke3d ------------------------------------ */

typedef struct {
  int seq_id,autoload;
  char *file;
  char *comp_file, *reg_file;
#ifdef pp_LIGHT
  int use_lighting_file;
  char *light_file;
#endif
  int case_number;
  int loaded, display, d_display;
  int soot_loaded,water_loaded,hrrpuv_loaded;
  int blocknumber;
  int type;
  int is1, is2, js1, js2, ks1, ks2;
  int version;
  flowlabels label;
  char menulabel[128];
  float *times;
  int fire_alpha;
  int *timeslist;
  int n_times,iframe,lastiframe,n_times_full;
  int nchars_uncompressed;

  int ncomp_smoke_total;
  int *nchars_compressed_smoke, *nchars_compressed_smoke_full;
  unsigned char *smokeframe_in, *smokeframe_out, **smokeframe_comp_list;
  unsigned char *smokeview_tmp;
  unsigned char *smoke_comp_all;
  smokedata smoke, light;
#ifdef pp_LIGHT
  int ncomp_light_total;
  int *nchars_compressed_light, *nchars_compressed_light_full;
  unsigned char *lightframe_in, *lightframe_out, **lightframe_comp_list;
  unsigned char *lightview_tmp;
  unsigned char *light_comp_all;
#endif
  unsigned char *hrrpuv_color, *water_color, *soot_color;
  int hrrpuv_index, water_index, soot_index;
  int dir;

} smoke3d;

/* --------------------------  patch ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file,*size_file;
  char *comp_file, *reg_file;
  int version;
  int type;
  int firstshort;
  int compression_type;
  int num_memblocks;
  int setvalmin, setvalmax;
  float valmin, valmax;
  int blocknumber,loaded,display;
  flowlabels label;
  char scale[31];
  char menulabel[128];

} patch;

/* --------------------------  plot3d ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file;
  float time;
  int num_memblocks;
  int u, v, w, nvars;
  int blocknumber,loaded,display;
  flowlabels label[6];
  char menulabel[256],longlabel[256];
} plot3d;

/* --------------------------  zone ------------------------------------ */

typedef struct {
  int seq_id, autoload;
  char *file,*basefile;
  int loaded,display;
  flowlabels label[4];
  int setvalmin, setvalmax;
  float valmin, valmax;
  char scale[31];
} zone;

/* --------------------------  roomdata ------------------------------------ */

typedef struct {
  int valid;
  float dx,dy,dz;
  float x0,y0,z0;
  float x1,y1,z1;
  float pfloor, ylay, tl, tu, rho_L, rho_U;
  int itl, itu;
} roomdata;

/* --------------------------  zvent ------------------------------------ */

typedef struct {
  float x1,x2,y1,y2,z1,z2,yy,zz;
  int dir;
  roomdata *room1, *room2;
  float area;
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
} f_unit;

/* --------------------------  f_units ------------------------------------ */

typedef struct {
  int nunits;
  int active,submenuid;
  char unitclass[30]; /* ie: velocity, temperature */
  char unittypes[10][256]; /* list of equivalent units ie: U-VEL, V-VEL, W-VEL */
  int nunittypes;
  f_unit *units;
} f_units;

/*
smokeview.ini format

UNIT
nunit
menulabel
colorbarlabel
scale[0] scale[1]
.
.
.
menulabel
colorbarlabel
scale[0] scale[1]   (nunit times)
*/

typedef struct {
  texture face[6];
} skyboxdata;


#endif
