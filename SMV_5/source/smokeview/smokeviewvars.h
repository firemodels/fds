// $Date$ 
// $Revision$
// $Author$

#include <time.h>

#include "csphere.h"

#ifndef CPP
#include "smokeviewdefs.h"
#include "smokeheaders.h"
#endif

#include "threader.h"
#include "string_util.h"

SVEXTERN int SVDECL(usemenu,1);
SVEXTERN float hrrpuv_iso_color[4];
SVEXTERN int show_slice_terrain;
SVEXTERN int npropinfo;
SVEXTERN propdata *propinfo;
SVEXTERN float right_green, right_blue;

SVEXTERN int levelset_colorbar, wallthickness_colorbar;
SVEXTERN float glui_time;
SVEXTERN int show_mode;
SVEXTERN int cellcenter_interp, cellcenter_slice_active, cellcenter_bound_active;
SVEXTERN int part5colorindex, show_tracers_always;
SVEXTERN int navatar_colors;
SVEXTERN int select_avatar, selected_avatar_tag, view_from_selected_avatar;
SVEXTERN int select_device, selected_device_tag;
SVEXTERN float selected_avatar_pos[3], selected_avatar_angle;
SVEXTERN unsigned char select_device_color[4], *select_device_color_ptr;
SVEXTERN float *avatar_colors;
SVEXTERN int script_render_flag, script_itime;

SVEXTERN int show_slice_in_obst, offset_slice;
SVEXTERN int skip_slice_in_embedded_mesh;
SVEXTERN int n_embedded_meshes;

SVEXTERN int force_isometric;
SVEXTERN int updategluiview;
SVEXTERN int render_double,render_double_state,render_double_menu,render_from_menu;
SVEXTERN int usetexturebar;
SVEXTERN int show_smokelighting;
SVEXTERN int sb_atstart;
#ifdef pp_CULL
SVEXTERN int cullactive, show_cullports, cull_portsize;
SVEXTERN int cullsmoke, ncullplaneinfo;
SVEXTERN cullplanedata *cullplaneinfo;
SVEXTERN cullplanedata **sort_cullplaneinfo;
SVEXTERN int have_setpixelcount,update_initcullplane;
#endif
#ifdef pp_GPU
SVEXTERN int usegpu,gpuactive;
SVEXTERN int GPU_aspectratio;
SVEXTERN int GPU_smoke3d_rthick, GPU_skip, GPU_hrrcutoff, GPU_hrr;
SVEXTERN int GPU_firecolor, GPU_is_smoke;
SVEXTERN int GPU_smokeshade,GPU_smokealpha;
SVEXTERN int GPU_blank;
SVEXTERN int GPU_adjustalphaflag;
SVEXTERN unsigned int GPU_depthtexture;
SVEXTERN int i_hrrcutoff;
#endif

SVEXTERN float partfacedir[3];
SVEXTERN int demo_option;
SVEXTERN int small_font_height, large_font_height;
SVEXTERN float cb_valmin, cb_valmax, cb_val;
SVEXTERN int cb_colorindex;
SVEXTERN float rgbterrain[4*MAXRGB];
SVEXTERN int terrain_rgba_zmin[3];
SVEXTERN int terrain_rgba_zmax[3];
SVEXTERN float vertical_factor;

SVEXTERN char inputfilename_ext[4];

SVEXTERN float percentile_level;
SVEXTERN float fire_line_min, fire_line_max;
SVEXTERN int update_fire_line;
SVEXTERN int fire_line_index;
SVEXTERN int slice_bounds_dialog;

SVEXTERN int dwinHbase;
SVEXTERN int dwinH;

SVEXTERN float xtemp;


SVEXTERN char TITLEBASE[1024];

SVEXTERN float set_view_xyz[3];
SVEXTERN char INIfile[1024];
SVEXTERN char WRITEINIfile[1024];

SVEXTERN spherepoints SVDECL(*sphereinfo,NULL), SVDECL(*wui_sphereinfo,NULL);

SVEXTERN float tourcol_selectedpathline[3];
SVEXTERN float tourcol_selectedpathlineknots[3];
SVEXTERN float tourcol_selectedknot[3];
SVEXTERN float tourcol_selectedview[3];

SVEXTERN float tourcol_pathline[3];
SVEXTERN float tourcol_pathknots[3];
SVEXTERN float tourcol_text[3];

SVEXTERN float tourcol_avatar[3];



SVEXTERN float mat_ambient_orig[4];
SVEXTERN float mat_specular_orig[4];
SVEXTERN float *mat_ambient2;
SVEXTERN float *mat_specular2;

SVEXTERN GLfloat iso_specular[4];
SVEXTERN GLfloat iso_shininess;

SVEXTERN float block_ambient_orig[4];
SVEXTERN float *block_ambient2;

SVEXTERN GLfloat block_shininess;

SVEXTERN GLfloat light_position0[4];
SVEXTERN GLfloat light_position1[4];

SVEXTERN GLfloat ambientlight[4];
SVEXTERN GLfloat diffuselight[4];

SVEXTERN GLint screenWidth2, screenHeight2;

SVEXTERN int list_p3_index,list_slice_index,list_patch_index,list_iso_index;
SVEXTERN int list_p3_index_old, list_slice_index_old, list_patch_index_old,list_iso_index_old;

SVEXTERN float glui_block_xmin, glui_block_ymin, glui_block_zmin;
SVEXTERN float glui_block_xmax, glui_block_ymax, glui_block_zmax;

SVEXTERN float zonelevels256[256];
SVEXTERN float boundarylevels256[256];
SVEXTERN float partlevels256[256];
SVEXTERN float *zonet, *zoneylay, *zonetl, *zonetu, *zonepr, *zoneqfire;
SVEXTERN unsigned char *hazardcolor;
SVEXTERN unsigned char *izonetu;
SVEXTERN int nzonet;
SVEXTERN float barright;
SVEXTERN float *tspr;
SVEXTERN float tmin, tmax;

SVEXTERN int videoSTEREO;
SVEXTERN float fzero;

SVEXTERN char blank[2];

SVEXTERN float *sphere_xyz;
SVEXTERN int demo_mode;
SVEXTERN int update_demo;
SVEXTERN int menu_view_number;
SVEXTERN int mxplot3dvars;
SVEXTERN int loadplot3dall;
SVEXTERN char *shortp3label[MAXPLOT3DVARS], *unitp3label[MAXPLOT3DVARS];
SVEXTERN char *LESsystem,*LESendian;

SVEXTERN int show3dsmoke;
SVEXTERN float frustum[6][4];
SVEXTERN int showtime, showtime2, showplot3d, showpatch, showslice, showvslice, showsmoke, showzone, showiso, showevac;
SVEXTERN int vis_slice_contours;
SVEXTERN int update_slicecontours;
SVEXTERN int showevac_colorbar;
SVEXTERN int showiso_colorbar;
SVEXTERN int visgridloc;

SVEXTERN int valindex;

SVEXTERN float *rgb2_ini;
SVEXTERN float rgb_full[MAXRGB][4];
SVEXTERN float rgb_full2[MAXRGB][4];
SVEXTERN float rgb_slice[4*MAXRGB];
SVEXTERN float rgb_iso[4*MAXRGB];
SVEXTERN float rgb_patch[4*MAXRGB];
SVEXTERN float rgb_plot3d[4*MAXRGB];
SVEXTERN float rgb_part[4*MAXRGB];
SVEXTERN float rgb_trans[4*MAXRGB];
SVEXTERN float rgb_cad[MAXRGB][4];

SVEXTERN float iso_ambient[12];
SVEXTERN float iso_transparency;
SVEXTERN int n_iso_ambient;
SVEXTERN float *iso_ambient_ini;
SVEXTERN int n_iso_ambient_ini;
SVEXTERN float *rgb_ini;
SVEXTERN float rgb[MAXRGB][4];
SVEXTERN float mouse_deltax, mouse_deltay;
SVEXTERN float **rgbptr, **rgb_plot3d_contour;
SVEXTERN float char_color[4];
SVEXTERN float *rgb_step[255];
SVEXTERN float movedir[3];
SVEXTERN float rgb_base[MAXRGB][4];
SVEXTERN float bw_base[MAXRGB][4];
SVEXTERN int nrgb2;
SVEXTERN float rgb2[MAXRGB][3];
SVEXTERN float rgbhazard[MAXRGB][4];
SVEXTERN float inverse_modelview_setup[16];
SVEXTERN float modelview_setup[16];
SVEXTERN float modelview_rotate_last[16],modelview_rotate_save[16];
SVEXTERN float modelview_current[16];
SVEXTERN float modelview_scratch[16];

SVEXTERN camera *camera_current, *camera_save, *camera_last;
SVEXTERN camera *camera_external, *camera_internal, *camera_ini, *camera_external_save;
SVEXTERN camera camera_list_first, camera_list_last, **camera_list;
SVEXTERN int ncamera_list,i_view_list,init_camera_list_flag;
SVEXTERN int camera_max_id;
SVEXTERN int startup,startup_view_ini,selected_view;
SVEXTERN char label_startup_view[256];
SVEXTERN char *camera_label, *colorbar_label;

SVEXTERN int visPatchType[7];
SVEXTERN int setp3min[MAXPLOT3DVARS],p3_extreme_min[MAXPLOT3DVARS],p3_extreme_max[MAXPLOT3DVARS];
SVEXTERN float p3min[MAXPLOT3DVARS];
SVEXTERN float p3chopmin[MAXPLOT3DVARS];
SVEXTERN int setp3max[MAXPLOT3DVARS];
SVEXTERN int setp3chopmin[6];
SVEXTERN int setp3chopmax[6];


SVEXTERN float p3max[MAXPLOT3DVARS];
SVEXTERN float p3chopmax[MAXPLOT3DVARS];

SVEXTERN int trainer_pause;
SVEXTERN int trainee_location;
SVEXTERN int trainer_inside;
SVEXTERN int from_glui_trainer;
SVEXTERN int trainer_path_old;
SVEXTERN int trainer_outline;
SVEXTERN int trainer_viewpoints,ntrainer_viewpoints;
SVEXTERN int trainer_realtime;
SVEXTERN int trainer_path;
SVEXTERN float trainer_xzy[3],trainer_ab[2];
SVEXTERN float motion_ab[2], motion_dir[2];
SVEXTERN int trainerload;
SVEXTERN int fontsize_save;
SVEXTERN int showtrainer;
SVEXTERN int trainer_mode;
SVEXTERN int trainer_active;
SVEXTERN int slice_average_flag,slice_turbprop_flag;
SVEXTERN int show_slice_average,vis_slice_average;
SVEXTERN float slice_average_interval;

SVEXTERN float angle, dang, tourangle;
SVEXTERN int maxtourframes;
SVEXTERN int blockageSelect;
SVEXTERN unsigned int ntourknots;
SVEXTERN int itourknots;
SVEXTERN int stretch_var_black, stretch_var_white, move_var;

SVEXTERN int showhide_option;
SVEXTERN int snifferrornumber;
SVEXTERN int xyz_dir;
SVEXTERN int which_face;
SVEXTERN int showfontmenu;
SVEXTERN int showlightmenu;

SVEXTERN float VECFRACTION,vecfactor,veclength;
SVEXTERN int iveclengths;

SVEXTERN int glui_active;

SVEXTERN int drawColorLabel,olddrawColorLabel;
SVEXTERN int staticframe0,visStaticSmoke;
SVEXTERN int vis3DSmoke3D;
SVEXTERN int smokeskip,smokeskipm1;
SVEXTERN int nrooms,nzone, nzvents, nfires;
SVEXTERN float ratio, aspect;
SVEXTERN int visLIGHT0, visLIGHT1, visLIGHTMENU, UpdateLIGHTS;

SVEXTERN int screenWidth, screenHeight;
SVEXTERN int renderW, renderH, render_option;
SVEXTERN int glui_screenWidth, glui_screenHeight;
SVEXTERN int windowsize_pointer;
SVEXTERN int sethazardcolor;
SVEXTERN int mxpoints,mxframes,mxframepoints;
SVEXTERN int mxpoints_orig,mxframes_orig;
SVEXTERN int mxpoints_comm, mxframes_comm;
SVEXTERN int timedrag,colordrag,colorsplitdrag;
SVEXTERN int isonormtype,showisonormals;
SVEXTERN int global_changecolorindex;
SVEXTERN int fontindex,fontWoffset,fontHoffset;

SVEXTERN float xcenGLOBAL, ycenGLOBAL, zcenGLOBAL;
SVEXTERN float xbar, ybar, zbar;
SVEXTERN float xbar0, ybar0, zbar0;
SVEXTERN float xbarORIG, ybarORIG, zbarORIG;
SVEXTERN float xbar0ORIG, ybar0ORIG, zbar0ORIG;
SVEXTERN int ReadPlot3dFile, ReadIsoFile;
SVEXTERN int ReadVolSlice;
SVEXTERN int Read3DSmoke3DFile;
SVEXTERN int ReadZoneFile, ReadPartFile, ReadEvacFile;

SVEXTERN int unload_qdata;

SVEXTERN int editwindow_status;
SVEXTERN int startup_pass;
SVEXTERN int ntargtimes;
SVEXTERN int showtitle1, showtitle2;

SVEXTERN int slicefilenumber;
SVEXTERN int exportdata;
SVEXTERN int count, lastcount;
SVEXTERN int nspr;
SVEXTERN int RenderGif, RenderSkip;
SVEXTERN int isoframestep;
SVEXTERN int isoframeskip;
SVEXTERN int smoke3dframestep;
SVEXTERN int smoke3dframeskip;
SVEXTERN int vectorskip;
SVEXTERN int iframe, iframebeg, izone;
SVEXTERN int eyeview,eyeview_level;
SVEXTERN int eyeview_old,eyeview_SAVE,eyeview_last;
SVEXTERN int frameratevalue;
SVEXTERN int setpartmin, setpartmax, setslicemin, setslicemax, endian;
SVEXTERN float slice_line_contour_min;
SVEXTERN float slice_line_contour_max;
SVEXTERN int slice_line_contour_num;
SVEXTERN int setpartmin_old, setpartmax_old;
SVEXTERN int setpatchmin, setpatchmax, setzonemin, setzonemax;
SVEXTERN int loadpatchbysteps;
SVEXTERN int settargetmin, settargetmax;
SVEXTERN int setpartchopmin, setpartchopmax;
SVEXTERN int SVDECL(setslicechopmin,0), SVDECL(setslicechopmax,0);
SVEXTERN int SVDECL(setpatchchopmin,0), SVDECL(setpatchchopmax,0);
SVEXTERN float partchopmin,  partchopmax;
SVEXTERN float slicechopmin, slicechopmax;
SVEXTERN float SVDECL(patchchopmin,1.0), SVDECL(patchchopmax,0.0);
SVEXTERN int setisomin, setisomax;
SVEXTERN float isomin, isomax;
SVEXTERN int setisochopmin, setisochopmax;
SVEXTERN float isochopmin, isochopmax;

SVEXTERN int vis_onlythreshold, vis_threshold, canshow_threshold, activate_threshold;
SVEXTERN int settmin_p, settmin_b, settmin_s, settmin_z, settmin_i;
SVEXTERN int settmax_p, settmax_b, settmax_s, settmax_z, settmax_i;
SVEXTERN int set_no_part;
SVEXTERN float tmin_p, tmin_b, tmin_s, tmin_z, tmin_i;
SVEXTERN float tmax_p, tmax_b, tmax_s, tmax_z, tmax_i;
SVEXTERN float patchmin, patchmax;
SVEXTERN float targetmin, targetmax;
SVEXTERN float partmin, partmax, slicemin, slicemax;
SVEXTERN float zonemin, zonemax;
SVEXTERN float speedmax;
SVEXTERN int axissmooth;
SVEXTERN propdata *prop_evacdefault;
SVEXTERN float hrrpuv_max_smv;
SVEXTERN int FlowDir,ClipDir;
SVEXTERN int plotn;
SVEXTERN int stept;
SVEXTERN int plotstate;
SVEXTERN int visVector;
SVEXTERN int visSmokePart, visSprinkPart, havesprinkpart;
SVEXTERN int visaxislabels;
SVEXTERN int numplot3dvars;
SVEXTERN int skip;
SVEXTERN int p3dsurfacesmooth;
SVEXTERN int p3dsurfacetype;
SVEXTERN int parttype;
SVEXTERN int allexterior,showexterior;
SVEXTERN int allinterior;
SVEXTERN int showbounds,showmotion,showedit, showclip, showgluistereo, showtour, showlabels, showcolorbar, showwui;
#ifdef pp_SHOOTER
SVEXTERN int showshooterDLG;
SVEXTERN float shooter_xyz[3], shooter_dxyz[3], shooter_uvw[3], shooterpointsize;
SVEXTERN float shooter_velx, shooter_vely, shooter_velz, shooter_time, shooter_time_max;
SVEXTERN int shooter_cont_update,shooter_firstframe;
SVEXTERN float shooter_u0, shooter_z0, shooter_p, shooter_v_inf;
SVEXTERN float shooter_velmag, shooter_veldir, shooter_duration, shooter_history;
SVEXTERN int shooter_active;
SVEXTERN int shooter_fps,shooter_vel_type, shooter_nparts, visShooter, showshooter, nshooter_frames, max_shooter_points;
SVEXTERN shootpointdata *shootpointinfo;
SVEXTERN shoottimedata *shoottimeinfo;
SVEXTERN int *shooter_timeslist;
SVEXTERN int shooter_itime;
#endif
SVEXTERN int showgluitrainer;
SVEXTERN int colorbarcycle;
SVEXTERN int colorbartype,colorbartype_ini,colorbartype_default;
SVEXTERN int colorbartype_save;
SVEXTERN int colorbarpoint;
SVEXTERN int vectorspresent;
SVEXTERN int cb_hidesv;

SVEXTERN int visTarg, ReadTargFile;
SVEXTERN int showtarget;
SVEXTERN int visAIso;
SVEXTERN int surfincrement,visiso;
SVEXTERN int  isotest;
SVEXTERN int isolevelindex, isolevelindex2;
SVEXTERN float pref,pamb,tamb;
SVEXTERN int ntc_total, nspr_total, nheat_total;
SVEXTERN int n_devices;

SVEXTERN int npart_files, nslice_files, nvslice, nslice2, npatch2, nplot3d_files, npatch_files;
SVEXTERN int nevac;
SVEXTERN int current_particle_type,last_particle_type;
SVEXTERN int nsmoke3d_files;
SVEXTERN int niso_files, niso_bounds;
SVEXTERN int ntrnx, ntrny, ntrnz,npdim,nmeshes,clip_mesh;
SVEXTERN int nobst,nvent,noffset;
SVEXTERN int nlabels,visLabels,nlabelssmv;
SVEXTERN int ntarg_files;
SVEXTERN int showallslicevectors;
SVEXTERN float framerate;
SVEXTERN int ntimes, itime, itime_save, itimeold, seqnum,RenderTime,RenderTimeOld;
SVEXTERN int npqq, nopart;
SVEXTERN int uindex, vindex, windex;

SVEXTERN int p3cont2d, p3cont3dsmooth;
SVEXTERN int cullfaces;
SVEXTERN int showonly_hiddenfaces;
SVEXTERN int blockage_index;


SVEXTERN int windowresized;

SVEXTERN int SVDECL(updatemenu,0), first_display;
SVEXTERN int updatezoommenu;
SVEXTERN int updatemenu_count;
SVEXTERN int no_graphics;

SVEXTERN int updatefaces,updatefacelists;
SVEXTERN int updateOpenSMVFile;

SVEXTERN int periodic_reloads;
SVEXTERN int periodic_value;

SVEXTERN int slicefilenum;
SVEXTERN int partfilenum,zonefilenum;
SVEXTERN int targfilenum;

SVEXTERN int update_makeiblank_smoke3d, update_initcull;
SVEXTERN int setPDIM;
SVEXTERN int menustatus;
SVEXTERN int visTimeZone, visTimeSmoke, visTimeSlice, visTimePatch, visTimeIso, visTimeEvac;
SVEXTERN int vishmsTimelabel, visTimeLabels, visColorLabels;
SVEXTERN int visTitle, visFullTitle, visFramerate, visFramelabel, visTimelabel;
SVEXTERN int visHRRlabel;
#ifdef pp_memstatus
SVEXTERN int visAvailmemory;
#endif
SVEXTERN slice *sd_shown;
SVEXTERN vslice *vd_shown;
SVEXTERN int show_all_slices;
SVEXTERN int autoterrain,manual_terrain;
SVEXTERN float zterrain_max, zterrain_min;
SVEXTERN int revision_smv, revision_fds;
SVEXTERN int visBlocklabel;
SVEXTERN int visOpenVents,visDummyVents,visOtherVents;
SVEXTERN int visOpenVentsAsOutline;
SVEXTERN int visTitle0, visTitle1, visTitle2;
SVEXTERN int ntitles,ititle;
SVEXTERN int visSmoke, visZone;
SVEXTERN int visEvac;
SVEXTERN int visBlocks;
SVEXTERN int smooth_block_solid;
SVEXTERN int visSmoothAsNormal,visTransparentBlockage;
SVEXTERN int visBlocksSave;
SVEXTERN int blocklocation;
SVEXTERN int ncadgeom;
SVEXTERN int visFloor, visFrame;
SVEXTERN int visNormalEditColors;
SVEXTERN int visWalls, visGrid, visCeiling, cursorPlot3D;
SVEXTERN int visVZone, visHZone, viszonefire;
SVEXTERN int visSensor, visSensorNorm, hasSensorNorm, visSprink, visHeat;
SVEXTERN int visVents;
SVEXTERN int partframestep, sliceframestep, boundframestep;
SVEXTERN int partframeskip, sliceframeskip, boundframeskip;
SVEXTERN int boundzipstep, boundzipskip;
SVEXTERN int smoke3dzipstep, smoke3dzipskip;
SVEXTERN int slicezipstep, slicezipskip;
SVEXTERN int isozipstep, isozipskip;
SVEXTERN int evacframeskip, evacframestep;
SVEXTERN int partpointstep;
SVEXTERN int partpointstep_old;
SVEXTERN int partpointskip;
SVEXTERN int viewoption;
SVEXTERN int xyz_clipplane,xyz_clipplane_last;
SVEXTERN int clip_x,clip_y,clip_z,clip_i,clip_j,clip_k;
SVEXTERN int clip_X,clip_Y,clip_Z,clip_I,clip_J,clip_K;
SVEXTERN float clip_x_val, clip_y_val, clip_z_val;
SVEXTERN float clip_X_val, clip_Y_val, clip_Z_val;
SVEXTERN int stepclip_x,stepclip_y,stepclip_z;
SVEXTERN int stepclip_X,stepclip_Y,stepclip_Z;
SVEXTERN float partpointsize,vectorpointsize,streaklinewidth;
SVEXTERN float isopointsize, isolinewidth;
SVEXTERN float plot3dpointsize, plot3dlinewidth;
SVEXTERN float vectorlinewidth;
SVEXTERN float sprinklerabssize, sensorabssize, heatabssize;
SVEXTERN float sensorrelsize,sensorrelsizeMIN;

SVEXTERN float linewidth, ventlinewidth, highlight_linewidth,solidlinewidth;
SVEXTERN float sliceoffset_factor, ventoffset_factor;
SVEXTERN int visBLOCKold;

SVEXTERN int selectedcolorbar_index,selectedcolorbar_index2;
SVEXTERN int planar_terrain_slice;
SVEXTERN int nrgb;
SVEXTERN int nrgb_ini;
SVEXTERN int nrgb2_ini;
SVEXTERN int rgb_white, rgb_yellow, rgb_blue, rgb_red;
SVEXTERN int rgb_green, rgb_magenta, rgb_cyan, rgb_black;
SVEXTERN int numColorbars;
SVEXTERN int setbw,colorbarflip,showcolorbarlines;
SVEXTERN int setbwSAVE;
SVEXTERN int background_flip;
SVEXTERN float transparentlevel;
SVEXTERN int SVDECL(transparentflag,1);
SVEXTERN int transparentflagSAVE;
SVEXTERN int antialiasflag;
SVEXTERN int nrgb_full;
SVEXTERN int nrgb_cad;
SVEXTERN float eyexfactor, eyeyfactor, eyezfactor;
SVEXTERN int transparent_state;
SVEXTERN float tload_begin, tload_end;
SVEXTERN int tload_skip;
SVEXTERN int use_tload_begin,use_tload_end,use_tload_skip;

SVEXTERN float frameinterval;

SVEXTERN int defaulttour_loaded;
SVEXTERN int blockages_dirty;
SVEXTERN int usetextures;
SVEXTERN int canrestorelastview;
SVEXTERN int ntargets;
SVEXTERN int endian_data, endian_native, setendian;

SVEXTERN int mainwindow_id,dwinWW;
SVEXTERN int rendertourcount;

SVEXTERN float vecyz[4];


SVEXTERN float ventcolor_orig[4];
SVEXTERN float *ventcolor;
SVEXTERN float static_color[4];
SVEXTERN float sensorcolor[4];
SVEXTERN float sensornormcolor[4];
SVEXTERN float sprinkoncolor[4];
SVEXTERN float sprinkoffcolor[4];
SVEXTERN float heatoncolor[4];
SVEXTERN float heatoffcolor[4];
SVEXTERN float backgroundbasecolor[4];
SVEXTERN float backgroundcolor[4];
SVEXTERN float foregroundbasecolor[4];
SVEXTERN float foregroundcolor[4];
SVEXTERN float boundcolor[4];
SVEXTERN float timebarcolor[4];
 
SVEXTERN float redcolor[4];

SVEXTERN int loadfiles_at_startup;

SVEXTERN char *smokeviewbindir;
SVEXTERN char *smokeviewtempdir;

SVEXTERN int nmenus;
SVEXTERN menudata menuinfo[10000];
SVEXTERN int showbuild;
SVEXTERN int max_screenWidth, max_screenHeight;
SVEXTERN int saveW, saveH;
SVEXTERN char *texturedir;
SVEXTERN char TITLE1[1024];
SVEXTERN char TITLE2[1024];//xxx check
SVEXTERN char TITLERELEASE[1024];
SVEXTERN char TITLE[1024];
SVEXTERN char FULLTITLE[1024];
SVEXTERN char *partshortlabel,*partunitlabel;
SVEXTERN char emptylabel[2];
SVEXTERN void *large_font;
SVEXTERN void *small_font;

SVEXTERN int nsmoothblocks,nopenvents,nopenvents_nonoutline,ndummyvents,ntransparentblocks,ntransparentvents;
SVEXTERN int ntotal_smooth_blockages;
SVEXTERN float veclengths[NVECLENGTHS];
SVEXTERN float texture_origin[3];

SVEXTERN unsigned char rgb_below_min[3], rgb_above_max[3];
SVEXTERN int colorbar_select_index,update_colorbar_select_index;
SVEXTERN float world_eyepos[3];
SVEXTERN int tour_usecurrent;
SVEXTERN int visVentLines, visVentSolid;
SVEXTERN int isZoneFireModel;
SVEXTERN int output_slicedata;
SVEXTERN f_units *unitclasses,*unitclasses_default,*unitclasses_ini;
SVEXTERN int nunitclasses,nunitclasses_default,nunitclasses_ini;
SVEXTERN mesh *meshinfo,*current_mesh, *mesh_save, *mesh_last, *loaded_isomesh;
SVEXTERN float devicenorm_length;
SVEXTERN int ndeviceinfo;
SVEXTERN int last_prop_display;
SVEXTERN device *deviceinfo;
SVEXTERN int show_smokesensors,active_smokesensors,test_smokesensors;
SVEXTERN float smoke3d_cvis;
SVEXTERN sv_object **object_defs;
SVEXTERN sv_object *heat_detector_object_backup;
SVEXTERN sv_object *target_object_backup;
SVEXTERN sv_object *sprinkler_upright_object_backup;
SVEXTERN sv_object *smoke_detector_object_backup;
SVEXTERN sv_object *thcp_object_backup;
SVEXTERN sv_object *missing_device;
SVEXTERN sv_object *error_device;
SVEXTERN sv_object object_def_first, object_def_last;
SVEXTERN char **device_texture_list;
SVEXTERN int ndevice_texture_list, *device_texture_list_index;
SVEXTERN int nobject_defs;
SVEXTERN int svofile_exists;
SVEXTERN treedata *treeinfo;
SVEXTERN terraindata *terraininfo;
SVEXTERN int ntreeinfo, nterraininfo, visTerrainType;
SVEXTERN float treecolor[4], treecharcolor[4], trunccolor[4];
SVEXTERN int showterrain;
SVEXTERN float rgb_terrain[10][4];
SVEXTERN tourdata *tourinfo;
SVEXTERN keyframe **tourknotskeylist;
SVEXTERN tourdata **tourknotstourlist;
SVEXTERN keyframe *selected_frame;
SVEXTERN tourdata *selected_tour;
SVEXTERN int callfrom_tourglui;
SVEXTERN int showtours_whenediting;

SVEXTERN int SVDECL(*slice_loadstack,NULL),  SVDECL(nslice_loadstack,0),  SVDECL(islice_loadstack,0);
SVEXTERN int SVDECL(*mslice_loadstack,NULL), SVDECL(nmslice_loadstack,0), SVDECL(imslice_loadstack,0);
SVEXTERN int SVDECL(*vslice_loadstack,NULL), SVDECL(nvslice_loadstack,0), SVDECL(ivslice_loadstack,0);
SVEXTERN int SVDECL(*mvslice_loadstack,NULL),SVDECL(nmvslice_loadstack,0),SVDECL(imvslice_loadstack,0);
SVEXTERN int SVDECL(*subslice_menuindex,NULL),SVDECL(*subvslice_menuindex,NULL);

SVEXTERN float xtimeleft, xtimeright;
SVEXTERN int showstereo, showstereoOLD, show_parallax, showstereo_frame;

SVEXTERN int SVDECL(show_hrrcutoff,1), SVDECL(show_hrrcutoff_active,0),SVDECL(hrrpuv_loaded,0);
SVEXTERN int showglui3dsmoke;
SVEXTERN int showgluitour;
SVEXTERN int showalert;
SVEXTERN int trainerview;
SVEXTERN int stereoactive;
SVEXTERN int apertureindex;
SVEXTERN int zoomindex;
SVEXTERN int projection_type;
SVEXTERN float apertures[5];
SVEXTERN float aperture,aperture_glui,aperture_default;
SVEXTERN float zooms[5];
SVEXTERN float zoom;
SVEXTERN int rgbmask[16];
SVEXTERN GLint nredbits, ngreenbits, nbluebits;
SVEXTERN int nredshift, ngreenshift, nblueshift;
SVEXTERN float xxmax,xyzbox;
SVEXTERN float xpltb[256], ypltb[256], zpltb[256];  
SVEXTERN float xplts[256*256], yplts[256*256], zplts[256*256];
SVEXTERN float *targtimes;
SVEXTERN int *targtimeslist;
SVEXTERN int *zonetlist;
SVEXTERN int delete_view_is_disabled;
SVEXTERN int old_listview;

SVEXTERN int sv_age;
SVEXTERN int titlesafe_offset;
SVEXTERN int titlesafe_offsetBASE;
SVEXTERN int   reset_frame;
SVEXTERN float reset_time,start_frametime,stop_frametime;
SVEXTERN int reset_time_flag;
SVEXTERN float velocity_range;
SVEXTERN int niso_compressed;
SVEXTERN int nslice_loaded, npatch_loaded;
SVEXTERN int *slice_loaded_list, *patch_loaded_list;
SVEXTERN int *render_frame;
SVEXTERN int RenderOnceNow, RenderOnceNowR, RenderOnceNowL;
SVEXTERN char *fdsprefix, *fdsprefix2;
SVEXTERN char *endianfilename;
SVEXTERN char *targfilename;

SVEXTERN int pass_through;
SVEXTERN int *sorted_surfidlist,*inv_sorted_surfidlist,nsorted_surfidlist;
SVEXTERN char *trainer_filename, *test_filename;
SVEXTERN char *filename_sb;
SVEXTERN int read_smoothobst;
SVEXTERN FILE *STREAM_SB;
#ifdef pp_MESSAGE
SVEXTERN int show_glui_warning,show_glui_error,show_glui_abort;
#endif
SVEXTERN time_t smv_modtime;
SVEXTERN float temp_threshold;
SVEXTERN char *smvfilename, *smvmenufile,*databasefilename,*smvprogdir,*smvisofilename;
SVEXTERN scriptfiledata first_scriptfile, last_scriptfile, *default_script;
SVEXTERN scriptdata *scriptinfo, *current_script_command;
SVEXTERN char *script_dir_path;
SVEXTERN int nscriptinfo;
SVEXTERN scriptfiledata *script_recording;
SVEXTERN int runscript, noexit;
SVEXTERN int script_multislice, script_multivslice, script_iso;
SVEXTERN FILE *scriptoutstream;
SVEXTERN char *scriptinifilename2;
SVEXTERN char *logfilename;
SVEXTERN char *flushfile, *chidfilebase;
SVEXTERN char *hrrfilename;
SVEXTERN hrrdata *hrrinfo;
SVEXTERN char *smokezippath;
SVEXTERN char *INI_fds_filein, *fds_filein, *fds_fileout,*fds_fileout2;
SVEXTERN char *caseinifilename;
SVEXTERN char *zonelonglabels, *zoneshortlabels, *zoneunits;
SVEXTERN char *smokeviewini;
SVEXTERN int overwrite_all,erase_all;
SVEXTERN int compress_autoloaded;
#ifdef WIN32
SVEXTERN   char openfilebuffer[1024];
SVEXTERN   int openfileflag;
#endif
SVEXTERN float xyzmaxdiff;
SVEXTERN char ext_png[5];
SVEXTERN char ext_jpg[5];
#ifdef pp_GDGIF
SVEXTERN char ext_gif[5];
#endif
SVEXTERN int renderfiletype;
SVEXTERN char part_ext[6];
SVEXTERN char ini_ext[5];

SVEXTERN int updatehiddenfaces;
SVEXTERN int nsurfids;
SVEXTERN surfid *surfids;
SVEXTERN int key_state;
SVEXTERN float starteyex, starteyey;
SVEXTERN float eye_xyz0[3];
SVEXTERN float start_xyz0[3];
SVEXTERN int glui_move_mode;

SVEXTERN float timeoffset;
SVEXTERN float motion_factor,speed_factor;
SVEXTERN float speed_desired;
SVEXTERN int speed_I;
SVEXTERN float speed_crawl,speed_walk,speed_now;
SVEXTERN int status_now,old_status_now;
SVEXTERN int motion_flag;
SVEXTERN int npartpoints, npartframes;
SVEXTERN float xslicemid, yslicemid, zslicemid;
SVEXTERN float delx;
SVEXTERN float delz;
SVEXTERN float d_eye_xyz[3],dsave_eye_xyz[3];
SVEXTERN float eyex0, eyey0, eyez0;
SVEXTERN float viewx, viewy, viewz;
SVEXTERN float anglexy0,direction_angle0;
SVEXTERN int xm0, ym0;
SVEXTERN int touring;
SVEXTERN int update_tourlist;
SVEXTERN float desired_view_height;
SVEXTERN int thistime, lasttime, resetclock,initialtime;
SVEXTERN int realtime_flag;
SVEXTERN char timelabel[30];
SVEXTERN char frameratelabel[30];
SVEXTERN char framelabel[30];
SVEXTERN float **p3levels, *zonelevels;
SVEXTERN float **p3levels256;
SVEXTERN char ***colorlabelp3,***colorlabeliso;
SVEXTERN char **scalep3, **scalep3copy;
SVEXTERN char *partscale;
SVEXTERN char a_partscale[31];
SVEXTERN char *zonescale;
SVEXTERN char a_zonescale[31];
SVEXTERN int islicetype,islicetype_save,ipatchtype;
SVEXTERN int iisotype,iisottype;
SVEXTERN char **colorlabelpart;
SVEXTERN char **colorlabelpatch;
SVEXTERN char **colorlabelzone;

SVEXTERN int minfill, maxfill;

//SVEXTERN int *xycolor, *xzcolor, *yzcolor;
//SVEXTERN float *xycolorf, *xzcolorf, *yzcolorf;

SVEXTERN int *plotiso;
SVEXTERN float *times,cputimes[20];
SVEXTERN int cpuframe;

SVEXTERN float eyexINI, eyeyINI, eyezINI;
SVEXTERN float anglexyINI, angleyzINI;
SVEXTERN float direction_angleINI;

SVEXTERN float eyexINI0_last, eyeyINI0_last, eyezINI0_last;
SVEXTERN float anglexyINI0_last, anglexzINI0_last, angleyzINI0_last;
SVEXTERN float direction_angleINI0_last;

SVEXTERN float eyexINI0, eyeyINI0, eyezINI0;
SVEXTERN float anglexyINI0, anglexzINI0, angleyzINI0;
SVEXTERN float direction_angleINI0;

SVEXTERN float xyzeyeorig[3],xeyedir[3], yeyedir[3], zeyedir[3];
SVEXTERN int adjustalphaflag;
SVEXTERN int colorband, show_extremedata;
SVEXTERN int show_extreme_above, show_extreme_below;

SVEXTERN int use_iblank,arg_iblank;

SVEXTERN int script_index, ini_index;
SVEXTERN char script_inifile_suffix[1024];
SVEXTERN char script_renderdir[1024], script_renderfilesuffix[1024], script_renderfile[1024];
SVEXTERN inifiledata first_inifile, last_inifile;
SVEXTERN char scriptinifilename[1024];
SVEXTERN int highlight_block, highlight_mesh, highlight_flag;
SVEXTERN int updatesmoothblocks,menusmooth,use_menusmooth;
SVEXTERN int smoothing_blocks;
SVEXTERN int blocksneedsmoothing;
SVEXTERN int updategetlabels;

SVEXTERN int pixel_skip;
SVEXTERN float smoke_extinct,smoke_dens,smoke_pathlength;
SVEXTERN int smoke_alpha;
#ifdef pp_SMOKETEST
SVEXTERN int smoketest,show_smoketest;
#else
SVEXTERN int smoketest,show_smoketest;
#endif
SVEXTERN int showall_textures;

SVEXTERN int ncolorbars;
SVEXTERN int ndefaultcolorbars;
SVEXTERN colorbardata *colorbarinfo,*current_colorbar;

SVEXTERN int update_load_startup;
SVEXTERN int do_threshold;
SVEXTERN int ntotal_blockages;
SVEXTERN int updateindexcolors;
SVEXTERN int show_path_knots;
SVEXTERN int keyframe_snap;
SVEXTERN int tourviewtype;
SVEXTERN int show_tourlocus;
SVEXTERN int tourlocus_type;
SVEXTERN int iavatar_types, navatar_types;
SVEXTERN int iavatar_evac;
SVEXTERN sv_object **avatar_types;
SVEXTERN int glui_avatar_index;
SVEXTERN sv_object *avatar_defs_backup[2];
SVEXTERN int device_sphere_segments;
SVEXTERN int ntexturestack;

SVEXTERN float tourrad_avatar;
SVEXTERN int dirtycircletour;
SVEXTERN float *tour_t, *tour_t2, *tour_dist, *tour_dist2, *tour_dist3;
SVEXTERN float view_tstart, view_tstop;
SVEXTERN int tour_constant_vel;
SVEXTERN float tour_bias,tour_continuity;
SVEXTERN int view_ntimes;
SVEXTERN int ntours,selectedtour_index,selectedtour_index_old,selectedtour_index_ini;
SVEXTERN int update_selectedtour_index;
SVEXTERN int viewtourfrompath,viewalltours,viewanytours,edittour;
SVEXTERN int rotation_index_OLD;
SVEXTERN selectdata *selectfaceinfo;
SVEXTERN blockagedata **selectblockinfo;
SVEXTERN tickdata *tickinfo;
SVEXTERN int nticks,ntickssmv;
SVEXTERN int visTicks;
SVEXTERN float user_tick_origin[3], user_tick_max[3], user_tick_min[3], user_tick_step[3], user_tick_length, user_tick_width;
SVEXTERN int user_tick_nxyz[3], user_tick_sub, user_tick_option, vis_user_ticks, auto_user_tick_placement;
SVEXTERN int user_tick_show_x, user_tick_show_y, user_tick_show_z;
SVEXTERN int visCadTextures, visTerrainTexture;
SVEXTERN int bw_colorbar_index;
SVEXTERN int viscolorbarpath;
SVEXTERN int *sortedblocklist,*changed_idlist,nchanged_idlist;
SVEXTERN int nselectblocks;
SVEXTERN surface *surfaceinfo,sdefault,v_surfacedefault,e_surfacedefault;
SVEXTERN int surface_indices[7],surface_indices_bak[7];
SVEXTERN int wall_case;
SVEXTERN surface *surfacedefault, *vent_surfacedefault, *exterior_surfacedefault;
SVEXTERN char surfacedefaultlabel[256];
SVEXTERN int nsurfaces;
SVEXTERN int ntotalfaces;
SVEXTERN colordata *firstcolor;
SVEXTERN texture *textureinfo, *terrain_texture;
SVEXTERN GLuint texture_colorbar_id, texture_slice_colorbar_id, texture_patch_colorbar_id, texture_plot3d_colorbar_id, texture_iso_colorbar_id;
SVEXTERN float mscale[3];
SVEXTERN float xclip_min, yclip_min, zclip_min;
SVEXTERN float xclip_max, yclip_max, zclip_max;
SVEXTERN float nearclip,farclip;
SVEXTERN int updateclipvals;
SVEXTERN int updateUpdateFrameRateMenu;
SVEXTERN int ntextures,ntextures_loaded_used;
SVEXTERN int nskyboxinfo;
SVEXTERN skyboxdata *skyboxinfo;
SVEXTERN firedata *fireinfo;
SVEXTERN roomdata *roominfo;
SVEXTERN zvent *zventinfo;
SVEXTERN zone *zoneinfo;
SVEXTERN zone *activezone;
SVEXTERN particle *partinfo;
SVEXTERN int update_screensize;
SVEXTERN int part5show;
SVEXTERN int streak5show,streak5value, streak5step, showstreakhead;
SVEXTERN int nstreak_value; // 5
SVEXTERN char *streak_values[5]; // "1","2","4","8","16"
SVEXTERN float streak_rvalue[7]; // 1.0, 2.0 4.0, 8.0, 16.0 
SVEXTERN int streak_index, update_streaks;       // 0
SVEXTERN float float_streak5value;// 1.0
SVEXTERN part5class *partclassinfo;
SVEXTERN int npartclassinfo;
SVEXTERN part5prop *part5propinfo, *current_property;
SVEXTERN int npart5prop,ipart5prop,ipart5prop_old;
SVEXTERN int prop_index;
SVEXTERN targ *targinfo;
SVEXTERN slice *sliceinfo;
SVEXTERN camdata *caminfo;
SVEXTERN multislice *multisliceinfo;
SVEXTERN multivslice *multivsliceinfo;
SVEXTERN outline *outlineinfo;
SVEXTERN int dummyvents;
SVEXTERN int noutlineinfo;
SVEXTERN int nmultislices;
SVEXTERN int nmultivslices;
SVEXTERN int *sliceorderindex,*vsliceorderindex,*partorderindex;
SVEXTERN int *patchorderindex,*isoorderindex,*plot3dorderindex;
SVEXTERN int showfiles;
SVEXTERN databounds SVDECL(*slicebounds,NULL), SVDECL(*isobounds,NULL), SVDECL(*patchbounds,NULL);
SVEXTERN vslice *vsliceinfo;
SVEXTERN int force_redisplay;
SVEXTERN int setp3min_temp, setp3max_temp;
SVEXTERN int setp3chopmin_temp, setp3chopmax_temp;
SVEXTERN float p3chopmin_temp, p3chopmax_temp;
SVEXTERN float p3min_temp, p3max_temp;

SVEXTERN int smoke3d_external;

SVEXTERN smoke3d *smoke3dinfo;
SVEXTERN int smoke_shade, fire_red, fire_green, fire_blue;
SVEXTERN float smoke_shade4[4];
SVEXTERN float fire_halfdepth;
SVEXTERN float hrrpuv_cutoff, global_hrrpuv_cutoff;

SVEXTERN int smokecullflag;
SVEXTERN int smokedrawtest,smokedrawtest2;
SVEXTERN int visMAINmenus;
SVEXTERN int smoke3d_thick;
#ifdef pp_GPU
SVEXTERN float smoke3d_rthick;
#endif
SVEXTERN int smokedrawtest_nummin;
SVEXTERN int smokedrawtest_nummax;
SVEXTERN int ijkbarmax;
SVEXTERN int blockage_as_input, blockage_snapped;
SVEXTERN int show_cad_and_grid;
SVEXTERN int use_nistlogo;
SVEXTERN int benchmark_flag;
SVEXTERN labeldata *labelinfo;
SVEXTERN int *slicetypes, *isotypes, *vslicetypes, *patchtypes;
SVEXTERN plot3d *plot3dinfo;
SVEXTERN float *plot3dtimelist;
SVEXTERN int nplot3dtimelist;
SVEXTERN patch *patchinfo;
SVEXTERN iso *isoinfo;
SVEXTERN targpos *target_positions;

SVEXTERN blockagedata *bchighlight,*bchighlight_old;
SVEXTERN cadgeom *cadgeominfo;

SVEXTERN int smokediff;
SVEXTERN int render_size_index;
SVEXTERN int render_skip_index;
SVEXTERN int buffertype;
SVEXTERN int benchmark;
SVEXTERN int opengldefined;
SVEXTERN float bench_starttime,bench_stoptime;
SVEXTERN int nslicetypes;
SVEXTERN int nvslicetypes;
SVEXTERN int nisotypes;
SVEXTERN int *isosubmenus, nisosubmenus;
SVEXTERN int *loadpatchsubmenus, nloadpatchsubmenus;
SVEXTERN int npatchtypes;
SVEXTERN char SVDECL(**patchlabellist,NULL);
SVEXTERN int SVDECL(*patchlabellist_index,NULL);
SVEXTERN int *sliceindex, *isoindex;

SVEXTERN int have_vents_int;
SVEXTERN int nface_outlines, nface_textures, nface_transparent;
SVEXTERN int nface_normals_single, nface_normals_double, nface_transparent_double, nvent_transparent;
SVEXTERN int show_transparent_vents;
SVEXTERN int show_bothsides_int, show_bothsides_ext;
SVEXTERN float transparency_level;
SVEXTERN int transparency_override;
SVEXTERN facedata **face_transparent;
#ifdef INMAIN
  SVEXTERN float rgb_baseBASE[MAXRGB][4]=
{
  {0.000000, 0.000000, 1.000000},
  {0.000000, 0.281732, 0.959493},
  {0.000000, 0.540640, 0.841254},
  {0.000000, 0.755749, 0.654861},
  {0.000000, 0.909632, 0.415416},
  {0.000000, 0.989821, 0.142316},
  {0.142316, 0.989821, 0.000000},
  {0.415416, 0.909632, 0.000000},
  {0.654861, 0.755749, 0.000000},
  {0.841254, 0.540640, 0.000000},
  {0.959493, 0.281732, 0.000000},
  {1.000000, 0.000000, 0.000000}
};
  SVEXTERN float bw_baseBASE[MAXRGB][4]={
  {1,            1,                1},           
  {0.909090909,  0.909090909,      0.909090909},
  {0.818181818,  0.818181818,      0.818181818}, 
  {0.727272727,  0.727272727,      0.727272727}, 
  {0.636363636,  0.636363636,      0.636363636}, 
  {0.545454545,  0.545454545,      0.545454545}, 
  {0.454545455,  0.454545455,      0.454545455}, 
  {0.363636364,  0.363636364,      0.363636364}, 
  {0.272727273,  0.272727273,      0.272727273}, 
  {0.181818182,  0.181818182,      0.181818182}, 
  {0.090909091,  0.090909091,      0.090909091}, 
  {  0,            0,                0}         
};
  SVEXTERN float rgb2BASE[MAXRGB][3]={
  {1.0f, 1.0f, 1.0f}, /* white */
  {1.0f, 1.0f, 0.0f}, /* yellow */
  {0.0f, 0.0f, 1.0f}, /* blue */
  {1.0f, 0.0f, 0.0f}, /* red */
  {0.0f, 1.0f, 0.0f}, /* green */
  {1.0f, 0.0f, 1.0f}, /* magenta */
  {0.0f, 1.0f, 1.0f}, /* cyan */
  {0.0f, 0.0f, 0.0f}  /* black */
  };

  SVEXTERN float rgbhazardBASE[MAXRGB][4]={
  {0.0f, 0.0f, 1.0f,1.0}, /* blue */
  {0.0f, 1.0f, 0.0f,1.0}, /* green */
  {1.0f, 1.0f, 0.0f,1.0}, /* yellow */
  {1.0f, 0.0f, 1.0f,1.0}, /* magenta */
  {1.0f, 0.0f, 0.0f,1.0}  /* red */
  };

#else
  SVEXTERN float rgb_baseBASE[MAXRGB][4];
  SVEXTERN float bw_baseBASE[MAXRGB][4];
  SVEXTERN float rgb2BASE[MAXRGB][3];
  SVEXTERN float rgbhazardBASE[MAXRGB][4];
#endif


