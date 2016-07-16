#ifndef SMOKEVIEWVARS_H_DEFINED
#define SMOKEVIEWVARS_H_DEFINED
#include <time.h>

#include "MALLOC.h"
#ifdef CPP
#include "glui.h"
#endif
#include "datadefs.h"
#include "translate.h"
#include "csphere.h"
#include "smokeviewdefs.h"
#include "isodefs.h"
#include "contourdefs.h"
#include "histogram.h"
#include "flowfiles.h"
#ifndef CPP
#include <zlib.h>
#endif
#include "smokeheaders.h"
#include "threader.h"

SVEXTERN int SVDECL(movie_bitrate, 5000);
SVEXTERN int SVDECL(disable_reshape, 0);
SVEXTERN int SVDECL(nscreeninfo,26);
#ifdef pp_RENDER360_DEBUG
SVEXTERN int SVDECL(screenview, 0);
SVEXTERN int SVDECL(*screenvis,NULL);
#endif
SVEXTERN screendata SVDECL(*screeninfo,NULL);
SVEXTERN int SVDECL(nwidth360,1024), SVDECL(nheight360,512);
SVEXTERN unsigned int SVDECL(*screenmap360, NULL);
SVEXTERN int SVDECL(render_360, 0);

SVEXTERN int SVDECL(highlight_vertexdup, 0);
SVEXTERN int SVDECL(highlight_edge0, 0);
SVEXTERN int SVDECL(highlight_edge1, 0);
SVEXTERN int SVDECL(highlight_edge2, 0);
SVEXTERN int SVDECL(highlight_edgeother, 0);

SVEXTERN colorbardata SVDECL(*split_colorbar, NULL);
SVEXTERN int split_colorbar_index;
#ifdef INMAIN
SVEXTERN float splitvals[3]={-1.0,0.0,1.0};
#else
SVEXTERN float splitvals[3];
#endif
#ifdef INMAIN
SVEXTERN int colorsplit[12] = {  0,0,0,  64,64,255,  0,192,0,  255,0,0 };
#else
SVEXTERN int colorsplit[12];
#endif

SVEXTERN int SVDECL(show_zlevel, 0);
SVEXTERN float terrain_zlevel;
SVEXTERN float terrain_zmin, terrain_zmax;
SVEXTERN int SVDECL(show_texture_1dimage, 0), SVDECL(show_texture_2dimage, 0);
SVEXTERN int SVDECL(force_update_histograms, 1);
SVEXTERN float SVDECL(geom_vert_exag, 1.0);
SVEXTERN float SVDECL(geom_vecfactor, .030);
SVEXTERN int SVDECL(geom_ivecfactor, 30);
SVEXTERN int SVDECL(geom_outline_ioffset, 5);
SVEXTERN float SVDECL(geom_outline_offset,0.005);
SVEXTERN float SVDECL(geom_max_angle, 30.0), cos_geom_max_angle;
SVEXTERN int SVDECL(use_max_angle, 1);
SVEXTERN int SVDECL(update_setvents, 0);
#ifdef pp_SLICECOLORDEFER
SVEXTERN int SVDECL(use_set_slicecolor, 1);
#else
SVEXTERN int SVDECL(use_set_slicecolor, 0);
#endif
SVEXTERN int SVDECL(cvents_defined, 0);
#ifdef pp_SLICEDUP
SVEXTERN int SVDECL(slicedup_option , SLICEDUP_KEEPFINE);
SVEXTERN int SVDECL(vectorslicedup_option, SLICEDUP_KEEPALL);
SVEXTERN int SVDECL(nslicedups, 0);
#endif
SVEXTERN int SVDECL(vis_xtree, 0), SVDECL(vis_ytree, 0), SVDECL(vis_ztree, 1);
SVEXTERN int SVDECL(max_device_tree,0);
#ifdef INMAIN
SVEXTERN float northangle_position[3] = {0.0, 0.0, 0.1};
#else
SVEXTERN float northangle_position[3];
#endif
SVEXTERN float SVDECL(northangle, 0.0);
SVEXTERN int SVDECL(vis_northangle, 0), SVDECL(have_northangle,0);
#ifdef pp_PILOT
SVEXTERN int SVDECL(npilot_buckets, 8);
#ifdef pp_WINDROSE
SVEXTERN int SVDECL(npilot_nr, 8);
SVEXTERN int SVDECL(npilot_ntheta, 12);
#endif
SVEXTERN int SVDECL(pilot_viewtype, 0);
#endif
SVEXTERN int SVDECL(ngeomdiaginfo, 0), SVDECL(show_geometry_diagnostics,0);
SVEXTERN geomdiagdata SVDECL(*geomdiaginfo,NULL);
SVEXTERN int SVDECL(zone_rho, 1);
SVEXTERN int SVDECL(visventslab, 0), SVDECL(visventprofile,1);
SVEXTERN int SVDECL(update_readiso_geom_wrapup, UPDATE_ISO_OFF);
SVEXTERN int SVDECL(nmemory_ids, 0);
SVEXTERN int SVDECL(update_playmovie, 0);
SVEXTERN int SVDECL(play_movie_now, 1);
SVEXTERN int SVDECL(update_makemovie, 0),SVDECL(movie_filetype,AVI);
SVEXTERN char movie_name[1024], movie_ext[10], render_file_base[1024];
SVEXTERN int SVDECL(movie_framerate, 10), SVDECL(have_ffmpeg, 0), SVDECL(have_ffplay, 0), SVDECL(overwrite_movie, 1);

SVEXTERN int SVDECL(show_missing_objects, 1),SVDECL(have_missing_objects,0);
SVEXTERN int SVDECL(toggle_dialogs, 1);
SVEXTERN int SVDECL(script_render_width, 320), SVDECL(script_render_height, 240);
SVEXTERN int SVDECL(show_tetratest_labels, 1);
SVEXTERN float SVDECL(tetra_line_thickness, 2.0);
SVEXTERN float SVDECL(tetra_point_size, 10.0);
SVEXTERN int SVDECL(use_data_extremes, 1);
SVEXTERN int SVDECL(extreme_data_offset,1), SVDECL(colorbar_offset,0), SVDECL(colorbarflip,0);

#ifdef INMAIN
SVEXTERN float gvecphys[3]={0.0,0.0,-9.8};
SVEXTERN float gvecunit[3]={0.0,0.0,-1.0};
#else
SVEXTERN float gvecphys[3];
SVEXTERN float gvecunit[3];
#endif
SVEXTERN int SVDECL(update_have_gvec,0),SVDECL(gvec_down,1),SVDECL(have_gvec,0),SVDECL(changed_zaxis,0),SVDECL(showgravity,0);
SVEXTERN float SVDECL(slice_line_contour_width,1.0);
SVEXTERN int SVDECL(slice_contour_type,0);
SVEXTERN int SVDECL(viscadopaque,0);
SVEXTERN int SVDECL(structured_isopen,0), SVDECL(unstructured_isopen,0);
SVEXTERN float SVDECL(patchout_tmin,1.0), SVDECL(patchout_tmax,-1.0);
SVEXTERN float SVDECL(patchout_xmin,1.0), SVDECL(patchout_xmax,-1.0);
SVEXTERN float SVDECL(patchout_ymin,1.0), SVDECL(patchout_ymax,-1.0);
SVEXTERN float SVDECL(patchout_zmin,1.0), SVDECL(patchout_zmax,-1.0);
SVEXTERN int SVDECL(showpatch_both,0);
SVEXTERN int SVDECL(show_triangle_count,0);
SVEXTERN int SVDECL(n_geom_triangles,0);
SVEXTERN int SVDECL(show_device_orientation,0);
SVEXTERN float SVDECL(orientation_scale,1.0);
SVEXTERN char SVDECL(*script_labelstring,NULL);
SVEXTERN char SVDECL(*loaded_file,NULL);
SVEXTERN int SVDECL(clipon,0);
SVEXTERN int SVDECL(vectortype,0);
SVEXTERN float tetra_xyz[3];
SVEXTERN int SVDECL(show_test_in_tetra,0);
SVEXTERN int SVDECL(show_cutcells,0);
SVEXTERN int b_state[7],SVDECL(*box_state,b_state+1);
SVEXTERN int face_id[200],face_vis[10], face_vis_old[10];
SVEXTERN int SVDECL(update_volbox_controls,0);
SVEXTERN float SVDECL(face_factor,0.01);
SVEXTERN int SVDECL(have_volume,0);
SVEXTERN int SVDECL(show_volumes_interior,0), SVDECL(show_volumes_exterior,1);
SVEXTERN int SVDECL(show_faces_interior,0), SVDECL(show_faces_exterior,1);
SVEXTERN int SVDECL(show_volumes_solid,1);
SVEXTERN int SVDECL(show_volumes_outline,0);
SVEXTERN int SVDECL(show_slices_and_vectors,0);
SVEXTERN int SVDECL(vispilot,0);
SVEXTERN int SVDECL(compute_fed,0);
SVEXTERN int SVDECL(is_fed_colorbar, 0);
SVEXTERN int SVDECL(tour_global_tension_flag,1);
SVEXTERN float SVDECL(tour_global_tension,0.0);

SVEXTERN float box_bounds[6],box_bounds2[6],box_translate[3],tetra_vertices[12];
SVEXTERN int tetrabox_vis[10];
SVEXTERN int SVDECL(geomtest_option,NO_TEST);

SVEXTERN int SVDECL(convert_ini,0), SVDECL(convert_ssf,0);
SVEXTERN int SVDECL(update_ssf,0);
SVEXTERN char SVDECL(*ini_from,NULL), SVDECL(*ini_to,NULL);
SVEXTERN char SVDECL(*ssf_from, NULL), SVDECL(*ssf_to, NULL);

SVEXTERN int SVDECL(cache_boundarydata, 0);
SVEXTERN int SVDECL(tour_antialias,0);
SVEXTERN int SVDECL(tour_drag,0);

SVEXTERN int SVDECL(update_gslice,0);
SVEXTERN int SVDECL(wc_flag,0);
SVEXTERN circdata cvent_circ, object_circ;
#ifdef pp_BETA
SVEXTERN int SVDECL(show_all_units,1);
#else
SVEXTERN int SVDECL(show_all_units,0);
#endif
SVEXTERN int SVDECL(circle_outline,0);
SVEXTERN int SVDECL(gversion,0);
SVEXTERN unsigned char SVDECL(*patchmin_unit,NULL),SVDECL(*patchmax_unit,NULL);
SVEXTERN unsigned char SVDECL(*slicemin_unit,NULL),SVDECL(*slicemax_unit,NULL);
SVEXTERN unsigned char SVDECL(*plot3dmin_unit,NULL),SVDECL(*plot3dmax_unit,NULL);
SVEXTERN unsigned char SVDECL(*partmin_unit,NULL),SVDECL(*partmax_unit,NULL);
SVEXTERN char degC[3], degF[3];
SVEXTERN float SVDECL(tmax_part,16.0);
SVEXTERN int SVDECL(redirect,0);
SVEXTERN int SVDECL(tempdir_flag,0),SVDECL(time_flag,0);
SVEXTERN int SVDECL(script_render,0);
SVEXTERN int SVDECL(make_volrender_script,0);
SVEXTERN char SVDECL(*volrender_scriptname,NULL);
SVEXTERN float SVDECL(nongpu_vol_factor,1.0);
SVEXTERN float SVDECL(gpu_vol_factor,1.0);
SVEXTERN int SVDECL(disable_gpu,0);
SVEXTERN int SVDECL(render_state,0);
SVEXTERN int SVDECL(script_startframe,-1), SVDECL(script_skipframe,-1);
SVEXTERN int SVDECL(vol_startframe0,-1), SVDECL(vol_skipframe0,-1);
SVEXTERN int SVDECL(startframe0,-1), SVDECL(skipframe0,-1);
SVEXTERN int SVDECL(skip_render_frames,0);
SVEXTERN int SVDECL(update_smokecolorbar,0);
SVEXTERN int SVDECL(combine_meshes,1);
SVEXTERN int colorbar_left_pos, colorbar_right_pos, colorbar_down_pos, colorbar_top_pos;
SVEXTERN float scale_2d_x, scale_2d_y;
SVEXTERN int SVDECL(colorbar_delta,35);
SVEXTERN int colorbar_label_width;

SVEXTERN int timebar_left_width, timebar_right_width;
SVEXTERN int SVDECL(h_space,2), SVDECL(v_space,2);
SVEXTERN portdata VP_fullscreen, VP_title, VP_timebar, VP_colorbar, VP_scene, VP_info;
SVEXTERN int SVDECL(in_external,0);
SVEXTERN int SVDECL(label_list_index,0);
SVEXTERN labeldata LABEL_local, SVDECL(*LABEL_global_ptr,NULL), LABEL_default;

SVEXTERN int SVDECL(renderdoublenow,0);
SVEXTERN int SVDECL(nrender_rows,2);
SVEXTERN int port_pixel_width, port_pixel_height;
SVEXTERN float port_unit_width, port_unit_height;
SVEXTERN int SVDECL(scaled_font2d_height,12);
SVEXTERN float SVDECL(scaled_font2d_height2width,1.0);
SVEXTERN int SVDECL(scaled_font3d_height,32);
SVEXTERN float SVDECL(scaled_font3d_height2width,1.0);
SVEXTERN float quat_general[4], quat_rotation[16];

SVEXTERN float modelview_identity[16];
SVEXTERN mousedata mouseinfo;
SVEXTERN int SVDECL(use_glui_rotate,0);
SVEXTERN int SVDECL(show_fed_area,1);
SVEXTERN int SVDECL(*fed_areas,NULL);
SVEXTERN char default_fed_colorbar[255];

SVEXTERN int SVDECL(*meshvisptr,NULL);
SVEXTERN smoke3ddata SVDECL(**smoke3dinfo_sorted,NULL);
SVEXTERN int SVDECL(from_commandline,0);
SVEXTERN filelistdata SVDECL(*ini_filelist,NULL);
SVEXTERN int SVDECL(nini_filelist,0);
SVEXTERN float this_mouse_time, SVDECL(last_mouse_time,0.0);
SVEXTERN int move_gslice;

SVEXTERN int SVDECL(visGeomTextures,0);
SVEXTERN int nplotx_all, nploty_all, nplotz_all;
SVEXTERN int iplotx_all, iploty_all, iplotz_all;
SVEXTERN int SVDECL(iplot_state,0);
SVEXTERN int SVDECL(visx_all,0),SVDECL(visy_all,1),SVDECL(visz_all,0);
SVEXTERN float SVDECL(*plotx_all,NULL), SVDECL(*ploty_all,NULL), SVDECL(*plotz_all,NULL);
SVEXTERN int SVDECL(defer_file_loading,0);
SVEXTERN int SVDECL(regenerate_fed,0);
SVEXTERN int SVDECL(debug_count,0);
SVEXTERN geomdata SVDECL(**geominfoptrs,NULL);
SVEXTERN int SVDECL(ngeominfoptrs,0);
SVEXTERN int SVDECL(update_glui_wui,0);
SVEXTERN int SVDECL(update_glui_stereo,0);
SVEXTERN int SVDECL(update_glui_trainer,0);
SVEXTERN int SVDECL(update_glui_alert,0);
SVEXTERN int SVDECL(update_glui_tour,0);
SVEXTERN int SVDECL(update_glui_motion,0);
SVEXTERN int SVDECL(update_glui_message,0);
SVEXTERN int SVDECL(update_glui_labels,0);
SVEXTERN int SVDECL(update_glui_device,0);
SVEXTERN int SVDECL(update_glui_clip,0);
SVEXTERN int SVDECL(update_glui_geometry,0);
SVEXTERN int SVDECL(update_glui_colorbar,0);
SVEXTERN int SVDECL(update_glui_bounds,0);
SVEXTERN int SVDECL(update_glui_shooter,0);

SVEXTERN int SVDECL(update_glui_dialogs,0);
#ifdef pp_LANG
SVEXTERN langlistdata SVDECL(*langlistinfo,NULL);
SVEXTERN int nlanglistinfo,SVDECL(show_lang_menu,1);
SVEXTERN char startup_lang_code[3];
#endif

#ifdef pp_GPUTHROTTLE
  SVEXTERN float SVDECL(thisGPUtime,0.0), SVDECL(lastGPUtime,0.0);
  SVEXTERN float SVDECL(thisMOTIONtime,0.0), SVDECL(lastMOTIONtime,0.0);
  SVEXTERN int SVDECL(GPUnframes,0),SVDECL(MOTIONnframes,0);
#endif
SVEXTERN int SVDECL(mouse_down,0);
SVEXTERN int SVDECL(show_volsmoke_moving,0);
SVEXTERN int SVDECL(freeze_volsmoke,0);
SVEXTERN int SVDECL(show_iso_solid,1),SVDECL(show_iso_outline,1),SVDECL(show_iso_verts,0);
SVEXTERN int SVDECL(show_patch_solid, 1), SVDECL(show_patch_outline, 0), SVDECL(show_patch_verts, 0);
SVEXTERN int SVDECL(show_patch_ingas, 1), SVDECL(show_patch_insolid, 1), SVDECL(show_patch_incutcell, 1);
SVEXTERN int SVDECL(show_iso_normal, 0), SVDECL(smooth_iso_normal, 1);
SVEXTERN int SVDECL(show_faces_solid, 1), SVDECL(show_faces_outline, 1), SVDECL(show_geom_verts, 0);
SVEXTERN int SVDECL(show_geom_normal, 0), SVDECL(smooth_geom_normal, 1);
SVEXTERN geomlistdata SVDECL(*geomlistinfo, NULL);
SVEXTERN int SVDECL(have_volcompressed,0);
SVEXTERN int SVDECL(glui_load_volcompressed,0),SVDECL(load_volcompressed,0);
SVEXTERN int SVDECL(use_multi_threading,1);
SVEXTERN int SVDECL(load_at_rendertimes,1);
SVEXTERN int nvolrenderinfo;
SVEXTERN int SVDECL(compress_volsmoke,0),SVDECL(glui_compress_volsmoke,0);
SVEXTERN int SVDECL(read_vol_mesh,VOL_READNONE);
SVEXTERN int SVDECL(trainer_temp_index,0),SVDECL(trainer_oxy_index,0);
SVEXTERN int SVDECL(*trainer_temp_indexes,NULL),SVDECL(*trainer_oxy_indexes,NULL);
SVEXTERN int SVDECL(trainer_showall_mslice,0),SVDECL(trainer_cycle_mslice,1);
SVEXTERN int SVDECL(trainer_temp_n,0),SVDECL(trainer_oxy_n,0);
SVEXTERN char SVDECL(*tr_name,NULL);
SVEXTERN int SVDECL(show_smoke_lighting,0),SVDECL(have_lighting,0);
SVEXTERN int SVDECL(showdeviceval,0),SVDECL(showvdeviceval,0),SVDECL(colordeviceval,0);
SVEXTERN int SVDECL(showdevicetype,1), SVDECL(showdeviceunit,1);
SVEXTERN float SVDECL(device_valmin,0.0), SVDECL(device_valmax,1.0);
SVEXTERN devicedata SVDECL(**devicetypes,NULL);
SVEXTERN int SVDECL(ndevicetypes,0);
SVEXTERN int SVDECL(sort_geometry,1),SVDECL(sort_transparent_faces,0);
SVEXTERN isotri SVDECL(***iso_trans_list,NULL),SVDECL(***iso_opaques_list,NULL);
SVEXTERN int SVDECL(*niso_trans_list,NULL),SVDECL(*niso_opaques_list,NULL);
SVEXTERN int SVDECL(niso_timesteps,0);
SVEXTERN isotri SVDECL(**iso_trans,NULL),SVDECL(**iso_opaques,NULL);
SVEXTERN int SVDECL(niso_trans,0),SVDECL(niso_opaques,0);
SVEXTERN int SVDECL(sort_iso_triangles,1);
SVEXTERN int SVDECL(object_outlines,0);
SVEXTERN int SVDECL(usemenu,1),SVDECL(show_evac_slices,0);
SVEXTERN float direction_color[4], SVDECL(*direction_color_ptr,NULL);
SVEXTERN int SVDECL(constant_evac_coloring,1),SVDECL(data_evac_coloring,1),SVDECL(show_evac_colorbar,0);
SVEXTERN float hrrpuv_iso_color[4];
SVEXTERN int show_slice_terrain;
SVEXTERN int npropinfo;
SVEXTERN propdata SVDECL(*propinfo,NULL);
SVEXTERN float right_green, right_blue;

SVEXTERN int SVDECL(levelset_colorbar,-1), SVDECL(wallthickness_colorbar,-1);
SVEXTERN colorbardata SVDECL(*fire_colorbar,NULL);
SVEXTERN float SVDECL(glui_time,0.0);
SVEXTERN int show_mode;
SVEXTERN int SVDECL(cellcenter_slice_active,0), SVDECL(cellcenter_bound_active,0), SVDECL(facecenter_slice_active,0);
SVEXTERN int SVDECL(part5colorindex,0), SVDECL(show_tracers_always,0);
SVEXTERN int navatar_colors;
SVEXTERN int select_avatar, selected_avatar_tag, view_from_selected_avatar;
SVEXTERN int select_device, selected_device_tag;
SVEXTERN float selected_avatar_pos[3], selected_avatar_angle;
SVEXTERN unsigned char select_device_color[4], SVDECL(*select_device_color_ptr,NULL);
SVEXTERN float SVDECL(*avatar_colors,NULL);
SVEXTERN int SVDECL(script_render_flag,0), SVDECL(script_itime,0);

SVEXTERN int SVDECL(show_slice_in_obst,0), offset_slice;
SVEXTERN int skip_slice_in_embedded_mesh;
SVEXTERN int n_embedded_meshes;

SVEXTERN geomdata SVDECL(*geominfo,NULL);
SVEXTERN int SVDECL(ngeominfo,0);

SVEXTERN int npartframes_max;
SVEXTERN int force_isometric;
SVEXTERN int SVDECL(update_startup_view,0);
SVEXTERN int SVDECL(render_multi,0);
SVEXTERN int SVDECL(render_multi_state,0);
SVEXTERN int SVDECL(render_multi_menu, 0);
SVEXTERN int SVDECL(render_from_menu,0);
SVEXTERN int SVDECL(usetexturebar,1);
SVEXTERN int show_smokelighting;
SVEXTERN int SVDECL(cullgeom_portsize,16);
SVEXTERN int SVDECL(update_initcullgeom,1),SVDECL(cullgeom,1);
#ifdef pp_CULL
SVEXTERN int cullactive, SVDECL(show_cullports,0), SVDECL(cull_portsize,16);
SVEXTERN int cullsmoke, ncullplaneinfo;
SVEXTERN cullplanedata SVDECL(*cullplaneinfo,NULL);
SVEXTERN cullplanedata SVDECL(**sort_cullplaneinfo,NULL);
SVEXTERN int have_setpixelcount,update_initcullplane;
#endif
SVEXTERN int opengl_version;
SVEXTERN char opengl_version_label[256];

SVEXTERN int SVDECL(usevolrender,1);
SVEXTERN int SVDECL(usegpu,0),SVDECL(gpuactive,0);
#ifdef pp_GPU
SVEXTERN int GPU_aspectratio;
SVEXTERN int GPU_smoke3d_rthick, GPU_skip, GPU_hrrcutoff, GPU_hrr, GPU_hrrpuv_max_smv, GPU_hrrpuv_cutoff, GPU_smoke_albedo;
SVEXTERN int GPU_fire_alpha, GPU_firecolor, GPU_have_smoke, GPU_smokecolormap;
SVEXTERN int GPU_smokeshade,GPU_smokealpha;
SVEXTERN int GPU_adjustalphaflag;

SVEXTERN int GPUzone_zonedir;
SVEXTERN int GPUzone_zoneinside;
SVEXTERN int GPUzone_eyepos;
SVEXTERN int GPUzone_xyzmaxdiff;
SVEXTERN int GPUzone_boxmin, GPUzone_boxmax;
SVEXTERN int GPUzone_zlay;
SVEXTERN int GPUzone_odl, GPUzone_odu;

SVEXTERN int GPUvol_inside, GPUvol_eyepos, GPUvol_xyzmaxdiff, GPUvol_slicetype,GPUvol_dcell3;
SVEXTERN int GPUvol_gpu_vol_factor;
SVEXTERN int GPUvol_soot_density, GPUvol_fire, GPUvol_blockage;
SVEXTERN int GPUvol_fire_opacity_factor,GPUvol_volbw,GPUvol_mass_extinct;
SVEXTERN int GPUvol_temperature_min,GPUvol_temperature_cutoff,GPUvol_temperature_max;
SVEXTERN int GPUvol_boxmin, GPUvol_boxmax, GPUvol_drawsides;
SVEXTERN int GPUvol_smokecolormap, GPUvol_dcell, GPUvol_havefire;

SVEXTERN int GPU3dslice_valtexture,GPU3dslice_colormap;
SVEXTERN int GPU3dslice_val_min,GPU3dslice_val_max;
SVEXTERN int GPU3dslice_boxmin, GPU3dslice_boxmax;
SVEXTERN int GPU3dslice_transparent_level;
SVEXTERN int GPUvol_block_volsmoke;

#ifdef pp_GPUDEPTH
SVEXTERN int GPUvol_depthtexture, GPUvol_screensize,GPUvol_nearfar;
SVEXTERN GLuint SVDECL(depthtexture_id,0);
#endif

#endif
SVEXTERN int SVDECL(ncsvinfo,0);
SVEXTERN csvdata SVDECL(*csvinfo,NULL);
SVEXTERN int smoke_render_option;
SVEXTERN float fnear, ffar;
SVEXTERN float partfacedir[3];
SVEXTERN int SVDECL(demo_option,0);
SVEXTERN int small_font_height, large_font_height;
SVEXTERN float cb_valmin, cb_valmax, cb_val;
SVEXTERN int cb_colorindex;
SVEXTERN float rgbterrain[4*MAXRGB];
SVEXTERN int terrain_rgba_zmin[3];
SVEXTERN int terrain_rgba_zmax[3];
SVEXTERN float vertical_factor;

SVEXTERN char input_filename_ext[4];

SVEXTERN float percentile_level;
SVEXTERN float SVDECL(fire_line_min,150.0), SVDECL(fire_line_max,200.0);
SVEXTERN int SVDECL(update_fire_line,0);
SVEXTERN int SVDECL(fire_line_index,-1);
SVEXTERN int SVDECL(slice_bounds_dialog,1);

SVEXTERN float xtemp;

SVEXTERN float set_view_xyz[3],user_zaxis[3];
#ifdef INMAIN
  SVEXTERN float zaxis_angles[3]={0.000000, 90.000000, 0.000000};
#else
  SVEXTERN float zaxis_angles[3];
#endif

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
SVEXTERN float SVDECL(*mat_ambient2,NULL), SVDECL(*mat_specular2,NULL);

SVEXTERN GLfloat iso_specular[4];
SVEXTERN GLfloat iso_shininess;

SVEXTERN float block_ambient_orig[4];
SVEXTERN float SVDECL(*block_ambient2,NULL);
SVEXTERN float block_specular_orig[4];
SVEXTERN float SVDECL(*block_specular2,NULL);
SVEXTERN GLfloat SVDECL(block_shininess,100.0);

SVEXTERN int SVDECL(light_enabled0, 1);
SVEXTERN int SVDECL(light_enabled1, 1);
SVEXTERN GLfloat light_position0[4];
SVEXTERN GLfloat light_position1[4];

SVEXTERN int SVDECL(lightmodel_localviewer,0);
SVEXTERN int SVDECL(lightmodel_separatespecularcolor,0);
SVEXTERN GLfloat ambientlight[4];
SVEXTERN GLfloat diffuselight[4];

SVEXTERN int list_p3_index,list_slice_index,list_patch_index,list_iso_index;
SVEXTERN int list_p3_index_old, list_slice_index_old, list_patch_index_old,list_iso_index_old;

SVEXTERN float glui_block_xmin, glui_block_ymin, glui_block_zmin;
SVEXTERN float glui_block_xmax, glui_block_ymax, glui_block_zmax;

SVEXTERN int SVDECL(nzonetotal,0);
SVEXTERN float SVDECL(zoneglobalmin,0.0), SVDECL(zoneglobalmax,0.0);
SVEXTERN float SVDECL(zoneusermin,0.0), SVDECL(zoneusermax,0.0);
SVEXTERN float zonelevels256[256];
SVEXTERN float boundarylevels256[256];
SVEXTERN float partlevels256[256];
SVEXTERN float SVDECL(*zone_times,NULL), SVDECL(*zoneylay,NULL), SVDECL(*zonetl,NULL), SVDECL(*zonetu,NULL), SVDECL(*zonepr,NULL);
SVEXTERN float SVDECL(*zonerhol, NULL), SVDECL(*zonerhou, NULL);
SVEXTERN float SVDECL(*zoneqfire,NULL), SVDECL(*zonefheight,NULL), SVDECL(*zonefbase,NULL), SVDECL(*zonefdiam,NULL);
SVEXTERN float SVDECL(*zoneodl,NULL), SVDECL(*zoneodu,NULL), SVDECL(*zonevents,NULL);
SVEXTERN float SVDECL(maxslabflow, 0.0);
SVEXTERN int SVDECL(have_ventslab_flow,0);
SVEXTERN float SVDECL(*zoneslab_T, NULL), SVDECL(*zoneslab_F, NULL), SVDECL(*zoneslab_YB, NULL), SVDECL(*zoneslab_YT, NULL);
SVEXTERN int SVDECL(*zoneslab_n, NULL);
SVEXTERN int SVDECL(zonecsv, 0), SVDECL(nzvents, 0), SVDECL(nzhvents, 0), SVDECL(nzvvents, 0), SVDECL(nzmvents, 0);
SVEXTERN float zone_maxventflow;
SVEXTERN unsigned char SVDECL(*hazardcolor,NULL);
SVEXTERN float SVDECL(zone_ventfactor,1.0);
SVEXTERN unsigned char SVDECL(*izonetu,NULL);
SVEXTERN int nzone_times;
SVEXTERN float barright;
SVEXTERN float SVDECL(*tspr,NULL);
SVEXTERN float tmin_global, tmax_global;

SVEXTERN int videoSTEREO;
SVEXTERN float fzero;

SVEXTERN char blank_global[2];

SVEXTERN float SVDECL(*sphere_xyz,NULL);
SVEXTERN int demo_mode;
SVEXTERN int update_demo;
SVEXTERN int mxplot3dvars;
SVEXTERN int loadplot3dall;
SVEXTERN char *shortp3label[MAXPLOT3DVARS], *unitp3label[MAXPLOT3DVARS];
SVEXTERN char SVDECL(*LESsystem,NULL),SVDECL(*LESendian,NULL);

SVEXTERN int show3dsmoke;
SVEXTERN float frustum[6][4];
SVEXTERN int showtime, showtime2, showplot3d, showpatch, showslice, showvslice, showsmoke, showzone, showiso, showevac;
SVEXTERN int SVDECL(showvolrender,0);
SVEXTERN int vis_slice_contours;
SVEXTERN int update_slicecontours;
SVEXTERN int showevac_colorbar;
SVEXTERN int showiso_colorbar;
SVEXTERN int SVDECL(visgridloc,0);
SVEXTERN int valindex;

SVEXTERN int fire_colorbar_index,SVDECL(fire_colorbar_index_save,-1);
SVEXTERN int SVDECL(update_fire_colorbar_index,0);
SVEXTERN int SVDECL(fire_colorbar_index_ini,0);
SVEXTERN float SVDECL(*rgb2_ini,NULL);
SVEXTERN float rgb_full[MAXRGB][4];
SVEXTERN float rgb_full2[MAXRGB][4];
SVEXTERN float rgb_terrain2[4 * MAXRGB];
SVEXTERN float rgb_slice[4 * MAXRGB];
SVEXTERN float rgb_volsmokecolormap[4*MAXSMOKERGB];
SVEXTERN float rgb_slicesmokecolormap[4*MAXSMOKERGB];
SVEXTERN float rgb_iso[4*MAXRGB];
SVEXTERN float rgb_patch[4*MAXRGB];
SVEXTERN float rgb_plot3d[4*MAXRGB];
SVEXTERN float rgb_part[4*MAXRGB];
SVEXTERN float rgb_trans[4*MAXRGB];
SVEXTERN float rgb_cad[MAXRGB][4];

SVEXTERN float iso_transparency, SVDECL(*iso_colors,NULL), SVDECL(*iso_colorsbw,NULL);
SVEXTERN int glui_iso_colors[4], SVDECL(glui_iso_level,1), glui_iso_transparency;

SVEXTERN float SVDECL(*rgb_ini,NULL);
SVEXTERN float rgb[MAXRGB][4];
SVEXTERN float mouse_deltax, mouse_deltay;
SVEXTERN float SVDECL(**rgbptr,NULL), SVDECL(**rgb_plot3d_contour,NULL);
SVEXTERN float char_color[4];
SVEXTERN float *rgb_step[255];
SVEXTERN float movedir[3];
SVEXTERN float rgb_base[MAXRGB][4];
SVEXTERN float bw_base[MAXRGB][4];
SVEXTERN int nrgb2;
SVEXTERN float rgb2[MAXRGB][3];
SVEXTERN float inverse_modelview_setup[16];
SVEXTERN float modelview_setup[16];
SVEXTERN float modelview_rotate_last[16],modelview_rotate_save[16];
SVEXTERN float modelview_current[16];
SVEXTERN float modelview_scratch[16];

SVEXTERN cameradata SVDECL(*camera_current,NULL), SVDECL(*camera_save,NULL), SVDECL(*camera_last,NULL);
SVEXTERN cameradata SVDECL(*camera_external,NULL), SVDECL(*camera_internal,NULL);
SVEXTERN cameradata SVDECL(*camera_ini,NULL), SVDECL(*camera_external_save,NULL);
SVEXTERN cameradata camera_list_first, camera_list_last, SVDECL(**camera_list,NULL);
SVEXTERN int ncamera_list,i_view_list,init_camera_list_flag;
SVEXTERN int camera_max_id;
SVEXTERN int startup,startup_view_ini,selected_view;
SVEXTERN char label_startup_view[256];
SVEXTERN char SVDECL(*camera_label,NULL), SVDECL(*colorbar_label,NULL);

SVEXTERN int visPatchType[7];
SVEXTERN int p3_extreme_min[MAXPLOT3DVARS], p3_extreme_max[MAXPLOT3DVARS];

SVEXTERN int setp3min[MAXPLOT3DVARS], setp3min_save[MAXPLOT3DVARS];
SVEXTERN float p3min[MAXPLOT3DVARS], p3min_save[MAXPLOT3DVARS];

SVEXTERN int setp3max[MAXPLOT3DVARS], setp3max_save[MAXPLOT3DVARS];
SVEXTERN float p3max[MAXPLOT3DVARS], p3max_save[MAXPLOT3DVARS];

SVEXTERN int setp3chopmin[MAXPLOT3DVARS], setp3chopmax[MAXPLOT3DVARS];
SVEXTERN float p3chopmin[MAXPLOT3DVARS], p3chopmax[MAXPLOT3DVARS];

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
SVEXTERN int SVDECL(trainerload,0),SVDECL(trainerload_old,0);
SVEXTERN int fontsize_save;
SVEXTERN int trainer_mode;
SVEXTERN int trainer_active;
SVEXTERN int SVDECL(slice_average_flag,0);
SVEXTERN int show_slice_average,vis_slice_average;
SVEXTERN float slice_average_interval;

SVEXTERN int maxtourframes;
SVEXTERN int blockageSelect;
SVEXTERN int SVDECL(ntourknots,0);
SVEXTERN int SVDECL(itourknots,-1);
SVEXTERN int stretch_var_black, stretch_var_white, move_var;

SVEXTERN int SVDECL(research_mode,0);
SVEXTERN int SVDECL(showhide_option,SHOWALL_FILES);
SVEXTERN int snifferrornumber;
SVEXTERN int xyz_dir;
SVEXTERN int which_face;
SVEXTERN int showfontmenu;

SVEXTERN float SVDECL(vecfactor,1.0),SVDECL(veclength,0.0);

SVEXTERN int glui_active;

SVEXTERN int drawColorLabel,olddrawColorLabel;
SVEXTERN int vis3DSmoke3D;
SVEXTERN int smokeskip,smokeskipm1;
SVEXTERN int nrooms,nzoneinfo, nfires;
SVEXTERN float SVDECL(window_aspect_ratio,1.0), SVDECL(scene_aspect_ratio,1.0);
SVEXTERN int UpdateLIGHTS;

SVEXTERN int SVDECL(screenWidth,640), SVDECL(screenHeight,480);
SVEXTERN int SVDECL(screenWidthINI,640), SVDECL(screenHeightINI,480);
SVEXTERN int SVDECL(renderW,640), SVDECL(renderH,480), render_option;
SVEXTERN int SVDECL(glui_screenWidth,640), SVDECL(glui_screenHeight,480);
SVEXTERN int windowsize_pointer;
SVEXTERN int SVDECL(zonecolortype, ZONETEMP_COLOR);
SVEXTERN int mxframepoints;
SVEXTERN int SVDECL(timebar_drag,0),SVDECL(colorbar_drag,0),SVDECL(colorbar_splitdrag,0);
SVEXTERN int SVDECL(global_colorbar_index,-1);
SVEXTERN int fontindex;

SVEXTERN int SVDECL(custom_worldcenter,0),SVDECL(show_rotation_center,0);
SVEXTERN float SVDECL(xcenGLOBAL,0.5), SVDECL(ycenGLOBAL,0.5), SVDECL(zcenGLOBAL,0.5);
SVEXTERN float SVDECL(xcenCUSTOM,0.5), SVDECL(ycenCUSTOM,0.5), SVDECL(zcenCUSTOM,0.5);
SVEXTERN float SVDECL(xcenCUSTOMsmv,0.5), SVDECL(ycenCUSTOMsmv,0.5), SVDECL(zcenCUSTOMsmv,0.5);
SVEXTERN int glui_rotation_index,SVDECL(update_rotation_center,0);
SVEXTERN int glui_rotation_index_ini,SVDECL(update_rotation_center_ini,0);

SVEXTERN float xbar, ybar, zbar;
SVEXTERN float xbar0, ybar0, zbar0;
SVEXTERN float xbarORIG, ybarORIG, zbarORIG;
SVEXTERN float xbar0ORIG, ybar0ORIG, zbar0ORIG;
SVEXTERN int ReadPlot3dFile, ReadIsoFile;
SVEXTERN int ReadVolSlice;
SVEXTERN int Read3DSmoke3DFile;
SVEXTERN int ReadZoneFile, ReadPartFile, ReadEvacFile;

SVEXTERN int SVDECL(cache_qdata,1);

SVEXTERN int editwindow_status;
SVEXTERN int startup_pass;

SVEXTERN int slicefilenumber;
SVEXTERN int exportdata;
SVEXTERN int SVDECL(frame_count,1), SVDECL(last_frame_count,1);
SVEXTERN int nspr;
SVEXTERN int SVDECL(RenderSkip,1);
SVEXTERN int SVDECL(isoframestep_global,1),SVDECL(isoframeskip_global,0);
SVEXTERN int smoke3dframestep;
SVEXTERN int smoke3dframeskip;
SVEXTERN int vectorskip;
SVEXTERN int SVDECL(frame_index,0), SVDECL(first_frame_index,0), SVDECL(izone,0);
SVEXTERN int rotation_type,eyeview_level;
SVEXTERN int rotation_type_old,eyeview_SAVE,eyeview_last;
SVEXTERN int frameratevalue;
SVEXTERN int setpartmin, setpartmax, SVDECL(endian_smv,0);
SVEXTERN int SVDECL(setslicemin,PERCENTILE_MIN), SVDECL(setslicemax,PERCENTILE_MAX);
SVEXTERN int SVDECL(setslicemin_save,PERCENTILE_MIN), SVDECL(setslicemax_save,PERCENTILE_MAX);
SVEXTERN int SVDECL(setpatchmin_save, PERCENTILE_MIN), SVDECL(setpatchmax_save, PERCENTILE_MAX);
SVEXTERN int SVDECL(setpartmin_save, PERCENTILE_MIN), SVDECL(setpartmax_save, PERCENTILE_MAX);

SVEXTERN float slice_line_contour_min;
SVEXTERN float slice_line_contour_max;
SVEXTERN int slice_line_contour_num;
SVEXTERN int setpartmin_old, setpartmax_old;
SVEXTERN int setpatchmin, setpatchmax, SVDECL(setzonemin,GLOBAL_MIN), SVDECL(setzonemax,GLOBAL_MAX);
SVEXTERN int SVDECL(loadpatchbysteps,UNCOMPRESSED );
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
SVEXTERN float tmin_p, tmin_b, tmin_s, tmin_z, tmin_i;
SVEXTERN float tmax_p, tmax_b, tmax_s, tmax_z, tmax_i;
SVEXTERN float patchmin, patchmax;
SVEXTERN float targetmin, targetmax;
SVEXTERN float partmin, partmax;
SVEXTERN float slicemin, slicemax;
SVEXTERN float slicemin_save, slicemax_save;
SVEXTERN float patchmin_save, patchmax_save;
SVEXTERN float partmin_save, partmax_save;

SVEXTERN float SVDECL(zonemin,1.0), SVDECL(zonemax,0.0);
SVEXTERN float speedmax;
SVEXTERN int SVDECL(axislabels_smooth,1),SVDECL(axislabels_smooth_save,1);
SVEXTERN propdata SVDECL(*prop_evacdefault,NULL);
SVEXTERN float hrrpuv_max_smv;
SVEXTERN int FlowDir,ClipDir;
SVEXTERN int plotn;
SVEXTERN int stept;
SVEXTERN int plotstate;
SVEXTERN int visVector;
SVEXTERN int visSmokePart, visSprinkPart, havesprinkpart;
SVEXTERN int visaxislabels;
SVEXTERN int numplot3dvars;
SVEXTERN int SVDECL(skip_global,1);
SVEXTERN int p3dsurfacesmooth;
SVEXTERN int SVDECL(p3dsurfacetype,SURFACE_SOLID);

SVEXTERN int parttype;
SVEXTERN int allexterior,showexterior;
SVEXTERN int allinterior;
SVEXTERN int SVDECL(showedit_dialog,0);
SVEXTERN int SVDECL(showcolorbar_dialog,0);
SVEXTERN int SVDECL(showtour_dialog,0),SVDECL(showtrainer_dialog,0);
SVEXTERN int SVDECL(showtours,0);

SVEXTERN float shooter_xyz[3], shooter_dxyz[3], shooter_uvw[3], SVDECL(shooterpointsize,4.0);
SVEXTERN float shooter_velx, shooter_vely, shooter_velz, shooter_time, shooter_time_max;
SVEXTERN int SVDECL(shooter_cont_update,0),SVDECL(shooter_firstframe,0);
SVEXTERN float SVDECL(shooter_u0,2.0), SVDECL(shooter_z0,1.0), SVDECL(shooter_p,1.0/7.0), SVDECL(shooter_v_inf,1.0);
SVEXTERN float shooter_velmag, shooter_veldir, shooter_duration, SVDECL(shooter_history,10.0);
SVEXTERN int SVDECL(shooter_active,0);
SVEXTERN int shooter_fps,shooter_vel_type, shooter_nparts, SVDECL(visShooter,0), showshooter, nshooter_frames, max_shooter_points;
SVEXTERN shootpointdata SVDECL(*shootpointinfo,NULL);
SVEXTERN shoottimedata SVDECL(*shoottimeinfo,NULL);
SVEXTERN int SVDECL(*shooter_timeslist,NULL);
SVEXTERN int SVDECL(shooter_itime,0);

SVEXTERN int showgluitrainer;
SVEXTERN int colorbartype,colorbartype_ini,colorbartype_default;
SVEXTERN char colorbarname[1024];
SVEXTERN int SVDECL(update_colorbartype,0);
SVEXTERN int colorbartype_save;
SVEXTERN int colorbarpoint;
SVEXTERN int vectorspresent;
SVEXTERN int SVDECL(colorbar_hidescene,0);

SVEXTERN int visAIso;
SVEXTERN int surfincrement,visiso;
SVEXTERN int  isotest;
SVEXTERN int isolevelindex, isolevelindex2;
SVEXTERN float pref,pamb,tamb;
SVEXTERN int ntc_total, nspr_total, nheat_total;
SVEXTERN int n_devices;

SVEXTERN int npartinfo, nsliceinfo, nvsliceinfo, nslice2, npatch2, nplot3dinfo, npatchinfo;
SVEXTERN int nfedinfo;
SVEXTERN int nevac;
SVEXTERN int SVDECL(nsmoke3dinfo,0);
SVEXTERN int nisoinfo, niso_bounds;
SVEXTERN int ntrnx, ntrny, ntrnz,npdim,nmeshes,clip_mesh;
SVEXTERN int SVDECL(nOBST,0),SVDECL(nVENT,0),SVDECL(nCVENT,0),SVDECL(ncvents,0),noffset;
SVEXTERN int visLabels;
SVEXTERN int showallslicevectors;
SVEXTERN float framerate;
SVEXTERN int nglobal_times, SVDECL(ntimes_old,0), itimes, itime_save, itimeold, seqnum,RenderTime,RenderTimeOld;
SVEXTERN int nopart;
SVEXTERN int uindex, vindex, windex;

SVEXTERN int SVDECL(contour_type,0), SVDECL(p3cont3dsmooth,0);
SVEXTERN int cullfaces;
SVEXTERN int showonly_hiddenfaces;

SVEXTERN int windowresized;

SVEXTERN int SVDECL(updatemenu,0), first_display;
SVEXTERN int updatezoommenu,SVDECL(updatezoomini,0);
SVEXTERN int updatemenu_count;
SVEXTERN int SVDECL(use_graphics,1);

SVEXTERN int updatefaces,updatefacelists;
SVEXTERN int updateOpenSMVFile;

SVEXTERN int periodic_reloads;
SVEXTERN int periodic_value;

SVEXTERN int slicefilenum;
SVEXTERN int partfilenum,zonefilenum;
SVEXTERN int targfilenum;

SVEXTERN float min_gridcell_size;

SVEXTERN volfacelistdata SVDECL(*volfacelistinfo,NULL),SVDECL(**volfacelistinfoptrs,NULL);
SVEXTERN int SVDECL(nvolfacelistinfo,0);
SVEXTERN int SVDECL(update_makeiblank_smoke3d,0), SVDECL(update_initcull,0);
SVEXTERN int setPDIM;
SVEXTERN int menustatus;
SVEXTERN int SVDECL(visTimeZone,1), SVDECL(visTimeParticles,1), SVDECL(visTimeSlice,1), SVDECL(visTimePatch,1);
SVEXTERN int SVDECL(visTimeIso,1), SVDECL(visTimeEvac,1);
SVEXTERN int SVDECL(vishmsTimelabel,0), SVDECL(visTimebar,1);
SVEXTERN int SVDECL(visColorbar,1), SVDECL(visColorbar_save,1);
SVEXTERN int SVDECL(visTitle,1), SVDECL(visFullTitle,1), SVDECL(visFramerate,0);
SVEXTERN int SVDECL(visFramelabel,1), SVDECL(visTimelabel,1);
SVEXTERN int SVDECL(visHRRlabel,0);
#ifdef pp_memstatus
SVEXTERN int visAvailmemory;
#endif
SVEXTERN int SVDECL(block_volsmoke,1),SVDECL(smoke3dVoldebug,0);
SVEXTERN slicedata SVDECL(*sd_shown,NULL);
SVEXTERN vslicedata SVDECL(*vd_shown,NULL);
SVEXTERN int SVDECL(show_all_slices,1);
SVEXTERN int SVDECL(autoterrain,0),SVDECL(manual_terrain,0);
SVEXTERN float zterrain_max, zterrain_min;
SVEXTERN char SVDECL(*fds_version, NULL), SVDECL(*fds_githash, NULL);
SVEXTERN char smv_githash[256], smv_gitdate[256];
SVEXTERN int SVDECL(visMeshlabel, 1);
SVEXTERN int SVDECL(visOpenVents,1),SVDECL(visDummyVents,1),SVDECL(visOtherVents,1),SVDECL(visOtherVentsSAVE,1),SVDECL(visCircularVents,VENT_CIRCLE);
SVEXTERN int SVDECL(visOpenVentsAsOutline,0);
SVEXTERN int SVDECL(visParticles,1), SVDECL(visZone,0);
SVEXTERN int SVDECL(visEvac,1);
SVEXTERN int visBlocks;
SVEXTERN int SVDECL(outline_color_flag,0);
SVEXTERN int SVDECL(solid_state,-1),SVDECL(outline_state,-1);
SVEXTERN int visTransparentBlockage;
SVEXTERN int visBlocksSave;
SVEXTERN int SVDECL(blocklocation,BLOCKlocation_grid);
SVEXTERN int ncadgeom;
SVEXTERN int visFloor, visFrame;
SVEXTERN int visNormalEditColors;
SVEXTERN int visWalls, visGrid, visCeiling, cursorPlot3D;
SVEXTERN int SVDECL(visVZone,1), SVDECL(visHZone,0), SVDECL(viszonefire,1), SVDECL(visSZone,0);
SVEXTERN int visSensor, visSensorNorm, hasSensorNorm;
SVEXTERN int SVDECL(visVents, 1), SVDECL(visVentFlow, 1),SVDECL(visVentHFlow, 1),SVDECL(visVentVFlow, 1),SVDECL(visVentMFlow, 1);
SVEXTERN int partframestep, sliceframestep, boundframestep;
SVEXTERN int partframeskip, sliceframeskip, boundframeskip;
SVEXTERN int boundzipstep, boundzipskip;
SVEXTERN int smoke3dzipstep, smoke3dzipskip;
SVEXTERN int slicezipstep, slicezipskip;
SVEXTERN int isozipstep, isozipskip;
SVEXTERN int evacframeskip, evacframestep;
SVEXTERN int viewoption;
SVEXTERN int SVDECL(clip_mode,CLIP_OFF),clip_mode_last;
SVEXTERN int clip_i,clip_j,clip_k;
SVEXTERN int clip_I,clip_J,clip_K;
SVEXTERN clipdata clipinfo,colorbar_clipinfo;
SVEXTERN int stepclip_xmin,stepclip_ymin,stepclip_zmin;
SVEXTERN int stepclip_xmax,stepclip_ymax,stepclip_zmax;
SVEXTERN float partpointsize,SVDECL(vectorpointsize,2.0),streaklinewidth;
SVEXTERN float isopointsize, isolinewidth;
SVEXTERN float plot3dpointsize, plot3dlinewidth;
SVEXTERN int SVDECL(scaled_font3d_thickness,1);
SVEXTERN int SVDECL(scaled_font2d_thickness,1);
SVEXTERN float SVDECL(vectorlinewidth,1.0);
SVEXTERN int SVDECL(cell_center_text,0);
SVEXTERN float SVDECL(gridlinewidth,2.0),SVDECL(ticklinewidth,2.0);
SVEXTERN int SVDECL(zone_highlight,0),SVDECL(zone_highlight_room,0);
SVEXTERN int SVDECL(script_step,0), SVDECL(script_step_now,0);
SVEXTERN int SVDECL(script_keystate,0);
SVEXTERN int SVDECL(render_clip_left,0);
SVEXTERN int SVDECL(render_clip_right,0);
SVEXTERN int SVDECL(render_clip_bottom,0);
SVEXTERN int SVDECL(render_clip_top,0);
SVEXTERN int SVDECL(clip_rendered_scene,0);

SVEXTERN float sprinklerabssize, sensorabssize, heatabssize;
SVEXTERN float SVDECL(sensorrelsize,1.0),SVDECL(sensorrelsizeMIN,0.0);
SVEXTERN float SVDECL(vector_baselength,1.0);
SVEXTERN float SVDECL(vector_basediameter,0.1);
SVEXTERN float SVDECL(vector_headlength,0.2);
SVEXTERN float SVDECL(vector_headdiameter,0.2);

SVEXTERN float linewidth, ventlinewidth, highlight_linewidth,solidlinewidth;
SVEXTERN float SVDECL(sliceoffset_factor,0.1), SVDECL(ventoffset_factor,0.01);
SVEXTERN int visBLOCKold;

SVEXTERN int selectedcolorbar_index,selectedcolorbar_index2;
SVEXTERN int planar_terrain_slice;
SVEXTERN int nrgb;
SVEXTERN int nrgb_ini;
SVEXTERN int nrgb2_ini;
SVEXTERN int rgb_white, rgb_yellow, rgb_blue, rgb_red;
SVEXTERN int rgb_green, rgb_magenta, rgb_cyan, rgb_black;
SVEXTERN int numColorbars;
SVEXTERN int setbw,SVDECL(setbwdata,0);
SVEXTERN int setbwSAVE;
SVEXTERN int background_flip;
SVEXTERN float SVDECL(transparent_level,0.8);
SVEXTERN int SVDECL(use_transparency_data,1);
SVEXTERN int antialiasflag;
SVEXTERN int nrgb_full;
SVEXTERN int nrgb_cad;
SVEXTERN float eyexfactor, eyeyfactor, eyezfactor;
SVEXTERN int SVDECL(transparent_state,ALL_SOLID);
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

SVEXTERN int mainwindow_id;
SVEXTERN int rendertourcount;

#ifdef pp_MEMDEBUG
SVEXTERN int list_memcheck_index;
SVEXTERN int SVDECL(visUsagememory,0);
#endif
SVEXTERN float gslice_norm[3];
#ifdef INMAIN
SVEXTERN float tour_xyz[3]={0.0,0.0,0.0};
SVEXTERN float gslice_xyz[3]={-1000001.0,-1000001.0,-1000001.0};
SVEXTERN float gslice_normal_xyz[3]={0.0,0.0,1.0};
SVEXTERN float gslice_normal_azelev[2]={0.0,90.0};
#else
SVEXTERN float tour_xyz[3];
SVEXTERN float gslice_xyz[3];
SVEXTERN float gslice_normal_xyz[3];
SVEXTERN float gslice_normal_azelev[3];
#endif

SVEXTERN float gslice_xyz0[3],gslice_normal_azelev0[2];
SVEXTERN int SVDECL(vis_gslice_data,0),SVDECL(SHOW_gslice_data,0),SVDECL(SHOW_gslice_data_old,0),SVDECL(show_gslice_triangles,0);
SVEXTERN int SVDECL(show_gslice_triangulation,0);
SVEXTERN int SVDECL(show_gslice_normal,0),SVDECL(show_gslice_normal_keyboard,0);


SVEXTERN float ventcolor_orig[4];
SVEXTERN float SVDECL(*ventcolor,NULL);
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

SVEXTERN char SVDECL(*smokeviewtempdir,NULL);

SVEXTERN int SVDECL(showall_boundary,1);
SVEXTERN int nmenus;
SVEXTERN menudata menuinfo[10000];
SVEXTERN int showbuild;
SVEXTERN int max_screenWidth, max_screenHeight;
SVEXTERN int saveW, saveH;
SVEXTERN char SVDECL(*texturedir,NULL);
SVEXTERN char release_title[1024];
SVEXTERN char plot3d_title[1024];
SVEXTERN char SVDECL(*partshortlabel,NULL),SVDECL(*partunitlabel,NULL);
SVEXTERN char emptylabel[2];
SVEXTERN void SVDECL(*large_font,NULL);
SVEXTERN void SVDECL(*small_font,NULL);

SVEXTERN int nopenvents,nopenvents_nonoutline,ndummyvents,ntransparentblocks,ntransparentvents;
SVEXTERN int nventcolors;
SVEXTERN float SVDECL(**ventcolors,NULL);
SVEXTERN float texture_origin[3];

SVEXTERN int vslicecolorbarflag;
SVEXTERN int SVDECL(use_new_drawface,0);
#ifdef INMAIN
  SVEXTERN unsigned char rgb_below_min[3]={255-64,255-64,255-64}, rgb_above_max[3]={0,0,0};
#else
  SVEXTERN unsigned char rgb_below_min[3], rgb_above_max[3];
#endif
SVEXTERN int SVDECL(colorbar_select_index,-1),SVDECL(update_colorbar_select_index,0);
SVEXTERN float world_eyepos[3],scaled_eyepos[3];
SVEXTERN int tour_usecurrent;
SVEXTERN int isZoneFireModel;
SVEXTERN int SVDECL(output_slicedata,0),SVDECL(output_patchdata,0);
SVEXTERN f_units SVDECL(*unitclasses,NULL),SVDECL(*unitclasses_default,NULL),SVDECL(*unitclasses_ini,NULL);
SVEXTERN int nunitclasses,nunitclasses_default,nunitclasses_ini;
SVEXTERN meshdata SVDECL(*meshinfo,NULL),SVDECL(*current_mesh,NULL), SVDECL(*mesh_save,NULL);
SVEXTERN supermeshdata SVDECL(*supermeshinfo,NULL);
SVEXTERN int SVDECL(nsupermeshinfo,0);
SVEXTERN meshdata SVDECL(*mesh_last,NULL), SVDECL(*loaded_isomesh,NULL);
SVEXTERN float devicenorm_length;
SVEXTERN int SVDECL(ndeviceinfo,0),nvdeviceinfo,ndeviceinfo_exp;
SVEXTERN float max_dev_vel;
SVEXTERN int SVDECL(last_prop_display,-1);
SVEXTERN int SVDECL(devicetypes_index,0);
SVEXTERN devicedata SVDECL(*deviceinfo,NULL);
SVEXTERN vdevicedata SVDECL(*vdeviceinfo, NULL);
SVEXTERN vdevicesortdata SVDECL(*vdevices_sorted, NULL);
SVEXTERN int SVDECL(ntreedeviceinfo, 0), SVDECL(mintreesize, 3);
SVEXTERN treedevicedata SVDECL(*treedeviceinfo,NULL);
SVEXTERN int SVDECL(show_smokesensors,SMOKESENSORS_0255),active_smokesensors,test_smokesensors;
SVEXTERN float smoke3d_cvis;
SVEXTERN sv_object SVDECL(**object_defs,NULL), SVDECL(*heat_detector_object_backup,NULL), SVDECL(*target_object_backup,NULL);
SVEXTERN sv_object SVDECL(*sprinkler_upright_object_backup,NULL), SVDECL(*smoke_detector_object_backup,NULL);
SVEXTERN sv_object SVDECL(*thcp_object_backup,NULL), SVDECL(*missing_device,NULL), SVDECL(*error_device,NULL);
SVEXTERN sv_object object_def_first, object_def_last;
SVEXTERN char SVDECL(**device_texture_list,NULL);
SVEXTERN int ndevice_texture_list, SVDECL(*device_texture_list_index,NULL);
SVEXTERN int SVDECL(nobject_defs,0);
SVEXTERN int svofile_exists;
SVEXTERN treedata SVDECL(*treeinfo,NULL);
SVEXTERN terraindata SVDECL(*terraininfo,NULL);
SVEXTERN int SVDECL(ntreeinfo,0), SVDECL(nterraininfo,0), SVDECL(visTerrainType,0);
SVEXTERN float treecolor[4], treecharcolor[4], trunccolor[4];
SVEXTERN int showterrain;
SVEXTERN float rgb_terrain[10][4];
SVEXTERN tourdata SVDECL(*tourinfo,NULL);
SVEXTERN keyframe SVDECL(**tourknotskeylist,NULL);
SVEXTERN tourdata SVDECL(**tourknotstourlist,NULL);
SVEXTERN keyframe SVDECL(*selected_frame,NULL);
SVEXTERN tourdata SVDECL(*selected_tour,NULL);
SVEXTERN int callfrom_tourglui;
SVEXTERN int showtours_whenediting;

SVEXTERN int SVDECL(*slice_loadstack,NULL),  SVDECL(nslice_loadstack,0),  SVDECL(islice_loadstack,0);
SVEXTERN int SVDECL(*mslice_loadstack,NULL), SVDECL(nmslice_loadstack,0), SVDECL(imslice_loadstack,0);
SVEXTERN int SVDECL(*vslice_loadstack,NULL), SVDECL(nvslice_loadstack,0), SVDECL(ivslice_loadstack,0);
SVEXTERN int SVDECL(*mvslice_loadstack,NULL),SVDECL(nmvslice_loadstack,0),SVDECL(imvslice_loadstack,0);
SVEXTERN int SVDECL(*subslice_menuindex,NULL),SVDECL(*subvslice_menuindex,NULL);

SVEXTERN float xtimeleft, xtimeright;

SVEXTERN int SVDECL(stereoactive,0);
SVEXTERN int SVDECL(stereotype,STEREO_NONE), SVDECL(stereotypeOLD, STEREO_NONE);
SVEXTERN int SVDECL(show_parallax,0), SVDECL(stereotype_frame, BOTH_EYES);

SVEXTERN int SVDECL(show_hrrcutoff,1), SVDECL(show_hrrcutoff_active,0),SVDECL(hrrpuv_loaded,0);
SVEXTERN int trainerview;
SVEXTERN int apertureindex;
SVEXTERN int zoomindex;
SVEXTERN int projection_type;
SVEXTERN float apertures[5];
SVEXTERN float aperture,aperture_glui,aperture_default;
SVEXTERN float zooms[5];
SVEXTERN float zoom,SVDECL(zoomini,1.0);
SVEXTERN int rgbmask[16];
SVEXTERN GLint nredbits, ngreenbits, nbluebits;
SVEXTERN int nredshift, ngreenshift, nblueshift;
SVEXTERN float xyzbox;
SVEXTERN float xplts[256*256], yplts[256*256], zplts[256*256];
SVEXTERN float SVDECL(*targtimes,NULL);
SVEXTERN int SVDECL(*targtimeslist,NULL);
SVEXTERN int SVDECL(*zone_timeslist,NULL);
SVEXTERN int delete_view_is_disabled;
SVEXTERN int old_listview;

SVEXTERN int sv_age;
SVEXTERN int titlesafe_offset;
SVEXTERN int titlesafe_offsetBASE;
SVEXTERN int   reset_frame;
SVEXTERN float reset_time,start_frametime,stop_frametime;
SVEXTERN int reset_time_flag;
SVEXTERN float SVDECL(velocity_range,0.0);
SVEXTERN int niso_compressed;
SVEXTERN int nslice_loaded, npatch_loaded;
SVEXTERN int SVDECL(*slice_loaded_list,NULL), SVDECL(*patch_loaded_list,NULL);
SVEXTERN int SVDECL(*render_frame,NULL);
SVEXTERN int RenderOnceNow;
SVEXTERN char SVDECL(*fdsprefix,NULL), SVDECL(*fdsprefix2,NULL);
SVEXTERN char SVDECL(*endian_filename,NULL);
SVEXTERN char SVDECL(*target_filename,NULL);

SVEXTERN int SVDECL(update_bounds,0);
SVEXTERN int SVDECL(*sorted_surfidlist,NULL),SVDECL(*inv_sorted_surfidlist,NULL),nsorted_surfidlist;
SVEXTERN char SVDECL(*trainer_filename,NULL), SVDECL(*test_filename,NULL);
SVEXTERN FILE SVDECL(*STREAM_SB,NULL);
SVEXTERN time_t smv_modtime;
SVEXTERN float temp_threshold;
SVEXTERN char SVDECL(*smv_filename,NULL),SVDECL(*fed_filename,NULL),fed_filename_base[1024],SVDECL(*stop_filename,NULL);
SVEXTERN char SVDECL(*sliceinfo_filename,NULL);
SVEXTERN char SVDECL(*database_filename,NULL),SVDECL(*smokeview_bindir,NULL),SVDECL(*iso_filename,NULL);
SVEXTERN scriptfiledata first_scriptfile, last_scriptfile, SVDECL(*default_script,NULL);
#ifdef pp_LUA
SVEXTERN luascriptfiledata first_luascriptfile, last_luascriptfile, SVDECL(*default_luascript,NULL);
SVEXTERN int SVDECL(luascript_loaded,0);
#endif
SVEXTERN scriptdata SVDECL(*scriptinfo,NULL), SVDECL(*current_script_command,NULL);
SVEXTERN char SVDECL(*script_dir_path,NULL);
SVEXTERN int SVDECL(nscriptinfo,0);
SVEXTERN scriptfiledata SVDECL(*script_recording,NULL);
SVEXTERN int SVDECL(runscript,0), SVDECL(noexit,0);
#ifdef pp_LUA
SVEXTERN int SVDECL(runluascript,0);
SVEXTERN int SVDECL(exit_on_script_crash,0);
#endif
SVEXTERN int SVDECL(script_multislice,0), SVDECL(script_multivslice,0), SVDECL(script_iso,0);
SVEXTERN FILE SVDECL(*scriptoutstream,NULL);
SVEXTERN char SVDECL(*log_filename,NULL);
SVEXTERN FILE SVDECL(*LOG_FILENAME,NULL);
SVEXTERN char SVDECL(*flushfile,NULL), SVDECL(*chidfilebase,NULL);
SVEXTERN char SVDECL(*hrr_csv_filename,NULL),SVDECL(*devc_csv_filename,NULL),SVDECL(*exp_csv_filename,NULL);
SVEXTERN hrrdata SVDECL(*hrrinfo,NULL);
SVEXTERN char SVDECL(*smokezippath,NULL),SVDECL(*smokeviewpath,NULL);
SVEXTERN char SVDECL(*INI_fds_filein,NULL), SVDECL(*fds_filein,NULL);
SVEXTERN char SVDECL(*caseini_filename,NULL),SVDECL(*boundini_filename,NULL);
SVEXTERN char SVDECL(*zonelonglabels,NULL), SVDECL(*zoneshortlabels,NULL), SVDECL(*zoneunits,NULL);
SVEXTERN char SVDECL(*smokeviewini,NULL);
SVEXTERN int overwrite_all,erase_all;
SVEXTERN int compress_autoloaded;
SVEXTERN tridata SVDECL(**opaque_triangles,NULL),SVDECL(**transparent_triangles,NULL),SVDECL(**alltriangles,NULL);
SVEXTERN int SVDECL(nopaque_triangles,0),SVDECL(ntransparent_triangles,0),SVDECL(nalltriangles,0);
#ifdef WIN32
SVEXTERN   char openfilebuffer[1024];
SVEXTERN   int openfileflag;
#endif
SVEXTERN float xyzmaxdiff;
SVEXTERN char ext_png[5];
SVEXTERN char ext_jpg[5];
SVEXTERN int render_filetype;
SVEXTERN int SVDECL(renderfilelabel,0);
SVEXTERN char part_ext[6];
SVEXTERN char ini_ext[5];

SVEXTERN int SVDECL(updatehiddenfaces,1),SVDECL(hide_overlaps,0);
SVEXTERN int SVDECL(nsurfids,0);
SVEXTERN surfid SVDECL(*surfids,NULL);
SVEXTERN int key_state;
SVEXTERN float starteyex, starteyey;
SVEXTERN float eye_xyz0[3];
SVEXTERN float start_xyz0[3];
SVEXTERN int glui_move_mode;

SVEXTERN float timeoffset;
SVEXTERN int npartpoints, npartframes;
SVEXTERN float xslicemid, yslicemid, zslicemid;
SVEXTERN float delx;
SVEXTERN float delz;
SVEXTERN float d_eye_xyz[3],dsave_eye_xyz[3];
SVEXTERN float eyex0, eyey0, eyez0;
SVEXTERN float viewx, viewy, viewz;
SVEXTERN float anglexy0,azimuth0;
SVEXTERN int mouse_down_xy0[2];
SVEXTERN int touring;
SVEXTERN int update_tourlist;
SVEXTERN float desired_view_height;
SVEXTERN int thistime, lasttime, resetclock,initialtime;
SVEXTERN int realtime_flag;
SVEXTERN char timelabel[30];
SVEXTERN char frameratelabel[30];
SVEXTERN char framelabel[30];
SVEXTERN float SVDECL(**p3levels,NULL), SVDECL(*zonelevels,NULL);
SVEXTERN float SVDECL(**p3levels256,NULL);
SVEXTERN char SVDECL(***colorlabelp3,NULL),SVDECL(***colorlabeliso,NULL);
SVEXTERN char SVDECL(**scalep3,NULL), SVDECL(**scalep3copy,NULL);
SVEXTERN float SVDECL(*fscalep3,NULL);
SVEXTERN char SVDECL(*partscale,NULL);
SVEXTERN char a_partscale[31];
SVEXTERN char SVDECL(*zonescale,NULL);
SVEXTERN char a_zonescale[31];
SVEXTERN int islicetype,islicetype_save,ipatchtype;
SVEXTERN int iisotype,iisottype;
SVEXTERN char SVDECL(**colorlabelpart,NULL), SVDECL(**colorlabelpatch,NULL),  SVDECL(**colorlabelzone,NULL);

SVEXTERN int SVDECL(hilight_skinny,0);

SVEXTERN int minfill, maxfill;

SVEXTERN int SVDECL(*plotiso,NULL);
SVEXTERN float SVDECL(*global_times,NULL),cputimes[20];
SVEXTERN int cpuframe;

SVEXTERN float xyzeyeorig[3],xeyedir[3], yeyedir[3], zeyedir[3];
SVEXTERN int adjustalphaflag;
SVEXTERN int SVDECL(colorband,5);
SVEXTERN int SVDECL(have_extreme_mindata,0), SVDECL(have_extreme_maxdata,0);
SVEXTERN int SVDECL(show_extreme_mindata,0), SVDECL(show_extreme_maxdata,0);
SVEXTERN int SVDECL(show_extreme_mindata_save,0), SVDECL(show_extreme_maxdata_save,0);

SVEXTERN int SVDECL(use_iblank,1),SVDECL(iblank_set_on_commandline,0);

SVEXTERN int script_index, ini_index;
SVEXTERN char script_inifile_suffix[1024], vol_prefix[1024];
SVEXTERN char script_renderdir[1024], script_renderfilesuffix[1024], script_renderfile[1024];
SVEXTERN inifiledata first_inifile, last_inifile;
SVEXTERN char script_filename[1024];
#ifdef pp_LUA
SVEXTERN char luascript_filename[1024];
#endif
SVEXTERN int highlight_block, highlight_mesh, highlight_flag;
SVEXTERN int SVDECL(updategetobstlabels,1);

SVEXTERN int pixel_skip;
SVEXTERN float smoke_extinct,smoke_dens,smoke_pathlength;
SVEXTERN int smoke_alpha;
SVEXTERN int smoketest,show_smoketest;
SVEXTERN int showall_textures;
SVEXTERN int SVDECL(enable_texture_lighting,0);

SVEXTERN int SVDECL(ncolorbars,0);
SVEXTERN int ndefaultcolorbars;
SVEXTERN colorbardata SVDECL(*colorbarinfo,NULL),SVDECL(*current_colorbar,NULL);

SVEXTERN int SVDECL(ncolortableinfo, 0);
SVEXTERN colortabledata SVDECL(*colortableinfo, NULL);
SVEXTERN int SVDECL(i_colortable_list,0);

SVEXTERN int SVDECL(update_load_Files, 0);
SVEXTERN int do_threshold;
SVEXTERN int ntotal_blockages;
SVEXTERN int updateindexcolors;
SVEXTERN int show_path_knots;
SVEXTERN int keyframe_snap;
SVEXTERN int tourviewtype;
SVEXTERN int SVDECL(show_tourlocus,1);
SVEXTERN int tourlocus_type;
SVEXTERN int iavatar_types, navatar_types;
SVEXTERN int iavatar_evac;
SVEXTERN sv_object SVDECL(**avatar_types,NULL);
SVEXTERN int glui_avatar_index;
SVEXTERN sv_object *avatar_defs_backup[2];
SVEXTERN int SVDECL(device_sphere_segments,6);
SVEXTERN int ntexturestack;

SVEXTERN float SVDECL(fire_opacity_factor,3.0),SVDECL(mass_extinct,8700.0);
SVEXTERN float SVDECL(temperature_min,20.0),SVDECL(temperature_cutoff,700.0),SVDECL(temperature_max,1200.0);
SVEXTERN float SVDECL(global_hrrpuv_min,0.0),SVDECL(global_hrrpuv_cutoff,200.0),SVDECL(global_hrrpuv_max,1200.0);
SVEXTERN int SVDECL(volbw,0);
SVEXTERN float tourrad_avatar;
SVEXTERN int dirtycircletour;
SVEXTERN float SVDECL(*tour_t,NULL), SVDECL(*tour_t2,NULL), SVDECL(*tour_dist,NULL), SVDECL(*tour_dist2,NULL), SVDECL(*tour_dist3,NULL);
SVEXTERN float view_tstart, view_tstop;
SVEXTERN int tour_constant_vel;
SVEXTERN float tour_bias,tour_continuity;
SVEXTERN int view_ntimes;

SVEXTERN int SVDECL(ntours, 0);
SVEXTERN int SVDECL(selectedtour_index, TOURINDEX_MANUAL), SVDECL(selectedtour_index_old, TOURINDEX_MANUAL), SVDECL(selectedtour_index_ini, TOURINDEX_MANUAL);
SVEXTERN int SVDECL(update_selectedtour_index,0);
SVEXTERN int viewtourfrompath,viewalltours,viewanytours,edittour;

SVEXTERN selectdata SVDECL(*selectfaceinfo,NULL);
SVEXTERN blockagedata SVDECL(**selectblockinfo,NULL);
SVEXTERN tickdata SVDECL(*tickinfo,NULL);
SVEXTERN int SVDECL(ntickinfo,0),SVDECL(ntickinfo_smv,0);
SVEXTERN int visFDSticks;
SVEXTERN float user_tick_origin[3], user_tick_max[3], user_tick_min[3], user_tick_step[3], user_tick_length, user_tick_width;
SVEXTERN int user_tick_nxyz[3], user_tick_sub, user_tick_option, visUSERticks, auto_user_tick_placement;
SVEXTERN int user_tick_show_x, user_tick_show_y, user_tick_show_z;
SVEXTERN int visCadTextures, visTerrainTexture;
SVEXTERN int bw_colorbar_index;
SVEXTERN int SVDECL(viscolorbarpath,0);
SVEXTERN int SVDECL(*sortedblocklist,NULL),SVDECL(*changed_idlist,NULL),SVDECL(nchanged_idlist,0);
SVEXTERN int nselectblocks;
SVEXTERN surfdata SVDECL(*surfinfo,NULL),sdefault,v_surfacedefault,e_surfacedefault;
SVEXTERN int nsurfinfo;
SVEXTERN matldata SVDECL(*matlinfo,NULL);
SVEXTERN int nmatlinfo;
SVEXTERN int surface_indices[7],surface_indices_bak[7];
SVEXTERN int wall_case;
SVEXTERN surfdata SVDECL(*surfacedefault,NULL), SVDECL(*vent_surfacedefault,NULL), SVDECL(*exterior_surfacedefault,NULL);
SVEXTERN char surfacedefaultlabel[256];
SVEXTERN int ntotalfaces;
SVEXTERN colordata SVDECL(*firstcolor,NULL);
SVEXTERN texturedata SVDECL(*textureinfo,NULL), SVDECL(*terrain_texture,NULL);
SVEXTERN GLuint texture_colorbar_id, texture_slice_colorbar_id, texture_patch_colorbar_id, texture_plot3d_colorbar_id, texture_iso_colorbar_id, terrain_colorbar_id;
SVEXTERN GLuint volsmoke_colormap_id,slice3d_colormap_id,slicesmoke_colormap_id;
SVEXTERN int SVDECL(volsmoke_colormap_id_defined,-1);
SVEXTERN int SVDECL(slice3d_colormap_id_defined,-1);
SVEXTERN float mscale[3];
SVEXTERN float xclip_min, yclip_min, zclip_min;
SVEXTERN float xclip_max, yclip_max, zclip_max;
SVEXTERN float nearclip,farclip;
SVEXTERN int updateclipvals;
SVEXTERN int updateUpdateFrameRateMenu;
SVEXTERN int ntextures,ntextures_loaded_used;
SVEXTERN int SVDECL(nskyboxinfo,0);
SVEXTERN skyboxdata SVDECL(*skyboxinfo,NULL);
SVEXTERN firedata SVDECL(*fireinfo,NULL);
SVEXTERN roomdata SVDECL(*roominfo,NULL);
SVEXTERN zventdata SVDECL(*zventinfo,NULL);
SVEXTERN zonedata SVDECL(*zoneinfo,NULL);
SVEXTERN zonedata SVDECL(*activezone,NULL);
SVEXTERN partdata SVDECL(*partinfo,NULL);
SVEXTERN int SVDECL(update_screensize,0);
SVEXTERN int part5show;
SVEXTERN int streak5show,streak5value, streak5step, showstreakhead;
SVEXTERN int nstreak_rvalue; // 5
SVEXTERN float streak_rvalue[8]; // 1.0, 2.0 4.0, 8.0, 16.0 twfin
SVEXTERN int streak_index, update_streaks;       // 0
SVEXTERN float float_streak5value;// 1.0
SVEXTERN partclassdata SVDECL(*partclassinfo,NULL);
SVEXTERN int npartclassinfo;
SVEXTERN partpropdata SVDECL(*part5propinfo,NULL), SVDECL(*current_property,NULL);
SVEXTERN int SVDECL(npart5prop,0),ipart5prop,ipart5prop_old;
SVEXTERN int prop_index;
SVEXTERN slicedata SVDECL(*sliceinfo,NULL),SVDECL(*slicexyzinfo,NULL);
SVEXTERN feddata SVDECL(*fedinfo,NULL);
SVEXTERN camdata SVDECL(*caminfo,NULL);
SVEXTERN multislicedata SVDECL(*multisliceinfo,NULL);
SVEXTERN multivslicedata SVDECL(*multivsliceinfo,NULL);
SVEXTERN outlinedata SVDECL(*outlineinfo,NULL);
SVEXTERN int noutlineinfo;
SVEXTERN int nmultisliceinfo;
SVEXTERN int nmultivsliceinfo;
SVEXTERN int SVDECL(*sliceorderindex,NULL),SVDECL(*vsliceorderindex,NULL),SVDECL(*partorderindex,NULL);
SVEXTERN int SVDECL(*patchorderindex,NULL),SVDECL(*isoorderindex,NULL),SVDECL(*plot3dorderindex,NULL);
SVEXTERN int showfiles;
SVEXTERN boundsdata SVDECL(*slicebounds,NULL), SVDECL(*isobounds,NULL), SVDECL(*patchbounds,NULL);
SVEXTERN vslicedata SVDECL(*vsliceinfo,NULL);
SVEXTERN int force_redisplay;
SVEXTERN int setp3min_temp, setp3max_temp;
SVEXTERN int setp3chopmin_temp, setp3chopmax_temp;
SVEXTERN float p3chopmin_temp, p3chopmax_temp;
SVEXTERN float p3min_temp, p3max_temp;

SVEXTERN smoke3ddata SVDECL(*smoke3dinfo,NULL);
SVEXTERN int SVDECL(fire_red,255), SVDECL(fire_green,128), SVDECL(fire_blue,0);

SVEXTERN int SVDECL(smoke_red, 0), SVDECL(smoke_green, 0), SVDECL(smoke_blue, 0);
SVEXTERN float SVDECL(smoke_albedo, 0.0), SVDECL(fire_halfdepth,2.0);

SVEXTERN int SVDECL(show_firecolormap,0);
SVEXTERN int SVDECL(firecolormap_type,FIRECOLORMAP_DIRECT),SVDECL(firecolormap_type_save,FIRECOLORMAP_DIRECT);
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
SVEXTERN labeldata label_first, label_last, *label_first_ptr, *label_last_ptr;
SVEXTERN int SVDECL(*slicetypes,NULL), SVDECL(*isotypes,NULL), SVDECL(*vslicetypes,NULL), SVDECL(*patchtypes,NULL);
SVEXTERN plot3ddata SVDECL(*plot3dinfo,NULL);
SVEXTERN float SVDECL(*plot3dtimelist,NULL);
SVEXTERN int nplot3dtimelist;
SVEXTERN patchdata SVDECL(*patchinfo,NULL);
SVEXTERN isodata SVDECL(*isoinfo,NULL);

SVEXTERN blockagedata SVDECL(*bchighlight,NULL),SVDECL(*bchighlight_old,NULL);
SVEXTERN cadgeomdata SVDECL(*cadgeominfo,NULL);

SVEXTERN int smokediff;
SVEXTERN int render_size_index;
SVEXTERN int render_skip_index;
SVEXTERN int buffertype;
SVEXTERN int opengldefined;
SVEXTERN int SVDECL(restart_time,0);
SVEXTERN int nslicetypes;
SVEXTERN int nvslicetypes;
SVEXTERN int nisotypes;
SVEXTERN int SVDECL(*isosubmenus,NULL), nisosubmenus;
SVEXTERN int SVDECL(*loadpatchsubmenus,NULL), nloadpatchsubmenus;
SVEXTERN int npatchtypes;
SVEXTERN char SVDECL(**patchlabellist,NULL);
SVEXTERN int SVDECL(*patchlabellist_index,NULL);
SVEXTERN int SVDECL(*isoindex,NULL);

SVEXTERN int have_vents_int;
SVEXTERN int nface_outlines, nface_textures, nface_transparent;
SVEXTERN int nface_normals_single, nface_normals_double, nface_transparent_double, nvent_transparent;
SVEXTERN int show_transparent_vents;
SVEXTERN int SVDECL(show_bothsides_blockages, 0);
SVEXTERN int show_bothsides_int, show_bothsides_ext;
SVEXTERN float SVDECL(transparency_geom,0.2);
SVEXTERN int SVDECL(use_transparency_geom,0);
SVEXTERN facedata SVDECL(**face_transparent,NULL);
SVEXTERN int SVDECL(hidepatchsurface,0);
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

  SVEXTERN float rgbhazard[MAXRGB][4]={
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
  SVEXTERN float rgbhazard[MAXRGB][4];
#endif
#endif
