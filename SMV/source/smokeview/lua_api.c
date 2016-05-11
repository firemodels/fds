
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include "options.h"
#include "smokeviewvars.h"
#include "c_api.h"
#include "lua_api.h"

#include GLUT_H

lua_State* L;
int lua_displayCB(lua_State *L);

#ifdef WIN32
#define snprintf _snprintf
#endif

/* ------------------ load_script ------------------------ */
// There are two options for scripting, Lua and SSF. Which is run is set here
// based on the commandline arguments. If either (exclusive) of these values
// are set to true, then that script will run from within the display callback
// (Display_CB, in callbacks.c). These two loading routines are included to
// load the scripts early in the piece, before the display callback.
// Both runluascript and runscript are global.
int load_script(char *filename) {
  fprintf(stderr, "load_script: %s\n", filename);
  if (runluascript == 1 && runscript == 1) {
    fprintf(stderr, "Both a Lua script and an SSF script cannot be run "
                    "simultaneously\n");
    exit(1);
  }
  if (runluascript == 1) {
    // Load the Lua script in order for it to be run later.
    if (loadLuaScript(filename) != LUA_OK) {
      fprintf(stderr, "There was an error loading the script, and so it "
                      "will not run.\n");
      if (exit_on_script_crash) {
          exit(1); // exit with an error code
      }
      runluascript = 0; // set this to false so that the smokeview no longer
                       // tries to run the script as it failed to load
      fprintf(stderr, "Running smokeview normally.\n");
    } else {
      fprintf(stderr, "%s successfully loaded\n", filename);
    }
  }
#ifdef pp_LUA_SSF
  if (runscript == 1) {
    // Load the ssf script in order for it to be run later
    // This still uses the Lua interpreter
    if (loadSSFScript(filename) != LUA_OK) {
      fprintf(stderr, "There was an error loading the script, and so it "
                      "will not run.\n");
      if (exit_on_script_crash) {
        exit(1); // exit with an error code
      }
      runluascript = 0; // set this to false so that the smokeview no longer
                       // tries to run the script as it failed to load
      fprintf(stderr, "Running smokeview normally.\n");
    }
  }
#endif
    return 1;
}

/*
  Load a .smv file. This is currently not used as it is dependent on Smokeview
  being able to run without a .smv file loaded.
*/
int lua_loadsmvall(lua_State *L) {
  // The first argument is taken from the stack as a string.
  const char *filepath = lua_tostring(L, 1);
  printf("lua_loadsmvall filepath: %s\n",filepath);
  // The function from the C api is called using this string.
  loadsmvall(filepath);
  // 0 arguments are returned.
  return 0;
}

/*
  Set render clipping.
*/
int lua_renderclip(lua_State *L) {
  int flag = lua_toboolean(L, 1);
  int left = lua_tonumber(L, 2);
  int right = lua_tonumber(L, 3);
  int bottom = lua_tonumber(L, 4);
  int top = lua_tonumber(L, 5);
  return 0;
}

/*
  Render the current frame to a file.
*/
int lua_render(lua_State *L) {
  lua_displayCB(L);
  printf("performing lua render\n");
  printf("rendering to: %s\n", script_dir_path);
  const char *basename = lua_tostring(L, 1);
  printf("basename(lua): %s\n", basename);
  int ret = render(basename);
  lua_pushnumber(L, ret);
  return 1;
}

int lua_gsliceview(lua_State *L) {
  int data = lua_tonumber(L, 1);
  int show_triangles = lua_toboolean(L, 2);
  int show_triangulation = lua_toboolean(L, 3);
  int show_normal = lua_toboolean(L, 4);
  gsliceview(data, show_triangles, show_triangulation, show_normal);
  return 0;
}

int lua_gslicepos(lua_State *L) {
  float x = lua_tonumber(L, 1);
  float y = lua_tonumber(L, 2);
  float z = lua_tonumber(L, 3);
  gslicepos(x, y, z);
  return 0;
}
int lua_gsliceorien(lua_State *L) {
  float az = lua_tonumber(L, 1);
  float elev = lua_tonumber(L, 2);
  gsliceorien(az, elev);
  return 0;
}

int lua_settourview(lua_State *L) {
  int edittourArg = lua_tonumber(L, 1);
  int mode = lua_tonumber(L, 2);
  int show_tourlocusArg = lua_toboolean(L, 3);
  float tour_global_tensionArg = lua_tonumber(L, 4);
  settourview(edittourArg, mode, show_tourlocusArg, tour_global_tensionArg);
  return 0;
}

int lua_settourkeyframe(lua_State *L) {
  float keyframe_time = lua_tonumber(L, 1);
  settourkeyframe(keyframe_time);
  return 0;
}

/*
  Trigger the display callback.
*/
int lua_displayCB(lua_State *L) {
  // runluascript=0;
  Display_CB();
  // runluascript=1;
  return 0;
}

/*
  Hide the smokeview window. This should not currently be used as it prevents
  the display callback being called, and therefore the script will not
  continue (the script is called as part of the display callback).
*/
int lua_hidewindow(lua_State *L) {
  printf("hiding window\n");
  glutHideWindow();
  //once we hide the window the display callback is never called
  return 0;
}

/*
  By calling yieldscript, the script is suspended and the smokeview display is
  updated. It is necessary to call this before producing any outputs (such as
  renderings).
*/
int lua_yieldscript(lua_State *L) {
  printf("yielding\n");
  lua_yield(L, 0 /*zero results*/);
  return 0;
}

/*
  As with lua_yieldscript, but immediately resumes the script after letting the
  display callback run.
*/
int lua_tempyieldscript(lua_State *L) {
  printf("tempyielding\n");
  runluascript=1;
  lua_yield(L, 0 /*zero results*/);
  return 0;
}

/*
  Return the current frame number which Smokeivew has loaded.
*/
int lua_getframe(lua_State *L) {
  int framenumber = getframe();
  // Push a return value to the Lua stack.
  lua_pushinteger(L, framenumber);
  // Tell Lua that there is a single return value left on the stack.
  return 1;
}

/*
  Shift to a specific frame number.
*/
int lua_setframe(lua_State *L) {
  int f = lua_tonumber(L, 1);
  printf("lua_api: setting frame to %d\n", f);
  setframe(f);
  return 0;
}

/*
  Get the time value of the currently loaded frame.
*/
int lua_gettime(lua_State *L) {
  if(global_times!=NULL&&nglobal_times>0){
    float time = gettime();
    lua_pushnumber(L, time);
    return 1;
  } else {
    return 0;
  }

}

/*
  Shift to the closest frame to given a time value.
*/
int lua_settime(lua_State *L) {
  lua_displayCB(L);
  float t = lua_tonumber(L, 1);
  int return_code = settime(t);
  lua_pushnumber(L, return_code);
  return 1;
}

/*
  Load an FDS data file directly (i.e. as a filepath).
*/
int lua_loaddatafile(lua_State *L) {
  const char *filename = lua_tostring(L, 1);
  int return_value = loadfile(filename);
  lua_pushnumber(L, return_value);
  return 1;
}

/*
  Load a Smokeview config (.ini) file.
*/
int lua_loadinifile(lua_State *L) {
  const char *filename = lua_tostring(L, 1);
  loadinifile(filename);
  return 0;
}

/*
  Load an FDS vector data file directly (i.e. as a filepath). This function
  handles the loading of any additional data files necessary to display vectors.
*/
int lua_loadvdatafile(lua_State *L) {
  const char *filename = lua_tostring(L, 1);
  int return_value = loadvfile(filename);
  lua_pushnumber(L, return_value);
  return 1;
}

/*
  Load an FDS boundary file directly (i.e. as a filepath). This is equivalent
  to lua_loadfile, but specialised for boundary files. This is included to
  reflect the underlying code.
*/
int lua_loadboundaryfile(lua_State *L) {
  const char *filename = lua_tostring(L, 1);
  loadboundaryfile(filename);
  return 0;
}

/*
  Print a label to stdout.
*/
int lua_label(lua_State *L) {
  const char *thelabel = lua_tostring(L, 1);
  label(thelabel);
  return 0;
}

/*
  Load a slice file given the type of slice, the axis along which it exists and
  its position along this axis.
*/
int lua_loadslice(lua_State *L) {
  const char *type = lua_tostring(L, 1);
  int axis = lua_tonumber(L, 2);
  float distance = lua_tonumber(L, 3);
  loadslice(type, axis, distance);
  return 0;
}

/*
  Set the clipping mode, which determines which parts of the model are clipped
  (based on the set clipping values). This function takes an int, which is one
  of:
    0: No clipping.
    1: Clip blockages and data.
    2: Clip blockages.
    3: Clip data.
*/
int lua_set_clipping_mode(lua_State *L) {
  int mode = lua_tonumber(L, 1);
  set_clipping_mode(mode);
  return 0;
}

int lua_set_sceneclip_x(lua_State *L) {
  int clipMin = lua_toboolean(L, 1);
  float min = lua_tonumber(L, 2);
  int clipMax = lua_toboolean(L, 3);
  float max = lua_tonumber(L, 4);
  set_sceneclip_x(clipMin, min, clipMax, max);
  return 0;
}

int lua_set_sceneclip_x_min(lua_State *L) {
  int flag = lua_toboolean(L, 1);
  float value = lua_tonumber(L, 2);
  set_sceneclip_x_min(flag, value);
  return 0;
}

int lua_set_sceneclip_x_max(lua_State *L) {
  int flag = lua_toboolean(L, 1);
  float value = lua_tonumber(L, 2);
  set_sceneclip_x_max(flag, value);
  return 0;
}

int lua_set_sceneclip_y(lua_State *L) {
  int clipMin = lua_toboolean(L, 1);
  float min = lua_tonumber(L, 2);
  int clipMax = lua_toboolean(L, 3);
  float max = lua_tonumber(L, 4);
  set_sceneclip_y(clipMin, min, clipMax, max);
  return 0;
}

int lua_set_sceneclip_y_min(lua_State *L) {
  int flag = lua_toboolean(L, 1);
  float value = lua_tonumber(L, 2);
  set_sceneclip_y_min(flag, value);
  return 0;
}

int lua_set_sceneclip_y_max(lua_State *L) {
  int flag = lua_toboolean(L, 1);
  float value = lua_tonumber(L, 2);
  set_sceneclip_y_max(flag, value);
  return 0;
}

int lua_set_sceneclip_z(lua_State *L) {
  int clipMin = lua_toboolean(L, 1);
  float min = lua_tonumber(L, 2);
  int clipMax = lua_toboolean(L, 3);
  float max = lua_tonumber(L, 4);
  set_sceneclip_z(clipMin, min, clipMax, max);
  return 0;
}

int lua_set_sceneclip_z_min(lua_State *L) {
  int flag = lua_toboolean(L, 1);
  float value = lua_tonumber(L, 2);
  set_sceneclip_z_min(flag, value);
  return 0;
}

int lua_set_sceneclip_z_max(lua_State *L) {
  int flag = lua_toboolean(L, 1);
  float value = lua_tonumber(L, 2);
  set_sceneclip_z_max(flag, value);
  return 0;
}

/*
  Return a table (an array) of the times available in Smokeview. They key of the
  table is an int representing the frame number, and the value of the table is
  a float representing the time.
*/
int lua_get_global_times(lua_State *L) {
  PRINTF("lua: initialising global time table\n");
  lua_createtable(L, 0, nglobal_times);
  int i;
    for (i = 0; i < nglobal_times; i++) {
      lua_pushnumber(L, i);
      lua_pushnumber(L, global_times[i]);
      lua_settable(L, -3);
    }
    return 1;
}

/*
  Get the number of (global) frames available to smokeview.
*/
int lua_get_nglobal_times(lua_State *L) {
  lua_pushnumber(L, nglobal_times);
  return 1;
}

/*
  Get the number of meshes in the loaded model.
*/
int lua_get_nmeshes(lua_State *L) {
  lua_pushnumber(L, nmeshes);
  return 1;
}

/*
  Build a Lua table with information on the meshes of the model. The key of the
  table is the mesh number.
*/
// TODO: provide more information via this interface.
int lua_get_meshes(lua_State *L) {
  int entries = nmeshes;
  meshdata *infotable = meshinfo;
  PRINTF("lua: initialising mesh table\n");
  lua_createtable(L, 0, entries);
  int i;
  for (i = 0; i < entries; i++) {
    lua_pushnumber(L, i);
    lua_createtable(L, 0, 2);

    lua_pushnumber(L, infotable[i].ibar);
    lua_setfield(L, -2, "ibar");

    lua_pushnumber(L, infotable[i].jbar);
    lua_setfield(L, -2, "ibar");

    lua_pushnumber(L, infotable[i].kbar);
    lua_setfield(L, -2, "ibar");

    lua_settable(L, -3);
  }
  PRINTF("lua: done initialising mesh table\n");
  // Leaves one returned value on the stack, the mesh table.
  return 1;
}

/*
  Get the number of meshes in the loaded model.
*/
int lua_get_ndevices(lua_State *L) {
    lua_pushnumber(L, ndeviceinfo);
    return 1;
}

/*
  Build a Lua table with information on the devices of the model.
*/
int lua_get_devices(lua_State *L) {
  int entries = ndeviceinfo;
  devicedata *infotable = deviceinfo;
  PRINTF("lua: initialising device table\n");
  lua_createtable(L, 0, entries);
  int i;
  for (i = 0; i < entries; i++) {
    lua_pushnumber(L, i);
    lua_createtable(L, 0, 2);

    lua_pushstring(L, infotable[i].label);
    lua_setfield(L, -2, "label");

    lua_settable(L, -3);
  }
  return 1;
}

/*
  Get the number of CSV files available to the model.
*/
int lua_get_ncsvinfo(lua_State*L) {
  lua_pushnumber(L, ncsvinfo);
  return 1;
}

/*
  Load data about the loaded module into the lua interpreter.
  This initsmvdata is necessary to bring some data into the Lua interpreter
  from the model. This is included here rather than doing in the Smokeview
  code to increase separation. This will likely be removed in future versions.
*/
// TODO: Consider converting most of these to userdata, rather than copying them
// into the lua interpreter.
int lua_initsmvdata(lua_State *L) {
  lua_get_nglobal_times(L);
  lua_setglobal(L, "nglobal_times");
  lua_get_global_times(L);
  lua_setglobal(L, "global_times");

  lua_get_nmeshes(L);
  lua_setglobal(L, "nmeshes");
  lua_get_meshes(L);
  lua_setglobal(L, "meshinfo");

  lua_get_ndevices(L);
  lua_setglobal(L, "ndevices");
  lua_get_devices(L);
  lua_setglobal(L, "deviceinfo");

  lua_get_sliceinfo(L);
  lua_setglobal(L, "sliceinfo");

  lua_get_csvinfo(L);
  lua_setglobal(L, "csvinfo");
  // lua_get_geomdata(L);
  // lua_setglobal(L, "geomdata");
  return 0;
}

/*
  As with lua_initsmvdata(), but for information relating to Smokeview itself.
*/
int lua_initsmvproginfo(lua_State *L) {
  char version[256];
  char githash[256];

  getPROGversion(version);
  // getGitHash(githash);

  lua_createtable(L, 0, 6);

  lua_pushstring(L, version);
  lua_setfield(L, -2, "version");

  // lua_pushstring(L, githash);
  // lua_setfield(L, -2, "githash");

  lua_pushstring(L, __DATE__);
  lua_setfield(L, -2, "builddate");

  lua_pushstring(L, fds_githash);
  lua_setfield(L, -2, "fdsgithash");

  lua_pushstring(L, smokeviewpath);
  lua_setfield(L, -2, "smokeviewpath");

  lua_pushstring(L, smokezippath);
  lua_setfield(L, -2, "smokezippath");

  lua_pushstring(L, texturedir);
  lua_setfield(L, -2, "texturedir");

  lua_setglobal(L, "smokeviewProgram");
  return 0;
}

/*
  Build a Lua table with information on the slices of the model.
*/
// TODO: provide more information via this interface.
int lua_get_sliceinfo(lua_State *L) {
  PRINTF("lua: initialising slice table\n");
  lua_createtable(L, 0, nsliceinfo);
  int i;
  for (i = 0; i < nsliceinfo; i++) {
    lua_pushnumber(L, i);
    lua_createtable(L, 0, 16);

    if(sliceinfo[i].slicelabel != NULL) {
      lua_pushstring(L, sliceinfo[i].slicelabel);
      lua_setfield(L, -2, "label");
    }

    if(sliceinfo[i].label.longlabel != NULL) {
      lua_pushstring(L, sliceinfo[i].label.longlabel);
      lua_setfield(L, -2, "longlabel");
    }

    if(sliceinfo[i].label.shortlabel != NULL) {
      lua_pushstring(L, sliceinfo[i].label.shortlabel);
      lua_setfield(L, -2, "shortlabel");
    }

    lua_pushstring(L, sliceinfo[i].file);
    lua_setfield(L, -2, "file");

    lua_pushnumber(L, sliceinfo[i].slicetype);
    lua_setfield(L, -2, "slicetype");

    lua_pushnumber(L, sliceinfo[i].idir);
    lua_setfield(L, -2, "idir");

    lua_pushnumber(L, sliceinfo[i].sliceoffset);
    lua_setfield(L, -2, "sliceoffset");

    lua_pushnumber(L, sliceinfo[i].ijk_min[0]);
    lua_setfield(L, -2, "imin");

    lua_pushnumber(L, sliceinfo[i].ijk_max[0]);
    lua_setfield(L, -2, "imax");

    lua_pushnumber(L, sliceinfo[i].ijk_min[1]);
    lua_setfield(L, -2, "jmin");

    lua_pushnumber(L, sliceinfo[i].ijk_max[1]);
    lua_setfield(L, -2, "jmax");

    lua_pushnumber(L, sliceinfo[i].ijk_min[2]);
    lua_setfield(L, -2, "kmin");

    lua_pushnumber(L, sliceinfo[i].ijk_max[2]);
    lua_setfield(L, -2, "kmax");

    lua_pushnumber(L, sliceinfo[i].blocknumber);
    lua_setfield(L, -2, "blocknumber");

    lua_pushnumber(L, sliceinfo[i].position_orig);
    lua_setfield(L, -2, "position_orig");


    lua_pushstring(L, sliceinfo[i].slicedir);
    lua_setfield(L, -2, "slicedir");

    lua_settable(L, -3);
  }
  return 1;
}

/*
  Build a Lua table with information on the CSV files available to the model.
*/
// TODO: provide more information via this interface.
int lua_get_csvinfo(lua_State *L) {
  PRINTF("lua: initialising csv table\n");
  lua_createtable(L, 0, ncsvinfo);
  int i;
  for (i = 0; i < ncsvinfo; i++) {
    lua_pushnumber(L, i);
    lua_createtable(L, 0, 4);

    lua_pushstring(L, csvinfo[i].file);
    lua_setfield(L, -2, "file");

    lua_pushboolean(L, csvinfo[i].loaded);
    lua_setfield(L, -2, "loaded");

    lua_pushboolean(L, csvinfo[i].display);
    lua_setfield(L, -2, "display");

    lua_pushnumber(L, csvinfo[i].type);
    lua_setfield(L, -2, "type");

    lua_settable(L, -3);
  }
  return 1;
}


int lua_loadvslice(lua_State *L) {
  const char *type = lua_tostring(L, 1);
  int axis = lua_tonumber(L, 2);
  float distance = lua_tonumber(L, 3);
  loadvslice(type, axis, distance);
  return 0;
}

int lua_loadiso(lua_State *L) {
  const char *type = lua_tostring(L, 1);
  loadiso(type);
  return 0;
}

int lua_load3dsmoke(lua_State *L) {
  const char *smoke_type = lua_tostring(L, 1);
  load3dsmoke(smoke_type);
  return 0;
}

int lua_loadvolsmoke(lua_State *L) {
  int meshnumber = lua_tonumber(L, 1);
  loadvolsmoke(meshnumber);
  return 0;
}

int lua_loadvolsmokeframe(lua_State *L) {
  int meshnumber = lua_tonumber(L, 1);
  int framenumber = lua_tonumber(L, 1);
  loadvolsmokeframe(meshnumber, framenumber, 1);
    // returnval = 1; // TODO: determine if this is the correct behaviour.
                    // this is what is done in the SSF code.
  return 0;
}

/*
  Set the format of images which will be exported. The value should be a string.
  The acceptable values are:
    "JPG"
    "PNG"
*/
int lua_set_rendertype(lua_State *L) {
  const char *type = lua_tostring(L, 1);
  rendertype(type);
  return 0;
}

int lua_get_rendertype(lua_State *L) {
  int render_type = get_rendertype();
  switch (render_type) {
    case JPEG:
      lua_pushstring(L, "JPG");
      break;
    case PNG:
      lua_pushstring(L, "PNG");
      break;
    default:
      lua_pushstring(L, NULL);
      break;
  }
  return 1;
}

/*
  Set the format of movies which will be exported. The value should be a string.
  The acceptable values are:
    "WMV"
    "MP4"
    "AVI"
*/
int lua_set_movietype(lua_State *L) {
  const char *type = lua_tostring(L, 1);
  set_movietype(type);
  return 0;
}

int lua_get_movietype(lua_State *L) {
  int movie_type = get_movietype();
  switch (movie_type) {
    case WMV:
        lua_pushstring(L, "WMV");
        break;
    case MP4:
        lua_pushstring(L, "MP4");
        break;
    case AVI:
        lua_pushstring(L, "AVI");
        break;
    default:
        lua_pushstring(L, NULL);
        break;
  }
  return 1;
}

int lua_makemovie(lua_State *L) {
  const char *name = lua_tostring(L, 1);
  const char *base = lua_tostring(L, 2);
  float framerate = lua_tonumber(L, 3);
  makemovie(name, base, framerate);
}

int lua_loadtour(lua_State *L) {
  const char *name = lua_tostring(L, 1);
  int error_code = loadtour(name);
  lua_pushnumber(L, error_code);
  return 1;
}

int lua_loadparticles(lua_State *L) {
  const char *name = lua_tostring(L, 1);
  loadparticles(name);
  return 0;
}

int lua_partclasscolor(lua_State *L) {
  const char *color = lua_tostring(L, 1);
  partclasscolor(color);
  return 0;
}

int lua_partclasstype(lua_State *L) {
  const char *type = lua_tostring(L, 1);
  partclasstype(type);
  return 0;
}

int lua_plot3dprops(lua_State *L) {
  int variable_index = lua_tonumber(L, 1);
  int showvector = lua_toboolean(L, 2);
  int vector_length_index = lua_tonumber(L, 3);
  int display_type = lua_tonumber(L, 4);
  float vector_length = lua_tonumber(L, 5);
  plot3dprops(variable_index, showvector, vector_length_index, display_type,
              vector_length);
  return 0;
}

int lua_loadplot3d(lua_State *L) {
  int meshnumber = lua_tonumber(L, 1);
  float time_local = lua_tonumber(L, 2);
  loadplot3d(meshnumber, time_local);
  return 0;
}

int lua_unloadall(lua_State *L) {
  unloadall();
  return 0;
}

int lua_unloadtour(lua_State *L) {
  unloadtour();
  return 0;
}

int lua_setrenderdir(lua_State *L) {
  const char *dir = lua_tostring(L, 1);
  int return_code = setrenderdir(dir);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_getrenderdir(lua_State *L) {
  lua_pushstring(L, script_dir_path);
  return 1;
}

int lua_setviewpoint(lua_State *L) {
  const char *viewpoint = lua_tostring(L, 1);
  int errorcode = setviewpoint(viewpoint);
  lua_pushnumber(L, errorcode);
  return 1;
}

int lua_getviewpoint(lua_State *L) {
  lua_pushstring(L, camera_current->name);
  return 1;
}

int lua_exit_smokeview(lua_State  *L) {
  exit_smokeview();
  return 0;
}

int lua_setwindowsize(lua_State *L) {
  int width = lua_tonumber(L, 1);
  int height = lua_tonumber(L, 2);
  setwindowsize(width, height);
  // Using the Display_CB is not sufficient in this case,
  // control must be temporarily returned to the main glut loop.
  lua_tempyieldscript(L);
  return 0;
}

int lua_setgridvisibility(lua_State *L) {
  int selection = lua_tonumber(L, 1);
  setgridvisibility(selection);
  return 0;
}

int lua_setgridparms(lua_State *L) {
  int x_vis = lua_tonumber(L, 1);
  int y_vis = lua_tonumber(L, 2);
  int z_vis = lua_tonumber(L, 3);

  int x_plot = lua_tonumber(L, 4);
  int y_plot = lua_tonumber(L, 5);
  int z_plot = lua_tonumber(L, 6);

  setgridparms(x_vis, y_vis, z_vis, x_plot, y_plot, z_plot);

  return 0;
}

int lua_setcolorbarflip(lua_State *L) {
  int flip = lua_toboolean(L, 1);
  setcolorbarflip(flip);
  return 0;
}

int lua_getcolorbarflip(lua_State *L) {
  int flip = getcolorbarflip();
  lua_pushboolean(L, flip);
  return 1;
}

int lua_setcolorbarindex(lua_State *L) {
  int chosen_index = lua_tonumber(L, 1);
  setcolorbarindex(chosen_index);
  return 0;
}

int lua_getcolorbarindex(lua_State *L) {
  int index = getcolorbarindex();
  lua_pushnumber(L, index);
  return 1;
}

int lua_set_slice_in_obst(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_slice_in_obst(setting);
  return 0;
}

int lua_get_slice_in_obst(lua_State *L) {
  int setting = get_slice_in_obst();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_setcolorbarvisibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_colorbar_visibility(setting);
  return 0;
}

int lua_getcolorbarvisibility(lua_State *L) {
  int setting = get_colorbar_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_togglecolorbarvisibility(lua_State *L) {
  toggle_colorbar_visibility();
  return 0;
}

int lua_settimebarvisibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_timebar_visibility(setting);
  return 0;
}

int lua_gettimebarvisibility(lua_State *L) {
  int setting = get_timebar_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggletimebarvisibility(lua_State *L) {
  toggle_timebar_visibility();
  return 0;
}

// HRR Label Visbility
int lua_set_hrrlabel_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_hrrlabel_visibility(setting);
  return 0;
}

int lua_get_hrrlabel_visibility(lua_State *L) {
  int setting = get_hrrlabel_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_hrrlabel_visibility(lua_State *L) {
  toggle_hrrlabel_visibility();
  return 0;
}

int lua_camera_mod_eyex(lua_State *L) {
  float delta = lua_tonumber(L, 1);
  camera_mod_eyex(delta);
  return 0;
}

int lua_camera_set_eyex(lua_State *L) {
  float eyex = lua_tonumber(L, 1);
  camera_set_eyex(eyex);
  return 0;
}

int lua_camera_mod_eyey(lua_State *L) {
  float delta = lua_tonumber(L, 1);
  camera_mod_eyey(delta);
  return 0;
}

int lua_camera_set_eyey(lua_State *L) {
  float eyey = lua_tonumber(L, 1);
  camera_set_eyey(eyey);
  return 0;
}

int lua_camera_mod_eyez(lua_State *L) {
  float delta = lua_tonumber(L, 1);
  camera_mod_eyez(delta);
  return 0;
}

int lua_camera_set_eyez(lua_State *L) {
  float eyez = lua_tonumber(L, 1);
  camera_set_eyez(eyez);
  return 0;
}

int lua_camera_mod_az(lua_State *L) {
  float delta = lua_tonumber(L, 1);
  camera_mod_az(delta);
  return 0;
}

int lua_camera_set_az(lua_State *L) {
  float az = lua_tonumber(L, 1);
  camera_set_az(az);
  return 0;
}

int lua_camera_get_az(lua_State *L) {
  lua_pushnumber(L, camera_get_az());
  return 1;
}

int lua_camera_mod_elev(lua_State *L) {
  float delta = lua_tonumber(L, 1);
  camera_mod_elev(delta);
  return 0;
}

int lua_camera_set_elev(lua_State *L) {
  float elev = lua_tonumber(L, 1);
  camera_set_elev(elev);
  return 0;
}

int lua_camera_get_elev(lua_State *L) {
  lua_pushnumber(L, camera_get_elev());
  return 1;
}
int lua_camera_get_projection_type(lua_State *L) {
  float projection_type = camera_get_projection_type();
  lua_pushnumber(L, projection_type);
  return 1;
}
int lua_camera_set_projection_type(lua_State *L) {
  float projection_type = lua_tonumber(L, 1);
  int return_value = camera_set_projection_type(projection_type);
  lua_pushnumber(L, return_value);
  return 1;
}

int lua_camera_get_rotation_type(lua_State *L) {
  float rotation_type = camera_get_rotation_type();
  lua_pushnumber(L, rotation_type);
  return 1;
}

int lua_camera_get_rotation_index(lua_State *L) {
  float rotation_index = camera_get_rotation_index();
  lua_pushnumber(L, rotation_index);
  return 1;
}

int lua_camera_set_rotation_type(lua_State *L) {
  float rotation_type = lua_tonumber(L, 1);
  camera_set_rotation_type(rotation_type);
  return 0;
}

int lua_camera_get_zoom(lua_State *L) {
  lua_pushnumber(L, zoom);
  return 1;
}

int lua_camera_set_zoom(lua_State *L) {
  float x = lua_tonumber(L, 1);
  zoom = x;
  return 0;
}

int lua_camera_get_eyex(lua_State *L) {
  float eyex = camera_get_eyex();
  lua_pushnumber(L, eyex);
  return 1;
}

int lua_camera_get_eyey(lua_State *L) {
  float eyey = camera_get_eyex();
  lua_pushnumber(L, eyey);
  return 1;
}

int lua_camera_get_eyez(lua_State *L) {
  float eyez = camera_get_eyez();
  lua_pushnumber(L, eyez);
  return 1;
}
int lua_camera_set_viewdir(lua_State *L) {
  float xcen = lua_tonumber(L, 1);
  float ycen = lua_tonumber(L, 2);
  float zcen = lua_tonumber(L, 3);
  printf("lua_api: Setting viewDir to %f %f %f\n", xcen, ycen, zcen);
  camera_set_viewdir(xcen, ycen, zcen);
  return 0;
}

int lua_camera_get_viewdir(lua_State *L) {
  float xcen = camera_get_xcen();
  float ycen = camera_get_ycen();
  float zcen = camera_get_zcen();

  lua_createtable(L, 0, 3);

  lua_pushstring(L, "x");
  lua_pushnumber(L, xcen);
  lua_settable(L, -3);

  lua_pushstring(L, "y");
  lua_pushnumber(L, ycen);
  lua_settable(L, -3);

  lua_pushstring(L, "z");
  lua_pushnumber(L, zcen);
  lua_settable(L, -3);

  return 1;
}

int lua_set_slice_bound_min(lua_State *L) {
  const char *slice_type = lua_tostring(L, 1);
  int set = lua_toboolean(L, 2);
  float value = lua_tonumber(L, 3);
  printf("lua: setting %s min bound ", slice_type);
  if(set) {printf("ON");} else {printf("OFF");}
  printf(" with value of %f\n", value);
  set_slice_bound_min(slice_type, set, value);
  return 0;
}

int lua_get_slice_bound_min(lua_State *L) {
  const char *slice_type = lua_tostring(L, 1);
  float value = get_slice_bound_min(slice_type);
  lua_pushnumber(L, value);
  return 1;
}

int lua_get_slice_bound_max(lua_State *L) {
  const char *slice_type = lua_tostring(L, 1);
  float value = get_slice_bound_max(slice_type);
  lua_pushnumber(L, value);
  return 1;
}

int lua_set_slice_bound_max(lua_State *L) {
  const char *slice_type = lua_tostring(L, 1);
  int set = lua_toboolean(L, 2);
  float value = lua_tonumber(L, 3);
  printf("lua: setting %s max bound ", slice_type);
  if(set) {printf("ON");} else {printf("OFF");}
  printf(" with value of %f\n", value);
  set_slice_bound_max(slice_type, set, value);
  return 0;
}

// add the smokeview bin directory to the Lua path variables
void addLuaPaths() {
  // package.path is a path variable where Lua scripts and modules may be
  // found, typiclly text based files with the .lua extension.
  lua_getglobal(L, "package");
  int pathType = lua_getfield(L, -1, "path");
  const char *oldPath = lua_tostring(L, -1);
  int newLength = strlen(oldPath) + 1 + strlen(smokeview_bindir) + 5 +1;
  char newPath[newLength];
  strcpy(newPath, oldPath);
  strcat(newPath,";");
  strcat(newPath,smokeview_bindir);
  strcat(newPath,"?.lua");
  lua_pushstring(L, newPath);
  lua_setfield(L, -3, "path");
  lua_pop(L, 1); // pop the now redundant "path" variable from the stack

  // package.cpath is a path variable where Lua modules may be found,
  // typically binary (C based) files such as .dll or .so.
  int cpathType = lua_getfield(L, -1, "cpath");
  const char *oldCPath = lua_tostring(L, -1);
  int newLengthC = strlen(oldCPath) + 1 + 2*strlen(smokeview_bindir) + 10 +1;
  char newCPath[newLengthC];
  strcpy(newCPath, oldCPath);
  strcat(newCPath,";");
  strcat(newCPath,smokeview_bindir);
  strcat(newCPath,"?.dll;");
  strcat(newCPath,smokeview_bindir);
  strcat(newCPath,"?.so");
  lua_pushstring(L, newCPath);
  lua_setfield(L, -3, "cpath");
  lua_pop(L, 1); // pop the now redundant "cpath" variable from the stack
  return;
}

void initLua() {
  L = luaL_newstate();

  luaL_openlibs(L);
  lua_initsmvproginfo(L);
  addLuaPaths();
  lua_register(L, "set_slice_bound_min", lua_set_slice_bound_min);
  lua_register(L, "set_slice_bound_max", lua_set_slice_bound_max);
  lua_register(L, "get_slice_bound_min", lua_get_slice_bound_min);
  lua_register(L, "get_slice_bound_max", lua_get_slice_bound_max);
  lua_register(L, "loadsmvall", lua_loadsmvall);
  lua_register(L, "hidewindow", lua_hidewindow);
  lua_register(L, "yieldscript", lua_yieldscript);
  lua_register(L, "tempyieldscript", lua_tempyieldscript);
  lua_register(L, "displayCB", lua_displayCB);
  lua_register(L, "renderclip", lua_renderclip);
  lua_register(L, "renderC", lua_render);
  lua_register(L, "gsliceview", lua_gsliceview);
  lua_register(L, "gslicepos", lua_gslicepos);
  lua_register(L, "gsliceorien", lua_gsliceorien);
  lua_register(L, "settourkeyframe", lua_settourkeyframe);
  lua_register(L, "settourview", lua_settourview);
  lua_register(L, "getframe", lua_getframe);
  lua_register(L, "setframe", lua_setframe);
  lua_register(L, "gettime", lua_gettime);
  lua_register(L, "settime", lua_settime);
  lua_register(L, "loaddatafile", lua_loaddatafile);
  lua_register(L, "loadinifile", lua_loadinifile);
  lua_register(L, "loadvdatafile", lua_loadvdatafile);
  lua_register(L, "loadboundaryfile", lua_loadboundaryfile);
  lua_register(L, "label", lua_label);
  lua_register(L, "load3dsmoke", lua_load3dsmoke);
  lua_register(L, "loadvolsmoke", lua_loadvolsmoke);
  lua_register(L, "loadvolsmokeframe", lua_loadvolsmokeframe);
  lua_register(L, "set_rendertype", lua_set_rendertype);
  lua_register(L, "get_rendertype", lua_get_rendertype);
  lua_register(L, "set_movietype", lua_set_movietype);
  lua_register(L, "get_movietype", lua_get_movietype);
  lua_register(L, "makemovie", lua_makemovie);
  lua_register(L, "loadtour", lua_loadtour);
  lua_register(L, "loadparticles", lua_loadparticles);
  lua_register(L, "partclasscolor", lua_partclasscolor);
  lua_register(L, "partclasstype", lua_partclasstype);
  lua_register(L, "plot3dprops", lua_plot3dprops);
  lua_register(L, "loadplot3d", lua_loadplot3d);
  lua_register(L, "loadslice", lua_loadslice);
  // lua_register(L, "loadnamedslice", lua_loadnamedslice);
  lua_register(L, "loadvslice", lua_loadvslice);
  lua_register(L, "loadiso", lua_loadiso);
  lua_register(L, "unloadall", lua_unloadall);
  lua_register(L, "unloadtour", lua_unloadtour);
  lua_register(L, "setrenderdir", lua_setrenderdir);
  lua_register(L, "getrenderdir", lua_getrenderdir);
  lua_register(L, "setviewpoint", lua_setviewpoint);
  lua_register(L, "getviewpoint", lua_getviewpoint);
  lua_register(L, "exit", lua_exit_smokeview);
  lua_register(L, "getcolorbarflip", lua_getcolorbarflip);
  lua_register(L, "setcolorbarflip", lua_setcolorbarflip);
  lua_register(L, "setwindowsize", lua_setwindowsize);
  lua_register(L, "setgridvisibility", lua_setgridvisibility);
  lua_register(L, "setgridparms", lua_setgridparms);
  lua_register(L, "setcolorbarindex", lua_setcolorbarindex);
  lua_register(L, "getcolorbarindex", lua_getcolorbarindex);

  lua_register(L, "set_slice_in_obst", lua_set_slice_in_obst);
  lua_register(L, "get_slice_in_obst", lua_get_slice_in_obst);

  lua_register(L, "setcolorbarvisibility", lua_setcolorbarvisibility);
  lua_register(L, "getcolorbarvisibility", lua_getcolorbarvisibility);
  lua_register(L, "togglecolorbarvisibility", lua_togglecolorbarvisibility);

  lua_register(L, "settimebarvisibility", lua_settimebarvisibility);
  lua_register(L, "gettimebarvisibility", lua_gettimebarvisibility);
  lua_register(L, "toggletimebarvisibility", lua_toggletimebarvisibility);

  lua_register(L, "set_hrrlabel_visibility", lua_set_hrrlabel_visibility);
  lua_register(L, "get_hrrlabel_visibility", lua_get_hrrlabel_visibility);
  lua_register(L, "toggle_hrrlabel_visibility", lua_toggle_hrrlabel_visibility);

  lua_register(L, "camera_mod_eyex", lua_camera_mod_eyex);
  lua_register(L, "camera_set_eyex", lua_camera_set_eyex);
  lua_register(L, "camera_get_eyex", lua_camera_get_eyex);

  lua_register(L, "camera_mod_eyey", lua_camera_mod_eyey);
  lua_register(L, "camera_set_eyey", lua_camera_set_eyey);
  lua_register(L, "camera_get_eyey", lua_camera_get_eyey);

  lua_register(L, "camera_mod_eyez", lua_camera_mod_eyez);
  lua_register(L, "camera_set_eyez", lua_camera_set_eyez);
  lua_register(L, "camera_get_eyez", lua_camera_get_eyez);

  lua_register(L, "camera_mod_az", lua_camera_mod_az);
  lua_register(L, "camera_set_az", lua_camera_set_az);
  lua_register(L, "camera_get_az", lua_camera_get_az);
  lua_register(L, "camera_mod_elev", lua_camera_mod_elev);
  lua_register(L, "camera_set_elev", lua_camera_set_elev);
  lua_register(L, "camera_get_elev", lua_camera_get_elev);

  lua_register(L, "camera_set_viewdir", lua_camera_set_viewdir);
  lua_register(L, "camera_get_viewdir", lua_camera_get_viewdir);

  lua_register(L, "camera_get_zoom", lua_camera_get_zoom);
  lua_register(L, "camera_set_zoom", lua_camera_set_zoom);

  lua_register(L, "camera_get_rotation_type"
      , lua_camera_get_rotation_type);
      lua_register(L, "camera_get_rotation_index"
          , lua_camera_get_rotation_index);
  lua_register(L, "camera_set_rotation_type"
      , lua_camera_set_rotation_type);
  lua_register(L, "camera_get_projection_type"
      , lua_camera_get_projection_type);
  lua_register(L, "camera_set_projection_type"
      , lua_camera_set_projection_type);

  // lua_register(L, "get_clipping_mode", lua_get_clipping_mode);
  lua_register(L, "set_clipping_mode", lua_set_clipping_mode);
  lua_register(L, "set_sceneclip_x", lua_set_sceneclip_x);
  lua_register(L, "set_sceneclip_x_min", lua_set_sceneclip_x_min);
  lua_register(L, "set_sceneclip_x_max", lua_set_sceneclip_x_max);
  lua_register(L, "set_sceneclip_y", lua_set_sceneclip_y);
  lua_register(L, "set_sceneclip_y_min", lua_set_sceneclip_y_min);
  lua_register(L, "set_sceneclip_y_max", lua_set_sceneclip_y_max);
  lua_register(L, "set_sceneclip_z", lua_set_sceneclip_z);
  lua_register(L, "set_sceneclip_z_min", lua_set_sceneclip_z_min);
  lua_register(L, "set_sceneclip_z_max", lua_set_sceneclip_z_max);


  lua_register(L, "get_nglobal_times", lua_get_nglobal_times);

  //add fdsprefix as a variable in the lua environment
  lua_pushstring(L, fdsprefix);
  lua_setglobal(L, "fdsprefix");

  //nglobal_times is the number of frames
  // this cannot be set as a global at init as it will change
  // on the loading of smokeview cases
  // TODO: possibly change this to initialise the value when a new
  // smokeview case is loaded, in which case we will have a lua_case_init
  // function, to initialise these types of variables
  lua_register(L, "initsmvdata", lua_initsmvdata);

  lua_pushstring(L, script_dir_path);
  lua_setglobal(L, "current_script_dir");
  //lua_pushstring(L, renderfile_dir);
  //lua_setglobal(L, "current_render_dir");

  // a boolean value that determines if lua is running in smokeview
  lua_pushboolean(L, 1);
  lua_setglobal(L, "smokeviewEmbedded");

  // luaL_requiref (L, "smv", luaopen_smv, 1)
  luaL_dostring(L, "require(\"smv\")");
  // luaL_loadfile(L, "smv.lua");
  // int luaL_1dofile (lua_State *L, const char *filename);

}

// int luaopen_smv(lua_State *L){
//  static const luaL_Reg Obj_lib[] = {
//    { "method", &Obj_method },
//    { NULL, NULL }
//  };
//
//  static const luaL_Reg MyLib_lib[] = {
//    { "MakeObj", &MyLib_MakeObj },
//    { NULL, NULL }
//  };
//
//  luaL_newlib(L, MyLib_lib);
//
//  // Stack: MyLib
//  luaL_newmetatable(L, Obj_typename); // Stack: MyLib meta
//  luaL_newlib(L, Obj_lib);
//  lua_setfield(L, -2, "__index"); // Stack: MyLib meta
//
//  lua_pushstring(L, "__gc");
//  lua_pushcfunction(L, Obj__gc); // Stack: MyLib meta "__gc" fptr
//  lua_settable(L, -3); // Stack: MyLib meta
//  lua_pop(L, 1); // Stack: MyLib
//
//  return 1;
// }

void runScriptString(char *string) {
  luaL_dostring(L, string);
//  return 0;
}

int loadLuaScript(char *filename) {
  printf("scriptfile: %s\n", filename);
  // The display callback needs to be run once initially.
  // PROBLEM: the display CB does not work without a loaded case.
  runluascript=0;
  lua_displayCB(L);
  runluascript=1;
  printf("loading: %s\n", filename);
  const char *err_msg;
  lua_Debug info;
  int level = 0;
  int return_code = luaL_loadfile(L, filename);
  switch (return_code) {
    case LUA_OK:
      printf("%s loaded ok\n", filename);
      break;
    case LUA_ERRSYNTAX:
      fprintf(stderr, "Syntax error loading %s\n", filename);
      err_msg = lua_tostring (L, -1);
      fprintf(stderr, "error:%s\n", err_msg);
      level = 0;
      while (lua_getstack(L, level, &info)) {
        lua_getinfo(L, "nSl", &info);
        fprintf(stderr, "  [%d] %s:%d -- %s [%s]\n",
          level, info.short_src, info.currentline,
          (info.name ? info.name : "<unknown>"), info.what);
        ++level;
      }
      break;
    case LUA_ERRMEM:
        break;
    case LUA_ERRGCMM:
        break;
    case LUA_ERRFILE:
      fprintf(stderr, "Could not load file %s\n", filename);
      err_msg = lua_tostring (L, -1);
      fprintf(stderr, "error:%s\n", err_msg);
      level = 0;
      while (lua_getstack(L, level, &info)) {
        lua_getinfo(L, "nSl", &info);
        fprintf(stderr, "  [%d] %s:%d -- %s [%s]\n",
            level, info.short_src, info.currentline,
            (info.name ? info.name : "<unknown>"), info.what);
        ++level;
      }
      break;
  }
  printf("after lua loadfile\n");
  return return_code;
}

int loadSSFScript(char *filename) {
  // char filename[1024];
  //   if (strlen(script_filename) == 0) {
  //       strncpy(filename, fdsprefix, 1020);
  //       strcat(filename, ".ssf");
  //   } else {
  //       strncpy(filename, script_filename, 1024);
  //   }
  printf("scriptfile: %s\n", filename);
    // The display callback needs to be run once initially.
    // PROBLEM: the display CB does not work without a loaded case.
    runscript=0;
    lua_displayCB(L);
    runscript=1;
    const char* err_msg;
    lua_Debug info;
    int level =  0;
    char lString[1024];
    snprintf(lString, 1024, "require(\"ssfparser\")\nrunSSF(\"%s.ssf\")", fdsprefix);
    luaL_dostring(L, "require \"ssfparser\"");
    int return_code = luaL_loadstring(L, lString);
    switch (return_code) {
      case LUA_OK:
        printf("%s loaded ok\n", filename);
        break;
      case LUA_ERRSYNTAX:
        fprintf(stderr, "Syntax error loading %s\n", filename);
        err_msg = lua_tostring (L, -1);
        fprintf(stderr, "error:%s\n", err_msg);
        level = 0;
        while (lua_getstack(L, level, &info)) {
            lua_getinfo(L, "nSl", &info);
            fprintf(stderr, "  [%d] %s:%d -- %s [%s]\n",
                level, info.short_src, info.currentline,
                (info.name ? info.name : "<unknown>"), info.what);
            ++level;
        }
        break;
      case LUA_ERRMEM:
        break;
      case LUA_ERRGCMM:
        break;
      case LUA_ERRFILE:
        fprintf(stderr, "Could not load file %s\n", filename);
        err_msg = lua_tostring (L, -1);
        fprintf(stderr, "error:%s\n", err_msg);
        level = 0;
        while (lua_getstack(L, level, &info)) {
          lua_getinfo(L, "nSl", &info);
          fprintf(stderr, "  [%d] %s:%d -- %s [%s]\n",
              level, info.short_src, info.currentline,
              (info.name ? info.name : "<unknown>"), info.what);
          ++level;
        }
        break;
    }
    printf("after lua loadfile\n");
}

int yieldOrOkSSF = LUA_YIELD;
int runSSFScript() {
  if (yieldOrOkSSF == LUA_YIELD) {
    printf("running ssf script\n");
    yieldOrOkSSF = lua_resume(L,NULL,0);
    printf("resume done\n");
    if (yieldOrOkSSF == LUA_YIELD) {
      printf("  LUA_YIELD\n");
    } else if (yieldOrOkSSF == LUA_OK) {
      printf("  LUA_OK\n");
    } else if (yieldOrOkSSF == LUA_ERRRUN) {
      printf("  LUA_ERRRUN\n");
      const char *err_msg;
      err_msg = lua_tostring (L, -1);
      fprintf(stderr, "error:%s\n", err_msg);
      lua_Debug info;
      int level = 0;
      while (lua_getstack(L, level, &info)) {
          lua_getinfo(L, "nSl", &info);
          fprintf(stderr, "  [%d] %s:%d -- %s [%s]\n",
              level, info.short_src, info.currentline,
              (info.name ? info.name : "<unknown>"), info.what);
          ++level;
      };
    } else if (yieldOrOkSSF == LUA_ERRMEM) {
      printf("  LUA_ERRMEM\n");
    } else if (yieldOrOkSSF == LUA_ERRGCMM) {
      printf("  LUA_ERRGCMM\n");
    } else {
      printf("  resume code: %i\n", yieldOrOkSSF);
    }
  } else {
    printf("script completed\n");
    lua_close(L);
    glutIdleFunc(NULL);
  }
  return yieldOrOkSSF;
}


int yieldOrOk = LUA_YIELD;
int runLuaScript() {
  if (yieldOrOk == LUA_YIELD) {
    printf("running lua script\n");
    yieldOrOk = lua_resume(L,NULL,0);
    printf("resume done\n");
    if (yieldOrOk == LUA_YIELD) {
      printf("  LUA_YIELD\n");
    } else if (yieldOrOk == LUA_OK) {
      printf("  LUA_OK\n");
    } else if (yieldOrOk == LUA_ERRRUN) {
      printf("  LUA_ERRRUN\n");
      const char *err_msg;
      err_msg = lua_tostring (L, -1);
      fprintf(stderr, "error:%s\n", err_msg);
      lua_Debug info;
      int level = 0;
      while (lua_getstack(L, level, &info)) {
        lua_getinfo(L, "nSl", &info);
        fprintf(stderr, "  [%d] %s:%d -- %s [%s]\n",
            level, info.short_src, info.currentline,
            (info.name ? info.name : "<unknown>"), info.what);
        ++level;
      };
    } else if (yieldOrOk == LUA_ERRMEM) {
      printf("  LUA_ERRMEM\n");
    } else if (yieldOrOk == LUA_ERRGCMM) {
      printf("  LUA_ERRGCMM\n");
    } else {
      printf("  resume code: %i\n", yieldOrOk);
    }
  } else {
    printf("script completed\n");
    lua_close(L);
    glutIdleFunc(NULL);
  }
  return yieldOrOk;
}
