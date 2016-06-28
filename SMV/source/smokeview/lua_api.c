
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
#include "gd.h"

#include <unistd.h>

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

int lua_render_var(lua_State *L) {
  gdImagePtr RENDERimage;
  int return_code;
  char *imageData;
  int imageSize;

  // render image to RENDERimage gd buffer
  return_code = RenderFrameLuaVar(VIEW_CENTER, &RENDERimage);
  // convert to a simpler byte-buffer
  imageData = gdImagePngPtr(RENDERimage, &imageSize);
  // push to stack
  lua_pushlstring(L, imageData, imageSize);
  // destroy C copy
  gdImageDestroy(RENDERimage);

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


int lua_get_clipping_mode(lua_State *L) {
  lua_pushnumber(L, get_clipping_mode());
  return 1;
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



//////////////////////

int lua_set_colorbar_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_colorbar_visibility(setting);
  return 0;
}

int lua_get_colorbar_visibility(lua_State *L) {
  int setting = get_colorbar_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_colorbar_visibility(lua_State *L) {
  toggle_colorbar_visibility();
  return 0;
}

int lua_set_timebar_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_timebar_visibility(setting);
  return 0;
}

int lua_get_timebar_visibility(lua_State *L) {
  int setting = get_timebar_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_timebar_visibility(lua_State *L) {
  toggle_timebar_visibility();
  return 0;
}

// title
int lua_set_title_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_title_visibility(setting);
  return 0;
}

int lua_get_title_visibility(lua_State *L) {
  int setting = get_title_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_title_visibility(lua_State *L) {
  toggle_title_visibility();
  return 0;
}

// axis
int lua_set_axis_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_axis_visibility(setting);
  return 0;
}

int lua_get_axis_visibility(lua_State *L) {
  int setting = get_axis_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_axis_visibility(lua_State *L) {
  toggle_axis_visibility();
  return 0;
}

// frame
int lua_set_framelabel_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_framelabel_visibility(setting);
  return 0;
}

int lua_get_framelabel_visibility(lua_State *L) {
  int setting = get_framelabel_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_framelabel_visibility(lua_State *L) {
  toggle_framelabel_visibility();
  return 0;
}

// framerate
int lua_set_framerate_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_framerate_visibility(setting);
  return 0;
}

int lua_get_framerate_visibility(lua_State *L) {
  int setting = get_framerate_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_framerate_visibility(lua_State *L) {
  toggle_framerate_visibility();
  return 0;
}

// grid locations
int lua_set_gridloc_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_gridloc_visibility(setting);
  return 0;
}

int lua_get_gridloc_visibility(lua_State *L) {
  int setting = get_gridloc_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_gridloc_visibility(lua_State *L) {
  toggle_gridloc_visibility();
  return 0;
}

// hrrpuv cutoff
int lua_set_hrrcutoff_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_hrrcutoff_visibility(setting);
  return 0;
}

int lua_get_hrrcutoff_visibility(lua_State *L) {
  int setting = get_hrrcutoff_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_hrrcutoff_visibility(lua_State *L) {
  toggle_hrrcutoff_visibility();
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
// memory load
#ifdef pp_memstatus
int lua_set_memload_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_memload_visibility(setting);
  return 0;
}

int lua_get_memload_visibility(lua_State *L) {
  int setting = get_memload_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_memload_visibility(lua_State *L) {
  toggle_memload_visibility();
  return 0;
}
#endif

// mesh label
int lua_set_meshlabel_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_meshlabel_visibility(setting);
  return 0;
}

int lua_get_meshlabel_visibility(lua_State *L) {
  int setting = get_meshlabel_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_meshlabel_visibility(lua_State *L) {
  toggle_meshlabel_visibility();
  return 0;
}

// slice average
int lua_set_slice_average_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_slice_average_visibility(setting);
  return 0;
}

int lua_get_slice_average_visibility(lua_State *L) {
  int setting = get_slice_average_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_slice_average_visibility(lua_State *L) {
  toggle_slice_average_visibility();
  return 0;
}

// time
int lua_set_time_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_time_visibility(setting);
  return 0;
}

int lua_get_time_visibility(lua_State *L) {
  int setting = get_time_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_time_visibility(lua_State *L) {
  toggle_time_visibility();
  return 0;
}

// user settable ticks
int lua_set_user_ticks_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_user_ticks_visibility(setting);
  return 0;
}

int lua_get_user_ticks_visibility(lua_State *L) {
  int setting = get_user_ticks_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_user_ticks_visibility(lua_State *L) {
  toggle_user_ticks_visibility();
  return 0;
}

// version info
int lua_set_version_info_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_version_info_visibility(setting);
  return 0;
}

int lua_get_version_info_visibility(lua_State *L) {
  int setting = get_version_info_visibility();
  lua_pushboolean(L, setting);
  return 1;
}

int lua_toggle_version_info_visibility(lua_State *L) {
  toggle_version_info_visibility();
  return 0;
}

// set all
int lua_set_all_label_visibility(lua_State *L) {
  int setting = lua_toboolean(L, 1);
  set_all_label_visibility(setting);
  return 0;
}

//////////////////////////////////////

int lua_blockage_view_method(lua_State *L) {
  int setting = lua_tonumber(L, 1);
  int return_code = blockage_view_method(setting);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_blockage_outline_color(lua_State *L) {
  int setting = lua_tonumber(L, 1);
  int return_code = blockage_outline_color(setting);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_blockage_locations(lua_State *L) {
  int setting = lua_tonumber(L, 1);
  int return_code = blockage_locations(setting);
  lua_pushnumber(L, return_code);
  return 1;
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

int lua_set_ambientlight(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_ambientlight(r, g, b);
  return 0;
}

int lua_set_backgroundcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_backgroundcolor(r, g, b);
  return 0;
}

int lua_set_blockcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_blockcolor(r, g, b);
  return 0;
}

int lua_set_blockshininess(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_blockshininess(v);
  return 0;
}

int lua_set_blockspecular(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_blockspecular(r, g, b);
  return 0;
}

int lua_set_boundcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_boundcolor(r, g, b);
  return 0;
}

int lua_set_colorbar_textureflag(lua_State *L) {
  int setting = lua_tonumber(L, 1);
  int return_code = set_colorbar_textureflag(setting);
  lua_pushnumber(L, 1);
  return 1;
}

float getcolorfield(lua_State *L, int stack_index,const char *key) {
  if (!lua_istable(L, stack_index)){
    fprintf(stderr, "stack is not a table at index, cannot use getcolorfield\n");
    exit(1);
  }
  // if stack index is relative (negative) convert to absolute (positive)
  if (stack_index < 0) {
    stack_index = lua_gettop(L) + stack_index + 1;
  }
  lua_pushstring(L, key);
  lua_gettable(L, stack_index);
  float result = lua_tonumber(L, -1);
  lua_pop(L, 1);
  return result;
}

int get_color(lua_State *L, int stack_index, float *color) {
  if (!lua_istable(L, stack_index)){
    fprintf(stderr, "color table is not present\n");
    return 1;
  }
  float r = getcolorfield(L, stack_index, "r");
  float g = getcolorfield(L, stack_index, "g");
  float b = getcolorfield(L, stack_index, "b");
  color[0] = r;
  color[1] = g;
  color[2] = b;
  return 0;
}
int lua_set_colorbar_colors(lua_State *L) {
  printf("running: lua_set_colorbar_colors\n");
  if (!lua_istable(L, 1)){
    fprintf(stderr, "colorbar table is not present\n");
    return 1;
  }
  int ncolors = 0;
  lua_pushnil(L);
  while(lua_next(L,1)!=0) {
    ncolors++;
    lua_pop(L, 1);
  }
  int i;
  float *color;
  float colors[ncolors][3];
  for (i = 1; i <= ncolors; i++) {
    lua_pushnumber(L, i);
    lua_gettable(L, 1);
    get_color(L, -1, colors[i-1]);
  }

  int return_code = set_colorbar_colors(ncolors, colors);
  return 0;
}

int lua_get_colorbar_colors(lua_State *L) {
  int i;
  float *rgb_ini_copy_p = rgb_ini;
  lua_createtable(L, 0, nrgb_ini);
  for (i = 0; i < nrgb_ini; i++) {
    lua_pushnumber(L, i+1);
    lua_createtable(L, 0, 2);

    lua_pushnumber(L, *rgb_ini_copy_p);
    lua_setfield(L, -2, "r");

    lua_pushnumber(L, *(rgb_ini_copy_p + 1));
    lua_setfield(L, -2, "g");

    lua_pushnumber(L, *(rgb_ini_copy_p + 2));
    lua_setfield(L, -2, "b");

    lua_settable(L, -3);
    rgb_ini_copy_p += 3;
  }
  PRINTF("lua: done creating colorbar table\n");
  // Leaves one returned value on the stack, the mesh table.
  return 1;
}

int lua_set_color2bar_colors(lua_State *L) {
  int ncolors = lua_tonumber(L, 1);
  if (!lua_istable(L, -1)){
    fprintf(stderr, "colorbar table is not present\n");
    return 1;
  }
  int i;
  float *color;
  float colors[ncolors][3];
  for (i = 1; i <= ncolors; i++) {
    lua_pushnumber(L, i);
    lua_gettable(L,-2);
    get_color(L, -1, colors[i-1]);
  }

  int return_code = set_color2bar_colors(ncolors, colors);
  return 0;
}

int lua_get_color2bar_colors(lua_State *L) {
  int i;
  float *rgb_ini_copy_p = rgb2_ini;
  lua_createtable(L, 0, nrgb2_ini);
  for (i = 0; i < nrgb2_ini; i++) {
    lua_pushnumber(L, i+1);
    lua_createtable(L, 0, 2);

    lua_pushnumber(L, *rgb_ini_copy_p);
    lua_setfield(L, -2, "r");

    lua_pushnumber(L, *(rgb_ini_copy_p + 1));
    lua_setfield(L, -2, "g");

    lua_pushnumber(L, *(rgb_ini_copy_p + 2));
    lua_setfield(L, -2, "b");

    lua_settable(L, -3);
    rgb_ini_copy_p += 3;
  }
  PRINTF("lua: done creating color2bar table\n");
  // Leaves one returned value on the stack, the mesh table.
  return 1;
}

int lua_set_diffuselight(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_diffuselight(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_directioncolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_directioncolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_flip(lua_State  *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_flip(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_foregroundcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_foregroundcolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_heatoffcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_heatoffcolor(r, g, b);
  return 0;
}

int lua_set_heatoncolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_heatoncolor(r, g, b);
  return 0;
}

int lua_set_isocolors(lua_State *L) {
  float shininess = lua_tonumber(L, 1);
  float default_opaqueness = lua_tonumber(L, 2);
  float specular[3];
  get_color(L, 3, specular);
  int ncolors = 0;
  int x;
  // count the number of colours
  lua_pushnil(L);  /* first key */
  while (lua_next(L, 4) != 0) {
    lua_pop(L, 1); // remove value (leave key for next iteration)
    ncolors++;
    x = lua_tonumber(L, -1);
  }
  int i;
  float *color;
  float colors[ncolors][3];
  for (i = 1; i <= ncolors; i++) {
    if (!lua_istable(L, 4)) {
      fprintf(stderr, "isocolor table is not present\n");
      return 1;
    }
    lua_pushnumber(L, i);
    lua_gettable(L, 4);
    get_color(L, -1, colors[i-1]);
  }
  // for (i = 0; i < ncolors; i++) {
  //   printf("%d: %f %f %f\n", i,
  //     colors[i][0], colors[i][1], colors[i][2]);
  // }
  // specular = lua_tonumber(L, 3);
  // int return_code = set_diffuselight(r, g, b);
  return 0;
}

int lua_set_colortable(lua_State *L) {
  // int ncolors = lua_tonumber(L, 1);
  int ncolors = 0;
  int i = 0;
  // count the number of colours
  lua_pushnil(L);  /* first key */
  while (lua_next(L, 1) != 0) {
    lua_pop(L, 1); // remove value (leave key for next iteration)
    ncolors++;
  }
  // initialise arrays using the above count info
  float colors[ncolors][3];
  char names[ncolors][255];
  /* table is in the stack at index 't' */
  lua_pushnil(L);  /* first key */
  while (lua_next(L, 1) != 0) {
    /* uses 'key' (at index -2) and 'value' (at index -1) */
    strncpy(names[i], lua_tostring(L, -2), 255);
    get_color(L, -1, colors[i]);
    /* removes 'value'; keeps 'key' for next iteration */
    lua_pop(L, 1);
    i++;
  }
  return 0;
}

int lua_set_light0(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_light0(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_light1(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_light1(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_lightpos0(lua_State *L) {
  float a = lua_tonumber(L, 1);
  float b = lua_tonumber(L, 2);
  float c = lua_tonumber(L, 3);
  float d = lua_tonumber(L, 4);
  int return_code = set_lightpos0(a, b, c, d);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_lightpos1(lua_State *L) {
  float a = lua_tonumber(L, 1);
  float b = lua_tonumber(L, 2);
  float c = lua_tonumber(L, 3);
  float d = lua_tonumber(L, 4);
  int return_code = set_lightpos1(a, b, c, d);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_lightmodellocalviewer(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_lightmodellocalviewer(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_lightmodelseparatespecularcolor(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_lightmodelseparatespecularcolor(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sensorcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_sensorcolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sensornormcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_sensornormcolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_bw(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int return_code = set_bw(a, b);
  lua_pushnumber(L, return_code);
  return 1;
}


int lua_set_sprinkleroffcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_sprinkleroffcolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sprinkleroncolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_sprinkleroncolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_staticpartcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_staticpartcolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_timebarcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_timebarcolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_ventcolor(lua_State *L) {
  float r = lua_tonumber(L, 1);
  float g = lua_tonumber(L, 2);
  float b = lua_tonumber(L, 3);
  int return_code = set_ventcolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_gridlinewidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_gridlinewidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_isolinewidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_isolinewidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_isopointsize(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_isopointsize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_linewidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_linewidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_partpointsize(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_partpointsize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_plot3dlinewidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_plot3dlinewidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_plot3dpointsize(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_plot3dpointsize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sensorabssize(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_sensorabssize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sensorrelsize(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_sensorrelsize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sliceoffset(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_sliceoffset(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_smoothlines(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_smoothlines(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_spheresegs(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_spheresegs(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sprinklerabssize(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_sprinklerabssize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_streaklinewidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_streaklinewidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_ticklinewidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_ticklinewidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_usenewdrawface(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_usenewdrawface(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_veccontours(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_veccontours(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_veclength(lua_State *L) {
  float a = lua_tonumber(L, 1);
  float b = lua_tonumber(L, 2);
  float c = lua_tonumber(L, 3);
  int return_code = set_veclength(a, b, c);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_vectorlinewidth(lua_State *L) {
  float a = lua_tonumber(L, 1);
  float b = lua_tonumber(L, 2);
  int return_code = set_vectorlinewidth(a, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_vectorpointsize(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_vectorpointsize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_ventlinewidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_ventlinewidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_ventoffset(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_ventoffset(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_windowoffset(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_windowoffset(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_windowwidth(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_windowwidth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_windowheight(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_windowheight(v);
  lua_pushnumber(L, return_code);
  return 1;
}

// --  *** DATA LOADING ***

int lua_set_boundzipstep(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_boundzipstep(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_fed(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_fed(v);
  lua_pushnumber(L, return_code);
  return 1;
}


int lua_set_fedcolorbar(lua_State *L) {
  const char *name = lua_tostring(L, 1);
  int return_code = set_fedcolorbar(name);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_isozipstep(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_isozipstep(v);
  lua_pushnumber(L, return_code);
  return 1;
}


int lua_set_nopart(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_nopart(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showfedarea(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showfedarea(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sliceaverage(lua_State *L) {
  int flag = lua_tonumber(L, 1);
  float interval = lua_tonumber(L, 2);
  int vis = lua_tonumber(L, 3);
  int return_code = set_sliceaverage(flag, interval, vis);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_slicedataout(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_slicedataout(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_slicezipstep(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_slicezipstep(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_smoke3dzipstep(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_smoke3dzipstep(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_userrotate(lua_State *L) {
  int index = lua_tonumber(L, 1);
  int show_center = lua_tonumber(L, 2);
  float x = lua_tonumber(L, 3);
  float y = lua_tonumber(L, 4);
  float z = lua_tonumber(L, 5);
  int return_code = set_userrotate(index, show_center, x, y, z);
  lua_pushnumber(L, return_code);
  return 1;
}

// --  *** VIEW PARAMETERS ***
int lua_set_aperture(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_aperture(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_axissmooth(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_axissmooth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_blocklocation(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_blocklocation(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_boundarytwoside(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_boundarytwoside(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_clip(lua_State *L) {
  float v_near = lua_tonumber(L, 1);
  float v_far = lua_tonumber(L, 2);
  int return_code = set_clip(v_near, v_far);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_contourtype(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_contourtype(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_cullfaces(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_cullfaces(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_texturelighting(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_texturelighting(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_eyeview(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_eyeview(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_eyex(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_eyex(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_eyey(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_eyey(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_eyez(lua_State *L) {
  float v = lua_tonumber(L, 1);
  int return_code = set_eyez(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_fontsize(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_fontsize(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_frameratevalue(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_frameratevalue(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_geomdiags(lua_State *L) {
  int structured = lua_tonumber(L, 1);
  int unstructured = lua_tonumber(L, 2);
  int diagnostics = lua_tonumber(L, 3);
  int return_code = set_geomdiags(structured, unstructured, diagnostics);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showfaces_interior(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showfaces_interior(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showfaces_exterior(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showfaces_exterior(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showfaces_solid(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showfaces_solid(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showfaces_outline(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showfaces_outline(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_smoothgeomnormal(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_smoothgeomnormal(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showvolumes_interior(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showvolumes_interior(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showvolumes_exterior(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showvolumes_exterior(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showvolumes_solid(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showvolumes_solid(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showvolumes_outline(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showvolumes_outline(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_geomvertexag(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_geomvertexag(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_geommaxangle(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_geommaxangle(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_gversion(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_gversion(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_isotran2(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_isotran2(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_meshvis(lua_State *L) {
  int n = 0;
  int i = 0;
  // count the number of values
  lua_pushnil(L);
  while (lua_next(L, 1) != 0) {
    lua_pop(L, 1); // remove value (leave key for next iteration)
    n++;
  }
  // initialise arrays using the above count info
  int vals[n];
  /* table is in the stack at index 't' */
  lua_pushnil(L);  /* first key */
  while (lua_next(L, 1) != 0) {
    vals[i] = lua_tonumber(L, -2);
    /* removes 'value'; keeps 'key' for next iteration */
    lua_pop(L, 1);
    i++;
  }
  int return_code = set_meshvis(n, vals);
  return 0;
}

int lua_set_meshoffset(lua_State *L) {
  int meshnum = lua_tonumber(L, 1);
  int value = lua_tonumber(L, 2);
  int return_code = set_meshoffset(meshnum, value);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_northangle(lua_State *L) {
  int vis = lua_tonumber(L, 1);
  float x = lua_tonumber(L, 2);
  float y = lua_tonumber(L, 3);
  float z = lua_tonumber(L, 4);
  int return_code = set_northangle(vis, x, y, z);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_offsetslice(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_offsetslice(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_outlinemode(lua_State *L) {
  int highlight = lua_tonumber(L, 1);
  int outline = lua_tonumber(L, 2);
  int return_code = set_outlinemode(highlight, outline);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_p3dsurfacetype(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_p3dsurfacetype(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_p3dsurfacesmooth(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_p3dsurfacesmooth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_projection(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_projection(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_scaledfont(lua_State *L) {
  int height2d = lua_tonumber(L, 1);
  int height2dwidth = lua_tonumber(L, 2);
  int thickness2d = lua_tonumber(L, 3);
  int height3d = lua_tonumber(L, 3);
  int height3dwidth = lua_tonumber(L, 5);
  int thickness3d = lua_tonumber(L, 6);
  int return_code = set_scaledfont(height2d, height2dwidth, thickness2d,
                                   height3d, height3dwidth, thickness3d);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showalltextures(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showalltextures(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showaxislabels(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showaxislabels(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showblocklabel(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showblocklabel(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showblocks(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showblocks(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showcadandgrid(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showcadandgrid(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showcadopaque(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showcadopaque(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showceiling(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showceiling(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showcolorbars(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showcolorbars(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showcvents(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 1);
  int return_code = set_showcvents(a, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showdummyvents(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showdummyvents(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showevacslices(lua_State *L) {
  int show_slices = lua_tonumber(L, 1);
  int constant_coloring = lua_tonumber(L, 2);
  int show_colorbar = lua_tonumber(L, 3);
  int return_code = set_showevacslices(show_slices, constant_coloring,
                                       show_colorbar);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showfloor(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showfloor(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showframe(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showframe(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showframelabel(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showframelabel(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showframerate(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showframerate(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showgrid(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showgrid(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showgridloc(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showgridloc(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showhmstimelabel(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showhmstimelabel(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showhrrcutoff(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showhrrcutoff(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showiso(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showiso(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showisonormals(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showisonormals(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showlabels(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showlabels(v);
  lua_pushnumber(L, return_code);
  return 1;
}

#ifdef pp_memstatus
int lua_set_showmemload(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showmemload(v);
  lua_pushnumber(L, return_code);
  return 1;
}
#endif

int lua_set_showopenvents(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 1);
  int return_code = set_showopenvents(a, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showothervents(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showothervents(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showsensors(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int return_code = set_showsensors(a, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showsliceinobst(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showsliceinobst(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showsmokepart(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showsmokepart(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showsprinkpart(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showsprinkpart(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showstreak(lua_State *L) {
  int show = lua_tonumber(L, 1);
  int step = lua_tonumber(L, 2);
  int showhead = lua_tonumber(L, 3);
  int index = lua_tonumber(L, 4);
  int return_code = set_showstreak(show, step, showhead, index);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showterrain(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showterrain(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtetras(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int return_code = set_showtetras(a, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showthreshold(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  float c = lua_tonumber(L, 3);
  int return_code = set_showthreshold(a, b, c);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showticks(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showticks(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtimebar(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtimebar(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtimelabel(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtimelabel(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtitle(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtitle(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtracersalways(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtracersalways(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtriangles(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int c = lua_tonumber(L, 3);
  int d = lua_tonumber(L, 4);
  int e = lua_tonumber(L, 5);
  int f = lua_tonumber(L, 6);
  int return_code = set_showtriangles(a, b, c, d, e, f);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtransparent(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtransparent(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtranparentvents(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtransparentvents(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtrianglecount(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtrianglecount(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showventflow(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int c = lua_tonumber(L, 3);
  int d = lua_tonumber(L, 4);
  int e = lua_tonumber(L, 5);
  int return_code = set_showventflow(a, b, c, d, e);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showvents(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showvents(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showwalls(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showwalls(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_skipembedslice(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_skipembedslice(v);
  lua_pushnumber(L, return_code);
  return 1;
}

#ifdef pp_SLICEUP
int lua_set_slicedup(lua_State *L) {
  int scalar = lua_tonumber(L, 1);
  int vector = lua_tonumber(L, 1);
  int return_code = set_slicedup(scalar, vector);
  lua_pushnumber(L, return_code);
  return 1;
}
#endif

int lua_set_smokesensors(lua_State *L) {
  int show = lua_tonumber(L, 1);
  int test = lua_tonumber(L, 2);
  int return_code = set_smokesensors(show, test);
  lua_pushnumber(L, return_code);
  return 1;
}


// int set_smoothblocksolid(int v); // SMOOTHBLOCKSOLID
#ifdef pp_LANG
int lua_set_startuplang(lua_State *L) {
  const char *lang = lua_tostring(L, 1);
  int return_code = set_startuplang(lang);
  lua_pushnumber(L, return_code);
  return 1;
}
#endif

int lua_set_stereo(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_stereo(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_surfinc(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_surfinc(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_terrainparams(lua_State *L) {
  int r_min = lua_tonumber(L, 1);
  int g_min = lua_tonumber(L, 2);
  int b_min = lua_tonumber(L, 3);
  int r_max = lua_tonumber(L, 4);
  int g_max = lua_tonumber(L, 5);
  int b_max = lua_tonumber(L, 6);
  int vert_factor = lua_tonumber(L, 7);
  int return_code = set_terrainparams(r_min, g_min, b_min, r_max, g_max, b_max,
                                      vert_factor);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_titlesafe(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_titlesafe(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_trainermode(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_trainermode(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_trainerview(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_trainerview(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_transparent(lua_State *L) {
  int use_flag = lua_tonumber(L, 1);
  float level = lua_tonumber(L, 2);
  int return_code = set_transparent(use_flag, level);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_treeparms(lua_State *L) {
  int minsize = lua_tonumber(L, 1);
  int visx = lua_tonumber(L, 2);
  int visy = lua_tonumber(L, 3);
  int visz = lua_tonumber(L, 4);
  int return_code = set_treeparms(minsize, visx, visy, visz);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_twosidedvents(lua_State *L) {
  int internal = lua_tonumber(L, 1);
  int external = lua_tonumber(L, 2);
  int return_code = set_twosidedvents(internal, external);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_vectorskip(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_vectorskip(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_volsmoke(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int c = lua_tonumber(L, 3);
  int d = lua_tonumber(L, 4);
  int e = lua_tonumber(L, 5);
  float f = lua_tonumber(L, 6);
  float g = lua_tonumber(L, 7);
  float h = lua_tonumber(L, 8);
  float i = lua_tonumber(L, 9);
  float j = lua_tonumber(L, 10);
  float k = lua_tonumber(L, 11);
  float l = lua_tonumber(L, 12);
  int return_code = set_volsmoke(a, b, c, d, e, f, g, h, i, j, k, l);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_zoom(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int return_code = set_zoom(a, b);
  lua_pushnumber(L, return_code);
  return 1;
}

// *** MISC ***
int lua_set_cellcentertext(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_cellcentertext(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_inputfile(lua_State *L) {
  const char *inputfile = lua_tostring(L, 1);
  int return_code = set_inputfile(inputfile);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_labelstartupview(lua_State *L) {
  const char *viewname = lua_tostring(L, 1);
  int return_code = set_labelstartupview(viewname);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_pixelskip(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_pixelskip(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_renderclip(lua_State *L) {
  int use_flag = lua_tonumber(L, 1);
  int left = lua_tonumber(L, 2);
  int right = lua_tonumber(L, 3);
  int bottom = lua_tonumber(L, 4);
  int top = lua_tonumber(L, 5);
  int return_code = set_renderclip(use_flag, left, right, bottom, top);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_renderfilelabel(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_renderfilelabel(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_renderfiletype(lua_State *L) {
  int render = lua_tonumber(L, 1);
  int movie = lua_tonumber(L, 2);
  int return_code = set_renderfiletype(render, movie);
  lua_pushnumber(L, return_code);
  return 1;
}


// int lua_set_skybox(lua_State *L){
//   return 0;
// }

int lua_set_renderoption(lua_State *L) {
  int opt = lua_tonumber(L, 1);
  int rows = lua_tonumber(L, 1);
  int return_code = set_renderoption(opt, rows);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_unitclasses(lua_State *L) {
  int i = 0;
  int n = 0;
  if (!lua_istable(L, -1)){
    fprintf(stderr, "stack is not a table at index\n");
    exit(1);
  }
  lua_pushnil(L);
  while(lua_next(L,-2)!=0) {
    lua_pop(L,1);
    n++;
  }
  int indices[n];
  lua_pushnil(L);
  while(lua_next(L,-2)!=0) {
    indices[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
    i++;
  }
  int return_code = set_unitclasses(n, indices);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_zaxisangles(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int c = lua_tonumber(L, 3);
  int return_code = set_zaxisangles(a, b, c);
  lua_pushnumber(L, return_code);
  return 1;
}

// *** 3D SMOKE INFO ***
int lua_set_adjustalpha(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_adjustalpha(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_colorbartype(lua_State *L) {
  int type = lua_tonumber(L, 1);
  const char *label = lua_tostring(L, 2);
  int return_code = set_colorbartype(type, label);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_extremecolors(lua_State *L) {
  int rmin = lua_tonumber(L, 1);
  int gmin = lua_tonumber(L, 2);
  int bmin = lua_tonumber(L, 3);
  int rmax = lua_tonumber(L, 4);
  int gmax = lua_tonumber(L, 5);
  int bmax = lua_tonumber(L, 6);
  int return_code = set_extremecolors(rmin, gmin, bmin, rmax, gmax, bmax);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_firecolor(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_firecolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_firecolormap(lua_State *L) {
  int type = lua_tonumber(L, 1);
  int index = lua_tonumber(L, 2);
  int return_code = set_firecolormap(type, index);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_firedepth(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_firedepth(v);
  lua_pushnumber(L, return_code);
  return 1;
}

// int set_gcolorbar(int ncolorbarini, ) {
//   colorbardata *cbi;
//   int r1, g1, b1;
//   int n;

//   initdefaultcolorbars();

//   ncolorbars = ndefaultcolorbars + ncolorbarini;
//   if(ncolorbarini>0)ResizeMemory((void **)&colorbarinfo, ncolorbars*sizeof(colorbardata));

//   for(n = ndefaultcolorbars; n<ncolorbars; n++){
//     char *cb_buffptr;

//     cbi = colorbarinfo + n;
//     fgets(buffer, 255, stream);
//     trim_back(buffer);
//     cb_buffptr = trim_front(buffer);
//     strcpy(cbi->label, cb_buffptr);

//     fgets(buffer, 255, stream);
//     sscanf(buffer, "%i %i", &cbi->nnodes, &cbi->nodehilight);
//     if(cbi->nnodes<0)cbi->nnodes = 0;
//     if(cbi->nodehilight<0 || cbi->nodehilight >= cbi->nnodes){
//       cbi->nodehilight = 0;
//     }

//     cbi->label_ptr = cbi->label;
//     for(i = 0; i<cbi->nnodes; i++){
//       int icbar;
//       int nn;

//       fgets(buffer, 255, stream);
//       r1 = -1; g1 = -1; b1 = -1;
//       sscanf(buffer, "%i %i %i %i", &icbar, &r1, &g1, &b1);
//       cbi->index_node[i] = icbar;
//       nn = 3 * i;
//       cbi->rgb_node[nn] = r1;
//       cbi->rgb_node[nn + 1] = g1;
//       cbi->rgb_node[nn + 2] = b1;
//     }
//     remapcolorbar(cbi);
//   }
//   return 0;
// } // GCOLORBAR

int lua_set_showextremedata(lua_State *L) {
  int show_extremedata = lua_tonumber(L, 1);
  int below = lua_tonumber(L, 2);
  int above = lua_tonumber(L, 3);
  int return_code = set_showextremedata(show_extremedata, below, above);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_smokecolor(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_smokecolor(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_smokecull(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_smokecull(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_smokeskip(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_smokeskip(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_smokealbedo(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_smokealbedo(v);
  lua_pushnumber(L, return_code);
  return 1;
}

#ifdef pp_GPU
int lua_set_smokerthick(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_smokerthick(v);
  lua_pushnumber(L, return_code);
  return 1;
}
#endif

int lua_set_smokethick(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_smokethick(v);
  lua_pushnumber(L, return_code);
  return 1;
}

#ifdef pp_GPU
int lua_set_usegpu(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_usegpu(v);
  lua_pushnumber(L, return_code);
  return 1;
}
#endif

// *** ZONE FIRE PARAMETRES ***
int lua_set_showhazardcolors(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showhazardcolors(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showhzone(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showhzone(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showszone(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showszone(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showvzone(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showvzone(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showzonefire(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showzonefire(v);
  lua_pushnumber(L, return_code);
  return 1;
}

// *** TOUR INFO ***
int lua_set_showpathnodes(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showpathnodes(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_showtourroute(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showtourroute(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_tourcolors_selectedpathline(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_tourcolors_selectedpathline(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}
int lua_set_tourcolors_selectedpathlineknots(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_tourcolors_selectedpathlineknots(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}
int lua_set_tourcolors_selectedknot(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_tourcolors_selectedknot(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}
int lua_set_tourcolors_pathline(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_tourcolors_selectedpathline(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}
int lua_set_tourcolors_pathknots(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_tourcolors_pathknots(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}
int lua_set_tourcolors_text(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_tourcolors_text(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}
int lua_set_tourcolors_avatar(lua_State *L) {
  int r = lua_tonumber(L, 1);
  int g = lua_tonumber(L, 2);
  int b = lua_tonumber(L, 3);
  int return_code = set_tourcolors_avatar(r, g, b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_tourconstantvel(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_tourconstantvel(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_viewalltours(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_viewalltours(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_viewtimes(lua_State *L) {
  float start = lua_tonumber(L, 1);
  float stop = lua_tonumber(L, 2);
  int ntimes = lua_tonumber(L, 3);
  int return_code = set_viewtimes(start, stop, ntimes);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_viewtourfrompath(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_viewtourfrompath(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_avatarevac(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_avatarevac(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int *lua_get_int_array(lua_State *L, int snumber) {
  // count the length of vals
  int nvals = 0;
  lua_pushnil(L);
  while(lua_next(L, snumber)!=0) {
    nvals++;
    lua_pop(L, 1);
  }
  // create array of vals
  int *vals;
  vals = (int *) calloc(nvals, sizeof(int));
  int i = 0;
  lua_pushnil(L);
  while(lua_next(L, snumber)!=0) {
    vals[i] = lua_tonumber(L, -1);
    i++;
    lua_pop(L,1);
  }
  return vals;
}

float *lua_get_float_array(lua_State *L, int snumber) {
  // count the length of vals
  int nvals = 0;
  lua_pushnil(L);
  while(lua_next(L, snumber)!=0) {
    nvals++;
    lua_pop(L, 1);
  }
  // create array of vals
  float *vals;
  vals = (float *) calloc(nvals, sizeof(float));
  int i = 0;
  lua_pushnil(L);
  while(lua_next(L, snumber)!=0) {
    vals[i] = lua_tonumber(L, -1);
    i++;
    lua_pop(L,1);
  }
  return vals;
}


int lua_set_geometrytest(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  float c = lua_tonumber(L, 3);
  float d = lua_tonumber(L, 4);
  // count the length of vals
  int *vals;
  float *b1Vals;
  float *b2Vals;
  float *b3Vals;
  // these arrays must be freed
  vals = lua_get_int_array(L, 5);
  b1Vals = lua_get_float_array(L, 6);
  b2Vals = lua_get_float_array(L, 7);
  b3Vals = lua_get_float_array(L, 8);

  int return_code = set_geometrytest(a, b, c, d, vals, b1Vals, b2Vals, b3Vals);
  free(vals);
  free(b1Vals);
  free(b2Vals);
  free(b3Vals);
  // int set_geometrytest(int a, int b, int c, int d, int vals[],
  //                    float b1Vals[], float b2Vals[], float b3Vals[]); // GEOMETRYTEST
  lua_pushnumber(L, return_code);
  return 1;
}

// int set_geometrytest(int a, int b, int c, int d, int vals[],
//                      float b1Vals[], float b2Vals[], float b3Vals[]) {
//   int *v;
//   int ii;
//   geomtest_option = a;
//   show_tetratest_labels = b;
//   tetra_line_thickness = c;
//   tetra_point_size = d;
//   v = tetrabox_vis;
//   ONEORZERO(show_tetratest_labels);
//   for(ii = 0; ii<10; ii++){
//     v[ii] = vals[ii];
//     ONEORZERO(v[ii]);
//   }
//   for(ii = 0; ii<6; ii++){
//     box_bounds2[ii] = b1Vals[ii];
//   }
//   for(ii = 0; ii<12; ii++){
//      tetra_vertices[ii] = b2Vals[ii];
//   }
//   for(ii = 0; ii<3; ii++){
//     box_translate[ii] = b3Vals[ii];
//   }
//   return 0;
// } //  GEOMETRYTEST

int lua_set_devicevectordimensions(lua_State *L) {
  float baselength = lua_tonumber(L, 1);
  float basediameter = lua_tonumber(L, 2);
  float headlength = lua_tonumber(L, 3);
  float headdiameter = lua_tonumber(L, 4);
  int return_code = set_devicevectordimensions(baselength, basediameter,
                                               headlength, headdiameter);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_devicebounds(lua_State *L) {
  float min = lua_tonumber(L, 1);
  float max = lua_tonumber(L, 2);
  int return_code = set_devicebounds(min, max);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_deviceorientation(lua_State *L) {
  int a = lua_tonumber(L, 1);
  float b = lua_tonumber(L, 2);
  int return_code = set_deviceorientation(a,b);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_gridparms(lua_State *L) {
  int vx = lua_tonumber(L, 1);
  int vy = lua_tonumber(L, 2);
  int vz = lua_tonumber(L, 3);
  int px = lua_tonumber(L, 4);
  int py = lua_tonumber(L, 5);
  int pz = lua_tonumber(L, 6);
  int return_code = set_gridparms(vx, vy, vz, px, py, pz);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_gsliceparms(lua_State *L) {
  int i;
  int vis_data = lua_tonumber(L, 1);
  int vis_triangles = lua_tonumber(L, 2);
  int vis_triangulation = lua_tonumber(L, 3);
  int vis_normal = lua_tonumber(L, 4);
  float xyz[3];
  // TODO: use named fields (e.g. xyz)
  for (i = 0; i < 3; i++) {
    lua_pushnumber(L, i);
    lua_gettable(L, 5);
    xyz[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
    i++;
  }
  float azelev[2];
  for (i = 0; i < 2; i++) {
    lua_pushnumber(L, i);
    lua_gettable(L, 6);
    azelev[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
    i++;
  }
  int return_code = set_gsliceparms(vis_data, vis_triangles, vis_triangulation,
                                    vis_normal, xyz, azelev);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_loadfilesatstartup(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_loadfilesatstartup(v);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_mscale(lua_State *L) {
  float a = lua_tonumber(L, 1);
  float b = lua_tonumber(L, 2);
  float c = lua_tonumber(L, 3);
  int return_code = set_mscale(a,b,c);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_sliceauto(lua_State *L) {
  lua_pushnil(L);
  int n = 0;
  while(lua_next(L, -2) != 0) {
    lua_pop(L,1);
    n++;
  }
  int i = 0;
  int vals[n];
  lua_pushnil(L);
  while(lua_next(L, -2) != 0) {
    vals[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
    i++;
  }
  int return_code = set_sliceauto(n, vals);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_msliceauto(lua_State *L) {
  lua_pushnil(L);
  int n = 0;
  while(lua_next(L, -2) != 0) {
    lua_pop(L,1);
    n++;
  }
  int i = 0;
  int vals[n];
  lua_pushnil(L);
  while(lua_next(L, -2) != 0) {
    vals[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
    i++;
  }
  int return_code = set_msliceauto(n, vals);
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_set_compressauto(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_compressauto(v);
  lua_pushnumber(L, return_code);
  return 1;
}


// int set_part5propdisp(int vals[]) {
//   char *token;

//   for(i = 0; i<npart5prop; i++){
//     partpropdata *propi;
//     int j;

//     propi = part5propinfo + i;
//     fgets(buffer, 255, stream);

//     trim_back(buffer);
//     token = strtok(buffer, " ");
//     j = 0;
//     while(token != NULL&&j<npartclassinfo){
//       int visval;

//       sscanf(token, "%i", &visval);
//       propi->class_vis[j] = visval;
//       token = strtok(NULL, " ");
//       j++;
//     }
//   }
//   CheckMemory;
//   continue;
// } // PART5PROPDISP

// int set_part5color(int n, int vals[]) {
//   int i;
//   for(i = 0; i<npart5prop; i++){
//     partpropdata *propi;

//     propi = part5propinfo + i;
//     propi->display = 0;
//   }
//   part5colorindex = 0;
//   i = n;
//   if(i >= 0 && i<npart5prop){
//     partpropdata *propi;

//     part5colorindex = i;
//     propi = part5propinfo + i;
//     propi->display = 1;
//   }
//   continue;
//   return 0;
// } // PART5COLOR

int lua_set_propindex(lua_State *L) {
  lua_pushnil(L);
  int n = 0;
  while(lua_next(L, -2) != 0) {
    lua_pop(L,1);
    n++;
  }
  int i = 0;
  int vals[n][2];
  lua_pushnil(L);
  while(lua_next(L, -2) != 0) {
    lua_pushnumber(L, 1);
    lua_gettable(L, -2);
    vals[i][0] = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_pushnumber(L, 1);
    lua_gettable(L, -2);
    vals[i][1] = lua_tonumber(L, -1);
    lua_pop(L, 1);

    lua_pop(L, 1);
    i++;
  }
  int return_code = set_propindex(n, vals);
  lua_pushnumber(L, return_code);
  return 1;
}



// int set_shooter(float xyz[], float dxyz[], float uvw[],
//                 float velmag, float veldir, float pointsize,
//                 int fps, int vel_type, int nparts, int vis, int cont_update,
//                 float duration, float v_inf) {
//   shooter_xyz[0] = xyz[0];
//   shooter_xyz[1] = xyz[1];
//   shooter_xyz[2] = xyz[2];

//   shooter_dxyz[0] = dxyz[0];
//   shooter_dxyz[1] = dxyz[1];
//   shooter_dxyz[2] = dxyz[2];

//   shooter_uvw[0] = uvw[0];
//   shooter_uvw[1] = uvw[1];
//   shooter_uvw[2] = uvw[2];

//   shooter_velmag = velmag;
//   shooter_veldir = veldir;
//   shooterpointsize = pointsize;

//   shooter_fps = fps;
//   shooter_vel_type = vel_type;
//   shooter_nparts = nparts;
//   visShooter = vis;
//   shooter_cont_update = cont_update;

//   shooter_duration = duration;
//   shooter_v_inf = v_inf;

//   return 0;
// } // SHOOTER

int lua_set_showdevices(lua_State *L) {
  lua_pushnil(L);
  int n = 0;
  while(lua_next(L, -2) != 0) {
    lua_pop(L,1);
    n++;
  }
  int i = 0;
  const char *names[n];
  lua_pushnil(L);
  while(lua_next(L, -2) != 0) {
    names[i] = lua_tostring(L, -1);
    lua_pop(L, 1);
    i++;
  }
  int return_code = set_showdevices(n, names);
  lua_pushnumber(L, return_code);
  return 1;
} // SHOWDEVICES

int lua_set_showdevicevals(lua_State *L) {
  int a = lua_tonumber(L, 1);
  int b = lua_tonumber(L, 2);
  int c = lua_tonumber(L, 3);
  int d = lua_tonumber(L, 4);
  int e = lua_tonumber(L, 5);
  int f = lua_tonumber(L, 6);
  int g = lua_tonumber(L, 7);
  int h = lua_tonumber(L, 8);
  int return_code = set_showdevicevals(a,b,c,d,e,f,g,h);
  lua_pushnumber(L, return_code);
  return 1;
} // SHOWDEVICEVALS

int lua_set_showmissingobjects(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_showmissingobjects(v);
  lua_pushnumber(L, return_code);
  return 1;
} // SHOWMISSINGOBJECTS

int lua_set_tourindex(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_tourindex(v);
  lua_pushnumber(L, return_code);
  return 1;
} // TOURINDEX

// int set_userticks(int vis, int auto_place, int sub, float origin[],
//                   float min[], float max[], float step[],
//                   int show_x, int show_y, int show_z) {
//   visUSERticks = vis;
//   auto_user_tick_placement = auto_place;
//   user_tick_sub = sub;

//   user_tick_origin[0] = origin[0];
//   user_tick_origin[1] = origin[1];
//   user_tick_origin[2] = origin[2];

//   user_tick_min[0] = min[0];
//   user_tick_min[1] = min[1];
//   user_tick_min[2] = min[2];

//   user_tick_max[0] = max[0];
//   user_tick_max[1] = max[1];
//   user_tick_max[2] = max[2];

//   user_tick_step[0] = step[0];
//   user_tick_step[1] = step[1];
//   user_tick_step[2] = step[2];

//   user_tick_show_x = show_x;
//   user_tick_show_y = show_y;
//   user_tick_show_z = show_z;

//   return 0;
// } // USERTICKS

int lua_set_c_particles(lua_State *L) {
  int minFlag = lua_tonumber(L, 1);
  float minValue = lua_tonumber(L, 2);
  int maxFlag = lua_tonumber(L, 3);
  float maxValue = lua_tonumber(L, 4);
  const char *label = NULL;
  if (lua_gettop(L) == 5) {
    label = lua_tostring(L, 5);
  }
  int return_code = set_c_particles(minFlag, minValue, maxFlag, maxValue,
                                    label);
  lua_pushnumber(L, return_code);
  return 1;
} // C_PARTICLES

int lua_set_c_slice(lua_State *L) {
  int minFlag = lua_tonumber(L, 1);
  float minValue = lua_tonumber(L, 2);
  int maxFlag = lua_tonumber(L, 3);
  float maxValue = lua_tonumber(L, 4);
  const char *label = NULL;
  if (lua_gettop(L) == 5) {
    label = lua_tostring(L, 5);
  }
  int return_code = set_c_slice(minFlag, minValue, maxFlag, maxValue,
                                    label);
  lua_pushnumber(L, return_code);
  return 1;
} // C_SLICE

int lua_set_cache_boundarydata(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_cache_boundarydata(v);
  lua_pushnumber(L, return_code);
  return 1;
} // CACHE_BOUNDARYDATA

int lua_set_cache_qdata(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_cache_qdata(v);
  lua_pushnumber(L, return_code);
  return 1;
} // CACHE_QDATA

int lua_set_percentilelevel(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_percentilelevel(v);
  lua_pushnumber(L, return_code);
  return 1;
} // PERCENTILELEVEL

int lua_set_timeoffset(lua_State *L) {
  int v = lua_tonumber(L, 1);
  int return_code = set_timeoffset(v);
  lua_pushnumber(L, return_code);
  return 1;
} // TIMEOFFSET

int lua_set_tload(lua_State *L) {
  int beginFlag = lua_tonumber(L, 1);
  float beginVal = lua_tonumber(L, 2);
  int endFlag = lua_tonumber(L, 3);
  float endVal = lua_tonumber(L, 4);
  int skipFlag = lua_tonumber(L, 5);
  float skipVal = lua_tonumber(L, 6);
  int return_code = set_tload(beginFlag, beginVal, endFlag, endVal,
                              skipFlag, skipVal);
  lua_pushnumber(L, return_code);
  return 1;
} // TLOAD

int lua_set_v_slice(lua_State *L) {
  int minFlag = lua_tonumber(L, 1);
  float minValue = lua_tonumber(L, 2);
  int maxFlag = lua_tonumber(L, 3);
  float maxValue = lua_tonumber(L, 4);
  const char *label = lua_tostring(L, 5);
  float lineMin = lua_tonumber(L, 6);
  float lineMax = lua_tonumber(L, 7);
  int lineNum = lua_tonumber(L, 8);
  int return_code = set_v_slice(minFlag, minValue, maxFlag, maxValue,
                                label, lineMin, lineMax, lineNum);
  lua_pushnumber(L, return_code);
  return 1;
}


int lua_set_patchdataout(lua_State *L) {
  int outputFlag = lua_tonumber(L, 1);
  int tmin = lua_tonumber(L, 1);
  int tmax = lua_tonumber(L, 2);
  int xmin = lua_tonumber(L, 3);
  int xmax = lua_tonumber(L, 4);
  int ymin = lua_tonumber(L, 5);
  int ymax = lua_tonumber(L, 6);
  int zmin = lua_tonumber(L, 7);
  int zmax = lua_tonumber(L, 8);
  int return_code = set_patchdataout(outputFlag, tmin, tmax, xmin, xmax, ymin,
                                     ymax, zmin, zmax);
  lua_pushnumber(L, return_code);
  return 1;
} // PATCHDATAOUT

int lua_show_smoke3d_showall(lua_State *L) {
  int return_code = show_smoke3d_showall();
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_show_smoke3d_hideall(lua_State *L) {
  int return_code = show_smoke3d_hideall();
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_show_slices_showall(lua_State *L) {
  int return_code = show_slices_showall();
  lua_pushnumber(L, return_code);
  return 1;
}

int lua_show_slices_hideall(lua_State *L) {
  int return_code = show_slices_hideall();
  lua_pushnumber(L, return_code);
  return 1;
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
  lua_register(L, "render_var", lua_render_var);
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

  // colorbar
  lua_register(L, "set_colorbar_visibility", lua_set_colorbar_visibility);
  lua_register(L, "get_colorbar_visibility", lua_get_colorbar_visibility);
  lua_register(L, "toggle_colorbar_visibility", lua_toggle_colorbar_visibility);

  // timebar
  lua_register(L, "set_timebar_visibility", lua_set_timebar_visibility);
  lua_register(L, "get_timebar_visibility", lua_get_timebar_visibility);
  lua_register(L, "toggle_timebar_visibility", lua_toggle_timebar_visibility);

  // title
  lua_register(L, "set_title_visibility", lua_set_title_visibility);
  lua_register(L, "get_title_visibility", lua_get_title_visibility);
  lua_register(L, "toggle_title_visibility", lua_toggle_title_visibility);

  // axis
  lua_register(L, "set_axis_visibility", lua_set_axis_visibility);
  lua_register(L, "get_axis_visibility", lua_get_axis_visibility);
  lua_register(L, "toggle_axis_visibility", lua_toggle_axis_visibility);

  // frame label
  lua_register(L, "set_framelabel_visibility", lua_set_framelabel_visibility);
  lua_register(L, "get_framelabel_visibility", lua_get_framelabel_visibility);
  lua_register(L, "toggle_framelabel_visibility", lua_toggle_framelabel_visibility);

  // framerate
  lua_register(L, "set_framerate_visibility", lua_set_framerate_visibility);
  lua_register(L, "get_framerate_visibility", lua_get_framerate_visibility);
  lua_register(L, "toggle_framerate_visibility", lua_toggle_framerate_visibility);

  // grid locations
  lua_register(L, "set_gridloc_visibility", lua_set_gridloc_visibility);
  lua_register(L, "get_gridloc_visibility", lua_get_gridloc_visibility);
  lua_register(L, "toggle_gridloc_visibility", lua_toggle_gridloc_visibility);

  // hrrpuv cutoff
  lua_register(L, "set_hrrcutoff_visibility", lua_set_hrrcutoff_visibility);
  lua_register(L, "get_hrrcutoff_visibility", lua_get_hrrcutoff_visibility);
  lua_register(L, "toggle_hrrcutoff_visibility", lua_toggle_hrrcutoff_visibility);

  // hrr label
  lua_register(L, "set_hrrlabel_visibility", lua_set_hrrlabel_visibility);
  lua_register(L, "get_hrrlabel_visibility", lua_get_hrrlabel_visibility);
  lua_register(L, "toggle_hrrlabel_visibility", lua_toggle_hrrlabel_visibility);

  // memory load
#ifdef pp_memstatus
  lua_register(L, "set_memload_visibility", lua_set_memload_visibility);
  lua_register(L, "get_memload_visibility", lua_get_memload_visibility);
  lua_register(L, "toggle_memload_visibility", lua_toggle_memload_visibility);
#endif

  // mesh label
  lua_register(L, "set_meshlabel_visibility", lua_set_meshlabel_visibility);
  lua_register(L, "get_meshlabel_visibility", lua_get_meshlabel_visibility);
  lua_register(L, "toggle_meshlabel_visibility", lua_toggle_meshlabel_visibility);

  // slice average
  lua_register(L, "set_slice_average_visibility", lua_set_slice_average_visibility);
  lua_register(L, "get_slice_average_visibility", lua_get_slice_average_visibility);
  lua_register(L, "toggle_slice_average_visibility", lua_toggle_slice_average_visibility);

  // time
  lua_register(L, "set_time_visibility", lua_set_time_visibility);
  lua_register(L, "get_time_visibility", lua_get_time_visibility);
  lua_register(L, "toggle_time_visibility", lua_toggle_time_visibility);

  // user settable ticks
  lua_register(L, "set_user_ticks_visibility", lua_set_user_ticks_visibility);
  lua_register(L, "get_user_ticks_visibility", lua_get_user_ticks_visibility);
  lua_register(L, "toggle_user_ticks_visibility", lua_toggle_user_ticks_visibility);

  // version info
  lua_register(L, "set_version_info_visibility", lua_set_version_info_visibility);
  lua_register(L, "get_version_info_visibility", lua_get_version_info_visibility);
  lua_register(L, "toggle_version_info_visibility", lua_toggle_version_info_visibility);

  // set all
  lua_register(L, "set_all_label_visibility", lua_set_all_label_visibility);

  // set the blockage view method
  lua_register(L, "blockage_view_method", lua_blockage_view_method);
  lua_register(L, "blockage_outline_color", lua_blockage_outline_color);
  lua_register(L, "blockage_locations", lua_blockage_locations);

  lua_register(L, "set_colorbar_colors", lua_set_colorbar_colors);
  lua_register(L, "get_colorbar_colors", lua_get_colorbar_colors);
  lua_register(L, "set_color2bar_colors", lua_set_color2bar_colors);
  lua_register(L, "get_color2bar_colors", lua_get_color2bar_colors);

  // Camera API
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

  lua_register(L, "camera_get_rotation_type" , lua_camera_get_rotation_type);
  lua_register(L, "camera_get_rotation_index", lua_camera_get_rotation_index);
  lua_register(L, "camera_set_rotation_type", lua_camera_set_rotation_type);
  lua_register(L, "camera_get_projection_type", lua_camera_get_projection_type);
  lua_register(L, "camera_set_projection_type", lua_camera_set_projection_type);

  lua_register(L, "get_clipping_mode", lua_get_clipping_mode);
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

  lua_register(L, "set_ambientlight", lua_set_ambientlight);
  lua_register(L, "set_backgroundcolor", lua_set_backgroundcolor);
  lua_register(L, "set_blockcolor", lua_set_blockcolor);
  lua_register(L, "set_blockshininess", lua_set_blockshininess);
  lua_register(L, "set_blockspecular", lua_set_blockspecular);
  lua_register(L, "set_boundcolor", lua_set_boundcolor);
  lua_register(L, "set_colorbar_textureflag", lua_set_colorbar_textureflag);
  lua_register(L, "set_diffuselight", lua_set_diffuselight);
  lua_register(L, "set_directioncolor", lua_set_directioncolor);
  lua_register(L, "set_flip", lua_set_flip);
  lua_register(L, "set_foregroundcolor", lua_set_foregroundcolor);
  lua_register(L, "set_heatoffcolor", lua_set_heatoffcolor);
  lua_register(L, "set_heatoncolor", lua_set_heatoncolor);
  lua_register(L, "set_isocolors", lua_set_isocolors);
  lua_register(L, "set_colortable", lua_set_colortable);
  lua_register(L, "set_light0", lua_set_light0);
  lua_register(L, "set_light1", lua_set_light1);
  lua_register(L, "set_lightpos0", lua_set_lightpos0);
  lua_register(L, "set_lightpos1", lua_set_lightpos1);
  lua_register(L, "set_lightmodellocalviewer", lua_set_lightmodellocalviewer);
  lua_register(L, "set_lightmodelseparatespecularcolor",
               lua_set_lightmodelseparatespecularcolor);
  lua_register(L, "set_sensorcolor", lua_set_sensorcolor);
  lua_register(L, "set_sensornormcolor", lua_set_sensornormcolor);
  lua_register(L, "set_bw", lua_set_bw);
  lua_register(L, "set_sprinkleroffcolor", lua_set_sprinkleroffcolor);
  lua_register(L, "set_sprinkleroncolor", lua_set_sprinkleroncolor);
  lua_register(L, "set_staticpartcolor", lua_set_staticpartcolor);
  lua_register(L, "set_timebarcolor", lua_set_timebarcolor);
  lua_register(L, "set_ventcolor", lua_set_ventcolor);
  lua_register(L, "set_gridlinewidth", lua_set_gridlinewidth);
  lua_register(L, "set_isolinewidth", lua_set_isolinewidth);
  lua_register(L, "set_isopointsize", lua_set_isopointsize);
  lua_register(L, "set_linewidth", lua_set_linewidth);
  lua_register(L, "set_partpointsize", lua_set_partpointsize);
  lua_register(L, "set_plot3dlinewidth", lua_set_plot3dlinewidth);
  lua_register(L, "set_plot3dpointsize", lua_set_plot3dpointsize);
  lua_register(L, "set_sensorabssize", lua_set_sensorabssize);
  lua_register(L, "set_sensorrelsize", lua_set_sensorrelsize);
  lua_register(L, "set_sliceoffset", lua_set_sliceoffset);
  lua_register(L, "set_smoothlines", lua_set_smoothlines);
  lua_register(L, "set_spheresegs", lua_set_spheresegs);
  lua_register(L, "set_sprinklerabssize", lua_set_sprinklerabssize);
  lua_register(L, "set_streaklinewidth", lua_set_streaklinewidth);
  lua_register(L, "set_ticklinewidth", lua_set_ticklinewidth);
  lua_register(L, "set_usenewdrawface", lua_set_usenewdrawface);
  lua_register(L, "set_veccontours", lua_set_veccontours);
  lua_register(L, "set_veclength", lua_set_veclength);
  lua_register(L, "set_vectorlinewidth", lua_set_vectorlinewidth);
  lua_register(L, "set_vectorpointsize", lua_set_vectorpointsize);
  lua_register(L, "set_ventlinewidth", lua_set_ventlinewidth);
  lua_register(L, "set_ventoffset", lua_set_ventoffset);
  lua_register(L, "set_windowoffset", lua_set_windowoffset);
  lua_register(L, "set_windowwidth", lua_set_windowwidth);
  lua_register(L, "set_windowheight", lua_set_windowheight);

  lua_register(L, "set_boundzipstep", lua_set_boundzipstep);
  lua_register(L, "set_fed", lua_set_fed);
  lua_register(L, "set_fedcolorbar", lua_set_fedcolorbar);
  lua_register(L, "set_isozipstep", lua_set_isozipstep);
  lua_register(L, "set_nopart", lua_set_nopart);
  lua_register(L, "set_showfedarea", lua_set_showfedarea);
  lua_register(L, "set_sliceaverage", lua_set_sliceaverage);
  lua_register(L, "set_slicedataout", lua_set_slicedataout);
  lua_register(L, "set_slicezipstep", lua_set_slicezipstep);
  lua_register(L, "set_smoke3dzipstep", lua_set_smoke3dzipstep);
  lua_register(L, "set_userrotate", lua_set_userrotate);

  lua_register(L, "set_aperture", lua_set_aperture);
  lua_register(L, "set_axissmooth", lua_set_axissmooth);
  lua_register(L, "set_blocklocation", lua_set_blocklocation);
  lua_register(L, "set_boundarytwoside", lua_set_boundarytwoside);
  lua_register(L, "set_clip", lua_set_clip);
  lua_register(L, "set_contourtype", lua_set_contourtype);
  lua_register(L, "set_cullfaces", lua_set_cullfaces);
  lua_register(L, "set_texturelighting", lua_set_texturelighting);
  lua_register(L, "set_eyeview", lua_set_eyeview);
  lua_register(L, "set_eyex", lua_set_eyex);
  lua_register(L, "set_eyey", lua_set_eyey);
  lua_register(L, "set_eyez", lua_set_eyez);
  lua_register(L, "set_fontsize", lua_set_fontsize);
  lua_register(L, "set_frameratevalue", lua_set_frameratevalue);
  lua_register(L, "set_geomdiags", lua_set_geomdiags);
  lua_register(L, "set_showfaces_interior", lua_set_showfaces_interior);
  lua_register(L, "set_showfaces_exterior", lua_set_showfaces_exterior);
  lua_register(L, "set_showfaces_solid", lua_set_showfaces_solid);
  lua_register(L, "set_showfaces_outline", lua_set_showfaces_outline);
  lua_register(L, "set_smoothgeomnormal", lua_set_smoothgeomnormal);
  lua_register(L, "set_showvolumes_interior", lua_set_showvolumes_interior);
  lua_register(L, "set_showvolumes_exterior", lua_set_showvolumes_exterior);
  lua_register(L, "set_showvolumes_solid", lua_set_showvolumes_solid);
  lua_register(L, "set_showvolumes_outline", lua_set_showvolumes_outline);
  lua_register(L, "set_geomvertexag", lua_set_geomvertexag);
  lua_register(L, "set_geommaxangle", lua_set_geommaxangle);
  lua_register(L, "set_gversion", lua_set_gversion);
  lua_register(L, "set_isotran2", lua_set_isotran2);
  lua_register(L, "set_meshvis", lua_set_meshvis);
  lua_register(L, "set_meshoffset", lua_set_meshoffset);

  lua_register(L, "set_northangle", lua_set_northangle);
  lua_register(L, "set_offsetslice", lua_set_offsetslice);
  lua_register(L, "set_outlinemode", lua_set_outlinemode);
  lua_register(L, "set_p3dsurfacetype", lua_set_p3dsurfacetype);
  lua_register(L, "set_p3dsurfacesmooth", lua_set_p3dsurfacesmooth);
  lua_register(L, "set_projection", lua_set_projection);
  lua_register(L, "set_scaledfont", lua_set_scaledfont);
  lua_register(L, "set_showalltextures", lua_set_showalltextures);
  lua_register(L, "set_showaxislabels", lua_set_showaxislabels);
  lua_register(L, "set_showblocklabel", lua_set_showblocklabel);
  lua_register(L, "set_showblocks", lua_set_showblocks);
  lua_register(L, "set_showcadandgrid", lua_set_showcadandgrid);
  lua_register(L, "set_showcadopaque", lua_set_showcadopaque);
  lua_register(L, "set_showceiling", lua_set_showceiling);
  lua_register(L, "set_showcolorbars", lua_set_showcolorbars);
  lua_register(L, "set_showcvents", lua_set_showcvents);
  lua_register(L, "set_showdummyvents", lua_set_showdummyvents);
  lua_register(L, "set_showevacslices", lua_set_showevacslices);
  lua_register(L, "set_showfloor", lua_set_showfloor);
  lua_register(L, "set_showframe", lua_set_showframe);
  lua_register(L, "set_showframelabel", lua_set_showframelabel);
  lua_register(L, "set_showframerate", lua_set_showframerate);
  lua_register(L, "set_showgrid", lua_set_showgrid);
  lua_register(L, "set_showgridloc", lua_set_showgridloc);
  lua_register(L, "set_showhmstimelabel", lua_set_showhmstimelabel);
  lua_register(L, "set_showhrrcutoff", lua_set_showhrrcutoff);
  lua_register(L, "set_showiso", lua_set_showiso);
  lua_register(L, "set_showisonormals", lua_set_showisonormals);
  lua_register(L, "set_showlabels", lua_set_showlabels);
#ifdef pp_memstatus
  lua_register(L, "set_showmemload", lua_set_showmemload);
#endif
  lua_register(L, "set_showopenvents", lua_set_showopenvents);
  lua_register(L, "set_showothervents", lua_set_showothervents);
  lua_register(L, "set_showsensors", lua_set_showsensors);
  lua_register(L, "set_showsliceinobst", lua_set_showsliceinobst);
  lua_register(L, "set_showsmokepart", lua_set_showsmokepart);
  lua_register(L, "set_showsprinkpart", lua_set_showsprinkpart);
  lua_register(L, "set_showstreak", lua_set_showstreak);
  lua_register(L, "set_showterrain", lua_set_showterrain);
  lua_register(L, "set_showtetras", lua_set_showterrain);
  lua_register(L, "set_showthreshold", lua_set_showthreshold);
  lua_register(L, "set_showticks", lua_set_showticks);
  lua_register(L, "set_showtimebar", lua_set_showtimebar);
  lua_register(L, "set_showtimelabel", lua_set_showtimelabel);
  lua_register(L, "set_showtitle", lua_set_showtitle);
  lua_register(L, "set_showtracersalways", lua_set_showtracersalways);
  lua_register(L, "set_showtriangles", lua_set_showtriangles);
  lua_register(L, "set_showtransparent", lua_set_showtransparent);
  lua_register(L, "set_showtransparentvents", lua_set_showtranparentvents);
  lua_register(L, "set_showtrianglecount", lua_set_showtrianglecount);
  lua_register(L, "set_showventflow", lua_set_showventflow);
  lua_register(L, "set_showvents", lua_set_showvents);
  lua_register(L, "set_showwalls", lua_set_showwalls);
  lua_register(L, "set_skipembedslice", lua_set_skipembedslice);
#ifdef pp_SLICEUP
  lua_register(L, "set_slicedup", lua_set_slicedup);
#endif
  lua_register(L, "set_smokesensors", lua_set_smokesensors);
#ifdef pp_LANG
  lua_register(L, "set_startuplang", lua_set_startuplang);
#endif
  lua_register(L, "set_stereo", lua_set_stereo);
  lua_register(L, "set_surfinc", lua_set_surfinc);
  lua_register(L, "set_terrainparams", lua_set_terrainparams);
  lua_register(L, "set_titlesafe", lua_set_titlesafe);
  lua_register(L, "set_trainermode", lua_set_trainermode);
  lua_register(L, "set_trainerview", lua_set_trainerview);
  lua_register(L, "set_transparent", lua_set_transparent);
  lua_register(L, "set_treeparms", lua_set_treeparms);
  lua_register(L, "set_twosidedvents", lua_set_twosidedvents);
  lua_register(L, "set_vectorskip", lua_set_vectorskip);
  lua_register(L, "set_volsmoke", lua_set_volsmoke);
  lua_register(L, "set_zoom", lua_set_zoom);
  lua_register(L, "set_cellcentertext", lua_set_cellcentertext);
  lua_register(L, "set_inputfile", lua_set_inputfile);
  lua_register(L, "set_labelstartupview", lua_set_labelstartupview);
  lua_register(L, "set_pixelskip", lua_set_pixelskip);
  lua_register(L, "set_renderclip", lua_set_renderclip);
  lua_register(L, "set_renderfilelabel", lua_set_renderfilelabel);
  lua_register(L, "set_renderfiletype", lua_set_renderfiletype);

  // lua_register(L, "set_skybox", lua_set_skybox);
  lua_register(L, "set_renderoption", lua_set_renderoption);
  lua_register(L, "set_unitclasses", lua_set_unitclasses);
  lua_register(L, "set_zaxisangles", lua_set_zaxisangles);
  lua_register(L, "set_adjustalpha", lua_set_adjustalpha);
  lua_register(L, "set_colorbartype", lua_set_colorbartype);
  lua_register(L, "set_extremecolors", lua_set_extremecolors);
  lua_register(L, "set_firecolor", lua_set_firecolor);
  lua_register(L, "set_firecolormap", lua_set_firecolormap);
  lua_register(L, "set_firedepth", lua_set_firedepth);
  // lua_register(L, "set_golorbar", lua_set_gcolorbar);
  lua_register(L, "set_showextremedata", lua_set_showextremedata);
  lua_register(L, "set_smokecolor", lua_set_smokecolor);
  lua_register(L, "set_smokecull", lua_set_smokecull);
  lua_register(L, "set_smokeskip", lua_set_smokeskip);
  lua_register(L, "set_smokealbedo", lua_set_smokealbedo);
#ifdef pp_GPU // TODO: register anyway, but tell user it is not available
  lua_register(L, "set_smokerthick", lua_set_smokerthick);
#endif
  lua_register(L, "set_smokethick", lua_set_smokethick);
#ifdef pp_GPU
  lua_register(L, "set_usegpu", lua_set_usegpu);
#endif
  lua_register(L, "set_showhazardcolors", lua_set_showhazardcolors);
  lua_register(L, "set_showhzone", lua_set_showhzone);
  lua_register(L, "set_showszone", lua_set_showszone);
  lua_register(L, "set_showvzone", lua_set_showvzone);
  lua_register(L, "set_showzonefire", lua_set_showzonefire);
  lua_register(L, "set_showpathnodes", lua_set_showpathnodes);
  lua_register(L, "set_showtourroute", lua_set_showtourroute);
  lua_register(L, "set_tourcolors_selectedpathline",
               lua_set_tourcolors_selectedpathline);
  lua_register(L, "set_tourcolors_selectedpathlineknots",
               lua_set_tourcolors_selectedpathlineknots);
  lua_register(L, "set_tourcolors_selectedknot",
               lua_set_tourcolors_selectedknot);
  lua_register(L, "set_tourcolors_pathline",
               lua_set_tourcolors_pathline);
  lua_register(L, "set_tourcolors_pathknots",
               lua_set_tourcolors_pathknots);
  lua_register(L, "set_tourcolors_text",
               lua_set_tourcolors_text);
  lua_register(L, "set_tourcolors_avatar",
               lua_set_tourcolors_avatar);
  lua_register(L, "set_tourconstantvel", lua_set_tourconstantvel);
  lua_register(L, "set_viewalltours", lua_set_viewalltours);
  lua_register(L, "set_viewtimes", lua_set_viewtimes);
  lua_register(L, "set_viewtourfrompath", lua_set_viewtourfrompath);
  lua_register(L, "set_avatarevac", lua_set_avatarevac);
  lua_register(L, "set_geometrytest", lua_set_geometrytest);
  lua_register(L, "set_devicevectordimensions",
               lua_set_devicevectordimensions);
  lua_register(L, "set_devicebounds", lua_set_devicebounds);
  lua_register(L, "set_deviceorientation", lua_set_deviceorientation);
  lua_register(L, "set_gridparms", lua_set_gridparms);
  lua_register(L, "set_gsliceparms", lua_set_gsliceparms);
  lua_register(L, "set_loadfilesatstartup", lua_set_loadfilesatstartup);
  lua_register(L, "set_mscale", lua_set_mscale);
  lua_register(L, "set_sliceauto", lua_set_sliceauto);
  lua_register(L, "set_msliceauto", lua_set_msliceauto);
  lua_register(L, "set_compressauto", lua_set_compressauto);
  // lua_register(L, "set_part5propdisp", lua_set_part5propdisp);
  // lua_register(L, "set_part5color", lua_set_part5color);
  lua_register(L, "set_propindex", lua_set_propindex);
  // lua_register(L, "set_shooter", lua_set_shooter);
  lua_register(L, "set_showdevices", lua_set_showdevices);
  lua_register(L, "set_showdevicevals", lua_set_showdevicevals);
  lua_register(L, "set_showmissingobjects", lua_set_showmissingobjects);
  lua_register(L, "set_tourindex", lua_set_tourindex);
  // lua_register(L, "set_userticks", lua_set_userticks);
  lua_register(L, "set_c_particles", lua_set_c_particles);
  lua_register(L, "set_c_slice", lua_set_c_slice);
  lua_register(L, "set_cache_boundarydata", lua_set_cache_boundarydata);
  lua_register(L, "set_cache_qdata", lua_set_cache_qdata);
  lua_register(L, "set_percentilelevel", lua_set_percentilelevel);
  lua_register(L, "set_timeoffset", lua_set_timeoffset);
  lua_register(L, "set_tload", lua_set_tload);
  lua_register(L, "set_v_slice", lua_set_v_slice);
  lua_register(L, "set_patchdataout", lua_set_patchdataout);

  lua_register(L, "show_smoke3d_showall", lua_show_smoke3d_showall);
  lua_register(L, "show_smoke3d_hideall", lua_show_smoke3d_hideall);
  lua_register(L, "show_slices_showall", lua_show_slices_showall);
  lua_register(L, "show_slices_hideall", lua_show_slices_hideall);

  lua_register(L, "get_nglobal_times", lua_get_nglobal_times);

  //add fdsprefix (the path plus  CHID) as a variable in the lua environment
  lua_pushstring(L, fdsprefix);
  lua_setglobal(L, "fdsprefix");

  lua_pushstring(L, chidfilebase);
  lua_setglobal(L, "chid");

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
  char cwd[1000];
  getcwd(cwd,1000);
  printf("cwd: %s\n", cwd);
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
