#define FDS_OBST 0
#define FDS_GEOM 1
#define LEN_BUFFER 1024
#define EARTH_RADIUS 6371000.0
#define INTERP1D(f,v1,v2) ((1.0-(f))*(v1)+(f)*(v2))
#define DONT_INTERPOLATE 0
#define INTERPOLATE 1

/* --------------------------  elevdata ------------------------------------ */

typedef struct {
  int ncols, nrows, nz, use_it;
  float xllcorner, yllcorner, cellsize;
  char *headerfile, *datafile;
  char filelabel[13];
  float lat_min, lat_max, long_min, long_max, dlong, dlat;
  float lat_min_orig, lat_max_orig, long_min_orig, long_max_orig;
  float val_min, val_max;
  float xmax, ymax, zmin, zmax;
  float xref, yref, longref, latref;
  float *valbuffer;
  gdImagePtr image;
} elevdata;

EXTERNCPP void GenerateFDSInputFile(char *casename, elevdata *fds_elevs, int option);
EXTERNCPP int GetElevations(char *elevfile, elevdata *fds_elevs);

SVEXTERN char libdir[1024];
SVEXTERN int SVDECL(border_buffer,300), SVDECL(show_maps,0);

