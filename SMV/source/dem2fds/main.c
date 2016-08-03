#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "string_util.h"
#include "file_util.h"
#include "datadefs.h"
#include "MALLOC.h"


#define FDS_OBST 0
#define FDS_GEOM 1
#define LEN_BUFFER 1024
#define EARTH_RADIUS 6371000.0
#define INTERP1D(f,v1,v2) ((1.0-(f))*(v1)+(f)*(v2))
#define DONT_INTERPOLATE 0
#define INTERPOLATE 1

/* --------------------------  elevdata ------------------------------------ */

typedef struct {
  int ncols, nrows, use_it;
  float xllcorner, yllcorner, cellsize;
  char *fileheader, *filedata;
  float lat_min, lat_max;
  float long_min, long_max;
  float val_min, val_max;
  float dlong, dlat;
  float *valbuffer;
  float xmax, ymax, zmin, zmax;
  float xref, yref;
  float longref, latref;
} elevdata;

/* ------------------ ShowExample ------------------------ */

void ShowExample(void){
  fprintf(stderr, " Example input file for generating an FDS input file from elevation data\n\n");
  fprintf(stderr, " // minimum longitude, maximum longitude, number of longitudes\n");
  fprintf(stderr, " LONGMINMAX\n");
  fprintf(stderr, "  -77.25 -77.20 100\n\n");
  fprintf(stderr, " // minimum latitude, maximum latitude, number of latitudes\n");
  fprintf(stderr, " LATMINMAX\n");
  fprintf(stderr, "  39.12 39.15 100\n");
}

  /* ------------------ Usage ------------------------ */

void Usage(char *prog){
 char githash[LEN_BUFFER];
 char gitdate[LEN_BUFFER];

  getGitInfo(githash,gitdate);    // get githash

  fprintf(stderr, "\n%s (%s) %s\n", prog, githash, __DATE__);
  fprintf(stderr, "Create an FDS input file using elevation data\n");
  fprintf(stderr, "  obtained from http://viewer.nationalmap.gov \n\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  dem2fds [-g|-o][-h][-v] casename.in\n");
  fprintf(stderr, "where\n");
  fprintf(stderr, "  -e - show an example input file\n");
  fprintf(stderr, "  -g - create an FDS input file using &GEOM keywords\n");
  fprintf(stderr, "  -o - create an FDS input file using &OBST keywords (default)\n");
  fprintf(stderr, "  -h - display this message\n");
  fprintf(stderr, "  -v - show version information\n");
}

/* ------------------ GetLongLats ------------------------ */

void GetLongLats(
  float longref, float latref, float xref, float yref,
  float xmax, float ymax, int nxmax, int nymax,
  float *longlats
) {
  int j;
  float dx, dy;

  dy = ymax / (float)(nymax - 1);
  dx = xmax / (float)(nxmax - 1);

  for (j = nymax-1; j >=0; j--) {
    int i;
    float dlat, dyval, coslat;

    dyval = (float)j*dy - yref;
    dlat = dyval / EARTH_RADIUS;
    coslat = cos(latref + dlat);
    for (i = 0; i < nxmax; i++) {
      float dxval, top, dlong;

      dxval = (float)i*dx - xref;
      top = sin(dxval / (2.0*EARTH_RADIUS));
      dlong = 2.0*asin(top / coslat);
      *longlats++ = longref + RAD2DEG*dlong;
      *longlats++ = latref + RAD2DEG*dlat;
    }
  }
}

/* ------------------ SphereDistance ------------------------ */

float SphereDistance(float llong1, float llat1, float llong2, float llat2){
  // https://en.wikipedia.org/wiki/Great-circle_distance
  // a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlong/2)^2
  // c = 2*asin(sqrt(a))
  // d = R*c
  // R = RAD_EARTH

  float deg2rad;
  float a, c;
  float dlat, dlong;

  deg2rad = 4.0*atan(1.0)/180.0;
  llat1 *= deg2rad;
  llat2 *= deg2rad;
  llong1 *= deg2rad;
  llong2 *= deg2rad;
  dlat = llat2 - llat1;
  dlong = llong2 - llong1;
  a = pow(sin(dlat / 2.0), 2) + cos(llat1)*cos(llat2)*pow(sin(dlong / 2.0), 2);
  c = 2.0 * asin(sqrt(ABS(a)));
  return EARTH_RADIUS*c;
}

/* ------------------ GenerateFDS ------------------------ */

void GenerateFDS(char *casename, elevdata *fds_elevs, int option){
  char fdsfile[LEN_BUFFER],*ext;
  char basename[LEN_BUFFER];
  int nlong, nlat,nz;
  int i,j;
  float llat1, llat2, llong1, llong2;
  float xmax, ymax, zmin, zmax;
  float *xgrid, *ygrid;
  int count;
  int ibar, jbar, kbar;
  float *vals, *valsp1;
  FILE *streamout = NULL;

  strcpy(basename, casename);
  ext = strrchr(basename, '.');
  if (ext != NULL)ext[0] = 0;

  strcpy(fdsfile,basename);
  strcat(fdsfile,".fds");
  streamout = fopen(fdsfile,"w");
  if(streamout==NULL){
    fprintf(stderr,"***error: unable to open %s for output\n",fdsfile);
    return;
  }

  llong1 = fds_elevs->long_min;
  llong2 = fds_elevs->long_max;
  nlong = fds_elevs->ncols;

  llat1 = fds_elevs->lat_min;
  llat2 = fds_elevs->lat_max;
  nlat = fds_elevs->nrows;

  zmin = fds_elevs->val_min;
  zmax = fds_elevs->val_max;
  nz = 30;

  vals = fds_elevs->valbuffer;

  xmax = fds_elevs->xmax;   
  ymax = fds_elevs->ymax;

  ibar = nlong - 1;
  jbar = nlat - 1;
  kbar = nz;

  NewMemory((void **)&xgrid, sizeof(float)*(ibar+1));
  for(i=0;i<ibar+1;i++){
    xgrid[i] = xmax*(float)i/(float)ibar;
  }

  NewMemory((void **)&ygrid, sizeof(float)*(jbar+1));
  for(i=0;i<jbar+1;i++){
    ygrid[i] = ymax*(float)i/(float)jbar;
  }


  fprintf(streamout,"&HEAD CHID='%s', TITLE='terrain' /\n",basename);
  fprintf(streamout,"&MESH IJK = %i, %i, %i, XB = 0.0, %f, 0.0, %f, %f, %f /\n",ibar,jbar,kbar,xmax,ymax,zmin,zmax);
  if(option==FDS_OBST){
    fprintf(streamout,"&MISC TERRAIN_CASE = .TRUE., TERRAIN_IMAGE = '%s.png' /\n", basename);
  }
  fprintf(streamout,"&TIME T_END = 0. /\n");
  fprintf(streamout,"&VENT XB = 0.0, 0.0, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", ymax, zmin, zmax);
  fprintf(streamout,"&VENT XB =  %f,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, xmax, ymax, zmin, zmax);
  fprintf(streamout,"&VENT XB = 0.0,  %f, 0.0, 0.0, %f, %f, SURF_ID = 'OPEN' /\n", xmax, zmin, zmax);
  fprintf(streamout,"&VENT XB = 0.0,  %f,  %f,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, ymax, ymax, zmin, zmax);
  fprintf(streamout,"&VENT XB = 0.0,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, ymax, zmax, zmax);
  fprintf(streamout,"&MATL ID = 'matl1', DENSITY = 1000., CONDUCTIVITY = 1., SPECIFIC_HEAT = 1., RGB = 122,117,48 /\n");
  fprintf(streamout,"&SURF ID = 'surf1', RGB = 122,117,48 TEXTURE_MAP='%s.png' /\n",basename);


  if(option==FDS_GEOM){
    fprintf(streamout,"&GEOM ID='terrain', SURF_ID='surf1',MATL_ID='matl1',\nIJK=%i,%i,XB=%f,%f,%f,%f,\nZVALS=\n",
                      nlong,nlat,0.0,xmax,0.0,ymax);
    count = 1;
    for(j = 0; j < jbar + 1; j++){
      for(i = 0; i < ibar + 1; i++){
        fprintf(streamout," %f,", vals[count-1]);
        if(count % 10 == 0)fprintf(streamout,"\n");
        count++;
      }
    }
    fprintf(streamout,"/\n");
  }
  if(option==FDS_OBST){
    count = 0;
    valsp1 = vals + nlong;
    for(j = 0; j < jbar; j++){
      float ycen;

      ycen = (ygrid[j] + ygrid[j + 1]) / 2.0;
      for(i = 0; i < ibar; i++){
        float vavg, xcen;

        xcen = (xgrid[i] + xgrid[i + 1]) / 2.0;
        vavg = (vals[count]+vals[count+1]+valsp1[count]+valsp1[count+1])/4.0;
        fprintf(streamout,"&OBST XB=%f,%f,%f,%f,0.0,%f SURF_ID='surf1'/\n",
                   xgrid[i],xgrid[i+1],ygrid[j],ygrid[j+1],vavg);
        count++;
      }
      count++;
    }
  }
  fprintf(streamout,"&TAIL /\n");
  
  fprintf(stderr,"FDS input file: %s\n",fdsfile);
  fprintf(stderr,"xmax: %f\n",xmax);
  fprintf(stderr,"ymax: %f\n",ymax);
  fprintf(stderr,"xref: %f\n",fds_elevs->xref);
  fprintf(stderr,"yref: %f\n",fds_elevs->yref);
  fprintf(stderr, "longref: %f\n", fds_elevs->longref);
  fprintf(stderr, "latref: %f\n", fds_elevs->latref);
}

/* ------------------ GetElevFile ------------------------ */

elevdata *GetElevFile(elevdata *elevinfo, int nelevinfo, float longval, float latval){
  int i;

  for(i = 0; i < nelevinfo; i++){
    elevdata *elevi;

    elevi = elevinfo + i;
    if(longval<elevi->long_min || longval>elevi->long_max)continue;
    if(latval<elevi->lat_min || latval>elevi->lat_max)continue;
    return elevi;
  }
  return NULL;
}

/* ------------------ GetElevation ------------------------ */

float GetElevation(elevdata *elevinfo, int nelevinfo, float longval, float latval, int interp_option, int *have_val){
  elevdata *elevi;
  int ival, jval;
  int ival2, jval2;
  int index11, index12, index21, index22;
  float val11, val12, val21, val22;
  float val1, val2;
  float factor_x, factor_y;
  float return_val;
  float ivalx, jvaly;

  *have_val = 0;
  elevi = GetElevFile(elevinfo, nelevinfo, longval, latval);
  if(elevi == NULL)return 0.0;
  if(elevi->valbuffer == NULL){
    FILE *stream;
    float *data_buffer;

    stream = fopen(elevi->filedata, "rb");
    if (stream == NULL)return 0.0;
    NewMemory((void **)&data_buffer, elevi->ncols*elevi->nrows * sizeof(float));
    elevi->valbuffer = data_buffer;
    fread(data_buffer, sizeof(float), elevi->ncols*elevi->nrows, stream);
    fclose(stream);
  }

  ivalx = (longval - elevi->long_min)/elevi->cellsize;
   ival = CLAMP(ivalx,0,elevi->ncols-1);

  jvaly = (elevi->lat_max - latval)/elevi->cellsize;
   jval = CLAMP(jvaly,0,elevi->nrows-1);

  index11 =  jval*elevi->ncols + ival;
  val11 = elevi->valbuffer[index11];
  if (interp_option == 0) {
    *have_val = 1;
    return val11;
  }

  ival2 = CLAMP(ival + 1, 0, elevi->ncols - 1);
  jval2 = CLAMP(jval - 1, 0, elevi->nrows - 1);

  index12 =  jval*elevi->ncols + ival2;
  val12 = elevi->valbuffer[index12];

  index21 = jval2*elevi->ncols + ival;
  val21 = elevi->valbuffer[index21];

  index22 = jval2*elevi->ncols + ival2;
  val22 = elevi->valbuffer[index22];

  factor_x = CLAMP(ivalx - ival,0.0,1.0);
  factor_y = CLAMP(jvaly - jval,0.0,1.0);

  val1 = INTERP1D(factor_x,val11,val12);
  val2 = INTERP1D(factor_x,val21,val22);

  return_val = INTERP1D(factor_y,val2,val1);
  *have_val = 1;
  return return_val;
}

/* ------------------ GenerateElevs ------------------------ */

int GenerateElevs(char *elevfile, elevdata *fds_elevs){
  int nelevinfo,i,j;
  filelistdata *fileheaders;
  FILE *stream_in;
  elevdata *elevinfo;
  int ibar, jbar, kbar;
  int longlat_defined = 0;
  float xmin_exclude, ymin_exclude, xmax_exclude, ymax_exclude;
  float longmin, longmax, latmin, latmax;
  int nlongs=100, nlats=100;
  float dlat, dlong;
  int count, *have_vals, have_data=0;
  float valmin, valmax, *vals;
  char *ext;
  float longref=-1000.0, latref=-1000.0;
  float xref=0.0, yref=0.0;
  float xmax = -1000.0, ymax = -1000.0, zmin=-1000.0, zmax=-1000.0;
  float *longlats = NULL, *longlatsorig;

  nelevinfo = get_nfilelist(".", "*.hdr");
  if(nelevinfo == 0)return 0;

  get_filelist(".","*.hdr", nelevinfo, &fileheaders);
  NewMemory((void **)&elevinfo, nelevinfo*sizeof(elevdata));
  for(i = 0; i < nelevinfo; i++){
    filelistdata *filei;
    elevdata *elevi;
    char file[LEN_BUFFER], *filedatai;
    int lenfile;

    filei = fileheaders + i;
    elevi = elevinfo + i;
    strcpy(file, filei->file);
    ext = strrchr(file, '.');
    if(ext!=NULL)ext[0] = 0;
    strcat(file, ".flt");
    lenfile = strlen(file);
    NewMemory((void **)&filedatai, (lenfile + 1) * sizeof(char));
    strcpy(filedatai, file);
    elevi->fileheader = filei->file;
    elevi->filedata = filedatai;
  }
  for(i = 0; i < nelevinfo; i++){
    elevdata *elevi;
    char buffer[LEN_BUFFER];

    elevi = elevinfo + i;
    elevi->use_it = 0;

    stream_in = fopen(elevi->fileheader, "r");
    if(stream_in == NULL)continue;

    if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)continue;
    trim_back(buffer);
    sscanf(buffer+5," %i", &elevi->ncols);

    if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)continue;
    trim_back(buffer);
    sscanf(buffer+5, " %i", &elevi->nrows);

    if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)continue;
    trim_back(buffer);
    sscanf(buffer+9, " %f", &elevi->xllcorner);

    if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)continue;
    trim_back(buffer);
    sscanf(buffer+9, " %f", &elevi->yllcorner);

    if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)continue;
    trim_back(buffer);
    sscanf(buffer+8, " %f", &elevi->cellsize);

    elevi->long_min = elevi->xllcorner;
    elevi->long_max = elevi->long_min + (float)elevi->ncols*elevi->cellsize;

    elevi->lat_min = elevi->yllcorner;
    elevi->lat_max = elevi->lat_min + (float)elevi->nrows*elevi->cellsize;

    elevi->valbuffer = NULL;

    elevi->use_it = 1;

    fclose(stream_in);
  }

  stream_in = fopen(elevfile, "r");
  if (stream_in == NULL) {
    fprintf(stderr,"***error: unable to open file %s for input\n", elevfile);
    return 0;
  }
  while(!feof(stream_in)){
    char buffer[LEN_BUFFER], *buffer2;

    CheckMemory;

    if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
    buffer2 = strstr(buffer, "//");
    if(buffer2 != NULL)buffer2[0] = 0;
    buffer2 = trim_frontback(buffer);
    if(strlen(buffer2) == 0)continue;

    if(match(buffer, "GRID") == 1){
      ibar = 10;
      jbar = 10;
      kbar = 10;
      if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%i %i %i", &ibar, &jbar, &kbar);
      continue;
    }

    if(match(buffer, "LONGLATREF") == 1){
      have_data = 1;
      longref = 1000.0;
      latref = 1000.0;
      if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f", &longref, &latref);
      continue;
    }

    if (match(buffer, "XYREF") == 1) {
      if (fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f", &xref, &yref);
      continue;
    }

    if (match(buffer, "XYMAX") == 1) {
      xmax = 1000.0;
      ymax = 1000.0;
      if (fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f", &xmax, &ymax);
      continue;
    }

    if (match(buffer, "ZMINMAX") == 1) {
      if (fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f", &zmin, &zmax);
      continue;
    }

    if(match(buffer, "LONGLATMINMAX") == 1){
      have_data = 1;
      longmin = -1000.0;
      longmax = -1000.0;
      nlongs = 50;
      ibar = nlongs;
      latmin = -1000.0;
      latmax = -1000.0;
      nlats = 50;
      kbar = nlats;
      if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f %i %f %f %i", &longmin, &longmax, &nlongs, &latmin, &latmax, &nlats);
      longref = (longmin + longmax) / 2.0;
      latref = (latmin + latmax) / 2.0;
      xmax = SphereDistance(longmin, latref, longmax, latref);
      ymax = SphereDistance(longref, latmin, longref, latmax);
      xref = xmax / 2.0;
      yref = ymax / 2.0;

      continue;
    }

    if(match(buffer, "EXCLUDE") == 1){
      if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f %f %f", &xmin_exclude, &ymin_exclude, &xmax_exclude, &ymax_exclude);
      longlat_defined = 1;
      continue;
    }
  }
  fclose(stream_in);

  NewMemory((void **)&longlatsorig, 2 * nlongs*nlats * sizeof(float));
  longlats = longlatsorig;
  GetLongLats(longref, latref, xref, yref,
    xmax, ymax, nlongs, nlats, longlats);

  for (i = 0; i < nlongs*nlats; i++) {
    float llong, llat;

    llong = *longlats++;
    llat = *longlats++;
    if (GetElevFile(elevinfo, nelevinfo, llong, llat) == NULL) {
      fprintf(stderr, "***error: elevation data not available for \n");
      fprintf(stderr, "    longitude/latitude: (%f %f) \n",llong, llat);
      for (i = 0; i < nelevinfo; i++) {
        elevdata *elevi;

        elevi = elevinfo + i;
        fprintf(stderr, " header file: %s bounds: %f %f %f %f\n",
          elevi->fileheader, elevi->long_min, elevi->lat_min, elevi->long_max, elevi->lat_max);
      }
      FREEMEMORY(longlatsorig);
      return 0;
    }
  }

  dlat = (latmax - latmin) / (float)(nlats-1);
  dlong = (longmax - longmin) / (float)(nlongs-1);
  NewMemory((void **)&vals, nlongs*nlats*sizeof(float));
  NewMemory((void **)&have_vals, nlongs*nlats*sizeof(int));
  longlats = longlatsorig;
  fds_elevs->valbuffer = vals;
  fds_elevs->nrows = nlats;
  fds_elevs->ncols = nlongs;
  fds_elevs->long_min = longmin;
  fds_elevs->long_max = longmax;
  fds_elevs->lat_min = latmin;
  fds_elevs->lat_max = latmax;
  fds_elevs->xmax = xmax;
  fds_elevs->ymax = ymax;
  fds_elevs->xref = xref;
  fds_elevs->yref = yref;
  fds_elevs->longref = longref;
  fds_elevs->latref = latref;

  count = 0;
  for (j = 0; j < nlats; j++) {
    for(i = 0; i < nlongs; i++){
      float llong, llat, elevij;
      int have_val;

      llong = *longlats++;
      llat = *longlats++;
      elevij = GetElevation(elevinfo, nelevinfo, llong, llat, INTERPOLATE, &have_val);
      vals[count]=elevij;
      have_vals[count] = have_val;
      if(have_val == 1){
        if(count == 0){
          valmin = elevij;
          valmax = elevij;
        }
        else{
          valmin = MIN(valmin, elevij);
          valmax = MAX(valmax, elevij);
        }
      }
      count++;
    }
  }
  fds_elevs->val_min = valmin;
  fds_elevs->val_max = valmax;
  FREEMEMORY(have_vals);
  FREEMEMORY(longlatsorig);
  return 1;
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  int gen_fds = FDS_OBST;
  char *casename = NULL;
  char file_default[LEN_BUFFER];
  elevdata fds_elevs;

  if(argc == 1){
    Usage("dem2fds");
    return 0;
  }

  strcpy(file_default, "terrain");

  initMALLOC();
  set_stdout(stdout);
  for(i = 1; i<argc; i++){
    int lenarg;
    char *arg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'e':
        ShowExample();
        exit(1);
        break;
      case 'h':
        Usage("dem2fds");
        exit(1);
        break;
      case 'o':
        gen_fds = FDS_OBST;
        break;
      case 'g':
        gen_fds = FDS_GEOM;
        break;
      case 'v':
        PRINTversion("dem2fds");
        exit(1);
        break;
      default:
        Usage("dem2fds");
        exit(1);
        break;
      }
    }
    else{
      if(casename == NULL)casename = argv[i];
    }
  }
  if(casename == NULL)casename = file_default;
  if (GenerateElevs(casename,&fds_elevs) == 1) {
    GenerateFDS(casename, &fds_elevs, gen_fds);
  }
  return 0;
}
