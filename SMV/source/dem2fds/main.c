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

/* --------------------------  elevdata ------------------------------------ */

typedef struct {
  int ncols, nrows, use_it;
  float xllcorner, yllcorner, cellsize;
  char *fileheader, *filedata;
  float lat_min, lat_max;
  float long_min, long_max;
  float dlong, dlat;
  float *valbuffer;
} elevdata;

/* ------------------ example ------------------------ */

void show_example(void){
  fprintf(stderr, " Example input file for generating an FDS input file from elevation data\n\n");
  fprintf(stderr, " // minimum longitude, maximum longitude, number of longitudes\n");
  fprintf(stderr, " LONGMINMAX\n");
  fprintf(stderr, "  -77.25 -77.20 100\n\n");
  fprintf(stderr, " // minimum latitude, maximum latitude, number of latitudes\n");
  fprintf(stderr, " LATMINMAX\n");
  fprintf(stderr, "  39.12 39.15 100\n");
}

  /* ------------------ usage ------------------------ */

void usage(char *prog){
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

/* ------------------ dist ------------------------ */

float dist(float llong1, float llong2, float llat1, float llat2){
  // https://en.wikipedia.org/wiki/Great-circle_distance
  // a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlong/2)^2
  // c = 2*asin(sqrt(a))
  // d = R*c
  // R = 6371000

  float deg2rad;
  float a, c, R, distance;
  float dlat, dlong;

  deg2rad = 4.0*atan(1.0)/180.0;
  llat1 *= deg2rad;
  llat2 *= deg2rad;
  llong1 *= deg2rad;
  llong2 *= deg2rad;
  dlat = llat2 - llat1;
  dlong = llong2 - llong1;
  a = pow(sin(dlat / 2.0), 2) + cos(llat1)*cos(llat2)*pow(sin(dlong / 2.0), 2);
  c = 2 * asin(sqrt(ABS(a)));
  R = 6371000;
  distance = R*c;
  return distance;
}

/* ------------------ generate_fds ------------------------ */

void generate_fds(char *casename, int option){
  char buffer[LEN_BUFFER], elevfile[LEN_BUFFER], fdsfile[LEN_BUFFER],*ext;
  char basename[LEN_BUFFER];
  int nlong, nlat,nz;
  int i,j;
  float llat1, llat2, llong1, llong2;
  float xmax, ymax, zmin, zmax;
  float xmin_clip = -1.0, xmax_clip = -1.0;
  float ymin_clip = -1.0, ymax_clip = -1.0;
  int clip_vals = 0;
  float *xgrid, *ygrid;
  int count;
  int ibar, jbar, kbar;
  float **valptrs;
  FILE *streamin = NULL, *streamout = NULL;

  strcpy(basename, casename);
  ext = strrchr(basename, '.');
  if (ext != NULL)ext[0] = 0;

  strcpy(elevfile,basename);
  strcat(elevfile,"_elevs.csv");
  streamin = fopen(elevfile,"r");
  if(streamin==NULL){
    fprintf(stderr,"***error: unable to open %s for input\n",elevfile);
    return;
  }

  strcpy(fdsfile,basename);
  strcat(fdsfile,".fds");
  streamout = fopen(fdsfile,"w");
  if(streamout==NULL){
    fprintf(stderr,"***error: unable to open %s for output\n",fdsfile);
    fclose(streamin);
    return;
  }

  fgets(buffer, LEN_BUFFER, streamin);
  trim_back(buffer);
  sscanf(buffer, "%f %f %i %f %f %i %f %f %i %f %f %f %f",
    &llong1, &llong2,&nlong,
    &llat1, &llat2, &nlat,
    &zmin,&zmax,&nz,
    &xmin_clip, &xmax_clip, &ymin_clip, &ymax_clip
  );
  if (xmin_clip > -0.5&&xmax_clip > -0.5&&xmax_clip > xmin_clip&&
    ymin_clip > -0.5&&ymax_clip > -0.5&&ymax_clip > ymin_clip) {
    clip_vals = 1;
  }

  xmax = (int)(dist(llong1, llong2, llat1, llat1)+0.5);

  ymax = (int)(dist(llong1, llong1, llat1, llat2)+0.5);

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
    fprintf(streamout,"&MISC TERRAIN_CASE = .TRUE., TERRAIN_IMAGE = '%s.png' /\n", casename);
  }
  fprintf(streamout,"&TIME T_END = 0. /\n");
  fprintf(streamout,"&VENT XB = 0.0, 0.0, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", ymax, zmin, zmax);
  fprintf(streamout,"&VENT XB =  %f,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, xmax, ymax, zmin, zmax);
  fprintf(streamout,"&VENT XB = 0.0,  %f, 0.0, 0.0, %f, %f, SURF_ID = 'OPEN' /\n", xmax, zmin, zmax);
  fprintf(streamout,"&VENT XB = 0.0,  %f,  %f,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, ymax, ymax, zmin, zmax);
  fprintf(streamout,"&VENT XB = 0.0,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, ymax, zmax, zmax);
  fprintf(streamout,"&MATL ID = 'matl1', DENSITY = 1000., CONDUCTIVITY = 1., SPECIFIC_HEAT = 1., RGB = 122,117,48 /\n");
  fprintf(streamout,"&SURF ID = 'surf1', RGB = 122,117,48 TEXTURE_MAP='%s.png' /\n",casename);


  NewMemory((void **)&valptrs, sizeof(float *)*nlat);
  for(j = 0; j < nlat; j++){
    int idummy;
    float dummy, llat, llong, elev;
    float *vals;

    NewMemory((void **)&vals, sizeof(float)*nlong);
    valptrs[j] = vals;
    for(i = 0; i < nlong; i++){
      fgets(buffer, LEN_BUFFER, streamin);
      sscanf(buffer, "%i,%f,%f,%f,%f", &idummy, &llat, &llong, &dummy, &elev);
      vals[i] = elev;
    }
  }

  if(option==FDS_GEOM){
    fprintf(streamout,"&GEOM ID='terrain', SURF_ID='surf1',MATL_ID='matl1',\nIJK=%i,%i,XB=%f,%f,%f,%f,\nZVALS=\n",
                      nlong,nlat,0.0,xmax,0.0,ymax);
    count = 1;
    for(j = 0; j < jbar + 1; j++){
      float *vals;

      vals = valptrs[jbar - j];
      for(i = 0; i < ibar + 1; i++){
        fprintf(streamout," %f,", vals[i]);
        if(count % 10 == 0)fprintf(streamout,"\n");
        count++;
      }
    }
    fprintf(streamout,"/\n");
  }
  if(option==FDS_OBST){
    for(j = 0; j < jbar; j++){
      float *vals, *valsp1, ycen;

      vals = valptrs[j];
      valsp1 = valptrs[j+1];
      ycen = (ygrid[j] + ygrid[j + 1]) / 2.0;
      for(i = 0; i < ibar; i++){
        float vavg, xcen;

        xcen = (xgrid[i] + xgrid[i + 1]) / 2.0;
        if (clip_vals == 1 && xcen > xmin_clip&&xcen<xmax_clip&&ycen>ymin_clip&&ycen < ymax_clip)continue;
        vavg = (vals[i]+vals[i+1]+valsp1[i]+valsp1[i+1])/4.0;
        fprintf(streamout,"&OBST XB=%f,%f,%f,%f,0.0,%f SURF_ID='surf1'/\n", xgrid[i],xgrid[i+1],ygrid[j],ygrid[j+1],vavg);
      }
    }
  }
  fprintf(streamout,"&TAIL /\n");
}

/* ------------------ get_elevfile ------------------------ */

elevdata *get_elevfile(elevdata *elevinfo, int nelevinfo, float longval, float latval){
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

/* ------------------ get_elev ------------------------ */

float get_elevation(elevdata *elevinfo, int nelevinfo, float longval, float latval, int *have_val){
  elevdata *elevi;
  int index, ival, jval;
  float return_val;

  *have_val = 0;
  elevi = get_elevfile(elevinfo, nelevinfo, longval, latval);
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
  ival = CLAMP((longval - elevi->long_min)/elevi->cellsize,0,elevi->ncols-1);
  jval = CLAMP((elevi->lat_max - latval)/elevi->cellsize,0,elevi->nrows-1);
  index = jval*elevi->ncols +ival;
  return_val = elevi->valbuffer[index];
  *have_val = 1;
  return return_val;
}

/* ------------------ generate_elevs ------------------------ */

int generate_elevs(char *elevfile){
  int nelevinfo,i,j;
  filelistdata *fileheaders;
  FILE *stream_in, *stream_out;
  elevdata *elevinfo;
  int ibar, jbar, kbar;
  float dx, dy;
  float longc, latc;
  int longlat_defined = 0;
  float xmin_exclude, ymin_exclude, xmax_exclude, ymax_exclude;
  float longmin, longmax, latmin, latmax;
  int nlongs, nlats;
  float dlat, dlong;
  int count=0, *have_vals;
  float valmin, valmax, *vals;
  char outfile[LEN_BUFFER], *ext;

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

    if(match(buffer, "DXDY") == 1){
      dx = 1000.0;
      dy = 1000.0;
      if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f", &dx, &dy);
      continue;
    }

    if(match(buffer, "LONGMINMAX") == 1){
      longc = 1000.0;
      latc = 1000.0;
      if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f %i", &longmin, &longmax, &nlongs);
      continue;
    }

    if(match(buffer, "LATMINMAX") == 1){
      longc = 1000.0;
      latc = 1000.0;
      if(fgets(buffer, LEN_BUFFER, stream_in) == NULL)break;
      sscanf(buffer, "%f %f %i", &latmin, &latmax, &nlats);
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

  strcpy(outfile,elevfile);
  ext = strrchr(outfile, '.');
  if (ext != NULL)ext[0] = 0;
  strcat(outfile, "_elevs.csv");
  stream_out = fopen(outfile, "w");

  if(
    get_elevfile(elevinfo, nelevinfo, longmin, latmin) == NULL ||
    get_elevfile(elevinfo, nelevinfo, longmin, latmax) == NULL ||
    get_elevfile(elevinfo, nelevinfo, longmax, latmin) == NULL ||
    get_elevfile(elevinfo, nelevinfo, longmax, latmax) == NULL
    ){
    fprintf(stderr,"***error: elevation data not available for all longitudes/latitude \n");
    fprintf(stderr,"          pairs within the rectangle (%f,%f) (%f %f)\n", longmin, latmin, longmax, latmax);
    for(i = 0; i < nelevinfo; i++){
      elevdata *elevi;

      elevi = elevinfo + i;
      fprintf(stream_out," header file: %s bounds: %f %f %f %f\n",
        elevi->fileheader, elevi->long_min, elevi->lat_min, elevi->long_max, elevi->lat_max);
    }
    return 0;
  }

  dlat = (latmax - latmin) / (float)(nlats-1);
  dlong = (longmax - longmin) / (float)(nlongs-1);
  NewMemory((void **)&vals, nlongs*nlats*sizeof(float));
  NewMemory((void **)&have_vals, nlongs*nlats*sizeof(int));
  for(j = 0; j < nlats; j++){
    float latj;

    latj = latmin + (float)j*dlat;
    for(i = 0; i < nlongs; i++){
      float longi;
      float elevij;
      int have_val;

      longi = longmin + (float)i*dlong;

      elevij = get_elevation(elevinfo, nelevinfo, longi, latj, &have_val);
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
  count = 0;
  fprintf(stream_out," %f %f %i %f %f %i %f %f %i\n", longmin, longmax, nlongs, latmin, latmax, nlats,valmin,valmax,30);
  for(j = 0; j < nlats; j++){
    float latj;

    latj = latmin + (float)j*dlat;
    for(i = 0; i < nlongs; i++){
      float longi;

      longi = longmin + (float)i*dlong;

      if(have_vals[count] == 1){
        float elevij;

        elevij = vals[count];
        fprintf(stream_out, " %i,%f,%f,%f,%f\n", count+1, longi, latj, 39.37*elevij/12.0, elevij);
      }
      count++;
    }
  }
  fclose(stream_out);
  FREEMEMORY(have_vals);
  FREEMEMORY(vals);
  return 1;
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  int gen_fds = FDS_OBST;
  char *casename = NULL;
  char file_default[LEN_BUFFER];

  if(argc == 1){
    usage("dem2fds");
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
        show_example();
        exit(1);
        break;
      case 'h':
        usage("dem2fds");
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
        usage("dem2fds");
        exit(1);
        break;
      }
    }
    else{
      if(casename == NULL)casename = argv[i];
    }
  }
  if(casename == NULL)casename = file_default;
  if (generate_elevs(casename) == 1) {
    generate_fds(casename, gen_fds);
  }
  return 0;
}
