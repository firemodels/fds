#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "string_util.h"
#include "file_util.h"
#include "MALLOC.h"

#define GENERATE_GEOM 1
#define GENERATE_OBSTS 2

int nelevsperfile=400;

/* ------------------ usage ------------------------ */

void usage(char *prog){
 char githash[256];
 char gitdate[256];

  getGitInfo(githash,gitdate);    // get githash

  fprintf(stderr, "\n%s (%s) %s\n", prog, githash, __DATE__);
  fprintf(stderr, "  Create an FDS input file using elevation data.\n\n");
  fprintf(stderr, "Usage:\n\n");
  fprintf(stderr, "  dem2fds [-e|-g|-o][-h][-n n][-v] casename\n\n");
  fprintf(stderr, "where\n\n");
  fprintf(stderr, "  -e - create elevation files (default)\n");
  fprintf(stderr, "  -g - create FDS input file from elevations using &GEOM keywords\n");
  fprintf(stderr, "  -o - create FDS input file from elevations using &OBST keywords\n");
  fprintf(stderr, "  -n n - number of elevations per file. default: %i \n",nelevsperfile);
  fprintf(stderr, "  -h - display this message\n");
  fprintf(stderr, "  -v - show versioning information\n\n");
  fprintf(stderr, " Usage examples:\n");
  fprintf(stderr, " 1. create elevation files:                dem2fds -e casename.in\n");
  fprintf(stderr, " 2. create FDS file using &GEOM keywords : dem2fds -g casename\n");
  fprintf(stderr, " 3. create FDS file using &OBST keywords : dem2fds -o casename\n");
}

#define ABS(a) ((a)>=0 ? (a) : (-(a)))
#define LENBUFFER 1024

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
  char buffer[LENBUFFER], elevfile[LENBUFFER], fdsfile[LENBUFFER];
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

  strcpy(elevfile,casename);
  strcat(elevfile,"_elevs.csv");
  streamin = fopen(elevfile,"r");
  if(streamin==NULL){
    printf("***error: unable to open %s for input\n",elevfile);
    return;
  }

  strcpy(fdsfile,casename);
  strcat(fdsfile,".fds");
  streamout = fopen(fdsfile,"w");
  if(streamout==NULL){
    printf("***error: unable to open %s for output\n",fdsfile);
    fclose(streamin);
    return;
  }
  
  fgets(buffer, LENBUFFER, streamin);
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


  fprintf(streamout,"&HEAD CHID='%s', TITLE='terrain' /\n",casename);
  fprintf(streamout,"&MESH IJK = %i, %i, %i, XB = 0.0, %f, 0.0, %f, %f, %f /\n",ibar,jbar,kbar,xmax,ymax,zmin,zmax);
  if(option==GENERATE_OBSTS){
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
      fgets(buffer, LENBUFFER, streamin);
      sscanf(buffer, "%i,%f,%f,%f,%f", &idummy, &llat, &llong, &dummy, &elev);
      vals[i] = elev;
    }
  }

  if(option==GENERATE_GEOM){
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
  if(option==GENERATE_OBSTS){
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

  /* ------------------ generate_latlongs ------------------------ */

void generate_longlats(char *elevfile){
  char buffer[LENBUFFER], casename[LENBUFFER], *ext;
  char fileout[LENBUFFER];
  float lat1, lat2, long1, long2;
  int nlat, nlong;
  int line_count, file_count;
  FILE *streamin=NULL, *streamout = NULL;
  int i;

  strcpy(casename, elevfile);
  ext=strrchr(casename, '.');
  if (ext != NULL)ext[0] = 0;

  streamin = fopen(elevfile, "r");
  if (streamin == NULL) {
    printf("***error: unable to open %s for input\n", elevfile);
    return;
  }

  fgets(buffer, LENBUFFER, streamin);
  sscanf(buffer, "%f %f %i %f %f %i", &long1, &long2, &nlong, &lat1, &lat2, &nlat);
  if (nlat < 2)nlat = 2;
  if (nlong < 2)nlong = 2;
  line_count = 1;
  file_count = 1;

  sprintf(fileout, "%s_longlats_%03i.csv", casename, file_count);
  streamout = fopen(fileout, "w");
  if (streamout == NULL) {
    printf("***error: unable to open %s for output\n", fileout);
    fclose(streamin);
    return;
  }
  for(i = 0; i < nlat; i++){
    int j;
    float llat;

    llat = (lat1*(float)(nlat - 1 - i) + lat2*(float)i) / (float)(nlat - 1);
    for(j = 0; j<nlong; j++){
      float llong;

      llong = (long1*(float)(nlong - 1 - j) + long2*(float)j) / (float)(nlong - 1);
      if(line_count>nelevsperfile){
        file_count++;
        fclose(streamout);
        sprintf(fileout, "%s_longlats_%03i.csv", casename, file_count);
        streamout = fopen(fileout, "w");
        if (streamout == NULL) {
          printf("***error: unable to open %s for output\n", fileout);
          fclose(streamin);
          return;
        }
        line_count = 1;
      }
      fprintf(streamout, "%f,%f\n", llong, llat);
      line_count++;
    }
  }
}

/* ------------------ main ------------------------ */


int main(int argc, char **argv){
  int i;
  int gen_fdsgeom = 0;
  int gen_fdsobst = 0;
  char *casename = NULL;
  char file_default[1000];

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
      case 'h':
        usage("dem2geom");
        exit(1);
        break;
      case 'o':
        gen_fdsobst = 1;
        break;
      case 'e':
        break;
      case 'n':
        ++i;
        if(i < argc)sscanf(argv[i], "%i", &nelevsperfile);
        printf("nelevsperfile=%i\n",nelevsperfile);
        break;
      case 'g':
        gen_fdsgeom = 1;
        break;
      case 'v':
        PRINTversion("dem2geom");
        exit(1);
        break;
      default:
        usage("dem2geom");
        exit(1);
        break;
      }
    }
    else{
      if(casename == NULL)casename = argv[i];
    }
  }
  if(casename == NULL)casename = file_default;
  if(gen_fdsgeom == 1){
    generate_fds(casename,GENERATE_GEOM);
  }
  else if(gen_fdsobst == 1){
    generate_fds(casename,GENERATE_OBSTS);
  }
  else{
    generate_longlats(casename);
  }
  return 0;
}
