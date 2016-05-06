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

/* ------------------ usage ------------------------ */

void usage(char *prog){
 char githash[256];
 char gitdate[256];

  getGitInfo(githash,gitdate);    // get githash

  fprintf(stderr, "\n%s (%s) %s\n", prog, githash, __DATE__);
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
  float a, c, R, dist;
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
  dist = R*c;
  return dist;
}

/* ------------------ generate_elevs ------------------------ */

void generate_fds(char *filebase, int option){
  char buffer[LENBUFFER];
  int nlong, nlat,nz;
  int i,j;
  float llat1, llat2, llong1, llong2;
  float xmax, ymax, zmin, zmax, dx, dy;
  float *xgrid, *ygrid;
  int count;
  int ibar, jbar, kbar;
  float **valptrs;


  fgets(buffer, LENBUFFER, stdin);
  trim_back(buffer);
  sscanf(buffer, "%f %f %i %f %f %i %f %f %i", &llong1, &llong2,&nlong,&llat1, &llat2, &nlat, &zmin,&zmax,&nz);

  xmax = (int)(dist(llong1, llong2, llat1, llat1)+0.5);
  dx = xmax / (float)nlong;

  ymax = (int)(dist(llong1, llong1, llat1, llat2)+0.5);
  dy = ymax / (float)nlat;

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


  printf("&HEAD CHID='%s', TITLE='terrain' /\n",filebase);
  printf("&MESH IJK = %i, %i, %i, XB = 0.0, %f, 0.0, %f, %f, %f /\n",ibar,jbar,kbar,xmax,ymax,zmin,zmax);
  if(option==GENERATE_OBSTS){
    printf("&MISC TERRAIN_CASE = .TRUE., TERRAIN_IMAGE = '%s.png' /\n", filebase);
  }
  printf("&TIME T_END = 0. /\n");
  printf("&VENT XB = 0.0, 0.0, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", ymax, zmin, zmax);
  printf("&VENT XB =  %f,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, xmax, ymax, zmin, zmax);
  printf("&VENT XB = 0.0,  %f, 0.0, 0.0, %f, %f, SURF_ID = 'OPEN' /\n", xmax, zmin, zmax);
  printf("&VENT XB = 0.0,  %f,  %f,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, ymax, ymax, zmin, zmax);
  printf("&VENT XB = 0.0,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", xmax, ymax, zmax, zmax);
  printf("&MATL ID = 'matl1', DENSITY = 1000., CONDUCTIVITY = 1., SPECIFIC_HEAT = 1., RGB = 122,117,48 /\n");
  printf("&SURF ID = 'surf1', RGB = 122,117,48 TEXTURE_MAP='%s.png' /\n",filebase);


  NewMemory((void **)&valptrs, sizeof(float *)*nlat);
  for(j = 0; j < nlat; j++){
    int idummy;
    float dummy, llat, llong, elev;
    float *vals;

    NewMemory((void **)&vals, sizeof(float)*nlong);
    valptrs[j] = vals;
    for(i = 0; i < nlong; i++){
      fgets(buffer, LENBUFFER, stdin);
      sscanf(buffer, "%i,%f,%f,%f,%f", &idummy, &llat, &llong, &dummy, &elev);
      vals[i] = elev;
    }
  }

  if(option==GENERATE_GEOM){
    printf("&GEOM ID='terrain', SURF_ID='surf1',MATL_ID='matl1',\nIJK=%i,%i,XB=%f,%f,%f,%f,\nZVALS=\n",
                      nlong,nlat,0.0,xmax,0.0,ymax);
    count = 1;
    for(j = 0; j < jbar + 1; j++){
      float *vals;

      vals = valptrs[jbar - j];
      for(i = 0; i < ibar + 1; i++){
        printf(" %f,", vals[i]);
        if(count % 10 == 0)printf("\n");
        count++;
      }
    }
    printf("/\n");
  }
  if(option==GENERATE_OBSTS){
    for(j = 0; j < jbar; j++){
      float *vals, *valsp1;

      vals = valptrs[j];
      valsp1 = valptrs[j+1];
      for(i = 0; i < ibar; i++){
        float vavg;

        vavg = (vals[i]+vals[i+1]+valsp1[i]+valsp1[i+1])/4.0;
        printf("&OBST XB=%f,%f,%f,%f,0.0,%f SURF_ID='surf1'/\n", xgrid[i],xgrid[i+1],ygrid[j],ygrid[j+1],vavg);
      }
    }
  }
  printf("&TAIL /\n");
}

  /* ------------------ generate_latlongs ------------------------ */

void generate_longlats(char *filebase){
  char buffer[LENBUFFER];
  char fileout[LENBUFFER];
  float lat1, lat2, long1, long2;
  int nlat, nlong;
  int line_count, file_count;
  FILE *streamin = NULL, *streamout = NULL;
  int i;

  sprintf(fileout, "%s_longlats_%03i.csv", filebase, 1);

  fgets(buffer, LENBUFFER, stdin);
  sscanf(buffer, "%f %f %i %f %f %i", &long1, &long2, &nlong, &lat1, &lat2, &nlat);
  line_count = 1;
  file_count = 1;
  streamout = fopen(fileout, "w");
  for(i = 0; i < nlat; i++){
    int j;
    float llat;

    llat = (lat1*(float)(nlat - 1 - i) + lat2*(float)i) / (float)(nlat - 1);

    for(j = 0; j<nlong; j++){
      float llong;

      llong = (long1*(float)(nlong - 1 - j) + long2*(float)j) / (float)(nlong - 1);
      if(line_count>400){
        file_count++;
        fclose(streamout);
        sprintf(fileout, "%s_longlats_%03i.csv", filebase, file_count);
        streamout = fopen(fileout, "w");
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
  int gen_elevs = 0;
  int gen_obsts = 0;
  char *filebase = NULL;
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
        gen_obsts = 1;
        break;
      case 'e':
        gen_elevs = 1;
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
      if(filebase == NULL){
        filebase = argv[i];
      }
    }
  }
  if(filebase == NULL)filebase = file_default;
  if(gen_elevs == 1){
    generate_fds(filebase,GENERATE_GEOM);
  }
  else if(gen_obsts == 1){
    generate_fds(filebase,GENERATE_OBSTS);
  }
  else{
    generate_longlats(filebase);
  }
}
