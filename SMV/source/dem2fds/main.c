#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "string_util.h"
#include "file_util.h"


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

void generate_elevs(void){
  char buffer[LENBUFFER];
  int nlong, nlat;
  int i;
  float llat1, llat2, llong1, llong2;
  float deltax, deltay, zmin, zmax;


  fgets(buffer, LENBUFFER, stdin);
  sscanf(buffer, "%i %i", &nlong, &nlat);
  fgets(buffer, LENBUFFER, stdin);
  sscanf(buffer, "%f %f %f %f %f %f", &llong1, &llong2,&llat1, &llat2, &zmin,&zmax);

  deltax = (int)(dist(llong1, llong2, llat1, llat1)+0.5);
  deltay = (int)(dist(llong1, llong1, llat1, llat2)+0.5);

  printf("&HEAD CHID='terrain', TITLE='terrain' /\n");
  printf("&MESH IJK = %i, %i, %i, XB = 0.0, %f, 0.0, %f, %f, %f /\n",nlong,nlat,30,deltax,deltay,zmin,zmax);
  printf("&TIME T_END = 0. /\n");
  printf("&VENT XB = 0.0, 0.0, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", deltay, zmin,   zmax);
  printf("&VENT XB =  %f,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", deltax, deltax, deltay, zmin, zmax);
  printf("&VENT XB = 0.0,  %f, 0.0, 0.0, %f, %f, SURF_ID = 'OPEN' /\n", deltax, zmin,   zmax);
  printf("&VENT XB = 0.0,  %f,  %f,  %f, %f, %f, SURF_ID = 'OPEN' /\n", deltax, deltay, deltay, zmin, zmax);
  printf("&VENT XB = 0.0,  %f, 0.0,  %f, %f, %f, SURF_ID = 'OPEN' /\n", deltax, deltay,   zmax, zmax);

  printf("&GEOM ID='terrain', VERTS=\n");
  for(i = 0; i < nlong; i++){
    int j;
    float xx;

    xx = 0.0*(float)(nlong - i) + deltax*(float)i;
    xx /= (float)(nlong - 1);

    for(j = 0; j < nlat; j++){
      int idummy;
      float dummy, llat, llong, elev;
      float yy;

      yy = 0.0*(float)(nlat - j) + deltay*(float)j;
      yy /= (float)(nlat - 1);

      fgets(buffer, LENBUFFER, stdin);
      sscanf(buffer, "%i,%f,%f,%f,%f", &idummy, &llat, &llong, &dummy, &elev);
      printf(" %f,%f,%f,\n", xx,yy, elev);
    }
  }
#define IJ(i,j) (nlat*(j-1)+i)
  printf(" FACES=\n");
  for(i = 1; i < nlong; i++){
    int j;

    for(j = 1; j < nlat; j++){
      //  j+1
      //  j
      //     i     i+1
      printf(" %i, %i, %i, %i, %i, %i,\n", IJ(i, j), IJ(i + 1, j+1), IJ(i + 1, j ), IJ(i, j), IJ(i , j + 1), IJ(i+1, j + 1));
    }
  }
  printf("/\n");
  printf("&TAIL /\n");

}

  /* ------------------ generate_latlongs ------------------------ */

void generate_latlongs(void){
  char buffer[LENBUFFER];
  char filebase[LENBUFFER], fileout[LENBUFFER];
  float lat1, lat2, long1, long2;
  int nlat, nlong;
  int line_count, file_count;
  FILE *streamin = NULL, *streamout = NULL;
  int i;

  strcpy(filebase, "elevations");
  sprintf(fileout, "%s%i.csv", filebase, 1);

  fgets(buffer, LENBUFFER, stdin);
  sscanf(buffer, "%f %f %i %f %f %i", &lat1, &lat2, &nlat, &long1, &long2, &nlong);
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
        sprintf(fileout, "%s%i.csv", filebase, file_count);
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
        break;
      case 'e':
        gen_elevs = 1;
        break;
      case 'v':
        version("dem2geom");
        exit(1);
        break;
      default:
        usage("dem2geom");
        exit(1);
        break;
      }
    }
  }
  if(gen_elevs == 1){
    generate_elevs();
  }
  else{
    generate_latlongs();
  }
}
