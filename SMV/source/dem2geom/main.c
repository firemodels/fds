#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "string_util.h"
#include "file_util.h"


/* ------------------ usage ------------------------ */

void usage(char *prog){
 char githash[256];
 char gitdate[256];

  getGitInfo(githash,gitdate);    // get githash

  fprintf(stderr, "\n%s (%s) %s\n", prog, githash, __DATE__);
}

/* ------------------ main ------------------------ */

#define LENBUFFER 1024

int main(int argc, char **argv){
  char buffer[LENBUFFER],*buffptr;
  int i;
  char filebase[LENBUFFER], fileout[LENBUFFER];
  FILE *streamin=NULL,*streamout=NULL;
  float lat1, lat2, long1, long2;
  int nlat, nlong;
  int line_count,file_count;

  set_stdout(stdout);
  buffptr=buffer;
  strcpy(filebase, "elevations");
  sprintf(fileout, "%s%i", filebase,1);
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
      case 'v':
        version("dem2geom");
        exit(1);
        break;
      default:
        usage("demo2geom");
        exit(1);
        break;
      }
    }
  }
  fgets(buffer, LENBUFFER, stdin);
  sscanf(buffer, "%f %f %i %f %f %i", &lat1, &lat2, &nlat, &long1, &long2, &nlong);
  line_count = 0;
  file_count = 1;
  streamout=fopen(fileout, "w");
  for(i = 0; i<nlat; i++){
    int j;
    float llat;

    llat = (lat1*(float)(nlat-1-i)+lat2*(float)i)/(float)(nlat-1);

    for(j = 0; j<nlong; j++){
      float llong;

      llong = (long1*(float)(nlong-1-j)+long2*(float)j)/(float)(nlong-1);
      if(line_count>500){
        file_count++;
        fclose(streamout);
        sprintf(fileout, "%s%i", filebase, file_count);
        streamout = fopen(fileout, "w");
        line_count = 0;
      }
      fprintf(streamout,"%f,%f\n", llat, llong);
      line_count++;

    }
  }
}
