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
#include "gd.h"
#include "dem_util.h"

/* ------------------ Usage ------------------------ */

void Usage(char *prog){
 char githash[LEN_BUFFER];
 char gitdate[LEN_BUFFER];

  getGitInfo(githash,gitdate);    // get githash

  fprintf(stderr, "\n%s (%s) %s\n", prog, githash, __DATE__);
  fprintf(stderr, "Create an FDS input file using elevation data\n");
  fprintf(stderr, "  obtained from http://viewer.nationalmap.gov \n\n");
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "  dem2fds [-dir dir][-geom|-obst][-help][-nobuffer][-version] casename.in\n");
  fprintf(stderr, "where\n");
  fprintf(stderr, "  -dir dir  - directory containing elevation and map files (default .)\n");
  fprintf(stderr, "  -geom     - create an FDS input file using &GEOM keywords\n");
  fprintf(stderr, "  -help     - display this message\n");
  fprintf(stderr, "  -nobuffer - create a terrain map assuming no buffer exists between maps.\n");
  fprintf(stderr, "              Otherwise assume that a 300 pixel buffer exists bewteen maps.\n");
  fprintf(stderr, "  -obst     - create an FDS input file using &OBST keywords\n");
  fprintf(stderr, "  -version  - show version information\n");
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
  strcpy(libdir, ".");

  initMALLOC();
  set_stdout(stdout);
  for(i = 1; i<argc; i++){
    int lenarg;
    char *arg,*libdirptr;


    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      if(strncmp(arg, "-dir", 4) == 0){
        i++;
        libdirptr = argv[i];
        if(file_exists(libdirptr) == 1){
          strcpy(libdir, libdirptr);
        }
      }
      else if(strncmp(arg, "-show", 5)==0|| strncmp(arg, "-s", 2) == 0){
        show_maps = 1;
      }
      else if(strncmp(arg, "-help", 5) == 0|| strncmp(arg, "-h", 2) == 0){
        Usage("dem2fds");
        return 1;
      }
      else if(strncmp(arg, "-nobuffer", 8) == 0|| strncmp(arg, "-n", 2) == 0){
        border_buffer = 0;
      }
      else if(strncmp(arg, "-obst", 5) == 0|| strncmp(arg, "-o", 2) == 0){
        gen_fds = FDS_OBST;
      }
      else if(strncmp(arg, "-geom", 5) == 0|| strncmp(arg, "-g", 2) == 0){
        gen_fds = FDS_GEOM;
      }
      else if(strncmp(arg, "-version", 8) == 0|| strncmp(arg, "-v", 2) == 0){
        PRINTversion("dem2fds");
        return 1;
      }
      else{
        Usage("dem2fds");
        return 1;
      }
    }
    else{
      if(casename == NULL)casename = argv[i];
    }
  }
  if(casename == NULL)casename = file_default;
  if (GetElevations(casename,&fds_elevs) == 1) {
    GenerateFDSInputFile(casename, &fds_elevs, gen_fds);
  }
  return 0;
}
