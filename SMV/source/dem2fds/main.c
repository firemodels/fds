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
  fprintf(stderr, "  dem2fds [-g|-o][-h][-v] casename.in\n");
  fprintf(stderr, "where\n");
  fprintf(stderr, "  -d dir - directory containing elevation files (default .)\n");
  fprintf(stderr, "  -g - create an FDS input file using &GEOM keywords\n");
  fprintf(stderr, "  -o - create an FDS input file using &OBST keywords (default)\n");
  fprintf(stderr, "  -h - display this message\n");
  fprintf(stderr, "  -v - show version information\n");
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
      switch(arg[1]){
      case 'd':
        i++;
        libdirptr = argv[i];
        if(file_exists(libdirptr) == 1){
          strcpy(libdir, libdirptr);
        }
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
