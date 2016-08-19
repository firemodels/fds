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

  GetGitInfo(githash,gitdate);    // get githash

  fprintf(stdout, "\n%s (%s) %s\n", prog, githash, __DATE__);
  fprintf(stdout, "Create an FDS input file using elevation and image\n");
  fprintf(stdout, "  data obtained from http://viewer.nationalmap.gov \n\n");
  fprintf(stdout, "Usage:\n");
  fprintf(stdout, "  dem2fds [options] casename.in\n");
  fprintf(stdout, "  -dir dir  - directory containing map files and elevation files if -elevdir \n");
  fprintf(stdout, "              is not specified (default: '.')\n");
  fprintf(stdout, "  -elevdir dir  - directory containing elevation files\n");
  fprintf(stdout, "  -elevs    - only output elevations, do not create a complete FDS input file\n");
  fprintf(stdout, "  -geom     - create an FDS input file using &GEOM keywords (experimental)\n");
  fprintf(stdout, "  -help     - display this message\n");
  fprintf(stdout, "  -nobuffer - create a terrain map assuming no buffer exists between maps.\n");
  fprintf(stdout, "              Otherwise assume that a 300 pixel buffer exists between maps.\n");
  fprintf(stdout, "  -obst     - create an FDS input file using &OBST keywords\n");
  fprintf(stdout, "  -show     - show image boundaries (black outline) and \n");
  fprintf(stdout, "              fds scenario boundary (red outline)\n");
  fprintf(stdout, "  -surf surf_id - specify surf ID for use by OBSTs or geometry \n");
  fprintf(stdout, "  -version  - show version information\n");
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  int gen_fds = FDS_OBST;
  char *casename = NULL;
  char file_default[LEN_BUFFER];
  elevdata fds_elevs;
  int fatal_error = 0;

  if(argc == 1){
    Usage("dem2fds");
    return 0;
  }

  strcpy(file_default, "terrain");
  strcpy(image_dir, ".");
  strcpy(elev_dir, "");
  strcpy(surf_id, "surf1");

  initMALLOC();
  set_stdout(stdout);

  for(i = 1; i<argc; i++){
    int lenarg;
    char *arg;


    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      if(strncmp(arg, "-dir", 4) == 0){
        i++;
        if(file_exists(argv[i]) != 1)fatal_error = 1;
      }
      else if (strncmp(arg, "-elevdir", 8) == 0) {
        i++;
        if (file_exists(argv[i]) != 1)fatal_error = 1;
      }
    }
    else{
      if(casename == NULL)casename = argv[i];
    }
  }

  if(casename==NULL){
    fprintf(stderr, "\n***error: input file not specified\n");
    return 1;
  }
  if(file_exists(casename)!=1){
    fprintf(stderr, "\n***error: input file %s does not exist\n",casename);
    return 1;
  }
  if (fatal_error == 1) {
    fprintf(stderr, "\ncase: %s\n", casename);
  }

  fatal_error=0;
  for(i = 1; i<argc; i++){
    int lenarg;
    char *arg;


    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      if(strncmp(arg, "-dir", 4) == 0){
        i++;
        if(file_exists(argv[i]) == 1){
          strcpy(image_dir, argv[i]);
          if (strlen(elev_dir) == 0) {
            strcpy(elev_dir, image_dir);
          }
        }
        else {
          fprintf(stderr, "***error: directory %s does not exist or cannot be accessed\n",argv[i]);
          fatal_error = 1;
        }
      }
      else if (strncmp(arg, "-elevdir", 8) == 0) {
        i++;
        if (file_exists(argv[i]) == 1) {
          strcpy(elev_dir, argv[i]);
        }
        else {
          fprintf(stderr, "***error: directory %s does not exist or cannot be accessed\n", argv[i]);
          fatal_error = 1;
        }
      }
      else if(strncmp(arg, "-elevs", 6) == 0 ) {
        elev_file = 1;
      }
      else if(strncmp(arg, "-geom", 5) == 0 || strncmp(arg, "-g", 2) == 0){
        gen_fds = FDS_GEOM;
      }
      else if(strncmp(arg, "-help", 5) == 0 || strncmp(arg, "-h", 2) == 0){
        Usage("dem2fds");
        return 1;
      }
      else if(strncmp(arg, "-nobuffer", 8) == 0 || strncmp(arg, "-n", 2) == 0){
        border_buffer = 0;
      }
      else if(strncmp(arg, "-obst", 5) == 0 || strncmp(arg, "-o", 2) == 0){
        gen_fds = FDS_OBST;
      }
      else if(strncmp(arg, "-show", 5) == 0){
        show_maps = 1;
      }
      else if(strncmp(arg, "-surf", 5) == 0){
        i++;
        strcpy(surf_id, argv[i]);
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

  if (fatal_error == 1) return 1;

  if (strlen(elev_dir) == 0) {
    strcpy(elev_dir, image_dir);
  }
  if(casename == NULL)casename = file_default;
  if(GetElevations(casename,&fds_elevs) == 1) {
    GenerateFDSInputFile(casename, &fds_elevs, gen_fds);
  }
  return 0;
}
