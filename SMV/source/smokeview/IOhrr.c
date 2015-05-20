// $Date: 2014-07-24 13:03:47 -0400 (Thu, 24 Jul 2014) $ 
// $Revision: 20000 $
// $Author: gforney $

// svn revision character string
char IOhrr_revision[]="$Revision: 20000 $";

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include GLUT_H

#include "smokeviewvars.h"

/* ------------------ readhrr ------------------------ */
#define LENBUFFER 1024
void readhrr(int flag, int *errorcode){
  FILE *HRRFILE;
  int ntimeshrr, nfirst;
  char buffer[LENBUFFER];
  float *hrrtime, *hrrval;
  int display=0;
  int ntimes_saved;

  *errorcode=0;
  if(hrrinfo!=NULL){
    display = hrrinfo->display;
    FREEMEMORY(hrrinfo->times_csv);
    FREEMEMORY(hrrinfo->times);
    FREEMEMORY(hrrinfo->hrrval_csv);
    FREEMEMORY(hrrinfo->hrrval);
    FREEMEMORY(hrrinfo->timeslist);
  }
  FREEMEMORY(hrrinfo);
  if(flag==UNLOAD)return;

  NewMemory((void **)&hrrinfo,sizeof(hrrdata));
  hrrinfo->file=hrr_csv_filename;
  hrrinfo->times_csv=NULL;
  hrrinfo->times=NULL;
  hrrinfo->timeslist=NULL;
  hrrinfo->hrrval_csv=NULL;
  hrrinfo->hrrval=NULL;
  hrrinfo->ntimes_csv=0;
  hrrinfo->loaded=1;
  hrrinfo->display=display;
  hrrinfo->itime=0;

  HRRFILE=fopen(hrrinfo->file,"r");
  if(HRRFILE==NULL){
    readhrr(UNLOAD,errorcode);
    return;
  }


// size data

  ntimeshrr=0;
  nfirst=-1;
  while(!feof(HRRFILE)){
    if(fgets(buffer,LENBUFFER,HRRFILE)==NULL)break;
    if(nfirst==-1&&strstr(buffer,".")!=NULL)nfirst=ntimeshrr;
    ntimeshrr++;
  }
  ntimes_saved=ntimeshrr;
  ntimeshrr-=nfirst;

  rewind(HRRFILE);
  NewMemory((void **)&hrrinfo->times_csv,ntimeshrr*sizeof(float));
  NewMemory((void **)&hrrinfo->hrrval_csv,ntimeshrr*sizeof(float));

// read data
  
  hrrtime=hrrinfo->times_csv;
  hrrval=hrrinfo->hrrval_csv;
  ntimeshrr=0;

// read no more than the number of lines found during first pass

  while(ntimeshrr<ntimes_saved&&!feof(HRRFILE)){
    if(fgets(buffer,LENBUFFER,HRRFILE)==NULL)break;
    if(ntimeshrr<nfirst){
      ntimeshrr++;
      continue;
    }
    stripcommas(buffer);
    sscanf(buffer,"%f %f",hrrtime,hrrval);
    hrrtime++;
    hrrval++;
    ntimeshrr++;
  }
  hrrinfo->ntimes_csv=ntimeshrr-nfirst;
  fclose(HRRFILE);
}
