// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_GPU
#include <GL/glew.h>
#endif
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "contourdefs.h"
#include "isodefs.h"

#include "flowfiles.h"
#include "smokeviewapi.h"
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char unit_revision[]="$Revision$";

/* ------------------ get_unit_class ------------------------ */

f_units *get_unit_class(char *label){
  int i;

  for(i=0;i<nunitclasses;i++){
    f_units *unitclass;
    int j;

    unitclass = unitclasses + i;
    for(j=0;j<unitclass->nunittypes;j++){
      if(strcmp(label,unitclass->unittypes[j])==0)return unitclass;
    }
  }
  return NULL;
}

/* ------------------ unit_type_match ------------------------ */

int unit_type_match(char *unit_type, f_units *unit_class){
  int i;

  for(i=0;i<unit_class->nunittypes;i++){
    if(strcmp(unit_class->unittypes[i],unit_type)==0)return 1;
  }
  return 0;
}

/* ------------------ add_unit_class ------------------------ */

void add_unit_class(flowlabels *label){
  f_units *ut;
  f_unit *units;
  int i;

  nunitclasses++;
  ResizeMemory((void **)&unitclasses,nunitclasses*sizeof(f_units));

  ut = unitclasses + nunitclasses - 1;

  ut->submenuid=0;

  ut->diff_index=1;
  ut->nunits=2;
  ut->active=0;
  strcpy(ut->unitclass,label->shortlabel);
  ut->nunittypes=1;
  strcpy(ut->unittypes[0],label->shortlabel);

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,label->unit);
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;

  strcpy(units[1].unit,"% diff");
  units[1].scale[0]=0.1667;
  units[1].scale[1]=0.0;
  units[1].rel_defined=0;

}

/* ------------------ init_unit_defs ------------------------ */

void init_unit_defs(void){
  int j;
  f_units *unitclass;

  if(unitclasses_ini!=NULL){
    unitclasses=unitclasses_ini;
    nunitclasses=nunitclasses_ini;
  }
  else{
    unitclasses=unitclasses_default;
    nunitclasses=nunitclasses_default;
  }
  if(smokediff==0)return;
  for(j=0;j<nslice_files;j++){
    slice *slicej;
    char *shortlabel;
      
    slicej = sliceinfo + j;
    shortlabel = slicej->label.shortlabel;
    unitclass = get_unit_class(shortlabel);
    if(unitclass==NULL){
      add_unit_class(&slicej->label);
    }
  }
  for(j=0;j<npatch_files;j++){
    patch *patchj;
    char *shortlabel;
      
    patchj = patchinfo + j;
    shortlabel = patchj->label.shortlabel;
    unitclass = get_unit_class(shortlabel);
    if(unitclass==NULL){
      add_unit_class(&patchj->label);
    }
  }
  for(j=0;j<nplot3d_files;j++){
    plot3d *plot3dj;
    char *shortlabel;
    int n;
      
    plot3dj = plot3dinfo + j;
    for(n=0;n<5;n++){
      shortlabel = plot3dj->label[n].shortlabel;
      unitclass = get_unit_class(shortlabel);
      if(unitclass==NULL){
        add_unit_class(&plot3dj->label[n]);
      }
    }
  }
  CheckMemory;
}

/* ------------------ update_unit_defs ------------------------ */

void update_unit_defs(void){
  int i, j;

  if(smokediff==0)return;
  for(i=0;i<nunitclasses;i++){
    float valmin, valmax, diff_maxmin;
    int firstslice, firstpatch, firstplot3d, diff_index;

    firstpatch=1;
    for(j=0;j<npatch_files;j++){
      patch *patchj;
      
      patchj = patchinfo + j;
      if(patchj->loaded==0||patchj->display==0)continue;
      if(unit_type_match(patchj->label.shortlabel,unitclasses+i)==0)continue;
      if(firstpatch==1){
        firstpatch=0;
        valmin=patchj->diff_valmin;
        valmax=patchj->diff_valmax;
      }
      else{
        if(patchj->diff_valmin<valmin)valmin=patchj->diff_valmin;
        if(patchj->diff_valmax>valmax)valmax=patchj->diff_valmax;
      }
    }

    firstslice=1;
    for(j=0;j<nslice_files;j++){
      slice *slicej;
      
      slicej = sliceinfo + j;
      if(slicej->loaded==0||slicej->display==0)continue;
      if(unit_type_match(slicej->label.shortlabel,unitclasses+i)==0)continue;
      if(firstslice==1){
        firstslice=0;
        valmin=slicej->diff_valmin;
        valmax=slicej->diff_valmax;
      }
      else{
        if(slicej->diff_valmin<valmin)valmin=slicej->diff_valmin;
        if(slicej->diff_valmax>valmax)valmax=slicej->diff_valmax;
      }
    }

    firstplot3d=1;
    for(j=0;j<nplot3d_files;j++){
      plot3d *plot3dj;
      int n;
      
      plot3dj = plot3dinfo + j;
      if(plot3dj->loaded==0||plot3dj->display==0)continue;
      for(n=0;n<5;n++){
        char *shortlabel;

        shortlabel = plot3dj->label[n].shortlabel;
        if(unit_type_match(shortlabel,unitclasses+i)==0)continue;
        if(firstplot3d==1){
          firstplot3d=0;
          valmin=plot3dj->diff_valmin[n];
          valmax=plot3dj->diff_valmax[n];
        }
         else{
          if(plot3dj->diff_valmin[n]<valmin)valmin=plot3dj->diff_valmin[n];
          if(plot3dj->diff_valmax[n]>valmax)valmax=plot3dj->diff_valmax[n];
        }
      }
    }

    diff_index=unitclasses[i].diff_index;
    if(diff_index!=-1&&(firstslice==0||firstpatch==0||firstplot3d==0)){
      int idiff;

      diff_maxmin=valmax-valmin;
      if(diff_maxmin!=0.0){
        unitclasses[i].units[diff_index].scale[0]=100.0/diff_maxmin;

        idiff = diff_maxmin + 0.5;
        diff_maxmin=idiff;
        sprintf(unitclasses[i].units[diff_index].rel_val,"%f",diff_maxmin);
        trimzeros(unitclasses[i].units[diff_index].rel_val);
        unitclasses[i].units[diff_index].rel_defined=1;
      }
    }
  }

}

/* ------------------ getunitinfo ------------------------ */

void getunitinfo(const char *shortlabel, int *unitclass, int *unittype){
  int i,j;
  *unitclass=-1;
  *unittype=-1;
  for(i=0;i<nunitclasses;i++){
    for(j=0;j<unitclasses[i].nunittypes;j++){
      if(STRCMP(unitclasses[i].unittypes[j],shortlabel)==0){
        *unitclass=i;
        *unittype=unitclasses[i].active;
        return;
      }
    }
  }
}

/* ------------------ Init_Units ------------------------ */

void InitUnits(void){
  int i;
  f_unit *units;
  f_units *ut;

  nunitclasses_default=2;
  NewMemory((void **)&unitclasses_default,nunitclasses_default*sizeof(f_units));
  
  for(i=0;i<nunitclasses_default;i++){
    unitclasses_default[i].submenuid=0;
  }
  ut=unitclasses_default;

  // temperature units

  if(smokediff==1){
    ut->diff_index=3;
    ut->nunits=4;
  }
  else{
    ut->diff_index=-1;
    ut->nunits=3;
  }
  ut->active=0;
  strcpy(ut->unitclass,"temperature");
  ut->nunittypes=2;
  strcpy(ut->unittypes[0],"temp");
  strcpy(ut->unittypes[1],"d_temp");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"C");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"F");
  units[1].scale[0]=1.8;
  if(smokediff==1){
    units[1].scale[1]=0.0;
    units[1].rel_defined=0;
  }
  else{
    units[1].scale[1]=32.0;
  }
  strcpy(units[2].unit,"K");
  units[2].scale[0]=1.0;
  if(smokediff==1){
    units[2].scale[1]=273.15;
  }
  else{
    units[2].scale[1]=0.0;
  }
  if(smokediff==1){
    strcpy(units[3].unit,"% diff");
    units[3].scale[0]=0.1667;
    units[3].scale[1]=0.0;
    units[3].rel_defined=0;
  }

  // velocity units

  ut = unitclasses_default+1;
  ut->active=0;
  if(smokediff==1){
    ut->diff_index=3;
    ut->nunits=4;
  }
  else{
    ut->diff_index=-1;
    ut->nunits=3;
  }
  strcpy(ut->unitclass,"velocity");
  ut->nunittypes=4;
  strcpy(ut->unittypes[0],"U-VEL");
  strcpy(ut->unittypes[1],"V-VEL");
  strcpy(ut->unittypes[2],"W-VEL");
  strcpy(ut->unittypes[3],"Speed");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"m/s");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"mph");
  units[1].scale[0]=2.236931818;
  units[1].scale[1]=0.0;
  strcpy(units[2].unit,"f/s");
  units[2].scale[0]=3.2808333333;
  units[2].scale[1]=0.0;
  if(smokediff==1){
    strcpy(units[3].unit,"% diff");
    units[3].scale[0]=0.16667;
    units[3].scale[1]=0.0;
    units[3].rel_defined=0;
  }
  CheckMemory;
}
