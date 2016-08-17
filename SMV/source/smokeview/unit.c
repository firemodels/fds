#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include GLUT_H

#include "smokeviewvars.h"

/* ------------------ GetUnitClass ------------------------ */

f_units *GetUnitClass(char *unit){
  int i;

  for(i=0;i<nunitclasses;i++){
    f_units *unitclass;

    unitclass = unitclasses + i;
    if(strlen(unit)!=strlen(unitclass->units->unit))continue;
    if(STRCMP(unit,unitclass->units->unit)==0)return unitclass;
  }
  return NULL;
}

/* ------------------ UnitTypeMatch ------------------------ */

int UnitTypeMatch(char *unit, f_units *unit_class){

  if(strlen(unit_class->units->unit)!=strlen(unit))return 1;
  if(strcmp(unit_class->units->unit,unit)!=0)return 1;
  return 0;
}

/* ------------------ IsUnitPresent ------------------------ */

int IsUnitPresent(char *label, char *unit){
  if(strlen(label)!=strlen(unit)||STRCMP(label,unit)!=0)return 0;
  return 1;
}

/* ------------------ SetUnitVis ------------------------ */

void SetUnitVis(void){
  int i;
  int j;

  for(i=0;i<nunitclasses;i++){
    f_units *uci;

    uci = unitclasses + i;
    uci->visible=0;

    for(j=0;j<nsliceinfo;j++){
      slicedata *slicej;

      slicej = sliceinfo + j;
      if(IsUnitPresent(slicej->label.unit,uci->units->unit)==1){
        uci->visible=1;
        break;
      }
    }
    if(uci->visible==1)continue;

    for(j=0;j<npatchinfo;j++){
      patchdata *patchj;

      patchj = patchinfo + j;
      if(IsUnitPresent(patchj->label.unit,uci->units->unit)==1){
        uci->visible=1;
        break;
      }
    }
    if(uci->visible==1)continue;

    for(j=0;j<nplot3dinfo;j++){
      plot3ddata *plot3dj;
      int n;

      plot3dj = plot3dinfo + j;
      for(n=0;n<5;n++){
        if(IsUnitPresent(plot3dj->label[n].unit,uci->units->unit)==1){
          uci->visible=1;
          break;
        }
      }
      if(uci->visible==1)break;
    }
  }
}

/* ------------------ InitUnitDefs ------------------------ */

void InitUnitDefs(void){
  if(unitclasses_ini!=NULL){
    unitclasses=unitclasses_ini;
    nunitclasses=nunitclasses_ini;
  }
  else{
    unitclasses=unitclasses_default;
    nunitclasses=nunitclasses_default;
  }
  CheckMemory;
}

/* ------------------ UpdateUnitDefs ------------------------ */

void UpdateUnitDefs(void){
  int i, j;

  if(smokediff==0)return;
  for(i=0;i<nunitclasses;i++){
    float valmin, valmax, diff_maxmin;
    int firstslice, firstpatch, firstplot3d, diff_index;

    firstpatch=1;
    for(j=0;j<npatchinfo;j++){
      patchdata *patchj;

      patchj = patchinfo + j;
      if(patchj->loaded==0||patchj->display==0)continue;
      if(UnitTypeMatch(patchj->label.unit,unitclasses+i)!=0)continue;
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
    for(j=0;j<nsliceinfo;j++){
      slicedata *slicej;

      slicej = sliceinfo + j;
      if(slicej->loaded==0||slicej->display==0)continue;
      if(UnitTypeMatch(slicej->label.unit,unitclasses+i)!=0)continue;
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
    for(j=0;j<nplot3dinfo;j++){
      plot3ddata *plot3dj;
      int n;

      plot3dj = plot3dinfo + j;
      if(plot3dj->loaded==0||plot3dj->display==0)continue;
      for(n=0;n<5;n++){
        if(UnitTypeMatch(plot3dj->label[n].unit,unitclasses+i)!=0)continue;
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
        TrimZeros(unitclasses[i].units[diff_index].rel_val);
        unitclasses[i].units[diff_index].rel_defined=1;
      }
    }
  }

}

/* ------------------ GetUnitVal ------------------------ */

float GetUnitVal(const char *unittype, float oldval){
  int i;

  for(i=0;i<nunitclasses;i++){
    if(STRCMP(unitclasses[i].unitclass,unittype)==0){
      int unit_index;
      float *scale;
      float val;

      unit_index = unitclasses[i].unit_index;
      scale = unitclasses[i].units[unit_index].scale;
      val =  scale[0]*oldval + scale[1];
      val = (float)((int)(val*10.0 + 0.5))/10.0;
      return val;
    }
  }
  return oldval;
}

/* ------------------ GetUnitInfo ------------------------ */

void GetUnitInfo(const char *unitlabel, int *unitclass, int *unittype){
  int i;

  *unitclass=-1;
  *unittype=-1;
  for(i=0;i<nunitclasses;i++){
    if(strlen(unitclasses[i].units->unit)!=strlen(unitlabel))continue;
    if(STRCMP(unitclasses[i].units->unit,unitlabel)==0){
      *unitclass=i;
      *unittype=unitclasses[i].unit_index;
      return;
    }
  }
}

/* ------------------ InitUnits ------------------------ */

void InitUnits(void){
  int i;
  f_unit *units;
  f_units *ut;

  nunitclasses_default=6;
  NewMemory((void **)&unitclasses_default,nunitclasses_default*sizeof(f_units));

  for(i=0;i<nunitclasses_default;i++){
    unitclasses_default[i].submenuid=0;
  }
  ut=unitclasses_default;

  // temperature units

  ut->diff_index=-1;
  ut->nunits=3;
  ut->unit_index=0;
  strcpy(ut->unitclass,"Temperature");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,degC);
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,degF);
  units[1].scale[0]=1.8;
  units[1].scale[1]=32.0;
  strcpy(units[2].unit,"K");
  units[2].scale[0]=1.0;
  units[2].scale[1]=273.15;

  // velocity units

  ut = unitclasses_default+1;
  ut->unit_index=0;
  ut->diff_index=-1;
  ut->nunits=3;
  strcpy(ut->unitclass,"Velocity");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"m/s");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"mph");
  units[1].scale[0]=2.236931818;
  units[1].scale[1]=0.0;
  strcpy(units[2].unit,"ft/s");
  units[2].scale[0]=3.2808333333;
  units[2].scale[1]=0.0;
  CheckMemory;

  // distance units

  ut = unitclasses_default + 2;
  ut->unit_index=0;
  ut->diff_index=-1;
  ut->nunits=2;
  strcpy(ut->unitclass,"Distance");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"m");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"ft");
  units[1].scale[0]=3.280833333;
  units[1].scale[1]=0.0;
  CheckMemory;

  // volume flow units

  ut = unitclasses_default + 3;
  ut->unit_index=0;
  ut->diff_index=-1;
  ut->nunits=2;
  strcpy(ut->unitclass,"Volume Flow");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"m^3/s");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"cfm");
  units[1].scale[0]=2118.86720;
  units[1].scale[1]=0.0;
  CheckMemory;

  // Mass fraction units

  ut = unitclasses_default + 4;
  ut->unit_index=0;
  ut->diff_index=-1;
  ut->nunits=2;
  strcpy(ut->unitclass,"Mass Fraction");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"kg/kg");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"g/kg");
  units[1].scale[0]=1000.0;
  units[1].scale[1]=0.0;
  CheckMemory;


  // Concentration

  ut = unitclasses_default + 5;
  ut->unit_index=0;
  ut->diff_index=-1;
  ut->nunits=2;
  strcpy(ut->unitclass,"Concentration");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"mol/mol");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"ppm");
  units[1].scale[0]=1000000.0;
  units[1].scale[1]=0.0;
  CheckMemory;
}
