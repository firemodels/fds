// $Date: 2008-03-05 16:31:28 -0500 (Wed, 05 Mar 2008) $ 
// $Revision: 1409 $
// $Author: gforney $

#include "options.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "transpose.h"
#include "MALLOC.h"
#include "svn_revision.h"

// svn revision character string
char csv_convert_revision[]="$Revision: 1409 $";

/* ------------------ handle_key ------------------------ */

int handle_key(char *key){
  int i;

  for(i=0;i<nkeys;i++){
    if(strcmp(keys[i],key)==0){
      return 1;
    }
  }
  return 0;
}

/* ------------------ get_key_val ------------------------ */

float get_key_val(char *key){
  int i;
  float return_val;

  for(i=0;i<nkeys;i++){
    if(strcmp(keys[i],key)==0){
      return_val = key_vals[i];
      return return_val;
    }
  }
  return 0.0;
}

/* ------------------ read_cfg ------------------------ */

int read_cfg(void){
  FILE *stream;
  int ikey=0;
  int i;
  char *token;

  stream=fopen(config_file,"r");
  if(stream==NULL)return -1;


  if(fgets(cfg_line0,1024,stream)==NULL){
    strcpy(cfg_line0,"");
    strcpy(cfg_line1,"");
    strcpy(cfg_line2,"");
    return -1;
  }
  if(fgets(cfg_line1,1024,stream)==NULL){
    strcpy(cfg_line0,"");
    strcpy(cfg_line1,"");
    strcpy(cfg_line2,"");
    return -1;
  }
  if(fgets(cfg_line2,1024,stream)==NULL){
    strcpy(cfg_line0,"");
    strcpy(cfg_line1,"");
    strcpy(cfg_line2,"");
    return -1;
  }
  fclose(stream);

  token=strtok(cfg_line0,",");
  if(token!=NULL){
    key_label=token;
    token=strtok(NULL,",");
    key_unit=token;
  }
  if(key_label!=NULL){
    key_label = trim_both(key_label);
  }
  if(key_unit!=NULL){
    key_unit = trim_both(key_unit);
  }
  nkeys=1;
  for(i=0;i<strlen(cfg_line1);i++){
    if(cfg_line1[i]==',')nkeys++;
  }
  NewMemory((void **)&keys,nkeys*sizeof(char *));
  token=strtok(cfg_line1,",");
  ikey=0;
  while(token!=NULL){
    keys[ikey++]=token;
    token=strtok(NULL,",");
  }
  for(i=0;i<ikey;i++){
    keys[i]=trim_both(keys[i]);
  }
  NewMemory((void **)&key_vals,nkeys*sizeof(float));
  token=strtok(cfg_line2,",");
  ikey=0;
  while(token!=NULL){
    sscanf(token,"%f",key_vals+ikey++);
    token=strtok(NULL,",");
  }
  return 0;
}

/* ------------------ convert_csv ------------------------ */

int convert_csv(char *csv_in, char *csv_out){
  FILE *stream;
  int i;
  int nlines;
  int ncommas, maxcommas;
  char line[1024];
  int icol, irow;
  float *times;
  int ntimes,ntimes2;
  float dt_max;

#define IJ(irow,icol) (icol) + (maxcommas+1)*(irow)

  read_cfg();

  stream = fopen(csv_in,"r");
  if(stream==NULL)return -1;

  // pass 1 - count lines, and fields

  ncommas=0;
  maxcommas=0;
  nlines=0;
  while(!feof(stream)){
    if(fgets(line,1024,stream)==NULL)break;
    nlines++;
    trim(line);
    ncommas=0;
    for(i=0;i<strlen(line);i++){
      if(line[i]==',')ncommas++;
    }
    if(ncommas>maxcommas)maxcommas=ncommas;
  }

  // allocate memory
  csv_fields=NULL;
  NewMemory((void **)&csv_fields,nlines*(maxcommas+1)*sizeof(char *));
  for(i=0;i<nlines*(maxcommas+1);i++){
    csv_fields[i]=NULL;
  }

  rewind(stream);
  irow=0;
  while(!feof(stream)){
    int nlength;
    char *line2;
    int ncommas=0;
    int icol;
    int ij;

    if(fgets(line,1024,stream)==NULL)break;

    trim(line);
    nlength = strlen(line);

    NewMemory((void **)&line2,nlength+1);
    strcpy(line2,line);


    ncommas=0;
    icol=0;
    ij = IJ(irow,icol);
    csv_fields[ij]=line2;
    for(i=0;i<nlength;i++){
      if(line2[i]=='"')line2[i]=' ';
      if(line2[i]==','){
        line2[i]=0;
        icol++;
        ij = IJ(irow,icol);
        csv_fields[ij]=line2+i+1;
      }
    }
    line2[nlength]=0;
    irow++;
  }
  fclose(stream);
  for(i=0;i<nlines*(maxcommas+1);i++){
    csv_fields[i]=trim_both(csv_fields[i]);
  }


  ntimes = nlines-2;
  NewMemory((void **)&times,ntimes*sizeof(float));
  for(i=0;i<ntimes;i++){
    char *field;
    int ij;

    ij=IJ(i+2,0);
    field = csv_fields[ij];
    sscanf(field,"%f",times+i);
  }
  dt_max=times[1]-times[0];
  for(i=2;i<ntimes;i++){
    float dt;

    dt = times[i]-times[i-1];
    if(dt>dt_max)dt_max=dt;
  }
  if(dt_skip<dt_max)dt_skip=dt_max+1;
  ntimes2=0;
  for(i=0;i<ntimes;i++){
    if(i*dt_skip<=times[ntimes-1])ntimes2++;
  }


  stream=fopen(csv_out,"w");
  if(stream==NULL)return -1;

  if(key_unit!=NULL){
    fprintf(stream,"%s,",key_unit);
  }
  else{
    fprintf(stream,",");
  }
  for(i=0;i<ntimes2;i++){
    fprintf(stream,"%i",(int)i*dt_skip);
    if(i!=ntimes2-1)fprintf(stream,",");
  }
  fprintf(stream,"\n");
  for(icol=0;icol<=maxcommas;icol++){
    if(icol!=0&&handle_key(csv_fields[IJ(1,icol)])==0)continue;
    for(irow=0;irow<nlines;irow++){
      char *field;
      int ij;

      //  times[i] int(f2)*dt_skip times[i+1]
      //   v1         vx              v2
      //   (vx-v1)/(v2-v1) = (ts-t1)/(t2-t1)
      //  vx = v1 + (v2-v1)*(ts-t1)/(t2-t1)
      //      = (v1*t2 - v1*t1 - v2*t1 + v1*t1 + ts*(v2-v1))/(t2-t1)
      //      = (v1*t2 - v2*t1 + ts*(v2-v1))/(t2-t1)
      //      = (v1*(t2-ts) + v2*(ts - t1) )/(t2-t1)
      //      = w1*v1 + w2*v2
      //   where w1 = (t2-ts)/(t2-t1) and w2 = 1 - w1

      if(irow==2){
        float val;

        if(icol==0){
          if(key_label!=NULL){
            fprintf(stream,"%s",key_label);
          }
          else{
            fprintf(stream,"");
          }
        }
        ij=IJ(1,icol);
        field = csv_fields[ij];
        if(handle_key(field)==1){
          val = get_key_val(field);
          fprintf(stream,"%f,",val);
        }
        else{
          fprintf(stream,",",val);
        }
      }
      if(irow>1&&irow<nlines-1){
        float t1, ts, t2, w1, w2, v1,v2;

        t1 = times[irow-2];
        t2 = times[irow-1];
        if(irow==2){
          ts = 0.0;
        }
        else{
          if((int)(t1/dt_skip)==(int)(t2/dt_skip))continue;
          ts = (int)(t2/dt_skip)*dt_skip;
        }
        w1 = (t2-ts)/(t2-t1);
        w2 = 1.0 - w1;
        ij=IJ(irow,icol);
        field = csv_fields[ij];
        sscanf(field,"%f",&v1);
        ij=IJ(irow+1,icol);
        field = csv_fields[ij];
        sscanf(field,"%f",&v2);
        if(icol==0){
          fprintf(stream,"DT_%04i,",(int)ts);
        }
        else{
          fprintf(stream,"%f,",w1*v1+w2*v2);
        }
        continue;
      }
      if(irow==nlines-1){
        fprintf(stream,"\n");
      }
      if(irow!=nlines-1&&irow>=2){
        fprintf(stream,",");
      }
    }
  }

  fclose(stream);
  return 0;
}
