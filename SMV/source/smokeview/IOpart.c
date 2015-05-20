#define XYZ_EXTRA 7

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include GLUT_H
#include <math.h>

#include "smv_endian.h"
#include "update.h"
#include "smokeviewvars.h"

int tagscompare( const void *arg1, const void *arg2 );
void copy_dep_vals(part5class *partclassi, part5data *datacopy, float *colorptr, propdata *prop, int j);

void draw_SVOBJECT(sv_object *object, int frame_index_local, propdata *prop, int recurse_level,float *valrgb, int vis_override);
void update_all_partvis(partdata *parti);
void update_partvis(int first_frame,partdata *parti, part5data *datacopy, int nclasses);
int get_tagindex(const partdata *parti, part5data **data, int tagval);
#define READPASS 1
#define READFAIL 0

#define FORTPART5READ(var,size) \
returncode=READPASS;\
FSEEK(PART5FILE,4,SEEK_CUR);if(ferror(PART5FILE)==1||feof(PART5FILE)==1)returncode=READFAIL;\
if(returncode==READPASS){\
  fread(var,4,size,PART5FILE);\
  if(ferror(PART5FILE)==1||feof(PART5FILE)==1)returncode=READFAIL;\
}\
if(returncode==READPASS){\
  if(endianswitch==1)endian_switch(var,size);\
  FSEEK(PART5FILE,4,SEEK_CUR);\
  if(ferror(PART5FILE)==1||feof(PART5FILE)==1)returncode=READFAIL;\
}

/* ------------------ drawpart_frame ------------------------ */

void drawpart_frame(void){
  partdata *parti;

  if(staticframe0==0||frame_index!=0){
    int i;

    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->loaded==0||parti->display==0)continue;
      if(parti->evac==1){
        drawEvac(parti);
        SNIFF_ERRORS("after drawEvac");
      }
      else{
        drawPart(parti);
        SNIFF_ERRORS("after drawPart");
      }
    }
  }
  if(visStaticSmoke==1&&staticframe0==1){
    int i;

    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->loaded==0||parti->display==0)continue;
      drawStaticPart(parti);
      SNIFF_ERRORS("after drawStaticPart");
    }
  }
}

/* ------------------ drawevac_frame ------------------------ */

void drawevac_frame(void){
  int i;

  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->loaded==0||parti->display==0||parti->evac==0)continue;
    drawEvac(parti);
  }
  SNIFF_ERRORS("after drawEvac 2");
}

/* ------------------ freepart5data ------------------------ */

void freepart5data(part5data *datacopy){
  FREEMEMORY(datacopy->cvals);
  FREEMEMORY(datacopy->sx);
  FREEMEMORY(datacopy->sy);
  FREEMEMORY(datacopy->sz);
  FREEMEMORY(datacopy->dsx);
  FREEMEMORY(datacopy->dsy);
  FREEMEMORY(datacopy->dsz);
  FREEMEMORY(datacopy->avatar_angle);
  FREEMEMORY(datacopy->avatar_width);
  FREEMEMORY(datacopy->avatar_height);
  FREEMEMORY(datacopy->avatar_depth);
  FREEMEMORY(datacopy->tags);
  FREEMEMORY(datacopy->sort_tags);
  FREEMEMORY(datacopy->vis_part);
  FREEMEMORY(datacopy->rvals);
  FREEMEMORY(datacopy->irvals);
}

/* ------------------ freeallpart5data ------------------------ */

void freeallpart5data(partdata *parti){
  int i;
  part5data *datacopy;

  if(parti->data5==NULL)return;
  datacopy = parti->data5;
  for(i=0;i<parti->ntimes*parti->nclasses;i++){
    freepart5data(datacopy);
    datacopy++;
  }
  FREEMEMORY(parti->data5);
}

/* ------------------ initpart5data ------------------------ */

void initpart5data(part5data *datacopy, part5class *partclassi){
  datacopy->cvals=NULL;
  datacopy->partclassbase=partclassi;
  datacopy->sx=NULL;
  datacopy->sy=NULL;
  datacopy->sz=NULL;
  datacopy->dsx=NULL;
  datacopy->dsy=NULL;
  datacopy->dsz=NULL;
  datacopy->avatar_angle=NULL;
  datacopy->avatar_width=NULL;
  datacopy->avatar_height=NULL;
  datacopy->avatar_depth=NULL;
  datacopy->tags=NULL;
  datacopy->vis_part=NULL;
  datacopy->sort_tags=NULL;
  datacopy->rvals=NULL;
  datacopy->irvals=NULL;
}


/* ------------------ getpart5data ------------------------ */

void getpart5data(partdata *parti, int partframestep_local, int partpointstep_local, int nf_all, float *delta_time, FILE_SIZE *file_size){
  FILE *PART5FILE;
  int one;
  int endianswitch=0;
  int version;
  int nclasses;
  int i;
  int skip_local;
  size_t returncode;
  float time_local;
  int nparts;
  int *numtypes=NULL,*numtypescopy, *numpoints=NULL;
  int numtypes_temp[2];
  char *reg_file;
  part5data *datacopy;
  int count;
  int count2;
  int first_frame=1;
  int local_starttime=0, local_stoptime=0;

  reg_file=parti->reg_file;

  PART5FILE=fopen(reg_file,"rb");
  if(PART5FILE==NULL)return;

  *file_size=get_filesize(reg_file);
  FSEEK(PART5FILE,4,SEEK_CUR);fread(&one,4,1,PART5FILE);FSEEK(PART5FILE,4,SEEK_CUR);
  if(one!=1)endianswitch=1;

  FORTPART5READ(&version,1);if(returncode==0)goto wrapup;
  FORTPART5READ(&nclasses,1);if(returncode==0)goto wrapup;
  NewMemory((void **)&numtypes,2*nclasses*sizeof(int));
  NewMemory((void **)&numpoints,nclasses*sizeof(int));
  numtypescopy=numtypes;
  numtypes_temp[0]=0;
  numtypes_temp[1]=0;
  CheckMemory;
  for(i=0;i<nclasses;i++){
    FORTPART5READ(numtypes_temp,2);
    if(returncode==0)goto wrapup;
    *numtypescopy++=numtypes_temp[0];
    *numtypescopy++=numtypes_temp[1];
    skip_local = 2*(numtypes_temp[0]+numtypes_temp[1])*(8 + 30);
    returncode=FSEEK(PART5FILE,skip_local,SEEK_CUR);
    if(returncode!=0)goto wrapup;
  }
  CheckMemory;

  datacopy = parti->data5;
  count=0;
  count2=-1;
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  for(;;){
    int doit;

    CheckMemory;
    if(count>=nf_all)break;
    FORTPART5READ(&time_local,1);
    
    if(count%partframestep_local!=0||(settmin_p==1&&time_local<tmin_p-TEPS)||(settmax_p==1&&time_local>tmax_p+TEPS)){
      doit=0;
    }
    else{
      count2++;
      doit=1;
    }
    count++;

    if(returncode==0)break;
    if(doit==1){
      PRINTF("particle time=%.2f",time_local);
      parti->times[count2]=time_local;
    }
    for(i=0;i<nclasses;i++){
      part5class *partclassi;
      int factor=256*128-1;

      partclassi = parti->partclassptr[i];
      FORTPART5READ(&nparts,1);
      if(returncode==0)goto wrapup;
      numpoints[i]=nparts;
      skip_local=0;
      CheckMemory;
      if(doit==1){
        short *sx, *sy, *sz;
        float *xyz;
        float *angle, *width, *depth, *height;
        int j;

        if(parti->evac==1){
          FORTPART5READ(partclassi->xyz,XYZ_EXTRA*nparts);
        }
        else{
          FORTPART5READ(partclassi->xyz,3*nparts);
        }
        CheckMemory;
        if(nparts>0){
          if(returncode==0)goto wrapup;
          xyz = partclassi->xyz;
          sx = datacopy->sx;
          sy = datacopy->sy;
          sz = datacopy->sz;
          if(parti->evac==1){
            angle=datacopy->avatar_angle;
            width=datacopy->avatar_width;
            depth=datacopy->avatar_depth;
            height=datacopy->avatar_height;
          }
          for(j=0;j<nparts;j++){
            float xx, yy, zz;

            xx = NORMALIZE_X(xyz[         j])/xbar;
            yy = NORMALIZE_Y(xyz[  nparts+j])/ybar;
            zz = NORMALIZE_Z(xyz[2*nparts+j])/zbar;

            sx[j] = factor*xx;
            sy[j] = factor*yy;
            sz[j] = factor*zz;
            if(parti->evac==1){
              angle[j] =xyz[j+3*nparts];
              width[j] =xyz[j+4*nparts];
              depth[j] =xyz[j+5*nparts];
              height[j]=xyz[j+6*nparts];
            }
          }
          CheckMemory;
        }
      }
      else{
        if(parti->evac==1){
          skip_local += 4 + XYZ_EXTRA*4*nparts + 4;  
        }
        else{
          skip_local += 4 + 3*4*nparts + 4;  
        }
      }
      CheckMemory;
      if(doit==1){
        int *sort_tags;
        int j;

        sort_tags=datacopy->sort_tags;
        FORTPART5READ(datacopy->tags,nparts);
        CheckMemory;
        if(nparts>0){
          if(returncode==0)goto wrapup;
          for(j=0;j<nparts;j++){
            sort_tags[2*j]=datacopy->tags[j];
            sort_tags[2*j+1]=j;
          }
          qsort( sort_tags, (size_t)nparts, 2*sizeof(int), tagscompare );
        }
      }
      else{
        skip_local += 4 + 4*nparts + 4;  // skip over tag for now
      }
      CheckMemory;
      if(doit==1){
        if(numtypes[2*i]>0){
          FORTPART5READ(datacopy->rvals,nparts*numtypes[2*i]);
          if(returncode==0)goto wrapup;
        }
      }
      else{
        if(numtypes[2*i]>0){
          skip_local += 4 + 4*nparts*numtypes[2*i] + 4;  // skip over vals for now
        }
      }
      CheckMemory;
      if(numtypes[2*i+1]>0){
        skip_local += 4 + 4*nparts*numtypes[2*i+1] + 4;
      }

      
      returncode=0;
      if(skip_local>0){
        returncode=FSEEK(PART5FILE,skip_local,SEEK_CUR);
        if(returncode!=0)goto wrapup;
      }
      CheckMemory;
      if(doit==1)datacopy++;
    }
    CheckMemory;
    if(first_frame==1)first_frame=0;
    if(doit==1)PRINTF(" completed\n");

  }
wrapup:
  local_stoptime = glutGet(GLUT_ELAPSED_TIME);
  *delta_time=(local_stoptime-local_starttime)/1000.0;

  CheckMemory;
  update_all_partvis(parti);
  FREEMEMORY(numtypes);
  FREEMEMORY(numpoints);
  fclose(PART5FILE);
}

/* ------------------ get_part5prop ------------------------ */

void print_part5prop(void){
  int i;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;
    PRINTF("label=%s min=%f max=%f\n",propi->label->longlabel,propi->valmin,propi->valmax);
    PRINTF("   glbmin=%f glbmax=%f\n",propi->global_min,propi->global_max);
    PRINTF("   permin=%f permax=%f\n",propi->percentile_min,propi->percentile_max);
    PRINTF("\n");
  }
}


/* ------------------ get_part5prop_index_s ------------------------ */

int get_part5prop_index_s(char *shortlabel){
  int i;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;
    if(strcmp(propi->label->shortlabel,shortlabel)==0)return i;
  }
  return -1;
}

/* ------------------ get_part5prop_index ------------------------ */

int get_part5prop_index(char *label){
  int i;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;
    if(strcmp(propi->label->longlabel,label)==0)return i;
  }
  return 0;
}

/* ------------------ get_part5prop_s ------------------------ */

part5prop *get_part5prop_s(char *label){
  int i;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;
    if(strcmp(propi->label->shortlabel,label)==0)return propi;
  }
  return NULL;
}

/* ------------------ get_part5prop ------------------------ */

part5prop *get_part5prop(char *label){
  int i;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;
    if(strcmp(propi->label->longlabel,label)==0)return propi;
  }
  return NULL;
}

/* ------------------ init_part5prop ------------------------ */

void init_part5prop(void){
  int i,j,k;

  // 0.  only needed if init_part5prop is called more than once
  // (and if so, need to also free memory of each component)

  FREEMEMORY(part5propinfo);
  npart5prop=0;

  // 1.  count max number of distinct variables

  for(i=0;i<npartclassinfo;i++){
    part5class *partclassi;

    partclassi = partclassinfo + i;
    npart5prop+=(partclassi->ntypes-1);
  }

  // 2. now count the exact amount and put labels into array just allocated

  if(npart5prop>0){
    NewMemory((void **)&part5propinfo,npart5prop*sizeof(part5prop));
    npart5prop=0;

    for(i=0;i<npartclassinfo;i++){
      int ii;
      part5class *partclassi;

      partclassi = partclassinfo + i;
      for(j=1;j<partclassi->ntypes;j++){
        flowlabels *flowlabel;
        int define_it;

        define_it = 1;
        flowlabel = partclassi->labels + j;
        for(k=0;k<npart5prop;k++){
          part5prop *propi;
          char *proplabel;

          propi = part5propinfo + k;
          proplabel = propi->label->longlabel;
          if(strcmp(proplabel,flowlabel->longlabel)==0){
            define_it=0;
            break;
          }
        }
        if(define_it==1){
          part5prop *propi;

          propi = part5propinfo + npart5prop;

          propi->human_property=0;
          propi->particle_property=0;
          propi->label=flowlabel;

          propi->setvalmin=0;
          propi->setvalmax=0;
          propi->set_global_bounds=1;
          propi->global_min=100000000.0;
          propi->global_max=-propi->global_min;
          propi->valmin=1.0;
          propi->valmax=0.0;
          propi->percentile_min=1.0;
          propi->percentile_max=0.0;
          propi->user_min=1.0;
          propi->user_max=0.0;
          propi->display=0;


          propi->setchopmin=0;
          propi->setchopmax=0;
          propi->chopmin=1.0;
          propi->chopmax=0.0;

          propi->buckets=NULL;
          propi->partlabels=NULL;
          NewMemory((void **)&propi->partlabels,256*sizeof(char *));
          for(ii=0;ii<256;ii++){
            char *labeli;

            labeli=NULL;
            NewMemory((void **)&labeli,11);
            propi->partlabels[ii]=labeli;
          }
          NewMemory((void **)&propi->scale,256);
          

          npart5prop++;
        }
      }

    }
  }
  for(i=0;i<npart5prop;i++){
    part5prop *propi;
    int ii;

    propi = part5propinfo + i;

    propi->class_present=NULL;
    propi->class_vis=NULL;
    propi->class_types=NULL;
    NewMemory((void **)&propi->class_types,npartclassinfo*sizeof(unsigned int));
    NewMemory((void **)&propi->class_present,npartclassinfo*sizeof(unsigned char));
    NewMemory((void **)&propi->class_vis,npartclassinfo*sizeof(unsigned char));
    for(ii=0;ii<npartclassinfo;ii++){
      propi->class_vis[ii]=1;
      propi->class_present[ii]=0;
      propi->class_types[ii]=0;
    }
  }
  for(i=0;i<npartclassinfo;i++){
    part5class *partclassi;
 
    partclassi = partclassinfo + i;
    for(j=1;j<partclassi->ntypes;j++){
      flowlabels *flowlabel;
      part5prop *classprop;

      flowlabel = partclassi->labels + j;
      classprop = get_part5prop(flowlabel->longlabel);
      if(classprop!=NULL){
        if(partclassi->kind==1){
          classprop->human_property=1;
        }
        else{
          classprop->particle_property=1;
        }
        classprop->class_present[i]=1;
        classprop->class_types[i]=j-2;
      }
    }
  }
}

/* ------------------ update_partvis ------------------------ */

void update_all_partvis(partdata *parti){
  part5data *datacopy;
  int i,j;
  int firstframe=1;

  datacopy = parti->data5;
  for(i=0;i<parti->ntimes;i++){
    for(j=0;j<parti->nclasses;j++){
      update_partvis(firstframe,parti,datacopy,parti->nclasses);
      datacopy++;
    }
    if(firstframe==1)firstframe=0;
  }
}


/* ------------------ update_partvis ------------------------ */

void update_partvis(int first_frame,partdata *parti, part5data *datacopy, int nclasses){
  int nparts;
  unsigned char *vis_part;

  nparts=datacopy->npoints;
  vis_part=datacopy->vis_part;

  if(first_frame==1){
    int ii;
    for(ii=0;ii<nparts;ii++){
      if(ii%partpointstep==0){
        vis_part[ii]=1;
      }
      else{
        vis_part[ii]=0;
      }
    }
  }
  else{
    int ii;
    part5data *datalast;
    int nvis=0,nleft;
      
    for(ii=0;ii<nparts;ii++){
      int tag_index;

      datalast = datacopy-nclasses;
      tag_index = get_tagindex(parti,&datalast,datacopy->tags[ii]);
      if(partpointstep==1||(tag_index!=-1&&datalast->vis_part[tag_index]==1)){
        datacopy->vis_part[ii]=1;
        nvis++;
      }
      else{
        datacopy->vis_part[ii]=0;
      }
    }

    nleft = nparts/partpointstep - nvis;
    if(nleft>0){
      for(ii=0;ii<nparts;ii++){
        if(datacopy->vis_part[ii]==1)continue;
        if(nleft>0){
          datacopy->vis_part[ii]=1;
          nleft--;
        }
      }
    }
  }
}

/* ------------------ get_tagindex ------------------------ */

int get_tagindex(const partdata *partin, part5data **datain, int tagval){
  int *returnval;
  part5data *data;
  int i;

  for(i=-1;i<npartinfo;i++){
    const partdata *parti;

    if(i==-1){
      parti=partin;
      data = *datain;
    }
    else{
      parti=partinfo + i;
      if(parti==partin)continue;
      if(parti->loaded==0||parti->display==0)continue;
      data = parti->data5+(*datain-partin->data5);
    }
    if(parti->loaded==0||parti->display==0)continue;

    if(data->npoints==0)continue;
    ASSERT(data->sort_tags!=NULL);
    returnval=bsearch(&tagval,data->sort_tags,data->npoints,2*sizeof(int),tagscompare);
    if(returnval==NULL)continue;
    *datain=data;
    return *(returnval+1);
  }
  return -1;
}

/* ------------------ getpart5nframes ------------------------ */

int getpart5nframes(partdata *parti){
  FILE *stream;
  char buffer[256];
  float time_local;
  //int count;
  char *reg_file, *size_file;
  int i;
  int stat_sizefile, stat_regfile;
  STRUCTSTAT stat_sizefile_buffer, stat_regfile_buffer;
  int nframes_all;

  reg_file=parti->reg_file;
  size_file=parti->size_file;

  // if size file doesn't exist then generate it

  stat_sizefile=STAT(size_file,&stat_sizefile_buffer);
  stat_regfile=STAT(reg_file,&stat_regfile_buffer);
  if(stat_regfile!=0)return -1;

  // create a size file if 1) the size does not exist
  //                       2) base file is newer than the size file
  if(stat_sizefile!=0||
    stat_regfile_buffer.st_mtime>stat_sizefile_buffer.st_mtime){
    //create_part5sizefile(reg_file,size_file);
      {
        int lenreg, lensize, error;
        int angle_flag=0;

        trim(reg_file);
        trim(size_file);
        lenreg=strlen(reg_file);
        lensize=strlen(size_file);
        if(parti->evac==1){
          angle_flag=1;
          FORTfcreate_part5sizefile(reg_file,size_file, &angle_flag, &redirect, &error, lenreg,lensize);
        }
        else{
          angle_flag=0;
          FORTfcreate_part5sizefile(reg_file,size_file, &angle_flag, &redirect, &error, lenreg,lensize);
        }
      }
  }
  
  stream=fopen(size_file,"r");
  if(stream==NULL)return -1;

  nframes_all=0;
  for(;;){
    int exitloop;

    if(fgets(buffer,255,stream)==NULL)break;
    sscanf(buffer,"%f",&time_local);
    exitloop=0;
    for(i=0;i<parti->nclasses;i++){
      if(fgets(buffer,255,stream)==NULL){
        exitloop=1;
        break;
      }
    }
    if(exitloop==1)break;
    nframes_all++;
  }
  fclose(stream);
  return nframes_all;
}

/* ------------------ get_partframes ------------------------ */

int get_min_partframes(void){
  int i;
  int min_frames=-1;

  for(i=0;i<npartinfo;i++){
    partdata *parti;
    int nframes;

    parti = partinfo + i;
    nframes = getpart5nframes(parti);
    if(nframes>0){
      if(min_frames==-1){
        min_frames=nframes;
      }
      else{
        if(nframes!=-1&&nframes<min_frames)min_frames=nframes;
      }
    }
  }
  return min_frames;
}

/* ------------------ getpart5header ------------------------ */

void getpart5header(partdata *parti, int partframestep_local, int *nf_all){
  FILE *stream;
  char buffer[256];
  float time_local;
  int count;
  char *reg_file, *size_file;
  int i;
  int stat_sizefile, stat_regfile;
  STRUCTSTAT stat_sizefile_buffer, stat_regfile_buffer;
  int nframes_all;

  reg_file=parti->reg_file;
  size_file=parti->size_file;

  // if size file doesn't exist then generate it

  parti->ntimes=0;

  stat_sizefile=STAT(size_file,&stat_sizefile_buffer);
  stat_regfile=STAT(reg_file,&stat_regfile_buffer);
  if(stat_regfile!=0)return;

  // create a size file if 1) the size does not exist
  //                       2) base file is newer than the size file
  if(stat_sizefile!=0||
    stat_regfile_buffer.st_mtime>stat_sizefile_buffer.st_mtime){
    //create_part5sizefile(reg_file,size_file);
      {
        int lenreg, lensize, error;
        int angle_flag=0;

        trim(reg_file);
        trim(size_file);
        lenreg=strlen(reg_file);
        lensize=strlen(size_file);
        if(parti->evac==1){
          angle_flag=1;
          FORTfcreate_part5sizefile(reg_file,size_file, &angle_flag, &redirect, &error, lenreg,lensize);
        }
        else{
          angle_flag=0;
          FORTfcreate_part5sizefile(reg_file,size_file, &angle_flag, &redirect, &error, lenreg,lensize);
        }
      }
  }
  
  stream=fopen(size_file,"r");
  if(stream==NULL)return;

    // pass 1: count frames
  
  nframes_all=0;
  for(;;){
    int exitloop;

    if(fgets(buffer,255,stream)==NULL)break;
    sscanf(buffer,"%f",&time_local);
    exitloop=0;
    for(i=0;i<parti->nclasses;i++){
      if(fgets(buffer,255,stream)==NULL||(npartinfo>1&&npartframes_max!=-1&&nframes_all+1>npartframes_max)){
        exitloop=1;
        break;
      }
    }
    if(exitloop==1)break;
    nframes_all++;
    if((nframes_all-1)%partframestep_local!=0||
       (settmin_p!=0&&time_local<tmin_p-TEPS)||
       (settmax_p!=0&&time_local>tmax_p+TEPS)){
       continue;
    }
    (parti->ntimes)++;
  }
  rewind(stream);
  *nf_all = nframes_all;

  // allocate memory for number of time steps * number of classes

  NewMemory((void **)&parti->data5,parti->nclasses*parti->ntimes*sizeof(part5data));
  NewMemory((void **)&parti->times,parti->ntimes*sizeof(float));

  // free memory for x, y, z frame data 

  for(i=0;i<parti->nclasses;i++){
    part5class *partclassi;

    partclassi = parti->partclassptr[i];
    FREEMEMORY(partclassi->xyz);
    partclassi->maxpoints=0;
  }

  // pass 2 - allocate memory for x, y, z frame data
  //          
  {
    part5data *datacopy;
    int fail;

    fail=0;
    count=-1;
    datacopy=parti->data5;
    for(i=0;i<nframes_all;i++){
      int j;

      count++;
      fail=0;
      if(fgets(buffer,255,stream)==NULL){
        fail=1;
        break;
      }
      sscanf(buffer,"%f",&time_local);
      if(count%partframestep_local!=0||
         (settmin_p!=0&&time_local<tmin_p-TEPS)||
         (settmax_p!=0&&time_local>tmax_p+TEPS)){
        for(j=0;j<parti->nclasses;j++){
          if(fgets(buffer,255,stream)==NULL){
            fail=1;
            break;
          }
        }
        if(fail==1)break;
        continue;
      }
      for(j=0;j<parti->nclasses;j++){
        int npoints ,ntypes;

        part5class *partclassj;

        datacopy->time = time_local;
        partclassj = parti->partclassptr[j];
        initpart5data(datacopy,partclassj);
        if(fgets(buffer,255,stream)==NULL){
          fail=1;
          break;
        }
        sscanf(buffer,"%i",&datacopy->npoints);
        npoints=datacopy->npoints;
        if(npoints>partclassj->maxpoints)partclassj->maxpoints=npoints;
        if(npoints>0){
          NewMemory((void **)&datacopy->tags,npoints*sizeof(int));
          NewMemory((void **)&datacopy->sort_tags,2*npoints*sizeof(int));
          NewMemory((void **)&datacopy->vis_part,npoints*sizeof(unsigned char));
          NewMemory((void **)&datacopy->sx,npoints*sizeof(short));
          NewMemory((void **)&datacopy->sy,npoints*sizeof(short));
          NewMemory((void **)&datacopy->sz,npoints*sizeof(short));
          NewMemory((void **)&datacopy->dsx,npoints*sizeof(float));
          NewMemory((void **)&datacopy->dsy,npoints*sizeof(float));
          NewMemory((void **)&datacopy->dsz,npoints*sizeof(float));
          if(parti->evac==1){
            NewMemory((void **)&datacopy->avatar_angle,npoints*sizeof(float));
            NewMemory((void **)&datacopy->avatar_width,npoints*sizeof(float));
            NewMemory((void **)&datacopy->avatar_depth,npoints*sizeof(float));
            NewMemory((void **)&datacopy->avatar_height,npoints*sizeof(float));
          }
          ntypes = datacopy->partclassbase->ntypes;
          if(ntypes>0){
            NewMemory((void **)&datacopy->rvals, ntypes*npoints*sizeof(float));
            NewMemory((void **)&datacopy->irvals,ntypes*npoints*sizeof(unsigned char));
          }
        }
        datacopy++;
      }
      if(fail==1)break;
    }
    if(fail==1)parti->ntimes=i;
    fclose(stream);
  }

  // allocate memory for x, y, z and tag for the maximum frame size
  //           don't need to allocate memory for all frames

  for(i=0;i<parti->nclasses;i++){
    part5class *partclassi;

    partclassi = parti->partclassptr[i];
    if(partclassi->maxpoints>0){
      if(parti->evac==1){
        NewMemory((void **)&partclassi->xyz,XYZ_EXTRA*partclassi->maxpoints*sizeof(float));
      }
      else{
        NewMemory((void **)&partclassi->xyz,3*partclassi->maxpoints*sizeof(float));
      }
    }
  }

}

/* ------------------ readpart5 ------------------------ */

void readpart5(char *file, int ifile, int flag, int *errorcode){
  size_t lenfile;
  int error=0;
  partdata *parti;
  int nf_all;
  int local_starttime0, local_stoptime0;
  float delta_time0, delta_time;
  FILE_SIZE file_size;

  local_starttime0 = glutGet(GLUT_ELAPSED_TIME);

  ASSERT(ifile>=0&&ifile<npartinfo);
  parti=partinfo+ifile;

  freeallpart5data(parti);

  if(parti->loaded==0&&flag==UNLOAD)return;


  *errorcode=0;
  partfilenum=ifile;
  if(parti->evac==0){
    ReadPartFile=0;
  }
  else{
    ReadEvacFile=0;
  }
  parti->loaded=0;
  parti->display=0;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatemenu=1;

  FREEMEMORY(parti->times); 

  if(colorlabelpart!=NULL){
    int n;

    for(n=0;n<MAXRGB;n++){
      FREEMEMORY(colorlabelpart[n]);
    }
    FREEMEMORY(colorlabelpart);
  }

  if(flag==UNLOAD){
    Update_Times();
    updatemenu=1;
    updatePart5extremes();
#ifdef pp_MEMPRINT
    PRINTF("After particle file unload: \n");
    PrintMemoryInfo;
#endif
    return;
  }

  lenfile = strlen(file);
  if(lenfile==0){
    readpart("",ifile,UNLOAD,&error);
    Update_Times();
    return;
  }
  
  PRINTF("Sizing particle data: %s\n",file);
  getpart5header(parti, partframestep, &nf_all);

  PRINTF("Loading particle data: %s\n",file);
  getpart5data(parti,partframestep,partpointstep, nf_all, &delta_time, &file_size);
  updateglui();

#ifdef pp_MEMPRINT
  PRINTF("After particle file load: \n");
  PrintMemoryInfo;
#endif
  if(parti->evac==0){
    ReadPartFile=1;
  }
  else{
    ReadEvacFile=1;
  }
  if(parti->evac==0){
    visParticles=1;
  }
  else{
    visEvac=1;
  }
  /* convert particle temperatures into integers pointing to an rgb color table */

  PRINTF("computing particle color levels \n");

  adjustpart5bounds(parti);
  NewMemory((void **)&colorlabelpart,MAXRGB*sizeof(char *));
  {
    int n;

    for(n=0;n<MAXRGB;n++){
      colorlabelpart[n]=NULL;
    }
    for(n=0;n<nrgb;n++){
      NewMemory((void **)&colorlabelpart[n],11);
    }
  }
  getPart5Colors(parti,nrgb);
  updateglui();
#ifdef pp_MEMPRINT
  PRINTF("After particle file load: \n");
  PrintMemoryInfo;
#endif
  if(parti->evac==0){
    visParticles=1;
    ReadPartFile=1;
  }
  else{
    visEvac=1;
    ReadEvacFile=1;
  }

  parttype=0;
  PART_CB_INIT();
  ParticlePropShowMenu(part5colorindex);
  parti->loaded=1;
  parti->display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  Update_Times();
  updatePart5extremes();
  updatemenu=1;
  Idle_CB();

  local_stoptime0 = glutGet(GLUT_ELAPSED_TIME);
  delta_time0 = (local_stoptime0-local_starttime0)/1000.0;

  if(file_size!=0&&delta_time>0.0){
    float loadrate;

    loadrate = ((float)file_size*8.0/1000000.0)/delta_time;
    PRINTF(" %.1f MB loaded in %.2f s - rate: %.1f Mb/s",
    (float)file_size/1000000.,delta_time,loadrate);
  }
  else{
    PRINTF(" %.1f MB downloaded in %.2f s",
    (float)file_size/1000000.,delta_time);
  }
  PRINTF(" (overhead: %.2f s)\n",delta_time0-delta_time);

  glutPostRedisplay();
}

/* ------------------ update_all_partvis2 ------------------------ */

void update_all_partvis2(void){
  partdata *parti;
  int i;
  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    if(parti->loaded==1)update_all_partvis(parti);
  }
}

/* ------------------ readpart ------------------------ */

void readpart(char *file, int ifile, int flag, int *errorcode){
  int nmax, n, i;
  FILE_SIZE lenfile;
  float *tcopy;
  unsigned char *isprinkcopy;
  int error=0;
  int bytesperpoint;
  int skip_local;
  int statfile,statfile2;
  STRUCTSTAT statbuffer,statbuffer2;
  char partsizefile[1024],buffer[1024];
  FILE *sizefile;
  int readpartsize=1;
  int partpointstepold, partframestepold;
  size_t return_code;
  int nb,nv;
  float xbox, ybox, zbox;
  partdata *parti;
  int blocknumber;
  mesh *meshi;
  float offset_x, offset_y, offset_z;
  int file_unit;

  ASSERT(ifile>=0&&ifile<=npartinfo);
  parti=partinfo+ifile;
  if(parti->version==1){
    readpart5(file,ifile,flag,errorcode);
    return;
  }
  blocknumber=parti->blocknumber;
  meshi=meshinfo+blocknumber;
  if(parti->loaded==0&&flag==UNLOAD)return;


  nb=meshi->nbptrs;
  nv=meshi->nvents;

  *errorcode=0;
  partfilenum=ifile;
  if(partinfo[ifile].evac==0){
    ReadPartFile=0;
  }
  else{
    ReadEvacFile=0;
  }
  partinfo[ifile].loaded=0;
  partinfo[ifile].display=0;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatemenu=1;

  FREEMEMORY(parti->times); 
  FREEMEMORY(parti->xpart);  FREEMEMORY(parti->ypart);  FREEMEMORY(parti->zpart);  
  FREEMEMORY(parti->xpartb); FREEMEMORY(parti->ypartb); FREEMEMORY(parti->zpartb); 
  FREEMEMORY(parti->xparts); FREEMEMORY(parti->yparts); FREEMEMORY(parti->zparts); 
  FREEMEMORY(parti->tpart);  FREEMEMORY(parti->itpart); 
  FREEMEMORY(parti->isprink);
  FREEMEMORY(parti->sframe); 
  FREEMEMORY(parti->bframe);
  FREEMEMORY(parti->sprframe)

  if(flag==UNLOAD){
    Update_Times();
    updatemenu=1;
#ifdef pp_MEMPRINT
    PRINTF("After particle file unload: \n");
    PrintMemoryInfo;
#endif
    return;
  }

  if(colorlabelpart!=NULL){
    PRINTF("freeing colorlabelpart\n");
    for(n=0;n<MAXRGB;n++){
      FREEMEMORY(colorlabelpart[n]);
    }
    FREEMEMORY(colorlabelpart);
  }

  lenfile = strlen(file);
  if(lenfile==0){
    readpart("",ifile,UNLOAD,&error);
    Update_Times();
    return;
  }
  
  PRINTF("Sizing particle data: %s\n",file);
  file_unit=15;
  FORTget_file_unit(&file_unit,&file_unit);
  FORTgetsizes(&file_unit,file,&nb,&nv,&nspr,&mxframepoints,&staticframe0,&error,lenfile);
  STRCPY(partsizefile,file);
  STRCAT(partsizefile,".sz");
  statfile=STAT(file,&statbuffer);
  statfile2=STAT(partsizefile,&statbuffer2);
  if(statfile==0&&statfile2==0&&difftime(statbuffer2.st_mtime,statbuffer.st_mtime)>0){
    sizefile=fopen(partsizefile,"r");
    if(sizefile!=NULL){
      if(fgets(buffer,255,sizefile)!=NULL){
        sscanf(buffer,"%i %i %i %i %i %i %i %i %i",
          &nb,&nv,&nspr,&mxframepoints,&staticframe0,&npartpoints,&npartframes,&partframestepold,&partpointstepold);
        fclose(sizefile);
        if(partframestepold==partframestep&&partpointstepold==partpointstep)readpartsize=0;
      }
    }
  }
  if(readpartsize==1){
    FORTgetdata1(&file_unit,&parttype,&error);
    FORTgetsizes2(&file_unit,&settmin_p,&tmin_p,&settmax_p,&tmax_p,
                     &nspr, &partframestep, &partpointstep, &npartpoints, &npartframes, &error);
    sizefile=fopen(partsizefile,"w");
    if(sizefile!=NULL){
      fprintf(sizefile,"%i %i %i %i %i %i %i %i %i",
          nb,nv,nspr,mxframepoints,staticframe0,npartpoints,npartframes,partframestep,partpointstep);
    }
    else{
      fprintf(stderr,"*** Error:  unable to write to %s\n",partsizefile);
    }
    file_unit=15;
    FORTget_file_unit(&file_unit,&file_unit);
    FORTgetsizes(&file_unit,file,&nb,&nv,&nspr,&mxframepoints,&staticframe0,&error,lenfile);
  }
  if(staticframe0==1)first_frame_index=1;
  if(error!=0){
    fprintf(stderr,"*** Error: problem reading %s\n",file);
    return;
  }
  if(npartpoints<=0){
    fprintf(stderr,"*** Warning: the particle file:%s is empty\n",file);
    return;
  }
  if(nspr>0){
    if(tspr==NULL){
      return_code=NewMemory((void **)&tspr,sizeof(float)*nspr);
    }
    else{
      return_code=ResizeMemory((void **)&tspr,sizeof(float)*nspr);
    }
    if(return_code==0){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpart("",ifile,UNLOAD,&error);
      return;
    }
  }

  if(NewMemory((void **)&parti->times,sizeof(float)*npartframes)==0){
    *errorcode=1;
    FORTclosefortranfile(&file_unit);
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  if(NewMemory((void **)&parti->xparts,npartpoints*sizeof(short))==0||
     NewMemory((void **)&parti->yparts,npartpoints*sizeof(short))==0||
     NewMemory((void **)&parti->zparts,npartpoints*sizeof(short))==0){
    *errorcode=1;
    FORTclosefortranfile(&file_unit);
    fprintf(stderr,"*** Error: memory allocation failed while attempting .\n");
    fprintf(stderr,"          to load %s\n",parti->file);
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  bytesperpoint=7;

  if(NewMemory((void **)&parti->tpart,npartpoints*sizeof(float))==0||
     NewMemory((void **)&parti->itpart,npartpoints*sizeof(unsigned char))==0||
     NewMemory((void **)&parti->isprink,npartpoints*sizeof(unsigned char))==0||
     NewMemory((void **)&parti->bframe,npartframes*sizeof(int))==0||
     NewMemory((void **)&parti->sframe,npartframes*sizeof(int))==0||
     NewMemory((void **)&parti->sprframe,npartframes*sizeof(int))==0){
      *errorcode=1;
    fprintf(stderr,"*** Error: memory allocation failed while attempting .\n");
    fprintf(stderr,"          to load %s\n",parti->file);
      FORTclosefortranfile(&file_unit);
      readpart("",ifile,UNLOAD,&error);
      return;
  }
  for(i=0;i<npartpoints;i++){
    parti->isprink[i]=0;
  }
  for(i=0;i<npartframes;i++){
    parti->sprframe[i]=0;
  }


  PRINTF("Loading particle data: %s\n",file);
  FORTgetdata1(&file_unit,&parttype,&error);
  if(partfilenum>=0&&partfilenum<npartinfo){
    partshortlabel=partinfo[partfilenum].label.shortlabel;
    partunitlabel=partinfo[partfilenum].label.unit;
  }
  else{
    partshortlabel=emptylabel;
    partunitlabel=emptylabel;
  }
  if(error!=0){
    *errorcode=1;
    fprintf(stderr,"*** Error: problem reading %s\n",file);
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  xbox=DENORMALIZE_X(xbar);
  ybox=DENORMALIZE_Y(ybar);
  zbox=DENORMALIZE_Z(zbar);
  offset_x=meshi->offset[XXX];
  offset_y=meshi->offset[YYY];
  offset_z=meshi->offset[ZZZ];
  FORTgetdata2(&file_unit,
    parti->xparts,parti->yparts,parti->zparts,
    parti->tpart,&parti->droplet_type,parti->isprink,
    tspr,parti->bframe,parti->sframe,parti->sprframe,parti->times,&nspr,&npartpoints,&npartframes,&parti->ntimes,
    &settmin_p,&settmax_p,&tmin_p,&tmax_p,&partframestep,&partpointstep, 
    &xbar0, &xbox, &ybar0, &ybox, &zbar0, &zbox,
    &offset_x, &offset_y, &offset_z, &redirect,
    &error,1);
  if(error!=0||parti->ntimes==0){
    if(error!=0)fprintf(stderr,"*** Error: problem reading %s\n",file);
    *errorcode=1;
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  
  nmax = parti->bframe[parti->ntimes-1]+parti->sframe[parti->ntimes-1];
  PRINTF("loaded: points=%i, size=%i KBytes, frames=%i\n",nmax,nmax*bytesperpoint/1024,parti->ntimes);

  if(parttype==-1||parttype==-3){
    parti->particle_type=1;  /*  only color temperature */
  }
  else{
    parti->particle_type=0;
  }
  havesprinkpart=0;
  skip_local=0;
  if(staticframe0==1)skip_local=parti->sframe[0];
  tcopy=parti->tpart+skip_local;
  tmin_global=1000000000.0;
  tmax_global=-tmin_global;
  isprinkcopy=parti->isprink;
  for(n=skip_local;n<nmax;n++){
    if(*isprinkcopy==0&&parti->particle_type==0){
      tcopy++; 
      isprinkcopy++;
      continue;
    }
    if(*isprinkcopy==1&&parti->droplet_type==0){
      tcopy++; 
      isprinkcopy++;
      havesprinkpart=1;
      continue;
    }
    if(*tcopy<tmin_global)tmin_global=*tcopy;
    if(*tcopy>tmax_global)tmax_global=*tcopy;
    if(*isprinkcopy==1){
      havesprinkpart=1;
    }
    tcopy++; 
    isprinkcopy++;
  }
  /* convert particle temperatures into integers pointing to an rgb color table */

  PRINTF("computing particle color levels \n");
  if(parti->particle_type!=0||parti->droplet_type!=0){
    adjustpartbounds(parti->tpart,parti->particle_type,parti->droplet_type,parti->isprink,
      skip_local,nmax,setpartmin,&tmin_global,setpartmax,&tmax_global);
  }
  if(setpartmin == SET_MIN){
    tmin_global = partmin;
  }
  if(setpartmax == SET_MAX){
    tmax_global = partmax;
  }
  partmin=tmin_global;
  partmax=tmax_global;
  if(NewMemory((void **)&colorlabelpart,MAXRGB*sizeof(char *))==0){
    FORTclosefortranfile(&file_unit);
    readpart("",ifile,UNLOAD,&error);
    *errorcode=1;
    return;
  }
  for(n=0;n<MAXRGB;n++){colorlabelpart[n]=NULL;}
  for(n=0;n<nrgb;n++){
    if(NewMemory((void **)&colorlabelpart[n],11)==0){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpart("",ifile,UNLOAD,&error);
      return;
    }
  }
  getPartColors(parti->tpart, skip_local, nmax, 
    parti->itpart,parti->isprink,parti->particle_type,parti->droplet_type,
    &tmin_global, &tmax_global, nrgb, colorlabelpart, partscale,partlevels256);
  for(n=0;n<skip_local;n++){
    parti->itpart[n]=0;
  }
  FREEMEMORY(parti->tpart);
  FREEMEMORY(parti->isprink);
  if(ResizeMemory((void **)&parti->xparts,nmax*sizeof(short))==0||
     ResizeMemory((void **)&parti->yparts,nmax*sizeof(short))==0||
     ResizeMemory((void **)&parti->zparts,nmax*sizeof(short))==0){
    FORTclosefortranfile(&file_unit);
    *errorcode=1;
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  if(ResizeMemory((void **)&parti->itpart,nmax*sizeof(unsigned char))==0){
    FORTclosefortranfile(&file_unit);
    *errorcode=1;
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  updateglui();

#ifdef pp_MEMPRINT
  PRINTF("After particle file load: \n");
  PrintMemoryInfo;
#endif
  if(partinfo[ifile].evac==0){
    ReadPartFile=1;
  }
  else{
    ReadEvacFile=1;
  }
  if(partinfo[ifile].evac==0){
    visParticles=1;
  }
  else{
    visEvac=1;
  }
  partinfo[ifile].loaded=1;
  partinfo[ifile].display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  Update_Times();
  updatemenu=1;
  Idle_CB();

  glutPostRedisplay();
}

/* ----------------------- drawselect_avatars ----------------------------- */

void drawselect_avatars(void){
  int i;

  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->loaded==0||parti->display==0)continue;
    if(parti->evac==1){
      drawEvac(parti);
      SNIFF_ERRORS("after drawEvac");
    }
  }
}

/* ------------------ drawEvac ------------------------ */

void drawEvac(const partdata *parti){
  drawPart(parti);
}

/* ------------------ get_evacpart_color ------------------------ */

int get_evacpart_color(float **color_handle,part5data *datacopy, int show_default, int j, int itype){
  int is_human_color;
  float *colorptr;
  unsigned char *color;
  int showcolor;

  showcolor=1;
  is_human_color=0;
  if(current_property!=NULL&&strcmp(current_property->label->longlabel,"HUMAN_COLOR")==0&&navatar_colors>0){
    is_human_color=1;
  }
  if(show_default==1){
    colorptr=datacopy->partclassbase->rgb;
  }
  else{
    color = datacopy->irvals+itype*datacopy->npoints;
   if(is_human_color==1){
      colorptr = avatar_colors + 3*color[j];
    }
    else{
      colorptr=rgb_full[color[j]];
    } 
    if(current_property!=NULL&&(color[j]>current_property->imax||color[j]<current_property->imin))showcolor=0;
  }
  *color_handle=colorptr;
  return showcolor;
}

/* ------------------ drawPart5 ------------------------ */

void drawPart5(const partdata *parti){
  int ipframe;
  part5data *datacopy,*datapast;
  int nclasses;
  int i,j;
  int offset_terrain;
  propdata *prop;

  if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
    offset_terrain=1;
  }
  else{
    offset_terrain=0;
  }

  if(current_property==NULL)return;
  ipframe=parti->itime;
  if(ipframe<0){
    ipframe=0;
    } //xxx need to check this - why is ipframe < 0 ???
  nclasses = parti->nclasses;
  datacopy = parti->data5+nclasses*ipframe;
  CheckMemory;
  if(part5show==1){
    if(streak5show==0||(streak5show==1&&showstreakhead==1)){
      for(i=0;i<parti->nclasses;i++){
        short *sx, *sy, *sz;
        float *angle, *width, *depth, *height;
        unsigned char *vis, *color;
        part5class *partclassi;
        int partclass_index, itype, vistype, class_vis;
        int show_default;

        partclassi = parti->partclassptr[i];
        partclass_index = partclassi - partclassinfo;

        vistype=current_property->class_present[partclass_index];
        class_vis=current_property->class_vis[partclass_index];


        if(vistype==0||datacopy->npoints<=0||(vistype==1&&class_vis==0)){
          if(show_tracers_always==0||partclassi->ntypes>2){
            datacopy++;
            continue;
          }
        }
        itype = current_property->class_types[partclass_index];

        show_default = 0;
        if(itype==-1||(show_tracers_always==1&&partclassi->ntypes<=2)){
          show_default=1;
        }

        sx = datacopy->sx;
        sy = datacopy->sy;
        sz = datacopy->sz;
        vis = datacopy->vis_part;
        if(parti->evac==1){
          int avatar_type=0;

          angle=datacopy->avatar_angle;
          width=datacopy->avatar_width;
          depth=datacopy->avatar_depth;
          height=datacopy->avatar_height;
          CheckMemory;

          avatar_type=0;
          prop=datacopy->partclassbase->prop;
          if(prop==NULL)prop=prop_evacdefault;
          if(iavatar_evac!=-1)avatar_type=iavatar_evac;
          for(j=0;j<datacopy->npoints;j++){
            float az_angle;
            float *colorptr;

            if(vis[j]==1){
              int save_use_displaylist;

              glPushMatrix();
              glTranslatef(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]-SCALE2SMV(parti->zoffset));
              if(select_avatar==1&&selected_avatar_tag>0&&selected_avatar_tag==datacopy->tags[j]){
                selected_avatar_pos[0]=xplts[sx[j]];
                selected_avatar_pos[1]=yplts[sy[j]];
                selected_avatar_pos[2]=zplts[sz[j]];
                selected_avatar_angle = datacopy->avatar_angle[j];
              }
              glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
                 
              az_angle=angle[j];
              glRotatef(az_angle,0.0,0.0,1.0);

              get_evacpart_color(&colorptr,datacopy,show_default,j,itype);
              
              //  :W :D :H1 :SX :SY :SZ :R :G :B :HX :HY :HZ
              //  class color: rgbobject[0], rgbobject[1], rgbobject[2]

              if(prop!=NULL){
                int n;
                sv_object_frame *obj_frame;
                tokendata **evac_tokens,*evac_token;

                obj_frame=prop->smv_object->obj_frames[0];
                evac_tokens = obj_frame->evac_tokens;
                obj_frame->nevac_tokens=12;
                
                n=0;

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=width[j]; //:W

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=depth[j]; //:D

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=1.0;//:H1

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=1.0;//:SX

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=1.0;//:SY

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=height[j];  //:SZ

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=255*colorptr[0]; //:R

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=255*colorptr[1];//:G

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=255*colorptr[2];//:B

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=0.0;//:HX

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=0.0;//:HY

                evac_token=evac_tokens[n++];
                if(evac_token!=NULL)evac_token->evac_var=height[j]/2.0; //:HZ
                prop->draw_evac=1;
              }

              save_use_displaylist=avatar_types[avatar_type]->use_displaylist;
              if(select_avatar==1&&show_mode==SELECTOBJECT){
                int tagval;

                avatar_types[avatar_type]->select_mode=1;
                select_device_color_ptr=select_device_color;
                tagval=datacopy->tags[j];
                select_device_color[0]=tagval>>(ngreenbits+nbluebits);
                select_device_color[1]=tagval>>nbluebits;
                select_device_color[2]=tagval&rgbmask[nbluebits-1];
                avatar_types[avatar_type]->use_displaylist=0;
              }
              else{
                if(selected_avatar_tag>0&&select_avatar==1&&datacopy->tags[j]==selected_avatar_tag){
                  select_device_color_ptr=select_device_color;
                  select_device_color[0]=255;
                  select_device_color[1]=0;
                  select_device_color[2]=0;
                  avatar_types[avatar_type]->use_displaylist=0;
                }
                else{
                  select_device_color_ptr=NULL;
                  avatar_types[avatar_type]->select_mode=0;
                }
              }
              copy_dep_vals(partclassi,datacopy,colorptr,prop,j);
              draw_SVOBJECT(avatar_types[avatar_type],0,prop,0,NULL,0);
              select_device_color_ptr=NULL;
              avatar_types[avatar_type]->use_displaylist=save_use_displaylist;
              glPopMatrix();
            }
          }
          SNIFF_ERRORS("after draw in Evac");
        }
        else{
          glPointSize(partpointsize);
          if(offset_terrain==0){

            // *** draw particles as points
           
            if(datacopy->partclassbase->vis_type==PART_POINTS){
              glBegin(GL_POINTS);
              if(show_default==1){
                glColor4fv(datacopy->partclassbase->rgb);
                for(j=0;j<datacopy->npoints;j++){
                  if(vis[j]==1){
                    glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
                  }
                }
              }
              else{
                color=datacopy->irvals+itype*datacopy->npoints;
                for(j=0;j<datacopy->npoints;j++){
                  if(vis[j]==1){
                    if(current_property!=NULL&&(color[j]>current_property->imax||color[j]<current_property->imin))continue;
                    glColor4fv(rgb_full[color[j]]);
                    glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
                  }
                }
              }
              glEnd();
            }

            // *** draw particles using smokeview object

            if(datacopy->partclassbase->vis_type==PART_SMV_DEVICE){
              for(j=0;j<datacopy->npoints;j++){
                float *colorptr;

                if(vis[j]!=1)continue;
                  
                glPushMatrix();
                glTranslatef(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);

                glRotatef(-datacopy->partclassbase->elevation,0.0,1.0,0.0);
                glRotatef( datacopy->partclassbase->azimuth,  0.0,0.0,1.0);

              //  0->2   color
              //  3      diameter
              //  4      length

                if(show_default==1){
                  colorptr=datacopy->partclassbase->rgb;
                }
                else{
                  color = datacopy->irvals+itype*datacopy->npoints;
                  colorptr=rgb_full[color[j]];
                }

                prop=datacopy->partclassbase->prop;
                copy_dep_vals(partclassi,datacopy,colorptr,prop,j);
                glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));

                partfacedir[0]=xbar0+SCALE2SMV(world_eyepos[0])-xplts[sx[j]];
                partfacedir[1]=ybar0+SCALE2SMV(world_eyepos[1])-yplts[sy[j]];
                partfacedir[2]=zbar0+SCALE2SMV(world_eyepos[2])-zplts[sz[j]];
                
                draw_SVOBJECT(prop->smv_object,0,prop,0,NULL,0);
                glPopMatrix();
              }
            }

            // *** draw particle as lines

            if(datacopy->partclassbase->vis_type==PART_LINES
              &&((datacopy->dsx!=NULL&&datacopy->dsy!=NULL&&datacopy->dsz!=NULL)||datacopy->partclassbase->device_name!=NULL)
              ){
              float *dxv, *dyv, *dzv;
              float dx, dy, dz;
              int flag=0;

              if(datacopy->dsx!=NULL&&datacopy->dsy!=NULL&&datacopy->dsz!=NULL){
                flag=1;
                dxv = datacopy->dsx;
                dyv = datacopy->dsy;
                dzv = datacopy->dsz;
              }
              else{
                dx = datacopy->partclassbase->dx;
                dy = datacopy->partclassbase->dy;
                dz = datacopy->partclassbase->dz;
              }
              glBegin(GL_LINES);
              if(show_default==1){
                glColor4fv(datacopy->partclassbase->rgb);
                for(j=0;j<datacopy->npoints;j++){
                  if(vis[j]==1){
                    if(flag==1){
                      dx = dxv[j];
                      dy = dyv[j];
                      dz = dzv[j];
                    }
                    glVertex3f(xplts[sx[j]]-dx,yplts[sy[j]]-dy,zplts[sz[j]]-dz);
                    glVertex3f(xplts[sx[j]]+dx,yplts[sy[j]]+dy,zplts[sz[j]]+dz);
                  }
                }
              }
              else{
                color=datacopy->irvals+itype*datacopy->npoints;
                for(j=0;j<datacopy->npoints;j++){
                  if(vis[j]==1){
                    glColor4fv(rgb_full[color[j]]);
                    if(flag==1){
                      dx = dxv[j];
                      dy = dyv[j];
                      dz = dzv[j];
                    }
                    glVertex3f(xplts[sx[j]]-dx,yplts[sy[j]]-dy,zplts[sz[j]]-dz);
                    glVertex3f(xplts[sx[j]]+dx,yplts[sy[j]]+dy,zplts[sz[j]]+dz);
                  }
                }
              }
              glEnd();
            }
          }
          else{
            glBegin(GL_POINTS);
            if(show_default==1){
              glColor4fv(datacopy->partclassbase->rgb);
              for(j=0;j<datacopy->npoints;j++){
                float zoffset;
                float xx, yy, zz;
                int loc;

                xx = xplts[sx[j]];
                yy = yplts[sy[j]];
                zz = zplts[sz[j]];

                zoffset = get_zcell_val_offset(meshinfo,xx,yy,&loc);
                if(vis[j]==1)glVertex3f(xx,yy,zz+zoffset);
              }
            }
            else{
              color=datacopy->irvals+itype*datacopy->npoints;
              for(j=0;j<datacopy->npoints;j++){
                if(vis[j]==1){
                  glColor4fv(rgb_full[color[j]]);
                  glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
                }
              }
            }
            glEnd();
          }
        }

        datacopy++;
      }
    }
  }

  // draw streak lines

  datacopy = parti->data5+nclasses*ipframe;

  if(streak5show==1){
  for(i=0;i<parti->nclasses;i++){
    short *sx, *sy, *sz;
    short *sxx, *syy, *szz;
    unsigned char *vis;
    int k;
    int show_default;
    float *colorptr;

    part5class *partclassi;
    int partclass_index, itype, vistype, class_vis;

    partclassi = parti->partclassptr[i];
    partclass_index = partclassi - partclassinfo;

    vistype=current_property->class_present[partclass_index];
    class_vis=current_property->class_vis[partclass_index];

    if(vistype==0||datacopy->npoints<=0||(vistype==1&&class_vis==0)){
      if(show_tracers_always==0||partclassi->ntypes>2){
        datacopy++;
        continue;
      }
    }
    itype = current_property->class_types[partclass_index];

    show_default=0;
    if(itype==-1||(show_tracers_always==1&&partclassi->ntypes<=2)){
      show_default=1;
    }

    sx = datacopy->sx;
    sy = datacopy->sy;
    sz = datacopy->sz;
    vis = datacopy->vis_part;

    if(show_default==1){

      // draw the streak line

      get_evacpart_color(&colorptr,datacopy,show_default,0,itype);
      glColor4fv(colorptr);

      glLineWidth(streaklinewidth);
      for(j=0;j<datacopy->npoints;j++){
        int tagval;
        
        tagval=datacopy->tags[j];
        if(vis[j]==0)continue;
        glBegin(GL_LINE_STRIP);
        glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
        for(k=1;k<streak5step;k++){
          int jj;

          if(ipframe-k<0)break;
          datapast = parti->data5+nclasses*(ipframe-k)+i;
          jj = get_tagindex(parti,&datapast,tagval);
          if(jj<0)break;
          sxx = datapast->sx;
          syy = datapast->sy;
          szz = datapast->sz;
          glVertex3f(xplts[sxx[jj]],yplts[syy[jj]],zplts[szz[jj]]);
        }
        glEnd();
      }

      // draw the dot at the end of the streak line
    }
    else{

      // draw the streak line

      for(j=0;j<datacopy->npoints;j++){
        int tagval;

        tagval=datacopy->tags[j];
        if(vis[j]==0)continue;
        if(get_evacpart_color(&colorptr,datacopy,show_default,j,itype)==0)continue;
        
        glBegin(GL_LINE_STRIP);
        glColor4fv(colorptr);
        glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
        for(k=1;k<streak5step;k++){
          int jj;

          if(ipframe-k<0)break;
          datapast = parti->data5+nclasses*(ipframe-k)+i;
          jj = get_tagindex(parti,&datapast,tagval);
          if(jj<0||datapast->irvals==NULL)break;
          sxx = datapast->sx;
          syy = datapast->sy;
          szz = datapast->sz;

          get_evacpart_color(&colorptr,datacopy,show_default,jj,itype);
          glColor4fv(colorptr);
          glVertex3f(xplts[sxx[jj]],yplts[syy[jj]],zplts[szz[jj]]);
        }
        glEnd();
      }
    }

    datacopy++;
  }
  }

}

/* ------------------ copy_dep_vals ------------------------ */

void copy_dep_vals(part5class *partclassi, part5data *datacopy, float *colorptr, propdata *prop, int j){
  int ii;
  int ndep_vals;
  float *dep_vals;

  if(prop==NULL)return;
  dep_vals=partclassi->fvars_dep;
  ndep_vals=partclassi->nvars_dep;
  for(ii=0;ii<partclassi->nvars_dep-3;ii++){

    unsigned char *var_type;
    unsigned char color_index;
    part5prop *varprop;
    float valmin, valmax;
    char *shortlabel;
    flowlabels *label;

    shortlabel=NULL;
    varprop=NULL;
    label = datacopy->partclassbase->labels+ii+2;
    if(label!=NULL)shortlabel=label->shortlabel;
    if(shortlabel!=NULL)varprop = get_part5prop_s(shortlabel);
    if(varprop!=NULL){
      var_type = datacopy->irvals + ii*datacopy->npoints;
      color_index=var_type[j];
      valmin=varprop->valmin;
      valmax=varprop->valmax;
      dep_vals[ii]=valmin + color_index*(valmax-valmin)/255.0;
    }
    else{
      dep_vals[ii]=1.0;
    }
  }
  
  dep_vals[ndep_vals-3]=colorptr[0]*255;
  dep_vals[ndep_vals-2]=colorptr[1]*255;
  dep_vals[ndep_vals-1]=colorptr[2]*255;
  prop->nvars_dep=partclassi->nvars_dep;
  prop->smv_object->visible=1;
  for(ii=0;ii<prop->nvars_dep;ii++){
    prop->fvars_dep[ii]=partclassi->fvars_dep[ii];
  }
  prop->nvars_dep=partclassi->nvars_dep;
  for(ii=0;ii<partclassi->nvars_dep;ii++){
    prop->vars_dep_index[ii]=partclassi->vars_dep_index[ii];
  }
  prop->tag_number = datacopy->tags[j];
}

/* ------------------ drawPart ------------------------ */

void drawPart(const partdata *parti){
  short *xpoints, *ypoints, *zpoints;
  unsigned char *itpoint=NULL;


  int n;
  int nsmokepoints, nsprpoints;
  int ipframe;
  int droplet_type, particle_type;
  float *rgb_smoke, *rgb_ismoke;

  ipframe=parti->itime;
  rgb_smoke = rgb_part;

  if(parti->times[0]>global_times[itimes])return;

  if(parti->version==1){
    drawPart5(parti);
    return;
  }

  droplet_type = parti->droplet_type;
  particle_type = parti->particle_type;

  /* define the data locations to look at */
  
    xpoints = parti->xparts + parti->bframe[ipframe];
    ypoints = parti->yparts + parti->bframe[ipframe];
    zpoints = parti->zparts + parti->bframe[ipframe];

  /* isprinkframe = isprink + bframe[ipframe];*/

  itpoint = parti->itpart + parti->bframe[ipframe];
  nsprpoints = parti->sprframe[ipframe];
  nsmokepoints = parti->sframe[ipframe]-nsprpoints;

  glPointSize(partpointsize);
  glBegin(GL_POINTS);
  if(parti->version==0){
    if(visSmokePart!=0){
      if(particle_type==0){
        for (n = 0; n < nsmokepoints; n++){
          glColor4fv(rgb[itpoint[n]]);
          glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
        }
      }
      else{
        for (n = 0; n < nsmokepoints; n++){
          rgb_ismoke = rgb_smoke + 4*itpoint[n];
          if(rgb_ismoke[3]>0.5){
            glColor4fv(rgb_ismoke);
            glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
          }
        }
      }
    }
    if(visSprinkPart==1){
      if(droplet_type==0){
        glColor4fv(rgb[rgb_blue]);
        for (n = nsmokepoints; n < nsmokepoints+nsprpoints; n++){
          glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
        }
      }
      else{
        for (n = nsmokepoints; n < nsmokepoints+nsprpoints; n++){
          glColor4fv(rgb_full[itpoint[n]]);
          glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
        }
      }
    }
  }

  glEnd();



}

/* ------------------ drawStaticPart ------------------------ */

void drawStaticPart(const partdata *parti){

  short *xpoints, *ypoints, *zpoints;

  int n;
  int nsmokepoints, nsprpoints;
  int ipframe;
  
  /* define the data locations to look at */

  if(parti->version!=0)return;
  ipframe=0;
  xpoints = parti->xparts + parti->bframe[ipframe];
  ypoints = parti->yparts + parti->bframe[ipframe];
  zpoints = parti->zparts + parti->bframe[ipframe];

  /* isprinkframe = isprink + bframe[ipframe];*/

  nsprpoints = parti->sprframe[ipframe];
  nsmokepoints = parti->sframe[ipframe]-nsprpoints;

  glPointSize(partpointsize);

  glColor4fv(static_color);
  glBegin(GL_POINTS);
  for(n=0;n<nsmokepoints+nsprpoints;n++){
    glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
  }
  glEnd();
}

/* ------------------ tagscompare ------------------------ */

int tagscompare( const void *arg1, const void *arg2 ){
  int i, j;

  i = *(int *)arg1;
  j = *(int *)arg2;
  if(i<j)return -1;
  if(i>j)return 1;
  return 0;
  
}

/* ------------------ partcompare ------------------------ */

int partcompare( const void *arg1, const void *arg2 ){
  partdata *parti, *partj;
  int i, j;

  i = *(int *)arg1;
  j = *(int *)arg2;
  
  parti = partinfo + i;
  partj = partinfo + j;

  if(parti->version==1){
    if(parti->blocknumber<partj->blocknumber)return -1;
    if(parti->blocknumber>partj->blocknumber)return 1;
  }
  else{
    if(strcmp(parti->label.longlabel,partj->label.longlabel)<0)return -1;
    if(strcmp(parti->label.longlabel,partj->label.longlabel)>0)return 1;
    if(parti->blocknumber<partj->blocknumber)return -1;
    if(parti->blocknumber>partj->blocknumber)return 1;
  }
  return 0;
}

/* ------------------ updatepartmenulabels ------------------------ */

void updatepartmenulabels(void){
  int i;
  partdata *parti;
  char label[128];
  int lenlabel;

  if(npartinfo>0){
    FREEMEMORY(partorderindex);
    NewMemory((void **)&partorderindex,sizeof(int)*npartinfo);
    for(i=0;i<npartinfo;i++){
      partorderindex[i]=i;
    }
    qsort( (int *)partorderindex, (size_t)npartinfo, sizeof(int), partcompare );

    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      STRCPY(parti->menulabel,"");
      if(parti->evac==1){
        STRCAT(parti->menulabel,"humans");
      }
      else{
        if(parti->version==1){
          STRCAT(parti->menulabel,"particles");
        }
        else{
          STRCAT(parti->menulabel,parti->label.longlabel);
        }
      }
      lenlabel=strlen(parti->menulabel);
      if(nmeshes>1){
	    mesh *partmesh;

		partmesh = meshinfo + parti->blocknumber;
        sprintf(label,"%s",partmesh->label);
        if(lenlabel>0)STRCAT(parti->menulabel,", ");
        STRCAT(parti->menulabel,label);
      }
      if(showfiles==1||lenlabel==0){
        if(nmeshes>1||lenlabel>0)STRCAT(parti->menulabel,", ");
        STRCAT(parti->menulabel,parti->file);
      }
    } 
  }


}
