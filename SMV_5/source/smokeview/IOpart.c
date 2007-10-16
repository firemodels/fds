// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <time.h>
#include <sys/stat.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char IOpart_revision[]="$Revision$";

int tagscompare( const void *arg1, const void *arg2 );

void ParticlePropShowMenu(int val);
void PART_CB_INIT(void);
void update_all_partvis(particle *parti);
void update_partvis(int first_frame,part5data *datacopy, int nclasses);
int get_tagindex(part5data *data, int tagval);
void endian_switch(void *val, int nval);
#define READPASS 1
#define READFAIL 0
#define SEEKPASS 0

#define FORTPART5READORIG(var,size) fseek(PART5FILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,PART5FILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           fseek(PART5FILE,4,SEEK_CUR)


#define FORTPART5READ(var,size) \
returncode=READPASS;\
fseek(PART5FILE,4,SEEK_CUR);if(ferror(PART5FILE)==1||feof(PART5FILE)==1)returncode=READFAIL;\
if(returncode==READPASS){\
  fread(var,4,size,PART5FILE);\
  if(ferror(PART5FILE)==1||feof(PART5FILE)==1)returncode=READFAIL;\
}\
if(returncode==READPASS){\
  if(endianswitch==1)endian_switch(var,size);\
  fseek(PART5FILE,4,SEEK_CUR);\
  if(ferror(PART5FILE)==1||feof(PART5FILE)==1)returncode=READFAIL;\
}



/* ------------------ freepart5data ------------------------ */

void freepart5data(part5data *datacopy){
  FREEMEMORY(datacopy->cvals);
  FREEMEMORY(datacopy->ivals);
  FREEMEMORY(datacopy->sx);
  FREEMEMORY(datacopy->sy);
  FREEMEMORY(datacopy->sz);
  FREEMEMORY(datacopy->tags);
  FREEMEMORY(datacopy->sort_tags);
  FREEMEMORY(datacopy->vis_part);
  FREEMEMORY(datacopy->rvals);
  FREEMEMORY(datacopy->irvals);
  FREEMEMORY(datacopy->ivals);
}

/* ------------------ freeallpart5data ------------------------ */

void freeallpart5data(particle *parti){
  int i;
  part5data *datacopy;

  if(parti->data5==NULL)return;
  datacopy = parti->data5;
  for(i=0;i<parti->nframes*parti->nclasses;i++){
    freepart5data(datacopy);
    datacopy++;
  }
  FREEMEMORY(parti->data5);
}

/* ------------------ initpart5data ------------------------ */

void initpart5data(part5data *datacopy, part5class *partclassi){
  datacopy->cvals=NULL;
  datacopy->ivals=NULL;
  datacopy->partclassbase=partclassi;
  datacopy->sx=NULL;
  datacopy->sy=NULL;
  datacopy->sz=NULL;
  datacopy->tags=NULL;
  datacopy->vis_part=NULL;
  datacopy->sort_tags=NULL;
  datacopy->rvals=NULL;
  datacopy->irvals=NULL;
  datacopy->ivals=NULL;
}


/* ------------------ getpart5data ------------------------ */

void getpart5data(particle *parti, int partframestep, int partpointstep){
  FILE *PART5FILE;
  int one;
  int endianswitch=0;
  int version;
  int nclasses;
  int i;
  int skip;
  size_t returncode;
  float time;
  int nparts;
  int *numtypes=NULL,*numtypescopy, *numpoints=NULL;
  int numtypes_temp[2];
  char *reg_file;
  part5data *datacopy;
  int count=-1;
  int count2=-1;
  int first_frame=1;

  reg_file=parti->reg_file;

  PART5FILE=fopen(reg_file,"rb");
  if(PART5FILE==NULL)return;

  fseek(PART5FILE,4,SEEK_CUR);fread(&one,4,1,PART5FILE);fseek(PART5FILE,4,SEEK_CUR);
  if(one!=1)endianswitch=1;

  FORTPART5READ(&version,1);if(returncode==0)goto wrapup;
  FORTPART5READ(&nclasses,1);if(returncode==0)goto wrapup;
  NewMemory((void **)&numtypes,2*nclasses*sizeof(int));
  NewMemory((void **)&numpoints,nclasses*sizeof(int));
  numtypescopy=numtypes;
  numtypes_temp[0]=0;
  numtypes_temp[1]=0;
  for(i=0;i<nclasses;i++){
    FORTPART5READ(numtypes_temp,2);
    if(returncode==0)goto wrapup;
    *numtypescopy++=numtypes_temp[0];
    *numtypescopy++=numtypes_temp[1];
    skip = 2*(numtypes_temp[0]+numtypes_temp[1])*(8 + 30);
    returncode=fseek(PART5FILE,skip,SEEK_CUR);
    if(returncode!=0)goto wrapup;
  }

  datacopy = parti->data5;
  for(;;){
    int doit;

    count++;
    if(count>=parti->nframes)break;
    if(count%partframestep==0){
      count2++;
      doit=1;
    }
    else{
      doit=0;
    }

    FORTPART5READ(&time,1);
    if(returncode==0)break;
    printf("particle time=%f",time);
    if(doit==1){
      parti->ptimes[count2]=time;
    }
    for(i=0;i<nclasses;i++){
      part5class *partclassi;
      int factor=256*128;

      partclassi = partclassinfo + i;
      FORTPART5READ(&nparts,1);
      if(returncode==0)goto wrapup;
      numpoints[i]=nparts;
      skip=0;
      if(doit==1){
        short *sx, *sy, *sz;
        float *xyz;
        int j;

        FORTPART5READ(partclassi->xyz,3*nparts);
        if(nparts>0){
          if(returncode==0)goto wrapup;
          xyz = partclassi->xyz;
          sx = datacopy->sx;
          sy = datacopy->sy;
          sz = datacopy->sz;


          for(j=0;j<nparts;j++){
            float xx, yy, zz;

            xx = (xyz[         j]-xbar0)/xyzmaxdiff;
            xx /= xbar;

            yy = (xyz[  nparts+j]-ybar0)/xyzmaxdiff;
            yy /= ybar;

            zz = (xyz[2*nparts+j]-zbar0)/xyzmaxdiff;
            zz /= zbar;

            sx[j] = factor*xx;
            sy[j] = factor*yy;
            sz[j] = factor*zz;
          }
        }
      }
      else{
        skip = 4 + 3*4*nparts + 4;  
      }
      if(doit==1){
        int *sort_tags;
        int j;

        sort_tags=datacopy->sort_tags;
        FORTPART5READ(datacopy->tags,nparts);
        if(nparts>0){
          if(returncode==0)goto wrapup;
          for(j=0;j<nparts;j++){
            sort_tags[2*j]=datacopy->tags[j];
            sort_tags[2*j+1]=j;
          }
          qsort( sort_tags, (size_t)nparts, 2*sizeof(int), tagscompare );
      //  update_partvis(first_frame,datacopy,nclasses);
        }
      }
      else{
        skip = 4 + 4*nparts + 4;  // skip over tag for now
      }
      if(numtypes[2*i]>0){
        //skip += 4 + 4*nparts*numtypes[2*i] + 4;  // skip over vals for now
        FORTPART5READ(datacopy->rvals,nparts*numtypes[2*i]);
        if(returncode==0)goto wrapup;
      }
      if(numtypes[2*i+1]>0){
        skip += 4 + 4*nparts*numtypes[2*i+1] + 4;
      }

      
      returncode=0;
      if(skip>0){
        returncode=fseek(PART5FILE,skip,SEEK_CUR);
        if(returncode!=0)goto wrapup;
      }
      datacopy++;
    }
    if(first_frame==1)first_frame=0;
    printf(" completed\n");

  }
wrapup:
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
    printf("label=%s min=%f max=%f\n",propi->label->longlabel,propi->valmin,propi->valmax);
    printf("   glbmin=%f glbmax=%f\n",propi->global_min,propi->global_max);
    printf("   permin=%f permax=%f\n",propi->percentile_min,propi->percentile_max);
    printf("\n");
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
          propi->label=flowlabel;

          propi->setvalmin=0;
          propi->setvalmax=0;
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
        classprop->class_present[i]=1;
        classprop->class_types[i]=j-2;
      }
    }
  }
}

/* ------------------ update_partvis ------------------------ */

void update_all_partvis(particle *parti){
  part5data *datacopy;
  int i,j;
  int firstframe=1;

  datacopy = parti->data5;
  for(i=0;i<parti->nframes;i++){
    for(j=0;j<parti->nclasses;j++){
      update_partvis(firstframe,datacopy,parti->nclasses);
      datacopy++;
    }
    if(firstframe==1)firstframe=0;
  }
}


/* ------------------ update_partvis ------------------------ */

void update_partvis(int first_frame,part5data *datacopy, int nclasses){
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

    datalast = datacopy-nclasses;
    for(ii=0;ii<nparts;ii++){
      int tag_index;

      tag_index = get_tagindex(datalast,datacopy->tags[ii]);
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

int get_tagindex(part5data *data, int tagval){
  int *returnval;

  if(data->npoints==0)return 1;
  ASSERT(data->sort_tags!=NULL);
  returnval=bsearch(&tagval,data->sort_tags,data->npoints,2*sizeof(int),tagscompare);
  if(returnval==NULL)return -1;
  return *(returnval+1);
}

/* ------------------ setpart5sizefile ------------------------ */

void create_part5sizefile(char *reg_file, char *size_file){
  FILE *size_stream, *PART5FILE;
  int one;
  int endianswitch=0;
  int version;
  int nclasses;
  int i;
  int skip;
  size_t returncode;
  float time;
  int nparts;
  int *numtypes=NULL,*numtypescopy, *numpoints=NULL;
  int numtypes_temp[2];
  char buffer_out[1024];

  PART5FILE=fopen(reg_file,"rb");
  if(PART5FILE==NULL)return;
  size_stream=fopen(size_file,"w");
  if(size_stream==NULL){
    printf("*** warning:  unable to write to %s\n",size_file);
    fclose(PART5FILE);
    return;
  }

  fseek(PART5FILE,4,SEEK_CUR);fread(&one,4,1,PART5FILE);fseek(PART5FILE,4,SEEK_CUR);
  if(one!=1)endianswitch=1;

  FORTPART5READ(&version,1);
  FORTPART5READ(&nclasses,1);
  NewMemory((void **)&numtypes,2*nclasses*sizeof(int));
  NewMemory((void **)&numpoints,nclasses*sizeof(int));
  numtypescopy=numtypes;
  numtypes_temp[0]=0;
  numtypes_temp[1]=0;
  for(i=0;i<nclasses;i++){
    FORTPART5READ(numtypes_temp,2);
    *numtypescopy++=numtypes_temp[0];
    *numtypescopy++=numtypes_temp[1];
    skip = 2*(numtypes_temp[0]+numtypes_temp[1])*(8 + 30);
    fseek(PART5FILE,skip,SEEK_CUR);
  }

  for(;;){
    FORTPART5READ(&time,1);
    if(returncode==0)break;
    sprintf(buffer_out,"%f ",time);
     for(i=0;i<nclasses;i++){
      FORTPART5READ(&nparts,1);
      if(returncode==0)goto wrapup;
      numpoints[i]=nparts;
      skip = 4 + 4*nparts*3 + 4;
      skip += 4 + 4*nparts + 4;
      if(numtypes[2*i]>0)skip += 4 + 4*nparts*numtypes[2*i] + 4;
      if(numtypes[2*i+1]>0)skip += 4 + 4*nparts*numtypes[2*i+1] + 4;
      
      returncode=fseek(PART5FILE,skip,SEEK_CUR);
      if(returncode!=0)goto wrapup;
    }

    fprintf(size_stream,"%f\n",time);
#ifdef _DEBUG
    printf("%f\n",time);
#endif
    for(i=0;i<nclasses;i++){
      fprintf(size_stream,"  %i\n",numpoints[i]);
#ifdef _DEBUG
      printf("  %i\n",numpoints[i]);
#endif
    }
//    fprintf(size_stream,"\n");
  }
wrapup:
  FREEMEMORY(numtypes);
  fclose(PART5FILE);
  fclose(size_stream);
}

/* ------------------ getpart5header ------------------------ */

void getpart5header(particle *parti, int partframestep){
  FILE *stream;
  char buffer[256];
  float time;
  int count=-1;
  char *reg_file, *size_file;
  int i,j;
  int stat_sizefile, stat_regfile;
  struct stat stat_sizefile_buffer, stat_regfile_buffer;

  reg_file=parti->reg_file;
  size_file=parti->size_file;

  // if size file doesn't exist then generate it

  parti->nframes=0;

  stat_sizefile=stat(size_file,&stat_sizefile_buffer);
  stat_regfile=stat(reg_file,&stat_regfile_buffer);
  if(stat_regfile!=0)return;

  // create a size file if 1) the size does not exist
  //                       2) base file is newer than the size file
  if(stat_sizefile!=0||
    stat_regfile_buffer.st_mtime>stat_sizefile_buffer.st_mtime){
    //create_part5sizefile(reg_file,size_file);
      {
        int lenreg, lensize, error;

        trim(reg_file);
        trim(size_file);
        lenreg=strlen(reg_file);
        lensize=strlen(size_file);

        FORTfcreate_part5sizefile(reg_file,size_file, &error, 
          lenreg,lensize);
      }
  }
  
  stream=fopen(size_file,"r");
  if(stream==NULL)return;

    // pass 1: count frames

  for(;;){
    if(fgets(buffer,255,stream)==NULL)break;
    count++;
    if(count%partframestep!=0)continue;
    sscanf(buffer,"%f",&time);
    (parti->nframes)++;
  }
  rewind(stream);

  // allocate memory for number of time steps * number of classes

  NewMemory((void **)&parti->data5,parti->nclasses*parti->nframes*sizeof(part5data));
  NewMemory((void **)&parti->ptimes,parti->nframes*sizeof(float));


  // free memory for x, y, z frame data 

  for(i=0;i<parti->nclasses;i++){
    part5class *partclassi;

    partclassi = partclassinfo + i;
    FREEMEMORY(partclassi->xyz);
    partclassi->maxpoints=0;
  }

  // pass 2 - allocate memory for x, y, z frame data
  //          
  {
    part5data *datacopy;
    int fail;

    fail=0;
    datacopy=parti->data5;
    for(i=0;i<parti->nframes;i++){
      if(fgets(buffer,255,stream)==NULL){
        fail=1;
        break;
      }
      sscanf(buffer,"%f",&datacopy->time);
      for(j=0;j<parti->nclasses;j++){
        int npoints ,ntypes;

        part5class *partclassj;

        partclassj = partclassinfo + j;
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
          ntypes = datacopy->partclassbase->ntypes;
          if(ntypes>0){ //xxx need to check this fix (was ntypes>2)
            NewMemory((void **)&datacopy->rvals,(ntypes-0)*npoints*sizeof(float));  //xxx need to check this fix
            NewMemory((void **)&datacopy->irvals,(ntypes-0)*npoints*sizeof(unsigned char)); //xxx need to check this fix
          }
        }
        datacopy++;
      }
      if(fail==1)break;
    }
    if(fail==1)parti->nframes=i;
    fclose(stream);
  }

  // allocate memory for x, y, z and tag for the maximum frame size
  //           don't need to allocate memory for all frames

  for(i=0;i<parti->nclasses;i++){
    part5class *partclassi;

    partclassi = partclassinfo + i;
    if(partclassi->maxpoints>0)NewMemory((void **)&partclassi->xyz,3*partclassi->maxpoints*sizeof(float));
  }

}

/* ------------------ readpart5 ------------------------ */

void readpart5(char *file, int ifile, int flag, int *errorcode){
  size_t lenfile;
  int error=0;
  int ibar,jbar,kbar;
  particle *parti;
  int blocknumber;
  mesh *meshi;

  parti=partinfo+ifile;

  freeallpart5data(parti);

  blocknumber=parti->blocknumber;
  meshi=meshinfo+blocknumber;
  if(parti->loaded==0&&flag==UNLOAD)return;


  ibar=meshi->ibar;
  jbar=meshi->jbar;
  kbar=meshi->kbar;

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

  FREEMEMORY(parti->ptimes); 

  if(colorlabelpart!=NULL){
    int n;

    for(n=0;n<MAXRGB;n++){
      FREEMEMORY(colorlabelpart[n]);
    }
    FREEMEMORY(colorlabelpart);
  }

  if(flag==UNLOAD){
    updatetimes();
    updatemenu=1;
#ifdef _DEBUG
    printf("After particle file unload: ");
    PrintMemoryInfo;
#endif
    return;
  }

  lenfile = strlen(file);
  if(lenfile==0){
    readpart("",ifile,UNLOAD,&error);
    updatetimes();
    return;
  }
  
  printf("Sizing particle data: %s\n",file);
  getpart5header(parti, partframestep);

  offsetmax=5;
  if(offsetmax>ibar/4)offsetmax=ibar/4;
  if(offsetmax>jbar/4)offsetmax=jbar/4;
  if(offsetmax>kbar/4)offsetmax=kbar/4;
  
  printf("Loading particle data: %s\n",file);
  getpart5data(parti,partframestep,partpointstep);
  updateglui();

#ifdef _DEBUG
  printf("After particle file load: ");
  PrintMemoryInfo;
#endif
  if(parti->evac==0){
    ReadPartFile=1;
  }
  else{
    ReadEvacFile=1;
  }
  if(parti->evac==0){
    visSmoke=1;
  }
  else{
    visEvac=1;
  }
  /* convert particle temperatures into integers pointing to an rgb color table */

  printf("computing particle color levels \n");

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
#ifdef _DEBUG
  printf("After particle file load: ");
  PrintMemoryInfo;
#endif
  if(parti->evac==0){
    visSmoke=1;
    ReadPartFile=1;
  }
  else{
    visEvac=1;
    ReadEvacFile=1;
  }

  parttype=0;
  PART_CB_INIT();
  ParticlePropShowMenu(0);
  parti->loaded=1;
  parti->display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatetimes();
  updatemenu=1;
  IDLE();

  glutPostRedisplay();
}

/* ------------------ update_all_partvis2 ------------------------ */

void update_all_partvis2(void){
  particle *parti;
  int i;
  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    if(parti->loaded==1)update_all_partvis(parti);
  }
}

/* ------------------ readpart ------------------------ */

void readpart(char *file, int ifile, int flag, int *errorcode){
  int nmax, n, i;
  size_t lenfile;
  float *tcopy;
  unsigned char *isprinkcopy;
  int error=0;
  int bytesperpoint;
  int skip;
  int statfile,statfile2;
  struct stat statbuffer,statbuffer2;
  char partsizefile[1024],buffer[1024];
  FILE *sizefile;
  int readpartsize=1;
  int partpointstepold, partframestepold;
  int npartframes2, npartpoints2;
  size_t return_code;
  int ibar,jbar,kbar;
  int nb,nv;
  float xbox, ybox, zbox;
  particle *parti;
  int blocknumber;
  mesh *meshi;
  float offset_x, offset_y, offset_z;

  parti=partinfo+ifile;
  if(parti->version==1){
    readpart5(file,ifile,flag,errorcode);
    return;
  }
  blocknumber=parti->blocknumber;
  meshi=meshinfo+blocknumber;
  if(parti->loaded==0&&flag==UNLOAD)return;


  ibar=meshi->ibar;
  jbar=meshi->jbar;
  kbar=meshi->kbar;
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

  FREEMEMORY(parti->ptimes); 
  FREEMEMORY(parti->xpart);  FREEMEMORY(parti->ypart);  FREEMEMORY(parti->zpart);  
  FREEMEMORY(parti->xpartb); FREEMEMORY(parti->ypartb); FREEMEMORY(parti->zpartb); 
  FREEMEMORY(parti->xparts); FREEMEMORY(parti->yparts); FREEMEMORY(parti->zparts); 
  FREEMEMORY(parti->tpart);  FREEMEMORY(parti->itpart); 
  FREEMEMORY(parti->isprink);
  FREEMEMORY(parti->sframe); 
  FREEMEMORY(parti->bframe);
  FREEMEMORY(parti->sprframe)

  if(colorlabelpart!=NULL){
    for(n=0;n<MAXRGB;n++){FREEMEMORY(colorlabelpart[n]);}
    FREEMEMORY(colorlabelpart);
  }

  if(flag==UNLOAD){
    updatetimes();
    updatemenu=1;
#ifdef _DEBUG
    printf("After particle file unload: ");
    PrintMemoryInfo;
#endif
    return;
  }

  lenfile = strlen(file);
  if(lenfile==0){
    readpart("",ifile,UNLOAD,&error);
    updatetimes();
    return;
  }
  
  printf("Sizing particle data: %s\n",file);
  FORTgetsizes(file,&ibar,&jbar,&kbar,&nb,&nv,&nspr,&mxframepoints,&endian,&staticframe0,&error,lenfile);
  STRCPY(partsizefile,file);
  STRCAT(partsizefile,".sz");
  statfile=stat(file,&statbuffer);
  statfile2=stat(partsizefile,&statbuffer2);
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
    FORTgetdata1(&parttype,&error);
    FORTgetsizes2(&settmin_p,&tmin_p,&settmax_p,&tmax_p,
                     &nspr, &partframestep, &partpointstep, &npartpoints, &npartframes, &error);
    sizefile=fopen(partsizefile,"w");
    if(sizefile!=NULL){
      fprintf(sizefile,"%i %i %i %i %i %i %i %i %i",
          nb,nv,nspr,mxframepoints,staticframe0,npartpoints,npartframes,partframestep,partpointstep);
    }
    else{
      printf("*** warning:  unable to write to %s\n",partsizefile);
    }
    FORTgetsizes(file,&ibar,&jbar,&kbar,&nb,&nv,&nspr,&mxframepoints,&endian,&staticframe0,&error,lenfile);
  }
  npartpoints2=npartpoints;
  npartframes2=npartframes;
  if(npartframes>mxframes||npartpoints>mxpoints){
    if(npartframes>mxframes)npartframes=mxframes;
    if(npartpoints>mxpoints)npartpoints=mxpoints;
  }
  iframebeg=0;
  if(staticframe0==1)iframebeg=1;
  if(error!=0){
    printf("*** warning: problem reading %s\n",file);
    return;
  }
  if(npartpoints<=0){
    printf("*** warning: the particle file:%s is empty\n",file);
    return;
  }
  offsetmax=5;
  if(offsetmax>ibar/4)offsetmax=ibar/4;
  if(offsetmax>jbar/4)offsetmax=jbar/4;
  if(offsetmax>kbar/4)offsetmax=kbar/4;
  if(nspr>0){
    if(tspr==NULL){
      return_code=NewMemory((void **)&tspr,sizeof(float)*nspr);
    }
     else{
      return_code=ResizeMemory((void **)&tspr,sizeof(float)*nspr);
    }
     if(return_code==0){
      *errorcode=1;
      FORTclosepart();
      readpart("",ifile,UNLOAD,&error);
      return;
    }
  }

  if(NewMemory((void **)&parti->ptimes,sizeof(float)*mxframes)==0){
    *errorcode=1;
    FORTclosepart();
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  if(NewMemory((void **)&parti->xparts,npartpoints*sizeof(short))==0||
     NewMemory((void **)&parti->yparts,npartpoints*sizeof(short))==0||
     NewMemory((void **)&parti->zparts,npartpoints*sizeof(short))==0){
    *errorcode=1;
    FORTclosepart();
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  bytesperpoint=7;

  if(NewMemory((void **)&parti->tpart,npartpoints*sizeof(float))==0||
     NewMemory((void **)&parti->itpart,npartpoints*sizeof(unsigned char))==0||
     NewMemory((void **)&parti->isprink,npartpoints*sizeof(unsigned char))==0||
     NewMemory((void **)&parti->bframe,mxframes*sizeof(int))==0||
     NewMemory((void **)&parti->sframe,mxframes*sizeof(int))==0||
     NewMemory((void **)&parti->sprframe,mxframes*sizeof(int))==0){
      *errorcode=1;
      FORTclosepart();
      readpart("",ifile,UNLOAD,&error);
      return;
  }
  for(i=0;i<npartpoints;i++){
    parti->isprink[i]=0;
  }
  for(i=0;i<mxframes;i++){
    parti->sprframe[i]=0;
  }


  printf("Loading particle data: %s\n",file);
  FORTgetdata1(&parttype,&error);
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
    printf("*** warning: problem reading %s\n",file);
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  xbox=xbar0+xbar*xyzmaxdiff;
  ybox=ybar0+ybar*xyzmaxdiff;
  zbox=zbar0+zbar*xyzmaxdiff;
  offset_x=meshi->offset[0];
  offset_y=meshi->offset[1];
  offset_z=meshi->offset[2];
  FORTgetdata2(
    parti->xparts,parti->yparts,parti->zparts,
    parti->tpart,&parti->droplet_type,parti->isprink,
    tspr,parti->bframe,parti->sframe,parti->sprframe,parti->ptimes,&nspr,&npartpoints,&mxframes,&parti->nframes,
    &settmin_p,&settmax_p,&tmin_p,&tmax_p,&partframestep,&partpointstep, 
    &xbar0, &xbox, &ybar0, &ybox, &zbar0, &zbox,
    &offset_x, &offset_y, &offset_z,
    &error,1);
  if(npartframes2>mxframes||npartpoints2>mxpoints){
    if(npartframes2>mxframes){
      printf("*** warning number of frames (%i) in particle file is greater than %i\n",npartframes2,mxframes);
      printf("use: smokeview -frames %i casename\n",npartframes2);
      printf("to view all particle frames\n\n");
    }
    if(npartpoints2>mxpoints){
      printf("*** warning number of particles (%i) in particle \n    file is greater than %i\n",npartpoints2,mxpoints);
      printf("use: smokeview -points %i casename\n",npartpoints2);
      printf("to view all particles\n\n");
    }
  }
  if(error!=0||parti->nframes==0){
    if(error!=0)printf("*** warning: problem reading %s\n",file);
    *errorcode=1;
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  
  nmax = parti->bframe[parti->nframes-1]+parti->sframe[parti->nframes-1];
  printf("loaded: points=%i, size=%i KBytes, frames=%i\n",nmax,nmax*bytesperpoint/1024,parti->nframes);

  if(parttype==-1||parttype==-3){
    parti->particle_type=1;  /*  only color temperature */
  }
  else{
    parti->particle_type=0;
  }
  havesprinkpart=0;
  skip=0;
  if(staticframe0==1)skip=parti->sframe[0];
  tcopy=parti->tpart+skip;
  tmin=1000000000.0;
  tmax=-tmin;
  isprinkcopy=parti->isprink;
  for(n=skip;n<nmax;n++){
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
    if(*tcopy<tmin)tmin=*tcopy;
    if(*tcopy>tmax)tmax=*tcopy;
    if(*isprinkcopy==1){
      havesprinkpart=1;
    }
    tcopy++; 
    isprinkcopy++;
  }
  /* convert particle temperatures into integers pointing to an rgb color table */

  printf("computing particle color levels \n");
  if(parti->particle_type!=0||parti->droplet_type!=0){
    adjustpartbounds(parti->tpart,parti->particle_type,parti->droplet_type,parti->isprink,
      skip,nmax,setpartmin,&tmin,setpartmax,&tmax);
  }
  if(setpartmin == SET_MIN){
    tmin = partmin;
  }
  if(setpartmax == SET_MAX){
    tmax = partmax;
  }
  partmin=tmin;
  partmax=tmax;
  if(NewMemory((void **)&colorlabelpart,MAXRGB*sizeof(char *))==0){
    FORTclosepart();
    readpart("",ifile,UNLOAD,&error);
    *errorcode=1;
    return;
  }
  for(n=0;n<MAXRGB;n++){colorlabelpart[n]=NULL;}
  for(n=0;n<nrgb;n++){
    if(NewMemory((void **)&colorlabelpart[n],11)==0){
      *errorcode=1;
      FORTclosepart();
      readpart("",ifile,UNLOAD,&error);
      return;
    }
  }
  getPartColors(parti->tpart, skip, nmax, 
    parti->itpart,parti->isprink,parti->particle_type,parti->droplet_type,
    &tmin, &tmax, nrgb, colorlabelpart, partscale,partlevels256);
  for(n=0;n<skip;n++){
    parti->itpart[n]=0;
  }
  FREEMEMORY(parti->tpart);
  FREEMEMORY(parti->isprink);
  if(ResizeMemory((void **)&parti->xparts,nmax*sizeof(short))==0||
     ResizeMemory((void **)&parti->yparts,nmax*sizeof(short))==0||
     ResizeMemory((void **)&parti->zparts,nmax*sizeof(short))==0){
    FORTclosepart();
    *errorcode=1;
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  if(ResizeMemory((void **)&parti->itpart,nmax*sizeof(unsigned char))==0){
    FORTclosepart();
    *errorcode=1;
    readpart("",ifile,UNLOAD,&error);
    return;
  }
  updateglui();

#ifdef _DEBUG
  printf("After particle file load: ");
  PrintMemoryInfo;
#endif
  if(partinfo[ifile].evac==0){
    ReadPartFile=1;
  }
  else{
    ReadEvacFile=1;
  }
  if(partinfo[ifile].evac==0){
    visSmoke=1;
  }
  else{
    visEvac=1;
  }
  partinfo[ifile].loaded=1;
  partinfo[ifile].display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatetimes();
  updatemenu=1;
  IDLE();

  glutPostRedisplay();
}

/* ------------------ drawEvac ------------------------ */

void drawEvac(const particle *parti){
  drawPart(parti);
}


/* ------------------ drawPart5 ------------------------ */

void drawPart5(const particle *parti){
  int ipframe;
  part5data *datacopy,*datapast;
  int nclasses;
  int i,j;

  if(current_property==NULL)return;
  ipframe=parti->iframe;
  nclasses = parti->nclasses;
  datacopy = parti->data5+nclasses*ipframe;
  CheckMemory;
  if(part5show==1){
    if(streak5show==0||(streak5show==1&&showstreakhead==1)){
      glPointSize(partpointsize);
      glBegin(GL_POINTS);
      for(i=0;i<parti->nclasses;i++){
        short *sx, *sy, *sz;
        unsigned char *vis, *color;
        part5class *partclassi;
        int partclass_index, itype, vistype, class_vis;

        partclassi = parti->partclassptr[i];
        partclass_index = partclassi - partclassinfo;

        vistype=current_property->class_present[partclass_index];
        class_vis=current_property->class_vis[partclass_index];


        if(vistype==0||datacopy->npoints<=0||(vistype==1&&class_vis==0)){
          datacopy++;
          continue;
        }
        itype = current_property->class_types[partclass_index];

        sx = datacopy->sx;
        sy = datacopy->sy;
        sz = datacopy->sz;
        vis = datacopy->vis_part;

        if(itype==-1){
          glColor4fv(datacopy->partclassbase->rgb);
          for(j=0;j<datacopy->npoints;j++){
            if(vis[j]==1)glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
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

        datacopy++;
      }
      glEnd();
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

    part5class *partclassi;
    int partclass_index, itype, vistype, class_vis;

    partclassi = parti->partclassptr[i];
    partclass_index = partclassi - partclassinfo;

    vistype=current_property->class_present[partclass_index];
    class_vis=current_property->class_vis[partclass_index];

    if(vistype==0||datacopy->npoints<=0||(vistype==1&&class_vis==0)){
      datacopy++;
      continue;
    }
    itype = current_property->class_types[partclass_index];

    sx = datacopy->sx;
    sy = datacopy->sy;
    sz = datacopy->sz;
    vis = datacopy->vis_part;

    if(itype==-1){

      // draw the streak line

      glColor4fv(datacopy->partclassbase->rgb);
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
          jj = get_tagindex(datapast,tagval);
          if(jj<0)break;
          sxx = datapast->sx;
          syy = datapast->sy;
          szz = datapast->sz;
          glVertex3f(xplts[sxx[jj]],yplts[syy[jj]],zplts[szz[jj]]);
        }
        glEnd();
      }

      // draw the dot at the end of the streak line

      if(showstreakhead==1){
        sx = datacopy->sx;
        sy = datacopy->sy;
        sz = datacopy->sz;
        vis = datacopy->vis_part;
        glPointSize(6.0);
        glColor4fv(datacopy->partclassbase->rgb);
        glBegin(GL_POINTS);
        for(j=0;j<datacopy->npoints;j++){
          if(vis[j]==0)continue;
          glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
        }
        glEnd();
      }
    }
    else{
      unsigned char *color;

      // draw the streak line

      color=datacopy->irvals+itype*datacopy->npoints;

      for(j=0;j<datacopy->npoints;j++){
        int tagval;
        tagval=datacopy->tags[j];
        if(vis[j]==0)continue;
        
        glBegin(GL_LINE_STRIP);
        glColor4fv(rgb_full[color[j]]);
        glVertex3f(xplts[sx[j]],yplts[sy[j]],zplts[sz[j]]);
        for(k=1;k<streak5step;k++){
          int jj;

          if(ipframe-k<0)break;
          datapast = parti->data5+nclasses*(ipframe-k)+i;
          jj = get_tagindex(datapast,tagval);
          if(jj<0||datapast->irvals==NULL)break;
          sxx = datapast->sx;
          syy = datapast->sy;
          szz = datapast->sz;
          color=datapast->irvals+itype*datapast->npoints;

          glColor4fv(rgb_full[color[jj]]);
          glVertex3f(xplts[sxx[jj]],yplts[syy[jj]],zplts[szz[jj]]);
        }
        glEnd();
      }
    }

    datacopy++;
  }
  }

}

/* ------------------ drawPart ------------------------ */

void drawPart(const particle *parti){
  short *xpoints, *ypoints, *zpoints;
  unsigned char *itpoint=NULL;


  int n;
  int nsmokepoints, nsprpoints;
  int ipframe;
  int droplet_type, particle_type;
  float *rgb_smoke, *rgb_ismoke;

  ipframe=parti->iframe;
  rgb_smoke = rgb_part;

  if(parti->ptimes[0]>times[itime])return;

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
        for (n = 0; n < nsmokepoints; n++) {
          glColor4fv(rgb[itpoint[n]]);
          glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
        }
      }
      else{
        for (n = 0; n < nsmokepoints; n++) {
      //  glColor4fv(rgb_full[itpoint[n]]);
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
        for (n = nsmokepoints; n < nsmokepoints+nsprpoints; n++) {
          glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
        }
      }
      else{
        for (n = nsmokepoints; n < nsmokepoints+nsprpoints; n++) {
          glColor4fv(rgb_full[itpoint[n]]);
          glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
        }
      }
    }
  }

  glEnd();



}

/* ------------------ drawStaticPart ------------------------ */

void drawStaticPart(const particle *parti){

  short *xpoints, *ypoints, *zpoints;

  int n;
  int nsmokepoints, nsprpoints;
  int ipframe;
  
  /* define the data locations to look at */

  ipframe=0;
  if(parti->version==0){
    xpoints = parti->xparts + parti->bframe[ipframe];
    ypoints = parti->yparts + parti->bframe[ipframe];
    zpoints = parti->zparts + parti->bframe[ipframe];
  }
  /* isprinkframe = isprink + bframe[ipframe];*/

  nsprpoints = parti->sprframe[ipframe];
  nsmokepoints = parti->sframe[ipframe]-nsprpoints;

  glPointSize(partpointsize);

  glColor4fv(static_color);
  glBegin(GL_POINTS);
  if(parti->version==0){
    for(n=0;n<nsmokepoints+nsprpoints;n++){
      glVertex3f(xplts[xpoints[n]],yplts[ypoints[n]],zplts[zpoints[n]]);
    }
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
  particle *parti, *partj;
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
  particle *parti;
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
        sprintf(label,"Mesh %i",1+parti->blocknumber);
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

/* ------------------ update_visSmokePart ------------------------ */

void update_visSmokePart(void){
//  int smoke_all=1;
//  int smoke_some=0;
//  int i;
//  particle *parti;

  /*
  for(i=0;i<npart;i++){
    parti = partinfo + i;
    if(parti->loaded==0||parti->evac==1)continue;
    if(parti->display_smoke==0)smoke_all=0;
    if(parti->display_smoke==1)smoke_some=1;
  }
  if(smoke_all==1){
    visSmokePart=2;
  }
  else{
    if(smoke_some==1){
      visSmokePart=1;
    }
    else{
      visSmokePart=0;
    }
  }
  */
}
