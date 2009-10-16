// $Date: 2009-10-15 13:47:38 -0400 (Thu, 15 Oct 2009) $ 
// $Revision: 4936 $
// $Author: gforney $

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"

// svn revision character string
char IOdboundary_revision[]="$Revision: 4936 $";

/* ------------------ setup_slice ------------------------ */

void setup_boundary(FILE *stream_out){
  casedata *case1, *case2;
  int i;

  case1 = caseinfo;
  case2 = caseinfo + 1;

  for(i=0;i<case1->nboundary_files;i++){
    boundary *boundaryi;

    boundaryi = case1->boundaryinfo + i;
    boundaryi->boundary2 = getboundary(boundaryi,case2);
    if(boundaryi->boundary2!=NULL&&stream_out!=NULL){
      char outfile[1024];

      fprintf(stream_out,"%s\n",boundaryi->keyword);
      make_outfile(outfile,NULL,boundaryi->file,".bf");
      fprintf(stream_out,"%s\n",outfile);
      fprintf(stream_out,"%s\n",boundaryi->label.longlabel);
      fprintf(stream_out,"%s\n",boundaryi->label.shortlabel);
      fprintf(stream_out,"%s\n",boundaryi->label.unit);
    }
  }
}

/* ------------------ getboundary ------------------------ */

boundary *getboundary(boundary *boundaryin, casedata *case2){
  int i;
  float dx, dy, dz;

  dx = boundaryin->boundarymesh->dx/2.0;
  dy = boundaryin->boundarymesh->dy/2.0;
  dz = boundaryin->boundarymesh->dz/2.0;

  for(i=0;i<case2->nboundary_files;i++){
    boundary *boundaryout;

    boundaryout = case2->boundaryinfo + i;
    if(boundaryin->boundarytype!=boundaryout->boundarytype)continue;
    if(strcmp(boundaryin->label.longlabel,boundaryout->label.longlabel)!=0)continue;
    return boundaryout;
  }
  return NULL;
}

/* ------------------ diff_slices ------------------------ */

void diff_boundarys(void){
  int j;

  for(j=0;j<caseinfo->nboundary_files;j++){
    char *file1, *file2;
    char fullfile1[1024], fullfile2[1024], outfile[1024];
    boundary *boundaryi, *boundary1, *boundary2;
    FILE *stream;
    int unit1, unit2, unit3;
    FILE_SIZE len1,len2;
    //int error1=0,error2a=0;

    boundaryi = caseinfo->boundaryinfo+j;
    boundary1 = boundaryi;
    if(boundaryi->boundary2==NULL)continue;
    boundary2 = boundaryi->boundary2;
    file1 = boundaryi->file;
    file2 = boundaryi->boundary2->file;
    fullfile(fullfile1,sourcedir1,file1);
    fullfile(fullfile2,sourcedir2,file2);

    stream=fopen(fullfile1,"r");
    if(stream==NULL)continue;
    fclose(stream);

    stream=fopen(fullfile2,"r");
    if(stream==NULL)continue;
    fclose(stream);

    make_outfile(outfile,destdir,file1,".bf");
    if(strlen(outfile)==0)continue;
    stream=fopen(outfile,"w");
    if(stream==NULL)continue;
    fclose(stream);

    unit1=11;
    len1=strlen(fullfile1);
   // FORTopenboundary(fullfile1,&unit1,&caseinfo->endian,&error1,len1);

    unit2=12;
    len2=strlen(fullfile2);
   // FORTopenboundary(fullfile2,&unit2,&caseinfo->endian,&error2a,len2);


    FORTclosefortranfile(&unit1);
    FORTclosefortranfile(&unit2);
    FORTclosefortranfile(&unit3);
  }
}
