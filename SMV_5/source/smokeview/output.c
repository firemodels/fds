// $Date: 2007-10-07 22:08:47 -0400 (Sun, 07 Oct 2007) $ 
// $Revision: 800 $
// $Author: gforney $

#include "options.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <stdio.h>
#include <string.h>
#include "flowfiles.h"
#include "smokeviewdefs.h"
#include "MALLOC.h"
#include "smokeviewvars.h"

// svn revision character string
char output_revision[]="$Revision: 800 $";

/* ------------------ outputAxisLabels ------------------------ */

void outputAxisLabels(){
  float x, y, z;
  float x0, y0;
  char XX[1]={'X'}, YY[1]={'Y'}, ZZ[1]={'Z'};
  int ibar,jbar,kbar;
  float *xplt,*yplt,*zplt;
  int i;
  mesh *meshi;

  glColor3fv(foregroundcolor);
  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    ibar=meshi->ibar;
    jbar=meshi->jbar;
    kbar=meshi->kbar;
    xplt=meshi->xplt;
    yplt=meshi->yplt;
    zplt=meshi->zplt;


    x0 = xplt[0]-0.02;
    y0 = yplt[0]-0.02;
    x = (xplt[0]+xplt[ibar])/2.0;
    y = (yplt[0]+yplt[jbar])/2.0;
    z = (zplt[0]+zplt[kbar])/2.0;
    glRasterPos3f(x,y0,(float)0.0);
    glutBitmapCharacter(large_font,XX[0]);
    glRasterPos3f(x0,y,(float)0.0);
    glutBitmapCharacter(large_font,YY[0]);
    glRasterPos3f(x0,y0,z);
    glutBitmapCharacter(large_font,ZZ[0]);
  }




}

/* ------------------ output3Text ------------------------ */

void output3Text(float *color, float x, float y, float z, const char *string){
  int len, i;

  if(string==NULL)return;
  glColor3fv(color);
  glRasterPos3f(x, y, z);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)glutBitmapCharacter(large_font,string[i]);
}

/* ------------------ outputLargeText ------------------------ */

void outputLargeText(float x, float y, const char *string){
  int len, i;
  if(string==NULL)return;
  glColor3fv(foregroundcolor);
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,string[i]);
}


/* ------------------ outputText ------------------------ */

void outputText(float x, float y, const char *string){
  int len, i;

  if(string==NULL)return;
  glColor3fv(foregroundcolor);
  glRasterPos2f(x, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++)glutBitmapCharacter(large_font,string[i]);
}

/* ------------------ outputBarText ------------------------ */

void outputBarText(float x, float y, const GLfloat *color, const char *string){
  int len, i;
  int length;
  float xlength;

  if(string==NULL)return;

  length=glutBitmapLength(small_font, (const unsigned char *)string); 
  xlength = length*barright/dwinWW+0.02;

  glColor3fv(color);
  glRasterPos2f(x-xlength, y);
  len = (int) strlen(string);
  for (i = 0; i < len; i++){
    glutBitmapCharacter(small_font,string[i]);
  }
}

/* ------------------ bench_out ------------------------ */

void bench_out(float localframerate){
  FILE *fileout;
  int i;
  particle *parti;
  slice *slicei,*u,*v,*w,*val;
  vslice *vd;
  smoke3d *smoke3di;
  patch *patchi;
  iso *isoi;
#ifdef pp_memstatus
  unsigned int memused, tmem;
#endif
  GLint nred, ngreen, nblue, nalpha, ndepth;

#ifdef pp_memstatus
  MEMSTATUS(0,NULL,&memused,&tmem);
#endif

  printf("*** benchmarking completed\n\n");
  fileout=fopen("svbenchmark.txt","w");
  if(fileout==NULL)fileout=stdout;

  fprintf(fileout,"If you would like to report the benchmark results recorded\n"); 
  fprintf(fileout,"below, please email them to glenn.forney@nist.gov .\n");
  fprintf(fileout,"Use \"Smokeview Benchmark\" for your subject line.\n\n");

  fprintf(fileout,"Fill in the blanks below as best you can.\n\n");
  
  fprintf(fileout," %s\n",TITLE);
  fprintf(fileout,"Average Frame Rate:%3.1f\n",localframerate);
  fprintf(fileout,"    CPU Type/Speed:\n");
  fprintf(fileout,"        Video Card:\n");
  glGetIntegerv(GL_RED_BITS,&nred);    
  glGetIntegerv(GL_GREEN_BITS,&ngreen);
  glGetIntegerv(GL_BLUE_BITS,&nblue); 
  glGetIntegerv(GL_DEPTH_BITS,&ndepth);
  glGetIntegerv(GL_ALPHA_BITS,&nalpha);
  fprintf(fileout,"          Red bits:%i\n",nred);
  fprintf(fileout,"        Green bits:%i\n",ngreen);
  fprintf(fileout,"         Blue bits:%i\n",nblue);
  fprintf(fileout,"        Alpha bits:%i\n",nalpha);
  fprintf(fileout,"        Depth bits:%i\n",ndepth);
  fprintf(fileout,"      Window width:%i\n",screenWidth);
  fprintf(fileout,"     Window height:%i\n",screenHeight);

#ifdef pp_memstatus
  fprintf(fileout,"       Memory Used:%iMB\n",(int)memused);
  fprintf(fileout,"      Total Memory:%iMB\n",(int)tmem);
#else
  fprintf(fileout,"       Memory Used:\n");
  fprintf(fileout,"      Total Memory:\n");
#endif
  if(fds_filein==NULL){
    fprintf(fileout,"          FDS Case:\n");
  }
  else{
    fprintf(fileout,"          FDS Case:%s\n",fds_filein);
  }
  fprintf(fileout,"   Files displayed:\n");

  for(i=0;i<npartinfo;i++){
    parti=partinfo+i;
    if(parti->loaded==0||parti->display==0)continue;
    fprintf(fileout,"      %s\n",parti->file);
  }
  for(i=0;i<nslice;i++){
    slicei=sliceinfo+i;
    slicei->benchvis=0;
    if(slicei->loaded==1&&slicei->display==1)slicei->benchvis=1;
  }
  for(i=0;i<nvslice;i++){
    vd = vsliceinfo + i;
    if(vd->loaded==0||vd->display==0||sliceinfo[vd->ival].type!=islicetype)continue;
    if(vd->val==NULL)continue;
    u = vd->u;
    v = vd->v;
    w = vd->w;
    val = vd->val;
    if(u!=NULL)u->benchvis=1;
    if(v!=NULL)v->benchvis=1;
    if(w!=NULL)w->benchvis=1;
    if(val!=NULL)val->benchvis=1;
  }
  for(i=0;i<nslice;i++){
    slicei=sliceinfo+i;
    if(slicei->benchvis==0)continue;
    fprintf(fileout,"      %s\n",slicei->file);
  }
  for(i=0;i<nsmoke3d;i++){
    smoke3di=smoke3dinfo+i;
    if(smoke3di->loaded==0||smoke3di->display==0)continue;
    fprintf(fileout,"      %s\n",smoke3di->file);
  }
  for(i=0;i<npatch_files;i++){
    patchi=patchinfo+i;
    if(patchi->loaded==0||patchi->display==0)continue;
    fprintf(fileout,"      %s\n",patchi->file);
  }
  for(i=0;i<niso;i++){
    isoi=isoinfo+i;
    if(isoi->loaded==0||isoi->display==0)continue;
    fprintf(fileout,"      %s\n",isoi->file);
  }

  benchmark_flag=0;
  if(fileout!=stdout)fclose(fileout);
}

/* ------------------ drawLabels ------------------------ */

void drawLabels(void){
  int i;

  for(i=0;i<nlabels;i++){
    labeldata *labelcopy;

    labelcopy=labelinfo+i;
    {
      float *labelcolor,*tstart_stop,*xyz;
      int drawlabel=0;

      tstart_stop=labelcopy->tstart_stop;
      xyz=labelcopy->xyz;
      if(labelcopy->useforegroundcolor==1){
        labelcolor=foregroundcolor;
      }
      else{
        labelcolor=labelcopy->rgb;
      }
      if(plotstate!=DYNAMIC_PLOTS||showtime==0)drawlabel=1;
      if(drawlabel==0&&plotstate==DYNAMIC_PLOTS&&showtime==1){
        if(tstart_stop[0]<0.0||tstart_stop[1]<0.0)drawlabel=1;
        if(drawlabel==0&&times[itime]>=tstart_stop[0]&&times[itime]<=tstart_stop[1])drawlabel=1;
      }
      if(drawlabel==1){
        {
          float xpos,ypos,zpos;
          xpos = (xyz[0]-xbar0)/xyzmaxdiff;
          ypos = (xyz[1]-ybar0)/xyzmaxdiff;
          zpos = (xyz[2]-zbar0)/xyzmaxdiff;
         output3Text(labelcolor,xpos,ypos,zpos,labelcopy->label);
        }
      }
    }
  }
}

