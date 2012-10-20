// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char output_revision[]="$Revision$";

#include "options.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "string_util.h"
#include "smokeviewvars.h"

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
  for(i=0;i<nmeshes;i++){
    meshi=meshinfo+i;
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


/* ------------------ outputSText3 ------------------------ */

void outputSText3(float x, float y, float z, char *string){ 
  char *c;
  float u[3]={0.0,0.0,1.0},v[3];
  float axis[3],angle,scale=0.001,theta;
  float quateye[4],quatz[4],quat2[4],rot[16];
  float scale_x, scale_y;


  if(string==NULL)return;
  scale_x = xyzmaxdiff*((float)scaled_font3d_size/(float)104.76)/(float)port_pixel_width;
  scale_y = xyzmaxdiff*((float)scaled_font3d_size/(float)152.38)/(float)port_pixel_height;
  glPushMatrix();
  glTranslatef(x,y,z);
  v[0]=world_eyepos[0]-x;
  v[1]=world_eyepos[1]-y;
  v[2]=world_eyepos[2]-z;
  rotateu2v(u,v,axis,&angle);
  theta=atan2(v[0],-v[1])*RAD2DEG;
  angleaxis2quat(theta*DEG2RAD,u,quatz);
  angleaxis2quat(angle,axis,quateye);
  mult_quat(quateye,quatz,quateye);
  quat2rot(quateye,rot);
  /*
  glLoadIdentity();
  glMultMatrixf(rot);*/
  glRotatef(90.0,cos(theta*DEG2RAD),sin(theta*DEG2RAD),0.0);
  glRotatef(theta,u[0],u[1],u[2]);
  glScalef(scale_x,scale_y,1.0);
  for (c=string; *c != '\0'; c++){
    glutStrokeCharacter(GLUT_STROKE_ROMAN,*c);
  }
  glPopMatrix();
}


/* ------------------ outputSText2r ------------------------ */

void outputSText2r(float x, float y, float z, char *string){ 
  char *c;
  int total_width=0;
  float scale_x, scale_y;

  if(string==NULL)return;
  total_width=0;
  for (c=string; *c != '\0'; c++){
    total_width+=glutStrokeWidth(GLUT_STROKE_ROMAN,*c);
  }
  glPushMatrix();
  scale_x = port_unit_width*((float)scaled_font2d_size/(float)104.76)/(float)port_pixel_width;
  scale_y = port_unit_height*((float)scaled_font2d_size/(float)152.38)/(float)port_pixel_height;
  glTranslatef(x-scale_x*total_width,y,z);
  glScalef(scale_x,scale_y,1.0);
  for (c=string; *c != '\0'; c++){
    glutStrokeCharacter(GLUT_STROKE_ROMAN,*c);
  }
  glPopMatrix();
}

/* ------------------ outputSText2 ------------------------ */

void outputSText2(float x, float y, float z, char *string){ 
  char *c;
  int total_width=0;
  float scale_x, scale_y;

  if(string==NULL)return;
  total_width=0;
  for (c=string; *c != '\0'; c++){
    total_width+=glutStrokeWidth(GLUT_STROKE_ROMAN,*c);
  }
  glPushMatrix();
  scale_x = port_unit_width*((float)scaled_font2d_size/(float)104.76)/(float)port_pixel_width;
  scale_y = port_unit_height*((float)scaled_font2d_size/(float)152.38)/(float)port_pixel_height;
  glTranslatef(x,y,z);
  glScalef(scale_x,scale_y,1.0);
  for (c=string; *c != '\0'; c++){
    glutStrokeCharacter(GLUT_STROKE_ROMAN,*c);
  }
  glPopMatrix();
}

/* ------------------ output3Val ------------------------ */

void output3Val(float x, float y, float z, float val){
  char *c,string[256];

  sprintf(string,"%f",val);
  trimzeros(string);
  output3Text(foregroundcolor,x,y,z,string);
}

/* ------------------ output3Text ------------------------ */

void output3Text(float *color, float x, float y, float z, char *string){
  char *c;

  if(string==NULL)return;
  glColor3fv(color);

  if(fontindex==SCALED_FONT){
    outputSText3(x,y,z,string);
  }
  else{
    glRasterPos3f(x, y, z);
    for (c=string; *c!='\0'; c++){
      glutBitmapCharacter(large_font,*c);
    }
  }
}

/* ------------------ outputLargeText ------------------------ */

void outputLargeText(float x, float y, char *string){
  char *c;

  if(string==NULL)return;
  glColor3fv(foregroundcolor);
  glRasterPos2f(x, y);
  for (c=string; *c!='\0'; c++){
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*c);
  }
}

/* ------------------ outputText ------------------------ */

void outputText(float x, float y, char *string){
  char *c;

  if(string==NULL)return;
  glColor3fv(foregroundcolor);
  if(fontindex==SCALED_FONT){
    outputSText2(x,y,0.0,string);
    return;
  }
  else{
    glRasterPos2f(x, y);
    for (c=string; *c!='\0'; c++){
      glutBitmapCharacter(large_font,*c);
    }  
  }
}

/* ------------------ outputBarText ------------------------ */

void outputBarText(float x, float y, const GLfloat *color, char *string){
  int length;
  float xlength;
  float xbeg;
  char *c;

  if(string==NULL)return;
  glColor3fv(color);

  if(fontindex==SCALED_FONT){
    outputSText2r(x,y,0.0,string);
  }
  else{
    length=glutBitmapLength(small_font, (const unsigned char *)string); 
    xlength = length*barright/dwinWW+0.02;

    xbeg=x-xlength;
    if(xbeg<0.0)xbeg=0.0;
    glRasterPos2f(xbeg, y);
    for (c=string; *c!='\0'; c++){
      glutBitmapCharacter(small_font,*c);
    }
  }
}

/* ------------------ bench_out ------------------------ */

void bench_out(float localframerate){
  FILE *fileout;
  int i;
  partdata *parti;
  slicedata *slicei,*u,*v,*w,*val;
  vslicedata *vd;
  smoke3ddata *smoke3di;
  patchdata *patchi;
  isodata *isoi;
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
  for(i=0;i<nsliceinfo;i++){
    slicei=sliceinfo+i;
    slicei->benchvis=0;
    if(slicei->loaded==1&&slicei->display==1)slicei->benchvis=1;
  }
  for(i=0;i<nvsliceinfo;i++){
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
  for(i=0;i<nsliceinfo;i++){
    slicei=sliceinfo+i;
    if(slicei->benchvis==0)continue;
    fprintf(fileout,"      %s\n",slicei->file);
  }
  for(i=0;i<nsmoke3dinfo;i++){
    smoke3di=smoke3dinfo+i;
    if(smoke3di->loaded==0||smoke3di->display==0)continue;
    fprintf(fileout,"      %s\n",smoke3di->file);
  }
  for(i=0;i<npatchinfo;i++){
    patchi=patchinfo+i;
    if(patchi->loaded==0||patchi->display==0)continue;
    fprintf(fileout,"      %s\n",patchi->file);
  }
  for(i=0;i<nisoinfo;i++){
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
  labeldata *thislabel;

  for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
    float *labelcolor,*tstart_stop,*xyz;
    int drawlabel=0;

    tstart_stop=thislabel->tstart_stop;
    xyz=thislabel->xyz;
    if(thislabel->useforegroundcolor==1){
      labelcolor=foregroundcolor;
    }
    else{
      labelcolor=thislabel->rgb;
    }
    if(plotstate!=DYNAMIC_PLOTS||showtime==0)drawlabel=1;
    if(drawlabel==0&&plotstate==DYNAMIC_PLOTS&&showtime==1){
      if(tstart_stop[0]<0.0||tstart_stop[1]<0.0)drawlabel=1;
      if(drawlabel==0&&global_times[itimes]>=tstart_stop[0]-0.05&&global_times[itimes]<=tstart_stop[1]+0.05)drawlabel=1;
    }
    if(drawlabel==1){
      float xyz_pos[3];

      normalize_xyz(xyz_pos,xyz);
      output3Text(labelcolor,xyz_pos[0],xyz_pos[1],xyz_pos[2],thislabel->name);
    }
  }
}
