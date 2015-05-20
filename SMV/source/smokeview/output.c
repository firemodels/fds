// $Date: 2014-07-24 13:03:47 -0400 (Thu, 24 Jul 2014) $ 
// $Revision: 20000 $
// $Author: gforney $

// svn revision character string
char output_revision[]="$Revision: 20000 $";

#include "options.h"
#include GLUT_H
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "smokeviewvars.h"

/* ------------------ outputAxisLabels ------------------------ */

void outputAxisLabels(){
  float x, y, z;
  float x0, y0, z0;

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);

  x = (xbar0+xbarORIG)/2.0;
  y = (ybar0+ybarORIG)/2.0;
  z = (zbar0+zbarORIG)/2.0;
  x0 = xbar0 - SCALE2FDS(0.02);
  y0 = ybar0 - SCALE2FDS(0.02);
  z0 = zbar0 - SCALE2FDS(0.02);

  output3Text(foregroundcolor,   x,y0, z0, "X");
  output3Text(foregroundcolor, x0,  y, z0, "Y");
  output3Text(foregroundcolor, x0,y0,   z, "Z");
  
  glPopMatrix();
}

/* ------------------ outputSText3 ------------------------ */

void outputSText3(float x, float y, float z, char *string){ 
  char *c;
  float u[3]={0.0,0.0,1.0},v[3];
  float axis[3],angle,theta;
  float quateye[4],quatz[4],rot[16];
  float scale_x, scale_y;


  if(string==NULL)return;
  scale_x = SCALE2FDS(scaled_font3d_height2width*(float)scaled_font3d_height/(float)104.76)/(float)port_pixel_width;
  scale_y = SCALE2FDS((float)scaled_font3d_height/(float)152.38)/(float)port_pixel_height;
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
  scale_x = port_unit_width*(scaled_font2d_height2width*(float)scaled_font2d_height/(float)104.76)/(float)port_pixel_width;
  scale_y = port_unit_height*((float)scaled_font2d_height/(float)152.38)/(float)port_pixel_height;
  if(renderdoublenow!=0){
    scale_x *= (float)nrender_rows;
    scale_y *= (float)nrender_rows;
  }
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
  scale_x = (25.0/36.0)*port_unit_width*(scaled_font2d_height2width*(float)scaled_font2d_height/(float)104.76)/(float)port_pixel_width;
  scale_y = (12.0/18.0)*(25.0/18.0)*port_unit_height*((float)scaled_font2d_height/(float)152.38)/(float)port_pixel_height;
  if(renderdoublenow!=0){
    scale_x *= (float)nrender_rows;
    scale_y *= (float)nrender_rows;
  }
  glTranslatef(x,y,z);
  glScalef(scale_x,scale_y,1.0);
  glTranslatef(0.0,25.0,0.0);
  for (c=string; *c != '\0'; c++){
    glutStrokeCharacter(GLUT_STROKE_ROMAN,*c);
  }
  glPopMatrix();
}

/* ------------------ output3Val ------------------------ */

void output3Val(float x, float y, float z, float val){
  char string[256];

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
    scale_3dfont();
    outputSText3(x,y,z,string);
  }
  else{
    glRasterPos3f(x, y, z);
    for (c=string; *c!='\0'; c++){
      glutBitmapCharacter(large_font,(unsigned char)*c);
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
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,(unsigned char)*c);
  }
}

/* ------------------ outputText ------------------------ */

void outputText(float x, float y, char *string){
  char *c;

  if(string==NULL)return;
  glColor3fv(foregroundcolor);
  if(fontindex==SCALED_FONT){
    scale_2dfont();
    outputSText2(x,y,0.0,string);
    return;
  }
  else{
    glRasterPos2f(x, y);
    for (c=string; *c!='\0'; c++){
      glutBitmapCharacter(large_font,(unsigned char)*c);
    }  
  }
}

/* ------------------ outputBarText ------------------------ */

void outputBarText(float x, float y, const GLfloat *color, char *string){
  char *c;

  if(string==NULL)return;
  glColor3fv(color);

  if(fontindex==SCALED_FONT){
    scale_2dfont();
    outputSText2(x,y,0.0,string);
  }
  else{
    glRasterPos2f(x, y);
    for (c=string; *c!='\0'; c++){
      glutBitmapCharacter(small_font,(unsigned char)(*c));
    }
  }
}

/* ------------------ drawLabels ------------------------ */

void drawLabels(void){
  labeldata *thislabel;

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);
  for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
    float *labelcolor,*tstart_stop,*xyz;
    int drawlabel;

    drawlabel=0;
    tstart_stop=thislabel->tstart_stop;
    xyz=thislabel->xyz;
    if(thislabel->useforegroundcolor==1){
      labelcolor=foregroundcolor;
    }
    else{
      labelcolor=thislabel->frgb;
    }
    if(plotstate!=DYNAMIC_PLOTS||thislabel->show_always==1||showtime==0)drawlabel=1;
    if(drawlabel==0&&plotstate==DYNAMIC_PLOTS&&showtime==1){
      if(tstart_stop[0]<0.0||tstart_stop[1]<0.0)drawlabel=1;
      if(drawlabel==0&&global_times[itimes]>=tstart_stop[0]-0.05&&global_times[itimes]<=tstart_stop[1]+0.05)drawlabel=1;
    }
    if(drawlabel==1){
      output3Text(labelcolor,xyz[0],xyz[1],xyz[2],thislabel->name);
    }
  }
  glPopMatrix();
}


/* ------------------ LABEL_Next ------------------------ */

labeldata *LABEL_Next(labeldata *label){
  labeldata *thislabel;

  if(label==NULL)return NULL;
  if(label_first_ptr->next->next==NULL)return NULL;
  for(thislabel=label->next;thislabel!=label;thislabel=thislabel->next){
    if(thislabel->next==NULL)thislabel=label_first_ptr->next;
    if(thislabel->labeltype==TYPE_SMV)continue;
    return thislabel;
  }
  return NULL;
}

/* ------------------ LABEL_Previous ------------------------ */

labeldata *LABEL_Previous(labeldata *label){
  labeldata *thislabel;

  if(label==NULL)return NULL;
  if(label_last_ptr->prev->prev==NULL)return NULL;
  for(thislabel=label->prev;thislabel!=label;thislabel=thislabel->prev){
    if(thislabel->prev==NULL)thislabel=label_last_ptr->prev;
    if(thislabel->labeltype==TYPE_SMV)continue;
    return thislabel;
  }
  return NULL;
}

/* ------------------ LABEL_Init ------------------------ */

int LABEL_Init(labeldata *gl){
  labeldata *thislabel;

  for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
    if(thislabel->labeltype==TYPE_SMV)continue;
    LABEL_copy(gl,thislabel);
    return 1;
  }
  return 0;
}

/* ------------------ LABEL_Get_Nuserlabels ------------------------ */

int LABEL_Get_Nuserlabels(void){
  int count=0;
  labeldata *thislabel;

  for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
    if(thislabel->labeltype==TYPE_INI)count++;
  }
  return count;
}

/* ------------------ LABEL_get ------------------------ */

labeldata *LABEL_get(char *name){
  labeldata *thislabel;

  if(name==NULL)return NULL;
  for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
    if(thislabel->name==NULL)return NULL;
    if(strcmp(thislabel->name,name)==0)return thislabel;
  }
  return NULL;
}

/* ------------------ LABEL_insert_before ------------------------ */

void LABEL_insert_before(labeldata *listlabel, labeldata *label){
  labeldata *prev, *next;

  next = listlabel;
  prev = listlabel->prev;
  prev->next = label;
  label->prev = prev;
  next->prev=label;
  label->next=next;
}

/* ------------------ LABEL_delete ------------------------ */

void LABEL_delete(labeldata *label){
  labeldata *prev, *next;

  prev = label->prev;
  next =label->next;
  CheckMemory;
  FREEMEMORY(label);
  prev->next=next;
  next->prev=prev;
}

/* ------------------ LABEL_copy ------------------------ */

void LABEL_copy(labeldata *label_to, labeldata *label_from){
  labeldata *prev, *next;

  prev=label_to->prev;
  next=label_to->next;
  memcpy(label_to,label_from,sizeof(labeldata));
  label_to->prev=prev;
  label_to->next=next;

}

/* ------------------ LABEL_resort ------------------------ */

void LABEL_resort(labeldata *label){
  labeldata labelcopy;

  CheckMemory;
  memcpy(&labelcopy,label,sizeof(labeldata));
  CheckMemory;
  LABEL_delete(label);
  LABEL_insert(&labelcopy);
}

/* ------------------ LABEL_insert_after ------------------------ */

void LABEL_insert_after(labeldata *listlabel, labeldata *label){
  labeldata *prev, *next;

  prev = listlabel;
  next = listlabel->next;
  prev->next = label;
  label->prev = prev;
  next->prev=label;
  label->next=next;
}

/* ------------------ LABEL_print ------------------------ */

void LABEL_print(void){
  labeldata *thislabel;
  float *xyz;

  for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
    xyz = thislabel->xyz;
    PRINTF("label: %s position: %f %f %f\n",thislabel->name,xyz[0],xyz[1],xyz[2]);
  }
}

/* ------------------ LABEL_insert ------------------------ */

labeldata *LABEL_insert(labeldata *labeltemp){
  labeldata *newlabel, *thislabel;
  labeldata *firstuserptr, *lastuserptr;

  NewMemory((void **)&newlabel,sizeof(labeldata));
  memcpy(newlabel,labeltemp,sizeof(labeldata));

  thislabel = LABEL_get(newlabel->name);
  if(thislabel!=NULL){
    LABEL_insert_after(thislabel->prev,newlabel);
    return newlabel;
  }

  firstuserptr=label_first_ptr->next;
  if(firstuserptr==label_last_ptr)firstuserptr=NULL;

  lastuserptr=label_last_ptr->prev;
  if(lastuserptr==label_first_ptr)lastuserptr=NULL;

  if(firstuserptr!=NULL&&strcmp(newlabel->name,firstuserptr->name)<0){
    LABEL_insert_before(firstuserptr,newlabel);
    return newlabel;
  }
  if(lastuserptr!=NULL&&strcmp(newlabel->name,lastuserptr->name)>0){
    LABEL_insert_after(lastuserptr,newlabel);
    return newlabel;
  }
  if(firstuserptr==NULL&&lastuserptr==NULL){
    LABEL_insert_after(label_first_ptr,newlabel);
    return newlabel;
  }
  for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
    labeldata *nextlabel;

    nextlabel=thislabel->next;
    if(strcmp(thislabel->name,newlabel->name)<0&&strcmp(newlabel->name,nextlabel->name)<0){
      LABEL_insert_after(thislabel,newlabel);
      return newlabel;
    }
  }
  return NULL;
}

/* ----------------------- scale_2dfont ----------------------------- */

void scale_2dfont(void){
  if(render_multi!=0){
    glLineWidth((float)nrender_rows*(float)scaled_font2d_thickness);
  }
  else{
    glLineWidth((float)scaled_font2d_thickness);
  }
}

/* ----------------------- scale_3dfont ----------------------------- */

void scale_3dfont(void){
  if(render_multi!=0){
    glLineWidth((float)nrender_rows*(float)scaled_font3d_thickness);
  }
  else{
    glLineWidth((float)scaled_font3d_thickness);
  }
}

