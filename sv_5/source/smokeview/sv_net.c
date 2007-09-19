#include "options.h"
#ifdef pp_SVNET
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "sv_net.h"
#include "MALLOC.h"
#include "smokeviewvars.h"

// svn revision character string
char sv_net_revision[]="$Revision: 614 $";

void remove_svcom(svcom *svc);
svcom *first_svcom, *last_svcom;
void send_command(char *command);

tcpdata *first_tcpdata, *last_tcpdata;
/*

typedef struct _svcom {
  unsigned int src_id, dest_id;
  unsigned char dest_ip[4];
  int command_type;
  char char_data[1024];
  float  rvals[16];
  int ivals[16];
  struct _svcom *prev, *next;
} svcom;

*/

#define SVC_KEYBOARD 0
#define SVC_KEYBOARD_UP 1
#define SVC_MOUSE 2
#define SVC_MOTION 3


/* ------------------ init_svcom ------------------------ */

void init_svcom(void){
  NewMemory((void **)&first_svcom,sizeof(svcom));
  NewMemory((void **)&last_svcom,sizeof(svcom));
  strcpy(first_svcom->label,"first");
  strcpy(last_svcom->label,"last");
  first_svcom->next=last_svcom;
  first_svcom->prev=NULL;
  last_svcom->prev=first_svcom;
  last_svcom->next=NULL;

  NewMemory((void **)&first_tcpdata,sizeof(tcpdata));
  NewMemory((void **)&last_tcpdata,sizeof(tcpdata));
  strcpy(first_tcpdata->label,"first");
  strcpy(last_tcpdata->label,"last");
  first_tcpdata->next=last_tcpdata;
  first_tcpdata->prev=NULL;
  last_tcpdata->prev=first_tcpdata;
  last_tcpdata->next=NULL;
  n_tcpdata = 0;
  tcp_test=1;
}

/* ------------------ new_svcom ------------------------ */

svcom *new_svcom(void){
  svcom *svc;
  svcom *prev, *next;

  NewMemory((void **)&svc,sizeof(svcom));
  strcpy(svc->label,"");
  next = last_svcom;
  prev = next->prev;

  prev->next=svc;
  next->prev=svc;

  svc->next=next;
  svc->prev=prev;
  return svc;
}


/* ------------------ remove_svcom ------------------------ */

void remove_svcom(svcom *svc){
  svcom *prev, *next;

  prev = svc->prev;
  next = svc->next;
  FREEMEMORY(svc);
  prev->next=next;
  next->prev=prev;
}

/* ------------------ new_tcpdata ------------------------ */

tcpdata *new_tcpdata(void){
  tcpdata *tcp;
  tcpdata *prev, *next;

  NewMemory((void **)&tcp,sizeof(tcpdata));
  strcpy(tcp->label,"");
  next = last_tcpdata;
  prev = next->prev;

  prev->next=tcp;
  next->prev=tcp;

  tcp->next=next;
  tcp->prev=prev;
  return tcp;
}

/* ------------------ remove_tcpdata ------------------------ */

void remove_tcpdata(tcpdata *tcp){
  tcpdata *prev, *next;

  prev = tcp->prev;
  next = tcp->next;
  FREEMEMORY(tcp);
  prev->next=next;
  next->prev=prev;
}

/* ------------------ put_keyboard ------------------------ */

void put_keyboard(unsigned char key, int x, int y, int modifier_state){
  svcom *svc;
  char command[1024];

  svc = new_svcom();
  svc->command_type=SVC_KEYBOARD;
  svc->char_data[0]=key;
  svc->ivals[0]=x;
  svc->ivals[1]=y;
  svc->ivals[2]=modifier_state;
  strcpy(svc->label,"keyboard");
  if(n_tcpdata>0){
    sprintf(command,"keyboard_down %s %i %i %i",key,x,y,modifier_state);
    send_command(command);
  }
}


/* ------------------ put_motion ------------------------ */

void put_motion(int xm, int ym){
  svcom *svc;
  char command[1024];

  svc = new_svcom();
  svc->command_type=SVC_MOTION;
  svc->ivals[0]=xm;
  svc->ivals[1]=ym;
  strcpy(svc->label,"motion");
  if(n_tcpdata>0){
    sprintf(command,"motion %i %i ",xm,ym);
    send_command(command);
  }
}

/* ------------------ put_mouse ------------------------ */

void put_mouse(int button, int state, int x, int y, int modifier_state){
  svcom *svc;
  char command[1024];

  svc = new_svcom();
  svc->command_type=SVC_MOUSE;
  svc->ivals[0]=button;
  svc->ivals[1]=state;
  svc->ivals[2]=x;
  svc->ivals[3]=y;
  svc->ivals[4]=modifier_state;
  strcpy(svc->label,"mouse");
  if(n_tcpdata>0){
    sprintf(command,"mouse %i %i %i %i %i",button,state,x,y,modifier_state);
    send_command(command);
  }
}

/* ------------------ put_keyboard_up ------------------------ */

void put_keyboard_up(unsigned char key, int x, int y){
  svcom *svc;
  char command[1024];

  svc = new_svcom();
  svc->command_type=SVC_KEYBOARD_UP;
  svc->char_data[0]=key;
  svc->ivals[0]=x;
  svc->ivals[1]=y;
  strcpy(svc->label,"keyboard_up");
  if(n_tcpdata>0){
    sprintf(command,"keyboard_up %s %i %i",key,x,y);
    send_command(command);
  }
}
/* ------------------ parse_keyboard ------------------------ */

void parse_svcom(void){
  svcom *svc,*svc_next;
  int x, y;
  unsigned char c;
  int button, state, modifier_state;
  int xm, ym;
  int didit=0;

  svc=first_svcom->next;
  for(svc=first_svcom->next;svc->next!=NULL;){
    didit=1;
    switch (svc->command_type){
      case SVC_KEYBOARD:
        x = svc->ivals[0];
        y = svc->ivals[1];
        modifier_state = svc->ivals[2];
        c = svc->char_data[0];
        svc_keyboard(c,x,y,modifier_state);
        svc_next=svc->next;
        remove_svcom(svc);
        svc=svc_next; 
        break;
      case SVC_KEYBOARD_UP:
        x = svc->ivals[0];
        y = svc->ivals[1];
        c = svc->char_data[0];
        svc_keyboard_up(c,x,y);
        svc_next=svc->next;
        remove_svcom(svc);
        svc=svc_next; 
        break;
      case SVC_MOUSE:
        button = svc->ivals[0];
        state = svc->ivals[1];
        x = svc->ivals[2];
        y = svc->ivals[3];
        modifier_state = svc->ivals[4];
        svc_mouse(button,state,x,y,modifier_state);
        svc_next=svc->next;
        remove_svcom(svc);
        svc=svc_next; 
        break;
      case SVC_MOTION:
        xm = svc->ivals[0];
        ym = svc->ivals[1];
        svc_motion(xm,ym);
        svc_next=svc->next;
        remove_svcom(svc);
        svc=svc_next; 
        break;
      default:
        svc=svc->next;
        break;
    }
  }
  if(didit==1)glutPostRedisplay();
}


/* ------------------ send_command ------------------------ */

void send_command(char *command){
  tcpdata *tcp;

  tcp=first_tcpdata->next;
  for(tcp=first_tcpdata->next;tcp->next!=NULL;){
  }
}

#endif
