// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $

#include "options.h"
int htonl(int token);
int ntohl(int token);
int getendian(void);
#define BIG_ENDIAN 1
#define LITTLE_ENDIAN 0
extern int endianswitch;
// svn revision character string

char endian_revision[]="$Revision: 12156 $";

/* ------------------ getendian ------------------------ */

int getendian(void){
  short val;
  char *cval;
  val=1;
  cval = (char *)&val+1;
  return (int)(*cval);
}

/* ------------------ htonl ------------------------ */

int htonl(int val){ // convert val from host to network format
  int *val2ptr;
  unsigned char *buffer;
  unsigned char buffer2[4];

  //  if the host is a BIG_ENDIAN computer then there is nothing to do
  //  (since network format == BIG_ENDIAN)

  if(getendian()==BIG_ENDIAN)return val;
  
  buffer=(unsigned char *)&val;

  buffer2[0]=*(buffer+3);
  buffer2[1]=*(buffer+2);
  buffer2[2]=*(buffer+1);
  buffer2[3]=*buffer;

  val2ptr=(int *)&buffer2;
  return *val2ptr;

}

/* ------------------ ntohl ------------------------ */

int ntohl(int val){ // convert val from network to host format
  int *val2ptr;
  unsigned char *buffer;
  unsigned char buffer2[4];


  //  if the host is a BIG_ENDIAN computer then there is nothing to do
  //  (since network format == BIG_ENDIAN)

  if(getendian()==BIG_ENDIAN)return val;
  
  buffer=(unsigned char *)&val;

  buffer2[0]=*(buffer+3);
  buffer2[1]=*(buffer+2);
  buffer2[2]=*(buffer+1);
  buffer2[3]=*buffer;

  val2ptr=(int *)&buffer2;
  return *val2ptr;

}

void endian_switch(void *val, int nval){
  unsigned char *ca, *cb, *cc, *cd;
  unsigned char c1, c2, c3, c4;
  int i;

  if(endianswitch==0)return;
  for(i=0;i<nval;i++){
    int ii;

    ii=4*i;
    ca=(unsigned char *)val+ii;
    cb=(unsigned char *)val+ii+1;
    cc=(unsigned char *)val+ii+2;
    cd=(unsigned char *)val+ii+3;
    c1=*ca;
    c2=*cb;
    c3=*cc;
    c4=*cd;
    *ca=c4;
    *cb=c3;
    *cc=c2;
    *cd=c1;
  }
}
