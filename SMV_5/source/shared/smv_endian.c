// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include "smv_endian.h"

// svn revision character string

char endian_revision[]="$Revision$";

/* ------------------ getendian ------------------------ */

int getendian(void){
  short val;
  char *cval;
  val=1;
  cval = (char *)&val+1;
  return (int)(*cval);
}

/* ------------------ endian_switch ------------------------ */

float float_switch(float val){ 
  float *val2ptr;
  unsigned char *buffer;
  unsigned char buffer2[4];

  buffer=(unsigned char *)&val;

  buffer2[0]=*(buffer+3);
  buffer2[1]=*(buffer+2);
  buffer2[2]=*(buffer+1);
  buffer2[3]=*buffer;

  val2ptr=(float *)&buffer2;
  return *val2ptr;

}

int int_switch(int val){ 
  int *val2ptr;
  unsigned char *buffer;
  unsigned char buffer2[4];

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
