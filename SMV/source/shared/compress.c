// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ASSERT.h"
#include "MALLOC.h"
#include "compress.h"

// svn revision character string
char compress_revision[]="$Revision$";

#define MARK 255

/* ------------------ rle ------------------------ */

unsigned int rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out){
  unsigned char lastchar=MARK, cmark=MARK, thischar, *buffer_start;
  unsigned char *buffer_in_end;
  int nrepeats=1;

  buffer_start=buffer_out;
  buffer_in_end = buffer_in + nchars_in;

  while(buffer_in<buffer_in_end){
    thischar=*buffer_in;
    if(thischar==lastchar){
      nrepeats++;
    }
    else{
      nrepeats=1;
    }
    switch (nrepeats){
    case 1:
    case 2:
    case 3:
      *buffer_out=thischar;
      lastchar=thischar;
      break;
    case 4:
      buffer_out-=3;
      *buffer_out++=cmark;
      *buffer_out++=thischar;
    default:
      if(nrepeats!=4)buffer_out--;
      *buffer_out=nrepeats;
      if(nrepeats==254){
        nrepeats=1;
        lastchar=MARK;
      }
      break;
    }
    buffer_in++;
    buffer_out++;

  }
  return buffer_out-buffer_start;
}

/* ------------------ irle ------------------------ */

unsigned int irle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out){
  int nrepeats,nn;
  unsigned char thischar, *buffer_in_end;

  nn=0;
  buffer_in_end  = buffer_in  + nchars_in;

  while(buffer_in<buffer_in_end){
    if(*buffer_in==MARK){
      if(buffer_in+2>=buffer_in_end)break;
      buffer_in++;
      thischar=*buffer_in++;
      nrepeats=*buffer_in++;
      nn+=nrepeats;
      memset(buffer_out,thischar,nrepeats);
      buffer_out+=nrepeats;
    }
    else{
      *buffer_out++=*buffer_in++;
      nn++;
    }


  }
  return nn;
}

/* ------------------ compress_volsliceframe ------------------------ */

void compress_volsliceframe(float *data, int n_data, 
                float timeval, float *vmin_in, float *vmax_in,
                unsigned char **compressed_data, uLongf *ncompressed_data
                ){
  float valmin, valmax;
  int i;
  uLongf n_data_compressed2;
  unsigned char *c_data, *c_data_compressed1, *c_data_compressed2;
  int one=1;
  int version=0;
  int nbytes=1;

  // determine bounds

  if(vmin_in==NULL){
    valmin=data[0];
    for(i=1;i<n_data;i++){
      if(data[i]<valmin)valmin=data[i];
    }
  }
  else{
    valmin=*vmin_in;
  }

  if(vmax_in==NULL){
    valmax=data[0];
    for(i=1;i<n_data;i++){
      if(data[i]>valmax)valmax=data[i];
    }
  }
  else{
    valmax=*vmax_in;
  }

  // allocate buffers

  NewMemory((void **)&c_data,n_data);
  n_data_compressed2 = 1.01*(n_data+32) + 600;
  NewMemory((void **)&c_data_compressed1,n_data_compressed2);
  NewMemory((void **)&c_data_compressed2,n_data_compressed2);

  // scale data

  if(valmax>valmin){
    for(i=0;i<n_data;i++){
      float scaled_val;
      
      scaled_val=(data[i]-valmin)/(valmax-valmin);
      c_data[i]=255*scaled_val;
    }
  }
  else{
    for(i=0;i<n_data;i++){
      c_data[i]=0;
    }
  }

  //  compress data

  compress(c_data_compressed2,&n_data_compressed2, c_data, n_data);

  memcpy(c_data_compressed1+0,&one,4);
  memcpy(c_data_compressed1+4,&version,4);
  memcpy(c_data_compressed1+8,&n_data_compressed2,4);
  memcpy(c_data_compressed1+12,&nbytes,4);
  memcpy(c_data_compressed1+16,&n_data,4);
  memcpy(c_data_compressed1+20,&timeval,4);
  memcpy(c_data_compressed1+24,&valmin,4);
  memcpy(c_data_compressed1+28,&valmax,4);
  if(n_data_compressed2>0){
    memcpy(c_data_compressed1+32,c_data_compressed2,n_data_compressed2);
  }

  // resize and deallocate buffers
  ResizeMemory((void **)&c_data_compressed1, n_data_compressed2+32);
  FREEMEMORY(c_data);
  FREEMEMORY(c_data_compressed2);
  *compressed_data=c_data_compressed1;
  *ncompressed_data=n_data_compressed2;
}

/* ------------------ uncompress_volsliceframe ------------------------ */

int uncompress_volsliceframe(unsigned char *compressed_data,
                           float *data, int n_data, float *timeval,
                           unsigned char *fullbuffer
                ){
  float valmin, valmax;
  int i,ndatafile;
  unsigned char *c_data, *c_data_compressed1, *c_data_compressed2;
  uLongf countin,countout;
  
  valmin=*(float *)(compressed_data+24);
  valmax=*(float *)(compressed_data+28);
  *timeval=*(float *)(compressed_data+20);
  countin = *(int *)(compressed_data+8);
  ndatafile = *(int *)(compressed_data+16);

  uncompress(fullbuffer,&countout,compressed_data+32,countin);

  if(countout==ndatafile&&n_data>=countout){
    for(i=0;i<countout;i++){
      data[i]=valmin+fullbuffer[i]*(valmax-valmin)/255.0;
    }
  }
  return countout;
}
