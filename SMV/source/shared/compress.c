#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "MALLOC.h"
#include "compress.h"

#define MARK 255

/* ------------------ compress_zlib ------------------------ */

int compress_zlib(unsigned char *dest, uLongf *destLen, unsigned char *source, int sourceLen){
  return compress(dest, destLen, source, sourceLen);
}

/* ------------------ uncompress_zlib ------------------------ */

int uncompress_zlib(unsigned char *dest, uLongf *destLen, unsigned char *source, int sourceLen){
  return uncompress(dest, destLen, source, sourceLen);
}

/* ------------------ compress_rle ------------------------ */

unsigned int compress_rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out){
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
    switch(nrepeats){
    case 1:
    case 2:
    case 3:
      *buffer_out=thischar;
      lastchar=thischar;
      break;
    default:
      if(nrepeats==4){
        buffer_out-=3;
        *buffer_out++=cmark;
        *buffer_out++=thischar;
      }
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

/* ------------------ uncompress_rle ------------------------ */

unsigned int uncompress_rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out){
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

void compress_volsliceframe(float *data_in, int n_data_in, float timeval_in, float *vmin_in, float *vmax_in,
                unsigned char **compressed_data_out, uLongf *ncompressed_data_out
                ){
  float valmin, valmax;
  int i;
  uLongf n_data_compressed,n_data_compressedm32;
  unsigned char *c_data, *c_data_compressed;
  int one=1;
  int version=0;
  int nbytes=1;

  // determine bounds
  CheckMemory;
  if(vmin_in==NULL){
    valmin=data_in[0];
    for(i=1;i<n_data_in;i++){
      if(data_in[i]<valmin)valmin=data_in[i];
    }
    CheckMemory;
  }
  else{
    valmin=*vmin_in;
  }
  CheckMemory;

  if(vmax_in==NULL){
    valmax=data_in[0];
    for(i=1;i<n_data_in;i++){
      if(data_in[i]>valmax)valmax=data_in[i];
    }
    CheckMemory;
  }
  else{
    valmax=*vmax_in;
    CheckMemory;
  }
  CheckMemory;

  // allocate buffers

  n_data_compressed = 1.1*(n_data_in+32) + 600;
  n_data_compressedm32=n_data_compressed-32;
  NewMemory((void **)&c_data,n_data_in);
  NewMemory((void **)&c_data_compressed,n_data_compressed*sizeof(unsigned char));

  // scale data

  if(valmax>valmin){
    for(i=0;i<n_data_in;i++){
      c_data[i]=255*(data_in[i]-valmin)/(valmax-valmin);
    }
    CheckMemory;
  }
  else{
    memset(c_data,0,n_data_in);
    CheckMemory;
  }

  //  compress data

  compress(c_data_compressed+32,&n_data_compressedm32, c_data, n_data_in);

  n_data_compressed=n_data_compressedm32+32;
  CheckMemory;
// 1,completion,version
// 1,version,n_data_compressed,nbytes,n_data_in,time,valmin,valmax,data ....
  memcpy(c_data_compressed,&one,4);
  memcpy(c_data_compressed+4,&version,4);
  memcpy(c_data_compressed+8,&n_data_compressed,4);
  memcpy(c_data_compressed+12,&nbytes,4);
  memcpy(c_data_compressed+16,&n_data_in,4);
  memcpy(c_data_compressed+20,&timeval_in,4);
  memcpy(c_data_compressed+24,&valmin,4);
  memcpy(c_data_compressed+28,&valmax,4);
  CheckMemory;

  // resize and deallocate buffers

  ResizeMemory((void **)&c_data_compressed, n_data_compressed);
  FREEMEMORY(c_data);
  *compressed_data_out=c_data_compressed;
  *ncompressed_data_out=n_data_compressed;
}

/* ------------------ uncompress_volsliceframe ------------------------ */

int uncompress_volsliceframe(unsigned char *compressed_data_in,
                           float *data_out, int n_data_in, float *timeval_out,
                           unsigned char *fullbuffer
                ){
  float valmin, valmax;
  int i,ndatafile;
  uLongf countin,countout;

  valmin=*(float *)(compressed_data_in+24);
  valmax=*(float *)(compressed_data_in+28);
  *timeval_out=*(float *)(compressed_data_in+20);
  countin = *(int *)(compressed_data_in+8)-32;
  ndatafile = *(int *)(compressed_data_in+16);

  uncompress(fullbuffer,&countout,compressed_data_in+32,countin);

  if(countout==ndatafile&&n_data_in>=countout){
    for(i=0;i<countout;i++){
      data_out[i]=valmin+fullbuffer[i]*(valmax-valmin)/255.0;
    }
  }
  return countout;
}
