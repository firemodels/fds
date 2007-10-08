#include "options.h"
#ifdef EGZ
#include <stdlib.h>
#include <string.h>
#include "egz_stdio.h"
#include "MALLOC.h"

// svn revision character string
char egz_stdio_revision[]="$Revision$";
// $Date$ $Author$


int getendian(void);

/* ------------------ EGZ_FCLOSE ------------------------ */

int EGZ_FCLOSE(EGZ_FILE *egz_stream){
  int return_code;

#ifdef USE_ZLIB
  if(egz_stream->compression==1){
    ASSERT(egz_stream->stream!=NULL);
    return_code=gzclose(egz_stream->stream);
  }
  else{
    ASSERT(egz_stream->stream2!=NULL);
    return_code=fclose(egz_stream->stream2);
  }
#else
  ASSERT(egz_stream->stream!=NULL);
  return_code=fclose(egz_stream->stream);
#endif
  FREEMEMORY(egz_stream);
  return return_code;
}


/* ------------------ EGZ_FOPEN ------------------------ */

EGZ_FILE *EGZ_FOPEN(const char *file, const char *mode, int compress, int endian){
  /*
  endian - format (big or little endian) used in file.
           endian=0 ==> little endian (LINUX, Windows PC's)
           endian=1 ==> big endian (Unix workstation)
           endian=2 ==> unknown, an integer 1 is assumed to exist in the first four bytes of file
  */

  EGZ_FILE *egz_stream;
  int testval;

#ifdef USE_ZLIB
  gzFile *stream;
  FILE *stream2;
  char gzfile[1024];
#else
  FILE *stream;
#endif
  int native_endian;
  int endianswitch_save;
  int one=1;

  egz_stream=NULL;
  if(NewMemory((void **)&egz_stream,sizeof(EGZ_FILE))==0)return NULL;

  native_endian = getendian();
  if((native_endian==1&&endian==0)||(native_endian==0&&endian==1)){
    egz_stream->endianswitch=1;
  }
  else{
    egz_stream->endianswitch=0;
  }

  switch (mode[0]){
  case 'r':
#ifdef USE_ZLIB
    strcpy(gzfile,file);
    strcat(gzfile,".gz");
    stream = gzopen(gzfile,mode);
    if(stream!=NULL){
      egz_stream->compression=1;
      egz_stream->stream = stream;
      egz_stream->stream2 = NULL;
    }
    else{
      stream2 = fopen(file,mode);
      egz_stream->compression=0;
      egz_stream->stream=NULL;
      egz_stream->stream2=stream2;
      if(stream2==NULL){
        FREEMEMORY(egz_stream);
        return NULL;
      }
    }
#else
    stream=fopen(file,mode);
    egz_stream->stream=stream;
    if(stream==NULL){
      FREEMEMORY(egz_stream);
      return NULL;
    }
#endif
    if(endian==2){
      egz_stream->endianswitch=0;
      testval=0;
      EGZ_FREAD( &testval, sizeof(int), 1, egz_stream );
      EGZ_REWIND(egz_stream);
      if(testval==1)return egz_stream;
      egz_stream->endianswitch=1;
      EGZ_FREAD( &testval, sizeof(int), 1, egz_stream );
      EGZ_REWIND(egz_stream);
      if(testval==1)return egz_stream;
      FREEMEMORY(egz_stream);
      return NULL;
    }
    break;
  case 'w':
  case 'a':
#ifdef USE_ZLIB
    if(compress==1){
      strcpy(gzfile,file);
      strcat(gzfile,".gz");
      stream = gzopen(gzfile,mode);
      if(stream==NULL){
        FREEMEMORY(egz_stream);
        return NULL;
      }
      if(stream!=NULL){
        egz_stream->compression=1;
        egz_stream->stream = stream;
        egz_stream->stream2 = NULL;
      }
    }
    else{
      stream2 = fopen(file,mode);
      egz_stream->compression=0;
      egz_stream->stream=NULL;
      egz_stream->stream2=stream2;
      if(stream2==NULL){
        FREEMEMORY(egz_stream);
        return NULL;
      }
    }
#else
    stream=fopen(file,mode);
    egz_stream->stream=stream;
    if(stream==NULL){
      FREEMEMORY(egz_stream);
      return NULL;
    }
#endif
    if(endian==2&&mode[0]=='w'){
      endianswitch_save=egz_stream->endianswitch;
      egz_stream->endianswitch=0;
      EGZ_FWRITE(&one,4,1,egz_stream);
      egz_stream->endianswitch=endianswitch_save;
    }
    break;
  default:
    ASSERT(0);
    break;
  }
  return egz_stream;
}
//  size_t fwrite( const void *buffer, size_t size, size_t count, FILE *stream );


/* ------------------ EGZ_FWRITE ------------------------ */

size_t EGZ_FWRITE( void *buffer, size_t size, size_t count, const EGZ_FILE *egz_stream ){
  int ii;
  size_t return_size;

  unsigned char c1, c2, c3, c4;
  unsigned char *ca, *cb, *cc, *cd;
  size_t i;

  if(egz_stream->endianswitch==1){
    switch (size){
    case 1:
      break;
    case 2:
      for(i=0;i<count;i++){
        ii = 2*i;
        ca=(unsigned char *)buffer+ii;
        cb=(unsigned char *)buffer+ii+1;
        c1=*ca;
        c2=*cb;
        *cb=c1;
        *ca=c2;
      }
      break;
    case 4:
      for(i=0;i<count;i++){
        ii=4*i;
        ca=(unsigned char *)buffer+ii;
        cb=(unsigned char *)buffer+ii+1;
        cc=(unsigned char *)buffer+ii+2;
        cd=(unsigned char *)buffer+ii+3;
        c1=*ca;
        c2=*cb;
        c3=*cc;
        c4=*cd;
        *ca=c4;
        *cb=c3;
        *cc=c2;
        *cd=c1;
      }
      break;
    default:
      ASSERT(0);
    }
  }
#ifdef USE_ZLIB
  if(egz_stream->compression==1){
    return_size = (size_t)gzwrite(egz_stream->stream,buffer,size*count);
  }
  else{
    return_size = fwrite(buffer,size,count,egz_stream->stream2);
  }
#else
  return_size = fwrite(buffer,size,count,egz_stream->stream);
#endif
  return return_size;
}

/* ------------------ EGZ_FREAD ------------------------ */

size_t EGZ_FREAD( void *buffer, size_t size, size_t count, const EGZ_FILE *egz_stream ){
  int ii;
  size_t return_size;
  unsigned char c1, c2, c3, c4;
  unsigned char *ca, *cb, *cc, *cd;
  size_t i;

#ifdef USE_ZLIB
  if(egz_stream->compression==1){
    return_size = (size_t)gzread(egz_stream->stream,buffer,size*count);
  }
  else{
    return_size = fread(buffer,size,count,egz_stream->stream2);
  }
#else
  return_size = fread(buffer,size,count,egz_stream->stream);
#endif
  if(egz_stream->endianswitch==0)return return_size;
  switch (size){
  case 1:
    return return_size;
  case 2:
    for(i=0;i<count;i++){
      ii = 2*i;
      ca=(unsigned char *)buffer+ii;
      cb=(unsigned char *)buffer+ii+1;
      c1=*ca;
      c2=*cb;
      *cb=c1;
      *ca=c2;
    }
    break;
  case 4:
    for(i=0;i<count;i++){
      ii=4*i;
      ca=(unsigned char *)buffer+ii;
      cb=(unsigned char *)buffer+ii+1;
      cc=(unsigned char *)buffer+ii+2;
      cd=(unsigned char *)buffer+ii+3;
      c1=*ca;
      c2=*cb;
      c3=*cc;
      c4=*cd;
      *ca=c4;
      *cb=c3;
      *cc=c2;
      *cd=c1;
    }
    break;
  default:
    ASSERT(0);
  }
  return return_size;
}

/* ------------------ EGZ_FEOF ------------------------ */

int EGZ_FEOF(const EGZ_FILE *egz_stream ){  
  int return_val;

#ifdef USE_ZLIB
  if(egz_stream->compression==1){
    return_val = gzeof(egz_stream->stream);
  }
  else{
    return_val = feof(egz_stream->stream2);
  }
#else
  return_val = feof(egz_stream->stream);
#endif

  return return_val;
}

void EGZ_REWIND(const EGZ_FILE *egz_stream){
#ifdef USE_ZLIB
  if(egz_stream->compression==1){
    gzrewind(egz_stream->stream);
  }
  else{
    rewind(egz_stream->stream2);
  }
#else
  rewind(egz_stream->stream);
#endif

}
/* ------------------ EGZ_FSEEK ------------------------ */

int EGZ_FSEEK( const EGZ_FILE *egz_stream, long offset, int origin ){
  int return_val;
  
#ifdef USE_ZLIB
  if(egz_stream->compression==1){
    return_val = gzseek(egz_stream->stream,offset,origin);
  }
  else{
    return_val = fseek(egz_stream->stream2,offset,origin);
  }
#else
  return_val = fseek(egz_stream->stream,offset,origin);
#endif

  return return_val;
}

/* ------------------ EGZ_FTELL ------------------------ */

long EGZ_FTELL( const EGZ_FILE *egz_stream ){
  long return_val;

#ifdef USE_ZLIB
  if(egz_stream->compression==1){
    return_val =  gztell(egz_stream->stream);
  }
  else{
    return_val =  ftell(egz_stream->stream2);
  }
#else
  return_val =  ftell(egz_stream->stream);
#endif
  return return_val;
}
#endif
