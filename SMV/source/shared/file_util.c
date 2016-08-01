#include "options.h"
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>
#include <stdlib.h>
#ifdef pp_OSX
#include <unistd.h>
#endif
#include <math.h>
#ifdef WIN32
#include <io.h>
#include <direct.h>
#include <dirent_win.h>
#else
#include <dirent.h>
#endif
#include "MALLOC.h"
#include "smv_endian.h"

FILE *alt_stdout=NULL;

/* ------------------ FFLUSH ------------------------ */

int FFLUSH(void){
  int return_val=0;

  if(alt_stdout!=NULL){
    return_val = fflush(alt_stdout);
  }
  return return_val;
}

/* ------------------ PRINTF ------------------------ */

int PRINTF(const char * format, ...){
  va_list args;
  int return_val=0;

  if(alt_stdout!=NULL){
    va_start(args, format);
    return_val=vfprintf(alt_stdout, format, args);
    va_end(args);
  }
  return return_val;
}

/* ------------------ set_stdout ------------------------ */

void set_stdout(FILE *stream){
  alt_stdout=stream;
}

#define FILE_BUFFER 1000

/* ------------------ copyfile ------------------------ */

void copyfile(char *destdir, char *file_in, char *file_out, int mode){
  char buffer[FILE_BUFFER];
  FILE *streamin=NULL, *streamout=NULL;
  char *full_file_out=NULL;
  size_t chars_in;

  if(destdir==NULL||file_in==NULL)return;
  streamin=fopen(file_in,"rb");
  if(streamin==NULL)return;

  full_file_out=NULL;
  NewMemory((void **)&full_file_out,strlen(file_out)+strlen(destdir)+1+1);
  strcpy(full_file_out,destdir);
  if(destdir[strlen(destdir)-1]!=*dirseparator){
    strcat(full_file_out,dirseparator);
  }
  strcat(full_file_out,file_out);

  if(mode==REPLACE_FILE){
    streamout=fopen(full_file_out,"wb");
  }
  else if(mode==APPEND_FILE){
    streamout=fopen(full_file_out,"ab");
  }
  else{
    ASSERT(0);
  }

  if(streamout==NULL){
    FREEMEMORY(full_file_out);
    fclose(streamin);
    return;
  }
  PRINTF("  Copying %s to %s\n",file_in,file_out);
  for(;;){
    int end_of_file;

    end_of_file=0;
    chars_in=fread(buffer,1,FILE_BUFFER,streamin);
    if(chars_in!=FILE_BUFFER)end_of_file=1;
    if(chars_in>0)fwrite(buffer,chars_in,1,streamout);
    if(end_of_file==1)break;
  }
  FREEMEMORY(full_file_out);
  fclose(streamin);
  fclose(streamout);
}

/* ------------------ have_prog ------------------------ */

int have_prog(char *prog){
  if(system(prog) == 0)return 1;
  return 0;
}

/* ------------------ get_smokezippath ------------------------ */

char *get_smokezippath(char *progdir){
  STRUCTSTAT statbuffer;
  char *zip_path;

  if(progdir!=NULL){
    NewMemory((void **)&zip_path,strlen(progdir)+20);
    strcpy(zip_path,progdir);
  }
  else{
    NewMemory((void **)&zip_path,2+20);
    strcpy(zip_path,".");
    strcat(zip_path,dirseparator);
  }

  strcat(zip_path,"smokezip");
#ifdef WIN32
  strcat(zip_path,".exe");
#endif
  if(STAT(zip_path,&statbuffer)==0)return zip_path;
  FREEMEMORY(zip_path);
  return NULL;
}

/* ------------------ setdir ------------------------ */

char *setdir(char *argdir){
  int lendir;
  char *dir;

  lendir=strlen(argdir);
  NewMemory((void **)&dir,lendir+2);
  strcpy(dir,argdir);
  if(dir[lendir-1]!=dirseparator[0]){
    strcat(dir,dirseparator);
  }
  return dir;
}

/* ------------------ get_basefilename ------------------------ */

char *get_basefilename(char *buffer,char *file){
  char *filebase,*ext;

  strcpy(buffer,file);
#ifdef WIN32
  filebase=strrchr(buffer,'\\');
#else
  filebase=strrchr(buffer,'/');
#endif
  if(filebase==NULL){
    filebase=buffer;
  }
  else{
    filebase++;
  }
  ext = strrchr(filebase,'.');
  if(ext!=NULL)*ext=0;
  return filebase;
}

/* ------------------ get_filename ------------------------ */

char *get_filename(char *temp_dir, char *file, int flag){
  char *file2;
  char *file_out=NULL;
  FILE *stream=NULL;

  trim_back(file);
  file2=trim_front(file);
  if(flag==0){
    stream=fopen(file2,"r");
    if(can_write_to_dir(".")==1||stream!=NULL){
      NewMemory((void **)&file_out,strlen(file2)+1);
      strcpy(file_out,file2);
    }
    if(stream!=NULL)fclose(stream);
  }
  if(file_out==NULL&&temp_dir!=NULL&&can_write_to_dir(temp_dir)==1){
    NewMemory((void **)&file_out,strlen(temp_dir)+1+strlen(file2)+1);
    strcpy(file_out,"");
    strcat(file_out,temp_dir);
    strcat(file_out,dirseparator);
    strcat(file_out,file);
  }
  return file_out;
}

/* ------------------ fullfile ------------------------ */

void fullfile(char *file_out, char *dir, char *file){
  char *file2;

  trim_back(file);
  file2=trim_front(file);
  strcpy(file_out,"");
  if(dir!=NULL)strcat(file_out,dir);
  strcat(file_out,file2);
}

/* ------------------ filecat ------------------------ */

int filecat(char *file_in1, char *file_in2, char *file_out){
  char buffer[FILE_BUFFER];
  FILE *stream_in1, *stream_in2, *stream_out;
  int chars_in;

  if(file_in1==NULL||file_in2==NULL)return -1;
  if(file_out==NULL)return -2;

  stream_in1=fopen(file_in1,"r");
  if(stream_in1==NULL)return -1;

  stream_in2=fopen(file_in2,"r");
  if(stream_in2==NULL){
    fclose(stream_in1);
    return -1;
  }

  stream_out=fopen(file_out,"w");
  if(stream_out==NULL){
    fclose(stream_in1);
    fclose(stream_in2);
    return -2;
  }

  for(;;){
    int end_of_file;

    end_of_file=0;
    chars_in=fread(buffer,1,FILE_BUFFER,stream_in1);
    if(chars_in!=FILE_BUFFER)end_of_file=1;
    if(chars_in>0)fwrite(buffer,chars_in,1,stream_out);
    if(end_of_file==1)break;
  }
  fclose(stream_in1);

  for(;;){
    int end_of_file;

    end_of_file=0;
    chars_in=fread(buffer,1,FILE_BUFFER,stream_in2);
    if(chars_in!=FILE_BUFFER)end_of_file=1;
    if(chars_in>0)fwrite(buffer,chars_in,1,stream_out);
    if(end_of_file==1)break;
  }
  fclose(stream_in2);
  fclose(stream_out);
  return 0;

}

/* ------------------ make_outfile ------------------------ */

void make_outfile(char *outfile, char *destdir, char *file1, char *ext){
  char filename_buffer[1024], *file1_noext;

  trim_back(file1);
  strcpy(filename_buffer,trim_front(file1));
  file1_noext=strstr(filename_buffer,ext);
  strcpy(outfile,"");
  if(file1_noext==NULL)return;
  file1_noext[0]='\0';
  if(destdir!=NULL){
    strcpy(outfile,destdir);
  }
  strcat(outfile,filename_buffer);
  strcat(outfile,"_diff");
  strcat(outfile,ext);
}

/* ------------------ can_write_to_dir ------------------------ */

int can_write_to_dir(char *dir){

// returns 1 if the directory can be written to, 0 otherwise

  char *full_name;
  char file_name[256], *file_name_ptr;
  FILE *stream;
  int len;
  int return_val=0;

  if(dir==NULL||strlen(dir)==0)return 0;

  file_name_ptr=randstr(file_name,20);
  if(file_name_ptr==NULL)return 0;

  len = strlen(dir) + 20 + 1 + 1;

  NewMemory((void **)&full_name,len);

  strcpy(full_name,"");
  if(strcmp(dir,".")!=0&&strlen(dir)>0){
    strcat(full_name,dir);
    strcat(full_name,dirseparator);
  }

  strcat(full_name,file_name_ptr);

  stream=fopen(full_name,"wb");
  if(stream!=NULL){
    fclose(stream);
    return_val=1;
  }
  UNLINK(full_name);
  FREEMEMORY(full_name);
  return return_val;
}

/* ------------------ is_file_newer ------------------------ */

int is_file_newer(char *file1, char *file2){

// returns 1 if file1 is newer than file2, 0 otherwise

  STRUCTSTAT statbuff1, statbuff2;
  int statfile1, statfile2;

  if(file1==NULL||file2==NULL)return -1;

  statfile1=STAT(file1,&statbuff1);
  statfile2=STAT(file2,&statbuff2);
  if(statfile1!=0||statfile2!=0)return -1;

  if(statbuff1.st_mtime>statbuff2.st_mtime)return 1;
  return 0;
}

  /* ------------------ getfileinfo ------------------------ */

int getfileinfo(char *filename, char *source_dir, FILE_SIZE *filesize){
  STRUCTSTAT statbuffer;
  int statfile;
  char buffer[1024];

  if(source_dir==NULL){
    strcpy(buffer,filename);
  }
  else{
    strcpy(buffer,source_dir);
    strcat(buffer,filename);
  }
  if(filesize!=NULL)*filesize=0;
  statfile=STAT(buffer,&statbuffer);
  if(statfile!=0)return statfile;
  if(filesize!=NULL)*filesize=statbuffer.st_size;
  return statfile;
}

/* ------------------ get_filesize ------------------------ */

FILE_SIZE get_filesize(const char *filename){
  STRUCTSTAT statbuffer;
  int statfile;
  FILE_SIZE return_val;

  return_val=0;
  if(filename==NULL)return return_val;
  statfile=STAT(filename,&statbuffer);
  if(statfile!=0)return return_val;
  return_val = statbuffer.st_size;
  return return_val;
}

  /* ------------------ file_exists ------------------------ */

int file_exists(char *filename){

// returns 1 if the file filename exists, 0 otherwise

#ifdef WIN32
  if(filename==NULL||_access(filename,0)==-1){
    return 0;
  }
  else{
    return 1;
  }
#else
  STRUCTSTAT statbuffer;

  if(filename==NULL||STAT(filename,&statbuffer)!=0){
    return 0;
  }
  else{
    return 1;
  }
#endif
}

/* ------------------ free_filelist ------------------------ */

void free_filelist(filelistdata *filelist, int *nfilelist){
  int i;

  for(i=0;i<*nfilelist;i++){
    FREEMEMORY(filelist[i].file);
  }
  FREEMEMORY(filelist);
  *nfilelist=0;
}

  /* ------------------ get_nfilelist ------------------------ */

int get_nfilelist(const char *path, char *key){
  struct dirent *entry;
  DIR *dp;
  int maxfiles=0;

  dp = opendir(path);
  if(dp == NULL){
    perror("opendir");
    return 0;
  }
  while( (entry = readdir(dp)) ){
  //  if((entry->d_type==DT_DIR&&entry->d_name[0]!='.')||(entry->d_type==DT_REG&&match_wild(entry->d_name,key)==1)){
    if((entry->d_type==DT_REG&&match_wild(entry->d_name,key)==1)){
      maxfiles++;
      continue;
    }
  }
  closedir(dp);
  return maxfiles;
}

 /* ------------------ get_filelist ------------------------ */

int get_filelist(const char *path, char *key, int maxfiles, filelistdata **filelist){
  struct dirent *entry;
  DIR *dp;
  int nfiles=0;
  filelistdata *flist;

  // DT_DIR - is a directory
  // DT_REG - is a regular file

  dp = opendir(path);
  if(dp == NULL){
    perror("opendir");
    *filelist=NULL;
    return 0;
  }
  if(maxfiles==0){
    closedir(dp);
    *filelist=NULL;
    return 0;
  }
  NewMemory((void **)&flist,maxfiles*sizeof(filelistdata));
  /*
  while( (entry = readdir(dp))&&nfiles<maxfiles ){
    if(entry->d_type==DT_DIR&&entry->d_name[0]!='.'){
      char *file;
      filelistdata *flisti;

      flisti = flist + nfiles;
      NewMemory((void **)&file,strlen(entry->d_name)+1);
      strcpy(file,entry->d_name);
      flisti->file=file;
      flisti->type=1;
      nfiles++;
    }
  }
  rewinddir(dp);
  */
  while( (entry = readdir(dp))&&nfiles<maxfiles ){
    if(entry->d_type==DT_REG&&match_wild(entry->d_name,key)==1){
      char *file;
      filelistdata *flisti;

      flisti = flist + nfiles;
      NewMemory((void **)&file,strlen(entry->d_name)+1);
      strcpy(file,entry->d_name);
      flisti->file=file;
      flisti->type=0;
      nfiles++;
    }
  }
  *filelist=flist;
  closedir(dp);
  return nfiles;
}

/* ------------------ getfilesizelabel ------------------------ */

void getfilesizelabel(int size, char *sizelabel){
  int leftsize,rightsize;

#define sizeGB   1000000000
#define size100MB 100000000
#define size10MB   10000000
#define sizeMB     1000000
#define size100KB    100000
#define size10KB      10000

  if(size>=sizeGB){
    size/=size10MB;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i GB",leftsize,rightsize);
  }
  else if(size>=size100MB&&size<sizeGB){
    size/=sizeMB;
    leftsize=size;
    sprintf(sizelabel,"%i MB",leftsize);
  }
  else if(size>=size10MB&&size<size100MB){
    size/=size100KB;
    leftsize=size/10;
    rightsize=size-10*leftsize;
    sprintf(sizelabel,"%i.%i MB",leftsize,rightsize);
  }
  else if(size>=sizeMB&&size<size10MB){
    size/=size10KB;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i MB",leftsize,rightsize);
  }
  else if(size>=size100KB&&size<sizeMB){
    size/=1000;
    leftsize=size;
    sprintf(sizelabel,"%i KB",leftsize);
  }
  else{
    size/=10;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i KB",leftsize,rightsize);
  }
}

/* ------------------ getprogdir ------------------------ */

char *getprogdir(char *progname, char **svpath){

// returns the directory containing the file progname

  char *progpath, *lastsep, *smokeviewpath2;

  lastsep=strrchr(progname,dirseparator[0]);
  if(lastsep==NULL){
    char *dir;

    dir = which(progname);
    if(dir==NULL){
      NewMemory((void **)&progpath,(unsigned int)3);
      strcpy(progpath,".");
      strcat(progpath,dirseparator);
    }
    else{
      int lendir;

      lendir=strlen(dir);
      NewMemory((void **)&progpath,(unsigned int)(lendir+2));
      strcpy(progpath,dir);
      if(progpath[lendir-1]!=dirseparator[0])strcat(progpath,dirseparator);
    }
    NewMemory((void **)&smokeviewpath2,(unsigned int)(strlen(progpath)+strlen(progname)+1));
    strcpy(smokeviewpath2,progpath);
  }
  else{
    int lendir;

    lendir=lastsep-progname+1;
    NewMemory((void **)&progpath,(unsigned int)(lendir+1));
    strncpy(progpath,progname,lendir);
    progpath[lendir]=0;
    NewMemory((void **)&smokeviewpath2,(unsigned int)(strlen(progname)+1));
    strcpy(smokeviewpath2,"");;
  }
  strcat(smokeviewpath2,progname);
  *svpath=smokeviewpath2;
  return progpath;
}

/* ------------------ lastname ------------------------ */

char *lastname(char *argi){

// returns the file name contained in the full path name argi

  char *lastdirsep;
  char *dir, *filename, cwdpath[1000];

#ifdef WIN32
#define CHDIR _chdir
#define GETCWD _getcwd
#define SEP '\\'
#else
#define CHDIR chdir
#define GETCWD getcwd
#define SEP '/'
#endif

  filename=argi;
  lastdirsep=strrchr(argi,SEP);
  if(lastdirsep!=NULL){
    dir=argi;
    filename=lastdirsep+1;
    lastdirsep[0]=0;
    GETCWD(cwdpath,1000);
    if(strcmp(cwdpath,dir)!=0){
      CHDIR(dir);
    }
  }
  return filename;
}

/* ------------------ get_zonefilename ------------------------ */

char *get_zonefilename(char *bufptr){
  char *full_name, *last_name, *filename;
  STRUCTSTAT statbuffer;

  full_name=bufptr;
  if(STAT(full_name,&statbuffer)!=0)full_name=NULL;

  last_name=lastname(bufptr);
  if(STAT(last_name,&statbuffer)!=0)last_name=NULL;

  if(last_name!=NULL&&full_name!=NULL){
    if(strcmp(last_name,full_name)==0){
      last_name=NULL;
    }
  }

  if(last_name!=NULL&&full_name!=NULL){
    filename=last_name;
  }
  else if(last_name==NULL&&full_name!=NULL){
    filename=full_name;
  }
  else if(last_name!=NULL&&full_name==NULL){
    filename=last_name;
  }
  else{
    filename=NULL;
  }
  return filename;
}

/* ------------------ file_modtime ------------------------ */

time_t file_modtime(char *filename){

// returns the modification time of the file named filename

  STRUCTSTAT statbuffer;
  time_t return_val;
  int statfile;

  return_val=0;
  if(filename==NULL)return return_val;
  statfile=STAT(filename,&statbuffer);
  if(statfile!=0)return return_val;
  return_val = statbuffer.st_mtime;
  return return_val;
}

/* ------------------ which ------------------------ */

char *which(char *progname){

// returns the PATH directory containing the file progname

  char *pathlist, *pathlistcopy, *fullprogname, *prognamecopy;
  char *dir,*pathentry;
  char pathsep[2], dirsep[2];

#ifdef WIN32
  strcpy(pathsep,";");
  strcpy(dirsep,"\\");
#else
  strcpy(pathsep,":");
  strcpy(dirsep,"/");
#endif

  pathlist = getenv("PATH");
  if(pathlist==NULL||strlen(pathlist)==0||progname==NULL||strlen(progname)==0)return NULL;

  NewMemory((void **)&prognamecopy, (unsigned int)(strlen(progname)+4+1));
  strcpy(prognamecopy, progname);

  NewMemory((void **)&pathlistcopy, (unsigned int)(strlen(pathlist)+1));
  strcpy(pathlistcopy, pathlist);

#ifdef WIN32
  {
    const char *ext;

    ext = prognamecopy+strlen(progname)-4;
    if(strlen(progname)<=4||STRCMP(ext,".exe")!=0)strcat(prognamecopy, ".exe");
  }
#endif

  NewMemory((void **)&fullprogname, (unsigned int)(strlen(progname)+4+strlen(dirsep)+strlen(pathlist)+1));

  dir=strtok(pathlistcopy,pathsep);
  while(dir!=NULL&&strlen(dir)>0){
    strcpy(fullprogname,dir);
    strcat(fullprogname,dirsep);
    strcat(fullprogname,prognamecopy);
    if(file_exists(fullprogname)==1){
      NewMemory((void **)&pathentry,(unsigned int)(strlen(dir)+2));
      strcpy(pathentry,dir);
      strcat(pathentry,dirsep);
      FREEMEMORY(pathlistcopy);
      FREEMEMORY(fullprogname);
      FREEMEMORY(prognamecopy);
      return pathentry;
    }
    dir=strtok(NULL,pathsep);
  }
  FREEMEMORY(pathlistcopy);
  FREEMEMORY(fullprogname);
  FREEMEMORY(prognamecopy);
  return NULL;
}
