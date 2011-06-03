// $Date$ 
// $Revision$
// $Author$
#define IN_STRING_UTIL
#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#ifdef WIN32
#include <direct.h>
#endif
#include "MALLOC.h"
#include "string_util.h"

// svn revision character string
char string_util_revision[]="$Revision$";

/* ----------------------- fparsecsv ----------------------------- */

void fparsecsv(char *buffer, float *vals, int ncols, int *ntokens){
  /*! \fn void fparsecsv(char *buffer, float *vals, int ncols, int *ntokens)
      \brief copy comma delimited values from buffer into floating point array vals
       returning number of values found in ntokens
  */
  int nt=0;
  char *token;

  token=strtok(buffer,",");
  while(token!=NULL&&nt<ncols){
    sscanf(token,"%f",&vals[nt++]);
    token=strtok(NULL,",");
  }
  *ntokens=nt;
}

/* ----------------------- parsecsv ----------------------------- */

void parsecsv(char *buffer, char **tokens, int ncols, int *ntokens){
  /*! \fn void parsecsv(char *buffer, char **tokens, int ncols, int *ntokens)
      \brief copy comma delimited values from buffer into character array tokens
       returning number of values found in ntokens
  */
  int nt=0;
  char *token;

  token=strtok(buffer,",");
  while(token!=NULL&&nt<ncols){
    tokens[nt++]=token;
    token=strtok(NULL,",");
  }
  *ntokens=nt;
}

/* ------------------ getrowcols ------------------------ */

void getrowcols(FILE *stream, int *nrows, int *ncols){
  /*! \fn void getrowcols(FILE *stream, int *nrows, int *ncols)
      \brief find number of rows (nrows) and number of columns (ncols) in
       comma delimited file pointed to by stream
  */
  char buffer[1024];
  int nnrows=0,nncols=0;

  while(!feof(stream)){
    if(fgets(buffer,1024,stream)==NULL)break;
    if(nnrows==0){
      nncols=getcols(buffer);
    }
    nnrows++;
  }
  *nrows=nnrows;
  *ncols=nncols;
}

/* ------------------ getcols ------------------------ */

int getcols(char *buffer){
  /*! \fn int getcols(char *buffer)
      \brief return the number of columns in the character string buffer
  */
  char buffer2[1024];
  int ncols=0;
  char *comma;

  strcpy(buffer2,buffer);
  trim(buffer2);
  comma=strtok(buffer2,",");

  while(comma!=NULL){
    ncols++;
    comma=strtok(NULL,",");
  }
  return ncols;
}

/* ------------------ stripquotes ------------------------ */

void stripquotes(char *buffer){
  /*! \fn void stripquotes(char *buffer)
      \brief replaces quotes (") with blanks in the character string buffer
  */
  int i;
  char *c;

  for(i=0;i<strlen(buffer);i++){
    c=buffer+i;
    if(*c=='"')*c=' ';
  }
}
/* ------------------ stripcommas ------------------------ */

void stripcommas(char *buffer){
  /*! \fn void stripcommas(char *buffer)
      \brief replaces commas (,) with blanks in the character string buffer
  */
  int i;
  char *c;

  for(i=0;i<strlen(buffer);i++){
    c=buffer+i;
    if(*c==',')*c=' ';
  }
}

/* ------------------ file_modtime ------------------------ */

time_t file_modtime(char *filename){
  /*! \fn time_t file_modtime(char *filename)
      \brief returns the modification time of the file named filename
  */
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

/* ------------------ randint ------------------------ */

int randint(int min, int max){
  /*! \fn int randint(int min, int max)
      \brief returns a random integer inclusively between min and max 
  */
  int return_val;

  if (min>max){
    return_val = max+((min-max+1)*rand()/(RAND_MAX+1.0));
  }
  else{
    return_val = min+((max-min+1)*rand()/(RAND_MAX+1.0));
  }
  return return_val;
}

/* ------------------ randstr ------------------------ */

char *randstr(char* str, int length){
  /*! \fn char *randstr(char* str, int length)
      \brief returns a random character string of length length
  */
    int i;

    if (str==NULL||length<=0)return NULL;

    for (i=0;i<length;i++){
      str[i]=(char)randint(65,90);
    }
    str[length]=0;
    return str;
}

/* ------------------ can_write_to_dir ------------------------ */

int can_write_to_dir(char *dir){
  /*! \fn int can_write_to_dir(char *dir)
      \brief returns 1 if the directory can be written to, 0 otherwise
  */
  char full_name[1024];
  char file_name[1024], *file_name_ptr;
  FILE *stream;

  file_name_ptr=randstr(file_name,20);
  if(file_name_ptr==NULL)return 0;
  strcpy(full_name,"");
  if(dir!=NULL&&strcmp(dir,".")!=0&&strlen(dir)>0){
    strcat(full_name,dir);
    strcat(full_name,dirseparator);
  }

  strcat(full_name,file_name_ptr);
  
  stream=fopen(full_name,"wb");
  if(stream==NULL)return 0;
  fclose(stream);
  remove(full_name);
  return 1;
}

/* ------------------ is_file_newer ------------------------ */

int is_file_newer(char *file1, char *file2){
  /*! \fn int is_file_newer(char *file1, char *file2)
      \brief returns 1 if file1 is newer than file2, 0 otherwise
  */
  STRUCTSTAT statbuff1, statbuff2;
  int statfile1, statfile2;

  if(file1==NULL||file2==NULL)return -1;

  statfile1=STAT(file1,&statbuff1);
  statfile2=STAT(file2,&statbuff2);
  if(statfile1!=0||statfile2!=0)return -1;

  if(statbuff1.st_mtime>statbuff2.st_mtime)return 1;
  return 0;
}

/* ------------------ rootdir ------------------------ */

char *getprogdir(char *progname){
  /*! \fn char *getprogdir(char *progname)
      \brief returns the directory containing the file progname
  */
  char *progpath, *lastsep;
#ifdef WIN32
  char cdirsep='\\';
#else
  char cdirsep='/';
#endif

  lastsep=strrchr(progname,cdirsep);
  if(lastsep==NULL){
    char *dir;

    dir = which(progname);
    if(dir==NULL){
      NewMemory((void **)&progpath,(unsigned int)3);
      strcpy(progpath,".");
      strcat(progpath,dirseparator);
      return progpath;
    }
    else{
      int lendir;

      lendir=strlen(dir);
      NewMemory((void **)&progpath,(unsigned int)(lendir+2));
      strcpy(progpath,dir);
      if(progpath[lendir-1]!=cdirsep)strcat(progpath,dirseparator);
    }
    return progpath;
  }
  else{
    int lendir;

    lendir=lastsep-progname+1;
    NewMemory((void **)&progpath,(unsigned int)(lendir+1));
    strncpy(progpath,progname,lendir);
    progpath[lendir]=0;
    return progpath;
  }
}

/* ------------------ lastname ------------------------ */

char *lastname(char *argi){
  /*! \fn char *lastname(char *argi)
      \brief returns the file name contained in the full path name argi
  */
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

/* ------------------ trim ------------------------ */

void trim(char *line){
  /*! \fn void trim(char *line)
      \brief removes trailing white space from the character string line
  */
  char *blank=" ";
  const char *c;
  const char *lf="\n", *cr="\r";
  size_t len, i;
  
  len = strlen(line);
  c = line+len-1;
  for(i=0; i<len; i++){
    if(strncmp(c,blank,1)!=0&&strncmp(c,lf,1)!=0&&strncmp(c,cr,1)!=0){
      c++; 
      line[c-line]='\0';
      return;
    }
    c--;
  }
  *line='\0';
}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  /*! \fn char *trim_front(char *line)
      \brief returns a pointer to the first non-blank character in the character string line
  */
  char *blank=" ";
  const char *c;
  size_t i,len;

  c = line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(strncmp(c++,blank,1)!=0)return line+i;
  }
  return line;
}


/* ------------------ trimzeros ------------------------ */

void trimzeros(char *line){
  /*! \fn void trimzeros(char *line)
      \brief removes trailing zeros in the floating point number found in line
  */
  size_t i,len;
  char *c;

  len = strlen(line);
  c = line + len-1;
  for(i=len-1;i>0;i--){
    if(*c=='0'){
      c--;
      if(*c=='.'){
        line[i+1]='\0';
        return;
      }
      continue;
    }
    line[i+1]='\0';
    return;
  }
  line[0]='\0';
}

/* ------------------ trimmzeros ------------------------ */

void trimmzeros(char *line){
  /*! \fn void trimmzeros(char *line)
      \brief removes trailing zeros in each floating point number found in line
  */
  char linecopy[1024];
  char *token;

  trim(line);
  strcpy(linecopy,line);
  token=strtok(linecopy," ");
  strcpy(line,"");
  while(token!=NULL){
    trimzeros(token);
    strcat(line,token);
    strcat(line," ");
    token=strtok(NULL," ");
  }
}

/* ------------------ STRSTR ------------------------ */

char *STRSTR(char *c, const char *key){
  /*! \fn char *STRSTR(char *c, const char *key)
      \brief looks for first occurence of the string key found in c ignoring case
  */
  char *C,*CCOPY,*CC,*cc,*result;
  char *KEY,*KEYCOPY,*KEY2;
  size_t i, len,len2;
  int diff;

  if(c==NULL||key==NULL)return NULL;
  len=strlen(c);
  len2=strlen(key);
  if(len<1||len2<1)return NULL;
  if(NewMemory((void **)&C,(unsigned int)(len+1))==0)return NULL;
  CC=C;
  cc=c;
  CCOPY=C;
  if(NewMemory((void **)&KEY,(unsigned int)(len2+1))==0){
    FreeMemory(C);
    return NULL;
  }
  KEY2=KEY;
  KEYCOPY=KEY;
  for(i=0;i<len;i++){
    *CC++=(char)toupper(*cc++);
  }
  for(i=0;i<len2;i++){
    *KEY2++=(char)toupper(*key++);
  }
  *CC='\0';
  *KEY2='\0';
  result = strstr(C,KEY);
  if(result!=NULL)diff = result - C;
  FREEMEMORY(CCOPY);
  FREEMEMORY(KEYCOPY);
  if(result==NULL)return NULL;
  return c + diff;

}

/* ------------------ scalestring ------------------------ */

void scalestring(const char *stringfrom, char *stringto, const float *scale, float range){
  float val;
  sscanf(stringfrom,"%f",&val);
  val = scale[0]*val+scale[1];
  num2string(stringto,val,scale[0]*range);
}

/* ------------------ scalefloat2string ------------------------ */

void scalefloat2string(float floatfrom, char *stringto, const float *scale, float range){
  float val;
  val = scale[0]*floatfrom+scale[1];
  num2string(stringto,val,scale[0]*range);
}

/* ------------------ num2string ------------------------ */

void num2string(char *string, float tval,float range){
  float tval2,mant10;
  int exp10;

  //if(fabs((double)tval)<fabs((double)range)/100.0f)tval=0.0f;
  tval2=tval; 
  if(tval2<0.0)tval2=-tval2;
  if(0.01<=tval2&&tval2<0.1){
    sprintf(string,"%3.2f",tval);
  }
  else if(0.1<=tval2&&tval2<1.0){
    sprintf(string,"%3.2f",tval);
  }
  else if(1.0<=tval2&&tval2<10.0){
    sprintf(string,"%3.2f",tval);
  }
  else if(10.0<=tval2&&tval2<100.0){
    sprintf(string,"%3.1f",tval);
  }
  else if(100.0<=tval2&&tval2<1000.0){
    sprintf(string,"%3.0f",tval);
  }
  else if(1000.0<=tval2&&tval2<10000.0){
    sprintf(string,"%4.0f",tval);
  }
  else if(10000.0<=tval2&&tval2<100000.0){
    sprintf(string,"%5.0f",tval);
    }
  else if(tval2==0.0){STRCPY(string,"0.00");}
  else{
    mant10 = frexp10(tval,&exp10);
    mant10 = (float)((int)(10.0f*mant10+0.5f))/10.0f;
    if(mant10>=10.0f){
      mant10/=10.0f;
      exp10++;
    }
    if(exp10<-99)STRCPY(string,"0.00");
    else if(exp10>=-99&&exp10<-9){sprintf(string,"%2.1f%i",mant10,exp10);}
    else if(exp10>99)STRCPY(string,"***");
    else{
      if(exp10==0){sprintf(string,"%2.1f",mant10);}
      else{sprintf(string,"%2.1fE%i",mant10,exp10);}
    }

    /*sprintf(string,"%1.1e",tval); */
  }
  if(strlen(string)>9)printf("ut oh - overwriting string\n");


}

/* ------------------ get_string ------------------------ */

char *trim_string(char *buffer){
  /*! \fn char *trim_string(char *buffer)
      \brief removes trailing blanks from buffer and returns a pointer to the first non-blank character
  */
  int len;
  char *bufptr;

  len=strlen(buffer);
  buffer[len-1]='\0';
  bufptr=trim_front(buffer);
  trim(bufptr);
  return bufptr;
}


/* ------------------ STRCMP ------------------------ */

int STRCMP(const char *s1, const char *s2){
  /*! \fn int STRCMP(const char *s1, const char *s2)
      \brief same as the standard function, strcmp, but ignores case
  */
  while (toupper(*s1) == toupper(*s2++)){
		if (*s1++ == 0)return (0);
  }
	return (toupper(*(const unsigned char *)s1) - toupper(*(const unsigned char *)(s2 - 1)));
}

/* ------------------ get_chid ------------------------ */

char *get_chid(char *file, char *buffer){
  FILE *stream;
  char *chidptr,*c;
  unsigned int i;
  int found1st, found2nd;

  if(file==NULL)return NULL;
  stream=fopen(file,"r");
  if(stream==NULL)return NULL;

  while(!feof(stream)){
    found1st=0;
    found2nd=0;

    if(fgets(buffer,255,stream)==NULL)break;
    chidptr=strstr(buffer,"CHID");
    if(chidptr==NULL)continue;

    chidptr+=5;
    for(i=0;i<strlen(chidptr);i++){
      c=chidptr+i;
      if(*c=='\''){
        found1st=1;
        chidptr=c+1;
        break;
      }
    }
    if(found1st==0)break;

    for(i=0;i<strlen(chidptr);i++){
      c=chidptr+i;
      if(*c=='\''){
        found2nd=1;
        *c=0;
        break;
      }
    }
    break;
  }
  fclose(stream);
  if(found1st==0||found2nd==0)chidptr=NULL;
  return chidptr;
}
#ifdef pp_GPU

/* ------------------ log_base2 ------------------------ */

int log_base2(float xx){
  /*! \fn int log_base2(float xx)
      \brief returns the log base 2 of the floating point number xx
  */
  int r = 0;
  int x;
  x=xx;
  while( (x >> r) != 0){
    r++;
  }
  return r-1; // returns -1 for x==0, floor(log2(x)) otherwise
}
#endif

/* ------------------ array2string ------------------------ */

void array2string(float *vals, int nvals, char *string){
  /*! \fn void array2string(float *vals, int nvals, char *string)
      \brief convert an array of floating point numbers to a character string
  */
  char cval[30];
  int i;

  strcpy(string,"");
  if(nvals==0)return;
  for(i=0;i<nvals-1;i++){
    sprintf(cval,"%f",vals[i]);
    trimzeros(cval);
    strcat(string,cval);
    strcat(string,", ");
  }
  sprintf(cval,"%f",vals[nvals-1]);
  trimzeros(cval);
  strcat(string,cval);
}

  /* ------------------ getfileinfo ------------------------ */

int file_exists(char *filename){
  /*! \fn int file_exists(char *filename)
      \brief returns 1 if the file filename exists, 0 otherwise
  */
  STRUCTSTAT statbuffer;

  if(STAT(filename,&statbuffer)==0){
    return 1;
  }
  else{
    return 0;
  }
}

/* ------------------ which ------------------------ */

char *which(char *progname){
  /*! \fn char *which(char *progname)
      \brief returns the PATH directory containing the file progname
  */
  char *pathlistptr, fullpath[4096], pathlist[4096], prog[4096];
  char *dir,*returndir;
  const char *ext;
  char pathsep[2], dirsep[2];
  int lendir,lenprog;

#ifdef WIN32
  strcpy(pathsep,";");
  strcpy(dirsep,"\\");
#else
  strcpy(pathsep,":");
  strcpy(dirsep,"/");
#endif

  if(progname==NULL)return NULL;
  strcpy(prog,progname);
  progname=prog;

  pathlistptr=getenv("PATH");
  if(pathlistptr==NULL)return NULL;
  strcpy(pathlist,pathlistptr);
  
#ifdef WIN32
  lenprog=strlen(prog);
  ext=progname+lenprog-4;
  if(lenprog<=4||STRCMP(ext,".exe")!=0){
    strcat(progname,".exe");
 }
#endif
        
  dir=strtok(pathlist,pathsep);
  while(dir!=NULL){
    strcpy(fullpath,dir);
    strcat(fullpath,dirsep);
    strcat(fullpath,prog);
    if(file_exists(fullpath)==1){
      lendir=strlen(dir);
      if(lendir<=0)continue;
      NewMemory((void **)&returndir,(unsigned int)(lendir+2));
      strcpy(returndir,dir);
      strcat(returndir,dirsep);
#ifdef pp_BETA
      printf("Using %s in %s\n\n",prog,dir);
#endif
      return returndir;
    }
    dir=strtok(NULL,pathsep);
  }
#ifdef pp_BETA
  printf("%s not found in any path directory\n",prog);
#endif
  return NULL;
}

/* ------------------ MIN ------------------------ */

float MIN(float x,float y){
  /*! \fn float MIN(float x,float y)
      \brief returns the minimum of x and y
  */
  if(x<y){
    return x;
  }
  else{
    return y;
  }
}

/* ------------------ MAX ------------------------ */

float MAX(float x,float y){
  /*! \fn float MAX(float x,float y)
      \brief returns the maximum of x and y
  */
  if(x>y){
    return x;
  }
  else{
    return y;
  }
}

/* ------------------ frexp10 ------------------------ */

float frexp10(float x, int *exp10){
  float xabs, mantissa;

  xabs = fabs((double)x);
  if(x==0.0f){
    *exp10=0;
    return 0.0f;
  }
  mantissa = log10((double)xabs);
  *exp10 = (int)floor((double)mantissa);
      
  mantissa = pow((double)10.0f,(double)mantissa-(double)*exp10);
  if(x<0)mantissa = -mantissa;
  return mantissa;
}

/* ------------------ getstring ------------------------ */

char *getstring(char *buffer){
  /*! \fn *getstring(char *buffer)
      \brief return pointer to string contained between a pair of double quotes
  */
  char *begin,*end;
  int i;

  // if buffer contains msgid "string"
  // return a pointer to s in string

  begin=strchr(buffer,'"');
  if(begin==NULL)return NULL;
  begin++;
  end=strrchr(begin,'"');
  if(end==NULL)return NULL;
  end[0]=0;
  for(i=0;i<strlen(begin);i++){
    if(begin[i]!=' ')return begin;
  }
  return NULL;
}



