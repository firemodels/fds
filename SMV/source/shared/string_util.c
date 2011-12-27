// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char string_util_revision[]="$Revision$";

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
#include <dirent_win.h>
#else
#include <dirent.h>
#endif
#include "MALLOC.h"
#include "string_util.h"

/* ----------------------- fparsecsv ----------------------------- */

void fparsecsv(char *buffer, float *vals, int *valids, int ncols, int *ntokens){
  /*! \fn void fparsecsv(char *buffer, float *vals, int ncols, int *ntokens)
      \brief copy comma delimited values from buffer into floating point array vals
       returning number of values found in ntokens
  */
  int nt=0;
  char *token;

  token=strtok(buffer,",");
  while(token!=NULL&&nt<ncols){
    if(STRCMP(token,"NULL")==0){
      valids[nt]=0;
      vals[nt]=0.0;
    }
    else{
      valids[nt]=1;
      sscanf(token,"%f",&vals[nt]);
    }
    nt++;
    token=strtok(NULL,",");
    if(token!=NULL)trim(token);
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
  int i;
  int lenbuffer;

  lenbuffer=strlen(buffer);
  for(i=0;i<lenbuffer;i++){
    if(buffer[i]==',')buffer[i]=0;
  }
  tokens[nt++]=buffer;
  for(i=1;i<lenbuffer;i++){
    if(buffer[i]==0){
      tokens[nt++]=buffer+i+1;
    }
  }
  *ntokens=nt;
}

/* ------------------ getrowcols ------------------------ */

int getrowcols(FILE *stream, int *nrows, int *ncols){
  /*! \fn int getrowcols(FILE *stream, int *nrows, int *ncols)
      \brief find number of rows (nrows) and number of columns (ncols) in
       comma delimited file pointed to by stream and returns the length
       of the longest line
  */
  int nnrows=0,nncols=1,maxcols=0,linelength=0,maxlinelength=0;

  while(!feof(stream)){
    char ch;

    ch = getc(stream);
    linelength++;
    if(ch == ',')nncols++;
    if(ch=='\n'){
      if(linelength>maxlinelength)maxlinelength=linelength;
      if(nncols>maxcols)maxcols=nncols;
      linelength=0;
      nncols=1;
      nnrows++;
    }
  }
  *nrows=nnrows;
  *ncols=maxcols;
  rewind(stream);
  return maxlinelength;
}

/* ------------------ getrevision ------------------------ */

int getrevision(char *svn){
  char svn_string[256];
  char *svn_ptr;
  int return_val;

  svn_ptr=svn_string;
  svn=strchr(svn,':');
  if(svn==NULL||strlen(svn)<=4)return 0;
  
  svn++;
  strcpy(svn_ptr,svn);
  svn_ptr=trim_front(svn_ptr);
  svn_ptr[strlen(svn_ptr)-1]=0;
  trim(svn_ptr);
  sscanf(svn_ptr,"%i",&return_val);
  return return_val;
}

/* ------------------ stripquotes ------------------------ */

void stripquotes(char *buffer){
  /*! \fn void stripquotes(char *buffer)
      \brief replaces quotes (") with blanks in the character string buffer
  */
  char *c;

  for(c=buffer;c<buffer+strlen(buffer);c++){
    if(*c=='"')*c=' ';
  }
}
/* ------------------ stripcommas ------------------------ */

void stripcommas(char *buffer){
  /*! \fn void stripcommas(char *buffer)
      \brief replaces commas (,) with blanks in the character string buffer
  */
  char *c;

  for(c=buffer;c<buffer+strlen(buffer);c++){
    if(*c==',')*c=' ';
  }
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

/* ------------------ trim_commas ------------------------ */

void trim_commas(char *line){
  char *c;
  
  for(c = line + strlen(line) - 1;c>=line;c--){
    if(isspace(*c))continue;
    if(strncmp(c,",",1)!=0)break;
    *c=' ';
  }
}

/* ------------------ trim ------------------------ */

void trim(char *line){
  /*! \fn void trim(char *line)
      \brief removes trailing white space from the character string line
  */
  char *c;
  size_t len;

  if(line==NULL)return;
  len = strlen(line);
  if(len==0)return;
  for(c=line+len-1; c>=line; c--){
    if(isspace(*c))continue;
    *(c+1)='\0';
    return;
  }
  *line='\0';
}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  /*! \fn char *trim_front(char *line)
      \brief returns a pointer to the first non-blank character in the character string line
  */
  char *c;

  for(c=line;c<=line+strlen(line)-1;c++){
    if(!isspace(*c))return c;
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

/* ------------------ STRSTR ------------------------ */

char *STRSTR(char *string, const char *key){
  char *k,*s,*ss;

  if(string==NULL||key==NULL)return NULL;

  for(s=string;*s!=0;s++){
    for(k=(char *)key;*k!=0;k++){
      ss = s + (k-key);
      if(*ss==0)return NULL;
      if(toupper(*ss)!=toupper(*k))break;
    }
    if(*k==0)return s;
  }
  return NULL;
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

/* ------------------ trim_string ------------------------ */

char *trim_string(char *buffer){
  /*! \fn char *trim_string(char *buffer)
      \brief removes trailing blanks from buffer and returns a pointer to the first non-blank character
  */
  int len;
  char *bufptr;

  if(buffer==NULL)return NULL;
  len=strlen(buffer);
  buffer[len-1]='\0';
  bufptr=trim_front(buffer);
  trim(bufptr);
  return bufptr;
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

  /* ------------------ time2timelabel ------------------------ */

char *time2timelabel(float time, float dt, char *timelabel){
  char *timelabelptr;

  if(dt<0.001){
    sprintf(timelabel,"%4.4f",time);
  }
  else if(dt>=0.001&&dt<0.01){
    sprintf(timelabel,"%4.3f",time);
  }
  else if(dt>=0.01&&dt<0.1){
    sprintf(timelabel,"%4.2f",time);
  }
  else{
    sprintf(timelabel,"%4.1f",time);
  }
  trimzeros(timelabel);
  trim(timelabel);
  timelabelptr=trim_front(timelabel);
  return timelabelptr;
}

/* ------------------ match_upper ------------------------ */

int match_upper(char *buffer, const char *key){
  size_t lenbuffer;
  size_t lenkey;
  int i;

  lenkey=strlen(key);
  trim(buffer);
  lenbuffer=strlen(buffer);

  if(lenbuffer<lenkey)return 0;
  for(i=0;i<lenkey;i++){
    if(toupper(buffer[i])!=toupper(key[i]))return 0;
  }
  if(lenbuffer>lenkey&&buffer[lenkey]==':')return 2;
  if(lenbuffer>lenkey&&!isspace(buffer[lenkey]))return 0;
  return 1;
}

/* ------------------ match ------------------------ */

int match(char *buffer, const char *key){
  size_t lenbuffer;
  size_t lenkey;

  lenkey=strlen(key);
  lenbuffer=strlen(buffer);
  if(lenbuffer<lenkey)return 0;
  if(strncmp(buffer,key,lenkey) != 0)return 0;
  if(lenbuffer>lenkey&&!isspace(buffer[lenkey]))return 0;
  return 1;
}

/* ----------------------- remove_comment ----------------------------- */

void remove_comment(char *buffer){
  char *comment;
  comment = strstr(buffer,"//");
  if(comment!=NULL)comment[0]=0;
  return;
}
