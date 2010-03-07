// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeheaders.h"
#include "smokeviewvars.h"

// svn revision character string
char string_util_revision[]="$Revision$";

// extern char dirseparator[];

/* ------------------ rootdir ------------------------ */

void getdir(char *argstart){
  size_t i,nargi;
  char *argi;

  nargi=strlen(argstart);
  argi = argstart+nargi-1;
  for(i=0;i<nargi;i++){
    if(strncmp(argi,dirseparator,1)==0){
      *(argi+1)='\0';
      return;
    }
    argi--;
  }
  *argstart='\0';
  return;
}


/* ------------------ lastname ------------------------ */

char *lastname(char *argi){
  size_t i,nargi;
  char *origargi;
  nargi=strlen(argi);
  origargi=argi;
  argi = argi+nargi-1;
  for(i=0;i<nargi;i++){
    if(strncmp(argi,dirseparator,1)==0){
      if(i!=0){return argi+1;}
      else{
        return NULL;
      }
    }
    argi--;
  }
  return origargi;
}

/* ------------------ trim ------------------------ */

void trim(char *line){
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
  size_t i,j,lenline;
  char *linecopy, c;
  char buffer[1024];
  size_t ibeg;
  size_t start=1;
  size_t bufstart=0;
  size_t lenbuf;

  lenline=strlen(line);
  linecopy=line;
  for(i=0;i<lenline;i++){
    c=*linecopy++;
    if(start==1){
      if(c==' ')continue;
      ibeg=i;
      start=0;
      buffer[0]=c;
      continue;
    }
    buffer[i-ibeg]=c;
    if(c==' '){
      buffer[i-ibeg]='\0';
      trimzeros(buffer);
      lenbuf = strlen(buffer);
      for(j=bufstart;j<bufstart+lenbuf;j++){
        line[j]=buffer[j-bufstart];
      }
      bufstart+=strlen(buffer);
      line[bufstart]=' ';
      bufstart++;
      start=1;
    }
  }
  line[bufstart]='\0';
}

/* ------------------ STRSTR ------------------------ */

char *STRSTR(char *c, const char *key){
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

  if(fabs((double)tval)<fabs((double)range)/100.0f)tval=0.0f;
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
  while (toupper(*s1) == toupper(*s2++)){
		if (*s1++ == 0)return (0);
  }
	return (toupper(*(const unsigned char *)s1) - toupper(*(const unsigned char *)(s2 - 1)));
}

/* ------------------ get_chid ------------------------ */

char *get_chid(char *file){
  FILE *stream;
  char *chidptr,*c;
  char buffer[1024];
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

/* ------------------ array2string ------------------------ */

void parse_string(char *string, char **tokens, int *ntokens){
  int i, len, in_quote, ntok, in_token, last_in_token;
  char *c;

  c=string;
  ntok=0;
  in_quote=0;
  in_token=0;
  last_in_token=0;
  len=strlen(string);
  for(i=0;i<=len;i++){
    switch (*c) {
      case '"':
        in_quote=1-in_quote;
        in_token=1;
        break;
      case ' ':
        if(in_quote==0){
          in_token=0;
        }
        break;
      case 0:
        in_token=0;
        break;
      default:
        in_token=1;
        break;
    }
    if(in_token>last_in_token){
      tokens[ntok]=c;
      ntok++;
    }
    if(in_token<last_in_token){
      *c=0;
    }
    last_in_token=in_token;
    c++;
  }
  *ntokens=ntok;
}

  /* ------------------ getfileinfo ------------------------ */

int fileexist(char *filename){
  STRUCTSTAT statbuffer;
  int statfile;
  char buffer[1024];

  statfile=STAT(filename,&statbuffer);
  return statfile;
}
/* ------------------ array2string ------------------------ */

char *which(char *prog, char *pathlist){
  char *pathlistptr;
  char fullpath[4096];
  char *dir;
  int lenpath;
  char pathsep[2];
  char dirsep[2];

#ifdef WIN32
  strcpy(pathsep,";");
  strcpy(dirsep,"\\");
#else
  strcpy(pathsep,":");
  strcpy(dirsep,"/";
#endif

  if(prog==NULL)return NULL;
  pathlistptr=getenv("PATH");
  if(pathlistptr==NULL)return NULL;
  lenpath=strlen(pathlistptr);
  strcpy(pathlist,pathlistptr);
        
  dir=strtok(pathlist,pathsep);
  while(dir!=NULL){
    printf("dir=%s\n",dir);
    strcpy(fullpath,dir);
    strcat(fullpath,dirsep);
    strcat(fullpath,prog);
    if(fileexist(fullpath)==0){
      return dir;
    }
    dir=strtok(NULL,pathsep);
  }
  return NULL;
}

