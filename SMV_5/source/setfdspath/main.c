// $Date: 2009-03-01 08:51:14 -0500 (Sun, 01 Mar 2009) $ 
// $Revision: 3426 $
// $Author: gforney $

#include "options.h"
#define INMAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "MALLOC.h"
#include "ASSERT.h"


char *trim_front(char *line);
void usage(char *program_name);
void trim(char *line);
int parse_path_key(int flag, char *newentry);
int STRCMP(const char *s1, const char *s2);
char *STRSTR(char *c, const char *key);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *arg;
  int i;
  int clean_old_fds=0;
  char *newentry,*program_name;
  
  program_name=argv[0];
  newentry=argv[1];

  if(argc==1){
    usage(program_name);
    return 1;
  }
 
  for(i=1;i<argc;i++){
    arg=argv[i];
    if(arg[0]!='-'||strlen(arg)<=1)break;

    switch(arg[1]){
      case 'c':
        clean_old_fds=1;
        break;
    }
  }

  // add newentry to end of local path

  system("reg query hkey_current_user\\Environment /v Path > local_path.txt"); 

  if(parse_path_key(0,newentry)==1)return 1;
  if(clean_old_fds==0)return 0;

  system("reg query \"hkey_local_machine\\SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment\" /v Path > local_path.txt"); 
  if(parse_path_key(1,NULL)==1)return 1;

  return 0;
}

/* ------------------ parse_path_key ------------------------ */

int parse_path_key(int flag, char *newentry){
  char buffer[1024],tokens[1024], command[1024];
  FILE *stream=NULL;
  char *fullpath,*token;

  stream=fopen("local_path.txt","r");
  if(stream==NULL)return 1;

  if(fgets(buffer,1024,stream)==NULL){
    fclose(stream);
    return 1;
  }

  while(STRSTR(buffer,"Path")==NULL){
    if(fgets(buffer,1024,stream)==NULL){
      fclose(stream);
      return 1;
    }
  }

  fullpath=STRSTR(buffer,"REG_EXPAND_SZ");
  if(fullpath==NULL){
    fclose(stream);
    return 1;
  }
  trim(fullpath);
  fullpath+=sizeof("REG_EXPAND_SZ");
  fullpath=trim_front(fullpath);
  strcpy(tokens,fullpath);
  if(flag==0&&newentry!=NULL){
    token=strtok(tokens,";");
    while(token!=NULL){
      if(newentry!=NULL&&STRCMP(token,newentry)==0){
        stream=NULL;
        return 1;
      }
      token=strtok(NULL,";");
    }
    strcat(fullpath,";");
    strcat(fullpath,newentry);
    strcpy(command,"reg add hkey_current_user\\Environment /v Path /t REG_EXPAND_SZ /f /d "); 
    strcat(command,"\"");
    strcat(command,fullpath);
    strcat(command,"\"");
    system(command);
  }
  if(flag==1&&newentry==NULL){
    token=strtok(tokens,";");
    strcpy(fullpath,"");
    while(token!=NULL){
      if(STRSTR(token,"NIST")!=NULL&&STRSTR(token,"FDS")!=NULL)continue;
      strcat(fullpath,";");
      strcat(fullpath,token);
      token=strtok(NULL,";");
    }
    strcpy(command,"reg add hkey_current_user\\Environment /v Path /t REG_EXPAND_SZ /f /d "); 
    strcat(command,"\"");
    strcat(command,fullpath);
    strcat(command,"\"");
    system(command);
  }

  unlink("local_path.txt");
  fclose(stream);
  return 0;
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
  char *tab="\t";
  const char *c;
  size_t i,len;

  c = line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(strncmp(c,blank,1)!=0&&strncmp(c,tab,1)!=0)return line+i;
    c++;
  }
  return line;
}

/* ------------------ usage ------------------------ */

void usage(char *program_name){
  printf("%s\n",program_name);
  printf("  Add the FDS bin directory to the system path and/or remove\n");
  printf("  previus FDS bin directories from the system path\n");
  printf("Usage:\n\n");
  printf("  %s [-r] path_entry\n\n",program_name);
  printf("  path_entry - location of FDS and Smokeview executables\n");
  printf("  -r  - remove pre FDS 5.4 directories from the system path\n");
}

/* ------------------ STRCMP ------------------------ */

int STRCMP(const char *s1, const char *s2){
  while (toupper(*s1) == toupper(*s2++)){
		if (*s1++ == 0)return (0);
  }
	return (toupper(*(const unsigned char *)s1) - toupper(*(const unsigned char *)(s2 - 1)));
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
