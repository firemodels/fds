// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "MALLOC.h"
#include "ASSERT.h"
#include "svn_revision.h"

char main_revision[]="$Revision$";

#define BUFFER_SIZE 2050

int clean_old_fds=0;
int show_user_path=0;
int show_system_path=0;
int show_debug=0;

char *trim_front(char *line);
void usage(void);
void version(void);
void trim(char *line);
char *parse_path_key(int flag, char *path_buffer, char *newentry);
int STRCMP(const char *s1, const char *s2);
char *STRSTR(char *c, const char *key);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *arg;
  int i;
  char *newentry=NULL,*program_name;
  char pathbuffer[BUFFER_SIZE], command[BUFFER_SIZE+100];
  char *path;

#define ADD_USER_PATH 0
#define CLEAN_SYSTEM_PATH 1
#define DELETE_FILE 2

  initMM();
  program_name=argv[0];

  if(argc==1){
    usage();
    return 1;
  }
 
  for(i=1;i<argc;i++){
    arg=argv[i];
    if(arg[0]!='-'||strlen(arg)<=1)break;

    switch(arg[1]){
      case 'a':
        newentry=argv[i+1];
        i++;
        break;
      case 'd':
        show_debug=1;
        break;
      case 'r':
        clean_old_fds=1;
        break;
      case 'u':
        show_user_path=1;
        break;
      case 'v':
        version();
        return 0;
        break;
      case 's':
        show_system_path=1;
        break;
      default:
        usage();
        return 1;
    }
  }

  // user path

  if(newentry!=NULL||show_user_path==1){
    strcpy(command,"reg query hkey_current_user\\Environment /v Path > local_path.txt"); 
    if(show_debug==1)printf("executing: %s\n\n",command);
    system(command);
    path=parse_path_key(ADD_USER_PATH,pathbuffer,newentry);
    _unlink("local_path.txt");
    if(path!=NULL&&(show_user_path==1||show_debug==1)){
      printf("User path: %s\n\n",path);
    }
  }

  // system path

  if(clean_old_fds==1||show_system_path==1){
    strcpy(command,"reg query \"hkey_local_machine\\SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment\" /v Path > local_path.txt"); 
    if(show_debug==1)printf("executing: %s\n\n",command);
    system(command);
    path=parse_path_key(CLEAN_SYSTEM_PATH,pathbuffer,NULL);
    _unlink("local_path.txt");
    if(path!=NULL&&(show_system_path==1||show_debug==1)){
      printf("System path: %s\n\n",path);
    }
  }

  return 0;
}

/* ------------------ parse_path_key ------------------------ */

char *parse_path_key(int flag, char *buffer, char *newentry){
  char tokens[BUFFER_SIZE], command[BUFFER_SIZE+100];
  FILE *stream=NULL;
  char *fullpath=NULL,*token;
  int newentry_present;
  int offset=0;
  int offset_type=0;

  stream=fopen("local_path.txt","r");
  if(stream==NULL)return NULL;

  if(fgets(buffer,BUFFER_SIZE,stream)==NULL){
    fclose(stream);
    return NULL;
  }

  while(STRSTR(buffer,"Path")==NULL){
    if(fgets(buffer,BUFFER_SIZE,stream)==NULL){
      fclose(stream);
      return NULL;
    }
  }

  fullpath=STRSTR(buffer,"REG_EXPAND_SZ");
  offset=sizeof("REG_EXPAND_SZ");
  if(fullpath!=NULL)offset_type=1;
  if(fullpath==NULL){
    fullpath=STRSTR(buffer,"REG_SZ");
    if(fullpath!=NULL)offset_type=2;
    offset=sizeof("REG_SZ");
  }
  if(fullpath==NULL){
    fclose(stream);
    return NULL;
  }
  trim(fullpath);
  fullpath+=offset;
  fullpath=trim_front(fullpath);
  strcpy(tokens,fullpath);
  switch (flag){
    case ADD_USER_PATH:
      token=strtok(tokens,";");
      newentry_present=0;
      while(token!=NULL){
        if(newentry!=NULL&&STRCMP(token,newentry)==0)newentry_present=1;
        token=strtok(NULL,";");
      }
      // don't add a path entry if the entry is already in the user path
      if(newentry!=NULL&&newentry_present==0){
        strcat(fullpath,";");
        strcat(fullpath,newentry);
        strcpy(command,"reg add hkey_current_user\\Environment /v Path /t ");
        if(offset_type==1){
          strcat(command,"REG_EXPAND_SZ /f /d "); 
        }
        if(offset_type==2){
          strcat(command,"REG_SZ /f /d "); 
        }
        strcat(command,"\"");
        strcat(command,fullpath);
        strcat(command,"\"");
        if(show_debug==1)printf("executing: %s\n\n",command);
        system(command);
      }
      break;
    case CLEAN_SYSTEM_PATH:
      if(clean_old_fds==1){
        token=strtok(tokens,";");
        strcpy(fullpath,"");
        while(token!=NULL){
          if(STRSTR(token,"NIST")!=NULL&&(
            STRSTR(token,"fds")!=NULL||
            STRSTR(token,"utilities")!=NULL||
            STRSTR(token,"smokeview")!=NULL)
            ){
            token=strtok(NULL,";");
            continue;
          }
          strcat(fullpath,token);
          token=strtok(NULL,";");
          if(token!=NULL)strcat(fullpath,";");
        }
        trim(fullpath);
        {
          int lenstr;

          lenstr = strlen(fullpath);
          if(fullpath[lenstr-1]==';'){
            fullpath[lenstr-1]=0;
          }
        }
        strcpy(command,"reg add \"hkey_local_machine\\SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment\" /v Path /t REG_EXPAND_SZ /f /d "); 
        strcat(command,"\"");
        strcat(command,fullpath);
        strcat(command,"\"");
        if(show_debug==1)printf("executing: %s\n\n",command);
        system(command);
      }
      break;
  }

  fclose(stream);
  return fullpath;
}

/* ------------------ version ------------------------ */

void version(void){
    int svn_num;

    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("set_path %s - %s\n\n",VERSION,__DATE__);
    printf("Version: %s\n",VERSION);
    printf("Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
}

/* ------------------ usage ------------------------ */

void usage(void){
  int max_revision;

  max_revision = getmaxrevision();

  printf("set_path SVN:revision:%i\n",max_revision);
  printf("  Add the FDS bin directory to the user path and/or remove previous FDS bin\n");
  printf("  directories from the system path.  All parameters are optional.\n\n");
  printf("Usage:\n\n");
  printf("  set_path [-a path_entry][-d][-r][-s][-u][-v]\n\n");
  printf("where\n\n");
  printf("  -a path_entry - add the directory, path_entry, to the user path\n");
  printf("  -d - turn on debug printing\n");
  printf("  -r - remove pre FDS 5.4 directories from the system path\n");
  printf("  -s - show the system path\n");
  printf("  -u - show the user path\n");
  printf("  -v - show versioning information\n");
}
