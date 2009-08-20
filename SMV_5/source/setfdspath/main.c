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

int add_user_path=0;
int remove_user_path=0;
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
#define REMOVE_USER_PATH 3

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
        if(remove_user_path==0){
          add_user_path=1;
          newentry=argv[i+1];
        }
        i++;
        break;
      case 'c':
        remove_user_path=1;
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

  // add or show user path

  if((newentry!=NULL&&add_user_path==1)||show_user_path==1){
    strcpy(command,"reg query hkey_current_user\\Environment /v Path > local_path.txt"); 
    if(show_debug==1){
      if(add_user_path==1)printf("*** Querying user path (add_user_path=1)\n");
      printf("executing: %s\n\n",command);
    }
    system(command);
    if(newentry!=NULL&&add_user_path==1){
      path=parse_path_key(ADD_USER_PATH,pathbuffer,newentry);
    }
    else{
      path=parse_path_key(ADD_USER_PATH,pathbuffer,NULL);
    }
    _unlink("local_path.txt");
    if(path!=NULL&&(show_user_path==1||show_debug==1)){
      printf("User path: %s\n\n",path);
    }
  }

  // remove user path

  if(newentry!=NULL&&remove_user_path==1){
    strcpy(command,"reg query hkey_current_user\\Environment /v Path > local_path.txt"); 
    if(show_debug==1){
      printf("*** Querying user path (remove_user_path=1)\n");
      printf("executing: %s\n\n",command);
    }
    system(command);
    path=parse_path_key(REMOVE_USER_PATH,pathbuffer,newentry);
    _unlink("local_path.txt");
  }

  // system path

  if(clean_old_fds==1||show_system_path==1){
    strcpy(command,"reg query \"hkey_local_machine\\SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment\" /v Path > local_path.txt"); 
    if(show_debug==1){
      printf("*** Querying system path\n");
      printf("executing: %s\n\n",command);
    }
    system(command);
    path=parse_path_key(CLEAN_SYSTEM_PATH,pathbuffer,NULL);
    _unlink("local_path.txt");
    if(path!=NULL&&(show_system_path==1||show_debug==1)){
      printf("System path: %s\n\n",path);
    }
  }

  return 0;
}

/* ------------------ add_percen ------------------------ */

int add_percen(char *fullpath){
  char fullpath_copy[BUFFER_SIZE];
  size_t i,i2;
  int found=0;

  if(fullpath==NULL||strlen(fullpath)==0)return 1;
  strcpy(fullpath_copy,fullpath);
  i2=0;
  for(i=0;i<strlen(fullpath_copy);i++){
    if(fullpath_copy[i]=='%'){
      fullpath[i2++]='%';
      found=1;
    }
    fullpath[i2++]=fullpath_copy[i];
  }
  fullpath[i2]=0;
  return found;
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
  if(stream==NULL){
    printf("*** Error: unable to open local_path.txt to retrieve path info\n");
    return NULL;
  }

  if(fgets(buffer,BUFFER_SIZE,stream)==NULL){
    printf("*** Error: local_path.txt is empty\n");
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
      if(show_debug==1){
        printf("ADD_USER_PATH\n");
      }
      token=strtok(tokens,";");
      newentry_present=0;
      while(token!=NULL){
        if(newentry!=NULL&&STRCMP(token,newentry)==0)newentry_present=1;
        token=strtok(NULL,";");
      }
      // don't add a path entry if the entry is already in the user path
      if(newentry!=NULL&&newentry_present==0){
        int percen;
        FILE *streamcom;

        if(show_debug==1){
          printf("%s not found in User Path - so add it\n",newentry);
        }
        strcat(fullpath,";");
        strcat(fullpath,newentry);
        strcpy(command,"reg add hkey_current_user\\Environment /v Path /t ");
        percen=add_percen(fullpath);
        if(offset_type==1||percen==1){
          strcat(command,"REG_EXPAND_SZ /f /d "); 
        }
        if(offset_type==2){
          strcat(command,"REG_SZ /f /d "); 
        }
        strcat(command,"\"");
        strcat(command,fullpath);
        strcat(command,"\"");
        if(show_debug==1)printf("executing: %s\n\n",command);
        streamcom=fopen("setpath_addnew.bat","w");
        if(streamcom!=NULL){
          fprintf(streamcom,"@echo off\n");
          fprintf(streamcom,"%s",command);
          fclose(streamcom);
          system("setpath_addnew.bat");
          _unlink("setpath_addnew.bat");
        }
        printf("  %s, was added to the User Path\n",newentry);
        printf("A re-boot is required after this installation completes.\n");
      }
      else{
        if(show_debug==1){
          if(newentry!=NULL){
            printf("%s was found in the User Path - so will not be added\n");
          }
          else{
            printf("*** error the path entry variable is NULL\n");
          }
        }
      }
      break;
    case REMOVE_USER_PATH:
      if(show_debug==1){
        printf("REMOVE_USER_PATH\n");
      }
      if(remove_user_path==1){
        int old_fds_found=0;

        token=strtok(tokens,";");
        strcpy(fullpath,"");
        while(token!=NULL){
          if(STRSTR(token,newentry)!=NULL){
            token=strtok(NULL,";");
            old_fds_found=1;
            continue;
          }
          strcat(fullpath,token);
          token=strtok(NULL,";");
          if(token!=NULL)strcat(fullpath,";");
        }
        trim(fullpath);
        if(old_fds_found==1){
          int lenstr, percen;
          FILE *streamcom;

          if(show_debug==1){
            printf("old FDS path found\n");
          }
          lenstr = strlen(fullpath);
          if(fullpath[lenstr-1]==';'){
            fullpath[lenstr-1]=0;
          }
          strcpy(command,"reg add hkey_current_user\\Environment /v Path /t ");
          percen=add_percen(fullpath);
          if(offset_type==1||percen==1){
            strcat(command,"REG_EXPAND_SZ /f /d "); 
          }
          if(offset_type==2){
            strcat(command,"REG_SZ /f /d "); 
          }
          strcat(command,"\"");
          strcat(command,fullpath);
          strcat(command,"\"");
          if(show_debug==1)printf("executing: %s\n\n",command);
          streamcom=fopen("setpath_removenew.bat","w");
          if(streamcom!=NULL){
            fprintf(streamcom,"@echo off\n");
            fprintf(streamcom,"%s",command);
            fclose(streamcom);
            system("setpath_removenew.bat");
            _unlink("setpath_removenew.bat");
          }
          else{
            printf("*** Error unable to create the file setpath_removenew.bat\n");
          }
          printf("  %s, was removed from the User Path\n",newentry);
          printf("A re-boot is required to complete the installation.\n");
        }
        if(old_fds_found==0){
          printf("  %s not found in the User Path,\n",newentry);
        }
      }
      break;
    case CLEAN_SYSTEM_PATH:
      if(clean_old_fds==1){
        int old_fds_found=0;

      if(show_debug==1){
        printf("CLEAN_SYSTEM_PATH\n");
      }
        token=strtok(tokens,";");
        strcpy(fullpath,"");
        while(token!=NULL){
          if(STRSTR(token,"NIST")!=NULL&&(
            STRSTR(token,"fds")!=NULL||
            STRSTR(token,"utilities")!=NULL||
            STRSTR(token,"smokeview")!=NULL)
            ){
            token=strtok(NULL,";");
            old_fds_found=1;
            continue;
          }
          strcat(fullpath,token);
          token=strtok(NULL,";");
          if(token!=NULL)strcat(fullpath,";");
        }
        trim(fullpath);
        if(old_fds_found==1){
          int lenstr;
          FILE *streamcom;

          if(show_debug==1){
            printf("Old FDS path enty was found in the system path\n");
          }
          lenstr = strlen(fullpath);
          if(fullpath[lenstr-1]==';'){
            fullpath[lenstr-1]=0;
          }
          strcpy(command,"reg add \"hkey_local_machine\\SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment\" /v Path /t REG_EXPAND_SZ /f /d "); 
          strcat(command,"\"");
          add_percen(fullpath);
          strcat(command,fullpath);
          strcat(command,"\"");
          printf("  Found and being removed.\n");
          printf("A re-boot is required after this installation completes.\n");
          if(show_debug==1)printf("executing: %s\n\n",command);
          streamcom=fopen("setpath_removeold.bat","w");
          if(streamcom!=NULL){
            fprintf(streamcom,"@echo off\n");
            fprintf(streamcom,"%s",command);
            fclose(streamcom);
            system("setpath_removeold.bat");
            _unlink("setpath_removeold.bat");
          }
          else{
            printf("***error unable to create the file setpath_removeold.bat\n");
          }
        }
        else{
          printf("  None found.\n");
        }
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
  printf("  set_path [-a path_entry][-c][-d][-r][-s][-u][-v]\n\n");
  printf("where\n\n");
  printf("  -a path_entry - add the directory, path_entry, to the user path\n");
  printf("  -c - cancel or remove path entry defined by current working directory\n");
  printf("  -d - turn on debug printing\n");
  printf("  -r - remove pre FDS 5.4 directories from the system path\n");
  printf("  -s - show the system path\n");
  printf("  -u - show the user path\n");
  printf("  -v - show versioning information\n");
}
