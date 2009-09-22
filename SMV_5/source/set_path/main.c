// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <tchar.h>
#include <math.h>
#include <ctype.h>
#include "svn_revision.h"

char main_revision[]="$Revision$";

#define BUFFER_SIZE 10050

int add_path=0;
int remove_path=0;
int display_path=0;
int act_on_user_path=1;
int act_on_system_path=0;
int prompt_user=0;
int test_mode=1;
char path_type[10];


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
  char *newentry=NULL;
  char pathbuffer[BUFFER_SIZE], *pathbuffer_ptr;
  char tokens[BUFFER_SIZE], newpath[BUFFER_SIZE], *token;
  int newentry_present;

#define ADD_USER_PATH 0
#define REMOVE_USER_PATH 1

#define ADD_SYSTEM_PATH 0
#define REMOVE_SYSTEM_PATH 1

  pathbuffer_ptr=pathbuffer;
  if(argc==1){
    usage();
    return 1;
  }
 
  for(i=1;i<argc;i++){
    arg=argv[i];
    if(arg[0]!='-'||strlen(arg)<=1)break;

    switch(arg[1]){
      case 'a':
        add_path=1;
        remove_path=0;
        newentry=argv[i+1];
        i++;
        break;
      case 'd':
        display_path=1;
        break;
      case 'r':
        remove_path=1;
        add_path=0;
        newentry=argv[i+1];
        i++;
        break;
      case 'u':
        act_on_user_path=1;
        act_on_system_path=0;
        strcpy(path_type,"User");
        break;
      case 's':
        act_on_user_path=0;
        act_on_system_path=1;
        strcpy(path_type,"System");
        break;
      case 'p':
        prompt_user=1;
        break;
      case 't':
        test_mode=1;
        break;
      case 'v':
        version();
        return 0;
        break;
      default:
        usage();
        return 1;
    }
  }

  if(test_mode==1)display_path=1;

  // get the path (user or system)

  if(act_on_user_path==1){
    if(reg_path(REG_GET,REG_USER_PATH,pathbuffer)==0)return 1;
    if(display_path==1){
      if(strlen(pathbuffer)>0){
        printf("User path: %s\n\n",pathbuffer);
      }
    }
  }
  else if(act_on_system_path==1){
    if(reg_path(REG_GET,REG_SYSTEM_PATH,pathbuffer)==0)return 1;
    if(display_path==1){
      if(strlen(pathbuffer)>0){
        printf("System path: %s\n\n",pathbuffer);
      }
    }
  }
  else{
    usage();
    return 0;
  }

  if(add_path==1&&newentry!=NULL){
    strcpy(tokens,pathbuffer);
    token=strtok(tokens,";");
    newentry_present=0;
    while(token!=NULL){
      if(STRCMP(token,newentry)==0){
        newentry_present=1;
        break;
      }
      token=strtok(NULL,";");
    }
    if(newentry_present==0){
      int answer=0;
      char c_answer[10], *c_answer_ptr;

      if(strlen(pathbuffer)>0)strcat(pathbuffer,";");
      strcat(pathbuffer,newentry);
      if(prompt_user==1){
        printf(" Set %s path to:\n",path_type);
        printf("%s\n",pathbuffer);
        printf("y=yes, n=no\n");
        scanf("%s",c_answer);
        c_answer_ptr=trim_front(c_answer);
        if(c_answer_ptr!=NULL&&strlen(c_answer_ptr)>0&&toupper(c_answer_ptr[0])=='Y')answer=1;
      }
      if(prompt_user==0||answer==1){
        if(act_on_user_path==1){
          if(reg_path(REG_SET,REG_USER_PATH,pathbuffer)==0)return 1;
        }
        else{
          if(reg_path(REG_SET,REG_SYSTEM_PATH,pathbuffer)==0)return 1;
        }
      }
      if(display_path==1){
        printf("%s path: %s\n",path_type,pathbuffer);
      }
    }
  }

  if(remove_path==1&&newentry!=NULL){
    strcpy(tokens,pathbuffer);
    token=strtok(tokens,";");
    newentry_present=0;
    strcpy(newpath,"");
    while(token!=NULL){
      if(STRSTR(token,newentry)!=NULL){
        newentry_present=1;
        token=strtok(NULL,";");
        continue;
      }
      strcat(newpath,token);
      strcat(newpath,";");
      token=strtok(NULL,";");
    }
    if(newentry_present==1){
      int answer=0;
      char c_answer[10], *c_answer_ptr;

      if(prompt_user==1){
        printf(" Set %s path to:\n",path_type);
        printf("%s\n",pathbuffer);
        printf("y=yes, n=no\n");
        scanf("%s",c_answer);
        c_answer_ptr=trim_front(c_answer);
        if(c_answer_ptr!=NULL&&strlen(c_answer_ptr)>0&&toupper(c_answer_ptr[0])=='Y')answer=1;
      }
      if(prompt_user==0||answer==1){
        if(act_on_user_path==1){
          if(reg_path(REG_SET,REG_USER_PATH,newpath)==0)return 1;
        }
        else{
          if(reg_path(REG_SET,REG_SYSTEM_PATH,newpath)==0)return 1;
        }
      }
    }
    if(display_path==1){
      printf("%s path: %s\n",path_type,newpath);
    }
  }
  return 0;
}

/* ------------------ reg_path ------------------------ */

int reg_path(int setget, int pathtype, char *path){

  HKEY hKey, hTree;
  long lRet;
  char temp[10000];
  DWORD dwBufLen;
  int lenpath;

  LPCTSTR reg_path;
  
  char creg_user_path[]="Environment";
  LPCTSTR reg_user_path=creg_user_path;
  
  char creg_system_path[]="SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment";
  LPCTSTR reg_system_path=(LPCTSTR)creg_system_path;

  char cPATH[]="Path";
  LPCTSTR PATH=(LPCTSTR)cPATH;

  switch (pathtype) {
    case REG_USER_PATH:
      reg_path=reg_user_path;
      hTree=HKEY_CURRENT_USER;
      break;
    case REG_SYSTEM_PATH:
      reg_path=reg_system_path;
      hTree=HKEY_LOCAL_MACHINE;
      break;
  }
  switch (setget) {
    case REG_GET:
      lRet = RegOpenKeyEx( hTree, reg_path, 0, KEY_QUERY_VALUE, &hKey );
      if(lRet!=ERROR_SUCCESS){
        printf("RegOpenKeyEx error: %i\n",(int)lRet);
        return 0;
      }
      dwBufLen=sizeof(temp);
      lRet = RegQueryValueEx( hKey, PATH, NULL, NULL, (BYTE*)&temp, &dwBufLen );
      switch (lRet){
        case ERROR_FILE_NOT_FOUND:
          printf("%s path not present\n",path_type);
          strcpy(path,"");
          dwBufLen=0;
          break;
        case ERROR_SUCCESS:
          strncpy(path,temp,dwBufLen);
          path[dwBufLen]=0;
          break;
        default:
          printf("RegQueryValueEx error: %i\n",(int)lRet);
          return 0;
      }
      lRet = RegCloseKey( hKey);
      if(lRet!=ERROR_SUCCESS){
        printf("RegCloseKey error: %i\n",(int)lRet);
        return 0;
      }
      break;
    case REG_SET:
      lRet = RegOpenKeyEx( hTree, reg_path, 0, KEY_QUERY_VALUE | KEY_SET_VALUE, &hKey );
      if(lRet!=ERROR_SUCCESS){
        printf("RegOpenKeyEx error: %i\n",(int)lRet);
        return 0;
      }
      lenpath=strlen(path);
      if(lenpath>0){
        if(path[lenpath-1]==';'){
          path[lenpath-1]='\0';
        }
      }
      if(display_path==1){
        printf("Setting %s path to: %s\n",path_type,path);
      }
      if(test_mode==0){
        lRet = RegSetValueEx(hKey,PATH,0,REG_EXPAND_SZ,(LPBYTE)path,strlen(path)+1);
      }
      if(lRet!=ERROR_SUCCESS){
        printf("RegSetValueEx error: %i\n",(int)lRet);
        return 0;
      }
      lRet = RegCloseKey( hKey);
      if(lRet!=ERROR_SUCCESS){
        printf("RegCloseKey error: %i\n",(int)lRet);
        return 0;
      }
      break;
  }
  return 1;
}

/* ------------------ version ------------------------ */

void version(void){
    int svn_num;

    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("set_path %s - %s\n\n",VERSION,__DATE__);
    printf("Version: %s\n",VERSION);
    printf("Revision Number: %i\n",svn_num);
    printf("Build Date: %s\n",__DATE__);
}

/* ------------------ usage ------------------------ */

void usage(void){
  int max_revision;

  max_revision = getmaxrevision();

  printf("set_path SVN:revision:%i\n",max_revision);
  printf("  Modify or display the User or System path environmental variables.\n\n");
  printf("Usage:\n\n");
  printf("  set_path [-s][-u] [-a path_entry] [-r path_entry] [-d][-p][-v]\n\n");
  printf("where\n\n");
  printf("  -a path_entry - add path_entry to the User or System path\n");
  printf("  -r path_entry - remove any entry from the Path containing path_entry\n");
  printf("  -d - display path \n");
  printf("  -p - prompt user before making changes\n");
  printf("  -s - perform action on the System path\n");
  printf("  -u - perform action on the User path\n");
  printf("  -t - test mode, show but do not change Path variables\n");
  printf("  -v - show versioning information\n");
}
