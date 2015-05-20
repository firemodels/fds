#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <ctype.h>
#include "ASSERT.h"
#include "string_util.h"
#include "datadefs.h"

// dummy change to bump version to 1.0

#define BUFFER_SIZE 10050

int add_path=0;
int at_end=1;
int remove_path=0;
int display_path=0;
int act_on_user_path=1;
int act_on_system_path=0;
int prompt_user_flag=0;
int test_mode=0;
int path_summary=0;
int batch_mode=0;

char path_type[10];


void usage(void);
void version(void);

/* ------------------ backup_path ------------------------ */

void backup_path(char *path_type_local, char *pathbuffer){
  FILE *stream;
  char file[256],filebase[256];
  int i;

  strcpy(filebase,path_type_local);
  strcat(filebase,"_path_backup");
  for(i=0;i<=100;i++){
    if(i==0){
      strcpy(file,filebase);
      strcat(file,".txt");
    }
    else{
      sprintf(file,"%s_%03i.txt",filebase,i);
    }
    stream=fopen(file,"r");
    if(stream==NULL)break;
    fclose(stream);
  }

  stream=fopen(file,"w");
  if(stream!=NULL){
    fprintf(stream,"%s Path\n",path_type_local);
    fprintf(stream,"%s\n",pathbuffer);
    fclose(stream);
  }
}

/* ------------------ prompt_user ------------------------ */

int prompt_user(char *path_type_local, char *pathbuffer){
  int answer=0;
  char c_answer[10], *c_answer_ptr;

  c_answer_ptr=c_answer;
  printf("\nSet %s path to:\n",path_type_local);
  printf("%s ?\n",pathbuffer);
  printf("y=yes, n=no\n");
  scanf("%s",c_answer);
  c_answer_ptr=trim_front(c_answer);
  if(c_answer_ptr!=NULL&&strlen(c_answer_ptr)>0&&toupper(c_answer_ptr[0])=='Y')answer=1;
  return answer;
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *arg;
  int i;
  char *newentry=NULL;
  char pathbuffer[BUFFER_SIZE];
  char pathbuffer2[BUFFER_SIZE];
  char tokens[BUFFER_SIZE], newpath[BUFFER_SIZE], *token;
  int newentry_present;

  strcpy(path_type,"User");
  if(argc==1){
    usage();
    return 1;
  }
 
  for(i=1;i<argc;i++){
    arg=argv[i];
    if(arg[0]!='-'||strlen(arg)<=1)break;

    switch(arg[1]){
      case 'a':
        at_end=1;
        add_path=1;
        remove_path=0;
        newentry=argv[i+1];
        i++;
        break;
      case 'f':
        at_end=0;
        add_path=1;
        remove_path=0;
        newentry=argv[i+1];
        i++;
        break;
      case 'd':
        display_path=1;
        break;
      case 'b':
        batch_mode=1;
        break;
      case 'm':
        path_summary=1;
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
        prompt_user_flag=1;
        break;
      case 't':
        test_mode=1;
        break;
      case 'v':
        version();
        return 0;
      default:
        usage();
        return 1;
    }
  }

  if(batch_mode==0&&remove_path==1){
    prompt_user_flag=1;
  }
  if(test_mode==1)display_path=1;
  if(remove_path==0&&add_path==0)display_path=1;

  // get the path (user or system)

  if(act_on_user_path==1){
    if(reg_path(REG_GET,REG_USER_PATH,pathbuffer)==0)return 1;
    if(display_path==1){
      if(strlen(pathbuffer)>0){
        printf("  User path:\n%s\n",pathbuffer);
      }
    }
  }
  else if(act_on_system_path==1){
    if(reg_path(REG_GET,REG_SYSTEM_PATH,pathbuffer)==0)return 1;
    if(display_path==1){
      if(strlen(pathbuffer)>0){
        printf("  System path:\n%s\n",pathbuffer);
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
    if(newentry_present==1){
      if(path_summary==1){
        printf("  %s already present in the %s path.\n",newentry,path_type);
      }
    }
    else{
      int answer=0;

      if(prompt_user_flag==1){
        answer = prompt_user(path_type,pathbuffer);
      }
      if(prompt_user_flag==0||answer==1){
        backup_path(path_type,pathbuffer);
        if(at_end==1){
          if(strlen(pathbuffer)>0)strcat(pathbuffer,";");
          strcat(pathbuffer,newentry);
        }
        else{
          strcpy(pathbuffer2,newentry);
          strcat(pathbuffer2,";");
          strcat(pathbuffer2,pathbuffer);
          strcpy(pathbuffer,pathbuffer2);
        }
        if(act_on_user_path==1){
          if(reg_path(REG_SET,REG_USER_PATH,pathbuffer)==0)return 1;
        }
        else{
          if(reg_path(REG_SET,REG_SYSTEM_PATH,pathbuffer)==0)return 1;
        }
        if(path_summary==1){
          printf("  %s was added to the %s path.\n",newentry,path_type);
          printf("  Reboot your computer so path changes may take effect.\n");
        }
        if(display_path==1){
          printf("\n  %s path set to:\n%s\n",path_type,pathbuffer);
        }
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

      if(prompt_user_flag==1){
        answer = prompt_user(path_type,newpath);
      }
      if(prompt_user_flag==0||answer==1){
        backup_path(path_type,pathbuffer);
        if(act_on_user_path==1){
          if(reg_path(REG_SET,REG_USER_PATH,newpath)==0)return 1;
        }
        else{
          if(reg_path(REG_SET,REG_SYSTEM_PATH,newpath)==0)return 1;
        }
        if(path_summary==1){
          printf("  All directories containing %s were removed from the %s path.\n",newentry,path_type);
          printf("  Reboot your computer so the path change may take effect.\n");
        }
      }
      if(display_path==1){
        printf("\n  %s path set to:\n%s\n",path_type,newpath);
      }
    }
    else{
      if(path_summary==1){
        printf("  %s not found in the %s path.\n",newentry,path_type);
      }
    }
  }
  return 0;
}

/* ------------------ reg_path ------------------------ */

int reg_path(int setget, int pathtype, char *path){

  HKEY hKey, hTree=NULL;
  long lRet;
  char temp[10000];
  DWORD dwBufLen;
  int lenpath;

  LPCTSTR reg_path_local=NULL;
  
  char creg_user_path[]="Environment";
  LPCTSTR reg_user_path=creg_user_path;
  
  char creg_system_path[]="SYSTEM\\CurrentControlSet\\Control\\Session Manager\\Environment";
  LPCTSTR reg_system_path=(LPCTSTR)creg_system_path;

  char cPATH[]="Path";
  LPCTSTR PATH=(LPCTSTR)cPATH;

  switch (pathtype) {
    case REG_USER_PATH:
      reg_path_local=reg_user_path;
      hTree=HKEY_CURRENT_USER;
      break;
    case REG_SYSTEM_PATH:
      reg_path_local=reg_system_path;
      hTree=HKEY_LOCAL_MACHINE;
      break;
    default:
      ASSERT(0);
      break;
  }
  switch (setget) {
    case REG_GET:
      lRet = RegOpenKeyEx( hTree, reg_path_local, 0, KEY_QUERY_VALUE, &hKey );
      switch (lRet){
        case ERROR_FILE_NOT_FOUND:
          printf("%s path not present\n",path_type);
          strcpy(path,"");
          dwBufLen=0;
          return 0;
        case ERROR_SUCCESS:
          break;
        default:
          printf("RegOpenKeyEx error: %i\n",(int)lRet);
          return 0;
      }
      dwBufLen=sizeof(temp);
      lRet = RegQueryValueEx( hKey, PATH, NULL, NULL, (BYTE*)&temp, &dwBufLen );
      switch (lRet){
        case ERROR_SUCCESS:
          strncpy(path,temp,dwBufLen);
          path[dwBufLen]=0;
          break;
        case 2:
          strcpy(path," ");
          dwBufLen=0;
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
      lRet = RegOpenKeyEx( hTree, reg_path_local, 0, KEY_QUERY_VALUE | KEY_SET_VALUE, &hKey );
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
      if(test_mode==1){
        printf("Setting %s path to: %s\n",path_type,path);
      }
      else{
        lRet = RegSetValueEx(hKey,PATH,0,REG_EXPAND_SZ,(LPBYTE)path,strlen(path)+1);
        lRet=ERROR_SUCCESS;
        if(lRet!=ERROR_SUCCESS){
          printf("RegSetValueEx error: %i\n",(int)lRet);
          return 0;
        }
      }
      lRet = RegCloseKey( hKey);
      if(lRet!=ERROR_SUCCESS){
        printf("RegCloseKey error: %i\n",(int)lRet);
        return 0;
      }
      break;
    default:
      ASSERT(0);
      break;
  }
  return 1;
}

/* ------------------ version ------------------------ */

void version(void){
  char revision[256];
  char version[256];

    getPROGversion(version);

    getRevision(revision);    // get svn revision number
    printf("\n");
    printf("set_path %s - %s\n\n",version,__DATE__);
    printf("Version: %s\n",version);
    printf("Revision: %s\n",revision);
    printf("Build Date: %s\n",__DATE__);
}

/* ------------------ usage ------------------------ */

void usage(void){
  char revision[100];

  getRevision(revision);

  printf("set_path Revision:%s\n",revision);
  printf("  Modify or display the User or System path environmental variables.\n\n");
  printf("Usage:\n\n");
  printf("  set_path [-s][-u] [-a path_entry] [-r path_entry] [-d][-p][-v]\n\n");
  printf("where\n\n");
  printf("  -a entry - append entry to the path variable being modified\n");
  printf("  -f entry - prepend entry to the path variable being modified\n");
  printf("  -r label - remove any entry containing label from the path\n");
  printf("             variable being modified\n");
  printf("  -s - add/remove/display entries in the System path\n");
  printf("  -u - add/remove/display entries in the User path (default)\n");
  printf("  -m - display a summary message changes made or made to the path\n");
  printf("  -d - display path before and after changes are made\n");
  printf("  -b - batch or script mode.  Override prompt option when set_path is run from a script\n");
  printf("  -p - prompt user before making any changes (default when path entries are being removed)\n");
  printf("  -t - test, show but do not change Path variables\n");
  printf("  -v - show versioning information\n");
}
