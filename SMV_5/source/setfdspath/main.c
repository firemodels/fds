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

#define BUFFER_SIZE 10050

int add_path=0;
int remove_path=0;
int display_path=0;
int act_on_user_path=1;
int act_on_system_path=0;
int prompt_user=0;


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
  char path_type[10];

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
      case 'v':
        version();
        return 0;
        break;
      default:
        usage();
        return 1;
    }
  }

  // get the path (user or system)

  if(act_on_user_path==1){
    if(reg_path(REG_GET,REG_USER_PATH,pathbuffer_ptr)==0){
      return 1;
    }
  }
  else if(act_on_system_path==1){
    if(reg_path(REG_GET,REG_SYSTEM_PATH,pathbuffer)==0)return 1;
  }
  else{
    if(display_path==1){
      if(reg_path(REG_GET,REG_USER_PATH,pathbuffer)!=0){
        printf("User path: %s\n\n",pathbuffer);
      }
      if(reg_path(REG_GET,REG_SYSTEM_PATH,pathbuffer)!=0){
        printf("User path: %s\n\n",pathbuffer);
      }
    }
    else{
      usage();
    }
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
  printf("  -v - show versioning information\n");
}
