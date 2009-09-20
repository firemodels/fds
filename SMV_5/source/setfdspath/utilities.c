// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include "svn_revision.h"
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>
#include <tchar.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "MALLOC.h"
#include "ASSERT.h"

char utilities_revision[]="$Revision$";

/* ------------------ imax ------------------------ */

int imax(int a, int b){
  if(a>b){
    return a;
  }
  else{
    return b;
  }
}

/* ------------------ trim ------------------------ */

void trim(char *line){
  size_t len, i;
  unsigned char *c;
  
  len = strlen(line);
  c = (unsigned char *)(line+len-1);
  for(i=0; i<len; i++){
    if(isspace(*c)==0){
      c++; 
      line[c-(unsigned char *)line]='\0';
      return;
    }
    c--;
  }
  *line='\0';
}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  unsigned char *c;
  size_t i,len;

  c = (unsigned char *)line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(isspace(*c)==0)return line+i;
    c++;
  }
  return line+len;
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
  char C[10000],*result, KEY[10000];
  size_t i, len_c,len_key;

  if(c==NULL||key==NULL)return NULL;
  len_c=strlen(c);
  len_key=strlen(key);
  if(len_c<1||len_key<1)return NULL;

  for(i=0;i<len_c;i++){
    C[i]=toupper(c[i]);
  }
  C[len_c]='\0';
  for(i=0;i<len_key;i++){
    KEY[i]=toupper(key[i]);
  }
  KEY[len_key]='\0';
  result = strstr(C,KEY);
  if(result==NULL){
    return NULL;
  }
  else{
    return result + (c - C); 
  }
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

/* ------------------ getmaxrev ------------------------ */

#define MAXREV(cval) max_revision=imax(getrevision(cval),max_revision)
int getmaxrevision(void){
  int max_revision=0;

  MAXREV(main_revision);
  MAXREV(utilities_revision);
  return max_revision;
}

/* ------------------ reg_path ------------------------ */

int reg_path(int setget, int pathtype, char *path){
  // reg_path(REG_USER_PATH,REG_GET,path);
  // reg_path(REG_USER_PATH,REG_SET,path);
  // reg_path(REG_SYSTEM_PATH,REG_GET,path);
  // reg_path(REG_SYSTEM_PATH,REG_SET,path);
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
      if(lRet!=ERROR_SUCCESS){
        printf("RegQueryValueEx error: %i\n",(int)lRet);
        return 0;
      }
      lRet = RegCloseKey( hKey);
      if(lRet!=ERROR_SUCCESS){
        printf("RegCloseKey error: %i\n",(int)lRet);
        return 0;
      }
      strncpy(path,temp,dwBufLen);
      path[dwBufLen]=0;
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
      lRet = RegSetValueEx(hKey,PATH,0,REG_EXPAND_SZ,(LPBYTE)path,strlen(path)+1);
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
