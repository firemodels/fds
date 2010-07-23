// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
#include <process.h>
#include <windows.h>
#endif
#include "svn_revision.h"
#include "background.h"

//dummy change to bump version number to 0.9

// svn revision character string
char main_revision[]="$Revision$";
void run_command(void);

void usage(char *prog);
#ifdef WIN32
void GetSystemTimesAddress(void);
#endif
unsigned char cpuusage(void);

#ifdef pp_LINUX
int get_ncores(void);
float get_load(void);
float get_host_load(char *host);
unsigned char cpuusage_host(char *host,int ncores);
#endif

#ifdef pp_LINUX
void Sleep(int ticks){
  float time;

  time = ticks/1000.0;
  sleep(time);
}
#endif

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char *prog;
  int i;
  int debug;
  int argstart=-1;
  float delay_time=0.0;
  int cpu_usage, cpu_usage_max=25;
#ifdef pp_LINUX  
  char command_buffer[1024];
  FILE *stream=NULL;
#endif  

  int itime;
  char *arg;
  char *command;
  char *command_arg;
  int n;

#ifdef pp_LINUX  
  hostlistfile=NULL;
  host=NULL;
#ifdef pp_LINUX
  sprintf(pid,"%i",getpid());
#endif
#endif

  debug=0;
  prog=argv[0];

  if(argc==1){
    version();
    return 1;
  }

  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'){
      if(lenarg>1){
        switch(arg[1]){
          case 'd':
            if(strlen(arg)<=2||strcmp(arg,"-debug")!=0){
              i++;
              if(i<argc){
                arg=argv[i];
                sscanf(arg,"%f",&delay_time);
                if(delay_time<0.0)delay_time=0.0;
              }
            }
            else{
		          debug=1;
		        }
            break;
          case 'h':
#ifdef pp_LINUX
            if(strlen(arg)<=2||strcmp(arg,"-hosts")!=0){
#endif
              usage(prog);
              return 1;
#ifdef pp_LINUX
            }
            else{
              i++;
              if(i<argc){
                arg=argv[i];
                hostlistfile=malloc(strlen(arg)+1);
                strcpy(hostlistfile,arg);
              }
            }
#endif
            break;
          case 'u':
            i++;
            if(i<argc){
              arg=argv[i];
              sscanf(arg,"%i",&cpu_usage_max);
              if(cpu_usage_max<25)cpu_usage_max=25;
              if(cpu_usage_max>100)cpu_usage_max=100;
            }
            break;
          case 'v':
            version();
            return 1;
            break;
          default:
            printf("Unknown option: %s\n",arg);
            usage(prog);
            return 1;
            break;
	      }
      }
    }
    else{
      argstart=i;
      break;

    }
  }
#ifdef pp_LINUX
  nhostinfo=0;
  if(hostlistfile!=NULL){
    stream=fopen(hostlistfile,"r");
  }
  if(hostlistfile!=NULL&&stream!=NULL){
    char buffer[255];

    while(!feof(stream)){
      if(fgets(buffer,255,stream)==NULL)break;
      nhostinfo++;
    }
    if(nhostinfo>0){
      hostdata *hd;

      hostinfo=malloc(nhostinfo*sizeof(hostdata));
      hd = hostinfo;
      rewind(stream);
      while(!feof(stream)){
        char *hostname;

        if(fgets(buffer,255,stream)==NULL)break;
        trim(buffer);
        hostname=malloc(strlen(buffer)+1);
        strcpy(hostname,buffer);
        hd->hostname=hostname;
        hd->ncores=get_host_ncores(hostname);
        printf("host: %s ncores=%i\n",hostname,hd->ncores);
        hd++;
      }

    }
  }
  if(stream!=NULL)fclose(stream);
#endif

  if(argstart<0)return 0;

  itime = delay_time*1000;
  if(itime>0){
    Sleep(itime);
  }

#ifdef WIN32
  GetSystemTimesAddress();
  cpu_usage=cpuusage();
  Sleep(200);
  cpu_usage=cpuusage();
  while(cpu_usage>cpu_usage_max){
    Sleep(2000);
    cpu_usage=cpuusage();
    Sleep(200);
    cpu_usage=cpuusage();
  }
  command=argv[argstart];
  _spawnvp(_P_NOWAIT,command, argv+argstart);
#endif
#ifdef pp_LINUX
  if(hostinfo==NULL){
    cpu_usage=cpuusage();
    Sleep(200);
    cpu_usage=cpuusage();
    host=NULL;
    while(cpu_usage>cpu_usage_max){
      cpu_usage=cpuusage();
      Sleep(1000);
    }
  }
  else{
    int doit=0;

    while(doit==0){
      for(i=0;i<nhostinfo;i++){
        hostdata *hd;
        float fusage;

        hd = hostinfo + i;
        host = hd->hostname;
        cpu_usage=cpuusage_host(host,hd->ncores);
        fusage=(float)cpu_usage/255.0;
        if(debug==1)printf("host: %s cpu_usage=%f\n",host,fusage);
        if(cpu_usage<cpu_usage_max){
          if(debug==1)printf(" host %s is now free\n",host);
          doit=1;
          break;
        }
      }
      if(doit==0)Sleep(1000);
    }
  }
  strcpy(command_buffer,"");
  for(i=argstart;i<argc;i++){
    arg=argv[i];
    strcat(command_buffer,arg);
    if(i==argc-1){
      strcat(command_buffer,"  &");    
    }
    else{
      strcat(command_buffer," ");    
    }
  }
  printf("submitting command: %s\n",command_buffer);
  system(command_buffer);
#endif
  return 0;
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
  char prog_version[100];
  int svn_num;
  char pp[]="%";

  getPROGversion(prog_version);  // get version (ie 5.x.z)
  svn_num=getmaxrevision();    // get svn revision number

  printf("\n");
  printf("background %s(%i) - %s\n",prog_version,svn_num,__DATE__);
  printf("  Runs a program in the background when resources are available\n\nUsage:\n\n");
  printf("  %s",prog);
  printf(" [-d delay time (s) -h -u max_usage -v] prog [arguments]\n\n");

  printf("where\n\n");

  printf("  -d dtime  - wait dtime seconds before running prog in the background\n");
  printf("  -debug    - display debug messages\n");
  printf("  -h        - display this message\n");
#ifdef pp_LINUX  
  printf("  -host hostfiles - file containing a list of host names to run jobs on\n");
#endif  
  printf("  -u max    - wait to run prog until cpu usage is less than max (25-100%s)\n",pp);
  printf("  -v        - display version information\n");
  printf("  prog      - program to run in the background\n");
  printf("  arguments - command line arguments of prog\n\n");
  printf("Example:\n");
  printf("  background -d 1.5 -u 50 prog arg1 arg2\n");
  printf("    runs prog (with arguments arg1 and arg2) after 1.5 seconds\n    and when the CPU usage drops below 50%s\n",pp);
}

#ifdef WIN32
typedef BOOL ( __stdcall * pfnGetSystemTimes)( LPFILETIME lpIdleTime, LPFILETIME lpKernelTime, LPFILETIME lpUserTime );
static pfnGetSystemTimes s_pfnGetSystemTimes = NULL;

static HMODULE s_hKernel = NULL;

/* ------------------ GetSystemTimesAddress ------------------------ */

void GetSystemTimesAddress(){
	if( s_hKernel == NULL )
	{   
		s_hKernel = LoadLibrary( "Kernel32.dll" );
		if( s_hKernel != NULL )
		{
			s_pfnGetSystemTimes = (pfnGetSystemTimes)GetProcAddress( s_hKernel, "GetSystemTimes" );
			if( s_pfnGetSystemTimes == NULL )
			{
				FreeLibrary( s_hKernel ); s_hKernel = NULL;
			}
		}
	}
}

/* ------------------ cpuusage ------------------------ */

unsigned char cpuusage()
{
	FILETIME               ft_sys_idle;
	FILETIME               ft_sys_kernel;
	FILETIME               ft_sys_user;

	ULARGE_INTEGER         ul_sys_idle;
	ULARGE_INTEGER         ul_sys_kernel;
	ULARGE_INTEGER         ul_sys_user;

	static ULARGE_INTEGER	 ul_sys_idle_old;
	static ULARGE_INTEGER  ul_sys_kernel_old;
	static ULARGE_INTEGER  ul_sys_user_old;

	unsigned char usage = 0;

	// we cannot directly use GetSystemTimes on C language
	/* add this line :: pfnGetSystemTimes */
	s_pfnGetSystemTimes(&ft_sys_idle,    /* System idle time */
		&ft_sys_kernel,  /* system kernel time */
		&ft_sys_user);   /* System user time */

	CopyMemory(&ul_sys_idle  , &ft_sys_idle  , sizeof(FILETIME)); // Could been optimized away...
	CopyMemory(&ul_sys_kernel, &ft_sys_kernel, sizeof(FILETIME)); // Could been optimized away...
	CopyMemory(&ul_sys_user  , &ft_sys_user  , sizeof(FILETIME)); // Could been optimized away...

	usage  =
		(
		(
		(
		(
		(ul_sys_kernel.QuadPart - ul_sys_kernel_old.QuadPart)+
		(ul_sys_user.QuadPart   - ul_sys_user_old.QuadPart)
		)
		-
		(ul_sys_idle.QuadPart-ul_sys_idle_old.QuadPart)
		)
		*
		(100)
		)
		/
		(
		(ul_sys_kernel.QuadPart - ul_sys_kernel_old.QuadPart)+
		(ul_sys_user.QuadPart   - ul_sys_user_old.QuadPart)
		)
		);

	ul_sys_idle_old.QuadPart   = ul_sys_idle.QuadPart;
	ul_sys_user_old.QuadPart   = ul_sys_user.QuadPart;
	ul_sys_kernel_old.QuadPart = ul_sys_kernel.QuadPart;

	return usage;
}
#endif
#ifdef pp_LINUX

/* ------------------ get_ncores ------------------------ */

int get_ncores(void){
  FILE *stream;
  int ncores=0;
  char buffer[255];

  stream=fopen("/proc/cpuinfo","r");
  if(stream==NULL)return 1;
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    if(strlen(buffer)<9)continue;
    buffer[9]=0;
    if(strcmp(buffer,"processor")==0)ncores++;
  }
  if(ncores==0)ncores=1;
  fclose(stream);
  return ncores;
}

/* ------------------ get_host_ncores ------------------------ */

int get_host_ncores(char *host){
  FILE *stream;
  char buffer[1024];
  char command[1024];
  char localfile[1024];
  int ncores=0;
  
  strcpy(localfile,"/tmp/cpuinfo.");
  strcat(localfile,host);
  strcat(localfile,".");
  strcat(localfile,pid);

  strcpy(command,"ssh ");
  strcat(command,host);
  strcat(command," cat /proc/cpuinfo >");
  strcat(command,localfile);

  system(command);

  stream=fopen(localfile,"r");
  if(stream==NULL){
    printf("unable to open %s\n",localfile);
    return 1;
  }
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    if(strlen(buffer)<9)continue;
    buffer[9]=0;
    if(strcmp(buffer,"processor")==0)ncores++;
  }
  if(ncores==0){
    printf("0 cores found in %s\n",localfile);
    ncores=1;
  }
  fclose(stream);
  unlink(localfile);
  return ncores;
}

/* ------------------ get_host_load ------------------------ */

float get_host_load(char *host){
  FILE *stream;
  char buffer[1024];
  char command[1024];
  char localfile[1024];
  float load1;
  
  strcpy(localfile,"/tmp/loadavg.");
  strcat(localfile,host);
  strcat(localfile,".");
  strcat(localfile,pid);

  strcpy(command,"ssh ");
  strcat(command,host);
  strcat(command," cat /proc/loadavg >");
  strcat(command,localfile);

  system(command);

  stream=fopen(localfile,"r");
  if(stream==NULL)return 1.0;
  if(fgets(buffer,255,stream)==NULL)return 1.0;
  sscanf(buffer,"%f",&load1);
  fclose(stream);
  unlink(localfile);
  return load1;
}

/* ------------------ get_load ------------------------ */

float get_load(void){
  FILE *stream;
  char buffer[255];
  float load1;

  stream=fopen("/proc/loadavg","r");
  if(stream==NULL)return 1.0;
  if(fgets(buffer,255,stream)==NULL)return 1.0;
  sscanf(buffer,"%f",&load1);
  fclose(stream);
  return load1;
}

/* ------------------ cpuusage ------------------------ */

unsigned char cpuusage_host(char *host, int ncores){
  float load;
  unsigned char usage;

  load = get_host_load(host);
  if(load>ncores)load=ncores;
  usage = 100*(load/(float)ncores);
  return usage;
}

/* ------------------ cpuusage ------------------------ */

unsigned char cpuusage(){
  unsigned char usage;
  float load;
  int ncores;

  ncores = get_ncores();
  load = get_load();
  if(load>ncores)load=ncores;
  usage = 100*(load/(float)ncores);
  return usage;
}
#endif
