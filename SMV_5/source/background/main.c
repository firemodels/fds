// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <process.h>
#include <windows.h>
#include "svn_revision.h"
#include "background.h"

//dummy change to bump version number to 0.9

// svn revision character string
char main_revision[]="$Revision$";
void run_command(void);

void usage(char *prog);
void GetSystemTimesAddress(void);
CHAR cpuusage(void);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char *prog;
  int i;
  int argstart=-1;
  float delay_time=0.0;
  int cpu_usage, cpu_usage_max=25;

  int itime;
  char *arg;
  char *command;
  char *command_arg;
  int n;

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
            i++;
            if(i<argc){
              arg=argv[i];
              sscanf(arg,"%f",&delay_time);
              if(delay_time<0.0)delay_time=0.0;
            }
		        break;
          case 'h':
            usage(prog);
            return 1;
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

  if(argstart<0)return 0;

  itime = delay_time*1000;
  if(itime>0){
    Sleep(itime);
  }

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
  printf("  Runs a windows program in the background\n\nUsage:\n\n");
  printf("  %s",prog);
  printf(" [-d delay time (s) -h -u max_usage -v] prog [arguments]\n\n");

  printf("where\n\n");

  printf("  -d dtime  - wait dtime seconds before running prog in the background\n");
  printf("  -h        - display this message\n");
  printf("  -u max    - wait to run prog until cpu usage is less than max (25-100%s)\n",pp);
  printf("  -v        - display version information\n");
  printf("  prog      - program to run in the background\n");
  printf("  arguments - command line arguments of prog\n\n");
  printf("Example:\n");
  printf("  background -d 1.5 -u 50 prog arg1 arg2\n");
  printf("    runs prog (with arguments arg1 and arg2) after 1.5 seconds\n    and when the CPU usage drops below 50%s\n",pp);
}

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

CHAR cpuusage()
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

	CHAR  usage = 0;

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
