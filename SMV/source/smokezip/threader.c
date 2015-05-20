#define INTHREADER
#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include "MALLOC.h"
#include "zlib.h"
#include "svzip.h"

/* ------------------ mt_compress_all ------------------------ */
#ifdef pp_THREAD
void mt_compress_all(void){
  int i;
  pthread_t *thread_ids;
  int *index;

  NewMemory((void **)&thread_ids,mt_nthreads*sizeof(pthread_t));
  NewMemory((void **)&index,mt_nthreads*sizeof(int));
  NewMemory((void **)&threadinfo,mt_nthreads*sizeof(threaddata));

  for(i=0;i<mt_nthreads;i++){
    index[i]=i;
    pthread_create(&thread_ids[i],NULL,compress_all,&index[i]);
    threadinfo[i].stat=-1;
  }

  for(i=0;i<mt_nthreads;i++){
    pthread_join(thread_ids[i],NULL);
  }

  print_summary();
  FREEMEMORY(thread_ids);
  FREEMEMORY(index);
  FREEMEMORY(threadinfo);
}
#endif

/* ------------------ init_all_threads ------------------------ */

void init_pthread_mutexes(void){
#ifdef pp_THREAD
  pthread_mutex_init(&mutexCOMPRESS,NULL);
  pthread_mutex_init(&mutexPATCH,NULL);
  pthread_mutex_init(&mutexPATCH_BOUND,NULL);
  pthread_mutex_init(&mutexSLICE,NULL);
  pthread_mutex_init(&mutexVOLSLICE,NULL);
  pthread_mutex_init(&mutexSLICE_BOUND,NULL);
  pthread_mutex_init(&mutexISOS,NULL);
  pthread_mutex_init(&mutexSMOKE,NULL);
  pthread_mutex_init(&mutexPLOT3D,NULL);
  pthread_mutex_init(&mutexPART2ISO,NULL);
  pthread_mutex_init(&mutexPRINT,NULL);
#endif
}

/* ------------------ print_thread_stats ------------------------ */

void print_thread_stats(void){
#ifdef pp_THREAD
  int i;
  int sum;
  int lastthread;

  sum=0;
  for(i=0;i<mt_nthreads;i++){
    if(threadinfo[i].stat>0){
      sum++;
      lastthread=i;
    }
  }
  if(sum>0){
    for(i=0;i<mt_nthreads;i++){
      threaddata *ti;

      ti = threadinfo+i;
      if(ti->stat>0){
        PRINTF(" %s(%i%s)",ti->label,ti->stat,GLOBpp);
        if(lastthread!=i)PRINTF(",");
      }
    }
    PRINTF("\n");
  }
#endif
}


