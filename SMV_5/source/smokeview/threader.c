// $Date$ 
// $Revision$
// $Author$
#define INTHREADER
#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"

// svn revision character string
char threader_revision[]="$Revision$";

void compress_svzip2(void);


/* ------------------ init_all_threads ------------------------ */

void init_multi_threading(void){
#ifdef pp_THREAD
  if(mt_compress==1)pthread_mutex_init(&mutexCOMPRESS,NULL);
#endif
}

// *************** multi-threaded compression **************** 

#ifdef pp_THREAD
 /* ------------------ mt_compress_svzip ------------------------ */

void *mt_compress_svzip(void *arg){

  LOCK_COMPRESS
  compress_svzip2();
  updatemenu=1;
  UNLOCK_COMPRESS
  pthread_exit(NULL);
  return NULL;

   }
#endif

/* ------------------ compress_svzip ------------------------ */
#ifdef pp_THREAD
void compress_svzip(void){
  if(mt_compress==1){
    pthread_create( 
      &compress_thread_id, 
      NULL, 
      mt_compress_svzip, 
      NULL
      );
  }
  else{
    compress_svzip2();
  }
}
#else
void compress_svzip(void){
  compress_svzip2();
}
#endif
// ************** multi threaded blockage smoothing **********************

/* ------------------ mt_update_smooth_blockages ------------------------ */
#ifdef pp_THREAD
void *mt_update_smooth_blockages(void *arg){

  if(ifsmoothblock()==1){
    printf("Smoothing blockages in the background\n");
    update_smooth_blockages();
    updatefacelists=1;
  }
  pthread_exit(NULL);
  return NULL;

   }
#endif

/* ------------------ smooth_blockages ------------------------ */
#ifdef pp_THREAD
void smooth_blockages(void){
  smoothing_blocks=1;
  if(mt_compress==1){
    pthread_create( 
      &smooth_block_thread_id,
      NULL, 
      mt_update_smooth_blockages, 
      NULL
      );
  }
  else{
    blocksneedsmoothing=ifsmoothblock();
    if(blocksneedsmoothing==1){
      update_smooth_blockages();
    }
  }
}
#else
void smooth_blockages(void){
  smoothing_blocks=1;
    blocksneedsmoothing=ifsmoothblock();
    if(blocksneedsmoothing==1){
      update_smooth_blockages();
    }
}
#endif
