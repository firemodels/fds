// $Date$ 
// $Revision$
// $Author$
#define INTHREADER
#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "zlib.h"
#include "ASSERT.h"
#include "svzip.h"

// svn revision character string
char threader_revision[]="$Revision$";

void part2iso(part *parti);

/* ------------------ init_all_threads ------------------------ */

void init_multi_threading(void){
#ifdef pp_THREAD
  if(mt_part2iso==1)pthread_mutex_init(&mutexPART2ISO,NULL);
  if(mt_smoke==1)pthread_mutex_init(&mutexSMOKE,NULL);
#endif
}

// ************** multi threaded blockage smoothing **********************

/* ------------------ MT_part2iso ------------------------ */
#ifdef pp_THREAD
void *MT_part2iso(void *arg){
  part *parti;

  parti = (part *)arg;
  part2iso(parti);
  return NULL;

}

void *MT_convert_3dsmoke(void *arg){
  smoke3d *smoke3di;
  int smoke_index;

  smoke3di = (smoke3d *)arg;
  smoke_index=smoke3di-smoke3dinfo;
  //LOCK_SMOKE;
  convert_3dsmoke(smoke3di);
  //UNLOCK_SMOKE;
 // printf("smoke3d compression for file: %i\n",smoke_index);
  return NULL;

}
#endif


