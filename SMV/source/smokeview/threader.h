// $Date$ 
// $Revision$
// $Author$

#ifndef CPP
#ifdef pp_THREAD
#include <pthread.h>
#endif
#endif

#ifdef INTHREADER
#define MT_EXTERN
#else
#define MT_EXTERN extern CCC
#endif

// setup LOCKS

#ifdef pp_THREAD
#define LOCK_COMPRESS pthread_mutex_lock(&mutexCOMPRESS);
#define UNLOCK_COMPRESS pthread_mutex_unlock(&mutexCOMPRESS);
#else
#define LOCK_COMPRESS
#define UNLOCK_COMPRESS
#endif

// define mutex's and thread_ids

#ifndef CPP
#ifdef pp_THREAD
MT_EXTERN pthread_mutex_t mutexCOMPRESS;
MT_EXTERN pthread_t smooth_block_thread_id;
MT_EXTERN pthread_t compress_thread_id;
MT_EXTERN pthread_t update_all_patch_bounds_id;
#endif
#endif
