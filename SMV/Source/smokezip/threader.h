#ifndef THREADER_H_DEFINED
#define THREADER_H_DEFINED
#ifndef CPP
#ifdef pp_THREAD
#include <pthread.h>
#endif
#endif

#ifdef CPP
#define CCC "C"
#else
#define CCC
#endif

#ifdef INTHREADER
#define MT_EXTERN
#else
#define MT_EXTERN extern CCC
#endif

typedef struct {
  int file_index, thread_index, stat;
  char label[256];
} threaddata;


// setup LOCKS

#ifdef pp_THREAD
#define LOCK_COMPRESS      pthread_mutex_lock(&mutexCOMPRESS);
#define UNLOCK_COMPRESS    pthread_mutex_unlock(&mutexCOMPRESS);
#define LOCK_PATCH         pthread_mutex_lock(&mutexPATCH);
#define UNLOCK_PATCH       pthread_mutex_unlock(&mutexPATCH);
#define LOCK_PATCH_BOUND   pthread_mutex_lock(&mutexPATCH_BOUND);
#define UNLOCK_PATCH_BOUND pthread_mutex_unlock(&mutexPATCH_BOUND);
#define LOCK_SLICE         pthread_mutex_lock(&mutexSLICE);
#define UNLOCK_SLICE       pthread_mutex_unlock(&mutexSLICE);
#define LOCK_VOLSLICE         pthread_mutex_lock(&mutexVOLSLICE);
#define UNLOCK_VOLSLICE       pthread_mutex_unlock(&mutexVOLSLICE);
#define LOCK_SLICE_BOUND   pthread_mutex_lock(&mutexSLICE_BOUND);
#define UNLOCK_SLICE_BOUND pthread_mutex_unlock(&mutexSLICE_BOUND);
#define LOCK_ISOS          pthread_mutex_lock(&mutexISOS);
#define UNLOCK_ISOS        pthread_mutex_unlock(&mutexISOS);
#define LOCK_SMOKE         pthread_mutex_lock(&mutexSMOKE);
#define UNLOCK_SMOKE       pthread_mutex_unlock(&mutexSMOKE);
#define LOCK_PLOT3D        pthread_mutex_lock(&mutexPLOT3D);
#define UNLOCK_PLOT3D      pthread_mutex_unlock(&mutexPLOT3D);
#define LOCK_PART2ISO      pthread_mutex_lock(&mutexPART2ISO);
#define UNLOCK_PART2ISO    pthread_mutex_unlock(&mutexPART2ISO);
#define LOCK_PRINT         pthread_mutex_lock(&mutexPRINT);
#define UNLOCK_PRINT       pthread_mutex_unlock(&mutexPRINT);
#else
#define LOCK_COMPRESS
#define UNLOCK_COMPRESS
#define LOCK_PATCH
#define UNLOCK_PATCH
#define LOCK_PATCH_BOUND
#define UNLOCK_PATCH_BOUND
#define LOCK_SLICE
#define UNLOCK_SLICE
#define LOCK_VOLSLICE
#define UNLOCK_VOLSLICE
#define LOCK_SLICE_BOUND
#define UNLOCK_SLICE_BOUND
#define LOCK_ISOS
#define UNLOCK_ISOS
#define LOCK_SMOKE
#define UNLOCK_SMOKE
#define LOCK_PLOT3D
#define UNLOCK_PLOT3D
#define   LOCK_PART2ISO
#define UNLOCK_PART2ISO
#define   LOCK_PRINT
#define UNLOCK_PRINT
#endif

// define mutex's and thread_ids

#ifdef pp_THREAD
MT_EXTERN int mt_compress;
MT_EXTERN int mt_nthreads;
#endif

#ifndef CPP
#ifdef pp_THREAD
MT_EXTERN pthread_mutex_t mutexCOMPRESS,mutexPATCH,mutexSLICE,mutexISOS,mutexSMOKE,mutexPLOT3D,mutexVOLSLICE;
MT_EXTERN pthread_mutex_t mutexSLICE_BOUND,mutexPATCH_BOUND,mutexPART2ISO,mutexPRINT;
#endif
#endif

void init_pthread_mutexes(void);
void print_thread_stats(void);

#define NTHREADS_MAX 16

#endif
