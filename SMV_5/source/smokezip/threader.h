// $Date$ 
// $Revision$
// $Author$

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

// setup LOCKS

#ifdef pp_THREAD
#define LOCK_PART2ISO if(mt_part2iso==1)pthread_mutex_lock(&mutexPART2ISO);
#define UNLOCK_PART2ISO if(mt_part2iso==1)pthread_mutex_unlock(&mutexPART2ISO);
#define LOCK_SMOKE if(mt_smoke==1)pthread_mutex_lock(&mutexSMOKE);
#define UNLOCK_SMOKE if(mt_smoke==1)pthread_mutex_unlock(&mutexSMOKE);
#else
#define LOCK_PART2ISO
#define UNLOCK_PART2ISO
#define LOCK_SMOKE
#define UNLOCK_SMOKE
#endif

// define mutex's and thread_ids

#ifdef pp_THREAD
MT_EXTERN int mt_part2iso, mt_smoke;
#endif

#ifndef CPP
#ifdef pp_THREAD
MT_EXTERN pthread_mutex_t mutexPART2ISO;
MT_EXTERN pthread_mutex_t mutexSMOKE;
#endif
#endif

void init_multi_threading(void);
void *MT_part2iso(void *arg);
void MT_convert_part2iso(void);
void MT_compress_smoke3ds(void);
void *MT_convert_3dsmoke(void *arg);

