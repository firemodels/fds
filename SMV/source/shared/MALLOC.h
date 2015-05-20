#ifndef MALLOC_H_DEFINED
#define MALLOC_H_DEFINED

#ifdef pp_THREAD
#include <pthread.h>
#endif

#ifndef pp_OSX
#include <malloc.h>
#endif
#include "ASSERT.h"
#include "string_util.h"
#include "file_util.h"
#define memGarbage 0xA3

typedef int mallocflag;
typedef char bbyte;

#ifdef CPP
#define MMCCC "C"
#else
#define MMCCC
#endif

#ifdef INDMALLOC
#define MMEXTERN
#else
#define MMEXTERN extern MMCCC
#endif

#ifdef pp_MEMDEBUG
#ifdef BIT64
#define MMsize unsigned long long
#else
#define MMsize unsigned int
#endif
#endif

typedef struct {
  unsigned char marker;
  void *prev, *next;
  int memory_id;
#ifdef pp_MEMPRINT
  size_t size;
#endif  
  
} MMdata;

MMEXTERN MMdata MMfirst, MMlast;
MMEXTERN MMdata *MMfirstptr, *MMlastptr;
#ifdef pp_MEMDEBUG
MMEXTERN MMsize MMtotalmemory,MMmaxmemory;
#endif

#define debugByte 0xE1
#define markerByte 0xE1
#ifdef pp_MEMDEBUG
  #define sizeofDebugByte 1
#else
  #define sizeofDebugByte 0
#endif

#ifdef pp_THREAD
MMEXTERN pthread_mutex_t mutexSLICE_BOUND,mutexPATCH_BOUND,mutexPART2ISO,mutexPRINT,mutexMEM;
#define LOCK_MEM           pthread_mutex_lock(&mutexMEM)
#define UNLOCK_MEM         pthread_mutex_unlock(&mutexMEM)
#else
#define LOCK_MEM
#define UNLOCK_MEM
#endif

#ifdef pp_MEMDEBUG
#define NewMemory(f,g)           __NewMemory((f),(g),0,(#f),__FILE__,__LINE__)
#define ResizeMemory(f,g)        __ResizeMemory((f),(g),0,(#f),__FILE__,__LINE__)

#define NewMemoryMemID(f,g,h)    __NewMemory((f),(g),(h),(#f),__FILE__,__LINE__)
#define ResizeMemoryMemID(f,g,h) __ResizeMemory((f),(g),(j),(#f),__FILE__,__LINE__)
#else
#define NewMemory(f,g)           _NewMemory((f),(g),0,(#f),__FILE__,__LINE__)
#define ResizeMemory(f,g)        _ResizeMemory((f),(g),0,(#f),__FILE__,__LINE__)

#define NewMemoryMemID(f,g,h)    _NewMemory((f),(g),(h),(#f),__FILE__,__LINE__)
#define ResizeMemoryMemID(f,g,h) _ResizeMemory((f),(g),(h),(#f),__FILE__,__LINE__)
#endif

#define NEWMEMORY(f,g) NewMemory((void **)&(f),(g))

#ifdef pp_memstatus
#ifdef WIN32
void _memorystatus(unsigned int size,unsigned int *availmem, unsigned int *memused, unsigned int *totalmem);
#define MEMSTATUS(f,g,h,i) _memorystatus(f,g,h,i)
#else
#define MEMSTATUS(f,g,h,i)
#endif
#endif

#ifndef pp_memstatus
#define MEMSTATUS(f,g,h,i)
#endif

#ifdef pp_MEMPRINT
void _PrintMemoryInfo(void);
#define PrintMemoryInfo _PrintMemoryInfo()
#else
#define PrintMemoryInfo
#endif
#ifdef pp_MEMDEBUG
void getMemusage(MMsize totalmemory,char *MEMlabel);
void _CheckMemory(void);
void _CheckMemoryNOTHREAD(void);
void _CheckMemoryOn(void);
void _CheckMemoryOff(void);
void _PrintAllMemoryInfo(void);
int _CountMemoryBlocks(void);
#define ValidPointer(pv,size) _ValidPointer(pv, size)
#define GetTotalMemory(size) size=_GetTotalMemory()
#define CheckMemory _CheckMemory()
#define CheckMemoryNOTHREAD _CheckMemoryNOTHREAD()
#define CheckMemoryOn _CheckMemoryOn()
#define CheckMemoryOff _CheckMemoryOff()
#define PrintAllMemoryInfo _PrintAllMemoryInfo()
#define CountMemoryBlocks(f,g)  f=_CountMemoryBlocks()-g
char *_strcpy(char *s1, const char *s2);
char *_strcat(char *s1, const char *s2);
#define STRCPY(f,g) _strcpy((f),(g))
#define STRCAT(f,g) _strcat((f),(g))
#else
#define ValidPointer(pv,size)
#define GetTotalMemory 
#define CheckMemory
#define CheckMemoryOn
#define CheckMemoryOff
#define PrintAllMemoryInfo
#define STRCPY(f,g) strcpy((f),(g))
#define STRCAT(f,g) strcat((f),(g))
#define CountMemoryBlocks(f,g)
#endif

#ifdef pp_MEMDEBUG
typedef struct BLOCKINFO {
  struct BLOCKINFO *pbiNext;
  bbyte *pb;
  size_t size;
  mallocflag  fref;
  char filename[256], varname[256];
  int linenumber;
} blockinfo;

MMEXTERN MMsize _GetTotalMemory(void);
mallocflag CreateBlockInfo(bbyte *pbNew, size_t sizeNew);
void FreeBlockInfo(bbyte *pb);
void UpdateBlockInfo(bbyte *pbOld, bbyte *pbNew, size_t sizeNew);
size_t sizeofBlock(bbyte *pv);
MMEXTERN mallocflag __ResizeMemory(void **ppv, size_t sizeNew, int memory_id, char *varname, char *file, int linenumber);
MMEXTERN mallocflag __NewMemory(void **ppv, size_t size, int memory_id, char *varname, char *file,int linenumber);
#endif
MMEXTERN mallocflag _ResizeMemory(void **ppv, size_t sizeNew, int memory_id, char *varname, char *file,int linenumber);
MMEXTERN mallocflag _NewMemory(void **ppv, size_t size, int memory_id, char *varname, char *file,int linenumber);
MMEXTERN void FreeMemory(void *pv);
MMEXTERN mallocflag _ResizeMemoryNOTHREAD(void **ppv, size_t sizeNew, int memory_id);
MMEXTERN mallocflag _NewMemoryNOTHREAD(void **ppv, size_t size, int memory_id);
MMEXTERN void FreeMemoryNOTHREAD(void *pv);
MMEXTERN int memusage(void);
void initMALLOC(void);
void FreeAllMemory(int memory_id);
mallocflag _ValidPointer(void *pv, size_t size);

#define FREEMEMORY(f) if((f)!=NULL){LOCK_MEM;FreeMemoryNOTHREAD((f));UNLOCK_MEM;(f)=NULL;}
#endif
