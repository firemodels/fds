// $Date$ 
// $Revision$
// $Author$

#ifndef MALLOC_DEFINED
#define MALLOC_DEFINED

#ifndef pp_OSX
#include <malloc.h>
#endif
#include "ASSERT.h"
#define memGarbage 0xA3

typedef int mallocflag;
typedef char bbyte;

#define debugByte 0xE1
#ifdef _DEBUG
  #define sizeofDebugByte 1
#else
  #define sizeofDebugByte 0
#endif

#ifdef _DEBUG
#define NewMemory(f,g) __NewMemory((f),(g),__FILE__,__LINE__)
#define ResizeMemory(f,g) __ResizeMemory((f),(g),__FILE__,__LINE__)
#else
#define NewMemory(f,g) _NewMemory((f),(g))
#define ResizeMemory(f,g) _ResizeMemory((f),(g))
#endif

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

#ifdef _DEBUG
void _CheckMemory(void);
void _CheckMemoryOn(void);
void _CheckMemoryOff(void);
void _PrintMemoryInfo(void);
int _GetMemoryInfo(void);
#define CheckMemory _CheckMemory()
#define CheckMemoryOn _CheckMemoryOn()
#define CheckMemoryOff _CheckMemoryOff()
#define PrintMemoryInfo _PrintMemoryInfo()
#define GetMemoryInfo(f,g)  f=_GetMemoryInfo()-g
char *_strcpy(char *s1, const char *s2);
char *_strcat(char *s1, const char *s2);
#define STRCPY(f,g) _strcpy((f),(g))
#define STRCAT(f,g) _strcat((f),(g))
#else
#define CheckMemory
#define CheckMemoryOn
#define CheckMemoryOff
#define PrintMemoryInfo
#define STRCPY(f,g) strcpy((f),(g))
#define STRCAT(f,g) strcat((f),(g))
#define GetMemoryInfo(f,g)
#endif

#ifdef _DEBUG
typedef struct BLOCKINFO {
  struct BLOCKINFO *pbiNext;
  bbyte *pb;
  size_t size;
  mallocflag  fref;
  char filename[256];
  int linenumber;
} blockinfo;

mallocflag CreateBlockInfo(bbyte *pbNew, size_t sizeNew);
void FreeBlockInfo(bbyte *pb);
void UpdateBlockInfo(bbyte *pbOld, bbyte *pbNew, size_t sizeNew);
size_t sizeofBlock(bbyte *pv);
mallocflag __ResizeMemory(void **ppv, size_t sizeNew,char *file,int linenumber);
mallocflag __NewMemory(void **ppv, size_t size,char *file,int linenumber);
#endif
mallocflag _ResizeMemory(void **ppv, size_t sizeNew);
mallocflag _NewMemory(void **ppv, size_t size);
void FreeMemory(void *pv);
mallocflag ValidPointer(void *pv, size_t size);

#endif
#define FREEMEMORY(f) if((void *)(f)!=NULL){FreeMemory((void *)(f));(f)=NULL;}
