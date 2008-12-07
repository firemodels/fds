// $Date: 2007-10-09 15:35:36 -0400 (Tue, 09 Oct 2007) $ 
// $Revision: 824 $
// $Author: gforney $

#include "options.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "MALLOC.h"
#include "ASSERT.h"
#ifdef _DEBUG
static int checkmemoryflag=1;
#endif
#ifdef WIN32
#include <windows.h>
#endif

// svn revision character string
char dmalloc_revision[]="$Revision: 824 $";

#ifdef pp_memstatus
#ifdef WIN32
void _memorystatus(unsigned int size,unsigned int *availmem,unsigned int *physmemused, unsigned int *totalmem){
  MEMORYSTATUS stat;

    GlobalMemoryStatus (&stat);
    if(availmem!=NULL)*availmem=stat.dwMemoryLoad;
    if(totalmem!=NULL)*totalmem=stat.dwTotalPhys/(1024*1024);
    if(physmemused!=NULL)*physmemused=(stat.dwTotalPhys-stat.dwAvailPhys)/(1024*1024);
#ifdef _DEBUG
    if(size!=0&&size<=stat.dwAvailPhys-0.1*stat.dwTotalPhys){
      printf("*** Available Memory: %i M \n",
            stat.dwAvailPhys/(1024*1024));
    }
#endif
    if(size!=0&&size>stat.dwAvailPhys-0.1*stat.dwTotalPhys){
      printf("*** Low Memory Warning. Only %i M available for viewing data.\n",
            stat.dwAvailPhys/(1024*1024));
      printf("    Unload datafiles or system performance may degrade.\n");
    }
}
#endif
#endif

/* ------------------ _NewMemory ------------------------ */

mallocflag _NewMemory(void **ppv, size_t size){
  void **ppb=(void **)ppv;
#ifdef _DEBUG
  char *c;
#endif

  ASSERT(ppv != NULL && size != 0);
  *ppb = (void *)malloc(size+sizeofDebugByte);

#ifdef _DEBUG
  {
    CheckMemory;
    if(*ppb != NULL){
      if(sizeofDebugByte!=0){
       c = (char *)(*ppb) + size;
       *c=(char)debugByte;
      }
      memset(*ppb, memGarbage, size);
      if(!CreateBlockInfo(*ppb, size)){
        free(*ppb);
        *ppb=NULL;
      }
    }
    ASSERT(*ppb !=NULL);
  }
#endif
  return (*ppb != NULL);
}

/* ------------------ FreeMemory ------------------------ */

void FreeMemory(void *pv){
#ifdef _DEBUG
  int len_memory;
#endif
  ASSERT(pv != NULL);
#ifdef _DEBUG
  {
    CheckMemory;
    len_memory=sizeofBlock(pv);
    memset(pv, memGarbage, len_memory);
    FreeBlockInfo(pv);
  }
#endif
  free(pv);
}

/* ------------------ ResizeMemory ------------------------ */

mallocflag _ResizeMemory(void **ppv, size_t sizeNew){
  bbyte **ppb = (bbyte **)ppv;
  bbyte *pbNew;
#ifdef _DEBUG
  char *c;
  size_t sizeOld;
#endif
  ASSERT(ppb != NULL && sizeNew != 0);
#ifdef _DEBUG
  {
    CheckMemory;
    sizeOld = sizeofBlock(*ppb);
    if(sizeNew<sizeOld){
      memset((*ppb)+sizeNew,memGarbage,sizeOld-sizeNew);
    }
    else if (sizeNew > sizeOld){
      void *pbForceNew;
      if(_NewMemory((void **)&pbForceNew, sizeNew)){
        memcpy(pbForceNew, *ppb, sizeOld);
        FreeMemory(*ppb);
        *ppb = pbForceNew;
      }
    }
  }
#endif

  pbNew = realloc(*ppb, sizeNew+sizeofDebugByte);
  if(pbNew != NULL){
#ifdef _DEBUG
    {
      if(sizeofDebugByte!=0){
        c = pbNew + sizeNew;
        *c=(char)debugByte;
      }
      UpdateBlockInfo(*ppb, pbNew, sizeNew);
      if(sizeNew>sizeOld){
        memset(pbNew+sizeOld,memGarbage,sizeNew-sizeOld);
      }
    }
#endif
    *ppb = pbNew;
  }
  return (pbNew != NULL);
}

#ifdef _DEBUG
/* ------------------ pointer comparison defines ------------------------ */

#define fPtrLess(pLeft, pRight)   ((pLeft) <  (pRight))
#define fPtrGrtr(pLeft, pRight)   ((pLeft) >  (pRight))
#define fPtrEqual(pLeft, pRight)  ((pLeft) == (pRight))
#define fPtrLessEq(pLeft, pRight) ((pLeft) <= (pRight))
#define fPtrGrtrEq(pLeft, pRight) ((pLeft) >= (pRight))

static blockinfo *GetBlockInfo(bbyte *pb);
mallocflag __NewMemory(void **ppv, size_t size, char *file, int linenumber){
  void **ppb=(void **)ppv;
  blockinfo *pbi;
  int return_code;

  return_code=_NewMemory(ppb,size);
  pbi=GetBlockInfo((bbyte *)*ppb);
  pbi->linenumber=linenumber;
  if(strlen(file)<256){
    strcpy(pbi->filename,file);
  }
  else{
    strncpy(pbi->filename,file,255);
    strcat(pbi->filename,(char)0);
  }
  return return_code;
}

/* ------------------ ResizeMemory ------------------------ */

mallocflag __ResizeMemory(void **ppv, size_t size, char *file, int linenumber){
  void **ppb=(void **)ppv;
  blockinfo *pbi;
  int return_code;

  return_code=_ResizeMemory(ppb,size);
  pbi=GetBlockInfo((bbyte *)*ppb);
  pbi->linenumber=linenumber;
  if(strlen(file)<256){
    strcpy(pbi->filename,file);
  }
  else{
    strncpy(pbi->filename,file,255);
    strcat(pbi->filename,(char)0);
  }
  return return_code;
}

static blockinfo *pbiHead = NULL;

/* ------------------ GetBlockInfo ------------------------ */

static blockinfo *GetBlockInfo(bbyte *pb){
  blockinfo *pbi;
  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext)
  {
    bbyte *pbStart = pbi->pb;
    bbyte *pbEnd   = pbi->pb + pbi->size - 1;

    if(fPtrGrtrEq(pb, pbStart) && fPtrLessEq(pb, pbEnd))
      break;
  }
  ASSERT(pbi != NULL);
  return (pbi);
}

int _GetMemoryInfo(void){
  blockinfo *pbi;
  int n=0,size=0;

  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext)
  {
    n++;
    size += pbi->size;
  }
  return n;
}

/* ------------------ PrintMemoryInfo ------------------------ */

void _PrintMemoryInfo(void){
  blockinfo *pbi;
  int n=0,size=0;

  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext)
  {
    n++;
    size += pbi->size;
/*    printf("Block %i allocated in file=%s, linenumber=%i\n",n,pbi->filename,pbi->linenumber);*/
  }
  printf("nblocks=%i sizeblocks=%i\n",n,size);
}

/* ------------------ GetBlockInfo_nofail ------------------------ */

static blockinfo *GetBlockInfo_nofail(bbyte *pb){
  blockinfo *pbi;
  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext)
  {
    bbyte *pbStart = pbi->pb;
    bbyte *pbEnd   = pbi->pb + pbi->size - 1;

    if(fPtrGrtrEq(pb, pbStart) && fPtrLessEq(pb, pbEnd))
      break;
  }
  return (pbi);
}

/* ------------------ _CheckMemoryOn ------------------------ */

void _CheckMemoryOn(void){
  checkmemoryflag=1;
}

/* ------------------ _CheckMemoryOff ------------------------ */

void _CheckMemoryOff(void){
  checkmemoryflag=0;
}

/* ------------------ _CheckMemory ------------------------ */

void _CheckMemory(void){
  blockinfo *pbi;
  if(checkmemoryflag==0)return;
  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext)
  {
  if(sizeofDebugByte!=0)ASSERT((char)*(pbi->pb+pbi->size)==(char)debugByte);
  }
  return;
}

/* ------------------ CreateBlockInfo ------------------------ */

mallocflag CreateBlockInfo(bbyte *pbNew, size_t sizeNew){
  blockinfo *pbi;

  ASSERT(pbNew != NULL && sizeNew != 0);

  pbi = (blockinfo *)malloc(sizeof(blockinfo));
  if( pbi != NULL){
    pbi->pb = pbNew;
    pbi->size = sizeNew;
    pbi->pbiNext = pbiHead;
    pbiHead = pbi;
  }
  return (mallocflag)(pbi != NULL);
}

/* ------------------ FreeBlockIfno ------------------------ */

void FreeBlockInfo(bbyte *pbToFree){
  blockinfo *pbi, *pbiPrev;

  pbiPrev = NULL;
  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext){
    if(fPtrEqual(pbi->pb, pbToFree)){
      if(pbiPrev == NULL){
        pbiHead = pbi->pbiNext;
      }
      else{
        pbiPrev->pbiNext = pbi->pbiNext;
      }
      break;
    }
    pbiPrev = pbi;
  }
  ASSERT(pbi != NULL);
  if(sizeofDebugByte!=0)ASSERT((char)*(pbi->pb+pbi->size)==(char)debugByte);
  free(pbi);
}

/* ------------------ UpdateBlockInfo ------------------------ */

void UpdateBlockInfo(bbyte *pbOld, bbyte *pbNew, size_t sizeNew){
  blockinfo *pbi;

  ASSERT(pbNew != NULL && sizeNew != 0);

  pbi = GetBlockInfo(pbOld);
  ASSERT(pbOld == pbi->pb);

  pbi->pb = pbNew;
  pbi->size = sizeNew;
}

/* ------------------ sizeofBlock ------------------------ */

size_t sizeofBlock(bbyte *pb){
  blockinfo *pbi;

  pbi = GetBlockInfo(pb);
  ASSERT(pb==pbi->pb);
  if(sizeofDebugByte!=0)ASSERT((char)*(pbi->pb+pbi->size)==(char)debugByte);
  return(pbi->size);
}

/* ------------------ ValidPointer ------------------------ */

mallocflag ValidPointer(void *pv, size_t size){
  blockinfo *pbi;
  bbyte *pb = (bbyte *)pv;

  ASSERT(pv != NULL && size != 0);

  pbi = GetBlockInfo(pb);
  ASSERT(pb==pbi->pb);

  ASSERT(fPtrLessEq(pb+size,pbi->pb + pbi->size));

  if(sizeofDebugByte!=0)ASSERT((char)*(pbi->pb+pbi->size)==(char)debugByte);
  return(1);
}

/* ------------------ strcpy ------------------------ */

char *_strcpy(char *s1, const char *s2){
  blockinfo *pbi;
  int offset;
  CheckMemory;
  pbi = GetBlockInfo_nofail(s1);
  if(pbi!=NULL){
    offset = s1 - pbi->pb;
    ASSERT(pbi->size - offset >= strlen(s2)+1);
  }

  return strcpy(s1,s2);
}

/* ------------------ strcat ------------------------ */

char *_strcat(char *s1, const char *s2){
  blockinfo *pbi;
  int offset;
  CheckMemory;
  pbi = GetBlockInfo_nofail(s1);
  if(pbi!=NULL){
    offset = s1 - pbi->pb;
    ASSERT(pbi->size - offset >= strlen(s1)+strlen(s2)+1);
  }

  return strcat(s1,s2);
}
#endif
