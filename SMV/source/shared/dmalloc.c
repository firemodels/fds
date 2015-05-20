#include "options.h"
#define INDMALLOC
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "MALLOC.h"
#ifdef pp_MEMDEBUG
static int checkmemoryflag=1;
#endif
#ifdef WIN32
#include <windows.h>
#endif

#ifdef pp_MEMDEBUG
static blockinfo *GetBlockInfo(bbyte *pb);
#endif

#ifdef WIN32

// return memory usage between 0% and 100%

int memusage(void){
  MEMORYSTATUS stat;
  int load;

  GlobalMemoryStatus(&stat);    
  load=stat.dwMemoryLoad;
  return load;
}
#else

// return memory usage between 0% and 100%

int memusage(void){
  return 0;
}
#endif

/* ------------------ _memorystatus ------------------------ */

#ifdef pp_memstatus
#ifdef WIN32
void _memorystatus(unsigned int size,unsigned int *availmem,unsigned int *physmemused, unsigned int *totalmem){
  MEMORYSTATUS stat;

    GlobalMemoryStatus(&stat);
    if(availmem!=NULL)*availmem=stat.dwMemoryLoad;
    if(totalmem!=NULL)*totalmem=stat.dwTotalPhys/(1024*1024);
    if(physmemused!=NULL)*physmemused=(stat.dwTotalPhys-stat.dwAvailPhys)/(1024*1024);
#ifdef pp_MEMDEBUG
    if(size!=0&&size<=stat.dwAvailPhys-0.1*stat.dwTotalPhys){
      int memsize;

      memsize = stat.dwAvailPhys/(1024*1024);
      fprintf(stderr,"*** Available Memory: %i M \n",memsize);
    }
#endif
    if(size!=0&&size>stat.dwAvailPhys-0.1*stat.dwTotalPhys){
      fprintf(stderr,"*** Warning: Low Memory. Only %i M available for viewing data.\n",
           (int)stat.dwAvailPhys/(1024*1024));
      fprintf(stderr,"    Unload datafiles or system performance may degrade.\n");
    }
}
#endif
#endif

/* ------------------ initMM ------------------------ */

void initMALLOC(void){
  
  MMfirstptr=&MMfirst;
  MMlastptr=&MMlast;

  MMfirstptr->prev=NULL;
  MMfirstptr->next=MMlastptr;

  MMlastptr->prev=MMfirstptr;
  MMlastptr->next=NULL;
#ifdef pp_THREAD
  pthread_mutex_init(&mutexMEM,NULL);
#endif
#ifdef pp_MEMDEBUG
  MMmaxmemory=0;
  MMtotalmemory=0;
#endif

}

/* ------------------ _NewMemory ------------------------ */

mallocflag _NewMemory(void **ppv, size_t size, int memory_id, char *varname, char *file, int linenumber){
  mallocflag returnval;

  LOCK_MEM;
  returnval=_NewMemoryNOTHREAD(ppv, size, memory_id);
  if(returnval!=1){
    fprintf(stderr,"*** Error: memory allocation request of size %llu failed\n",(unsigned long long)size);
    if(varname!=NULL){
      fprintf(stderr,"             variable: %s\n",varname);
      fprintf(stderr,"                 size: %u\n",(unsigned int)size);
    }
    if(file!=NULL){
      fprintf(stderr,"                 file: %s\n",file);
      fprintf(stderr,"          line number: %i\n",linenumber);
    }
  }
  UNLOCK_MEM;
  return returnval;
}

/* ------------------ _NewMemoryNOTHREAD ------------------------ */

mallocflag _NewMemoryNOTHREAD(void **ppv, size_t size, int memory_id){
  void **ppb=(void **)ppv;
#ifdef pp_MEMDEBUG
  char *c;
#endif
  int infoblocksize;
  MMdata *this_ptr, *prev_ptr, *next_ptr;

  ASSERT(ppv != NULL && size != 0);
  infoblocksize=(sizeof(MMdata)+3)/4;
  infoblocksize*=4;

#ifdef pp_MEMDEBUG
  if(MMmaxmemory==0||MMtotalmemory+size<=MMmaxmemory){
    this_ptr = (void *)malloc(infoblocksize+size+sizeofDebugByte);
  }
  else{
    this_ptr = NULL;
  }
#else
  this_ptr = (void *)malloc(infoblocksize+size+sizeofDebugByte);
#endif
  if(this_ptr!=NULL){
    prev_ptr=MMfirstptr;
    next_ptr=MMfirstptr->next;

    prev_ptr->next=this_ptr;
    next_ptr->prev=this_ptr;

#ifdef pp_MEMPRINT
    this_ptr->size = size;
#endif
    this_ptr->memory_id = memory_id;
    this_ptr->prev=prev_ptr;
    this_ptr->next=next_ptr;
    this_ptr->marker=markerByte;

    *ppb=(char *)this_ptr+infoblocksize;
  }
  else{
    *ppb=NULL;
  }

#ifdef pp_MEMDEBUG
  {
    CheckMemoryNOTHREAD;
    if(*ppb != NULL){
      if(sizeofDebugByte!=0){
       c = (char *)(*ppb) + size;
       *c=(char)debugByte;
      }
      memset(*ppb, memGarbage, size);
      if(!CreateBlockInfo(*ppb, size)){
        free((char *)*ppb-infoblocksize);
        *ppb=NULL;
      }
    }
    MMtotalmemory+=size;
    ASSERT(*ppb !=NULL);
  }
#endif
  return (*ppb != NULL);
}

/* ------------------ FreeAllMemory ------------------------ */

void FreeAllMemory(int memory_id){
  MMdata *thisptr, *nextptr;
  int infoblocksize;
#ifdef _DEBUG
  int count = 0, count2 = 0;
  int nblocks = 0;
#endif

  LOCK_MEM;
  infoblocksize=(sizeof(MMdata)+3)/4;
  infoblocksize*=4;

#ifdef _DEBUG
  thisptr = MMfirstptr->next;
  for(;;){
    // if the 'thisptr' memory block is freed then thisptr is no longer valid.
    // so, nextptr (which is thisptr->next) must be defined before it is freed
    nextptr = thisptr->next;
    if(thisptr->next == NULL || thisptr->marker != markerByte)break;
    nblocks++;
    thisptr = nextptr;
  }
#endif

  thisptr = MMfirstptr->next;
  for(;;){
    // if the 'thisptr' memory block is freed then thisptr is no longer valid.
    // so, nextptr (which is thisptr->next) must be defined before it is freed
    nextptr = thisptr->next;
    if(thisptr->next == NULL || thisptr->marker != markerByte)break;
    if(memory_id == 0 || thisptr->memory_id == memory_id){
      FreeMemoryNOTHREAD((char *)thisptr + infoblocksize);
#ifdef _DEBUG
      count2++;
#endif
    }
#ifdef _DEBUG
    count++;
    if(count % 1000 == 0)printf("unloading %i blocks out of %i %i\n", count2, count,nblocks);
#endif
    thisptr = nextptr;
  }
  UNLOCK_MEM;
}

#ifdef pp_MEMPRINT
/* ------------------ _PrintMemoryInfo ------------------------ */

void _PrintMemoryInfo(void){
  MMdata *thisptr;
  int n = 0;
  LINT size = 0;

  for(thisptr = MMfirstptr->next;thisptr->next!=NULL;thisptr=thisptr->next){
    size += thisptr->size;
    n++;
  }
  PRINTF("nblocks=%i sizeblocks=%llu\n", n, size);
}
#endif

/* ------------------ FreeMemory ------------------------ */

void FreeMemory(void *pv){
  LOCK_MEM;
  FreeMemoryNOTHREAD(pv);
  UNLOCK_MEM;
}

/* ------------------ FreeMemoryNOTHREAD ------------------------ */

void FreeMemoryNOTHREAD(void *pv){
#ifdef pp_MEMDEBUG
  int len_memory;
#endif
  int infoblocksize;
  MMdata *this_ptr, *prev_ptr, *next_ptr;

  ASSERT(pv != NULL);
  infoblocksize=(sizeof(MMdata)+3)/4;
  infoblocksize*=4;
#ifdef pp_MEMDEBUG
  {
    blockinfo *meminfoblock;

    CheckMemoryNOTHREAD;
    meminfoblock = GetBlockInfo(pv);
    MMtotalmemory-=meminfoblock->size;
    len_memory=sizeofBlock((char *)pv);
    memset((char *)pv, memGarbage, len_memory);
    FreeBlockInfo((char *)pv);
  }
#endif
  this_ptr=(MMdata *)((char *)pv-infoblocksize);
  ASSERT(this_ptr->marker==markerByte);
  prev_ptr=this_ptr->prev;
  next_ptr=this_ptr->next;

  prev_ptr->next=next_ptr;
  next_ptr->prev=prev_ptr;
  free((char *)pv-infoblocksize);
}

/* ------------------ _ResizeMemory ------------------------ */

mallocflag _ResizeMemory(void **ppv, size_t sizeNew, int memory_id, char *varname, char *file, int linenumber){
  mallocflag returnval;

  LOCK_MEM;
  returnval=_ResizeMemoryNOTHREAD(ppv, sizeNew, memory_id);
  if(returnval!=1){
    fprintf(stderr,"*** Error: memory allocation request of size %llu failed\n",(unsigned long long)sizeNew);
    if(varname!=NULL){
      fprintf(stderr,"             variable: %s\n",varname);
      fprintf(stderr,"                 size: %u\n",(unsigned int)sizeNew);
    }
    if(file!=NULL){
      fprintf(stderr,"                 file: %s\n",file);
      fprintf(stderr,"          line number: %i\n",linenumber);
    }
  }
  UNLOCK_MEM;
  return returnval;
}

/* ------------------ _ResizeMemoryNOTHREAD ------------------------ */

mallocflag _ResizeMemoryNOTHREAD(void **ppv, size_t sizeNew, int memory_id){
  bbyte **ppold, *pbNew;
  int infoblocksize;
  MMdata *this_ptr, *prev_ptr, *next_ptr;
#ifdef pp_MEMDEBUG
  char *c;
  size_t sizeOld;
#endif

  infoblocksize=(sizeof(MMdata)+3)/4;
  infoblocksize*=4;

  ppold=(bbyte **)ppv;
  ASSERT(ppold != NULL && sizeNew != 0);
#ifdef pp_MEMDEBUG
  {
    CheckMemoryNOTHREAD;
    sizeOld = sizeofBlock(*ppold);
    if(sizeNew<sizeOld){
      memset((*ppold)+sizeNew,memGarbage,sizeOld-sizeNew);
    }
    else if(sizeNew > sizeOld){
      void *pbForceNew;

      if(_NewMemoryNOTHREAD((void **)&pbForceNew, sizeNew, 0)){
        memcpy(pbForceNew, *ppold, sizeOld);
        FreeMemoryNOTHREAD(*ppold);
        *ppold = pbForceNew;
      }
    }
  }
#endif

  this_ptr=(MMdata *)((char *)(*ppold)-infoblocksize);
  prev_ptr=this_ptr->prev;
  next_ptr=this_ptr->next;
  pbNew = realloc((char *)(*ppold)-infoblocksize, infoblocksize+sizeNew+sizeofDebugByte);
  if(pbNew != NULL){
    if(pbNew!=(char *)(*ppold)-infoblocksize){
      this_ptr=(MMdata *)pbNew;
      
      prev_ptr->next=this_ptr;
      next_ptr->prev=this_ptr;
      
#ifdef pp_MEMPRINT
      this_ptr->size = sizeNew;
#endif
      this_ptr->memory_id = memory_id;
      this_ptr->next=next_ptr;
      this_ptr->prev=prev_ptr;
      this_ptr->marker=markerByte;
    }
#ifdef pp_MEMDEBUG
    {
      if(sizeofDebugByte!=0){
        c = pbNew + infoblocksize + sizeNew;
        *c=(char)debugByte;
      }
      UpdateBlockInfo(*ppold, (char *)pbNew+infoblocksize, sizeNew);
      if(sizeNew>sizeOld){
        memset(pbNew+infoblocksize+sizeOld,memGarbage,sizeNew-sizeOld);
      }
    }
#endif
    *ppold = pbNew+infoblocksize;
  }
  return (pbNew != NULL);
}

#ifdef pp_MEMDEBUG
/* ------------------ pointer comparison defines ------------------------ */

#define fPtrEqual(pLeft, pRight)  ((pLeft) == (pRight))
#define fPtrLessEq(pLeft, pRight) ((pLeft) <= (pRight))
#define fPtrGrtrEq(pLeft, pRight) ((pLeft) >= (pRight))

/* ------------------ __NewMemory ------------------------ */

mallocflag __NewMemory(void **ppv, size_t size, int memory_id, char *varname, char *file, int linenumber){
  void **ppb=(void **)ppv;
  blockinfo *pbi;
  int return_code;
  char *varname2;
  char *file2;
  char ampersand='&';
#ifdef WIN32
  char dirsep='\\';
#else
  char dirsep='/';
#endif

  LOCK_MEM;
  return_code=_NewMemoryNOTHREAD(ppb,size,memory_id);
  pbi=GetBlockInfo((bbyte *)*ppb);
  pbi->linenumber=linenumber;

  file2=strrchr(file,dirsep);
  if(file2==NULL){
    file2=file;
  }
  else{
    file2++;
  }
  if(strlen(file2)<256){
    strcpy(pbi->filename,file2);
  }
  else{
    strncpy(pbi->filename,file2,255);
    strcat(pbi->filename,"\0");
  }

  varname2 = strchr(varname,ampersand);
  if(varname2==NULL){
    varname2=varname;
  }
  else{
    varname2++;
  }
  if(strlen(varname2)<256){
    strcpy(pbi->varname,varname2);
  }
  else{
    strncpy(pbi->varname,varname2,255);
    strcat(pbi->varname,"\0");
  }
  if(return_code!=1){
    fprintf(stderr,"*** Error: memory allocation request of size %llu failed\n",(unsigned long long)size);
    if(varname!=NULL){
      fprintf(stderr,"             variable: %s\n",varname);
      fprintf(stderr,"                 size: %u\n",(unsigned int)size);
    }
    if(file!=NULL){
      fprintf(stderr,"                 file: %s\n",file);
      fprintf(stderr,"          line number: %i\n",linenumber);
    }
  }
  UNLOCK_MEM;
  return return_code;
}

/* ------------------ __ResizeMemory ------------------------ */

mallocflag __ResizeMemory(void **ppv, size_t size, int memory_id, char *varname, char *file, int linenumber){
  void **ppb=(void **)ppv;
  blockinfo *pbi;
  int return_code;

  LOCK_MEM;
  return_code=_ResizeMemoryNOTHREAD(ppb,size,memory_id);
  pbi=GetBlockInfo((bbyte *)*ppb);
  pbi->linenumber=linenumber;
  if(strlen(file)<256){
    strcpy(pbi->filename,file);
  }
  else{
    strncpy(pbi->filename,file,255);
    strcat(pbi->filename,"\0");
  }
  if(strlen(varname)<256){
    strcpy(pbi->varname,varname);
  }
  else{
    strncpy(pbi->varname,varname,255);
    strcat(pbi->varname,"\0");
  }
  if(return_code!=1){
    fprintf(stderr,"*** Error: memory allocation request of size %llu failed\n",(unsigned long long)size);
    if(varname!=NULL){
      fprintf(stderr,"             variable: %s\n",varname);
      fprintf(stderr,"                 size: %u\n",(unsigned int)size);
    }
    if(file!=NULL){
      fprintf(stderr,"                 file: %s\n",file);
      fprintf(stderr,"          line number: %i\n",linenumber);
    }
  }
  UNLOCK_MEM;
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

    if(fPtrGrtrEq(pb, pbStart) && fPtrLessEq(pb, pbEnd))break;
  }
  ASSERT(pbi != NULL);
  return (pbi);
}

/* ------------------ CountMemoryBlocks ------------------------ */

int _CountMemoryBlocks(void){
  blockinfo *pbi;
  int n=0;

  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext){
    n++;
  }
  return n;
}

/* ------------------ GetTotalMemory ------------------------ */

MMsize _GetTotalMemory(void){
  return MMtotalmemory;
}

/* ------------------ PrintAllMemoryInfo ------------------------ */

void _PrintAllMemoryInfo(void){
  blockinfo *pbi;
  int n=0,size=0;

  PRINTF("\n\n");
  PRINTF("********************************************\n");
  PRINTF("********************************************\n");
  PRINTF("********************************************\n");
  for (pbi = pbiHead; pbi != NULL; pbi = pbi->pbiNext)
  {
    n++;
    size += pbi->size;
    PRINTF("%s allocated in %s at line %i\n",pbi->varname,pbi->filename,pbi->linenumber);
  }
  PRINTF("nblocks=%i sizeblocks=%i\n",n,size);
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
  LOCK_MEM;
  _CheckMemoryNOTHREAD();
  UNLOCK_MEM;
}

/* ------------------ _CheckMemory ------------------------ */

void _CheckMemoryNOTHREAD(void){
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

mallocflag _ValidPointer(void *pv, size_t size){
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

  LOCK_MEM;
  CheckMemoryNOTHREAD;
  pbi = GetBlockInfo_nofail(s1);
  if(pbi!=NULL){
    offset = s1 - pbi->pb;
    ASSERT(pbi->size - offset >= strlen(s2)+1);
  }
  UNLOCK_MEM;

  return strcpy(s1,s2);
}

/* ------------------ strcat ------------------------ */

char *_strcat(char *s1, const char *s2){
  blockinfo *pbi;
  int offset;

  LOCK_MEM;
  CheckMemoryNOTHREAD;
  pbi = GetBlockInfo_nofail(s1);
  if(pbi!=NULL){
    offset = s1 - pbi->pb;
    ASSERT(pbi->size - offset >= strlen(s1)+strlen(s2)+1);
  }
  UNLOCK_MEM;

  return strcat(s1,s2);
}
#endif
#ifdef pp_MEMDEBUG

/* ------------------ set_memcheck ------------------------ */

void set_memcheck(int index){
  switch(index){
  case 0:
    MMmaxmemory=0;
    break;
  case 1:
    MMmaxmemory=1000000000;
    break;
  case 2:
    MMmaxmemory=2000000000;
    break;
#ifdef BIT64
  case 3:
    MMmaxmemory=4000000000;
    break;
  case 4:
    MMmaxmemory=8000000000;
    break;
  default:
    ASSERT(0);
    break;
#endif
  }
}

/* ------------------ getMemusage ------------------------ */

void getMemusage(MMsize totalmemory,char *MEMlabel){
  int size;
  float rsize;

  if(totalmemory<1000000000){
    size = totalmemory/1000000;
    sprintf(MEMlabel,"%i MB",size);
  }
  else{
    rsize = totalmemory/1000000000.0;
    sprintf(MEMlabel,"%4.2f GB",rsize);
  }

}
#endif

