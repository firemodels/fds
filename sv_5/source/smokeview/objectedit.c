// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char objectedit_revision[]="$Revision$";

/* ------------------ movefblockage ------------------------ */

void movefblockage(float *xmin, float *xmax, float *ymin, float *ymax, float *zmin, float *zmax){
  blockagedata *bc;
  int imin, imax, jmin, jmax, kmin, kmax;
  int ibar,jbar,kbar;
  float *xplt_orig,*yplt_orig,*zplt_orig;
  ibar=current_mesh->ibar;
  jbar=current_mesh->jbar;
  kbar=current_mesh->kbar;
  xplt_orig=current_mesh->xplt_orig;
  yplt_orig=current_mesh->yplt_orig;
  zplt_orig=current_mesh->zplt_orig;

  bc=bchighlight;
  if(bc==NULL)return;
  imin = getxyzindex(*xmin,xplt_orig,ibar);
  imax = getxyzindex(*xmax,xplt_orig,ibar);
  if(imax<=imin){
    imax = imin;
    if(imax>ibar){
      imax=ibar;
      imin=ibar;
    }
  }
  *xmin = xplt_orig[imin];
  *xmax = xplt_orig[imax];

  jmin = getxyzindex(*ymin,yplt_orig,jbar);
  jmax = getxyzindex(*ymax,yplt_orig,jbar);
  if(jmax<=jmin){
    jmax = jmin;
    if(jmax>jbar){
      jmax=jbar;
      jmin=jbar;
    }
  }
  *ymin = yplt_orig[jmin];
  *ymax = yplt_orig[jmax];

  kmin = getxyzindex(*zmin,zplt_orig,kbar);
  kmax = getxyzindex(*zmax,zplt_orig,kbar);
  if(kmax<=kmin){
    kmax = kmin;
    if(kmax>kbar){
      kmax=kbar;
      kmin=kbar;
    }
  }
  *zmin = zplt_orig[kmin];
  *zmax = zplt_orig[kmax];

  bc->ijk[IMIN]=imin;
  bc->ijk[IMAX]=imax;
  bc->ijk[JMIN]=jmin;
  bc->ijk[JMAX]=jmax;
  bc->ijk[KMIN]=kmin;
  bc->ijk[KMAX]=kmax;
  bc->changed=1;
  if(bc->id>0&&bc->id<=nchanged_idlist){
    changed_idlist[bc->id]=1;
  }
  blockages_dirty=1;

}

/* ------------------ moveiblockage ------------------------ */

void moveiblockage(int ival){
  blockagedata *bc;
  int imin, imax, jmin, jmax, kmin, kmax;
  int diff;
  int ibar,jbar,kbar;

  ibar=current_mesh->ibar;
  jbar=current_mesh->jbar;
  kbar=current_mesh->kbar;

  bc=bchighlight;
  if(bc==NULL)return;
  switch (xyz_dir){
  case 0:
    imin = bc->ijk[IMIN];
    diff = bc->ijk[IMAX] - bc->ijk[IMIN];
    imax = ival + diff;
    if(imax>ibar){
      imax = ibar;
    }
    imin = imax - diff;
    bc->ijk[IMIN] = imin;
    bc->ijk[IMAX] = imax;
    break;
  case 1:
    jmin = bc->ijk[JMIN];
    diff = bc->ijk[JMAX] - bc->ijk[JMIN];
    jmax = ival + diff;
    if(jmax>jbar){
      jmax = jbar;
    }
    jmin = jmax - diff;
    bc->ijk[JMIN] = jmin;
    bc->ijk[JMAX] = jmax;
    break;
  case 2:
    kmin = bc->ijk[KMIN];
    diff = bc->ijk[KMAX] - bc->ijk[KMIN];
    kmax = ival + diff;
    if(kmax>kbar){
      kmax = kbar;
    }
    kmin = kmax - diff;
    bc->ijk[KMIN] = kmin;
    bc->ijk[KMAX] = kmax;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  bc->changed=1;
  if(bc->id>0&&bc->id<=nchanged_idlist){
    changed_idlist[bc->id]=1;
  }
  blockages_dirty=1;


}


/* ------------------ stretchiblockage ------------------------ */

void stretchiblockage(int ival){
  blockagedata *bc;

  bc = bchighlight;
  if(bc==NULL)return;
  switch (xyz_dir){
  case 0:
    if(which_face==-1){
      bc->ijk[IMIN] = ival;
      if(bc->ijk[IMIN]>=bc->ijk[IMAX])bc->ijk[IMIN]=bc->ijk[IMAX];
    }
    if(which_face==1){
      bc->ijk[IMAX] = ival;
      if(bc->ijk[IMAX]<=bc->ijk[IMIN])bc->ijk[IMAX]=bc->ijk[IMIN];
    }
    break;
  case 1:
    if(which_face==-1){
      bc->ijk[JMIN] = ival;
      if(bc->ijk[JMIN]>=bc->ijk[JMAX])bc->ijk[JMIN]=bc->ijk[JMAX];
    }
    if(which_face==1){
      bc->ijk[JMAX] = ival;
      if(bc->ijk[JMAX]<=bc->ijk[JMIN])bc->ijk[JMAX]=bc->ijk[JMIN];
    }
    break;
  case 2:
    if(which_face==-1){
      bc->ijk[KMIN] = ival;
      if(bc->ijk[KMIN]>=bc->ijk[KMAX])bc->ijk[KMIN]=bc->ijk[KMAX];
    }
    if(which_face==1){
      bc->ijk[KMAX] = ival;
      if(bc->ijk[KMAX]<=bc->ijk[KMIN])bc->ijk[KMAX]=bc->ijk[KMIN];
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  bc->changed=1;
  if(bc->id>0&&bc->id<=nchanged_idlist){
    changed_idlist[bc->id]=1;
  }
  blockages_dirty=1;

}

/* ------------------ getxyzindex ------------------------ */

int getxyzindex(float x,const float *xplt_orig,int ibar){
  int i;
  if(x<xplt_orig[0])return 0;
  if(x>=xplt_orig[ibar])return ibar;
  for(i=0;i<ibar;i++){
    if(xplt_orig[i]<=x&&x<xplt_orig[i+1])return i;
  }
  return 0;
}

/* ------------------ outputblockage ------------------------ */

void outputblockage(void){
  blockagedata *bc;
  float *xplttemp,*yplttemp,*zplttemp;

  xplttemp=current_mesh->xplt;
  yplttemp=current_mesh->yplt;
  zplttemp=current_mesh->zplt;

  bc = bchighlight;
  printf("&OBST XB=%f,%f,%f,%f,%f,%f / \n",
    xbar0+xyzmaxdiff*xplttemp[bc->ijk[IMIN]],xbar0+xyzmaxdiff*xplttemp[bc->ijk[IMAX]],
    ybar0+xyzmaxdiff*yplttemp[bc->ijk[JMIN]],ybar0+xyzmaxdiff*yplttemp[bc->ijk[JMAX]],
    zbar0+xyzmaxdiff*zplttemp[bc->ijk[KMIN]],zbar0+xyzmaxdiff*zplttemp[bc->ijk[KMAX]]
    );
}

/* ------------------ addnewobject ------------------------ */

void addnewobject(void){
  blockagedata *bc;
  int return_code1;
  mesh *meshi;
  int meshindex;

  meshi=current_mesh;
  meshindex = meshi - selected_case->meshinfo;
  meshi->nbptrs++;
  return_code1=NewMemory((void **)&bc,sizeof(blockagedata));
  if(meshi->nbptrs==1){
    NewMemory((void **)&meshi->blockageinfoptrs,sizeof(blockagedata *)*meshi->nbptrs);
    NewMemory((void **)&meshi->deletelist,sizeof(blockagedata *)*meshi->nbptrs);
  }
  else{
    ResizeMemory((void **)&meshi->blockageinfoptrs,sizeof(blockagedata *)*meshi->nbptrs);
    ResizeMemory((void **)&meshi->deletelist,sizeof(blockagedata *)*meshi->nbptrs);
  }
  meshi->blockageinfoptrs[meshi->nbptrs-1]=bc;

  if(return_code1==0){
    printf("*** Warning: unable to allocate needed memory in addnewobject\n");
    return;
  }
  initobst(bc,surfacedefault,-1,meshindex);

  bc->xmin=meshi->xplt_orig[0];
  bc->xmax=meshi->xplt_orig[1];
  bc->ymin=meshi->yplt_orig[0];
  bc->ymax=meshi->yplt_orig[1];
  bc->zmin=meshi->zplt_orig[0];
  bc->zmax=meshi->zplt_orig[1];

  bc->xmin += meshi->offset[0];
  bc->xmax += meshi->offset[0];
  bc->ymin += meshi->offset[1];
  bc->ymax += meshi->offset[1];
  bc->zmin += meshi->offset[2];
  bc->zmax += meshi->offset[2];

  bc->xmin = (bc->xmin-xbar0)/xyzmaxdiff;
  bc->xmax = (bc->xmax-xbar0)/xyzmaxdiff;
  bc->ymin = (bc->ymin-ybar0)/xyzmaxdiff;
  bc->ymax = (bc->ymax-ybar0)/xyzmaxdiff;
  bc->zmin = (bc->zmin-zbar0)/xyzmaxdiff;
  bc->zmax = (bc->zmax-zbar0)/xyzmaxdiff;

  bc->walltype=WALL_1;
  bc->walltypeORIG=WALL_1;
  backup_blockage(bc);
  CheckMemory;
  update_faces();
  highlight_block=meshi->nbptrs-1;
  highlight_mesh=meshi-selected_case->meshinfo;
  update_highlight_mesh();
  xyz_dir=0;
  which_face=0;
  update_xyzdir(xyz_dir);
  bchighlight_old=bchighlight;
  bchighlight=bc;
  update_blockvals(1);
}

/* ------------------ obstlabelcopy ------------------------ */

void obstlabelcopy(char **dest, const char *source){
  size_t len;

  len=0;
  if(source!=NULL){
    len = strlen(source);
    FREEMEMORY(*dest);
    if(len>0){
      NewMemory((void **)dest,(unsigned int)(len+1));
      strcpy(*dest,source);
    }
  }
  if(source==NULL||len==0){
    NewMemory((void **)dest,2);
    strcpy(*dest,"");
  }
}

/* ------------------ deleteblockage ------------------------ */

void deleteobject(void){
  int ndelete;

  if(bchighlight!=NULL){
    ndelete=current_mesh->ndeletelist;
    current_mesh->deletelist[ndelete]=bchighlight;
    current_mesh->ndeletelist++;
    bchighlight->del=1;
    update_faces();
    highlight_block=-1;
    xyz_dir=0;
    which_face=0;
    bchighlight_old=bchighlight;
    bchighlight=NULL;
    update_xyzdir(xyz_dir);
    update_blockvals(1);
  }
}

/* ------------------ undeleteblockage ------------------------ */

void undeleteobject(void){
  int ndelete;
  int i;
  mesh *meshi;

  meshi = current_mesh;
  if(meshi->ndeletelist>0){
    meshi->ndeletelist--;
    ndelete=meshi->ndeletelist;
    bchighlight_old=bchighlight;
    bchighlight=meshi->deletelist[ndelete];
    bchighlight->del=0;

    update_faces();
    for(i=0;i<meshi->nbptrs;i++){
      if(bchighlight==meshi->blockageinfoptrs[i]){
        highlight_block=i;
        break;
      }
    }
    xyz_dir=0;
    which_face=0;
    update_xyzdir(xyz_dir);
    update_blockvals(1);
  }
  return;
}
