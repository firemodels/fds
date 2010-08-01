// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#define INMAIN
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"

#define MARK 255

#define FORTREAD(read) fseek(BOUNDARYFILE,4,SEEK_CUR);returncode=read;fseek(BOUNDARYFILE,4,SEEK_CUR);

//dummy change to bump version number to 1.2.6

// svn revision character string
char main_revision[]="$Revision$";


/* ------------------ main ------------------------ */

int main(int argc, char **argv){

//  Bytef *source,*sourcecheck,*dest;
//  int sourceLen, destLen;
//  int returncode;

  char *arg;

  char *filebase;
  int filelen;
  char smvfile[1024];
  char smvfilebase[1024];
  char *ext;
  int addext=1;
  char inifile[1024];
  char inifilebase[1024];
  char *prog;
  int i;
  int endian_fds;
  int endian_info;

#ifdef pp_LIGHT
  nphotons=NPHOTONS;
#endif
  frameskip=-1;
  no_chop=0;
  autozip=0;
  make_demo=0;
  endf=0;
  syst=0;
  endianfile=NULL;
  destdir=NULL;
  sourcedir=NULL;
  lensourcedir=0;
  lendestdir=0;
  endianswitch=-1;
  overwrite_b=0;
  overwrite_s=0;
  overwrite_iso=0;
  get_bounds=0;
  get_slice_bounds=0;
  get_plot3d_bounds=0;
  get_boundary_bounds=0;
#ifdef pp_PART
  get_part_bounds=0;
#endif
#ifdef pp_LIGHT
  make_lighting_file=0;
  albedo=(float)0.3;
  nlightinfo=0;
  lightinfo=NULL;
  light_cdf=NULL;
  light_delta = 0.3;
  light_min=1.0;
  light_max=5.0;
#endif
  overwrite_slice=0;
  overwrite_plot3d=0;
  endian_info=0;
  cleanfiles=0;
  smoke3dzipstep=1;
  boundzipstep=1;
  slicezipstep=1;
  isozipstep=1;
  filesremoved=0;
#ifdef WIN32
  strcpy(dirseparator,"\\");
#else
  strcpy(dirseparator,"/");
#endif

  npatch_files=0;
  nsmoke3d_files=0;
#ifdef pp_PART
  npart_files=0;
  npartclassinfo=0;
  partinfo=NULL;
  partclassinfo=NULL;
  maxpart5propinfo=0;
  npart5propinfo=0;
#endif
  nslice_files=0;
  sliceinfo=NULL;
  nmeshes=0;
  niso_files=0;
  isoinfo=NULL;

  patchinfo=NULL;
  smoke3dinfo=NULL;
  strcpy(pp,"%");

  prog=argv[0];
  filebase=NULL;
  if(argc==1){
    version();
    return 1;
  }

  for(i=1;i<argc;i++){
    int lenarg;
    int lenarg2;
    char *arg2;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
#ifdef pp_LIGHT
      case 'l':
        make_lighting_file=1;
        break;
      case 'a':
        arg2=argv[i+1];
        if(strcmp(arg2,"auto")==0){
          autozip=1;
        }
        else{
          sscanf(arg2,"%f",&albedo);
          if(albedo<0.0)albedo=0.0;
          if(albedo>1.0)albedo=1.0;
          i++;
        }
        break;
#else
      case 'a':
        autozip=1;
        break;
#endif
      case 'b':
        if(strcmp(arg,"-bounds")==0){
          get_bounds=1;
          get_slice_bounds=1;
          get_plot3d_bounds=1;
          get_boundary_bounds=1;
#ifdef pp_PART
          get_part_bounds=1;
#endif
        }
        else if(strcmp(arg,"-bb")==0){
          get_boundary_bounds=1;
        }
        else if(strcmp(arg,"-bs")==0){
          get_slice_bounds=1;
        }
        else if(strcmp(arg,"-bp")==0){
          get_plot3d_bounds=1;
        }
#ifdef pp_PART
        else if(strcmp(arg,"-bP")==0){
          get_part_bounds=1;
        }
#endif
        else{
          overwrite_b=1;
        }
        break;
      case '2':
        overwrite_slice=1;
        break;
      case '3':
        overwrite_s=1;
        break;
#ifdef  pp_PART
      case 'P':
        overwrite_part=1;
        break;
#endif
      case 'p':
        overwrite_plot3d=1;
        break;
      case 'f':
        overwrite_b=1;
        overwrite_s=1;
        overwrite_slice=1;
        overwrite_iso=1;
        overwrite_plot3d=1;
#ifdef pp_PART
        overwrite_part=1;
#endif
        break;
      case 'i':
        overwrite_iso=1;
        break;
      case 'c':
        cleanfiles=1;
        break;
      case 'e':
        endian_info=1;
        break;
      case 'n':
        if(lenarg==2){
        }
        else if(strcmp(arg,"-no_chop")==0){
          no_chop=1;
        }
        i++;
        break;
      case 's':
        if(i+1>=argc)break;
        if(lenarg==2){
            lenarg2=strlen(argv[i+1]);
            NewMemory((void **)&sourcedir,lenarg2+2);
            strcpy(sourcedir,argv[i+1]);
            if(sourcedir[lenarg2-1]!=dirseparator[0]){
              strcat(sourcedir,dirseparator);
            }
            lensourcedir=strlen(sourcedir);
            if(getfileinfo(sourcedir,NULL,NULL)!=0){
              printf("The source directory specified, %s, does not exist or cannot be accessed\n",sourcedir);
              return 1;
            }
           i++;
        }
        else if(strcmp(arg,"-skip")==0){
          frameskip=-1;
          arg2=argv[i+1];
          sscanf(arg2,"%i",&frameskip);
          if(frameskip>0){
            slicezipstep=frameskip;
            isozipstep=frameskip;
            smoke3dzipstep=frameskip;
            boundzipstep=frameskip;
          }
          i++;
        }
        break;
      case 'd':
        if(strcmp(arg,"-demo")==0){
          autozip=1;
          make_demo=1;
          break;
        }
        if(i+1<argc){
          lenarg2=strlen(argv[i+1]);
          NewMemory((void **)&destdir,lenarg2+2);
          strcpy(destdir,argv[i+1]);
          if(destdir[lenarg2-1]!=dirseparator[0]){
            strcat(destdir,dirseparator);
          }
          lendestdir=strlen(destdir);
 //         if(getfileinfo(destdir,NULL,NULL)!=0){
 //           printf("The destination directory %s does not exist or cannot be accessed\n",destdir);
 //           return 1;
 //         }
          i++;
        }
        break;
      case 'h':
        usage(prog);
        return 1;
        break;
      case 'v':
        version();
        return 1;
        break;
      default:
        usage(prog);
        return 1;
      }
    }
    else{
      if(filebase==NULL){
        filebase=argv[i];
      }
    }
  }

  // construct smv filename
  
  if(filebase==NULL){
    usage(prog);
    return 1;
  }
  if(sourcedir==NULL){
    strcpy(smvfile,filebase);
  }
  else{
    strcpy(smvfile,sourcedir);
    strcat(smvfile,filebase);
  }
  strcpy(smvfilebase,filebase);
  filelen=strlen(filebase);
  if(filelen>4){
    ext=filebase+filelen-4;
    if(strcmp(ext,".smv")==0)addext=0;
  }
  if(addext==1)strcat(smvfile,".smv");
  
  // construct ini file name

  strcpy(inifile,smvfile);
  inifile[strlen(inifile)-4]=0;
  strcat(inifile,".ini");
  strcpy(inifilebase,filebase);
  strcat(inifilebase,".ini");

  strcpy(endianfilebase,"");

  // make sure smv file name exists

  if(getfileinfo(smvfile,NULL,NULL)!=0){
    printf("file: %s does not exist\n",smvfile);
    return 1;
  }

  // make sure smv file can be opened

  if(readsmv(smvfile)!=0)return 1;

  if(doiso==1&&niso_files>0){
    for(i=0;i<niso_files;i++){
      iso *isoi;

      isoi = isoinfo + i;
      if(isoi->blocknumber<0||isoi->blocknumber>=nmeshes){
        doiso=0;
        break;
      }
    }
  }
  if(nplot3d_files>0){
    plot3dinfo[0].dup=0;
    for(i=1;i<nplot3d_files;i++){
      plot3d *plot3di; 

      plot3di = plot3dinfo + i;

      plot3di->dup=0;
      plot3ddup(plot3di,i);
    }
  }
  if(npatch_files>0){
    patchinfo[0].dup=0;
    for(i=1;i<npatch_files;i++){
      patch *patchi; 

      patchi = patchinfo + i;

      patchi->dup=0;
      patchdup(patchi,i);
    }
  }
  if(nslice_files>0){
    sliceinfo[0].dup=0;
    for(i=1;i<nslice_files;i++){
      slice *slicei; 

      slicei = sliceinfo + i;

      slicei->dup=0;
      slicedup(slicei,i);
    }
  }

  if(getendian()==1){
      printf("Smokezip running on a big endian computer.\n");
  }
  else{
      printf("Smokezip running on a little endian computer.\n");
  }
  if(endf==0&&syst==0){
    printf("Warning: casename.end file is missing.  Endianness of\n");
    printf("         FDS boundary file data is unknown.\n");
    if(getendian()==1){
      printf("         Assuming FDS boundary data is big endian - \n");
    }
    if(getendian()==0){
      printf("         Assuming FDS boundary data is little endian - \n");
    }
    printf("         or equivalently assuming FDS and Smokezip are\n");
    printf("         being run on the same type of computer\n");
    endianswitch=0;
  }
  else{
    endian_fds=getendian()+endianswitch;
    if(endian_fds==2)endian_fds=0;
    if(endian_fds==1){
      printf("FDS was run on a big endian computer. \n");
    }
    else{
      printf("FDS was run on a little endian computer.\n");
    }
  }
  if(endian_info==1)return 0;

  readini(inifile);

  compress_patches();
  compress_smoke3ds();
  compress_slices();
  compress_plot3ds();
#ifdef pp_PART
  compress_parts();
#endif
  if(doiso==1)compress_isos();

  if(cleanfiles==0&&destdir!=NULL){
    printf("Copying .smv, .ini and .end files to %s directory\n",destdir);
    filecopy(destdir,smvfile,smvfilebase);
    filecopy(destdir,inifile,inifilebase);
    filecopy(destdir,endianfile,endianfilebase);
  }
  if(cleanfiles==1&&filesremoved==0){
    printf("No compressed files were removed\n");
  }
  if(make_demo==1){
    makesvd(destdir,smvfile);
  }

  return 0;
}
       
/* ------------------ filecopy ------------------------ */

#define SIZEBUFFER 1000000
void filecopy(char *destdir, char *file, char *filebase){
  char buffer[SIZEBUFFER];
  FILE *streamin;
  FILE *streamout;
  char *fileout=NULL;
  size_t chars_in;

  if(destdir==NULL||file==NULL)return;
  streamin=fopen(file,"rb");
  if(streamin==NULL)return;

  fileout=NULL;
  NewMemory((void **)&fileout,strlen(filebase)+strlen(destdir)+1);
  strcpy(fileout,destdir);
  strcat(fileout,filebase);

  streamout=fopen(fileout,"rb");
  if(streamout!=NULL){
    printf("  Warning: will not overwrite %s%s\n",destdir,file);
    fclose(streamout);
    fclose(streamin);
    return;
  }
  streamout=fopen(fileout,"wb");
  if(streamout==NULL){
    fclose(streamin);
    return;
  }
  printf("  Copying %s to %s\n",file,destdir);
  for(;;){
    int eof;
       
    eof=0;
    chars_in=fread(buffer,1,SIZEBUFFER,streamin);
    if(chars_in!=SIZEBUFFER)eof=1;
    if(chars_in>0)fwrite(buffer,chars_in,1,streamout);
    if(eof==1)break;
  }
}
       
       
/* ------------------ makesvd ------------------------ */

void makesvd(char *destdir, char *smvfile){
  char buffer[SIZEBUFFER];
  FILE *streamin;
  FILE *streamout;
  char *fileout=NULL;
  size_t chars_in;
  char *svd;

  if(smvfile==NULL)return;
  streamin=fopen(smvfile,"rb");
  if(streamin==NULL)return;

  fileout=NULL;
  if(destdir==NULL){
    NewMemory((void **)&fileout,strlen(smvfile)+2+1);
    strcpy(fileout,".");
    strcat(fileout,dirseparator);
  }
  else{
    NewMemory((void **)&fileout,strlen(smvfile)+strlen(destdir)+1);
    strcpy(fileout,destdir);
  }
  strcat(fileout,smvfile);
  svd = fileout + strlen(fileout) - 4;
  strcpy(svd,".svd");

  streamout=fopen(fileout,"wb");
  if(streamout==NULL){
    fclose(streamin);
    return;
  }
  printf("  Copying %s to %s\n",smvfile,fileout);
  for(;;){
    int eof;
       
    eof=0;
    chars_in=fread(buffer,1,SIZEBUFFER,streamin);
    if(chars_in!=SIZEBUFFER)eof=1;
    if(chars_in>0)fwrite(buffer,chars_in,1,streamout);
    if(eof==1)break;
  }
  fclose(streamout);
  fclose(streamin);
}
       
/* ------------------ usage ------------------------ */

void usage(char *prog){
  char pp[2];
  char smv_version[100];
  int svn_num;

  getSMZversion(smv_version);  // get Smokeview version (ie 5.x.z)
  svn_num=getmaxrevision();    // get svn revision number

  strcpy(pp,"%");
  printf("\n");
  printf("  smokezip %s(%i) - %s\n\n",smv_version,svn_num,__DATE__);
  printf("  Compresses Smokeview 3D smoke, slice, iso-surface and boundary files\n\n");
  printf("  %s",prog);
  printf(" [-c -f -3 -b -s -i -skip skipval -no_chop]");
#ifdef pp_LIGHT
  printf("[-a val -l]");
#endif
  printf("  [-c -d destdir -s sourcedir] casename\n\n");
  printf("  -2  - overwrites 2d slice compressed files\n");
  printf("  -3  - overwrites 3d smoke files\n");
  printf("  -b  - overwrites boundary compressed files\n");
  printf("  -i  - overwrites iso-surface compressed files\n");
  printf("  -p  - overwrites PLOT3D files\n");
#ifdef pp_PART
  printf("  -P  - overwrites particle files\n");
#endif
  printf("  -f  - overwrites all compressed files\n");
  printf("  -bounds - estimate data bounds for all file types\n");
  printf("  -bb - estimate data bounds for boundary files\n");
  printf("  -bs - estimate data bounds for slice files\n");
  printf("  -bp - estimate data bounds for plot3d files\n");
#ifdef pp_PART
  printf("  -bP - estimate data bounds for particle files\n");
#endif
  printf("  -d destdir - copies compressed files (and files needed by Smokeview\n");
  printf("               to view the case) to the directory destdir\n"); 
  printf("  -demo - Creates the files (compressed and .svd ) needed by the\n");
  printf("          Smokeview demonstrator mode.  Compresses files that are autoloaded, \n");
  printf("          uses (20.0,620.0) and (0.0,0.23) for temperature and oxygen bounds and\n");
  printf("          creates the .svd file which activates the Smokeview demonstrator mode.\n");
  printf("  -s sourcedir - specifies directory containing source files\n");
  printf("  -skip skipval - skip frames when compressing files\n");
  printf("  -no_chop - do not chop or truncate slice data.  Smokezip by default will compress\n");
  printf("             slice data truncating data above and below values specified in the .ini file\n");
  printf("  -auto - compress only files that are auto-loaded by Smokeview\n");
#ifdef pp_LIGHT
  printf("  -l  - create lighting file used with 3d smoke\n");
  printf("  -a albedo  - specifies albedo (0.0 to 1.0) of 3d smoke.\n");
  printf("               Used with the -l option\n");
#endif
  printf("  -c  - cleans or removes all compressed files\n");
  printf("  -h  - display this message\n\n");
  printf("  casename - Smokeview .smv file\n");
  printf("  Min and max bounds used to compress boundary files are obtained\n");
  printf("  from the casename.ini file or calculated by %s if casename.ini \n",prog);
  printf("  does not exist.  See http://fire.nist.gov/fds for more information.\n");
}
