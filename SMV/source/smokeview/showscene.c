// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "contourdefs.h"
#include "isodefs.h"

#include "flowfiles.h"
#include "smokeviewapi.h"
#include "MALLOC.h"
#include "smokeviewvars.h"
#include "viewports.h"
#include "update.h"

// svn revision character string
char showscene_revision[]="$Revision$";

/* ------------------ ShowScene ------------------------ */

void ShowScene(int mode, int view_mode, int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  CheckMemory;

  show_mode=mode;

  if(xyz_clipplane==1){
    if(clip_x==1)glDisable(GL_CLIP_PLANE0);
    if(clip_y==1)glDisable(GL_CLIP_PLANE1);
    if(clip_z==1)glDisable(GL_CLIP_PLANE2);
    if(clip_X==1)glDisable(GL_CLIP_PLANE3);
    if(clip_Y==1)glDisable(GL_CLIP_PLANE4);
    if(clip_Z==1)glDisable(GL_CLIP_PLANE5);
  }

/* ++++++++++++++++++++++++ update variables as needed +++++++++++++++++++++++++ */

  if(restart_time==1){
    restart_time=0;
    reset_itimes0();
  }
  if(loadfiles_at_startup&&update_load_startup==1){
    load_startup_smoke();
  }
//  if(updategluiview==1&&updateclipvals==0){
  if(updategluiview==1){
    camera *ca;

    ca = get_camera(label_startup_view);
    if(ca!=NULL){
      startup_view_ini = ca->view_id;
    }

    reset_glui_view(startup_view_ini);
    updategluiview=0;
  }
  if(menusmooth==1&&smoothing_blocks==0&&updatesmoothblocks==1){
    smooth_blockages();
  }
  if(update_tourlist==1){
    Update_Tourlist();
  }
  if(camera_current->dirty==1){
    update_camera(camera_current);
  }
  if(updateclipvals==1){
    clip2cam(camera_current);
    update_clip_all();
    updateclipvals=0;
    updategluiview=0;
  }
  if(update_selectedtour_index==1){
    update_tourindex();
  }
  if(trainer_mode==1&&fontindex!=LARGE_FONT)FontMenu(1);
  if(updateindexcolors==1){
    UpdateIndexColors();
  }
  if(force_isometric==1){
    force_isometric=0;
    projection_type=1;
    camera_current->projection_type=projection_type;
    ZoomMenu(-2);
  }

  updateShow();
  if(times!=NULL&&updateUpdateFrameRateMenu==1)sv_FrameRateMenu(frameratevalue);
  if(updatefaces==1)update_faces();
  if(updatefacelists==1)update_facelists();
  if(showstereo==0||showstereo==1)ClearBuffers(mode);

/* ++++++++++++++++++++++++ draw viewports +++++++++++++++++++++++++ */

  if(mode==RENDER){
    BLOCK_viewport(quad,          s_left,s_down,s_width,s_height);
    sniffErrors("after BLOCK_viewport");

    TIMEBAR_viewport(quad,          s_left,s_down,s_width,s_height);
    sniffErrors("after TIMEBAR_viewport");

    COLORBAR_viewport(quad,          s_left,s_down,s_width,s_height);
    sniffErrors("after COLORBAR_viewport");

    LOGO_viewport(quad,          s_left,s_down,s_width,s_height);
    sniffErrors("after LOGO_viewport");
    
    TITLE_viewport(quad,          s_left,s_down,s_width,s_height);
    sniffErrors("after TITLE_viewport");

    Scene_viewport(quad,view_mode,s_left,s_down,s_width,s_height);
    sniffErrors("after Scene_viewport");
  }

  

/* ++++++++++++++++++++++++ draw "fancy" colorbar +++++++++++++++++++++++++ */

  if(viscolorbarpath==1){
    if(cb_hidesv==1){
      setColorbarClipPlanes(0);
    }
    drawcolorbarpath();
    if(cb_hidesv==1){
      setColorbarClipPlanes(1);
    }
    sniffErrors("after setColorbarClipPlanes 1");
  }

  if(eyeview==1&&nskyboxinfo>0)draw_skybox();

  if(UpdateLIGHTS==1)updateLights(0);

  if(mode!=RENDER||viscolorbarpath!=1){
    setClipPlanes(0);
  }
  if(mode==RENDER){
    if(viscolorbarpath==1){
      if(cb_hidesv==1){
        setColorbarClipPlanes(1);
      }
      else{
        setColorbarClipPlanes(0);
      }
      sniffErrors("after setColorbarClipPlanes 2");
    }
    glPointSize((float)1.0);


    /* ++++++++++++++++++++++++ draw trees +++++++++++++++++++++++++ */

    if(ntreeinfo>0){
      drawtrees();
      sniffErrors("after drawtrees");
    }

/* ++++++++++++++++++++++++ draw particles +++++++++++++++++++++++++ */

    if(showsmoke==1){
      drawpart_frame();
    }

/* ++++++++++++++++++++++++ draw evacuation +++++++++++++++++++++++++ */

    if(showevac==1){
      drawevac_frame();
    }

/* ++++++++++++++++++++++++ draw targets +++++++++++++++++++++++++ */

    if(showtarget==1){
      drawTargets();
    }

/* ++++++++++++++++++++++++ draw sensors/sprinklers/heat detectors +++++++++++++++++++++++++ */

    if(xyz_clipplane==2){
      setClipPlanes(1);
    }
    draw_devices();
    if(xyz_clipplane==2){
      unsetClipPlanes();
    }
    sniffErrors("after draw_devices");

    if(visaxislabels==1||showedit==1){
      outputAxisLabels();
      sniffErrors("after outputAxisLables");
    }


 /* ++++++++++++++++++++++++ draw user ticks +++++++++++++++++++++++++ */

    if(vis_user_ticks==1){
      antialias(1);
      glDisable(GL_CLIP_PLANE0);
      glDisable(GL_CLIP_PLANE1);
      glDisable(GL_CLIP_PLANE2);
      glDisable(GL_CLIP_PLANE3);
      glDisable(GL_CLIP_PLANE4);
      glDisable(GL_CLIP_PLANE5);
      draw_user_ticks();
      if(mode!=RENDER||viscolorbarpath!=1){
        setClipPlanes(0);
      }
      antialias(0);
      sniffErrors("after drawticks");
    }

 /* ++++++++++++++++++++++++ draw ticks +++++++++++++++++++++++++ */

    if(visTicks==1&&nticks>0){
      drawticks();
      sniffErrors("after drawticks");
    }

    /* draw the box framing the simulation (corners at (0,0,0) (xbar,ybar,zbar) */


/* ++++++++++++++++++++++++ draw simulation frame (corners at (0,0,0) and (xbar,ybar,zbar) +++++++++++++++++++++++++ */

    if(isZoneFireModel==0&&visFrame==1&&highlight_flag==2){
      drawoutlines();
      sniffErrors("after drawoutlines");
    }


/* ++++++++++++++++++++++++ draw mesh +++++++++++++++++++++++++ */

    if(setPDIM==1){
      if(visGrid==1){
        int igrid;
        mesh *meshi;

        for(igrid=0;igrid<nmeshes;igrid++){
          meshi=meshinfo+igrid;
          drawgrid(meshi);
          sniffErrors("drawgrid");
        }
      }
    }
  } /* end of if(mode==RENDER) code segment */


/* ++++++++++++++++++++++++ draw selected devices +++++++++++++++++++++++++ */

  if(mode==SELECT){
    if(select_device==1){
     draw_devices();
      sniffErrors("after drawselect_devices");
      return;
    }
  }

/* ++++++++++++++++++++++++ draw selected avatars +++++++++++++++++++++++++ */

  if(mode==SELECT){
    if(select_avatar==1){
      drawselect_avatars();
      sniffErrors("after drawselect_avatars");
      return;
    }
  }

/* ++++++++++++++++++++++++ draw selected tours +++++++++++++++++++++++++ */

  if(mode==SELECT){
    if(edittour==1&&ntours>0){
      drawselect_tours();
      sniffErrors("after drawselect_tours");
      return;
    }
  }


/* ++++++++++++++++++++++++ draw tours +++++++++++++++++++++++++ */

  if(showtour==1){
    drawtours();
    sniffErrors("after drawTours");
  }

  /* ++++++++++++++++++++++++ draw stereo parallax indicator +++++++++++++++++++++++++ */
  
  if(show_parallax==1){
    antialias(1);
    glLineWidth(linewidth);
    glBegin(GL_LINES);
    glColor3fv(foregroundcolor);
    glVertex3f(0.75,0.0,0.25);
    glVertex3f(0.75,1.0,0.25);
    glEnd();
    antialias(0);
  }

  /* ++++++++++++++++++++++++ draw blockages +++++++++++++++++++++++++ */

  if(xyz_clipplane==2){
    setClipPlanes(1);
  }
  drawBlockages(mode,DRAW_OPAQUE);
  if(xyz_clipplane==2){
    unsetClipPlanes();
  }
  sniffErrors("drawBlockages");

  /* ++++++++++++++++++++++++ draw triangles +++++++++++++++++++++++++ */
  
  if(ngeominfo>0){
    draw_geom(DRAW_OPAQUE);
  }

#ifdef pp_SHOOTER
/* ++++++++++++++++++++++++ draw shooter points +++++++++++++++++++++++++ */

  if(showshooter!=0&&shooter_active==1){
    draw_shooter();
  }
#endif

/* ++++++++++++++++++++++++ draw terrain +++++++++++++++++++++++++ */

  if(visTerrainType!=4){
    int i;
    
    //shaded 17 0
    //stepped 18 1
    //line    19 2
    //texture 20 3
    //hidden 20 4

    for(i=0;i<nterraininfo;i++){
      terraindata *terri;
      int only_geom;

      terri = terraininfo + i;
      if(terri->loaded==1){
        only_geom=0;
      }
      else{
        only_geom=1;
      }
      switch (visTerrainType){
        case 1:
        case 2:
          if(cullfaces==1)glDisable(GL_CULL_FACE);
          DrawContours(&meshinfo[i].terrain_contour);
          if(cullfaces==1)glEnable(GL_CULL_FACE);
          break;
        case 0:
        case 3:
          if(visTerrainType==3&&terrain_texture!=NULL&&terrain_texture->loaded==1){
            drawterrain_texture(terri,only_geom);
          }
          else{
            drawterrain(terri,only_geom);
          }
          break;
      }
    }
  }

/* ++++++++++++++++++++++++ draw slice files +++++++++++++++++++++++++ */

  if(showslice==1&&use_transparency_data==0){
    drawslice_frame();
  } 

  /* ++++++++++++++++++++++++ draw boundary files +++++++++++++++++++++++++ */

  if(showpatch==1){
    drawpatch_frame();
  }

/* ++++++++++++++++++++++++ draw labels +++++++++++++++++++++++++ */

  if(visLabels==1){
    drawLabels();
  }

/* ++++++++++++++++++++++++ draw animated isosurfaces +++++++++++++++++++++++++ */

    //if(isoinfo!=NULL)drawspherepoints(sphereinfo);
  if(showiso==1){
    drawiso(DRAW_OPAQUE);
  }

/* ++++++++++++++++++++++++ draw zone fire modeling info +++++++++++++++++++++++++ */

  if(nrooms>0){
    drawroomgeom();
    sniffErrors("after drawroomgeom");
    if(showzone==1){
      drawfiredata();
      sniffErrors("after drawroomdata");
      if(ReadZoneFile==1&&nzvents>0){
        drawventdata();
        sniffErrors("after drawventdata");
      }
    }
  }


//**********************************************************************************
//**********************************************************************************
//**********************************************************************************
//    nothing transparent should be drawn before this portion of the code
//    (ie draw all opaque objects first then draw transparent objects
//**********************************************************************************
//**********************************************************************************
//**********************************************************************************

  /* ++++++++++++++++++++++++ draw triangles +++++++++++++++++++++++++ */
  
  if(ngeominfo>0){
    draw_geom(DRAW_TRANSPARENT);
  }

  if(showiso==1){
    drawiso(DRAW_TRANSPARENT);
  }

/* ++++++++++++++++++++++++ draw transparent faces +++++++++++++++++++++++++ */

  if(xyz_clipplane==2){
    setClipPlanes(1);
  }
  draw_transparent_faces();
  if(xyz_clipplane==2){
    unsetClipPlanes();
  }

/* ++++++++++++++++++++++++ draw 3D smoke +++++++++++++++++++++++++ */

  if(show3dsmoke==1||showvolrender==1){
    drawsmoke_frame();
  }

  if(active_smokesensors==1&&show_smokesensors!=0){
    getsmokesensors();
    draw_devices_val();
  }

/* ++++++++++++++++++++++++ draw zone fire modeling info +++++++++++++++++++++++++ */

  if(nrooms>0&&showzone==1){
    drawroomdata();
    sniffErrors("after drawroomdata");
  }

/* ++++++++++++++++++++++++ draw slice files +++++++++++++++++++++++++ */

  if(showslice==1&&use_transparency_data==1){
    drawslice_frame();
  } 
  sniffErrors("after drawslice");

/* ++++++++++++++++++++++++ draw transparent blockages +++++++++++++++++++++++++ */

//  draw_demo(20,20);
//  draw_demo2(1);
  drawBlockages(mode,DRAW_TRANSPARENT);
  sniffErrors("after drawBlokcages");

/* ++++++++++++++++++++++++ draw vector slice files +++++++++++++++++++++++++ */

  if(showvslice==1){
    drawvslice_frame();
  }
  sniffErrors("after drawvslice");

/* ++++++++++++++++++++++++ draw plot3d files +++++++++++++++++++++++++ */

  if(showplot3d==1){
    drawplot3d_frame();
  }
  sniffErrors("after drawplot3d");

  /* ++++++++++++++++++++++++ draw cross hairs +++++++++++++++++++++++++ */

#ifdef _DEBUG
    if(eyeview==EYE_CENTERED)drawMovedir();
#endif

/* ++++++++++++++++++++++++ render scene +++++++++++++++++++++++++ */

  Render(view_mode);

 /* ++++++++++++++++++++++++ draw "fancy" colorbar +++++++++++++++++++++++++ */

  if(viscolorbarpath==1){
    if(cb_hidesv==1){
      setColorbarClipPlanes(0);
    }
  }

  sniffErrors("end of loop");

}
