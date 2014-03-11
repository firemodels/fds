.SUFFIXES: .c .cpp .o .obj

SIZE=-m64

PLATFORM=

GLUT_INC_LOCATION=../include2

CFLAGS = -D pp_NOCLOSE

RM=/bin/rm
COMPILER=icc
COMPILER2=icc
CC=${COMPILER}
CPP=${COMPILER2}
CPPFLAGS=$(PLATFORM) -I${GLUT_INC_LOCATION} -I ../include2/glut ${SIZE}
LPATH=${GLUTPATH} 
INCS= 

.c.obj:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCS) $(FILTERC) $<
.c.o:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCS) $(FILTERC) $<
.cpp.obj:
	$(CPP) -c $(CFLAGS) $(CPPFLAGS) $(INCS) $(FILTERC) $<
.cpp.o:
	$(CPP) -c $(CFLAGS) $(CPPFLAGS) $(INCS) $(FILTERC) $<


LIBGLUI= -Llib -lglui

All: libglui.a 

GLUI_OBJS = glui_img_checkbox_0.o glui_img_checkbox_0_dis.o glui_img_checkbox_1.o glui_img_checkbox_1_dis.o\
glui_img_downarrow.o glui_img_leftarrow.o glui_img_listbox_down.o glui_img_listbox_up.o glui_img_listbox_up_dis.o\
glui_img_radiobutton_0.o glui_img_radiobutton_0_dis.o glui_img_radiobutton_1.o glui_img_radiobutton_1_dis.o \
glui_img_rightarrow.o glui_img_spindown_0.o glui_img_spindown_1.o glui_img_spindown_dis.o glui_img_spinup_0.o\
glui_img_spinup_1.o glui_img_spinup_dis.o glui_img_uparrow.o algebra3.o arcball.o glui.o glui_add_controls.o\
glui_bitmaps.o glui_bitmap_img_data.o glui_button.o glui_checkbox.o glui_column.o glui_commandline.o glui_control.o\
glui_edittext.o glui_filebrowser.o glui_list.o glui_listbox.o glui_mouse_iaction.o glui_node.o glui_panel.o \
glui_radio.o glui_rollout.o glui_rotation.o glui_scrollbar.o glui_separator.o glui_spinner.o glui_statictext.o\
glui_string.o glui_textbox.o glui_translation.o glui_tree.o glui_treepanel.o glui_window.o quaternion.o viewmodel.o
GLUI_OBJSWIN = $(GLUI_OBJS:.o=.obj)

libglui.a: $(GLUI_OBJS)
	@echo "Creating library"
	ar rcs libglui.a $(GLUI_OBJS)

libglui.lib: $(GLUI_OBJSWIN)
	lib $(GLUI_OBJSWIN) -OUT:libglui.lib 

