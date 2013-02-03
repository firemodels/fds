.SUFFIXES: .c .cpp

SIZE=-m64

PLATFORM=

GLUT_INC_LOCATION=../include2

CFLAGS = 

RM=/bin/rm
COMPILER=icc
CC=${COMPILER}
CPP=${COMPILER2}
CPPFLAGS=$(PLATFORM) -fpermissive -I${GLUT_INC_LOCATION} -I ../include2/glut ${SIZE}
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

GLUI_OBJS = glui_add_controls.o glui.o glui_bitmap_img_data.o glui_bitmaps.o glui_button.o glui_edittext.o glui_checkbox.o glui_node.o glui_radio.o glui_statictext.o glui_panel.o glui_separator.o glui_spinner.o glui_control.o glui_column.o glui_translation.o glui_rotation.o glui_mouse_iaction.o glui_listbox.o glui_rollout.o arcball.o algebra3.o quaternion.o

libglui.a: $(GLUI_OBJS)
	@echo "Creating library"
	ar rcs libglui.a $(GLUI_OBJS)

#.cpp.o:
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCS) $*.cpp

#.c.o:
#	$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCS) $*.c
clean:
	$(RM) *.o
