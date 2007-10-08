// $Date$ $Author$
#define DLL
#define svWINAPI

DLL void svWINAPI sv_init0(void);
DLL void svWINAPI sv_startup(char *file, int showpart);

DLL void svWINAPI sv_update(void); /* glut main event loop; never returns !! */
DLL void svWINAPI sv_unload(void);  /* unload menus and memory associated with loaded data files. 
                        then hide the window. */
/* menu calls */

DLL void svWINAPI sv_grid(int value); /* 0=hide grid, 1=show grid */
DLL void svWINAPI sv_BlockageMenu(int value);  /* 0=hide 1=showall 2=show outline */
DLL void svWINAPI sv_MainMenu(int value);      /* 1=tour 2=reset view 3=quit */
DLL void svWINAPI sv_FrameRateMenu(int value); /* value > 1000 => unlimited,
                                         value < 0 => step 
                                         else maxframe rate=value */


DLL void svWINAPI sv_startup_c(int argc, char **argv);
