The Linux version of Smokeview links to shared libraries. As a result,
the environmental variable, LD_LIBRARY_PATH, needs to be created or modified
(if it already exists). The following command (if you are using the csh or
tcsh shell) points to 32 bit shared libraries assuming you install FDS in your
home directory.  Change 32 to 64 if you are using the 64 bit version of
Smokeview.


setenv LD_LIBRARY_PATH ~/FDS/FDS5/bin/INTEL_LIB32
