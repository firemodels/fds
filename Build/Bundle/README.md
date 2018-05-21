# Bundle FDS

The scripts make_bundle.sh and make_bundle.bat are used for creating FDS/Smokeview installation bundles.  make_bundle.sh is used for  Linux and OSX platforms and make_bundle.bat is used for Windows.
Environment variables required by make_bundle.sh are defined in ~/.bundle/FDS_SMV_ENV.sh.  

To use this script, type ./make_bundle.sh in a bash shell (in this directory). To see usage information type ./make_bundle.sh -h . 
Likewise to use make_bundle.bat, type make_bundle.bat in a windows command shell (again in this directory).

These scripts build all the software (by calling other scripts), collect manuals, example files and libraries that are needed then builds the bundle script.  The created bundles are placed in the directory Utilities/uploads. 

Note, this is as work in progress, these scripts are not yet complete.  
