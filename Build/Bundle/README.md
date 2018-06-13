# Bundle FDS

The scripts make_bundle.sh and make_bundle.bat are used for creating FDS/Smokeview installation bundles.  The bash script make_bundle.sh is used to create bundles for Linux and OSX systmes and the batch file make_bundle.bat is used to create a bundle for Windows systems.
Environment variables required by make_bundle.sh are defined in ~/.bundle/FDS_SMV_ENV.sh.  

To use this script, open a bash shell, cd to Build/Bundle and type `./make_bundle.sh`.
To see usage information type: `./make_bundle.sh -h` . 
Likewise to use make_bundle.bat, type `make_bundle.bat` in a windows command shell.

These scripts build all the software (by calling other scripts), collect manuals, example files and libraries that are needed then builds the bundle script.  The created bundles are placed in the directory Build/Bundle/uploads. 

Note, this is as work in progress, these scripts are not yet complete.  
