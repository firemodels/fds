## Code Documentation using Doxygen

This directory contains a configuration file, `Doxyfile`, for use with the code documentation program [Doxygen](http://www.doxygen.org). Assuming that Doxygen is installed on your system, type:
```
doxygen Doxyfile
```
and a new directory called `html` will be created containing an HTML version of documentation for the FDS code.  The main page is `html/index.html`.

If you want to create your own configuration file, type
```
doxygen -g myDoxyfile
```
