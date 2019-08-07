## Intel Advisor

This folder contains the optimal build for Intel Advisor, which can assist in determining effective optimization locations for threading and vectorization. 

### Required Setup

1. Source advixe-vars.sh, from the installation folder for Inspector. For example:
```
source /opt/intel19/advisor/advixe-vars.sh
```
2. Compile the `advise` version of FDS using the script `make_fds.sh` in this directory. The relevant compiler options are listed [here](https://software.intel.com/en-us/get-started-with-advisor-before-you-begin).

### Recommended Steps

For useful analysis capabilities, using advixe-gui on the platform used for collection is ideal, or at least having access to the used source and executable on the machine where Advisor is intalled.

To setup, use `advixe-gui` before collecting data. X11 forwarding is necessary if logging in to a remote cluster. Inside the GUI, use the new project option to create a project for FDS, selecting the 'advise' version's executable as the target. 

### Collection

The base command used on one's platform to run FDS is a proper starting point. `mpiexec fds [test case]` is a common input, and thus we'll use it for the model here.

Before the fds executable in the run command, such as between mpiexec and the executable, place `advixe-cl -collect <analysis-type> [Optional actions] -- `. Most frequently, analysis-type can be 'survey'. This obtains a basic look at the program, and can be followed up with 'suitability' after annotation. Note that the others have not been incredibly useful/stable in FDS development. You can learn about annotating [here](https://software.intel.com/en-us/advisor-user-guide-annotating-code-for-deeper-analysis)

Thus, a possible input could be:
```
mpiexec -np 1 advixe-cl -collect survey -- $HOME/firemodels/fds/Build/impi_intel_linux_64_advise/fds_impi_intel_linux_64_advise simple_test.fds
```


### Analysis

#### Graphical User Interface (Recommended)

Here, the command line interface has not been explored, so advixe-gui is the offering that can be recommended. Open the project file with the GUI to analyze results.

### IMPORTANT NOTE

Generally, unless directly interested in testing the suitability of several locations at once, or in vectorization, Amplifier is a more useful tool for FDS development. It can be found in the 'vtune' build folder.