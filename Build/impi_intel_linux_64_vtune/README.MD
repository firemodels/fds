#Intel VTune Amplifier

This folder contains the optimal build for Intel VTune Amplifier, which can locate hotspots and the performance of both parallel and serial code. This helps for targeting optimizations.

### Recommended Setup

For useful analysis capabilities, using amplxe-gui on the platform used for collection is ideal, or at least having access to the used source and executable on the machine where Amplifier is installed.

To setup, use amplxe-gui before collecting data. X11 forwarding is necessary if logging in to a remote cluster. Inside the GUI, use the new project option to create a project for FDS, selecting the 'vtune' version's executable as the target. When prompted in analysis, select the 'vtune' version's folder for binary files, and firemodels/fds/Source for source files. To use the project, place collection results, obtained later, into the generated project folder.

### Collection

1. Compile the `vtune` version of FDS using the script `make_fds.sh` in this directory. The relevant compiler options are listed [here](https://software.intel.com/en-us/vtune-amplifier-help-compiler-switches-for-performance-analysis-on-linux-targets).
2. Run fds with 'qfds' using -a [result directory] - this analyzes the code in terms of high performance computing, focusing on serial vs parallel code. You can read more [here](https://software.intel.com/en-us/vtune-amplifier-help-collect)

What these commands do is run the command line version of Amplifier (amplxe-cl) on FDS, which returns an Amplifier folder that needs further analysis by the developer.

### Analysis

Due to the large amount of data collected in any Amplifier run, using the amplxe-gui is highly recommended With results placed/generated in the project folder, Amplifier can display the results in a convenient GUI, instead of necessitating parsing the output on the command line.

In the GUI, there are several display options under 'Bottom-up' in a results screen. When optimizing existing OpenMP directives, focus on using options that highlight OpenMP regions. Otherwise, focus on function or source file views. Hovering over data labels normally gives a cursory explanation of what is represented.