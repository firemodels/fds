## Intel Inspector

This folder contains a few configuration files that work with the [Intel Inspector](https://software.intel.com/en-us/node/622387), a tool that can help detect improperly coded OpenMP directives.

### Recommended Setup

For useful analysis capabilities, using inspxe-gui on the platform used for collection is ideal. By doing so, Inspector can have direct access to relevant source files, and thus give more useful results that do not require investigation with a secondary code viewer. 

To setup, use inspxe-gui before collecting data. X11 forwarding is necessary if logging in to a remote cluster. Inside the GUI, use the new project option to create a project for FDS, selecting the `inspect` version's executable as the target, the 'inspect' version's folder for binary files, and firemodels/fds/Source for source files. To use the project, place collection results, obtained later, into the generated project folder.

### Collection

1. Compile the `inspect` version of FDS using the script `make_fds.sh` in this directory. The relevant compiler options are listed [here](https://software.intel.com/en-us/inspector-user-guide-linux-building-applications).
2. Run fds with 'qfds' using -x [result directory] - this runs in Intel's ti2 form with a stack depth of 32. This specifically targets thread errors, and requires at least 2 OpenMP threads to be in use to collect results. You can read more [here](https://software.intel.com/en-us/inspector-user-guide-linux-threading-error-analysis-types)
3. From firemodels/fds/Verification/Thread_Check, you can get inspect_report.sh. This script can handle the multiple files generated in any Inspector run with more than one process. 

What these commands do is run the command line version of Inspector (inspxe-cl) on FDS, which returns an Inspector file that needs further analysis by the developer.

### Analysis

From firemodels/fds/Verification/Thread_Check, you can get inspect_report.sh. This script can handle the multiple files generated in any Inspector run with more than one process.

Alternatively, you can use the recommended inspxe-gui. With results placed/generated in the project folder, Inspector can display the results in a convenient GUI, instead of necessitating parsing the output on the command line.

