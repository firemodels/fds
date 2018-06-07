# Intel Trace Analyzer Configuration

The files in this folder are example .itarc files. Segments of these files should be pasted into a configuration file to be loaded into Intel Trace Analyzer as described below. In order to make these settings the default for Intel Trace Analyzer, paste into the hidden .itarc file in your home directory (the directory with your username for Linux and C:\Users\%username%\ for Windows).

**function_colors.itarc:**

Copy all of the "File" tag with name "fds_impi_intel_linux_64_db.stf". Paste into another file beneath the first closing "File" tag (```</File>```) and replace the name and path to match the file you would like to have colored functions for. Intel Trace Analyzer only allows function coloring settings on a per file basis, so there is no way to set these settings for all files to be opened. Note that to see function colors, it is necessary to right click on the Event Timeline and click "Ungroup Group Application".
