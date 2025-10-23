# FDS Input File Tools

This directory contains a few handy scripts to help in writing input files.

## Mass Producing Input Files from a Template

The Python script called `swaps.py` reads a comma-delimited spreadsheet file called `paramfile.csv` with the following format:

| Template   | param1 | param2 | param3 |
|------------|--------|--------|--------|
| file01.fds | 12.3   | text1  | 6      |
| file02.fds | 17.4   | text2  | 8      |
| file03.fds | 19.3   | text3  | 98     |

The script then reads the file called `Template` which is an FDS input file except that it contains the character strings `param1`, `param2`, and `param3`. It then creates the three FDS input files `file01.fds`, `file02.fds`, and `file03.fds` in which the character strings in the first row are replaced by the character strings in the subsequent rows. There is no number or string typing; the strings in the first row are replaced exactly as is by the strings in the subsequent rows. There is no need to add quotes or anything else.

## Creating a simple pan using OBST lines

The script `pan_with_lip.py` creates a simple circular pan with a lip. Input your desired dimensions directly in the script.
