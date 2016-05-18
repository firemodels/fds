Creating an FDS Terrain Case (preliminary)

Preliminaries
  Install or verify following is installed
  1.  Google Earth
  2.  latest test smokeview
  3.  screen capture program such as snagit (just a suggestion, not a recommendation)

Create FDS Input File
  1.  Open Google Earth
  2.  Identify a rectangular region of interest, noting the lat/lon coordinates of two opposite corners
  3.  Capture this region using your screen capture program.  Save the image as casename.jpg where casename is the 'CHID' of your FDS input file
  4.  creat a file named casename.in containing:
    longitude_begin  longitude_end n_longitudes  latitude_begin latitude_end n_latitudes

     a. note google earth gives latitudes and longitudes as d=degrees m=minutes s=seconds
        dem2fds, the program that uses casename.in, assumes that latitudes and longitudes are in decimal ie d + m/60 + s/3600
     b. also note longitudes for USA are negative
    5. type the command:
     dem2fds casename < casename.in

   this command generated a series of files named casename_longlats_xxx.csv .  each file contains 400 longitude/latitude pairs

6. Go to http://viewer.nationalmap.gov/theme/elevation/ and click on bulk

   a.  click on Choose Files and select a file for each casename_longlats_xxx.csv generated in step 5
   b.  click on Get Elevations
   c.  after that web site retrieved all the data requested save the file

7. Concatenate all elevation files generated in step 6 removing header line.  Call this file casename_elevs.in

8. Create an FDS input file using
   dem2fds -o casename


