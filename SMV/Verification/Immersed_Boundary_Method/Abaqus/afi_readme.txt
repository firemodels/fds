To run the Abaqus demo examples, please run the FDS job using the updated FDS first, and then type the following in the command line as you run a normal abaqus job:

abaqus -j job.inp -user afi.obj int

where job.inp is the abaqus input file. At this moment, the abaqus job name is assumed to be the same as the FDS job name, i.e. job=CHID; afi.obj is a library file compiled in 32-bit Windows OS. Libraries for other operation systems will be provided later on.

If you have any question regarding the FDS-Coupled simulations, please email to Charles Luo cluo@gem-innovation.com.