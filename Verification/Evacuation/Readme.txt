File started by T.K., May 4, 2009

==============================
Evacuation test cases for FDS5
==============================

===========================
TEST CASES FOR MEMORY USAGE
===========================

These test cases consume quite much RAM.  If one needs more memory
usage to test something, then the MESH IJK could be doubled (and so
on).  Note that the memory (fds5.4.0, svn3895) is mainly used for the
evacuation meshes.

   fds5.4.1, svn 4746:
   fire + no evac:  65 340 K RAM (Peak Memory usage, winXP32, testcase0)
   fire + evac:    217 928 K RAM (Peak Memory usage, winXP32, testcase1)
   no fire + evac: 156 972 K RAM (Peak Memory usage, winXP32, testcase2)

evac_memory_test0.fds

   One fire mesh, no evacuation meshes

evac_memory_test1.fds

   One fire mesh, one main evacuation mesh, four additional evacuation
   meshes (for each door)

evac_memory_test2.fds

   This consumes quite much RAM. If one needs more memory usage to
   test something, then the MESH IJK could be doubled (and so on).

   No fire mesh, one main evacuation mesh, four additional evacuation
   meshes (for each door)

========================
TEST CASES FOR SMOKEVIEW
========================

Next cases have two floors, which are connected, i.e., the agents are
going from the second floor to the first floor.  These test cases
use all the main features of the evacuation module, e.g., they move
agents from one mesh to some other mesh.

evac_smv_testcase1.fds

   Fire + evacuation meshes

evac_smv_testcase2.fds

   Just evacuation meshes

