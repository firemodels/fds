Latest changes: T.K., March 2, 2010, SVN $Revision: 4582 $

==============================
Evacuation test cases for FDS5
==============================

===========================
TEST CASES FOR MEMORY USAGE
===========================

These test cases consume quite much RAM.  If one needs more memory
usage to test something, then the MESH IJK could be doubled (and so
on).

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

   fds5.4.0, svn3895:
   fire + evac:    616 472 K RAM (Peak Memory usage, winXP32, test1)
   no fire + evac: 556 404 K RAM (Peak Memory usage, winXP32, test2)
   fire + no evac:  63 596 K RAM (Peak Memory usage, winXP32, test0)

   fds5.4.2, svn5104:
   fire + evac:    224 324 K RAM (Peak Memory usage, winXP32, test1)
   no fire + evac: 162 808 K RAM (Peak Memory usage, winXP32, test2)
   fire + no evac:  65 952 K RAM (Peak Memory usage, winXP32, test0)

========================
TEST CASES FOR SMOKEVIEW
========================

Next cases have two floors, which are connected, i.e., the agents are
going from the second floor to the first floor.  These test cases
use all the main features of the evacuation module, e.g., they move
agents from one mesh to some other mesh.

evac_smv_testcase0.fds

   Just fire mesh

evac_smv_testcase1.fds

   Fire + evacuation meshes

evac_smv_testcase2.fds

   Just evacuation meshes

