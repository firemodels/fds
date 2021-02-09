Latest changes: T.K., January 2, 2021
Git Revision FDS6.7.5-857-g94af16e64

=======================================================
Evacuation test cases for FDS 6.7.5, Evac version 2.6.0
=======================================================

Evac 2.6.0: A new fire+evacuation strategy, where fire and evacuation
            meshes are separated: fire+evacuation has 3 phases

Phase 1) &MISC EVACUATION_INITIALIZATION=T

         - Does an evacuation drill calculation according to the
           inputs

         - Writes CHID_evac.eff (guiding flow fields for evacuation
           movement) to the hard drive (could be used in the "MC" mode
           of the evacuation calculation

         - Writes CHID_evac.xyz file that is used, if the (separate)
           fire calculation is wanted to save the fire information for
           evacuation calculation (Phase 3)
           

Phase 2) &MISC EVACUATION_WRITE_FED=T

         - Does an ordinary fire calculation according to the inputs

         - Reads CHID_evac.xyz file (Phase 1) and writes CHID_evac.fed
           file that is used in Phase 3, where the evacuation
           calculation uses fire data

Phase 3) &MISC EVACUATION_MC_MODE=T

         - Does an evacuation calculation that uses fire information

         - Reads the CHID_evac.eff and CHID_evac.fed files (might
           recalculate CHID_evac.eff, if this is not found: This is
           still under construction, i.e., this is not yet tested for
           the Evac 2.6.0 strategy)

         - If &MISC EVACUATION=DRILL=T, then does not (try to) read
           the CHID_evac.fed fire data file, does an evacuation drill
           simulation in the "MC mode", i.e., reads in the
           CHID_evac.eff file, which will make the calculation faster.


===========================
TEST CASES FOR MEMORY USAGE
===========================

These test cases consume quite much RAM.  If one needs more memory
usage to test something, then the MESH IJK could be doubled (and so
on).

evac_memory_test0.fds

   One fire mesh, no evacuation meshes

   VERIFICATION OK IF: "STOP: FDS completed successfully"

evac_memory_test1a.fds,evac_memory_test1b.fds,evac_memory_test1c.fds

   One fire mesh, one main evacuation mesh
   a: Phase 1, writes the .eff and .xyz files
      copy evac_memory_test1a_evac.xyz => evac_memory_test1b_evac.xyz 
      copy evac_memory_test1a_evac.eff => evac_memory_test1c_evac.eff 
   b: Phase 2, reads the .xyz file and write the .fed file 
      copy evac_memory_test1b_evac.fed => evac_memory_test1c_evac.fed 
   c: Phase 3, reads the .eff and .fed files   

   VERIFICATION OK IF: a,b,c, "STOP: FDS completed successfully"

========================
TEST CASES FOR SMOKEVIEW
========================

Next cases have two floors, which are connected, i.e., the agents are
going from the second floor to the first floor.  These test cases
use all the main features of the evacuation module, e.g., they move
agents from one mesh to some other mesh.

evac_smv_testcase0.fds

   Just fire calculation (NO_EVACUATION=.TRUE.)

   VERIFICATION OK IF: "STOP: FDS completed successfully"

evac_smv_testcase1a.fds,evac_smv_testcase1b.fds,evac_smv_testcase1c.fds

   Fire + evacuation meshes
   a: Phase 1, writes the .eff and .xyz files
      copy evac_smv_testcase1a_evac.xyz => evac_smv_test1case1b_evac.xyz 
      copy evac_smv_testcase1a_evac.eff => evac_smv_test1case1c_evac.eff 
   b: Phase 2, reads the .xyz file and write the .fed file 
      copy evac_smv_testcase1b_evac.fed => evac_smv_test1case1c_evac.fed 
   c: Phase 3, reads the .eff and .fed files   

   VERIFICATION OK IF: a,b,c, "STOP: FDS completed successfully"

