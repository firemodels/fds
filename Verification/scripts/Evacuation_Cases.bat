%RUNFDS% Evacuation evac_memory_test1             fire41
copy Evacuation\evac_memory_test1_evac.fed Evacuation\evac_memory_test2_evac.fed
copy Evacuation\evac_memory_test1_evac.eff Evacuation\evac_memory_test2_evac.eff
%RUNFDS% Evacuation evac_memory_test2          fire42
%RUNFDS% Evacuation evac_smv_testcase1             fire41
copy Evacuation\evac_smv_testcase1_evac.fed Evacuation\evac_smv_testcase2_evac.fed
copy Evacuation\evac_smv_testcase1_evac.eff Evacuation\evac_smv_testcase2_evac.eff
%RUNFDS% Evacuation evac_smv_testcase2          fire42
