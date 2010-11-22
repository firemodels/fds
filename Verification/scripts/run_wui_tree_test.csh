#!/bin/csh -f

$RUNFDS Wui onetree_surf_1mesh.fds fire47
$SMOKEZIP -part2iso onetree_surf_1mesh
