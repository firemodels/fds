#!/bin/bash

coupler=../Utilities/coupling_emulator/intel_linux_64/coupling_emulator_linux_64
fds=../Build/intel_linux_64/fds_linux_64
$coupler $fds w
