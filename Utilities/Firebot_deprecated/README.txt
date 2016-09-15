=====================================
= Firebot configuration information =
=====================================

=========
= About =
=========

Firebot is an automatic verification and validation test bot that is run at a regular interval (nightly).
More details on the Firebot build stages can be found in the FDS Configuration Management Plan.

============================
= Running Firebot Manually =
============================



======================
= Installing Firebot =
======================

1. clone FDS-SMV repository

2. Ensure that the following software packages are installed on the system:

    Intel compilers and Intel Inspector
    LaTeX (TeX Live distribution), be sure to make this the default LaTeX in the system-wide PATH
    Matlab (test the command 'matlab')

3. Firebot uses email notifications for build status updates. Ensure that outbound emails can be sent using the 'mail' command.

4. Install libraries for Smokeview. On CentOS, you can use the following command:

    yum install mesa-libGL-devel mesa-libGLU-devel libXmu-devel libXi-devel xorg-x11-server-Xvfb

5. Add the following lines to firebot's ~/.bashrc file:

    . /usr/local/Modules/3.2.10/init/bash
    module load null modules torque-maui mpi/openmpi-1.8.1-gnu-ib

    export IFORT_COMPILER=/opt/intel/composerxe
    export IFORT_COMPILER_LIB=/opt/intel/composerxe/lib

    #FDS
    source ~/.bashrc_fds intel64

    # Set unlimited stack size
    ulimit -s unlimited

6. Setup passwordless SSH on the firebot account. Generate SSH keys and ensure that the head node can SSH into all of the compute nodes. Also, make sure that firebot's account information is propagated across all compute nodes (e.g., with the passsync or authcopy command).

7. Ensure that a queue named 'firebot' is created, enabled, and started in the torque queueing system and that nodes are defined for this queue. Test the 'qstat' command on firebot's account.


=========================
= Firebot files/scripts =
=========================

# firebot_linux_wrapper.sh or firebot_mac_wrapper.sh

    The _wrapper script uses a semaphore file that ensures multiple instances of Firebot do not run, which would cause file conflicts.
    This script should be called from crontab to start firebot.

# firebot_linux.sh

    This is the primary Firebot automated test script used on Linux with a queueing system (e.g., TORQUE).
    This script is invoked via crontab (details below).
    By default, Firebot is run in verification mode, but the -v option can be used to run Firebot in validation mode.
    Note: In validation mode, the file /FDS-SMV/Validation/Process_All_Output.sh is used to obtain a list of
    active (uncommended) validation cases.

===========
= Crontab =
===========

---------------------------------------------------------------------------
#### The following information is in the Linux firebot user's crontab: ####
---------------------------------------------------------------------------

PATH=/bin:/usr/bin:/usr/local/bin:/home2/smokevis2/firebot/firebot:$PATH
MAILTO=""

#  ========================
#  = Firebot build script =
#  ========================

# Update and run Firebot at 9:56 PM every night
# If no SVN argument is specified, then the latest SVN revision is used
56 21 * * * cd ~/firebot ; svn revert * ; svn up ; bash -lc "./firebot_linux_wrapper.sh"
