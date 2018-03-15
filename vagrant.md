# Compiling FDS with GNU Fortran, Open MPI and Intel Performance Libraries in Ubuntu Linux

The basic idea is to create a build environment to compile FDS source code in an automated way - reproducible, version controlled and fast.
Vagrant can provide this to every developer. A file `Vagrantfile` to create and provision a virtual machine is now added to the source code repository in top-level folder.

With regard to running FDS it is a good idea to isolate an end-user environment used for production jobs from an environment capable to compile FDS Fortran source code. End-users just install the official bundle and usually only care about system administration when they have issues. When support end-users reconfigure their environmentto compile nstall software dependencies required


## Why vagrant

If you are a developer, Vagrant will isolate dependencies and their configuration within a single disposable, consistent environment, without sacrificing any of the tools you are used to working with (editors, browsers, debuggers, etc.). Once you or someone else creates a single `Vagrantfile`, you just need to vagrant up and everything is installed and configured for you to work. Other members of your team create their development environments from the same configuration, so whether you are working on Linux, Mac OS X, or Windows, all your team members are running code in the same environment, against the same dependencies, all configured the same way. Say goodbye to "works on my machine" bugs.
The Vagrantfile is meant to be committed to version control with your project, if you use version control. This way, every person working with that project can benefit from Vagrant without any upfront work.

[Vagrant might be for everyone](https://www.vagrantup.com/intro/index.html#for-everyone)

## Install Vagrant 

Download Vagrant
- [Downloads for latest version of Vagrant](https://www.vagrantup.com/downloads.html)
    - Please download the proper binary package for your operating system and architecture. 

Verify installation

```bash
$ vagrant --version
Vagrant 1.9.8
```

>The installer will automatically add vagrant to your system path so that it is available in terminals.

## Install a provider

Vagrant comes with support out of the box for VirtualBox, a free, cross-platform consumer virtualization product.

VirtualBox must be installed on its own prior to using the provider, or the provider will display an error message asking you to install it. VirtualBox can be installed by downloading a package or installer for your operating system and using standard procedures to install that package.

- [Download a package or installer for your operating system](https://www.virtualbox.org/wiki/Downloads)

## Spin up Vagrant box 

The build environment uses Ubuntu 14.04.5 LTS (Trusty Tar) to demonstrate compiling FDS using free and open source software (FOSS). Boot your Vagrant environment. 

    $ vagrant up

## Provisioning

Provisioners in Vagrant allow you to automatically install software, alter configurations, and more on the machine as part of the `vagrant up` process.

Create a build environment
- Install GNU Fortran 7.2
- Build Open MPI 3.0.0
- Install Intel Math Kernel Library
- Compile FDS sources
  - build target `mpi_gnu_linux_64`

## Connect to Vagrant box

 Vagrant runs the virtual machine without a UI. To prove that it is running, you can SSH into the machine

    $ vagrant ssh

## Run FDS executable on Vagrant box

By default, Vagrant shares your project directory (remember, that is the one with the Vagrantfile) to the `/vagrant` directory in your guest machine.

Note that when you `vagrant ssh` into your machine, you're in `/home/vagrant`. `/home/vagrant` is a different directory from the synced `/vagrant` directory.
    
    $ pwd
    /home/vagrant

Run FDS in your home directory

    $ /vagrant/Build/mpi_gnu_linux_64/fds_mpi_gnu_linux_64 case.fds 

## Support for `GLMAT` pressure solver

To use `GLMAT` pressure solver in your models FDS must be compiled with Intel MKL library. Support for automated setup of Intel MKL in this build environment was added just recently.

Actually there is an issue with the automated way to enable MKL environment. Please use the steps below to set environment for Intel MKL manually.

Connect to Vagrant box

    $ vagrant ssh

Inside Vagrant box

    $ source /opt/intel/mkl/bin/mklvars.sh intel64
    $ source /vagrant/open-mpi-vars.sh
    $ cd /vagrant/Build/mpi_gnu_linux_64
    $ rm -rf *.mod *.o
    $ ./make_fds.sh


>It is possible that you need to delete existing `.mod` and `.o` files in build directory compiled without linking to MKL library.


 









