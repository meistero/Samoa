[![Logo](https://raw.githubusercontent.com/meistero/Samoa/master/logo_small.png)]( https://raw.githubusercontent.com/meistero/Samoa/master/logo.png)

sam(oa)² 
=======

Space-Filling Curves and Adaptive Meshes for Oceanic And Other Applications. <br>
Github repository: [https://github.com/meistero/Samoa](https://github.com/meistero/Samoa)

## Contents

1. Prerequisites
2. Installation
 1. Local systems
 2. SuperMUC
 3. Linux Cluster and MAC Cluster
3. Compilation
4. Execution

## Prerequisites

The following prerequisites are necessary in order to install and run sam(oa)²:
* [git](http://git-scm.com/)
* [scons](http://www.scons.org/)
* gfortran 4.7 or higher OR Intel Fortran Compiler 13.0 or higher
* (Optional) [ASAGI](https://github.com/tum-i5/ASAGI) v0.5.0 or higher for external geodata
* (Optional) Netcdf data files for ASAGI: For porous media flow, download the SPE10 data files from [SPE10](http://www.spe.org/web/csp/datasets/set02.htm#download). A script is included in the data directory that converts them to netcdf files. For the tsunami scenario the netcdf files can be generated from our [Tsunami repository](https://github.com/TUM-I5/tsunami)

## Installation

### Local Systems

Create a directory (named samoa_dir here) and execute the following steps:

    cd <samoa_dir>
    git clone https://github.com/meistero/samoa .

This will download the source files for samoa into samoa_dir. 

### SuperMUC

Create a directory on the SuperMUC (named samoa_dir here). SuperMUC restricts access to outside sources and thus does not allow connections to https servers. However, there are
two methods to clone sam(oa)² from github on the SuperMUC:
* By accessing the SuperMUC file system as a remote directory. git can then be executed locally:

        nohup sshfs <login>@supermuc.lrz.de:<samoa_dir> <local_dir>
        cd <local_dir>
        git clone https://github.com/meistero/samoa .

* By login with remote port forwarding. In this case an alternative URL must be used to clone the git repository:


        ssh -X <login>@supermuc.lrz.de -R <port>:github.com:9418
        cd <samoa_dir>
        git clone git://localhost:<port>/meistero/Samoa .

This will download the source files for samoa into samoa_dir. 
Additionally, in order to compile and run ASAGI and sam(oa)² on the SuperMUC, we must add the netcdf library to the CMAKE prefix path and load the following modules:

    module unload gcc
    module load git scons gcc/4.7 cmake/4.1 netcdf
    export CMAKE_PREFIX_PATH=$NETCDF_BASE

At this point, you should be able to compile ASAGI and sam(oa)².

### Linux Cluster and MAC Cluster

Create a directory (named samoa_dir here) and execute the following steps:

    cd <samoa_dir>
    git clone https://github.com/meistero/samoa .

This will download the source files for samoa into samoa_dir.
The following modules should be loaded before compiling ASAGI and sam(oa)² on the Linux and MAC clusters

    module unload gcc python
    module load git cmake scons netcdf gcc/4.7
    module load gnuplot

sam(oa)² supports both multithreaded and single-threaded MPI. Both ASAGI and sam(oa)² must link to the same respective libraries, thus it is necessary to compile ASAGI twice:
once without MT support and once with MT support. Rename the single-threaded library to "libasagi_nomt.so" and the multi-threaded library to "libasagi.so".

At this point, you should be able to compile ASAGI and sam(oa)².

## Compilation

In order to view all the compilation options sam(oa)² provides, you can execute the following command now:

    scons --help

Typical settings are:

    scons asagi_dir=<asagi_dir> compiler=gnu scenario=darcy -j<threads>
    scons asagi_dir=<asagi_dir> compiler=intel target=debug scenario=swe -j<threads>

If you wish to simulate simple scenarios that do not require data files you can also disable asagi with the flag

    scons asagi=No ...

Executables will be created in the directory samoa_dir/bin and should be run from samoa_dir.

## Execution

For execution parameters refer to the online help by calling the executable with '-h' or '--help'.

## Build Status

[![Build Status](https://travis-ci.org/meistero/Samoa.svg)](https://travis-ci.org/meistero/Samoa)

