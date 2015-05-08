Samoa
=====

Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications.
Website: [Samoa](https://github.com/meistero/Samoa)

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
* (Optional) [ASAGI](https://github.com/tum-i5/ASAGI) and netcdf data files for production runs (not publicly available, contact the developers)

## Installation

### Local systems

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
scons must be installed manually on the SuperMUC as the current version is not compatible with sam(oa)².
Assuming scons has been installed in the folder scons_dir, set the following variables:

    module load python
    export PATH=<scons_dir>/build/scripts/:$PATH
    export PYTHONPATH=<scons_dir>/build/lib/:$PYTHONPATH
    export SCONS_LIB_DIR=<scons_dir>/build/lib/

Additionally, in order to compile and run ASAGI and sam(oa)² on the SuperMUC the following modules must be loaded:

    module load git gcc/4.7 cmake netcdf
    module switch ccomp ccomp/intel/13.1
    module switch fortran fortran/intel/13.1
    module switch mpi.ibm mpi.intel

sam(oa)² supports both multithreaded and single-threaded MPI. Both ASAGI and sam(oa)² must link to the same respective libraries, thus it is necessary to compile ASAGI twice:
once without MT support and once with MT support. Rename the single-threaded library to "libasagi_nomt.so" and the multi-threaded library to "libasagi.so".

At this point, you should be able to compile ASAGI and sam(oa)².

### Linux Cluster and MAC Cluster

Create a directory (named samoa_dir here) and execute the following steps:

    cd <samoa_dir>
    git clone https://github.com/meistero/samoa .

This will download the source files for samoa into samoa_dir.
The following modules should be loaded before compiling ASAGI and sam(oa)² on the Linux and MAC clusters

    module unload gcc python
    module load git cmake netcdf gcc/4.7 python/2.7.5
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

    scons asagi=noasagi ...

Executables will be created in the directory samoa_dir/bin and should be run from samoa_dir.

## Execution

For execution parameters refer to the online help by calling the executable with '-h' or '--help'.

##Build status

![Travis build status not found](https://travis-ci.org/meistero/Samoa.svg)

