Samoa
=====

Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications.
Website: [Samoa](https://github.com/meistero/Samoa)

Contents:

1. Prerequisites
2. Installation
2a. SuperMUC
2b. Linux Cluster
3. Running

1. Prerequisites:
-----------------

The following prerequisites are necessary in order to install and run sam(oa)²:
* [git](http://git-scm.com/)
* [scons](http://www.scons.org/)
* [ASAGI](https://github.com/tum-i5/ASAGI)
* gfortran >= 4.7 OR Intel Fortran Compiler >= 13.0

2. Installation:
-----------------

After installing all the prerequisites, create a directory <samoa_dir> and execute the following steps:

    cd <samoa_dir>
    git clone https://github.com/meistero/samoa .

This will download the source files for samoa into the current directory. You can execute

    scons --help

now, in order to view all the compilation options sam(oa)² provides. Typical settings are:

    scons asagi_dir=<asagi_dir> compiler=gnu scenario=darcy -j<threads>
    scons asagi_dir=<asagi_dir> compiler=intel target=debug scenario=swe -j<threads>

Executables will be created in the directory <samoa_dir>/bin and should be run from <samoa_dir>

2a. SuperMUC:
-------------


