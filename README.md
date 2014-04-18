###FALCON###
====

Flux Assignment (with) LAD Convex Objectives (and) Normalization.
LAD, or Least Absolute Deviations, is the 1-norm equivalent of least squares.

Currently, two algorithms are implemented. The first algorithm estimates
enzyme complex abundance, which is fed in to the second algorithm that
also takes a genome-scale metabolic model as input. The output of the
second algorithm are estimated fluxes for the environment that the
expression and model were taken from.


### Installation ###
====

FALCON should work on all major platforms (Windows, Linux, OS X). 
At this time it is required that the user is familiar with the UNIX
command-line and has some experience installing UNIX software. It
is possible that it could work on other platforms as well, but some of
the MATLAB code would need to be ported to GNU Octave. In principle
should be straightforward.

#### Required Prerequisites ####

There are several essential prerequisites:

* **MATLAB**: currently required for running the flux-fitting algorithm and 
 running most of the analyses. This is the only non-free requirement. 

* **COBRA Toolbox**: A MATLAB package developed for analysis of genome-scale
 metabolic models. We assume the COBRA standards and use many of its functions
 in our code.  See the [official page](http://opencobra.sourceforge.net/) 
 for documenation and the protocol. We have used a fork of the project available 
 at [Ronan Fleming's repository](https://github.com/rmtfleming/cobratoolbox).


* **Boehm GC:** This should be available
 on most systems, though in Cygwin it will have to be explicitly installed
 using the setup program. The headers may need to be installed separately, which
 is necessary if one wishes to compile the minimum disjunction program.

* **libgcc** This is the GCC runtime library, and should be installed on most
 systems by default. Sometimes the appropriate version may need to be installed
 if you are using an executable of minDisj compiled by someone else 
 (if this is the case it is best you compile your own copy of minDisj).

#### Optional Packages ####
* **C compiler**: We have only tested recent versions of GCC. This is
 needed for building the implementation of the first algorithm in the
 pipeline (minimum-disjunction). For Windows, we suggest using GCC in
 [Cygwin](http://cygwin.com/), as this is also currently needed for an
 optional but desirable prerequisite: ATS.

* **ATS:** The first algorithm is written in ATS2, and while it should
 be possible to compile a version of it from the C sources we have
 generated, if one runs in to problems the easiest thing to do will be
 to try to compile from ATS. Nonetheless, please report any bugs
 encountered in compiling from C sources in the GitHub issue tracker.
 ATS2 can be installed by following directions on the [ATS
 site](http://www.ats-lang.org/DOWNLOAD/) and [the ATS2
 wiki](http://sourceforge.net/p/ats2-lang/wiki/Building%20and%20installing/).
 Users interested in altering or extending minDisj or GPR-related work
 file will need to install the ATS language, which has a style similar
 to StandardML.


#### Generic installation instructions ####
After having acquired MATLAB and having installed the COBRA Toolbox,
all that is needed is to install the minimum-disjunction program. It
may be worth trying out a [provided
binary](https://app.box.com/s/ujwvl9vsbzsjyq7qpjao), which are
available for some platforms. These are not guaranteed to be up to
date or to work on all distributions or versiosn of systems. If it
does not work, please try building it instead, using either the
provided C sources or the ATS sources.  Should building fail, please
submit a report to the GitHub issue tracker.

Copy or create a link to the file with the name minDisj, and
make sure it is in your PATH.

Once minDisj is built, consider adding the following to your startup.m
file (in UNIX, this should be in ~/Documents/MATLAB/), or run it
before each use of FALCON within MATLAB:

    addpath(genpath('MY_FALCON_DIRECTORY'));

This will be especially useful if you plan on running analyses 
similar to those found in the publication; otherwise, the top-level
directory should suffice.

#### Compiling the minimum-disjunction program: C-source method ####
Use this method if the provided binaries do not work for you
or there are none matching your operating system, and if you
aren't interested in experimenting with GPR-related aspects
of the algorithm.

[Download](http://sourceforge.net/projects/ats2-lang/files/ats2-lang/)
the latest ATS includes (click on the latest version,
e.g. ats2-postiats-x.y.z, then download the file
(ATS2-Postiats-include-x.y.z.tgz). Extract the archive, e.g.

    tar zxvf ATS2-Postiats-include-x.y.z.tgz

Set the environmental variable PATSHOME to the new directory:

    cd ATS2-Postiats-include-x.y.z
    export PATSHOME=`pwd`

Go in to the GPR subdirectory and build falcon_cmdline

    cd <MY_FALCON_DIRECTORY>/GPR
    make -C .. MYPORTDIR
    make PATSCC=gcc MYCCRULE=PORTABLE -f ../Makefile
    
And as stated above, copy or link to this file in your path as the file 'minDisj'.


#### Compiling the minimum-disjunction program: ATS source method ####
After installing ATS2 (see above), download the latest source
from this repository.

Go in to the GPR subdirectory and run

    make falcon_cmdline
    
And as stated above, copy or link to this file in your path as the file 'minDisj'.



### Additional information and related analyses ###

For other data and information related to FALCON and and associated
pubications, please see:

http://openwetware.org/wiki/Barker:Notebook/FALCON

In particular, data for testing the minimum-disjunction program [can be found there](http://openwetware.org/wiki/Barker:Notebook/FALCON#Yeast_RNA-Seq_analysis).
