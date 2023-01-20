GammaLib information
====================
* Version: 2.1.0.dev (20 January 2023)

[![Build Status](https://cta-jenkins.irap.omp.eu/buildStatus/icon?job=gammalib-integrate-os)](https://cta-jenkins.irap.omp.eu/job/gammalib-integrate-os/)

[![Quality Gate](https://cta-sonar.irap.omp.eu/api/badges/gate?key=gammalib)](https://cta-sonar.irap.omp.eu/dashboard/index/gammalib)


License information
===================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


What's new in this release?
===========================
See the files [NEWS](NEWS) and [ChangeLog](ChangeLog).


What is the GammaLib anyway?
============================
The GammaLib is a versatile toolbox for the high-level analysis of
astronomical gamma-ray data.  It is implemented as a C++ library
that is fully scriptable in the Python scripting language.  The library
provides core functionalities such as data input and output, interfaces
for parameter specifications, and a reporting and logging interface.
It implements instruments specific functionalities such as instrument
response functions and data formats.  Instrument specific functionalities
share a common interface to allow for extension of the GammaLib to
include new gamma-ray instruments.  The GammaLib provides an abstract
data analysis framework that enables simultaneous multi-mission analysis.


Web sites
=========
* http://cta.irap.omp.eu/gammalib/                  - for end users
* https://cta-redmine.irap.omp.eu/projects/gammalib - for developers


Prerequisites
=============
GammaLib should compile on every modern Unix system without any need
to install other libraries. 

To enable support for FITS file handling, however, the cfitsio library
from HEASARC needs to be installed.  cfitsio can be downloaded from 
http://heasarc.gsfc.nasa.gov/fitsio and detailed installation instructions
can be found there.  If cfitsio does not already exist on your system,
we recommend installation of cfitsio in the default GammaLib install
directory as a shared library by typing:

     $ ./configure --prefix=/usr/local/gamma
     $ make shared
     $ make install

GammaLib can also benefit from the presence of the readline library
that provides line-editing and history capabilities for text input
(GammaLib offers however also full functionality without having 
readline installed).  readline requires the ncurses library.  Both
libraries can be downloaded from http://ftp.gnu.org/gnu/.


Conda Installation
==================
The easiest is to install GammaLib via conda.  This also takes care of the
installation of cfitsio.  Assuming that you have installed anaconda, type
the following:

     $ conda config --append channels conda-forge
     $ conda config --append channels cta-observatory
     $ conda install gammalib


Linux Installation
==================
To build, verify and install GammaLib, simply type the following:

     $ ./configure
     $ make
     $ make check
     $ make install

If the folder does not contain any `configure` file, please run

     $ ./autogen.sh 

before invoking `configure`.

By default GammaLib installs itself in `/usr/local/gamma`.  If you need to
install GammaLib in a different location or in your home directory, use
the `--prefix` option to `./configure`.  For example:

     $ ./configure --prefix=/home/yourname/projects
     $ make
     $ make check
     $ make install

The file INSTALL details more about using configure. Also try

     $ ./configure --help.

The `make check` command will run an extensive unit test to verify that
GammaLib was correctly built.  Make sure that all tests were successful. 


Macintosh OS Installation
==========================
GammaLib is known to work on various flavors of Mac OS.  To cope with
different system versions and architectures, there are two Mac
specific configure options:

     $ ./configure --enable-universalsdk[=PATH]

creates a universal build of GammaLib.  The optional argument specifies
which MacOS SDK should be used to perform the build.  This defaults to
`/Developer/SDKs/MacOSX.10.4u.sdk`.  Specify `/` when building on a 10.5
system or higher, especially when building 64-bit code.

     $ ./configure --with-univeral-archs=VALUE

specifies the kind of universal build that should be created.  Possible
values are: `32-bit`, `3-way`, `intel` or `all`.  By default, a `32-bit` 
build will be made.  This option is only valid when
`--enable-universalsdk` is specified.

These options are in particular needed if your Python architecture differs
from the default architecture of your system.  To examine the Python
architecture you may type:

     $ file `which python`

which will return the architectures that are compiled in the Mach-0
executable:

     i386    32-bit intel
     ppc     32-bit powerpc
     ppc64   64-bit powerpc
     x86_64  64-bit intel

If Python is 32-bit (ppc, i386) but the compiler produces by default
64-bit code (ppc64, x86_64), the Python module will not work.  Using

     $ ./configure --enable-universalsdk=/

will force a universal 32-bit build which creates code for ppc and 
i386.  If on the other hand Python is 64-bit (ppc64, x86_64) but the
compiler produces by default 32-bit code (ppc, i386), the option

     $ ./configure --enable-universalsdk=/ --with-univeral-archs=3-way

will generate a universal build which contains 32-bit and 64-bit code.


BSD Installation
================
GammaLib has been tested on FreeBSD successfully.  Follow the Linux
installation instructions above.


Solaris Installation
====================
GammaLib compiles on Solaris and openSolaris using the GNU compiler
collection.  Follow the Linux installation instructions above.
GammaLib also compiles using the Sun Studio compiler, but due to
some unexpected handling of global variables, it is not fully
working (see "Known problems" below).
Please note that these is also a cfitsio problem on Solaris, please
refer to the section "Known problems" for more information.


Windows Installation
====================
On Windows GammaLib needs to be installed into a virtual machine running
a Linux distribution.


Testing
=======
After building and before installing GammaLib, you should run the
extensive unit test by typing:

    $ make check

If everything works successfully you should see

    ===================
    All 22 tests passed
    ===================

or

    ============================================================================
    Testsuite summary for gammalib 2.1.0.dev
    ============================================================================
    # TOTAL: 22
    # PASS:  22
    # SKIP:  0
    # XFAIL: 0
    # FAIL:  0
    # XPASS: 0
    # ERROR: 0
    ============================================================================

at the end of the test (depending on your automake version).  If no Python
support was compiled in, the number of tests performed will be reduced by one.


Setting up your environment
===========================
Before using GammaLib you have to setup some environment variables.
This will be done automatically by an initialisation script that will
be installed in the bin directory.

Assuming that you have installed GammaLib in the default directory 
`/usr/local/gamma` you need to add the following to your `$HOME/.bashrc` or 
`$HOME/.profile` script on a Linux machine:

    export GAMMALIB=/usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your
`$HOME/.cshrc` or `$HOME/.tcshrc` script:

    setenv GAMMALIB /usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.csh


Getting started
===============
The easiest way to start with GammaLib is by using the python interface.
To start, type the following:

    $ python
    >>> import gammalib
    >>> models=gammalib.GModels()
    >>> print models
    === GModels ===
    Number of models ..........: 0
    Number of parameters ......: 0

This examples allocates an empty GModels object that holds a collection
of models.

For examples, inspect the test directory, and the test directories
of the instrument specific interfaces, e.g. `inst/mwl/test`, `inst/lat/test`,
and `inst/cat/test`.


Documentation
=============
The doc directory (usually at `/usr/local/gamma/share/doc/gammalib`)
contains the most recent set of updated documentation for this release.
A detailed documentation can be created by typing:

    $ make doc

before installing the library. Two types of documentation exist: 
* code documentation
* user documentation

Code documentation is created using Doxygen.  You need Doxygen on your
system to generate the code documentation.  This includes man pages.
Doxygen can be obtained from http://www.stack.nl/~dimitri/doxygen/.

User documentation is created using Sphinx.  You need Sphinx on your
system to generate the user documentation.  Sphinx can be obtained
from http://sphinx-doc.org/install.html.


Bug reports
===========
To report or search for bugs, please use the GammaLib Bug Tracker at
https://cta-redmine.irap.omp.eu/projects/gammalib.  Before using the
tracker, please read 
https://cta-redmine.irap.omp.eu/projects/gammalib/wiki/Submission_guidelines


Known problems
==============

Python support
--------------

GammaLib comes with Python wrappers so that all classes can be
directly used from Python.  To compile-in Python support, GammaLib
needs the Python.h header file, which on many distributions is not
installed by default.
To make Python.h available, install the Python developer package
in your distribution.  Otherwise you will not be able to use GammaLib
from Python.

Readline support
----------------
Many distributions do not have the readline header files and symbolic
link to the shared library set up by default.  The same is true for the
ncurses library that is needed by readline.  To properly compile-in
readline support in GammaLib, make sure that the redline header files
exist and that the symbolic links are set.  In many distributions,
the shared redline and ncurses are located in `/lib`, while the symbolic
links are in `/usr/lib`.  Make sure that you have symbolic links with
names

     libreadline.so
     libncurses.so

Symbolic links with an additional number attached, such as
`libreadline.so.6` or `libncurses.so.5` (which are found in many
distributions in the `/lib` directory) are not sufficient.
In many distributions, the appropriate headers and symbolic links
are installed if the proper readline and ncurses developer
packages are installed.

Solaris
-------

Although GammaLib builds on Solaris using the Sun compiler, there are
problems with global symbols in the shared library that prevent the
model registry to work correctly.  Furthermore, GammaLib is not able
to catch its own exceptions, which prevents the FITS interface to work
correctly.
GammaLib has however been built and tested successfully using the GNU
compiler, and this is the only build method that is currently supported.
Problems have also been encountered when compiling cfitsio versions
more recent than 3.250.  The problems have been reported to the cfitsio
developer team, and are likely to be solved in the future.  For the
time being, it is recommended to use cfitsio version 3.250 on Solaris.

OpenSolaris
-----------
On OpenSolaris, the same problems concerning the SunStudio compiler
occur as for Solaris, and also here, the GNU compiler is the recommended
tool to build GammaLib.  Also here, cfitsio version 3.250 is the
recommended library as more recent version feature relocation
problems. GammaLib has been tested using gcc 4.3.2 on 
OpenSolaris 2009.06.  Make sure to create the symbolic links

     ln -s /usr/bin/gcc4.3.2 /usr/bin/gcc
     ln -s /usr/bin/g++4.3.2 /usr/bin/g++

which are not there by default to avoid excess warnings during
compilation. 
For OpenSolaris 2009.06, no readline package is available, hence 
readline support is not readily available.  Readline support can however
be enabled by directly installing the GNU ncurses and readline packages
from source (http://ftp.gnu.org/gnu).  This has been done successfully
on OpenSolaris 2009.06 for ncurses-5.9 and readline-6.2.  Make sure
to specify the `--with-shared` option when configuring ncurses, so
that the shared library will be built and installed.  Otherwise,
relocation errors may occur during compilation of GammaLib.


Contact
=======
To get in touch with the GammaLib developers and to contribute to the
project please contact Juergen Knoedlseder <jurgen.knodlseder@irap.omp.eu>.
