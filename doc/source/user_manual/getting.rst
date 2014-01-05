Getting GammaLib
================

Before you start
----------------

The procedure for building and installing GammaLib is modeled on GNU software
distributions. You do not need to have system administrator privileges
to compile and to install GammaLib.

You will need the following to build the software:

-  About 100 MB of free disk space.

-  An ANSI C++ compiler. We recommend building GammaLib with the GNU g++
   compiler.

-  GNU make

-  The cfitsio library for FITS file support together with the developer
   package that includes the ``cfitsio.h`` header.

Note that GammaLib compiles also in the absence of the cfitsio library, yet
without cfitsio, FITS file reading or writing is not supported.

Furthermore, the following optional packages are supported but are not
required to compile :

-  Python, including the Python developer package that includes the
   ``Python.h`` header file. If Python is present, the GammaLib Python module will
   be built and installed, allowing to script all GammaLib functionalities from
   within Python.

-  readline, including the readline developer package that provides the
   ``readline.h`` header file. If readline is present, the packages are used
   to enhance the user interface when entering parameters for ftools
   applications (see section [sec:app]).

If you plan to modify or to extend the GammaLib source code, the following
software is also required on your system:

-  GNU `autoconf <http://www.gnu.org/software/autoconf/>`_ and `automake
   <http://www.gnu.org/software/automake/>`_ is needed to rebuild the
   configure script and ``Makefile.am`` resource files following
   configuration modifications.

-  `swig <http://www.swig.org/>`_ is needed to rebuild the Python wrappers
   following Python interface modifications. Make sure to install the
   latest swig version () to guarantee the largest possible
   compatibility of the Python wrappers.

-  `Doxygen <http://www.doxygen.org/>`_ is needed to rebuild the software
   reference manual following code modifications.

The following sections provide some information about the installation
of cfitsio and readline.

.. _sec_cfitsio:

Installing cfitsio
~~~~~~~~~~~~~~~~~~

HEASARC's cfitsio library comes on many Linux distributions as
pre-compiled binary, and there are good chances that the package is
already installed on your system. For Mac OS X, cfitsio can be installed
from Mac Ports. If you use a pre-compiled binary, make sure that also
the developer package is installed on your system. The developer package
provides the ``cfitsio.h`` header file which is needed to compile in FITS
file support in GammaLib. Please refer to the documentation of your Linux
distribution to learn how to install pre-compiled binary packages (note
that the installation of pre-compiled binary packages usually requires
system administrator privileges).

If you need (or prefer) to install cfitsio from source, you can download
the latest source code from http://heasarc.gsfc.nasa.gov/fitsio.
Detailed installation instructions can also be found on this site. We
recommend that you install cfitsio as a shared library in the same
directory in which you will install , so that cfitsio is automatically
found by the GammaLib configure script. By default, GammaLib gets installed
into the directory ``/usr/local/gamma``.

You can install version 3.290 of cfitsio (the latest version that was
available during writing this manual) by executing the following command
sequence ($ denotes the UNIX shell prompt)::

    $ wget ftp://heasarc.gsfc.nasa.gov/software/fitsio/c/cfitsio3290.tar.gz
    $ tar xfz cfitsio3290.tar.gz
    $ cd cfitsio
    $ ./configure --prefix=/usr/local/gamma
    $ make shared
    $ sudo make install

The ``--prefix=/usr/local/gamma`` option specifies the directory into which
cfitsio gets installed. We choose here the default GammaLib installation
directory ``/usr/local/gamma``. As this directory is a system directory, we
need to use sudo for installation. If you decide to install cfitsio into
a local directory which is owned by yourself, it is sufficient to type
make install to install the library.

Installing readline
~~~~~~~~~~~~~~~~~~~

The readline package comes on all Linux distributions that are known to
us as pre-compiled binary, and it is almost certain that readline is already
installed on your system. Very often, however, the readline developer package
that provides the readline.h header file is not installed, and you need to
install this package yourself to enable readline support for GammaLib. Please
refer to the documentation of your Linux distribution to learn how to
install pre-compiled binary packages (note that the installation of
pre-compiled binary packages usually requires system administrator
privileges).

If you need (or prefer) to install readline from source, you need also
to install the ncurses library that is required by readline. Here is the
command line sequence that will install ncurses (version 5.9) and
readline (version 6.2) in the GammaLib default install directory
``/usr/local/gamma`` from source::

    $ wget http://ftp.gnu.org/gnu/ncurses/ncurses-5.9.tar.gz
    $ tar xfz ncurses-5.9.tar.gz
    $ cd ncurses-5.9
    $ ./configure --prefix=/usr/local/gamma
    $ make
    $ sudo make install
    $ cd ..
    $ wget http://ftp.gnu.org/gnu/readline/readline-6.2.tar.gz
    $ tar xfz readline-6.2.tar.gz
    $ cd readline-6.2
    $ ./configure --prefix=/usr/local/gamma
    $ make
    $ sudo make install

Note that sudo is only needed if you are not the owner of the install
directory.

Installing
----------

Downloading
~~~~~~~~~~~

To get the latest version of GammaLib, please visit the site
https://sourceforge.net/projects/gammalib/. The code can be downloaded
from this site by clicking on the download button. Alternatively, the
code can be downloaded and unpacked from the UNIX prompt using::

    $ wget --no-check-certificate https://downloads.sourceforge.net/project/gammalib/
    gammalib/gammalib-00-08-00.tar.gz
    $ tar xfz gammalib-00-08-00.tar.gz

The GammaLib source code can also be cloned using git. This method is
recommended if you plan to contribute to the development of the GammaLib
library. Assuming that git is installed on your system, you may clone
GammaLib using::

    $ git clone https://cta-git.irap.omp.eu/gammalib

In case that you get::

    error: SSL certificate problem, verify that the CA cert is OK.

you may add::

    $ export GIT_SSL_NO_VERIFY=true

before retrieving the code.

.. _sec_configure:

Configuring
~~~~~~~~~~~

Once you've downloaded and uncompressed GammaLib, step into the GammaLib
source code directory and type ::

    $ ./configure

to configure the library for compilation. Make sure that you type
``./configure`` and not simply configure to ensure that the configuration
script in the current directory is invoked and not some other
system-wide configuration script.

If you would like to install GammaLib in a different directory, use the optional
``--prefix`` argument during the configuration step. For example ::

    $ ./configure --prefix=/home/myname/gamma

installs GammaLib in the gamma directory that will be located in the user's
myname home directory. You can obtain a full list of configuration
options using ::

    $ ./configure --help

If configuration was successful, the script will terminate with printing
information about the configuration. This information is important in
case that you encounter installation problems, and may help you to
diagnose the problems. The typical output that you may see is as
follows::

      GammaLib configuration summary
      ==============================
      * FITS I/O support             (yes)   /usr/local/gamma/lib /usr/local/gamma/include
      * Readline support             (yes)    
      * Ncurses support              (yes)   
      * Python                       (yes)
      * Python.h                     (yes)
      * swig                         (yes)
      * Make Python bindings         (yes)
      * Multiwavelength interface    (yes)
      * Fermi-LAT interface          (yes)
      * CTA interface                (yes)
      * Doxygen                      (yes)   /usr/local/bin/doxygen
      * Perform NaN/Inf checks       (yes)   (default)
      * Perform range checking       (yes)   (default)
      * Optimize memory usage        (yes)   (default)
      - Compile in debug code        (no)    (default)
      - Enable code for profiling    (no)    (default)

The script informs whether cfitsio has been found (and eventually also
gives the directories in which the cfitsio library and the header file
resides), whether readline and ncurses have been found, and whether
Python including the Python.h header file is available. Although none of
these items is mandatory, we highly recommend to install cfitsio to
support FITS file reading and writing (see section [sec:cfitsio]), and
to install Python to enable GammaLib scripting.

If cfitsio is installed on your system but not found by the configure
script, it may be located in a directory that is not known to the
configure script. By default, configure will search for cfitsio (in the
given order) in the GammaLib install directory, in all standard paths (e.g.
``/usr/lib``, ``/usr/local/lib``, ...), and in some system specific locations,
including ``/opt/local/lib`` for Mac OS X. Assuming that you installed
cfitsio on your system in the directory ``/home/myname/cfitsio``, you may
explicitly specify this location to configure using the ``LDFLAGS`` and
``CPPFLAGS`` environment variables::

    $ ./configure LDFLAGS=-L/home/myname/cfitsio/lib CPPFLAGS=-I/home/myname/cfitsio/include

Here, ``LDFLAGS`` specifies the path where the shared cfitsio library is
located, while ``CPPFLAGS`` specifies the path where the ``cfitsio.h`` header
file is located. Note that ``-L`` has to prefix the library path and that ``-I``
has to prefix the header file path. With the same method, you may
specify any non-standard location for the readline and ncurses
libraries.

The configuration script also checks for the presence of swig, which is
used for building the Python wrapper files. Normally, swig is not needed
to create the Python bindings as the necessary wrapper files are shipped
with the GammaLib source code. If you plan, however, to modify or to extend the
Python interface, you will need swig to rebuild the Python wrappers
following changes to the interface.

The configuration summary informs also about all instrument dependent
interfaces that will be compiled into the GammaLib library. By default, all
available interfaces (multi-wavelength, *Fermi*-LAT, COMPTEL and CTA) will be
compiled into GammaLib. If you wish to disable a particular interface, you may
use the configure options ``--without-mwl``, ``--without-lat``, 
``--without-com`` or ``--without-cta``.
For example, ::

    $ ./configure --without-mwl --without-lat --without-com

will compile GammaLib without the multi-wavelength, the *Fermi*-LAT and 
the COMPTEL interfaces. In this case, only CTA data analysis will be supported.

GammaLib uses `Doxygen <http://www.doxygen.org/>`_ for code documentation, 
and the latest GammaLib reference manual can be found at
http://gammalib.sourceforge.net/doxygen/. In case that you want to
install the reference manual also locally on your machine, Doxygen is
needed to create the reference manual from the source code. Doxygen is
also needed if you plan to modify or extend the GammaLib library to allow
rebuilding the reference documentation after changes. Please read see
section [sec:doxygen] to learn how to build and to install the reference
manual locally.

Finally, there exist a number of options that define how exactly GammaLib
will be compiled.

By default, GammaLib makes use of OpenMP for multi-core processing. If you 
want to disable the multi-core processing, you may specify the 
``--disable-openmp`` option during configuration.

Several methods are able to detect invalid floating point values (either
``NaN`` or ``Inf``), and by default, these checks will be compiled in the
library to track numerical problems. If you want to disable these
checks, you may specify the ``--disable-nan-check`` option during
configuration.

Range checking is performed by default on all indices that are provided
to methods or operators (such as vector or matrix element indices, sky
pixels, event indices, etc.), at the expense of a small speed penalty
that arises from these verifications. You may disable these range
checkings by specifying the ``--disable-range-check`` option during
configuration.

In a few places there exists a trade-off between speed and memory
requirements, and a choice has to be made whether faster execution or
smaller memory allocation should be preferred. By default, smaller
memory allocation is preferred by GammaLib, but if you are not concerned about
memory allocation you may specify the ``--disable-small-memory`` option
during configuration to speed up the code.

If you develop code for GammaLib you may be interested in adding some special
debugging code, and this debugging code can be compiled in the library
by specifying the ``--enable-debug`` option during configuration. By default,
no debugging code will be added to GammaLib.

Another developer option concerns profiling, which may be of interest to
optimize the execution time of your code. If you would like to add
profiling information to the code (which will be at the expense of
execution time), you may specify the ``--enable-profiling`` option during
configuration, which adds the ``-pg`` flags to the compiler. By default,
profiling is disabled for GammaLib.

Mac OS X options
^^^^^^^^^^^^^^^^

The Mac OS X environment is special in that it supports different CPU
architectures (intel, ppc) and different addressing schemes (32-bit and
64-bit). To cope with different system versions and architectures, you
can build a universal binary by using the option ::

    $ ./configure --enable-universalsdk[=PATH]

The optional argument ``PATH`` specifies which OSX SDK should be used to
perform the build. By default, the SDK ``/Developer/SDKs/MacOSX.10.4u.sdk``
is used. If you want to build a universal binary on Mac OS X 10.5 or
higher, and in particular if you build 64-bit code, you have to specify
``--enable-universalsdk=/``.

A second option (which is only valid in combination with the
``--enable-universalsdk``) allows to specify the kind of universal build that
should be created::

    $ ./configure --enable-universalsdk[=PATH] --with-univeral-archs=VALUE

Possible options for ``VALUE`` are: ``32-bit``, ``3-way``, ``intel``, or ``all``. By
default, a 32-bit build will be made.

These options are in particular needed if your Python architecture
differs from the default architecture of your system. To examine the
Python architecture you may type::

    $ file `which python`

which will return the architectures that are compiled in the Python
executable::

      i386     32-bit intel
      ppc      32-bit powerpc
      ppc64    64-bit powerpc
      x86_64   64-bit intel

If Python is 32-bit (``ppc``, ``i386``) but the compiler produces by default
64-bit code (``ppc64``, ``x86_64``), the Python module will not work. Using ::

    $ ./configure --enable-universalsdk=/

will force a universal 32-bit build which creates code for ``ppc`` and ``i386``
architectures. If on the other hand Python is 64-bit (``ppc64``, ``x86_64``)
but the compiler produces by default 32-bit code (``ppc``, ``i386``), the option ::

    $ ./configure --enable-universalsdk=/ --with-univeral-archs=3-way

will generate a universal build which contains 32-bit and 64-bit code.

Building
~~~~~~~~

Once configured you can build GammaLib by typing ::

    $ make

This compiles all GammaLib code, including the Python wrappers, and builds the
dynamic library and Python module.

GammaLib building can profit from multi-processor or multi-core machines by
performing parallel compilation of source code within the modules. You
can enable this feature by typing ::

    $ make -j<n>

where ``<n>`` is a number that should be twice the number of cores or
processors that are available on your machine.

In case that you rebuild GammaLib after changing the configuration, we recommend
to clean the directory from any former build by typing ::

    $ make clean

prior to make. This will remove all existing object and library files
from the source code directory, allowing for a fresh clean build of the
library.

Testing
~~~~~~~

GammaLib comes with an extensive unit test that allows to validate the library
prior to installation. **We highly recommend to run this unit test
before installing the library (see section [sec:install]).**

To run the unit test type::

    $ make check

This will start a test of all GammaLib modules by using dedicated executables
which will print some progress and success information into the
terminal. After completion of all tests (and assuming that all
instrument dependent modules are enabled), you should see the following
message in your terminal::

    ===================
    All 19 tests passed
    ===================

(Note that the exact number of tests that is conducted depends on the 
configuration options).

.. _sec_install:

Installing
~~~~~~~~~~

GammaLib is finally installed by typing ::

    $ [sudo] make install

By default, GammaLib is installed in the system directory ``/usr/local/gamma``,
hence sudo needs to be prepended to enable writing in a system-level
directory. If you install GammaLib, however, in a local directory of which you
are the owner, or if you install GammaLib under root, you may simply specify
make install to initiate the installation process.

The installation step will copy all necessary files into the
installation directory. Information will be copied in the following
subdirectories:

-  ``bin`` contains GammaLib environment configuration scripts (see section
   [sec:environment])

-  ``include`` contains GammaLib header files (subdirectory gammalib)

-  ``lib`` contains the GammaLib library and Python module

-  ``share`` contains addition GammaLib information, such as a calibration database
   (subdirectory ``caldb``), documentation (subdirectory ``doc``), and Python
   interface definition files (subdirectory ``gammalib/swig``)

.. _sec_environment:

Setting up the GammaLib environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before using GammaLib you have to setup some environment variables. This will be
done automatically by an initialisation script that has been installed
in the bin subdirectory of the install directory. Assuming that you have
installed GammaLib in the default directory ``/usr/local/gamma`` you need to
add the following to your ``$HOME/.bashrc`` or ``$HOME/.profile`` script on a
Linux machine::

    export GAMMALIB=/usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your
``$HOME/.cshrc`` or ``$HOME/.tcshrc`` script::

    setenv GAMMALIB /usr/local/gamma
    source $GAMMALIB/bin/gammalib-init.csh

You then have to source your initialisation script by typing (for
example) ::

    $ source $HOME/.bashrc

and all environment variables are set correctly to use GammaLib properly.

.. _sec_doxygen:

Generating the reference documentation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reference documentation for GammaLib is generated directly from the source
code using the `Doxygen <http://www.doxygen.org/>`_ documentation system.
The latest GammaLib reference manual can be found at
http://gammalib.sourceforge.net/doxygen/.

The reference documentation is not shipped together with the source code
as this would considerably increase the size of the tarball. In case
that you want to install the reference manual also locally on your
machine, you first have to create the documentation using Doxygen.

Assuming that Doxygen is available on your machine (see section
:ref:`sec_configure`) you can create the reference documentation by typing ::

    $ make doxygen

Once created, you can install the reference manual by typing ::

    $ [sudo] make doxygen-install

By default, GammaLib is installed in the system directory ``/usr/local/gamma``,
hence sudo needs to be prepended to enable writing in a system-level
directory. If you install , however, in a local directory of which you
are the owner, or if you install GammaLib under root, you may simply specify
make install to initiate the installation process.

The reference manual will be installed in form of web-browsable HTML
files into the folder ::

      /usr/local/gamma/share/doc/gammalib/html/doxygen

You can access all web-based GammaLib documentation locally using
``file:///usr/local/gamma/share/doc/gammalib/html/index.html`` (assuming
that the GammaLib library has been installed in the default directory
``/usr/local/gamma``).

In addition, the reference manual will also be available as man pages
that will be installed into ::

      /usr/local/gamma/share/doc/gammalib/man

To access for example the information for the ``GApplication`` class, you
can type ::

    $ man GApplication

which then returns the synopsis and detailed documentation for the
requested class.

Getting support
---------------

Any question, bug report, or suggested enhancement related to GammaLib should be
submitted via the Tracker on https://cta-redmine.irap.omp.eu/projects/gammalib
or by sending an e-mail to the mailing list.

.. _sec_known_problems:

Known problems
--------------

Solaris (TBW)
