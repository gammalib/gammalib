Installation
============


.. note ::

   These installation instructions apply to gammalib-00-04-10 and
   later. gammalib-00-04-10 has been built, installed and tested
   successfully on Debian, Ubuntu, Mandriva, OpenSUSE, Scientific Linux,
   CentOS, gentoo, Mac OS X, FreeBSD, and OpenSolaris (using gcc), so
   hopefully it also installs on your distribution. In case you encounter
   problems, please post a report on the bug tracker or send an e-mail to the
   `mailing list <mailto:gammalib-users@lists.soureforge.net>`_.

Getting GammaLib
----------------

The latest version of the GammaLib source code, documentation, and
example programs are available on the World-Wide Web from
`https://sourceforge.net/projects/gammalib/ <https://sourceforge.net/projects/gammalib/>`_.

Prerequisites
-------------

GammaLib should compile on every modern Unix system without any need to
install other libraries. To enable support for FITS file handling,
however, the cfitsio library from HEASARC needs to be installed. cfitsio
exists on many distributions as a prebuilt library, so in general you
can directly install it using your package manager. Make sure that you
install the library and the development package, as the cfitsio header
file (which usually comes in the development package) is needed to
compiling GammaLib.

If cfitsio is not available as a prebuilt package, or if you encounter
some problems with the prebuilt package, cfitsio can be downloaded from
`http://heasarc.gsfc.nasa.gov/fitsio <http://heasarc.gsfc.nasa.gov/fitsio>`_
and be installed from the source (for details, please refer to the
installation instructions on the
`HEASARC <http://heasarc.gsfc.nasa.gov/fitsio>`_ site). We recommend
installation of cfitsio in the default GammaLib install directory as a
shared library by typing ::

   $ ./configure --prefix=/usr/local/gamma
   $ make shared
   $ make install

GammaLib can also benefit from the presence of the readline library that
provides line-editing and history capabilities for text input (GammaLib
offers however also full functionality without having readline
installed). readline (which depends on ncurses) is available on most
system as a prebuilt library, so also here we recommend to use your
package manager to install the libraries. Also here, the readline (and
ncurses) development packages are required, so that the header files
become available and the symbolic links are set correctly.

If readline and ncurses are not available as a prebuilt package, or if
you encounter some problems with the prebuilt packages, both libraries
can be downloaded from
`http://ftp.gnu.org/gnu/ <http://ftp.gnu.org/gnu/>`_.

Building GammaLib
-----------------

Download the latest source tarball, save it in an appropriate location,
and type ::

   $ tar xvfz GammaLib-XX-XX-XX.tar.gz

where ``XX-XX-XX`` is the version number of the library. Step in the created
directory using ::

   $ cd gammalib-XX-XX-XX

and build GammaLib by typing ::

   $ ./configure
   $ make

at the operating system prompt. The configuration command customizes the
Makefiles for the particular system, the make command compiles the
source files and builds the library, and the make install command
installs the library in the install directory. Type ``./configure`` and not
simply configure to ensure that the configuration script in the current
directory is run and not some other system-wide configuration script. By
default, the install directory is set to ``/usr/local/gamma``. To change the
install directory an optional ``--prefix`` argument should be given, for
example ::

   $ ./configure --prefix=/usr/local

If cfitsio and/or readline is not installed in a standard location for
libraries (such as ``/usr/lib`` or ``/usr/local/lib``), the appropriate
location(s) can be specified using the ``LDFLAGS`` (for the library
directory) and ``CPPFLAGS`` (for the include directory) options::

   $ ./configure LDFLAGS='-L/opt/local/lib' CPPFLAGS='-I/opt/local/include'

A full list of configuration options can be found using ::

   $ ./configure --help

Testing GammaLib
----------------

GammaLib should be tested by typing::

   $ make check

This will execute an extensive testing suite that should terminate with ::

   ===================
   All 17 tests passed
   ===================

Eventually, loading the shared cfitsio and/or readline libraries may
fail during the test if the libraries are not located in standard
locations. In this case, add the library directories to the
``LD_LIBRARY_PATH`` environment variables (``DYLD_LIBRARY_PATH`` on Mac OS
X), e.g. ::

   export LD_LIBRARY_PATH=/opt/local/lib:$LD_LIBRARY_PATH

or ::

   setenv LD_LIBRARY_PATH /opt/local/lib:$LD_LIBRARY_PATH

on C shell variants.

Installing GammaLib
-------------------

Install GammaLib by typing ::

   $ make install

at the operating system prompt.

Setting up your environment
---------------------------

Before using GammaLib you have to setup some environment variables. This
will be done automatically by an initialisation script that will be
installed in the bin directory. Assuming that you have installed
GammaLib in the default directory ``/usr/local/gamma`` you need to add the
following to your ``$HOME/.bashrc`` or ``$HOME/.profile`` script on a Linux
machine::

   export GAMMALIB=/usr/local/gamma
   source $GAMMALIB/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your
``$HOME/.cshrc`` or ``$HOME/.tcshrc`` script::

   setenv GAMMALIB /usr/local/gamma
   source $GAMMALIB/bin/gammalib-init.csh

Installing documentation
------------------------

The GammaLib web pages are shipped together with the source code and
will be installed in the directory ``$(prefix)/share/doc/gammalib/html``,
where ``$(prefix)`` is the installation base path, by default
``/usr/local/gamma``.

For developers, a more detailed documentation can be built from the
source code using Doxygen. Doxygen will scan the source files for code
annotations, and compiles a complete documentation of the implemented
C++ classes in a set of www pages. In addition, it will create a set of
man files that can be used by the man command. If you have not installed
Doxygen on your machine, but you need the development documentation,
please visit
`http://dogygen.org <http://doxygen.org>`_
to download and install the Doxygen package.

Once Doxygen is install, type the following to build and to install the
Doxygen documentation::

   $ ./configure
   $ make doxygen
   $ make doxygen-install

The Doxygen documentation will be installed into
``$(prefix)/share/doc/gammalib/html/doxygen`` and can be browsed using a
regular web browser.

To check man support, type for example ::

   $ man GObservations

and you should see the documentation for the GObservations C++ class.

Known Problems
--------------

* Python support
   GammaLib comes with Python wrappers so that all classes can be directly
   used from Python. To compile-in Python support, GammaLib needs the
   Python.h header file, which on many distributions is not installed by
   default. To make Python.h available, install the Python developer
   package in your distribution using the package manager. Otherwise you
   will not be able to use GammaLib from Python.

* Mac OS X
   The Python development package is not installed as default on Mac OS X,
   and consequently, the Python.h header file is missing that is needed to
   compile in the Python bindings. The configure script recognises this
   fact and adjust the build procedure accordingly, but you will not be
   able to use GammaLib from Python. So better install the Python
   development package before installing GammaLib (see above).

* Solaris
   Although GammaLib builds on Solaris using the Sun compiler, there are
   problems with global symbols in shared libraries and exception catching,
   which prevents the FITS interface to work correctly. GammaLib has
   however been built and tested successfully using the GNU compiler, and
   this is the only build method that is currently supported. Problems have
   also been encountered when compiling cfitsio versions more recent than
   3.250. The problems have been reported to the cfitsio developer team,
   and are likely to be solved in the future. For the time being, it is
   recommended to use cfitsio version 3.250 on Solaris.

* OpenSolaris
   On OpenSolaris, the same problems concerning the SunStudio compiler
   occur as for Solaris, and also here, the GNU compiler is the recommended
   tool to build GammaLib. Also here, cfitsio version 3.250 is the
   recommended library as more recent version feature relocation
   problems. GammaLib has been tested using gcc 4.3.2 on OpenSolaris
   2009.06. Make sure to create the symbolic links ::

      $ ln -s /usr/bin/gcc4.3.2 /usr/bin/gcc
      $ ln -s /usr/bin/g++4.3.2 /usr/bin/g++

   which are not there by default to avoid excess warnings during
   compilation.

Getting Help
------------

Any questions, bug reports, or suggested enhancements related to
GammaLib should be submitted via the
`issue tracker <https://cta-redmine.irap.omp.eu/projects/gammalib>`_
or the
`mailing list <mailto:gammalib-users@lists.soureforge.net>`_.