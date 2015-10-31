.. _installation:

Installation
============

.. note ::

   These installation instructions apply to gammalib-0.11.0 and
   later. gammalib-0.11.0 has been built, installed and tested
   successfully on Debian, Ubuntu, Mandriva, OpenSUSE, Scientific Linux,
   CentOS, Mac OS X, FreeBSD, and OpenSolaris (using gcc), so
   hopefully it also installs on your distribution. In case you encounter
   problems, please post a report on the bug tracker or send an e-mail to the
   `mailing list <mailto:ctools@irap.omp.eu>`_.

.. _getting:

Getting GammaLib
----------------

The latest version of the GammaLib source code, documentation, and
example programs are available on the World-Wide Web from
`http://cta.irap.omp.eu/ctools/download.html <http://cta.irap.omp.eu/ctools/download.html>`_.

.. _prerequisits:

Prerequisites
-------------

GammaLib should compile on every modern Unix system and on Mac OS X,
Windows is not supported. To enable support for FITS file handling
the cfitsio library from HEASARC needs to be installed. cfitsio exists
on many Linux distributions as a prebuilt library, so in general you
can directly install it using your package manager. Make sure that you 
install the library and the development package, as the cfitsio header
file (which usually comes only in the development package) is needed 
to compiling GammaLib.

If cfitsio is not available as a prebuilt package, or if you encounter
some problems with the prebuilt package, cfitsio can be downloaded from
`http://heasarc.gsfc.nasa.gov/fitsio <http://heasarc.gsfc.nasa.gov/fitsio>`_
and installed from the source files (for details, please refer to the
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
package manager to install the libraries if they are not already there.
Also here, the readline (and ncurses) development packages are required,
so that the header files become available.

If readline and ncurses are not available as a prebuilt package, or if
you encounter some problems with the prebuilt packages, both libraries
can be downloaded from
`http://ftp.gnu.org/gnu/ <http://ftp.gnu.org/gnu/>`_
and installed from the source files. We recommend to put also those in the
default GammaLib install directory ``/usr/local/gamma``.

.. _build:

Building GammaLib
-----------------

.. note ::

   You will now learn how to build the GammaLib library from the source files.
   For Mac OS X there exist also binary packages that are bundled together 
   with ctools, and we recommend using directly these binary packages. The
   Mac OS X binary packages can be downloaded 
   `here <http://cta.irap.omp.eu/ctools/download.html>`_.

.. note ::

   If you need the latest source code from the git repository that has not 
   yet been released you may also clone the gammalib git repository using ::

   $ git clone https://cta-git.irap.omp.eu/gammalib

   In case you encounter an SSL certificat problem, type ::

   $ export GIT_SSL_NO_VERIFY=true

   before cloning from git. Step then in the cloned ``gammalib`` directory
   and type ::

   $ git checkout devel
   $ ./autogen.sh

   This will make sure that you use the validated development trunk and 
   creates the configuration file that is needed to configure 
   GammaLib before compilation. Continue then with the procedure below
   starting with ``./configure``.

To build GammaLib from the source files, download the latest release
tarball, save it in an appropriate location,
and type ::

   $ tar xvfz gammalib-x.y.z.tar.gz

where ``x.y.z`` is the version number of the library. Step in the created
directory using ::

   $ cd gammalib-x.y.z

and build GammaLib by typing ::

   $ ./configure
   $ make

at the operating system prompt. The configuration command customizes the
Makefiles for your particular system, the make command compiles the
source files and builds the C++ library and Python module. Type
``./configure`` and not simply ``configure`` to ensure that the configuration
script in the current directory is used and not some other system-wide
configuration script. By default, the install directory is set to 
``/usr/local/gamma``. To change the install directory an optional 
``--prefix`` argument should be given, for example ::

   $ ./configure --prefix=/usr/local

If cfitsio and/or readline is not installed in a standard location for
libraries (such as ``/usr/lib`` or ``/usr/local/lib``), the appropriate
location(s) can be specified using the ``LDFLAGS`` (for the library
directory) and ``CPPFLAGS`` (for the include directory) options::

   $ ./configure LDFLAGS='-L/opt/local/lib' CPPFLAGS='-I/opt/local/include'

A full list of configuration options can be found using ::

   $ ./configure --help

.. _test:

Testing GammaLib
----------------

Before you install GammaLib you should test the C++ library and Python 
module by typing::

   $ make check

This will execute an extensive testing suite that should terminate with ::

   ===================
   All 20 tests passed
   ===================

Eventually, loading the shared cfitsio and/or readline libraries may
fail during the test if the libraries are not located in standard
locations. In this case, add the library directories to the
``LD_LIBRARY_PATH`` environment variables (``DYLD_LIBRARY_PATH`` on Mac OS
X), e.g. ::

   $ export LD_LIBRARY_PATH=/opt/local/lib:$LD_LIBRARY_PATH

.. _install:

Installing GammaLib
-------------------

Now you are ready to install GammaLib by typing ::

   $ make install

at the operating system prompt. You may need to prepend a ``sudo`` in
case that you need administrator privileges to access the install
directory. If you do not have such privileges, chose an install directory
that you can access using the ``--prefix`` option.

.. _setup:

Setting up your environment
---------------------------

Before using GammaLib you have to setup some environment variables. This
will be done automatically by an initialisation script that will be
installed in the bin directory. Assuming that you have installed
GammaLib in the default directory ``/usr/local/gamma`` you need to add the
following to your ``$HOME/.bashrc`` or ``$HOME/.profile`` script on a Linux
machine:

.. code-block:: bash

   export GAMMALIB=/usr/local/gamma
   source $GAMMALIB/bin/gammalib-init.sh

If you use C shell or a variant then add the following to your
``$HOME/.cshrc`` or ``$HOME/.tcshrc`` script:

.. code-block:: csh

   setenv GAMMALIB /usr/local/gamma
   source $GAMMALIB/bin/gammalib-init.csh

.. _documentation:

Installing documentation
------------------------

.. note ::

   The documentation of the latest GammaLib release can be found at
   `http://cta.irap.omp.eu/gammalib/ <http://cta.irap.omp.eu/gammalib/>`_.
   The documentation corresponding to the git development branch can be
   found at
   `http://cta.irap.omp.eu/gammalib-devel/ <http://cta.irap.omp.eu/gammalib-devel/>`_.
   You can however also install the GammaLib documentation locally on your
   machine, and this section describes how to do that.

The GammaLib documentation is shipped together with the source code and
will be installed in the directory ``$(prefix)/share/doc/gammalib/html``,
where ``$(prefix)`` is the installation base path, by default
``/usr/local/gamma``. This comprises user documentation and code
documentation.

To build the user documentation you need the Sphinx reStructuredText
documentation generator installed
(see `http://sphinx-doc.org/rest.html <http://sphinx-doc.org/rest.html>`_
for more information).
Code documentation is based on Doxygen, which also needs to be installed
on your system
(see `http://dogygen.org <http://doxygen.org>`_ to download and install
the Doxygen package).
Doxygen will scan the source files for code annotations, and compiles a
complete documentation of the implemented C++ classes in a set of html
pages. In addition, it will create a set of man files that can be accessed 
using the ``man`` command.

To build and install all documentation, type the following::

   $ ./configure
   $ make doc
   $ make install

To build only user documentation, type::

   $ ./configure
   $ make sphinx
   $ make install

and to build only Doxygen documentation, type::

   $ ./configure
   $ make doxygen
   $ make install

The Doxygen documentation will be installed into
``$(prefix)/share/doc/gammalib/html/doxygen`` and can be browsed using a
regular web browser.

To check man support, type for example ::

   $ man GObservations

and you should see the documentation for the GObservations C++ class.

.. _issues:

Known issues
------------

.. topic:: Many (but not all) unit tests fails

   Occasionally it may happen that the cfitsio library is not found when
   configuring GammaLib. The library will compile successfully without
   cfitsio, but in that case FITS I/O will not be supported. Consequently,
   many unit tests will fail. If you are sure that cfitsio is installed,
   but the path where the library and the path where the ``fitsio.h`` header
   reside are non-standard, you may add the paths explicitly during
   configuration using::

      $ ./configure LDFLAGS='-L/path/to/cftsio/library' CPPFLAGS='-I/path/to/fitsio.h/header'

   The same logic applies for finding the readline and ncurses libraries,
   although these libraries are not manadatory for getting the full
   GammaLib functionalities.

.. topic:: Python module does not work

   GammaLib includes a Python module that is built from so called wrapper 
   files that are autogenerated using the `swig <http://www.swig.org/>`_
   tool. These wrapper files are shipped with a GammaLib release, but if
   you use the code from git you need `swig <http://www.swig.org/>`_
   to generate the wrapper files during the build step. In any case,
   to compile the Python module GammaLib needs the ``Python.h`` header file
   which may not necessarily be installed on your system. Check the output
   of ``./configure`` to examine the configuration that GammaLib has
   detected. You may see the following::

    * Python                       (yes)
    * Python.h                     (yes)
    * Python wrappers              (yes)
    * swig                         (yes)

   Recall, if the wrappers exist you do not need `swig <http://www.swig.org/>`_,
   but if the wrappers don't exist you need `swig <http://www.swig.org/>`_.
   If the ``Python.h`` header file does not exist then install the Python
   development package.

.. topic:: Python unit test fail on El Capitan

    On Mac OS X 10.11 (El Capitan) the Python unit tests fail when GammaLib 
    is not installed due to the way the ``DYLD_LIBRARY_PATH`` environment
    variable is handled by the operating system. A (not very elegant)
    workaround is to install GammaLib before running the unit test using 
    ``make check``.


.. topic:: Mac OS X

   Dependent on the Mac OS X version you are using, not everything that
   is needed to install GammaLib may be available (e.g. automake, libtool, 
   etc...). The easiest way to get the needed software is using a package 
   management system such as `MacPorts <https://www.macports.org/>`_ or 
   `Homebrew <http://brew.sh/>`_. On a fresh El Capitan install you need 
   for example the following from `Homebrew <http://brew.sh/>`_::

     $ brew install automake
     $ brew install libtool
     $ brew install cfitsio
     $ brew install swig

   `swig <http://www.swig.org/>`_ is only necessary if you installed the 
   code from git. Dependening on your system, you also may need to install
   the Python development package.

.. topic:: Solaris

   Although GammaLib builds on Solaris using the Sun compiler, there are
   problems with global symbols in shared libraries and exception catching,
   which prevents the FITS interface to work correctly. GammaLib has
   however been built and tested successfully using the GNU compiler, and
   this is the only build method that is currently supported. Problems have
   also been encountered when compiling cfitsio versions more recent than
   3.250. The problems have been reported to the cfitsio developer team,
   and are likely to be solved in the future. For the time being, it is
   recommended to use cfitsio version 3.250 on Solaris.

.. topic:: OpenSolaris

   On OpenSolaris, the same problems concerning the SunStudio compiler
   occur as for Solaris, and also here, the GNU compiler is the recommended
   tool to build GammaLib. Also here, cfitsio version 3.250 is the
   recommended library as more recent version feature relocation
   problems. GammaLib has been tested using gcc 4.3.2 on OpenSolaris
   2009.06. Make sure to create the following symbolic links if they do
   not yet exist on your system ::

      $ ln -s /usr/bin/gcc4.3.2 /usr/bin/gcc
      $ ln -s /usr/bin/g++4.3.2 /usr/bin/g++

   They avoid excess warnings during compilation.

.. _help:

Getting Help
------------

Any questions, bug reports, or suggested enhancements related to
GammaLib should be submitted via the
`issue tracker <https://cta-redmine.irap.omp.eu/projects/gammalib/issues/new>`_
or the
`mailing list <mailto:ctools@irap.omp.eu>`_.