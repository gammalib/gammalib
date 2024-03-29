.. _install_source:

Installation from source code
=============================

Compiling GammaLib from source is the recommended method to install the package
in case that you do not use Anaconda (if you use Anaconda, please follow
the instructions at :ref:`install_conda`, and better do not try to compile
the source code against your Anaconda Python installation). You also need to
compile from source if you want to use the development version of GammaLib.

.. note ::
   These installation instructions apply to gammalib-2.0.0 and
   later. gammalib-2.0.0 was built, installed and tested
   successfully on Debian, Ubuntu, Mandriva, OpenSUSE, Scientific Linux,
   CentOS, Mac OS X, FreeBSD, and OpenSolaris (using gcc), so
   hopefully it also installs on your distribution.
   **Windows is not supported.**
   In case you encounter problems, please post a report on the bug tracker or
   send an e-mail to the `mailing list <mailto:ctools@irap.omp.eu>`_.


Getting GammaLib
----------------

The latest version of the GammaLib source code, documentation, and
example programs are available on the World-Wide Web from
`http://cta.irap.omp.eu/gammalib/admin/download.html <http://cta.irap.omp.eu/gammalib/admin/download.html>`_.


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
`http://heasarc.gsfc.nasa.gov/fitsio <https://heasarc.gsfc.nasa.gov/fitsio/>`_
and installed from the source files (for details, please refer to the
installation instructions on the
`HEASARC <https://heasarc.gsfc.nasa.gov/fitsio/>`_ site). We recommend
installation of cfitsio in the default GammaLib install directory as a
shared library by typing

.. code-block:: bash

   $ ./configure --prefix=/usr/local/gamma
   $ make shared
   $ make install

GammaLib can also benefit from the presence of the readline library that
provides line-editing and history capabilities for text input (GammaLib
offers however full functionality without having readline
installed). readline (which depends on ncurses) is available on most
system as a prebuilt library, so also here we recommend to use your
package manager to install the libraries if they are not already there.
Also here, the readline (and ncurses) development packages are required,
so that the header files become available.

If readline and ncurses are not available as a prebuilt package, or if
you encounter some problems with the prebuilt packages, both libraries
can be downloaded from
`https://ftp.gnu.org/gnu/ <https://ftp.gnu.org/gnu/>`_
and installed from the source files. We recommend to put also those in the
default GammaLib install directory ``/usr/local/gamma``.

GammaLib includes a Python module that is compatible with Python 2 (version
2.3 or higher) and Python 3 (all versions). To generate the Python module,
the Python development package including the ``Python.h`` header file needs
to be installed on your system.


Building GammaLib
-----------------

To build GammaLib from the source files, download the latest release
tarball, save it in an appropriate location,
and type

.. code-block:: bash

   $ tar xvfz gammalib-x.y.z.tar.gz

where ``x.y.z`` is the version number of the library. Step in the created
directory using

.. code-block:: bash

   $ cd gammalib-x.y.z

and build GammaLib by typing

.. code-block:: bash

   $ ./configure
   $ make

at the operating system prompt. The configuration command customizes the
Makefiles for your particular system, the make command compiles the
source files and builds the C++ library and Python module. Type
``./configure`` and not simply ``configure`` to ensure that the configuration
script in the current directory is used and not some other system-wide
configuration script.

.. note::
   To accelerate compilation you may type

   .. code-block:: bash

      $ make -j8

   which will instruct the automake system to perform a parallelized compilation,
   where ``8`` specifies the maximum number of parallel threads (set this number
   to about twice the number of cores you have available).

By default, the install directory is set to ``/usr/local/gamma``. To change the
install directory, provide an optional ``--prefix`` argument, for example

.. code-block:: bash

   $ ./configure --prefix=/usr/local

If cfitsio and/or readline is not installed in a standard location for
libraries (such as ``/usr/lib`` or ``/usr/local/lib``), you can specify
the appropriate location(s) using the ``LDFLAGS`` (for the library
directory) and ``CPPFLAGS`` (for the include directory) options

.. code-block:: bash

   $ ./configure LDFLAGS='-L/opt/local/lib' CPPFLAGS='-I/opt/local/include'

You can find a full list of configuration options using

.. code-block:: bash

   $ ./configure --help


Testing GammaLib
----------------

Before you install GammaLib you should test the C++ library and Python 
module by typing

.. code-block:: bash

   $ make check

This will execute an extensive testing suite that should terminate with ::

   ============================================================================
   Testsuite summary for gammalib 2.0.0
   ============================================================================
   # TOTAL: 22
   # PASS:  22

or ::

   ===================
   All 22 tests passed
   ===================

dependent on your autotools version.

Eventually, loading the shared cfitsio and/or readline libraries may
fail during the test if the libraries are not located in standard
locations. In this case, add the library directories to the
``LD_LIBRARY_PATH`` environment variables (``DYLD_LIBRARY_PATH`` on Mac OS
X), e.g.

.. code-block:: bash

   $ export LD_LIBRARY_PATH=/opt/local/lib:$LD_LIBRARY_PATH


Installing GammaLib
-------------------

Now you are ready to install GammaLib by typing

.. code-block:: bash

   $ make install

at the operating system prompt. You may need to prepend a ``sudo`` in
case that you need administrator privileges to access the install
directory. If you do not have such privileges, chose an install directory
that you can access using the ``--prefix`` option.


.. _setup_env:

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

To build and install all documentation, type the following

.. code-block:: bash

   $ ./configure
   $ make doc
   $ make install

To build only user documentation, type

.. code-block:: bash

   $ ./configure
   $ make sphinx
   $ make install

and to build only Doxygen documentation, type

.. code-block:: bash

   $ ./configure
   $ make doxygen
   $ make install

The Doxygen documentation will be installed into
``$(prefix)/share/doc/gammalib/html/doxygen`` and can be browsed using a
regular web browser.

To check man support, type for example

.. code-block:: bash

   $ man GObservations

and you should see the documentation for the GObservations C++ class.
