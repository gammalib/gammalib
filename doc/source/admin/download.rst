.. _sec_download:

Download
========

GammaLib can be obtained in form of releases or directly from the git 
development repository. Prefer a release if you intend using GammaLib
for production (and publications). Clone the code from git if you need
the most recent code that implements new features and corrects known
bugs.


Releases
--------

The latest GammaLib release is
`gammalib-1.4.3 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.4.3.tar.gz>`_
(31 October 2017).

Below a table of older GammaLib releases. Please read the :ref:`sec_release` to
learn more about new features and corrected bugs in a given release.

To download an older release, click on the corresponding release number:
`1.4.2 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.4.2.tar.gz>`_
`1.4.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.4.1.tar.gz>`_
`1.4.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.4.0.tar.gz>`_
`1.3.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.3.1.tar.gz>`_
`1.3.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.3.0.tar.gz>`_
`1.2.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.2.0.tar.gz>`_
`1.1.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.1.0.tar.gz>`_
`1.0.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.0.1.tar.gz>`_
`1.0.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.0.0.tar.gz>`_
`0.11.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-0.11.0.tar.gz>`_
`0.10.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-0.10.0.tar.gz>`_
`0.9.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-09-01.tar.gz>`_
`0.9.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-09-00.tar.gz>`_
`0.8.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-08-01.tar.gz>`_
`0.8.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-08-00.tar.gz>`_
`0.7.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-07-00.tar.gz>`_
`0.6.2 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-06-02.tar.gz>`_
`0.6.1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-06-01.tar.gz>`_
`0.5.0 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-05-00.tar.gz>`_
`0.4.2 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-02.tar.gz>`_
`0.4.11 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-11.tar.gz>`_
`0.4.10 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-10.tar.gz>`_
`0.4.9 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-09.tar.gz>`_
`0.4.7 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-00-04-07.tar.gz>`_


Development release
-------------------

The current GammaLib development release is
`gammalib-1.5.0.dev1 <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.5.0.dev1.tar.gz>`_.
This release reflects the status of the current ``devel`` branch of
the GammaLib git repository.


Git repository
--------------

To clone the gammalib source code, type

.. code-block:: bash

   $ git clone https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git
  
This will create a gammalib directory under the current working directory
that will contain the gammalib source code. In case that the cloning does
not work you may try adding

.. code-block:: bash

   $ export GIT_SSL_NO_VERIFY=true

or

.. code-block:: bash

   $ git config --global http.sslverify "false"

before retrieving the code.
Before you will be able to compile the code you need to generate the
configuration file using the ``autogen.sh`` script.
Also make sure that you're actually on the devel branch of the git
repository. GammaLib can be compiled and configured using
the following command sequence (the code will be installed into the 
``/usr/local/gamma`` directory):

.. code-block:: bash

   $ cd gammalib
   $ git checkout devel
   $ ./autogen.sh
   $ ./configure
   $ make
   $ make check
   $ sudo make install
   $ export GAMMALIB=/usr/local/gamma
   $ source $GAMMALIB/bin/gammalib-init.sh

Please read the :ref:`sec_installation` section if you need more information on
how to install GammaLib.

.. note::

  You need `swig <http://www.swig.org/>`_ on your system to build the
  Python wrappers when you get the code from Git. Python wrappers are
  not stored in the Git repository but are built using
  `swig <http://www.swig.org/>`_ from interface files located in the
  pyext folder. However, you do not need `swig <http://www.swig.org/>`_
  when fetching a release as the Python wrappers are bundled with the
  release tarballs.
