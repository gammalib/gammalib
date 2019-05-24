.. _install_devel:

Installing the development version
==================================

The current development version of the code can be downloaded as source
tarballs or directly from the `Git <https://git-scm.com/>`_ repository.

The current GammaLib development release is ``gammalib-1.7.0.dev``.


Source tarball
--------------

Download the source tarball from the following link

- `GammaLib <http://cta.irap.omp.eu/ctools/releases/gammalib/gammalib-1.7.0.dev.tar.gz>`_

and follow the instructions on :ref:`install_source`.


Git repository
--------------

To clone the gammalib source code, type

.. code-block:: bash

   $ git clone https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git

This will create a directory named gammalib under the current
working directory that will contain the gammalib source code.
In case that the cloning does not work you may try adding

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

Please read the :ref:`install_source` section if you need more information
on how to install GammaLib.

.. note::
   You need `swig <http://www.swig.org/>`_ on your system to build the
   Python wrappers when you get the code from Git. Python wrappers are
   not stored in the Git repository but are built using
   `swig <http://www.swig.org/>`_ from interface files located in the
   pyext folder. However, you do not need `swig <http://www.swig.org/>`_
   when fetching a release as the Python wrappers are bundled with the
   release tarballs.


Making a conda package
----------------------

To create a conda package for the development version, clone the gammalib
source code using

.. code-block:: bash

   $ git clone https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git

and type the following to create the conda packages (make sure that anaconda is
included in your ``$PATH`` environment):

.. code-block:: bash

   $ cd gammalib
   $ git checkout devel
   $ ./autogen.sh
   $ ./configure
   $ conda-build dev/conda.recipe

Once this is done, you can create a conda environment using the development
version as follows:

.. code-block:: bash

   $ conda create -n gammalib-devel python=3.6
   $ source activate gammalib-devel
   $ conda install --use-local gammalib=1.7.0.dev

Note that you can chose between Python 2.7, 3.5, 3.6 and 3.7 for your
conda environment.

