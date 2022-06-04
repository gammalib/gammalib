.. _install_conda:

Installing via Anaconda
=======================

`Anaconda <https://www.anaconda.com/download/>`_ (you can also use
`Miniconda <https://conda.io/miniconda.html>`_ for a smaller footprint) is a
scientific Python installation shipping with essentially all needed packages.
Install it according to the instructions on their homepage. You can use any
Anaconda Python version. An Anaconda Python installation is completely separate
from any existing system wide or user space Python installation, so Anaconda
can be tested without the fear of breaking an existing installation.

**Anaconda GammaLib packages exist for Mac OS X and Linux distributions. Windows
is not supported.**


Installing
----------

- `Install Anaconda following the instructions on their site <https://www.anaconda.com/download/>`_

- Add the ``conda-forge`` and ``cta-observatory`` channels to your Anaconda
  configuration

  .. code-block:: bash

     $ conda config --append channels conda-forge
     $ conda config --append channels cta-observatory

- We strongly recommend to work with separate Anaconda environments, and
  especially not use the special root environment (that is used for all conda
  commands and environment manipulations, package installations etc.) for
  anything besides updating the conda package itself. For example create
  a ``myenv`` environment as follows:

  .. code-block:: bash

     $ conda create -n myenv python=3.8  # or one of the following Python versions: 2.7, 3.5, 3.6, 3.7, 3.8, 3.9, 3.10
     $ conda activate myenv
     (myenv) $

- Install pre-compiled GammaLib conda package from Anaconda cloud with:

  .. code-block:: bash

     (myenv) $ conda install gammalib

  This will also install the required dependencies, and in particular cfitsio.


Testing
-------

Type the following to test the GammaLib package

.. code-block:: bash

   (myenv) $ python -c 'import gammalib; gammalib.test()'


Updating
--------

Type the following to update the GammaLib package

.. code-block:: bash

   (myenv) $ conda update gammalib

or (if you did not add the ``cta-observatory`` channel to your default
channels):

.. code-block:: bash

   (myenv) $ conda update -c cta-observatory gammalib
