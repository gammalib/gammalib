.. _dev_howto_module_codegen:

Add module using the code generator
===================================

  .. admonition:: What you will learn

     You will learn how you **add template code for an instrument module
     to GammaLib** using the code generator.

You should use the code generator that is included in GammaLib to add
a template for a new instrument module to the code base. Let's assume
that you want to create a module for the
`DAMPE <http://dpnc.unige.ch/dampe/>`_
satellite. Before starting, make sure that you get rid of any existing
compiled code:

.. code-block:: bash

   $ make clean

Now you should create a new branch that will contain the new module:

.. code-block:: bash

   $ git checkout -b add-dampe-module

To this new branch you should add now the new module using the code
generator. The code generator in invoked as follows from the root
of the GammaLib directory:

.. code-block:: bash

   $ dev/codegen.py
   GammaLib code generator
   =======================

   [1] Add instrument module
   [q] Quit
   Enter your choice:

..

  .. note::

     Note that you need gammalib installed and configured to use the code
     generator. You can only use the code generator from the root of the
     GammaLib code directory.

Now select ``1`` to add a new instrument module. This will ask you a few
questions:

.. code-block:: bash

   Add instrument module
   ---------------------
   Please enter a 3 digit-name for the module: dam
   Please enter the instrument name (e.g. "Fermi/LAT"): DAMPE
   Please enter your name (e.g. "Joe Public"): Joe Public

   All right. Have now:
   Module name .....: "dam"
   Instrument name .: "DAMPE"
   Your name .......: "Joe Public"
   Is this correct? (y/n):

Once you are satisfied with the answers, press ``y``. Now you have to answer
two more questions about the event format that should be supported. In most
cases your data will be in form of event list. If this is your case, answer
``y`` to the first question. In many cases you want to be able to bin the
data into a counts cube. If this is your case, answer also ``y`` to the second
question:

.. code-block:: bash

  All right. You want a new "dam" instrument module.
  Do you want event list support? (y/n): y
  Do you want binned event data support? (y/n): y
..

  .. note::

     The code generator also allows to add <missing `dev_howto_module_missing`>
     classes later, so if you are not yet sure you may also answer ``n``
     for now.

Now you are done. Enter ``q`` to exit the code generator.

You put the new code under
`Git <https://git-scm.com/>`_
control by typing

.. code-block:: bash

   $ git add .
   $ git commit -m "Add code templates for DAMPE instrument module"
   [add-dampe-module c573eff] Add code templates for DAMPE instrument module
    39 files changed, 5288 insertions(+), 13 deletions(-)
    create mode 100644 inst/dam/Makefile.am
    create mode 100644 inst/dam/README.md
    create mode 100644 inst/dam/include/GDAMEventAtom.hpp
    create mode 100644 inst/dam/include/GDAMEventBin.hpp
    create mode 100644 inst/dam/include/GDAMEventCube.hpp
    create mode 100644 inst/dam/include/GDAMEventList.hpp
    create mode 100644 inst/dam/include/GDAMInstDir.hpp
    create mode 100644 inst/dam/include/GDAMLib.hpp
    create mode 100644 inst/dam/include/GDAMObservation.hpp
    create mode 100644 inst/dam/include/GDAMResponse.hpp
    create mode 100644 inst/dam/include/GDAMRoi.hpp
    create mode 100644 inst/dam/pyext/GDAMEventAtom.i
    create mode 100644 inst/dam/pyext/GDAMEventBin.i
    create mode 100644 inst/dam/pyext/GDAMEventCube.i
    create mode 100644 inst/dam/pyext/GDAMEventList.i
    create mode 100644 inst/dam/pyext/GDAMInstDir.i
    create mode 100644 inst/dam/pyext/GDAMObservation.i
    create mode 100644 inst/dam/pyext/GDAMResponse.i
    create mode 100644 inst/dam/pyext/GDAMRoi.i
    create mode 100644 inst/dam/pyext/dam.i
    create mode 100644 inst/dam/src/GDAMEventAtom.cpp
    create mode 100644 inst/dam/src/GDAMEventBin.cpp
    create mode 100644 inst/dam/src/GDAMEventCube.cpp
    create mode 100644 inst/dam/src/GDAMEventList.cpp
    create mode 100644 inst/dam/src/GDAMInstDir.cpp
    create mode 100644 inst/dam/src/GDAMObservation.cpp
    create mode 100644 inst/dam/src/GDAMResponse.cpp
    create mode 100644 inst/dam/src/GDAMRoi.cpp
    create mode 100644 inst/dam/test/Makefile.am
    create mode 100644 inst/dam/test/test_DAM.cpp
    create mode 100644 inst/dam/test/test_DAM.hpp
    create mode 100644 inst/dam/test/test_DAM.py
  $ git push origin add-dampe-module

Finally, configure and compile your code, including the new module, by
typing

.. code-block:: bash

   $ autoconf
   $ automake
   $ ./configure
   $ make
   $ make check
