.. _dev_profiling_install:

Installing valgrind
===================

Mac OS X
--------

The tools necessary to profile GammaLib on Mac OS X are best installed using
Homebrew. Below the sequence of things that need to be installed:

.. code-block:: bash

   $ brew install qcachegrind
   $ brew install valgrind
   $ brew install GraphViz

Now you should have ``qcachegrind`` and ``valgrind`` installed:

.. code-block:: bash

   $ which qcachegrind
   /usr/local/bin/qcachegrind
   $ which valgrind
   /usr/local/bin/valgrind
