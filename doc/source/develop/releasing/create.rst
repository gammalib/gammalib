.. _dev_releasing_create:

Create the release
==================

The GammaLib release is created by invoking the release manager:

.. code-block:: bash

   $ dev/release.py
   GammaLib release manager
   ========================

   [1] Make a new release
   [2] Set the package version in current branch
   [3] Set the libtool version in current branch
   [4] Commit changes in current branch
   [5] Create tarball from current branch
   [6] Check tarball
   [q] Quit
   Enter your choice:

To create a new release type ``1``. You will then be guided through a number
of questions. In the example below we show how you go from the development
version ``1.2.0.dev`` to the release ``1.2.0``. The first step is to create
a ``release`` branch.

.. code-block:: bash

   Make a new release
   ------------------
   Step 1: Do you want to create a 'release' branch? (y/n): y
   Created new branch 'release'
   Switched to a new branch 'release'

Then you have to set the version number of the release.

.. code-block:: bash

   Step 2: Current GammaLib version is '1.2.0.dev'. Do you want to change the package version? (y/n): y
   Current GammaLib version is '1.2.0.dev'. Please enter new GammaLib version: 1.2.0
   Change version to '1.2.0'? (y/n): y
   GammaLib version changed to '1.2.0'

After this the libtool version number will be set. **The libtool version number
will in general be different from the release version number.** The new libtool
version number is determined from the source code and interface changes with
respect to the last release. You therefore have to answer a couple of
questions about the changes you did.

.. code-block:: bash

   Step 3: Current Libtool version is '2:0:0'. Do you want to change the libtool version? (y/n): y
   Current libtool version is '2:0:0'. Do you want to change the version? (y/n): y
    Has the source code changed since last release? (y/n): y
    Have interfaces been added since last release? (y/n): y
    Have interfaces been removed since last release? (y/n): y
    Have interfaces changed since last release? (y/n): y
   Libtool version changed to '3:0:0'

Then you have to commit all changes and push them into the git repository.

.. code-block:: bash

   Step 4: Commit changes? (y/n): y
   The following files have been changed:
    M README.md
    M configure.ac
    M doc/Doxyfile
    M doc/source/conf.py
    M gammalib.pc.in
    M sonar-project.properties
   Commit all changes? (y/n): y
   [release addb070] GammaLib package version set to '1.2.0' and libtool version set to '3:0:0'
    6 files changed, 8 insertions(+), 8 deletions(-)

   Step 5: Push changes? (y/n): y
   Counting objects: 10, done.
   Delta compression using up to 8 threads.
   Compressing objects: 100% (10/10), done.
   Writing objects: 100% (10/10), 880 bytes | 0 bytes/s, done.
   Total 10 (delta 8), reused 0 (delta 0)
   remote: To https://github.com/gammalib/gammalib.git
   remote:  * [new branch]      release -> release
   To https://cta-gitlab.irap.omp.eu/gammalib/gammalib.git
    * [new branch]      release -> release

You then should build and check the source tarball to verify that everything
went fine. Note that this step is not formally needed for a release since
the software release will start from the code in the ``release`` branch.

.. code-block:: bash

   Step 6: Build tarball? (y/n): y
   Log actions in logfile? (y/n): y
   Configure package
   Package configuration successful
   Create tarball
   Tarball creation successful

   Step 7: Check tarball? (y/n): y
   Log check in logfile? (y/n): y
   Configure package
   Package configuration successful
   Check tarball
   Tarball checking successful

Now you are done and can quite the release manager.

.. code-block:: bash

   [1] Make a new release
   [2] Set the package version in current branch
   [3] Set the libtool version in current branch
   [4] Commit changes in current branch
   [5] Create tarball from current branch
   [6] Check tarball
   [q] Quit
   Enter your choice: q


