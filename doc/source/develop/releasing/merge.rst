.. _dev_releasing_merge:

Merge release branch into master branch
=======================================

The master branch should always contain the last release of the GammaLib
package. You therefore have to merge now the release branch into the master
branch by typing:

.. code-block:: bash

   $ git checkout master
   $ git merge release
   $ git push origin
