Miscellaneous
=============

GammaLib Version Numbering
--------------------------

GammaLib applies a three-number version numbering scheme::

     major.minor.patch

A ``major`` revision of 0 indicates that the GammaLib design is not yet
frozen. At this level, external interfaces of GammaLib may change
without notification, and no interface control system is implemented.
The ``minor`` revision tag will be incremented for each new release,
signaling that new features have become available. The ``patch`` tag will be
incremented after correcting bugs that were reported on releases.

Once the GammaLib design is frozen, the ``major`` revision number will be
incremented to 1. From this moment on, external interfaces of GammaLib
will be under configuration control. If existing external interfaces
will be modified, the ``major`` revision number will be incremented. At the
same time, the libtool version number of the GammaLib will also be
incremented. The ``minor`` revision number will be incremented if
modifications and extensions are backward compatible. As before, the
``patch`` number will be incremented after correcting bugs that were
reported on releases.
