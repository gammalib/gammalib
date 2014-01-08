What is in this section?
------------------------

This section provides an overview over all GammaLib modules and their C++
classes, with particular emphasis on the relation between the classes
and their basic functionalities. It describes the purpose of all C++
classes and their primary usage, as well as their underlying
arithmetics. However, we do not provide a detailed description of the
interface and the inner workings of each C++ class. This information is
provided in the reference documentation, which can be found online at
http://gammalib.sourceforge.net/doxygen/, or which can be installed
locally on your machine (see section :ref:`sec_doxygen`).

Each GammaLib module is presented in a dedicated section, following the overview
shown in :ref:`fig_structure` from the top-left to the bottom right.
Instrument specific modules are described in a dedicated chapter (see 
:ref:`sec_inst`). All C++ classes of a module and their relations are
illustrated using a UML diagram.

To explain how to read such a diagram, we show an example for five
fictive classes in :ref:`fig_uml_template`. Our example shows a container
class that contains an arbitrary number of elements which are realized
by an abstract base class. Names of abstract base classes are indicated
in *italic* to highlight the fact that such classes can not be
instantiated. The possible number of elements that may be held by the
container (in this case any number) is indicated by the cardinality
0..\* situated next to the abstract base class. Our example shows also
two derived classes that inherit from the abstract base class. The
second derived class is associated with a single element of some other
class, indicated by the cardinality 1 next to the class. This other
class is not part of the actual module, and is thus shown in grey with a
dotted boundary.

.. _fig_uml_template:

.. figure:: uml_template.png
   :width: 100%

   UML usage
