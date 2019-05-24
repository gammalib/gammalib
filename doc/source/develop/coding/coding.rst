Coding conventions
==================

This section summarizes coding conventions for the GammaLib development.

C++ classes
-----------

General rules
^^^^^^^^^^^^^

Each class should be defined in a pair of individual files:

-  a header file that defines the class interface
   (with filename suffix ``.hpp``)

-  a source code file that implements the class
   (with filename suffix ``.cpp``)

In addition, a SWIG interface file should be provided for the Python
bindings (with filename suffix ``.i``).

In the header file, prefer forward declarations instead of ``#include``
directives to minimise the dependencies between the files. Specify only
the ``#include`` directives that are absolutely needed to compile the
code.

In the source file, put all ``#include`` directives that are necessary for
compilation of the specific file.

The C++ style header files should be used instead of the C style header
files to ensure maximum portability. The following table provides the
correspondence between C++ and C header files for headers commonly used
in GammaLib:

============= ============== =================
C++           C              Function examples
============= ============== =================
``<cctype>``                 
``<cmath>``   ``<math.h>``   ``std::abs``, ``std::cos``
``<cfloat>``                 
``<cstdio>``  ``<stdio.h>``  ``std::fopen``, ``std::fgets``, ``std::fclose``, ``std::fprintf``, ``std::sprintf``
``<cstdlib>``                
``<cstring>`` ``<string.h>`` ``std::strlen``
============= ============== =================

Note that functions and types should be prefixed by ``std::``. For example,
``cos`` becomes ``std::cos``, ``time_t`` becomes ``std::time_t``, etc.
One significant change between C and C++ is that ``fabs`` becomes ``std::abs``
since the C style ``abs`` function only applies to integers. Here, the
``std::`` prefix is crucial to distinguish the C++ function (which is also
defined for doubles) from the C function.

.. _sec_header:

Header file structure
^^^^^^^^^^^^^^^^^^^^^

The header file defines the interface of the class. Here is an example
of a header file.

.. code-block:: cpp

    /***************************************************************************
     *                        GClass.hpp  -  My nice class                     *
     * ----------------------------------------------------------------------- *
     *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
     * ----------------------------------------------------------------------- *
     *                                                                         *
     *  This program is free software: you can redistribute it and/or modify   *
     *  it under the terms of the GNU General Public License as published by   *
     *  the Free Software Foundation, either version 3 of the License, or      *
     *  (at your option) any later version.                                    *
     *                                                                         *
     *  This program is distributed in the hope that it will be useful,        *
     *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
     *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
     *  GNU General Public License for more details.                           *
     *                                                                         *
     *  You should have received a copy of the GNU General Public License      *
     *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
     *                                                                         *
     ***************************************************************************/
    /**
     * @file GClass.hpp
     * @brief Definition of my nice class interface
     * @author Juergen Knoedlseder
     */

    #ifndef GCLASS_HPP
    #define GCLASS_HPP

    /* __ Includes ___________________________________________________________ */
    #include <string>
    #include "GBase.hpp"


    /***********************************************************************//**
     * @class GClass
     *
     * @brief Illustration of a GammaLib class
     *
     * My nice class illustrates how a GammaLib class should be defined.
     ***************************************************************************/
    class GClass : public GBase {

    public:
        // Constructors and destructors
        GClass(void);
        GClass(const GClass& c);
        virtual ~GClass(void);
     
        // Operators
        GClass& operator=(const GClass& c);

        // Methods
        void        clear(void);
        GClass*     clone(void) const;
        std::string print(const GChatter& chatter = NORMAL) const;
      
    protected:
        // Protected methods
        void init_members(void);
        void copy_members(const GClass& c);
        void free_members(void);

        // Protected data members
        std::string     m_name;          //!< Name
    };

    #endif /* GCLASS_HPP */

The header file starts with a comment containing the file name and class
purpose, the copyright information and the license text. The years in
the copyright information should cover the years over which the file has
been modified, the author is the person who initially created the file.

Following the header comment is a comment that provides file information
to the Doxygen documentation system.

The subsequent

.. code-block:: cpp

        #ifndef GCLASS_HPP
        #define GCLASS_HPP

declarations together with the

.. code-block:: cpp

        #endif /* GCLASS_HPP */

declaration at the end protect the file from multiple inclusions of the
header. This is a crucial feature needed for proper compilation of the
code.

Now all header files are included. Standard header files are included
using the ``< >`` brackets, GammaLib header files are included using ``" "``. A
80 character long separator precedes the header inclusion. Further 80
character long separators may be added for additional sections, such as
constants, type definitions, forward declarations, etc. Use one
separator to precede each additional section.

The class definition is preceded by a comment block that will be used by
the Doxygen system to extract the class definition. Provide here the
class name, a brief one line description of the class, and an extended
detailed description of the class purpose.

The class definition is structured in several sections:

- Definition of public constructors
- Definition of public operators
- Definition of public methods
- Definition of protected methods
- Definition of protected members

Note that most classes will derive from the abstract interface class
:doxy:`GBase` which imposes the implementation of the :doxy:`GBase::clear`, the
:doxy:`GBase::clone` and the :doxy:`GBase::print` methods.

The definition of pure virtual methods should be done in a section that
is separate from the methods that are implemented.

Here an illustration of the expected structure, based on the
:doxy:`GObservation` class:

.. code-block:: cpp

    class GObservation : public GBase {

    public:
        // Constructors and destructors
        GObservation(void);
        GObservation(const GObservation& obs);
        virtual ~GObservation(void);

        // Operators
        virtual GObservation& operator=(const GObservation& obs);

        // Pure virtual methods
        virtual void          clear(void) = 0;
        virtual GObservation* clone(void) const = 0;
        virtual std::string   print(const GChatter& chatter = NORMAL) const = 0;

        // Virtual methods
        virtual double        model(const GModels& models, const GEvent& event, GVector* gradient = NULL) const;
        virtual double        npred(const GModels& models, GVector* gradient = NULL) const;

        // Implemented methods
        void                  name(const std::string& name);
        void                  id(const std::string& id);

    protected:
        // Protected methods
        void init_members(void);
        void copy_members(const GObservation& obs);
        void free_members(void);

        // Protected data area
        std::string m_name;         //!< Name of observation
        std::string m_id;           //!< Observation identifier
        std::string m_statistics;   //!< Optimizer statistics (default=poisson)
        GEvents*    m_events;       //!< Pointer to event container
    };

.. _sec_sourcecode:

Source code file structure
^^^^^^^^^^^^^^^^^^^^^^^^^^

The source code file implements the code of the class. Here is an
example of the start of a source code file.

.. code-block:: cpp

    /***************************************************************************
     *                        GClass.cpp  -  My nice class                     *
     * ----------------------------------------------------------------------- *
     *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
     * ----------------------------------------------------------------------- *
     *                                                                         *
     *  This program is free software: you can redistribute it and/or modify   *
     *  it under the terms of the GNU General Public License as published by   *
     *  the Free Software Foundation, either version 3 of the License, or      *
     *  (at your option) any later version.                                    *
     *                                                                         *
     *  This program is distributed in the hope that it will be useful,        *
     *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
     *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
     *  GNU General Public License for more details.                           *
     *                                                                         *
     *  You should have received a copy of the GNU General Public License      *
     *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
     *                                                                         *
     ***************************************************************************/
    /**
     * @file GClass.cpp
     * @brief Implementation of my nice class
     * @author Juergen Knoedlseder
     */

    /* __ Includes ___________________________________________________________ */
    #ifdef HAVE_CONFIG_H
    #include <config.h>
    #endif
    #include "GClass.hpp"
    #include "GTools.hpp"

    /* __ Method name definitions ____________________________________________ */
    #define G_CLEAR                                             "GClass::clear()"
    #define G_CLONE                                       "GClass::clone() const"
    #define G_PRINT                              "GClass::print(GChatter&) const"

    /* __ Compile options ____________________________________________________ */
    #define G_USE_MY_OPTION

    /* __ Debug options ______________________________________________________ */
    #define G_DEBUG_PRINT

    /* __ Constants __________________________________________________________ */
    const double pi = 3.14;

The include section starts with a conditional include of the code
configuration header file (see :ref:`sec_configure`). This makes
GammaLib compile options available to the source code.

The include section is followed by the declaration of method names.
These method names will be used in exceptions 
(see :ref:`sec_exceptions`).
Define the method names at the top of the file eases
the maintainability of the code, as changes in method names or
interfaces need only to be implemented in a single place. Method names
need only be defined for methods throwing exceptions.

Compile options are used to control which parts of the code should be
compiled. Such options may be used, for example, to compare different
algorithms or computation methods. They can also be used during
development, allowing an easy switch between the new and the old code
for comparison.

Debug options are compile options that are used to add additional code
for debugging. Often, these are print statements that allow to trace the
execution of the code. For code checked into the repository, all debug
options should be commented out.

.. _sec_python:

Python interface for C++ classes
--------------------------------

The Python interface for C++ classes is defined by a so-called `SWIG <http://www.swig.org/>`_
interface file. SWIG uses these interface files to build Python wrapper
files, which are C files that define the interface between GammaLib and
Python. The structure of the SWIG interface file follows closely that of
the header file, with a few exceptions. Here an example:

.. code-block:: cpp

    /***************************************************************************
     *                         GClass.i  -  My nice class                      *
     * ----------------------------------------------------------------------- *
     *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
     * ----------------------------------------------------------------------- *
     *                                                                         *
     *  This program is free software: you can redistribute it and/or modify   *
     *  it under the terms of the GNU General Public License as published by   *
     *  the Free Software Foundation, either version 3 of the License, or      *
     *  (at your option) any later version.                                    *
     *                                                                         *
     *  This program is distributed in the hope that it will be useful,        *
     *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
     *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
     *  GNU General Public License for more details.                           *
     *                                                                         *
     *  You should have received a copy of the GNU General Public License      *
     *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
     *                                                                         *
     ***************************************************************************/
    /**
     * @file GClass.i
     * @brief Python interface of my nice class
     * @author Juergen Knoedlseder
     */
    %{
    /* Put headers and other declarations here that are needed for compilation */
    #include "GClass.hpp"
    %}

    /***********************************************************************//**
     * @class GClass
     *
     * @brief Illustration of a GammaLib class
     *
     * My nice class illustrates how a GammaLib class should be defined.
     ***************************************************************************/
    class GClass : public GBase {
    public:
        // Constructors and destructors
        GClass(void);
        GClass(const GClass& c);
        virtual ~GClass(void);

        // Methods
        void        clear(void);
        GClass*     clone(void) const;
    };


    /***********************************************************************//**
     * @brief GClass class extension
     ***************************************************************************/
    %extend GClass {
        GClass copy() {
            return (*self);
        }
    };

The code starts with a section that is enclosed in ``%{ %}`` brackets. In
this section, all header files are specified that are needed to compile
the SWIG wrapper file.

Then follows the class definition, with the following differences with
respect to the definition in the header file:

-  it does not include the assignment operator
-  it does not include any access operator (these have to be implemented
   specifically, see below)
-  it does not include the ``print()`` method (see below)
-  it does not include protected or private members

Finally, there is a section with extension to the C++ class. Here,
methods are implemented that do not exist in the actual C++ class, but
that will exist in the Python interface.

In case that an access operator needs to be implemented, the
``__getitem__()`` and ``__setitem__()`` methods need to be added to the
class extensions. Here an example:

.. code-block:: cpp

  /***********************************************************************//**
   * @brief GObservations class extension
   ***************************************************************************/
  %extend GObservations {
      GObservation* __getitem__(const int& index) {
          if (index >= 0 && index < self->size()) {
              return (*self)[index];
          }
          else {
              throw GException::out_of_range("__getitem__(int)", index, self->size());
          }
      }
      void __setitem__(const int& index, const GObservation& val) {
          if (index >= 0 && index < self->size()) {
              self->set(index, val);
              return;
          }
          else {
              throw GException::out_of_range("__setitem__(int)", index, self->size());
          }
      }
  };

Note that the access operators perform explicit range checking because the
``[]`` operators used do not any range checking. In addition, in this was
the Python operator names ``__getitem__()`` and ``__setitem__()`` can be
specified in the exception.
