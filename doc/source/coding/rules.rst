General coding rules
====================

This section provides lists of general coding rules for the GammaLib
development.

C++ rules
---------

All GammaLib C++ code must be compatible with ISO standard C++98.
Do not use C++11 features.
(see `here <http://en.wikipedia.org/wiki/C%2B%2B#Standardization>`__
if you don't know what C++98 and C++11 is and note that C++03 is
just a minor correction of C++98).

The reason for this restriction is that C++ compilers / standard libraries
such as `gcc <http://gcc.gnu.org>`_ / `libstdc++ <http://gcc.gnu.org/libstdc++/>`_
or `clang <http://clang.llvm.org>`_ / `libc++ <http://libcxx.llvm.org>`_
only finished implementing C++11 in 2013 and using C++11 features in GammaLib would
make it impossible to use GammaLib on user and server machines with an old compiler.


Code format
^^^^^^^^^^^

-  Blocks are indented by 4 characters.

-  Do not use tabs (code formatting should be independent of editor
   configurations).

-  Do not exceed a line length of 80 characters (a few more characters
   are acceptable in exceptional cases).

-  Put a blank line at the end of each file (this is required by some
   compilers).

-  Each function starts with a curly bracket in the line following the
   function name. The return type is in the same line as the function
   name. Example:

   .. code-block:: cpp

       int function(void)
       {
           int i = 0;
           ...
           return i;
       }

-  Code blocks should be encompassed in curly brackets, even if the block
   consists only of a single line.

-  The opening curly bracket of a block starts in the same line as the
   related statement. Example:

   .. code-block:: cpp

       for (int i = 0; i < 10; ++i) {
           sum += i;
       }

-  Separate code elements by spaces; for example, do not write ``int i=0;``
   but write ``int i = 0;``.

-  Align successive similar lines on common elements. Here a common example
   illustrating the alignment on the ``=`` symbol:

   .. code-block:: cpp

       m_max    = par.m_max;
       m_prompt = par.m_prompt;
       sum     += par.m_sum;

   And here an example illustrating the alignment in a class definition on the
   member function name:

   .. code-block:: cpp

       void        log10GeV(const double& eng);
       void        log10TeV(const double& eng);
       std::string print(void) const;

Code semantics
^^^^^^^^^^^^^^

-  Each function and/or method terminates with a return statement.

-  Each function and/or method has only a single exit point (i.e. a
   single return statement).

-  Use ``explicit`` for constructors with single arguments to prevent unintended
   type conversions. The only exception to this rule is the copy
   constructor or type conversion constructors.

-  Specify void for function definitions without arguments, i.e. use
   ``void function(void)`` instead of ``void function()`` in the
   function declaration.

-  Use pre-incrementation in loops (pre-incrementation is faster than
   post-incrementation). Example:

   .. code-block:: cpp

       for (int i = 0; i < 10; ++i) {
           sum += i;
       }

-  Where possible (and appropriate), use ``std::vector`` containers instead
   of allocating memory. In other words: avoid direct memory allocation with ``new``.

-  Use the ``std::`` namespace prefix where possible; write for example

   .. code-block:: cpp

       std::sin(angle);
       std::cos(angle);

   You may not believe it, but droping the ``std::`` may on some systems
   lead to considerably slower code for trigonometric functions!

-  Provide comments, comments, comments!!!

Language features
^^^^^^^^^^^^^^^^^

-  Do not use macros.

-  Do not use ``#define`` directives for the declaration of constants. Use
   ``const`` instead.

-  Do not use ``std::strncpy``, ``std::memcpy`` or similar as these functions
   are corrupted on some systems.

-  If possible, pass arguments by reference.

-  Output arguments should be passed as pointers.

-  Use C++ (``std::string``) instead of C-style (``char*``) strings.

-  Use C++ casts instead of C-style casts.

-  Avoid using templates.


Python rules
------------

All code must be compatible with Python 2.6, 2.7 as well as 3.2 or later.
As of 2014 most users still have Python 2, even though Python 3 has been
out for five years, so using Python 3 features is out of the question.
Also GammaLib is mostly a C++ library with only some high-level analysis
scripts (e.g. running a complete simulation or analysis or plotting some results)
that are simple Python code and wouldn't profit much from using new Python 3 features.

Code format
^^^^^^^^^^^

-  Python code should follow the official 
   `PEP8 Python style <http://www.python.org/dev/peps/pep-0008/>`_:
   It says (among many other things) that you should indent with four
   spaces, not tabs. Following PEP8 is simple, because there's a
   `pep8 tool <https://github.com/jcrocholl/pep8>`_
   that you should run on your Python code before committing, e.g. like so:

   .. code-block:: bash

       $ pep8 test/test_python.py
       ...
       test/test_python.py:156:1: W191 indentation contains tabs
       test/test_python.py:156:1: W391 blank line at end of file
       test/test_python.py:156:1: W293 blank line contains whitespace
       $

   If you want you can even use the
   `autopep8 tool <https://github.com/hhatto/autopep8>`_
   which can automatically fix the
   formatting for almost all cases. Run ``pep8 -h`` and ``autopep8 -h`` to see
   the options you can use. PEP8 compliance is automatically checked by
   the continuous integration system.

Code semantics
^^^^^^^^^^^^^^

- When writing scripts, the following import conventions should be used::

     >>> import gammalib
     >>> energy = gammalib.GEnergy(10, 'TeV')
     >>> import numpy as np
     >>> import matplotlib.pyplot as plt  

  In interactive sessions (`IPython <http://ipython.org/>`_ is recommended) to save
  typing you can use ``import *`` to import everything into your local namespace::
  
     >>> from gammalib import *
     >>> energy = GEnergy(10, 'TeV')
     >>> from numpy import *
     >>> from matplotlib.pyplot import *
  
  Using ``import *`` in Python scripts is not recommended. Some reason are given in
  `this section <http://docs.python.org/3/howto/doanddont.html#from-module-import>`__
  in the Python docs.
  
- Most GammaLib code can be written equivalently in C++ and / or Python,
  however some things have been made `"Pythonic" <https://www.google.de/search?q=pythonic>`__
  in the SWIG Python wrapper.
  GammaLib objects can be printed by passing them to the Python ``print`` function.
  GammaLib container classes feature a Pythonic list-like interface, i.e. they
  can be iterated over and elements can be accessed with ``[]``. Examples:
  
  .. code-block:: python

     # Create a container with two elements
     >>> from gammalib import GEnergy, GEnergies
     >>> energies = GEnergies()
     >>> energies.append(GEnergy(1, 'TeV')) 
     >>> energies.append(GEnergy(3, 'TeV'))
     # Print GammaLib object Python-style.
     # C++ style would be calling a member function: energies[1].print()
     >>> print(energies[1])
     2 TeV
     # Use the container like a Python list 
     >>> len(energies)
     2
     >>> for energy in energies:
             print(energy)
     1000 GeV
     2 TeV

Language features
^^^^^^^^^^^^^^^^^

You can find lots of information on how to write Python code that works with
Python 2 and 3 on the web, e.g.
`here <http://docs.python.org/dev/howto/pyporting.html#use-same-source>`__
and `here <http://python3porting.com/noconv.html>`__.
If you don't have a Python 2 or Python 3 interpreter on your development
machine, don't worry. This will be checked during code review or during
continuous integration testing.

The most common Python 2 / 3 issues are ``print`` and integer division.

-  You should use the ``print()`` function, because the ``print`` statement
   was removed in Python 3::
   
      >>> print('hello')  # print function. Good. Works in Python 2 and 3.
      hello
      >>> print 'hello'   # print statement. Bad. Only works in Python 2.
      SyntaxError: invalid syntax

-  Integer division behaves differently in Python 2 and 3.
   
      >>> 5 / 2  # Python 3: integer division results in float
      2.5
      >>> 5 / 2  # Python 2: integer division results in int
      2

   If you do need integer division in your Python code, you should
   use one of these explicit forms that work in Python 2 and 3::
   
      >>> int(5 / 2)  # Python 2 and 3: truncate int part 
      2
      >>> 5 // 2  # Python 2 and 3: integer division
      2


A good way to enforce Python 3 compatible behaviour in your script
concerning ``print`` and integer division even when running it in Python 2
is to put the following line at the top of your Python script: 
   
.. code-block:: python

   from __future__ import print_function, division
