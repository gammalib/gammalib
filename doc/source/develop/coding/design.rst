Design conventions
==================

.. _sec_configuration:

Code configuration
------------------

The code configuration is controlled via the ``config.h`` header file that
is created during the configuration step of the GammaLib build. To make 
configuration
options available the following code has to be added to the source code
file:

.. code-block:: cpp

    #ifdef HAVE_CONFIG_H
    #include <config.h>
    #endif

.. warning::

   **The config.h file should never be included in header files.**
   The header files are used by the outside world without any knowledge
   about ``config.h``.

One of the most commonly used configuration options used throughout
GammaLib is related to range checking. Range checks are usually
performed when accessing array elements, assuring that no elements
outside the valid range are accessed. Range checking, however, is time
consuming, in particular if many elements are accessed subsequently.
GammaLib therefore allows to disable range checks. This can be done by
specifying ``./configure --disable-range-check`` when configuring GammaLib
for compilation. The ``--disable-range-check`` option undefines
G_RANGE_CHECK, and optional range checking is thus achieved by adding
for example

.. code-block:: cpp

    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_num) {
        throw GException::out_of_range("GVector::operator(int)", inx, m_num);
    }
    #endif

to the code.

The following table gives a list of important configuration options that
are available in config.h and that can be used to tune the code:

==================== =========================== ================================================= 
Definition           Option                      Usage
==================== =========================== ================================================= 
``G_DEBUG``          ``--enable-debug``          Code debugging
``G_PROFILE``        ``--enable-profiling``      Code profiling
``G_RANGE_CHECK``    ``--enable-range-check``    Performs range checking
``G_NAN_CHECK``      ``--enable-nan-check``      Check for ``NaN`` and ``Inf`` values
``G_SMALL_MEMORY``   ``--enable-small-memory``   Optimizes for small memory
``HAVE_OPENMP``      ``--enable-openmp``         Has OpenMP multi-threading support
``HAVE_LIBREADLINE`` ``--with-readline``         Has readline library
``HAVE_LIBCFITSIO``  ``--with-cfitsio``          Has cfitsio library
``HAVE_PYTHON``      ``--enable-python-binding`` Has Python bindings
``PACKAGE``          n.a.                        gammalib
``PACKAGE_PREFIX``   n.a.                        Installation location (e.g. ``/usr/local/gamma``)
``PACKAGE_STRING``   n.a.                        Full name and version of GammaLib
``PACKAGE_VERSION``  n.a.                        Version of GammaLib (format: ``x.y.z``)
==================== =========================== ================================================= 

Note that ``enable`` may be replaced by ``disable`` and ``with`` by ``without`` for
switching off an option.

C++ classes
-----------

Members
^^^^^^^

Class members should be either private or protected, the latter being
generally used when a derived class should be able to access base class
data.

Members should be prefixed by ``m_`` and should be in lower case. For long
member names, additional underscores may be added. Examples of valid
member names are

::

    m_num
    m_response
    m_grid_length
    m_axis_dir_qual

Initialisation, copying and deleting of class members should be gathered
in a single place to avoid memory leaks. For this purpose, each C++
class should have the following private or protected methods for memory
management:

-  ``init_members()``

   initializes all member variables and pointers to well defined initial
   values. The class should be fully operational and consistent with
   these initial values. All pointers that will hold dynamically
   allocated memory should be initialised to ``NULL``.

-  ``copy_members(const &A a)``

   copies all members from an instance ``a`` into the class.

-  ``free_members()``

   frees all memory that has been allocated by the class. Memory
   pointers should be set to ``NULL`` after the memory was deleted to signal
   that no valid memory is associated to the pointer. This allows for
   checking if memory has been allocated before it is accessed.

(in the above notation, ``A`` is the class name and ``a`` is an instance of the
class).

An example for valid ``init_members()``, ``copy_members(const &A a)`` and
``free_members()`` methods is:

.. code-block:: cpp

    void GEbounds::init_members(void)
    {
        m_num = 0;
        m_min = NULL;
        m_max = NULL;
        return;
    }

    void GEbounds::copy_members(const GEbounds& ebds)
    {
        m_num  = ebds.m_num;
        if (m_num > 0) {
            m_min = new GEnergy[m_num];
            m_max = new GEnergy[m_num];
            for (int i = 0; i < m_num; ++i) {
                m_min[i] = ebds.m_min[i];
                m_max[i] = ebds.m_max[i];
            }
        }
        return;
    }

    void GEbounds::free_members(void)
    {
        if (m_min != NULL) delete [] m_min;
        if (m_max != NULL) delete [] m_max;
        m_min = NULL;
        m_max = NULL;
        return;
    }

In this example, one may probably want to add a ``alloc_members()`` method
for memory allocation:

.. code-block:: cpp

    void GEbounds::alloc_members(const int& num)
    {
        if (num > 0) {
            m_min = new GEnergy[num];
            m_max = new GEnergy[num];
            for (int i = 0; i < num; ++i) {
                m_min[i] = 0.0;
                m_max[i] = 0.0;
            }
            m_num = num;
        }
        return;
    }

This example illustrates several design conventions:

-  Always check if a pointer is not ``NULL`` before de-allocating memory.

-  After de-allocation, always set the pointer immediately to ``NULL``.

-  Never allocate zero elements (check if the number of elements to be
   allocated is positive).

-  Always initialise allocated memory to well-defined values (do not
   expect that the compiler will do this for you).

Constructors, destructors and operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each class should have at least a void constructor, a copy constructor,
a destructor and an assignment operator. Additional constructors and
operators can be implemented as required. The following example shows
the basic implementation for these 4 methods. Due to the usage of the
``init_members()``, ``copy_members(const &A a)`` and the ``free_members()``
methods, most classes will have exactly this kind of syntax:

.. code-block:: cpp

    GEbounds::GEbounds(void)
    {
        init_members();
        return;
    }

    GEbounds::GEbounds(const GEbounds& ebds)
    {
        init_members();
        copy_members(ebds);
        return;
    }

    GEbounds::~GEbounds(void)
    {
        free_members();
        return;
    }

    GEbounds& GEbounds::operator= (const GEbounds& ebds)
    {
        if (this != &ebds) {
            free_members();
            init_members();
            copy_members(ebds);
        }
        return *this;
    }

Inheritance
^^^^^^^^^^^

Class inheritance is central feature of the C++ language, and is largely
used throughout GammaLib. Multiple inheritance is not used at the moment in GammaLib.
Because of the added complexity of multiple inheritance in C++ and in Python
there would have to be very good reasons to use it in GammaLib.

Although the inheritance philosophy may differ from class to class, the
following guidelines should be respected as far as possible:

-  The base class and derived class destructors should be declared
   virtual.

-  Avoid overloading of base class methods by derived class methods.
   Preferentially, define base class methods as pure virtual.

-  All base class methods that should be implemented in the derived
   class should be declared virtual. Exceptions are the ``init_members()``,
   the ``copy_members()`` and the ``free_members()`` methods that will be
   implemented in the base class and the derived class.

-  Base classes manage base class members, derived classes manage
   derived class members. By managing we mean here in particular memory
   allocation and de-allocation, but also proper initialization.

-  Derived class constructors should invoke base class constructors for
   proper base class initialization. A void constructor should look like

   .. code-block:: cpp

       GEventList::GEventList(void) : GEvents()
       {
           init_members();
           return;
       }

   and a copy constructor should look like

   .. code-block:: cpp

       GEventList::GEventList(const GEventList& list) : GEvents(list)
       {
           init_members();
           copy_members(list);
           return;
       }

-  Derived class operators should invoke base class operators, as
   illustrated by the following example:

   .. code-block:: cpp

       GEventList& GEventList::operator=(const GEventList& list)
       {
           if (this != &list) {
               this->GEvents::operator=(list);
               free_members();
               init_members();
               copy_members(list);
           }
           return *this;
       }

-  The ``clear()`` method of a derived class show invoke the ``free_members()``
   method of the base class, as illustrated by the following example:

   .. code-block:: cpp

       void GCTAEventList::clear(void)
       {
           free_members();
           this->GEventList::free_members();
           this->GEvents::free_members();
           this->GEvents::init_members();
           this->GEventList::init_members();
           init_members();
           return;
       }

-  Avoid as far as possible methods that are only defined in the derived
   class.

.. warning::

   **For a derived class, init_members(), copy_members(const &A a) and
   free_members() should only act on derived class members but not on base
   class members**. Any exception to this rule needs very careful
   documentation since it can easily be the source of memory leaks.


Method naming conventions
^^^^^^^^^^^^^^^^^^^^^^^^^

Uniform public method names should be provided throughout GammaLib for
all classes. Unless the public method names are very long (which should
be avoided), names should not comprise underscores as separators. Public
method names are all lowercase.

Private or protected method names may differ from this since they are hidden
within the class. Yet also here, all method names should be lowercase,
and the use of underscores should be limited.

Methods that set or retrieve class attributes should be named after the
attribute. Here an example for the attribute ``m_name``:

.. code-block:: cpp

    public:
        void        name(const std::string& name);
        std::string name(void) const;
    protected:
        m_name;

A method name that is used in multiple classes should always perform an
equivalent action. Here is a list of method names that are widely used
in GammaLib, together with their typical usage. The last column specifies where
these methods are used. Note that **the ``clear()``, ``clone()``, and ``print()``
methods should be implemented for all classes**.

============ ============================================= ==============
Method       Usage                                         Implementation
============ ============================================= ==============
``clear``    Set object to initial empty state             all classes
``clone``    Provides a deep copy of the class             all classes
``print``    Print object into string                      all classes (see :ref:`sec_output`)
``is_empty`` Checks for emptiness of object                :ref:`sec_containers`
``append``   Append element to list of elements            :ref:`sec_containers`
``extend``   Append container elements to list of elements :ref:`sec_containers`
``insert``   Insert element to list of elements            :ref:`sec_containers`
``remove``   Remove element from list of elements          :ref:`sec_containers`
``reserve``  Reserve memory for a number of elements       :ref:`sec_containers`
``load``     Load data from file (open, read, close)       if applicable
``save``     Save data into file (open, write, close)      if applicable
``open``     Open file                                     if applicable
``read``     Read data from open file                      if applicable
``write``    Write data into open file                     if applicable
``close``    Close file                                    if applicable
``name``     Name of object                                if applicable
``type``     Type of object                                if applicable
``size``     Size of object                                if applicable
``real``     Returns double precision value                if applicable
``integer``  Returns ``int`` value                         if applicable
``string``   Returns ``std::string`` value                 if applicable
============ ============================================= ==============

Note the difference between ``load()`` and ``read()`` and between ``save()`` and
``write()``. The ``load()`` and ``save()`` methods should take as arguments a file
name, and they will open the file, read or write some data, and then
close the file. In contrast, ``read()`` and ``write()`` will operate on files
that are already open, and after the read or write operation the files
will remain open. Typically, these methods take a ``GFits*`` or a
``GFitsHDU*`` pointer as argument.

Methods that perform checks should return a bool type and should start
with the prefix ``is_`` or ``has_``, e.g. ``is_empty()`` or ``has_min()``.

Method const declarations
^^^^^^^^^^^^^^^^^^^^^^^^^

All methods that do not alter accessible class members should be
declared ``const``. With accessible we mean here class members that can be
read or written in some way by one of the methods. Non-accessible class
members would be members that are only used internally, and for which no
consistent state has to be preserved for the outside world. These could
for example be members that hold pre-computed values.

Methods that do not alter accessible members, but that modify
non-accessible members, should also be declared ``const``. The
non-accessible class members need then to be declared ``mutable`` to avoid
compiler errors. Alternatively, the ``const_cast`` declaration can be used
to allow member modifications within a ``const`` method.

As example we show here part of the definition of :doxy:`GModelSpectralPlawPhotonFlux`:

.. code-block:: cpp

    class GModelSpectralPlawPhotonFlux : public GModelSpectral {
    public:
        virtual double eval(const GEnergy& srcEng) const;
        virtual void   read(const GXmlElement& xml);
    protected:
        // Protected members
        GModelPar       m_integral;        //!< Integral flux
        GModelPar       m_index;           //!< Spectral index
        GModelPar       m_emin;            //!< Lower energy limit (MeV)
        GModelPar       m_emax;            //!< Upper energy limit (MeV)

        // Cached members used for pre-computations
        mutable double  m_log_emin;        //!< Log(emin)
        mutable double  m_log_emax;        //!< Log(emax)
        mutable double  m_pow_emin;        //!< emin^(index+1)
        mutable double  m_pow_emax;        //!< emax^(index+1)
        mutable double  m_norm;            //!< Power-law normalization (for pivot energy 1 MeV)
        mutable double  m_g_norm;          //!< Power-law normalization gradient
        mutable double  m_power;           //!< Power-law factor
        mutable double  m_last_integral;   //!< Last integral flux
        mutable double  m_last_index;      //!< Last spectral index (MeV)
        mutable double  m_last_emin;       //!< Last lower energy limit (MeV)
        mutable double  m_last_emax;       //!< Last upper energy limit (MeV)
        mutable GEnergy m_last_energy;     //!< Last source energy
        mutable double  m_last_value;      //!< Last function value
        mutable double  m_last_g_integral; //!< Last integral flux gradient
        mutable double  m_last_g_index;    //!< Last spectral index gradient

This class has an internal cache for precomputation, which is
potentially updated when :doxy:`GModelSpectralPlawPhotonFlux::eval` is called. Here the corresponding code:

.. code-block:: cpp

    double GModelSpectralPlawPhotonFlux::eval(const GEnergy& srcEng) const
    {
        // Update precomputed values
        update(srcEng);

        // Compute function value
        double value = integral() * m_norm * m_power;

        // Return
        return value;
    }

As the pre-computation cache is not exposed to the external world but
fully handled within the class, ``eval()`` is declared const as it does not
modify any of the model parameters (which are ``m_integral``, ``m_index``,
``m_emin``, and ``m_emax``). It may however modified some of the cache
members, that’s why these members are declared ``mutable``. As there is
however no way to access these cache values from the outside (no method
exists to access them), the ``eval()`` method does not modify any
*observable* property of the class, hence it is declared ``const``.

Method arguments and return values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If possible, method arguments should always be passed by reference. To
protect references from changes by the method, **arguments passed by
reference should always be declared const**. Pointers should only be
used as arguments if ``NULL`` should be a possible argument value. Also
pointers should always be declared ``const``. Here an example based on the
definition of :doxy:`GObservation`:

.. code-block:: cpp

    class GObservation {
    public:
        void events(const GEvents* events);
        void statistics(const std::string& statistics);
    protected:
        std::string m_statistics;   //!< Optimizer statistics (default=poisson)
        GEvents*    m_events;       //!< Pointer to event container
    };

The ``statictics`` value is passed by reference because the class will hold
the actual value, while ``events`` is passed as a pointer because the class
will hold the pointer.

Numeric argument types should be either ``int`` or ``double``. Unless absolutely
necessary, avoid ``short int``, ``long``, and ``float``.

If a method returns a class member, the return value should be passed by
reference. Unless we explicitly want to modify a class member through
the method call, return values passed by reference should be declared
``const``.

If a method returns a base class object, a pointer should be returned.
**Do never return base class objects by reference, as this will lead to
code slicing if the method is used for object assignment.** Unless we
explicitly want to modify a class member through the method call, the
returned pointer should be declared ``const``.

Here an example based on the definition of :doxy:`GObservation`:

.. code-block:: cpp

    class GObservation {
    public:
        virtual double     ontime(void) const = 0;
        const GEvents*     events(void) const;
        const std::string& statistics(void) const { return m_statistics; }
    protected:
        std::string m_statistics;   //!< Optimizer statistics (default=poisson)
        GEvents*    m_events;       //!< Pointer to event container
    };

The :doxy:`GObservation::ontime` method does return a ``double`` by value as the ``ontime`` property
is not stored explicitly in the class (hence no reference can be
returned to it). On the other hand, the ``statistics`` method returns by
reference as the ``statistics`` property is stored as a data member (hence a
reference can be returned). Although we could have returned a reference
to the event container, this would lead to code slicing. Therefore, the
``events`` method returns a pointer. All returned references or pointers
are declared ``const`` to prevent modification of class members.

.. _sec_containers:

Container classes
^^^^^^^^^^^^^^^^^

Container classes are classes that contain list of elements. Two cases
are distinguished here: containers holding objects, and containers
holding pointers to objects.

Containers holding objects
~~~~~~~~~~~~~~~~~~~~~~~~~~

Containers holding objects should have element access operators
``operator[]`` implemented that return container elements by reference. A
non-const and a const version of the operator should exist. Eventually,
``at()`` methods could be added that always perform range checking. Here is
a list of mandatory methods for container classes holding objects:

========================================= ==========================================
Method                                    Usage
========================================= ==========================================
``operator[](const int&)``                Element access operator
``const operator[](const int&) const``    Element access operator (const version)
``void clear()``                          Delete all objects in container
``bool is_empty()``                       Checks whether container is empty
``int size()``                            Return number of elements in container
``void append(const e&)``                 Append an element to the container
``void insert(const int&, const e&)``     Insert an element into the container
``void extend(const C&)``                 Append another container to the container
``void remove(const int&)``               Removes an element from the container
``void reserve(const int&)``              Reserve memory space in a container
``std::string print()``                   Print container (see :ref:`sec_output`)
========================================= ==========================================

Containers holding pointers
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Containers holding pointers are different from those holding objects in that
their ``operator[]`` operators return a pointer, and in that they implement
a ``set()`` method for value setting. Here is a list of mandatory methods
for container classes holding pointers:

========================================= =====
Method                                    Usage
========================================= =====
``e* operator[](const int&)``             Element access operator
``const e* operator[](const int&) const`` Element access operator (const version)
``void set(const int&, const e&)``        Set an element of the container
``void clear()``                          Delete all objects in container
``bool is_empty()``                       Checks whether container is empty
``int size()``                            Return number of elements in container
``void append(const e&)``                 Append an element to the container
``void insert(const int&, const e&)``     Insert an element into the container
``void extend(const C&)``                 Append another container to the container
``void remove(const int&)``               Removes an element from the container
``std::string print()``                   Print container (see :ref:`sec_output`)
========================================= =====

.. _sec_output:

Output
^^^^^^

Output stream and logging operators should be implemented for every
class as friend operators (see :ref:`sec_header`). 
The usage of friend operators (instead of member operators) allows for correct
handling of code such as

.. code-block:: cpp

    log << std::endl << "This is a text" << std::endl;

To support these friend operators (and to support also the Python
interface), each class should have a ``print()`` method:

.. code-block:: cpp

    std::string print(const GChatter& chatter = NORMAL) const;

In case that the class derives from one of the standard interface classes 
:doxy:`GBase`, :doxy:`GContainer` and :doxy:`GRegistry`,
the output stream and logging operators are automatically implemented on
the level of the base class. 
In all other cases, the developer needs to implement these operators on the 
class level.
The operators will use the ``print()`` method to enable printing in the
following form:

.. code-block:: cpp

    std::ostream& operator<< (std::ostream& os, const GFits& fits)
    {
        os << fits.print();
        return os;
    }
    GLog& operator<< (GLog& log, const GFits& fits)
    {
        log << fits.print(log.chatter());
        return log;
    }

.. _sec_exceptions:

Exceptions
^^^^^^^^^^

Exceptions are largely used in GammaLib to handle the occurrence of
unexpected events. GammaLib exceptions are implemented by the :doxy:`GException`
class. For each new exception type, a new exception subclass is added.

Each exception returns the method name in which the exception occurs and
an exception message. The exception message is generally built from
values that are passed as arguments to the exception constructor.

Below a list of conventions for implementing and using exceptions:

-  Re-use existing exceptions as far as possible.
-  Pass exception arguments by reference.
-  Use exceptions only for events that cannot be handled by a method. Do
   not use exceptions to check a value or a state. Implement appropriate
   methods instead.
-  Never use exceptions in a destructor.
-  De-allocate all memory that is not de-allocated by the destructor
   before throwing an exception.
-  Always catch exceptions by reference.

Python interface for C++ classes
--------------------------------

Container classes
^^^^^^^^^^^^^^^^^

The Python interface for container classes should implement the
following class extensions:

=============== ============== =====
Extension       Method         Usage
=============== ============== =====
``__str__``     ``print``      Convert object to string
``__getitem__`` ``operator[]`` Element access (get)
``__setitem__`` ``operator[]`` Element access (set)
``__len__``     ``size``       Container size
=============== ============== =====

While ``__str__`` and ``__len__`` are implemented in the :doxy:`GContainer`
base class, ``__getitem__`` and ``__setitem__`` need to be implemented
as class extensions in the SWIG file (see :ref:`sec_python`).
