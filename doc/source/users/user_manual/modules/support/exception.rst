Exception handling
~~~~~~~~~~~~~~~~~~

GammaLib uses C++ exceptions to handle any kind of errors or unexpected
values that are encountered. Here an example of how an exception is thrown 
in the :doxy:`GModels::at` method:

**C++**

.. code-block:: cpp
   :linenos:

    #define G_AT "GModels::at(int&)"
    GModel* GModels::at(const int& index)
    {
        if (index < 0 || index >= size()) {
            throw GException::out_of_range(G_AT, "Model index", index, size());
        }
        return m_models[index];
    }

Line 1 defines the name of the method that is always passed to the 
exception handle to track the origin of the exception. The definition 
should always provide the class name, the method name and the argument
types (as several methods with the same name but a different set of 
arguments may exist). Variable names or const declarations are omitted 
from the definition.

The method is implemented in lines 2-8. Before accessing a model in line 
7, the method checks whether the provided index is in the valid range. 
Note that lower and upper boundary of the index value is systematically 
checked in all GammaLib methods that perform index checking. If one of the 
boundary conditions is violated, the ``throw`` statement is used to throw 
an object of type :doxy:`GException::out_of_range`. The object is constructed 
by passing the method name (defined by ``G_AT``), a text string that 
describes the parameter that is out of the valid range, the value of the 
parameter, and the maximum number of elements that are expected in the 
range. This specific instance of the :doxy:`GException::out_of_range` class
assumes that the lower boundary of the valid range is 0, hence it does not 
need to be specified explicitely as an argument.

The actual GammaLib code implements a wealth of possible exceptions, yet 
in a future version of the code, this wealth should be reduced to a 
limited set of standard exceptions. The first class of exceptions are
logic exceptions, which are those that the client could in principle have 
tested before calling the method. These comprise:

======================== =====
Logic exceptions         Usage
======================== =====
:doxy:`invalid_value`    An invalid value has been encountered in the method.
:doxy:`invalid_argument` One of the arguments passed to the method is invalid.
:doxy:`out_of_range`     An index is outside the expected range.
:doxy:`fits_error`       An error has occured in FITS file handling.
======================== =====

The second class of exceptions are runtime exceptions, which are those 
that are not testable by the client. Typical runtime exceptions are 
underflow or overflow errors. So far, only one runtime exception is 
implemented in GammaLib:

=============================== =====
Runtime exceptions              Usage
=============================== =====
:doxy:`feature_not_implemented` The method has not been implemented.
=============================== =====
