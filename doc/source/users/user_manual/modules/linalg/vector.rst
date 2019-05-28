Vectors
=======

The :doxy:`GVector` class implements a vector as a one-dimensional array of
successive double precision floating point values. You may allocate
an vectors as follows

**C++**

.. code-block:: cpp
   :linenos:

   GVector vector1;                           // Allocates empty vector
   GVector vector2(10);                       // Allocates vector with 10 elements
   GVector vector3 = vector1;                 // Allocates vector by copying

**Python**

.. code-block:: python
   :linenos:

   vector1 = gammalib.GVector()               # Allocates empty vector
   vector2 = gammalib.GVector(10)             # Allocates vector with 10 elements
   vector3 = vector1.copy()                   # Allocates vector by copying

On allocation, all elements of a vector are set to 0.

There are three special constructors which you can use to allocate vectors
with 1, 2 or 3 elements:

**C++**

.. code-block:: cpp
   :linenos:

   GVector vector1(1.0);                      // Allocates vector [1.0]
   GVector vector2(1.0, 2.0);                 // Allocates vector [1.0, 2.0]
   GVector vector3(1.0, 2.0, 3.0);            // Allocates vector [1.0, 2.0, 3.0]

**Python**

.. code-block:: python
   :linenos:

   vector1 = gammalib.GVector(1.0)            # Allocates vector [1.0]
   vector2 = gammalib.GVector(1.0, 2.0)       # Allocates vector [1.0, 2.0]
   vector3 = gammalib.GVector(1.0, 2.0, 3.0)  # Allocates vector [1.0, 2.0, 3.0]

You can access vector elements using the :doxy:`GVector::operator[]` operator.
As usual in C++ and Python, indices start with 0:

**C++**

.. code-block:: cpp
   :linenos:

   std::cout << vector2[1] << std::endl;      // Print 2nd element
   2.0                                        // Result

**Python**

.. code-block:: python
   :linenos:

   print(vector2[1])                          # Print 2nd element
   2.0                                        # Result

Note that in C++ the :doxy:`GVector::operator[]` operator does not check the
validity of the element index. If index checking is needed in C++, use the
:doxy:`GVector::at` method:

**C++**

.. code-block:: cpp
   :linenos:

   std::cout << vector2.at(1) << std::endl;   // Print 2nd element with index checking
   2.0                                        // Result

In Python, index checking is always performed.

All GammaLib objects are printable, and vectors can be printed using

**C++**

.. code-block:: cpp
   :linenos:

   std::cout << vector3 << std::endl;         // Print vector
   (1, 2, 3)                                  // Result

**Python**

.. code-block:: python
   :linenos:

   print(vector3)                              # Print result
   (1, 2, 3)                                   # Result


You can handle vectors pretty much the same way you handle floating point 
variables. :doxy:`GVector` supports various arithmetic operations:

**C++**

.. code-block:: cpp
   :linenos:

   GVector a(1.0, 2.0, 3.0);                  // A vector
   GVector b(2.0, 4.0, 6.0);                  // Another vector
   GVector c;                                 // Result vector
   double  s = 2.0;                           // A floating point value
   c  = a + b;                                // Vector + Vector addition
   c += a;                                    // Unary vector addition
   c  = a + s;                                // Vector + Scalar addition
   c += s;                                    // Unary scalar addition
   c  = s + b;                                // Scalar + Vector addition (only C++)
   c  = a - b;                                // Vector - Vector subtraction
   c -= a;                                    // Unary vector subtraction
   c  = a - s;                                // Vector - Scalar subtraction
   c -= s;                                    // Unary scalar subtraction
   c  = s - b;                                // Scalar - Vector subtraction (only C++)
   s  = a * b;                                // Vector * Vector multiplication (dot product)
   c  = a * s;                                // Vector * Scalar multiplication
   c *= s;                                    // Unary scalar multiplication
   c  = s * b;                                // Scalar * Vector multiplication (only C++)
   c  = a / s;                                // Vector * Scalar division
   c /= s;                                    // Unary scalar division
   c  = s;                                    // Assigns scalar to all vector elements
   c  = -a;                                   // Vector negation

**Python**

.. code-block:: python
   :linenos:

   a  = gammalib.GVector(1.0, 2.0, 3.0)       # A vector
   b  = gammalib.GVector(2.0, 4.0, 6.0)       # Another vector
   c  = gammalib.GVector()                    # Result vector
   s  = 2.0                                   # A floating point value
   c  = a + b                                 # Vector + Vector addition
   c += a                                     # Unary vector addition
   c  = a + s                                 # Vector + Scalar addition
   c += s                                     # Unary scalar addition
   c  = a - b                                 # Vector - Vector subtraction
   c -= a                                     # Unary vector subtraction
   c  = a - s                                 # Vector - Scalar subtraction
   c -= s                                     # Unary scalar subtraction
   s  = a * b                                 # Vector * Vector multiplication (dot product)
   c  = a * s                                 # Vector * Scalar multiplication
   c *= s                                     # Unary scalar multiplication
   c  = a / s                                 # Vector * Scalar division
   c /= s                                     # Unary scalar division
   c  = -a                                    # Vector negation

Most of these operations operate element-wise. For example, scalar 
additions or subtractions add or subtract a given scalar value from every 
vector element. And scalar multiplications and divisions multiply or 
divide every vector element by a given value. The dot product implements 
the usual formula

.. math::
   s = \sum_{i=0}^{N-1} a_i b_i

(where :math:`N` is the number of vector elements).
It is obvious that the dot product, as well as vector addition and 
subtraction, require vectors of identical dimensions. If vectors are not 
identical, an :doxy:`GException::vector_mismatch` exception will be thrown.

You can compare vectors using

Finally, you can use the comparison operators

**C++**

.. code-block:: cpp
   :linenos:

   int equal   = (a == b);                    // True if all elements equal
   int unequal = (a != b);                    // True if at least one elements unequal

**Python**

.. code-block:: python
   :linenos:

   equal   = (a == b)                         # True if all elements equal
   unequal = (a != b)                         # True if at least one elements unequal

In addition to the operators, you can apply the following mathematical functions
to vectors

**C++**

.. code-block:: cpp
   :linenos:

   c = acos(a);                               // Vector with acos of all vector elements
   c = acosh(a);                              // Vector with acosh of all vector elements
   c = asin(a);                               // Vector with asin of all vector elements
   c = asinh(a);                              // Vector with asinh of all vector elements
   c = atan(a);                               // Vector with atan of all vector elements
   c = atanh(a);                              // Vector with atanh of all vector elements
   c = cos(a);                                // Vector with cos of all vector elements
   c = cosh(a);                               // Vector with cosh of all vector elements
   c = exp(a);                                // Vector with exp of all vector elements
   c = abs(a);                                // Vector with absolute values of all vector elements
   c = log(a);                                // Vector with natural logarithm of all vector elements
   c = log10(a);                              // Vector with base 10 logarithm of all vector elements
   c = sin(a);                                // Vector with sin of all vector elements
   c = sinh(a);                               // Vector with sinh of all vector elements
   c = tan(a);                                // Vector with tan of all vector elements
   c = tanh(a);                               // Vector with tanh of all vector elements
   c = pow(a, 2.7);                           // Vector with power of 2.7 of all vector elements (only C++)
   c = cross(a, b);                           // Vector with cross product of 3 element vectors
   s = norm(a);                               // Vector norm |a|
   s = min(a);                                // Minimum element of vector
   s = max(a);                                // Maximum element of vector
   s = sum(a);                                // Sum of vector elements

**Python**

.. code-block:: python
   :linenos:

   c = a.acos()                               # Vector with acos of all vector elements
   c = a.acosh()                              # Vector with acosh of all vector elements
   c = a.asin()                               # Vector with asin of all vector elements
   c = a.asinh()                              # Vector with asinh of all vector elements
   c = a.atan()                               # Vector with atan of all vector elements
   c = a.atanh()                              # Vector with atanh of all vector elements
   c = a.cos()                                # Vector with cos of all vector elements
   c = a.cosh()                               # Vector with cosh of all vector elements
   c = a.exp()                                # Vector with exp of all vector elements
   c = a.abs()                                # Vector with absolute values of all vector elements
   c = a.log()                                # Vector with natural logarithm of all vector elements
   c = a.log10()                              # Vector with base 10 logarithm of all vector elements
   c = a.sin()                                # Vector with sin of all vector elements
   c = a.sinh()                               # Vector with sinh of all vector elements
   c = a.tan()                                # Vector with tan of all vector elements
   c = a.tanh()                               # Vector with tanh of all vector elements
   c = a.cross(b)                             # Vector with cross product of 3 element vectors
   s = a.norm()                               # Vector norm |a|
   s = a.min()                                # Minimum element of vector
   s = a.max()                                # Maximum element of vector
   s = a.sum()                                # Sum of vector elements

Finally, the following methods exist:

**C++**

.. code-block:: cpp
   :linenos:

   int size     = a.size();                   // Returns number of vector elements
   int nonzeros = a.non_zeros();              // Returns number of non-zero vector elements

**Python**

.. code-block:: python
   :linenos:

   size     = a.size(a)                       # Returns number of vector elements
   nonzeros = a.non_zeros()                   # Returns number of non-zero vector elements
