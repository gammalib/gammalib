Matrix arithmetics
==================

The following description of matrix arithmetics applies to all storage classes.
The following matrix operators have been implemented:

**C++**

.. code-block:: cpp
   :linenos:

   GMatrix A(10,10);                          // A 10 x 10 matrix
   GMatrix B(10,10);                          // Another 10 x 10 matrix
   GMatrix C;                                 // Result matrix
   GVector v(10);                             // Vector with 10 elements
   double  s = 2.0;                           // Floating point value
   C  = A + B;                                // Matrix addition
   C  = A - B;                                // Matrix subtraction
   C  = A * B;                                // Matrix multiplication
   C  = A * v;                                // Vector multiplication
   C  = A * s;                                // Right-handed scalar multiplication
   C  = s * A;                                // Left-handed scalar Matrix multiplication (only C++)
   C  = A / s;                                // Scalar division
   C  = -A;                                   // Negation
   A += B;                                    // Unary matrix addition
   A -= B;                                    // Unary matrix subtraction
   A *= B;                                    // Unary matrix multiplications
   A *= s;                                    // Unary matrix scalar multiplication
   A /= s;                                    // Unary matrix scalar division

**Python**

.. code-block:: python
   :linenos:

   A  = gammalib.GMatrix(10,10)               # A 10 x 10 matrix
   B  = gammalib.GMatrix(10,10)               # Another 10 x 10 matrix
   v  = gammalib.GVector(10)                  # Vector with 10 elements
   s  = 2.0                                   # Floating point value
   C  = A + B                                 # Matrix addition
   C  = A - B                                 # Matrix subtraction
   C  = A * B                                 # Matrix multiplication
   C  = A * v                                 # Vector multiplication
   C  = A * s                                 # Scalar multiplication
   C  = A / s                                 # Scalar division
   C  = -A                                    # Negation
   A += B                                     # Unary matrix addition
   A -= B                                     # Unary matrix subtraction
   A *= B                                     # Unary matrix multiplications
   A *= s                                     # Unary matrix scalar multiplication
   A /= s                                     # Unary matrix scalar division

You can use the comparison operators

**C++**

.. code-block:: cpp
   :linenos:

   int equal   = (A == B);                    // True if all elements equal
   int unequal = (A != B);                    // True if at least one elements unequal

**Python**

.. code-block:: python
   :linenos:

   equal   = (A == B)                         # True if all elements equal
   unequal = (A != B)                         # True if at least one elements unequal

In addition to the operators, you can apply the following mathematical functions
to a matrix:

**C++**

.. code-block:: cpp
   :linenos:

   C = A.abs();                               // Matrix with absolute values of all matrix elements
   C = A.transpose();                         // Transpose matrix
   C = A.invert();                            // Invert matrix
   v = A.solve(v);                            // Solve matrix equation x = M x v
   s = A.min();                               // Minimum element of matrix
   s = A.max();                               // Maximum element of matrix
   s = A.sum();                               // Sum of matrix elements

**Python**

.. code-block:: python
   :linenos:

   C = A.abs()                                # Matrix with absolute values of all matrix elements
   C = A.transpose()                          # Transpose matrix
   C = A.invert()                             # Invert matrix
   v = A.solve(v)                             # Solve matrix equation x = M x v
   s = A.min()                                # Minimum element of matrix
   s = A.max()                                # Maximum element of matrix
   s = A.sum()                                # Sum of matrix elements

.. warning::
   The ``invert()`` and ``solve()`` methods are so far only implemented for
   sparse matrices.
