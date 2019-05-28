Matrix storage classes
======================

A matrix is a two-dimensional array of double precision floating point 
values, arranged in rows and columns. Three matrix storage classes are 
implemented in GammaLib:

- :doxy:`GMatrix` which explicitly stores all elements of a matrix

- :doxy:`GMatrixSymmetric` which implements a symmetric matrix and only stores 
  the lower-left triangle of matrix elements

- :doxy:`GMatrixSparse` which implements a sparse matrix and only stores the 
  non-zero elements of a matrix

All matrix classes derive from the abstract :doxy:`GMatrixBase` class.

In the most general case, the rows and columns of a matrix are stored in
a continuous array of :math:`{\tt rows} \times {\tt columns}` memory
locations. This storage type is referred to as a *full matrix*, and is
implemented by the class :doxy:`GMatrix`. Operations on full matrixes are in
general relatively fast, but memory requirements may be important to
hold all the elements. In general matrixes are stored column-wise
(or in column-major format). For example, the matrix ::

        1  2  3  4  5
        6  7  8  9 10
       11 12 13 14 15 

is stored in memory as ::

        |  1  6 11 |  2  7 12 |  3  8 13 |  4  9 14 |  5 10 15 |

Many physical or mathematical problems treat with a subclass of matrixes
that is symmetric, i.e. for which the element :math:`(row,col)` is identical to
the element :math:`(col,row)`. In this case, the duplicated elements need not to
be stored. The class :doxy:`GMatrixSymmetric` implements such a storage type.
:doxy:`GMatrixSymmetric` stores the lower-left triangle of the matrix in 
column-major format. For illustration, the matrix ::

        1  2  3  4
        2  5  6  7
        3  6  8  9
        4  7  9 10

is stored in memory as ::

        |  1  2  3  4 |  5  6  7 |  8  9 | 10 |

This divides the storage requirements to hold the matrix elements by
almost a factor of two.

Finally, quite often one has to deal with matrixes that contain a large
number of zeros. Such matrixes are called *sparse matrixes*. If only the
non-zero elements of a sparse matrix are stored the memory requirements
are considerably reduced. This goes however at the expense of matrix
element access, which has become now more complex. In particular,
filling efficiently a sparse matrix is a non-trivial problem (see
:ref:`sec_matrix_filling`). Sparse matrix storage is implemented by
the :doxy:`GMatrixSparse` class. A :doxy:`GMatrixSparse` object contains
three one-dimensional arrays to store the matrix elements: a double type
array that contains in continuous column-major order all non-zero
elements, an int type array that contains for each non-zero element the
row number of its location, and an int type array that contains the
storage location of the first non-zero element for each matrix column.
To illustrate this storage format, the matrix ::

        1  0  0  7
        2  5  0  0
        3  0  6  0
        4  0  0  8

is stored in memory as ::

        |  1  2  3  4 |  5 |  6 |  7  8 |  Matrix elements
        |  0  1  2  3 |  1 |  2 |  0  3 |  Row indices for all elements
        |  0          |  4 |  5 |  6    |  Storage location of first element of each column

This example is of course not very economic, since the total number of
Bytes used to store the matrix is
:math:`8 \times 8 + (8 + 4) \times 4 = 112` Bytes, while a full
:math:`4 \times 4` matrix is stored in
:math:`(4 \times 4) \times 8 = 128` Bytes (recall: a double type values
takes 8 Bytes, an int type value takes 4 Bytes). For realistic large
systems, however, the gain in memory space can be dramatical.

The usage of the :doxy:`GMatrix`, :doxy:`GMatrixSymmetric` and :doxy:`GMatrixSparse`
classes is analoguous in that they implement basically all functions and 
methods in an identical way. So from the semantics the user has not to worry 
about the storage class. However, matrix element access speeds are not
identical for all storage types, and if performance is an issue (as it
certainly always will be), the user has to consider matrix access more
carefully (see :ref:`sec_matrix_filling`).

You allocate a matrix using the constructors

**C++**

.. code-block:: cpp
   :linenos:

   GMatrix          A(10,20);                 // Full 10 x 20 matrix
   GMatrixSymmetric B(10,10);                 // Symmetric 10 x 10 matrix
   GMatrixSparse    C(1000,10000);            // Sparse 1000 x 10000 matrix

   GMatrix          D(0,0);                   // WRONG: empty matrix not allowed
   GMatrixSymmetric E(20,22);                 // WRONG: symmetric matrix requested

**Python**

.. code-block:: python
   :linenos:

   A = gammalib.GMatrix(10,20)                # Full 10 x 20 matrix
   B = gammalib.GMatrixSymmetric(10,10)       # Symmetric 10 x 10 matrix
   C = gammalib.GMatrixSparse(1000,10000)     # Sparse 1000 x 10000 matrix

   D = gammalib.GMatrix(0,0)                  # WRONG: empty matrix not allowed
   E = gammalib.GMatrixSymmetric(20,22)       # WRONG: symmetric matrix requested


In the constructor, the first argument specifies the number of rows, the
second the number of columns: ``A(row,column)``. A symmetric matrix needs of
course an equal number of rows and columns. And an empty matrix is not
allowed. All matrix elements are initialised to 0 by the matrix
allocation.

Storage class conversion constructors exist for all three classes
to transform one storage class into another:

**C++**

.. code-block:: cpp
   :linenos:

   GMatrix          full(10,10);
   GMatrixSymmetric symmetric(plain);
   GMatrixSparse    sparse(symmetric);

**Python**

.. code-block:: python
   :linenos:

   full      = gammalib.GMatrix(10,10)
   symmetric = gammalib.GMatrixSymmetric(full)
   sparse    = gammalib.GMatrixSparse(symmetric)

Matrix elements are accessed using

**C++**

.. code-block:: cpp
   :linenos:

   std::cout << B(5,7) << std::endl;          // Matrix access in C++ through () operator

**Python**

.. code-block:: python
   :linenos:

   print(B[5,7])                              # Matrix access in Python through [] operator

Filling of sparse matrix elements using this operator is possible, although
this can be time consuming due to internal memory management. In general,
each fill of a new non-zero element needs to shift all elements
that are located after that element in memory.

To reduce the memory management overhead in the filling of a sparse matrix,
methods have been implemented that allow to fill a matrix column wise:

**C++**

.. code-block:: cpp
   :linenos:

   GMatrixSparse sparse(10,5);
   GVector       column(10);
   column[0] = 1.0;
   column[1] = 2.0;
   column[5] = 8.0;
   sparse.column(0, column);
   sparse.add_to_column(0, column);

**Python**

.. code-block:: cpp
   :linenos:

   sparse    = gammalib.GMatrixSparse(10,5)
   column    = gammalib.GVector(10)
   column[0] = 1.0
   column[1] = 2.0
   column[5] = 8.0
   sparse.column(0, column)
   sparse.add_to_column(0, column)

Line 1 allocates a sparse matrix with 10 rows and 5 columns, line 2
instantiates a vector with 10 elements. In lines 3-5, 3 elements of
the vector are set to specific values, all other elements will default
to 0. In line 6, the elements of the vector are set as the elements
of the first matrix column (column 0). Line 7 differs from line 6 in
that the elements are now not set but added to the existing matrix
elements.

To further reduce the memory management overhead for the column-wise
fill of a sparse matrix, a "fill-stack" has been implemented. The
"fill-stack" is a buffer that implements a queue for columns that are
to be set or added to the matrix. The columns will be stored in this
"fill-stack" in the order they are provided, and only once the 
"fill-stack" is full, or upon request, the "fill-stack" will be flushed
into memory. The "fill-stack" is used as follows:

**C++**

.. code-block:: cpp
   :linenos:

   sparse.stack_init(size, entries);
   ...
   sparse.column(0, column);
   ...
   sparse.stack_flush();
   ...
   sparse.stack_destroy();

The ``stack_init(size, entries)`` method initialises the "fill-stack",
where ``size`` is the size of the allocated memory buffer and ``entries``
is the maximum number of columns that will be held by the buffer.
If ``size`` is set to 0 (the default value), a default ``size`` value of
512 is used. If ``entries`` is set to 0 (the default value), the number of
matrix columns is taken as default ``entries`` value. Note that a too large
number of elements will produce some overhead due to "fill-stack"
management, hence ``entries`` should not exceeed a value of the order of
10-100.

The ``stack_flush()`` method flushes the stack, which is mandatory
before any usage of the matrix. Note that the "fill-stack" **is not
inserted automatically** before any matrix operation, hence manual stack
flushing is needed to make all filled matrix elements available for usage.
The ``stack_destroy()`` method will flush the stack and free all stack
elements. This method should be called once no filling is required anymore.
If ``stack_destroy()`` is called immediately after filling, no call to 
``stack_flush()`` is needed as the ``stack_destroy()`` method flushes the
stack before destroying it. The matrix stack is also destroyed by the
sparse matrix destructor, hence manual stack destruction is not
mandatory.
