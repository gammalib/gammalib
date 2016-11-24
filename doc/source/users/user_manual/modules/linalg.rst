.. _sec_linalg:

Linear algebra
--------------

Overview
~~~~~~~~

The following figure presents an overview over the C++ classes of
the linear algebra module and their relations.

.. _fig_uml_linalg:

.. figure:: uml_linalg.png
   :width: 55%
   :align: center

   *Linear algebra module*

The linear algebra module provides classes for vector and matrix 
manipulation. The :doxy:`GVector` class implements a double precision
floating point vector, the classes :doxy:`GMatrix`, :doxy:`GMatrixSymmetric`
and :doxy:`GMatrixSparse` implement a matrix of double precision
floating point values for a generic matrix, a symmetric matrix and
a sparse matrix, respectively. All matrix classes derive from the
abstract :doxy:`GMatrixBase` base class that provides common functionalities
to all matrix classes.

The classes :doxy:`GSparseNumeric` and :doxy:`GSparseSymbolic` are classes
that are used by :doxy:`GMatrixSparse` but these classes are not
exposed to the outside world (i.e. the class definitions are not
part of the GammaLib interface).


Matrix Storage classes
~~~~~~~~~~~~~~~~~~~~~~

Three matrix storage classes are implemented:

* :doxy:`GMatrix` implements a plain matrix storage class which stores
  all elements of the matrix in memory;
* :doxy:`GMatrixSymmetric` implements a symmetric matrix storage class
  which stores only a triangle of matrix elements, imposing thus
  strict matrix symmetry;
* :doxy:`GMatrixSparse` implements a sparse matrix storage class which
  stores only non-zero matrix elements in a column-wise organisation.

Storage class conversion constructors exist for all three classes
to transform one storage class into another:

.. code-block:: cpp
   :linenos:

    GMatrix          plain(10,10);
    GMatrixSymmetric symmetric(plain);
    GMatrixSparse    sparse(symmetric);

In the above example, a plain matrix is instantiated in line 1, the
plain matrix is converted into a symmetric matrix in line 2, and the
symmetric matrix is converted into a sparse matrix in line 3.
Additional complementary storage class conversion constructors exist,
but conversion to a symmetric matrix is of course only possible if the
matrix is indeed symmetric.

Matrix elements are accessed using the ``operator()``. Filling of
sparse matrix elements using this operator is possible, although this
can be time consuming due to internal memory management. In general,
each fill of a new non-zero element needs to shift all elements
that are located after that element in memory.

To reduce the memory management overhead in the filling of a sparse 
matrix, methods have been implemented that allow to fill a
matrix column wise:

.. code-block:: cpp
   :linenos:

    GMatrixSparse sparse(10,5);
    GVector       column(10);
    column[0] = 1.0;
    column[1] = 2.0;
    column[5] = 8.0;
    sparse.column(0, column);
    sparse.add_to_column(0, column);

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

.. code-block:: cpp

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
