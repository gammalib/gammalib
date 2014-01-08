/***************************************************************************
 *                  GMatrixSparse.hpp - Sparse matrix class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2013 by Juergen Knoedlseder                         *
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
 * @file GMatrixSparse.hpp
 * @brief Sparse matrix class definition
 * @author Juergen Knoedlseder
 */

#ifndef GMATRIXSPARSE_HPP
#define GMATRIXSPARSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GMatrixBase.hpp"

/* __ Definitions ________________________________________________________ */
#define G_SPARSE_MATRIX_DEFAULT_MEM_BLOCK      512  // Default Memory block size
#define G_SPARSE_MATRIX_DEFAULT_STACK_ENTRIES 1000  // Max # of stack entries
#define G_SPARSE_MATRIX_DEFAULT_STACK_SIZE     512  // Max # of stack elements

/* __ Forward declarations _______________________________________________ */
class GMatrix;
class GMatrixSymmetric;
class GMatrixSparse;
class GVector;
class GSparseSymbolic;
class GSparseNumeric;


/***********************************************************************//**
 * @class GMatrixSparse
 *
 * @brief Sparse matrix class interface definition
 *
 * This class implements a sparse matrix. The class only stores non-zero
 * elements, which can considerably reduce the memory requirements for large
 * systems that are sparsely filled. Also the computations will be faster
 * because only non-zero elements will be used in the operations.
 *
 * For a description of the common matrix methods, please refer to the
 * GMatrixBase class.
 *
 * The class has been inspired from the code CSparse that can be downloaded
 * from http://www.cise.ufl.edu/research/sparse/CSparse/ which has been
 * written by Timothy A. Davis to accompany his book "Direct Methods for
 * Sparse Linear Systems". In the CSparse code, a sparse matrix is
 * implemented by the structure cs_sparse. The members of this
 * structure map as follows to the members of the GMatrixSparse class:
 *
 * <table border> 
 *  <tr> 
 *     <td><b>CSparse</b></td> 
 *     <td><b>GMatrixSparse</b></td> 
 *     <td><b>Dimension</b></td> 
 *     <td><b>Description</b></td> 
 *  </tr>
 *  <tr>
 *     <td>nzmax</td>
 *     <td>m_elements</td>
 *     <td>n.a.</td>
 *     <td>Maximum number of entries</td>
 *  </tr>
 *  <tr>
 *     <td>m</td>
 *     <td>m_rows</td>
 *     <td>n.a.</td>
 *     <td>Number of rows</td>
 *  </tr>
 *  <tr>
 *     <td>n</td>
 *     <td>m_cols</td>
 *     <td>n.a.</td>
 *     <td>Number of columns</td>
 *  </tr>
 *  <tr>
 *     <td>*p</td>
 *     <td>*m_colstart</td>
 *     <td>m_cols+1</td>
 *     <td>Column pointers (index of first element of each column in 
 *         m_data array)</td>
 *  </tr>
 *  <tr>
 *     <td>*i</td>
 *     <td>*m_rowinx</td>
 *     <td>m_elements</td>
 *     <td>Row indices for all non-zero elements</td>
 *  </tr>
 *  <tr>
 *     <td>*x</td>
 *     <td>*m_data</td>
 *     <td>m_elements</td>
 *     <td>Numerical values of all non-zero elements</td>
 *  </tr>
 * </table>
 *
 * Matrix allocation is done using the constructors
 *
 *     GMatrixSparse sparsematrix(rows, columns);
 *     GMatrixSparse sparsematrix(rows, columns, elements);
 *     GMatrixSparse sparsematrix(matrix);
 *     GMatrixSparse sparsematrix(sparsematrix);
 *     GMatrixSparse sparsematrix(symmetricmatrix);
 *
 * where @p rows and @p columns specify the number of rows and columns of
 * the matrix and @p elements is the size of the memory that should be
 * pre-allocated for matrix storage (optional). Storage conversion
 * constructors exist that allow allocating a sparse matrix by copying from
 * a general matrix of type GMatrix and a symmetric matrix of type
 * GMatrixSymmetric.
 *
 * Because the location of a specific matrix element is not necessarily known
 * in advance, the class implements a so called pending element for matrix
 * element setting. A pending element is an element that is scheduled for
 * insertion in the matrix, but which has not yet been inserted. The matrix
 * element access operators () and at() return the reference to this pending
 * element if a memory location for a given row and column does not yet
 * exist. The protected method fill_pending() is then called before all
 * operations that are done on the matrix to insert the element before the
 * operation is executed.
 *
 * Another mechanism that has been implemented to speed up the fill of the
 * sparse matrix is a "fill stack" that is used for insertion or addition
 * of matrix columns by the
 *
 *     matrix.column(index, vector);
 *     matrix.column(index, values, rows, number);
 *     matrix.add_to_column(index, vector);
 *     matrix.add_to_column(index, values, rows, number);
 *
 * methods. The fill stack is a buffer that implements a queue for columns
 * that are to be inserted into the matrix. This is more efficient than
 * insertion of elements one by one, as every insertion of an element may
 * require to shift all elements back by one memory location.
 *
 * The fill stack is used as follows:
 *
 *     matrix.stack_init(size, entries);
 *     ...
 *     matrix.column(index, vector); // (or one of the other column methods)
 *     ...
 *     matrix.stack_flush();
 *     ...
 *     matrix.stack_destroy();
 *
 * The stack_init(size, entries) method initialises the fill stack, where
 * @p size is the size of the allocated memory buffer and @p entries is the
 * maximum number of entries (i.e. columns) that will be held by the buffer.
 * If @p size is set to 0 (the default value), a default @p size value of
 * 512 is used. If @p entries is set to 0 (the default value), the number of
 * matrix columns is taken as default @p entries value. Note that a too large
 * number of elements will produce some overhead due to fill stack
 * management, hence @p entries should not exceed a value of the order of
 * 10-100.
 *
 * The stack_flush() method flushes the stack, which is mandatory before any
 * usage of the matrix. Note that the fill stack IS NOT INSERTED AUTOMATICALLY
 * before any matrix operation, hence manual stack flushing is needed to make
 * all filled matrix elements available for usage. The stack_destroy() method
 * will flush the stack and free all stack elements. This method
 * should be called once no filling is required anymore. If stack_destroy()
 * is called immediately after filling, no call to stack_flush() is needed as
 * the stack_destroy() method flushes the stack before destroying it. The
 * matrix stack is also destroyed by the sparse matrix destructor, hence
 * manual stack destruction is not mandatory.
 * 
 * Except for *m_rowinx which is implemented on the level of GMatrixSparse,
 * all other members are implemented by the base class GMatrixBase.
 ***************************************************************************/
class GMatrixSparse : public GMatrixBase {

    // Friend classes
    friend class GSparseSymbolic;
    friend class GSparseNumeric;

    // Some friend functions that we should not expose ... but I don't
    // know what to do with them as they are also needed by GSparseSymbolic,
    // yet they have to access low-level stuff from the matrix class, hence
    // they need to be friends ...
    friend GMatrixSparse cs_symperm(const GMatrixSparse& matrix, const int* pinv);
    friend GMatrixSparse cs_transpose(const GMatrixSparse& matrix, int values);
    friend double        cs_cumsum(int* p, int* c, int n);

public:
    // Constructors and destructors
    GMatrixSparse(void);
    explicit GMatrixSparse(const int& rows, const int& columns,
                           const int& elements = 0);
    GMatrixSparse(const GMatrix& matrix);
    GMatrixSparse(const GMatrixSparse& matrix);
    GMatrixSparse(const GMatrixSymmetric& matrix);
    virtual ~GMatrixSparse(void);

    // Implemented pure virtual base class operators
    virtual double&        operator()(const int& row, const int& column);
    virtual const double&  operator()(const int& row, const int& column) const;
    virtual GVector        operator*(const GVector& vector) const;

    // Overloaded virtual base class operators
    virtual bool           operator==(const GMatrixSparse& matrix) const;
    virtual bool           operator!=(const GMatrixSparse& matrix) const;

    // Other operators
    virtual GMatrixSparse& operator=(const GMatrixSparse& matrix);
    virtual GMatrixSparse& operator=(const double& value);
    virtual GMatrixSparse  operator+(const GMatrixSparse& matrix) const;
    virtual GMatrixSparse  operator-(const GMatrixSparse& matrix) const;
    virtual GMatrixSparse  operator*(const GMatrixSparse& matrix) const;
    virtual GMatrixSparse  operator-(void) const;
    virtual GMatrixSparse& operator+=(const GMatrixSparse& matrix);
    virtual GMatrixSparse& operator-=(const GMatrixSparse& matrix);
    virtual GMatrixSparse& operator*=(const GMatrixSparse& matrix);
    virtual GMatrixSparse& operator*=(const double& scalar);
    virtual GMatrixSparse& operator/=(const double& scalar);

    // Implemented pure virtual base class methods
    virtual void           clear(void);
    virtual GMatrixSparse* clone(void) const;
    virtual double&        at(const int& row, const int& column);
    virtual const double&  at(const int& row, const int& column) const;
    virtual GVector        row(const int& row) const;
    virtual void           row(const int& row, const GVector& vector);
    virtual GVector        column(const int& column) const;
    virtual void           column(const int& column, const GVector& vector);
    virtual void           add_to_row(const int& row, const GVector& vector);
    virtual void           add_to_column(const int& column, const GVector& vector);
    virtual double         fill(void) const;
    virtual double         min(void) const;
    virtual double         max(void) const;
    virtual double         sum(void) const;
    virtual std::string    print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void          column(const int& column, const double* values,
                         const int* rows, int number);
    void          add_to_column(const int& column, const double* values,
                                const int* rows, int number);
    GMatrixSparse transpose(void) const;
    GMatrixSparse invert(void) const;
    GVector       solve(const GVector& vector) const;
    GMatrixSparse abs(void) const;
    GMatrixSparse cholesky_decompose(const bool& compress = true) const;
    GVector       cholesky_solver(const GVector& vector, const bool& compress = true) const;
    GMatrixSparse cholesky_invert(const bool& compress = true) const;
    void          set_mem_block(const int& block);
    void          stack_init(const int& size = 0, const int& entries = 0);
    int           stack_push_column(const GVector& vector, const int& col);
    int           stack_push_column(const double* values, const int* rows,
                                    const int& number, const int& col);
    void          stack_flush(void);
    void          stack_destroy(void);

private:
    // Private methods
    void init_members(void);
    void copy_members(const GMatrixSparse& m);
    void free_members(void);
    void alloc_members(const int& rows, const int& cols, const int& elements = 0);
    void init_stack_members(void);
    void free_stack_members(void);
    int  get_index(const int& row, const int& column) const;
    void fill_pending(void);
    void alloc_elements(int start, const int& num);
    void free_elements(const int& start, const int& num);
    void remove_zero_row_col(void);
    void insert_zero_row_col(const int& rows, const int& cols);
    void mix_column_prepare(const int* src1_row, int src1_num,
                            const int* src2_row, int src2_num,
                            int* num_1, int* num_2, int* num_mix);
    void mix_column(const double* src1_data, const int* src1_row, int src1_num,
                    const double* src2_data, const int* src2_row, int src2_num,
                    double* dst_data, int* dst_row, int* dst_num);

    // Private data members
    int*             m_rowinx;       //!< Row-indices of all elements
    double           m_zero;         //!< The zero element (needed for data access)
    double           m_fill_val;     //!< Element to be filled
    int              m_fill_row;     //!< Row into which element needs to be filled
    int              m_fill_col;     //!< Column into which element needs to be filled
    int              m_mem_block;    //!< Memory block to be allocated at once
    GSparseSymbolic* m_symbolic;     //!< Holds GSparseSymbolic object after decomposition
    GSparseNumeric*  m_numeric;      //!< Holds GSparseNumeric object after decomposition

    // Fill-stack
    int     m_stack_max_entries;  //!< Maximum number of entries in the stack
    int     m_stack_size;         //!< Maximum number of elements in the stack
    int     m_stack_entries;      //!< Number of entries in the stack
    int*    m_stack_colinx;       //!< Column index for each entry [m_stack_entries]
    int*    m_stack_start;        //!< Start in stack for each entry [m_stack_entries+1]
    double* m_stack_data;         //!< Stack data [m_stack_size]
    int*    m_stack_rowinx;       //!< Stack row indices [m_stack_size]
    int*    m_stack_work;         //!< Stack flush integer working array [m_cols]
    int*    m_stack_rows;         //!< Stack push integer working array [m_cols]
    double* m_stack_values;       //!< Stack push double buffer [m_cols]
};


/***********************************************************************//**
 * @brief Set memory block size
 *
 * @param[in] block Memory block size.
 *
 * Sets the size of the memory block that will be allocated at once.
 ***************************************************************************/
inline
void GMatrixSparse::set_mem_block(const int& block)
{
    m_mem_block = (block > 0) ? block : 1;
    return;
}


/***********************************************************************//**
 * @brief Binary matrix addition
 *
 * @param[in] matrix Matrix.
 * @return Result of matrix addition.
 *
 * Returns the sum of two matrices. The method makes use of the unary
 * addition operator.
 ***************************************************************************/
inline
GMatrixSparse GMatrixSparse::operator+(const GMatrixSparse& matrix) const
{
    GMatrixSparse result = *this;
    result += matrix;
    return result;
}


/***********************************************************************//**
 * @brief Binary matrix subtraction
 *
 * @param[in] matrix Matrix.
 * @return Result of matrix subtraction.
 *
 * Returns the difference between two matrices. The method makes use of the
 * unary subtraction operator.
 ***************************************************************************/
inline
GMatrixSparse GMatrixSparse::operator-(const GMatrixSparse& matrix) const
{
    GMatrixSparse result = *this;
    result -= matrix;
    return result;
}


/***********************************************************************//**
 * @brief Binary matrix multiplication
 *
 * @param[in] matrix Matrix.
 * @return Result of matrix multiplication.
 *
 * Returns the product of two matrices. The method makes use of the unary
 * multiplication operator.
 ***************************************************************************/
inline
GMatrixSparse GMatrixSparse::operator*(const GMatrixSparse& matrix) const
{
    GMatrixSparse result = *this;
    result *= matrix;
    return result;
}


/***********************************************************************//**
 * @brief Scale matrix elements
 *
 * @param[in] scalar Scale factor.
 * @return Matrix with elements multiplied by @p scalar.
 *
 * Returns a matrix where all elements have been multiplied by the specified
 * @p scalar value.
 ***************************************************************************/
inline
GMatrixSparse& GMatrixSparse::operator*=(const double& scalar)
{
    fill_pending();
    scale_elements(scalar);
    return *this;
}


/***********************************************************************//**
 * @brief Divide matrix elements
 *
 * @param[in] scalar Scalar.
 * @return Matrix with elements divided by @p scalar.
 *
 * Returns a matrix where all elements have been divided by the specified
 * @p scalar value.
 ***************************************************************************/
inline
GMatrixSparse& GMatrixSparse::operator/=(const double& scalar)
{
    *this *= 1.0/scalar;
    return *this;
}


/***********************************************************************//**
 * @brief Multiply matrix by scalar
 *
 * @param[in] matrix Matrix.
 * @param[in] scalar Scalar.
 * @return Matrix divided by @p scalar.
 ***************************************************************************/
inline
GMatrixSparse operator*(const GMatrixSparse& matrix, const double& scalar)
{
    GMatrixSparse result = matrix;
    //result.fill_pending();
    result *= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Multiply matrix by scalar
 *
 * @param[in] scalar Scalar.
 * @param[in] matrix Matrix.
 * @return Matrix divided by @p scalar.
 ***************************************************************************/
inline
GMatrixSparse operator*(const double& scalar, const GMatrixSparse& matrix)
{
    GMatrixSparse result = matrix;
    //result.fill_pending();
    result *= scalar;
    return result;
}


/***********************************************************************//**
 * @brief Divide matrix by scalar
 *
 * @param[in] matrix Matrix.
 * @param[in] scalar Scalar.
 * @return Matrix divided by @p scalar.
 ***************************************************************************/
inline
GMatrixSparse operator/(const GMatrixSparse& matrix, const double& scalar)
{
    GMatrixSparse result = matrix;
    //result.fill_pending();
    result /= scalar;
    return result;
}

#endif /* GMATRIXSPARSE_HPP */
