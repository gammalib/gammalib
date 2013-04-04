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
class GSparseSymbolic;
class GSparseNumeric;


/***********************************************************************//**
 * @class GMatrixSparse
 *
 * @brief Sparse matrix class interface defintion
 *
 * This class implements a sparse matrix class. The class only stores
 * non-zero elements, which can considerably reduce the memory requirements
 * for large systems that are sparsly filled.
 *
 * The class has been inspired from the code CSparse that can be downloaded
 * from http://www.cise.ufl.edu/research/sparse/CSparse/ 
 * CSparse has been written by Timothy A. Davis to accompany his book
 * "Direct Methods for Sparse Linear Systems". In the CSparse code, a sparse
 * matrix is implemented by the structure cs_sparse. The members of this
 * structure map as follows to the member of the GMatrixSparse class:
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
    virtual GVector        row(const int& row) const;
    virtual void           row(const int& row, const GVector& vector);
    virtual GVector        column(const int& column) const;
    virtual void           column(const int& column, const GVector& vector);
    virtual void           add_to_row(const int& row, const GVector& vector);
    virtual void           add_to_column(const int& column, const GVector& vector);
    virtual void           transpose(void);
    virtual void           invert(void);
    virtual void           negate(void);
    virtual void           abs(void);
    virtual double         fill(void) const;
    virtual double         min(void) const;
    virtual double         max(void) const;
    virtual double         sum(void) const;
    virtual std::string    print(void) const;

    // Other methods
    void    column(const double* values, const int* rows,
                   int number, const int& col);
    void    add_to_column(const double* values, const int* rows,
                          int number, const int& col);
    void    cholesky_decompose(bool compress = true);
    GVector cholesky_solver(const GVector& vector, bool compress = true);
    void    cholesky_invert(bool compress = true);
    void    set_mem_block(const int& block);
    void    stack_init(const int& size = 0, const int& entries = 0);
    int     stack_push_column(const GVector& vector, const int& col);
    int     stack_push_column(const double* values, const int* rows,
                              const int& number, const int& col);
    void    stack_flush(void);
    void    stack_destroy(void);

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
 * @brief Negate matrix elements
 *
 * @return Matrix with negated elements.
 ***************************************************************************/
inline
GMatrixSparse GMatrixSparse::operator-(void) const
{
    GMatrixSparse result = *this;
    result.negate();
    return result;
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


/***********************************************************************//**
 * @brief Return transpose of matrix
 *
 * @param[in] matrix Matrix.
 * @return Transpose of matrix.
 ***************************************************************************/
inline
GMatrixSparse transpose(const GMatrixSparse& matrix)
{
    GMatrixSparse result = matrix;
    result.transpose();
    return result;
}


/***********************************************************************//**
 * @brief Return inverse of matrix
 *
 * @param[in] matrix Matrix.
 * @return Inverse of matrix.
 ***************************************************************************/
inline
GMatrixSparse invert(const GMatrixSparse& matrix)
{
    GMatrixSparse result = matrix;
    result.invert();
    return result;
}


/***********************************************************************//**
 * @brief Return matrix with absolute values of all elements
 *
 * @param[in] matrix Matrix.
 * @return Matrix with elements being the absolute elements of the input
 *         matrix.
 ***************************************************************************/
inline
GMatrixSparse abs(const GMatrixSparse& matrix)
{
    GMatrixSparse result = matrix;
    result.abs();
    return result;
}


/***********************************************************************//**
 * @brief Return Cholesky decomposition of matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] compress Use matrix compression (defaults to true).
 * @return Cholesky decomposition of matrix.
 ***************************************************************************/
inline
GMatrixSparse cholesky_decompose(const GMatrixSparse& matrix, bool compress)
{
    GMatrixSparse result = matrix;
    result.cholesky_decompose(compress);
    return result;
}


/***********************************************************************//**
 * @brief Return inverse matrix using Cholesky decomposition
 *
 * @param[in] matrix Matrix.
 * @param[in] compress Use matrix compression (defaults to true).
 * @return Inverse of matrix.
 ***************************************************************************/
inline
GMatrixSparse cholesky_invert(const GMatrixSparse& matrix, bool compress)
{
    GMatrixSparse result = matrix;
    result.cholesky_invert(compress);
    return result;
}

#endif /* GMATRIXSPARSE_HPP */
