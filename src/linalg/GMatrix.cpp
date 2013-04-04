/***************************************************************************
 *                    GMatrix.cpp - General matrix class                   *
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
 * @file GMatrix.cpp
 * @brief General matrix class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GException.hpp"
#include "GTools.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GMatrixSparse.hpp"
#include "GMatrixSymmetric.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CONSTRUCTOR                          "GMatrix::GMatrix(int&, int&)"
#define G_ACCESS                              "GMatrix::operator(int&, int&)"
#define G_OP_MUL_VEC                           "GMatrix::operator*(GVector&)"
#define G_OP_ADD                              "GMatrix::operator+=(GMatrix&)"
#define G_OP_SUB                              "GMatrix::operator-=(GMatrix&)"
#define G_OP_MAT_MUL                          "GMatrix::operator*=(GMatrix&)"
#define G_EXTRACT_ROW                                    "GMatrix::row(int&)"
#define G_SET_ROW                              "GMatrix::row(int&, GVector&)"
#define G_EXTRACT_COLUMN                              "GMatrix::column(int&)"
#define G_SET_COLUMN                        "GMatrix::column(int&, GVector&)"
#define G_ADD_TO_ROW                    "GMatrix::add_to_row(int&, GVector&)"
#define G_ADD_TO_COLUMN              "GMatrix::add_to_column(int&, GVector&)"
#define G_INVERT                                          "GMatrix::invert()"
#define G_EXTRACT_LOWER                   "GMatrix::extract_lower_triangle()"
#define G_EXTRACT_UPPER                   "GMatrix::extract_upper_triangle()"


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void matrix constructor
 *
 * Constructs an empty matrix (i.e. a matrix without any elements; the number
 * of rows and columns of the matrix will be zero).
 *
 * @todo Verify that the class is save against empty matrix objects.
 ***************************************************************************/
GMatrix::GMatrix(void) : GMatrixBase()
{
    // Initialise class members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Matrix constructor
 *
 * @param[in] rows Number of rows [>0].
 * @param[in] columns Number of columns [>0].
 *
 * @exception GException::empty
 *            Specified number of rows or columns is not valid.
 *
 * Constructs a matrix with the specified number of @p rows and @p columns.
 * The memory for all matrix elements will be allocated and all matrix
 * elements will be initialized to 0.
 ***************************************************************************/
GMatrix::GMatrix(const int& rows, const int& columns) : GMatrixBase()
{
    // Continue only if matrix is valid
    if (rows > 0 && columns > 0) {

        // Initialise class members for clean destruction
        init_members();

        // Allocate matrix memory
        alloc_members(rows, columns);

    }
    else {
        throw GException::empty(G_CONSTRUCTOR);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
GMatrix::GMatrix(const GMatrix& matrix) : GMatrixBase(matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(matrix);

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrixSymmetric to GMatrix storage class convertor
 *
 * @param[in] matrix Symmetric matrix (GMatrixSymmetric).
 *
 * Constructs a GMatrix object from a symmetric matrix object of type
 * GMatrixSymmetric. Since the result is a generic matrix, the constructor
 * will succeed in all cases.
 ***************************************************************************/
GMatrix::GMatrix(const GMatrixSymmetric& matrix) : GMatrixBase(matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(matrix.rows(), matrix.columns());

    // Fill matrix. We benefit here from the symmetry of the matrix and
    // loop only over the lower triangle of the symmetric matrix to perform
    // the fill.
    for (int col = 0; col < m_cols; ++col) {
        for (int row = col; row < m_rows; ++row) {
            double value      = matrix(row, col);
            (*this)(row, col) = value;
            (*this)(col, row) = value;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrixSparse to GMatrix storage class convertor
 *
 * @param[in] matrix Sparse matrix (GMatrixSparse).
 *
 * Constructs a GMatrix object from a sparse matrix object of type
 * GMatrixSparse. Since the result is a generic matrix, the constructor
 * will succeed in all cases.
 ***************************************************************************/
GMatrix::GMatrix(const GMatrixSparse& matrix) : GMatrixBase(matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct matrix
    alloc_members(matrix.rows(), matrix.columns());

    // Fill matrix
    for (int col = 0; col < m_cols; ++col) {
        for (int row = 0; row < m_rows; ++row) {
            (*this)(row, col) = matrix(row, col);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GMatrix::~GMatrix(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] matrix Matrix.
 * @return Matrix.
 ***************************************************************************/
GMatrix& GMatrix::operator=(const GMatrix& matrix)
{
    // Execute only if object is not identical
    if (this != &matrix) {

        // Assign base class members. Note that this method will also perform
        // the allocation of the matrix memory and the copy of the matrix
        // attributes.
        this->GMatrixBase::operator=(matrix);

        // Free derived class members
        free_members();

        // Initialise derived class members
        init_members();

        // Copy derived class members
        copy_members(matrix);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return reference to matrix element
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] column Matrix column [0,...,columns()-1].
 *
 * @exception GException::out_of_range
 *            Row or column index out of range.
 *
 * Returns a reference to the matrix element at @p row and @p column.
 ***************************************************************************/
double& GMatrix::operator()(const int& row, const int& column)
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_ACCESS, row, column, m_rows, m_cols);
    }
    #endif

    // Return element
    return (m_data[m_colstart[column]+row]);
}


/***********************************************************************//**
 * @brief Return reference to matrix element (const version)
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] column Matrix column [0,...,columns()-1].
 *
 * @exception GException::out_of_range
 *            Row or column index out of range.
 *
 * Returns a const reference to the matrix element at @p row and @p column.
 ***************************************************************************/
const double& GMatrix::operator()(const int& row, const int& column) const
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || column < 0 || column >= m_cols) {
      throw GException::out_of_range(G_ACCESS, row, column, m_rows, m_cols);
    }
    #endif

    // Return element
    return (m_data[m_colstart[column]+row]);
}



/***********************************************************************//**
 * @brief Vector multiplication
 *
 * @param[in] vector Vector.
 *
 * @exception GException::matrix_vector_mismatch
 *            Vector length differs from number of columns in matrix.
 *
 * This method performs a vector multiplication of a matrix. The vector
 * multiplication will produce a vector. The matrix multiplication can only
 * be performed when the number of matrix columns is identical to the length
 * of the vector.
 ***************************************************************************/
GVector GMatrix::operator*(const GVector& vector) const
{
    // Raise an exception if the matrix and vector dimensions are not compatible
    if (m_cols != vector.size()) {
        throw GException::matrix_vector_mismatch(G_OP_MUL_VEC, vector.size(),
                                                 m_rows, m_cols);
    }

    // Perform vector multiplication
    GVector result(m_rows);
    for (int row = 0; row < m_rows; ++row) {
        double sum = 0.0;
        for (int col = 0; col < m_cols; ++col) {
            sum += (*this)(row,col) * vector[col];
        }
        result[row] = sum;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Unary matrix addition operator
 *
 * @param[in] matrix Matrix.
 *
 * @exception GException::matrix_mismatch
 *            Incompatible matrix size.
 *
 * This method performs a matrix addition. The operation can only succeed
 * when the dimensions of both matrices are identical.
 ***************************************************************************/
GMatrix& GMatrix::operator+=(const GMatrix& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows || m_cols != matrix.m_cols) {
        throw GException::matrix_mismatch(G_OP_ADD,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Add all matrix elements
    const double* src = matrix.m_data;
    double*       dst = m_data;
    for (int i = 0; i < m_elements; ++i) {
        *dst++ += *src++;
    }

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Unary matrix subtraction operator
 *
 * @param[in] matrix Matrix.
 *
 * @exception GException::matrix_mismatch
 *            Incompatible matrix size.
 *
 * This method performs a matrix addition. The operation can only succeed
 * when the dimensions of both matrices are identical.
 ***************************************************************************/
GMatrix& GMatrix::operator-=(const GMatrix& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows || m_cols != matrix.m_cols) {
        throw GException::matrix_mismatch(G_OP_SUB,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Subtract all matrix elements
    const double* src = matrix.m_data;
    double*       dst = m_data;
    for (int i = 0; i < m_elements; ++i) {
        *dst++ -= *src++;
    }

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Unary matrix multiplication operator
 *
 * @param[in] matrix Matrix.
 *
 * @exception GException::matrix_mismatch
 *            Incompatible matrix size.
 *
 * This method performs a matrix multiplication. The operation can only
 * succeed when the dimensions of both matrices are compatible.
 *
 * In case of rectangular matrices the result matrix does not change and
 * the operation is performed inplace. For the general case the result
 * matrix changes in size and for simplicity a new matrix is allocated to
 * hold the result.
 ***************************************************************************/
GMatrix& GMatrix::operator*=(const GMatrix& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_cols != matrix.m_rows) {
        throw GException::matrix_mismatch(G_OP_MAT_MUL,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Case A: Matrices are rectangular, so perform 'inplace' multiplication
    if (m_rows == m_cols) {
        for (int row = 0; row < m_rows; ++row) {
            GVector v_row = this->row(row);
            for (int col = 0; col < m_cols; ++col) {
                double sum = 0.0;
                for (int i = 0; i < m_cols; ++i) {
                    sum += v_row[i] * matrix(i,col);
                }
                (*this)(row,col) = sum;
            }
        }
    }

    // Case B: Matrices are not rectangular, so we cannot work inplace
    else {

        // Allocate result matrix
        GMatrix result(m_rows, matrix.m_cols);

        // Loop over all elements of result matrix
        for (int row = 0; row < m_rows; ++row) {
            for (int col = 0; col < matrix.m_cols; ++col) {
                double sum = 0.0;
                for (int i = 0; i < m_cols; ++i) {
                    sum += (*this)(row,i) * matrix(i,col);
                }
                result(row,col) = sum;
            }
        }

        // Assign result
        *this = result;

    } // endelse: matrices were not rectangular

    // Return result
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                            Public methods                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear matrix
 ***************************************************************************/
void GMatrix::clear(void)
{
    // Free members
    free_members();

    // Initialise private members
    init_members();
    
    // Return
    return; 
}


/***********************************************************************//**
 * @brief Clone matrix
 *
 * @return Pointer to deep copy of matrix.
 ***************************************************************************/
GMatrix* GMatrix::clone(void) const
{
    // Clone matrix
    return new GMatrix(*this);
}


/***********************************************************************//**
 * @brief Extract row as vector from matrix
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @return Vector of matrix row elements (columns() elements).
 *
 * @exception GException::out_of_range
 *            Invalid row index specified.
 *
 * Extracts one @p row of the matrix into a vector. The vector will contain
 * columns() elements.
 ***************************************************************************/
GVector GMatrix::row(const int& row) const
{
    // Raise an exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_EXTRACT_ROW, row, 0, m_rows-1);
    }
    #endif

    // Create result vector
    GVector result(m_cols);

    // Extract row into vector. Recall that the matrix elements are stored
    // column wise, hence we have to increment the element counter i by
    // the number of rows to go into the next matrix column.
    for (int col = 0, i = row; col < m_cols; ++col, i += m_rows) {
        result[col] = m_data[i];
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Set row in matrix
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] vector Vector of matrix row elements (columns() elements).
 *
 * @exception GException::out_of_range
 *            Invalid row index specified.
 * @exception GException::matrix_vector_mismatch
 *            Vector does not match the matrix dimensions.
 *
 * Sets the elements from a vector as the elements of a matrix row. The
 * length of the vector must be identical to the number of columns in the
 * matrix.
 ***************************************************************************/
void GMatrix::row(const int& row, const GVector& vector)
{
    // Raise an exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_SET_ROW, row, 0, m_rows-1);
    }
    #endif

    // Raise an exception if the matrix and vector dimensions are not
    // compatible
    if (m_cols != vector.size()) {
        throw GException::matrix_vector_mismatch(G_SET_ROW, vector.size(),
                                                 m_rows, m_cols);
    }

    // Set row from vector. Recall that the matrix elements are stored
    // column wise, hence we have to increment the element counter i by
    // the number of rows to go into the next matrix column.
    for (int col = 0, i = row; col < m_cols; ++col, i += m_rows) {
         m_data[i] = vector[col];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract column as vector from matrix
 *
 * @param[in] column Column index [0,...,columns()-1].
 *
 * @exception GException::out_of_range
 *            Invalid row index specified.
 *
 * This method extracts a matrix column into a vector.
 ***************************************************************************/
GVector GMatrix::column(const int& column) const
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_EXTRACT_COLUMN, column, 0, m_cols-1);
    }
    #endif

    // Create result vector
    GVector result(m_rows);

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start of data in matrix
        int i = m_colstart[column];

        // Extract column into vector
        for (int row = 0; row < m_rows; ++row) {
            result[row] = m_data[i++];
        }

    } // endif: matrix had elements

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Set matrix column from vector
 *
 * @param[in] column Column index [0,...,columns()-1].
 * @param[in] vector Vector.
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * Inserts the content of a vector into a matrix column. Any previous
 * content in the matrix column will be overwritted.
 ***************************************************************************/
void GMatrix::column(const int& column, const GVector& vector)
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_SET_COLUMN, column, 0, m_cols-1);
    }
    #endif

    // Raise an exception if the matrix and vector dimensions are not
    // compatible
    if (m_rows != vector.size()) {
        throw GException::matrix_vector_mismatch(G_SET_COLUMN, vector.size(),
                                                 m_rows, m_cols);
    }

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start index of data in matrix
        int i = m_colstart[column];

        // Insert column into vector
        for (int row = 0; row < m_rows; ++row) {
            m_data[i++] = vector[row];
        }

    } // endif: matrix had elements

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add row to matrix elements
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] vector Vector of matrix row elements (columns() elements).
 *
 * @exception GException::out_of_range
 *            Invalid row index specified.
 * @exception GException::matrix_vector_mismatch
 *            Vector does not match the matrix dimensions.
 ***************************************************************************/
void GMatrix::add_to_row(const int& row, const GVector& vector)
{
    // Raise an exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_ADD_TO_ROW, row, 0, m_rows-1);
    }
    #endif

    // Raise an exception if the matrix and vector dimensions are not
    // compatible
    if (m_cols != vector.size()) {
        throw GException::matrix_vector_mismatch(G_ADD_TO_ROW, vector.size(),
                                                 m_rows, m_cols);
    }

    // Add row from vector. Recall that the matrix elements are stored
    // column wise, hence we have to increment the element counter i by
    // the number of rows to go into the next matrix column.
    for (int col = 0, i = row; col < m_cols; ++col, i += m_rows) {
         m_data[i] += vector[col];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add vector column into matrix
 *
 * @param[in] column Column index [0,...,columns()-1].
 * @param[in] vector Vector.
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * Adds the content of a vector to a matrix column.
 ***************************************************************************/
void GMatrix::add_to_column(const int& column, const GVector& vector)
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_ADD_TO_COLUMN, column, 0, m_cols-1);
    }
    #endif

    // Raise an exception if the matrix and vector dimensions are not
    // compatible
    if (m_rows != vector.size()) {
        throw GException::matrix_vector_mismatch(G_ADD_TO_COLUMN, vector.size(),
                                                 m_rows, m_cols);
    }

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start index of data in matrix
        int i = m_colstart[column];

        // Add column into vector
        for (int row = 0; row < m_rows; ++row) {
            m_data[i++] += vector[row];
        }

    } // endif: matrix had elements

    // Return
    return;
}


/***********************************************************************//**
 * @brief Transpose matrix
 *
 * The transpose operation exchanges the number of rows against the number
 * of columns. For a square matrix the exchange is done inplace. Otherwise
 * a copy of the matrix is made.
 ***************************************************************************/
void GMatrix::transpose(void)
{
    // Case A: Matrix is square then simply swap the elements
    if (m_rows == m_cols) {
        double  swap;
        double* ptr_dst;
        double* ptr_src;
        for (int row = 0; row < m_rows; ++row) {
            for (int col = row; col < m_cols; ++col) {
                ptr_dst  = m_data + m_cols*row + col;
                ptr_src  = m_data + m_rows*col + row;
                swap     = *ptr_dst;
                *ptr_dst = *ptr_src;
                *ptr_src = swap;
            }
        }
    }

    // Case B: Non-rectangular transpose. No code optimization has been
    // done so far.
    else {
        GMatrix result(m_cols, m_rows);
        for (int row = 0; row < m_rows; ++row) {
            for (int col = 0; col < m_cols; ++col) {
                result(col, row) = (*this)(row, col);
            }
        }
        *this = result;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invert matrix
 *
 * @exception GException::feature_not_implemented
 *            Feature not yet implemented.
 *
 * @todo Needs to be implemented.
 ***************************************************************************/
void GMatrix::invert(void)
{
    // Throw exception
    throw GException::feature_not_implemented(G_INVERT);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Negate all matrix elements
 *
 * Negates all matrix elements.
 ***************************************************************************/
void GMatrix::negate(void)
{
    // Negate all matrix elements
    double* ptr = m_data;
    for (int i = 0; i < m_elements; ++i, ++ptr) {
        *ptr = -(*ptr);
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Take absolute value of matrix elements
 *
 * Replaces all elements of the matrix by their absolute values.
 ***************************************************************************/
void GMatrix::abs(void)
{
    // Take the absolute value of all matrix elements
    double* ptr = m_data;
    for (int i = 0; i < m_elements; ++i, ++ptr) {
        *ptr = std::abs(*ptr);
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return fill of matrix
 *
 * @return Matrix fill (between 0 and 1).
 *
 * Returns the fill of the matrix. The fill of a matrix is defined as the
 * number non-zero elements devided by the number of total elements. By
 * definiton, the fill is comprised in the interval [0,..,1].
 *
 * The fill of a matrix with zero elements will be set to 0.
 ***************************************************************************/
double GMatrix::fill(void) const
{
    // Initialise result
    double result = 0.0;

    // Continue only if matrix has elements
    if (m_elements > 0) {

        // Determine the number of zero elements
        int zero = 0;
        for (int i = 0; i < m_elements; ++i) {
            if (m_data[i] == 0.0) {
                zero++;
            }
        }

        // Compute fill
        result = 1.0-double(zero)/double(m_elements);

    } // endif: there were elements in the matrix

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Extract lower triangle of matrix
 *
 * @exception GException::matrix_not_square
 *            Matrix is not a square matrix.
 *
 * This method extracts the lower triangle of a matrix into another matrix.
 * (including the diagonal elements).
 * All remaining matrix elements will be zero.
 *
 * Triangular extraction only works for square matrixes.
 ***************************************************************************/
GMatrix GMatrix::extract_lower_triangle(void) const
{
    // Raise an exception if matrix is not square
    if (m_rows != m_cols) {
        throw GException::matrix_not_square(G_EXTRACT_LOWER, m_rows, m_cols);
    }

    // Define result matrix
    GMatrix result(m_rows, m_cols);

    // Extract all elements
    for (int col = 0; col < m_cols; ++col) {
        int i = m_colstart[col] + col;
        for (int row = col; row < m_rows; ++row, ++i) {
            result.m_data[i] = m_data[i];
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Extract upper triangle of matrix
 *
 * @exception GException::matrix_not_square
 *            Matrix is not a square matrix.
 *
 * This method extracts the upper triangle of a matrix into another matrix.
 * (including the diagonal elements).
 * All remaining matrix elements will be zero.
 *
 * Triangular extraction only works for square matrixes.
 ***************************************************************************/
GMatrix GMatrix::extract_upper_triangle(void) const
{
    // Raise an exception if matrix is not squared
    if (m_rows != m_cols) {
        throw GException::matrix_not_square(G_EXTRACT_UPPER, m_rows, m_cols);
    }

    // Define result matrix
    GMatrix result(m_rows, m_cols);

    // Extract all elements
    for (int col = 0; col < m_cols; ++col) {
        int i = m_colstart[col];
        for (int row = 0; row <= col; ++row, ++i) {
            result.m_data[i] = m_data[i];
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Set Euler rotation matrix around x axis
 *
 * @param[in] angle Rotation angle (degrees)
 ***************************************************************************/
void GMatrix::eulerx(const double& angle)
{
    // Free members
    GMatrixBase::free_members();

    // Initialise members
    GMatrixBase::init_members();

    // Construct 3*3 matrix
    alloc_members(3,3);

    // Compute angles
    double cosangle = std::cos(angle * deg2rad);
    double sinangle = std::sin(angle * deg2rad);

    // Set matrix elements
    (*this)(0,0) =       1.0;
    (*this)(0,1) =       0.0;
    (*this)(0,2) =       0.0;
    (*this)(1,0) =       0.0;
    (*this)(1,1) =  cosangle;
    (*this)(1,2) = -sinangle;
    (*this)(2,0) =       0.0;
    (*this)(2,1) =  sinangle;
    (*this)(2,2) =  cosangle;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set Euler rotation matrix around y axis
 *
 * @param[in] angle Rotation angle (degrees)
 ***************************************************************************/
void GMatrix::eulery(const double& angle)
{
    // Free members
    GMatrixBase::free_members();

    // Initialise members
    GMatrixBase::init_members();

    // Construct 3*3 matrix
    alloc_members(3,3);

    // Compute angles
    double cosangle = std::cos(angle * deg2rad);
    double sinangle = std::sin(angle * deg2rad);

    // Set matrix elements
    (*this)(0,0) =  cosangle;
    (*this)(0,1) =       0.0;
    (*this)(0,2) =  sinangle;
    (*this)(1,0) =       0.0;
    (*this)(1,1) =       1.0;
    (*this)(1,2) =       0.0;
    (*this)(2,0) = -sinangle;
    (*this)(2,1) =       0.0;
    (*this)(2,2) =  cosangle;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set Euler rotation matrix around z axis
 *
 * @param[in] angle Rotation angle (degrees)
 ***************************************************************************/
void GMatrix::eulerz(const double& angle)
{
    // Free members
    GMatrixBase::free_members();

    // Initialise members
    GMatrixBase::init_members();

    // Construct 3*3 matrix
    alloc_members(3,3);

    // Compute angles
    double cosangle = std::cos(angle * deg2rad);
    double sinangle = std::sin(angle * deg2rad);

    // Set matrix elements
    (*this)(0,0) =  cosangle;
    (*this)(0,1) = -sinangle;
    (*this)(0,2) =       0.0;
    (*this)(1,0) =  sinangle;
    (*this)(1,1) =  cosangle;
    (*this)(1,2) =       0.0;
    (*this)(2,0) =       0.0;
    (*this)(2,1) =       0.0;
    (*this)(2,2) =       1.0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print matrix
 *
 * @return String containing matrix information.
 ***************************************************************************/
std::string GMatrix::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GMatrix ===");

    // Append information
    result.append("\n"+parformat("Number of rows")+str(m_rows));
    if (m_rowsel != NULL) {
        result.append(" (compressed "+str(m_num_rowsel)+")");
    }
    result.append("\n"+parformat("Number of columns")+str(m_cols));
    if (m_colsel != NULL) {
        result.append(" (compressed "+str(m_num_colsel)+")");
    }
    result.append("\n"+parformat("Number of elements")+str(m_elements));
    result.append("\n"+parformat("Number of allocated cells")+str(m_alloc));

    // Append elements and compression schemes
    result.append(print_elements());
    result.append(print_row_compression());
    result.append(print_col_compression());

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GMatrix::init_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
void GMatrix::copy_members(const GMatrix& matrix)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMatrix::free_members(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate matrix memory
 *
 * @param[in] rows Number of matrix rows.
 * @param[in] columns Number of matrix columns.
 *
 * Allocates memory for the matrix elements. The method assumes that no
 * memory has been allocated for the matrix elements and the column start
 * index array. The method allocates the memory for matrix elements and the
 * column start indices, sets all matrix elements to 0, and sets the column
 * start indices.
 *
 * If either the number of rows or the number of columns is non-positive, the
 * method does nothing.
 ***************************************************************************/
void GMatrix::alloc_members(const int& rows, const int& columns)
{
    // Determine number of elements to store in matrix
    int elements = rows * columns;

    // Continue only if number of elements is positive
    if (elements > 0) {

        // Allocate matrix array and column start index array.
        m_data     = new double[elements];
        m_colstart = new int[columns+1];

        // Store matrix size (logical, storage, allocated)
        m_rows     = rows;
        m_cols     = columns;
        m_elements = elements;
        m_alloc    = elements;

        // Set-up column start indices
        m_colstart[0] = 0;
        for (int col = 1; col <= m_cols; ++col) {
            m_colstart[col] = m_colstart[col-1] + m_rows;
        }
        
        // Initialise matrix elements to 0
        for (int i = 0; i < m_elements; ++i) {
            m_data[i] = 0.0;
        }

    } // endif: number of elements was positive

    // Return
    return;
}
