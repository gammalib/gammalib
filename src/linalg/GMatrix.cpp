/***************************************************************************
 *                   GMatrix.cpp  -  General matrix class                  *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2012 by Juergen Knoedlseder                         *
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
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS1                              "GMatrix::operator(int&,int&)"
#define G_ACCESS2                        "GMatrix::operator(int&,int&) const"
#define G_OP_MUL_VEC                           "GMatrix::operator*(GVector&)"
#define G_OP_ADD                              "GMatrix::operator+=(GMatrix&)"
#define G_OP_SUB                              "GMatrix::operator-=(GMatrix&)"
#define G_OP_MAT_MUL                          "GMatrix::operator*=(GMatrix&)"
#define G_INVERT                                          "GMatrix::invert()"
#define G_ADD_COL                           "GMatrix::add_col(GVector&,int&)"
#define G_EXTRACT_ROW                            "GMatrix::extract_row(int&)"
#define G_EXTRACT_COL                            "GMatrix::extract_col(int&)"
#define G_EXTRACT_LOWER               "GMatrix::extract_lower_triangle(void)"
#define G_EXTRACT_UPPER               "GMatrix::extract_upper_triangle(void)"
#define G_INSERT_COL                     "GMatrix::insert_col(GVector&,int&)"
#define G_ALLOC_MEMBERS                   "GMatrix::alloc_members(int&,int&)"


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void matrix constructor
 *
 * This method will allocate a generic matrix object without any elements.
 * The number of rows and columns of the matrix will be zero.
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
 * @param[in] rows Number of rows.
 * @param[in] cols Number of columns.
 *
 * This method allocates a generic matrix object with the specified number
 * of rows and columns.
 ***************************************************************************/
GMatrix::GMatrix(const int& rows, const int& cols) : GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(rows, cols);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] matrix Matrix.
 *
 * This method copies a matrix object.
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
 * @brief GSymMatrix to GMatrix storage class convertor
 *
 * @param[in] matrix Symmetric matrix (GSymMatrix).
 *
 * This constructor converts a symmetric matrix (of type GSymMatrix) into a
 * generic matrix. As the result is generic, the conversion will succeed in
 * all cases. 
 ***************************************************************************/
GMatrix::GMatrix(const GSymMatrix& matrix) : GMatrixBase(matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(matrix.rows(), matrix.cols());

    // Fill matrix. We benefit here from the symmetry of the matrix and
    // loop only over the lower triangle of the symmetric matrix to perform
    // the fill.
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = col; row < matrix.rows(); ++row) {
            double value      = matrix(row, col);
            (*this)(row, col) = value;
            (*this)(col, row) = value;
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief GSparseMatrix to GMatrix storage class convertor
 *
 * @param[in] matrix Sparse matrix (GSparseMatrix).
 *
 * This constructor converts a sparse matrix (of type GSparseMatrix) into a
 * generic matrix. As the result is generic, the conversion will succeed in
 * all cases. 
 ***************************************************************************/
GMatrix::GMatrix(const GSparseMatrix& matrix) : GMatrixBase(matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct matrix
    alloc_members(matrix.rows(), matrix.cols());

    // Fill matrix
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = 0; row < matrix.rows(); ++row) {
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
 * @brief Access operator
 *
 * @param[in] row Matrix row.
 * @param[in] col Matrix column.
 *
 * @exception GException::out_of_range
 *            Row or column index out of range.
 ***************************************************************************/
double& GMatrix::operator()(const int& row, const int& col)
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_ACCESS1, row, col, m_rows, m_cols);
    }
    #endif

    // Return element
    return m_data[m_colstart[col]+row];
}


/***********************************************************************//**
 * @brief Access operator (const version)
 *
 * @param[in] row Matrix row.
 * @param[in] col Matrix column.
 *
 * @exception GException::out_of_range
 *            Row or column index out of range.
 ***************************************************************************/
const double& GMatrix::operator()(const int& row, const int& col) const
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols) {
      throw GException::out_of_range(G_ACCESS2, row, col, m_rows, m_cols);
    }
    #endif

    // Return element
    return m_data[m_colstart[col]+row];
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
 * be performed when the number of matrix columns is equal to the length of
 * the vector.
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

    // Add matrices using base class method
    addition(matrix);

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

    // Subtract matrices using base class method
    subtraction(matrix);

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
            GVector v_row = extract_row(row);
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
 * @brief Clear instance
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
 * @brief Clone object
 ***************************************************************************/
GMatrix* GMatrix::clone(void) const
{
    // Clone this image
    return new GMatrix(*this);
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
 * @brief Add vector column into matrix
 *
 * @param[in] vector Vector.
 * @param[in] col Column index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * Adds the content of a vector to a matrix column.
 ***************************************************************************/
void GMatrix::add_col(const GVector& vector, const int& col)
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_ADD_COL, col, 0, m_cols-1);
    }
    #endif

    // Raise an exception if the matrix and vector dimensions are not
    // compatible
    if (m_rows != vector.size()) {
        throw GException::matrix_vector_mismatch(G_ADD_COL, vector.size(),
                                                 m_rows, m_cols);
    }

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start index of data in matrix
        int i = m_colstart[col];

        // Add column into vector
        for (int row = 0; row < m_rows; ++row) {
            m_data[i++] += vector[row];
        }

    } // endif: matrix had elements

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert vector column into matrix
 *
 * @param[in] vector Vector.
 * @param[in] col Column index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * Inserts the content of a vector into a matrix column. Any previous
 * content in the matrix column will be overwritted.
 ***************************************************************************/
void GMatrix::insert_col(const GVector& vector, const int& col)
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_INSERT_COL, col, 0, m_cols-1);
    }
    #endif

    // Raise an exception if the matrix and vector dimensions are not
    // compatible
    if (m_rows != vector.size()) {
        throw GException::matrix_vector_mismatch(G_INSERT_COL, vector.size(),
                                                 m_rows, m_cols);
    }

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start index of data in matrix
        int i = m_colstart[col];

        // Insert column into vector
        for (int row = 0; row < m_rows; ++row) {
            m_data[i++] = vector[row];
        }

    } // endif: matrix had elements

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract row as vector from matrix
 *
 * @param[in] row Row to be extracted (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid row index specified.
 *
 * This method extracts a matrix row into a vector.
 ***************************************************************************/
GVector GMatrix::extract_row(const int& row) const
{
    // Raise an exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_EXTRACT_ROW, row, 0, m_rows-1);
    }
    #endif

    // Create result vector
    GVector result(m_cols);

    // Extract row into vector
    for (int col = 0, i = row; col < m_cols; ++col, i+=m_rows) {
        result[col] = m_data[i];
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Extract column as vector from matrix
 *
 * @param[in] col Column to be extracted (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid row index specified.
 *
 * This method extracts a matrix column into a vector.
 ***************************************************************************/
GVector GMatrix::extract_col(const int& col) const
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_EXTRACT_COL, col, 0, m_cols-1);
    }
    #endif

    // Create result vector
    GVector result(m_rows);

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start of data in matrix
        int i = m_colstart[col];

        // Extract column into vector
        for (int row = 0; row < m_rows; ++row) {
            result[row] = m_data[i++];
        }

    } // endif: matrix had elements

    // Return vector
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
 * @brief Return fill of matrix
 *
 * This method returns the fill of a matrix.
 * The fill of a matrix is defined as the number non-zero elements devided
 * by the number of total elements. By definiton, the fill is comprised
 * in the interval [0,..,1].
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
    double cosangle = cos(angle * deg2rad);
    double sinangle = sin(angle * deg2rad);

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
    double cosangle = cos(angle * deg2rad);
    double sinangle = sin(angle * deg2rad);

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
    double cosangle = cos(angle * deg2rad);
    double sinangle = sin(angle * deg2rad);

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
 ***************************************************************************/
std::string GMatrix::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GMatrix ===");
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
 * @brief Allocates matrix memory
 *
 * @param[in] rows Number of rows (>1).
 * @param[in] cols Number of columns (>1).
 *
 * @exception GException::empty
 *            Attempt to allocate zero size matrix
 *
 * This method is the main constructor code that allocates and initialises
 * memory for matrix elements. The method assumes that no memory has been
 * allocated for the matrix elements and the column start index array.
 * The method allocates the memory for matrix elements and the column start
 * indices, sets all matrix elements to 0.0, and sets the column start
 * indices.
 ***************************************************************************/
void GMatrix::alloc_members(const int& rows, const int& cols)
{
    // Determine number of elements to store in matrix
    int elements = rows*cols;

    // Throw exception if requested matrix size is zero
    if (elements == 0) {
        throw GException::empty(G_ALLOC_MEMBERS);
    }

    // Allocate matrix array and column start index array.
    m_data     = new double[elements];
    m_colstart = new int[cols+1];

    // Store matrix size (logical, storage, allocated)
    m_rows     = rows;
    m_cols     = cols;
    m_elements = elements;
    m_alloc    = elements;

    // Set-up column start indices
    m_colstart[0] = 0;
    for (int col = 1; col <= m_cols; ++col) {
        m_colstart[col] = m_colstart[col-1] + m_rows;
    }
        
    // Initialise matrix elements to 0.0
    for (int i = 0; i < m_elements; ++i) {
        m_data[i] = 0.0;
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           Friend functions                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Return matrix with absolute values of all elements
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
GMatrix abs(const GMatrix& matrix)
{
    // Define result matrix
    GMatrix result = matrix;

    // Convert all elements to absolute values  
    for (int i = 0; i < result.m_elements; ++i) {
        result.m_data[i] = std::abs(result.m_data[i]);
    }

    // Return result
    return result;
}
