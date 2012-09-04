/***************************************************************************
 *                       GMatrix.cpp  -  Matrix class                      *
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
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GSymMatrix.hpp"
#include "GSparseMatrix.hpp"
#include "GTools.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS1                                "GMatrix::operator(int,int)"
#define G_ACCESS2                          "GMatrix::operator(int,int) const"
#define G_OP_MUL_VEC                           "GMatrix::operator* (GVector)"
#define G_OP_ADD                              "GMatrix::operator+= (GMatrix)"
#define G_OP_SUB                              "GMatrix::operator-= (GMatrix)"
#define G_OP_MAT_MUL                          "GMatrix::operator*= (GMatrix)"
#define G_ADD_COL                            "GMatrix::add_col(GVector, int)"
#define G_EXTRACT_ROW                             "GMatrix::extract_row(int)"
#define G_EXTRACT_COL                             "GMatrix::extract_col(int)"
#define G_EXTRACT_LOWER               "GMatrix::extract_lower_triangle(void)"
#define G_EXTRACT_UPPER               "GMatrix::extract_upper_triangle(void)"
#define G_INSERT_COL                      "GMatrix::insert_col(GVector, int)"
#define G_CONSTRUCTOR                        "GMatrix::constructor(int, int)"


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void matrix constructor (contains no elements)
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
 * @param[in] rows Number of rows in matrix
 * @param[in] cols Number of columns in matrix
 *
 * First calls the void base class constructor that initialises all base
 * class members, then call the GMatrix initialisation method and finally
 * construct the object.
 ***************************************************************************/
GMatrix::GMatrix(int rows, int cols) : GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Construct full matrix
    constructor(rows, cols);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] matrix General matrix (GMatrix).
 *
 * First calls the void base class copy constructor that copys all base
 * class members, then call the GMatrix member copy method.
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
 ***************************************************************************/
GMatrix::GMatrix(const GSymMatrix& matrix) : GMatrixBase(matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct matrix
    constructor(matrix.rows(), matrix.cols());

    // Fill matrix
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = col; row < matrix.rows(); ++row) {
            (*this)(row, col) = matrix(row, col);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief GSparseMatrix to GMatrix storage class convertor
 *
 * @param[in] matrix Sparse matrix (GSparseMatrix).
 ***************************************************************************/
GMatrix::GMatrix(const GSparseMatrix& matrix) : GMatrixBase(matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Construct matrix
    constructor(matrix.rows(), matrix.cols());

    // Fill matrix
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = col; row < matrix.rows(); ++row) {
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
 =                             GMatrix operators                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] m GMatrix instance to be assigned
 ***************************************************************************/
GMatrix& GMatrix::operator= (const GMatrix& m)
{
    // Execute only if object is not identical
    if (this != &m) {

        // Copy base class members
        this->GMatrixBase::operator=(m);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(m);

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
double& GMatrix::operator() (int row, int col)
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
        throw GException::out_of_range(G_ACCESS1, row, col, m_rows, m_cols);
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
const double& GMatrix::operator() (int row, int col) const
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
      throw GException::out_of_range(G_ACCESS2, row, col, m_rows, m_cols);
    #endif

    // Return element
    return m_data[m_colstart[col]+row];
}



/***********************************************************************//**
 * @brief Vector multiplication
 *
 * @param[in] v GVector with which object is to be multiplied
 ***************************************************************************/
GVector GMatrix::operator* (const GVector& v) const
{
    // Raise an exception if the matrix and vector dimensions are not compatible
    if (m_cols != v.size())
        throw GException::matrix_vector_mismatch(G_OP_MUL_VEC, v.size(),
                                                 m_rows, m_cols);

    // Perform vector multiplication
    GVector result(m_rows);
    for (int row = 0; row < m_rows; ++row) {
        double sum = 0.0;
        for (int col = 0; col < m_cols; ++col)
            sum += (*this)(row,col) * v[col];
        result[row] = sum;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Unary matrix addition operator
 *
 * @param[in] m GMatrix to be added
 ***************************************************************************/
GMatrix& GMatrix::operator+= (const GMatrix& m)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != m.m_rows || m_cols != m.m_cols)
        throw GException::matrix_mismatch(G_OP_ADD, m_rows, m_cols,
                                          m.m_rows, m.m_cols);

    // Add matrices
    addition(m);

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Unary matrix subtraction operator
 *
 * @param[in] m GMatrix to be subtracted
 ***************************************************************************/
GMatrix& GMatrix::operator-= (const GMatrix& m)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != m.m_rows || m_cols != m.m_cols)
        throw GException::matrix_mismatch(G_OP_SUB, m_rows, m_cols,
                                          m.m_rows, m.m_cols);

    // Subtract matrices
    subtraction(m);

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Matrix multiplication operator
 *
 * @param[in] m GMatrix to be multiplied
 *
 * In case of rectangular matrices the result matrix does not change and
 * the operations is performed inplace. For the general case the result
 * matrix changes in size and for simplicity a new matrix is allocated to
 * hold the result.
 ***************************************************************************/
GMatrix& GMatrix::operator*= (const GMatrix& m)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_cols != m.m_rows)
        throw GException::matrix_mismatch(G_OP_MAT_MUL, m_rows, m_cols,
                                          m.m_rows, m.m_cols);

    // Case A: Matrices are rectangular, so perform 'inplace' multiplication
    if (m_rows == m_cols) {
        for (int row = 0; row < m_rows; ++row) {
            GVector v_row = extract_row(row);
            for (int col = 0; col < m_cols; ++col) {
                double sum = 0.0;
                for (int i = 0; i < m_cols; ++i)
                    sum += v_row[i] * m(i,col);
                (*this)(row,col) = sum;
            }
        }
    }

    // Case B: Matrices are not rectangular, so we cannot work inplace
    else {

        // Allocate result matrix
        GMatrix result(m_rows, m.m_cols);

        // Loop over all elements of result matrix
        for (int row = 0; row < m_rows; ++row) {
            for (int col = 0; col < m.m_cols; ++col) {
                double sum = 0.0;
                for (int i = 0; i < m_cols; ++i)
                    sum += (*this)(row,i) * m(i,col);
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
 =                            GMatrix methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Transpose matrix
 *
 * The transpose operation exchanges the number of rows against the number
 * of columns. For a square matrix the exchange is done inplace. Otherwise
 * a copy of the matrix is made.
 ***************************************************************************/
void GMatrix::transpose(void)
{
    // Case A: matrix is square then simply swap the elements
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

    // Case B: Dummy for non-rectangular transpose
    else {
        GMatrix result(m_cols, m_rows);
        for (int row = 0; row < m_rows; ++row) {
            for (int col = 0; col < m_cols; ++col)
                result(col, row) = (*this)(row, col);
        }
        *this = result;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Invert matrix
 *
 * TODO: Needs to be implemented, is just a dummy for now.
 ***************************************************************************/
void GMatrix::invert(void)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Add vector column into matrix
 *
 * @param[in] v GVector to be added into column
 * @param[in] col Column into which vector should be added (starting from 0)
 ***************************************************************************/
void GMatrix::add_col(const GVector& v, int col)
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col >= m_cols)
        throw GException::out_of_range(G_ADD_COL, 0, col, m_rows, m_cols);
    #endif

    // Raise an exception if the matrix and vector dimensions are not compatible
    if (m_rows != v.size())
        throw GException::matrix_vector_mismatch(G_ADD_COL, v.size(),
                                                 m_rows, m_cols);

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start index of data in matrix
        int i = m_colstart[col];

        // Add column into vector
        for (int row = 0; row < m_rows; ++row)
            m_data[i++] += v[row];

    } // endif: matrix had elements

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert vector column into matrix
 *
 * @param[in] v Vector to be inserted into matrix.
 * @param[in] col Column at which vector is to be inserted (starting from 0).
 ***************************************************************************/
void GMatrix::insert_col(const GVector& v, int col)
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col >= m_cols)
        throw GException::out_of_range(G_INSERT_COL, 0, col, m_rows, m_cols);
    #endif

    // Raise an exception if the matrix and vector dimensions are not compatible
    if (m_rows != v.size())
         throw GException::matrix_vector_mismatch(G_INSERT_COL, v.size(),
                                                  m_rows, m_cols);

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start index of data in matrix
        int i = m_colstart[col];

        // Insert column into vector
        for (int row = 0; row < m_rows; ++row)
            m_data[i++] = v[row];

    } // endif: matrix had elements

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract row as vector from matrix
 *
 * @param[in] row Row to be extracted (starting from 0)
 ***************************************************************************/
GVector GMatrix::extract_row(int row) const
{
    // Raise an exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row >= m_rows)
        throw GException::out_of_range(G_EXTRACT_ROW, row, 0, m_rows, m_cols);
    #endif

    // Create result vector
    GVector result(m_cols);

    // Extract row into vector
    for (int col = 0, i = row; col < m_cols; ++col, i+=m_rows)
        result[col] = m_data[i];

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Extract column as vector from matrix
 *
 * @param[in] col Column to be extracted (starting from 0)
 ***************************************************************************/
GVector GMatrix::extract_col(int col) const
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col >= m_cols)
        throw GException::out_of_range(G_EXTRACT_COL, 0, col, m_rows, m_cols);
    #endif

    // Create result vector
    GVector result(m_rows);

    // Continue only if we have elements
    if (m_elements > 0) {

        // Get start of data in matrix
        int i = m_colstart[col];

        // Extract column into vector
        for (int row = 0; row < m_rows; ++row)
            result[row] = m_data[i++];

    } // endif: matrix had elements

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Extract lower triangle of matrix as matrix
 *
 * Triangular extraction only works for square matrixes.
 ***************************************************************************/
GMatrix GMatrix::extract_lower_triangle(void) const
{
    // Raise an exception if matrix is not squared
    if (m_rows != m_cols)
        throw GException::matrix_not_square(G_EXTRACT_LOWER, m_rows, m_cols);

    // Define result matrix
    GMatrix result(m_rows, m_cols);

    // Extract all elements
    for (int col = 0; col < m_cols; ++col) {
        int i = m_colstart[col] + col;
            for (int row = col; row < m_rows; ++row, ++i)
                result.m_data[i] = m_data[i];
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Extract upper triangle of matrix as matrix
 *
 * Triangular extraction only works for square matrixes.
 ***************************************************************************/
GMatrix GMatrix::extract_upper_triangle(void) const
{
    // Raise an exception if matrix is not squared
    if (m_rows != m_cols)
        throw GException::matrix_not_square(G_EXTRACT_UPPER, m_rows, m_cols);

    // Define result matrix
    GMatrix result(m_rows, m_cols);

    // Extract all elements
    for (int col = 0; col < m_cols; ++col) {
        int i = m_colstart[col];
        for (int row = 0; row <= col; ++row, ++i)
            result.m_data[i] = m_data[i];
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return fill of matrix
 *
 * The fill of a matrix is defined as the number non-zero elements devided
 * by the number of total elements. By definiton, the fill is comprised
 * in the interval [0,..,1]. The fill of an undefined matrix is defined to
 * be 0.
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
            if (m_data[i] == 0.0)
                zero++;
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
    constructor(3,3);

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
    constructor(3,3);

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
    constructor(3,3);

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


/*==========================================================================
 =                                                                         =
 =                          GMatrix private methods                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Allocate matrix
 *
 * @param[in] rows Number of rows to be allocated.
 * @param[in] cols Number of columns to be allocated.
 *
 * This is the main constructor code that allocates and initialises memory
 * for matrix elements.
 ***************************************************************************/
void GMatrix::constructor(int rows, int cols)
{
    // Determine number of elements to store in matrix
    int elements = rows*cols;

    // Continue only if we have elements to store
    if (elements > 0) {

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
        for (int col = 1; col <= m_cols; ++col)
            m_colstart[col] = m_colstart[col-1] + m_rows;

        // Initialise matrix elements to 0.0
        for (int i = 0; i < m_elements; ++i)
            m_data[i] = 0.0;

    } // endif: we had elements

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             GMatrix friends                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] m Matrix to put in output stream.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GMatrix& m)
{
    // Put header in stream
    os << "=== GMatrix ===" << std::endl;
    if (m.m_rowsel != NULL)
        os << " Number of rows ............: " << m.m_rows
           << " (compressed " << m.m_num_rowsel << ")" << std::endl;
    else
        os << " Number of rows ............: " << m.m_rows << std::endl;
    if (m.m_colsel != NULL)
        os << " Number of columns .........: " << m.m_cols
           << " (compressed " << m.m_num_colsel << ")" << std::endl;
    else
        os << " Number of columns .........: " << m.m_cols << std::endl;
    os << " Number of elements ........: " << m.m_elements << std::endl;
    os << " Number of allocated cells .: " << m.m_alloc << std::endl;

    // Dump elements and compression schemes
    m.dump_elements(os);
    m.dump_row_comp(os);
    m.dump_col_comp(os);

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Return matrix with absolute values of all elements
 *
 * @param[in] m Matrix for which absolute values are to be returned.
 ***************************************************************************/
GMatrix abs(const GMatrix& m)
{
    // Define result matrix
    GMatrix result = m;

    // Convert all elements to absolute values  
    for (int i = 0; i < result.m_elements; ++i)
        result.m_data[i] = std::abs(result.m_data[i]);

    // Return result
    return result;
}
