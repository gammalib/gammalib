/***************************************************************************
 *                   GMatrixSparse.cpp - Sparse matrix class               *
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
 * ----------------------------------------------------------------------- *
 * This class implements the compressed sparse column format. The          *
 * following arrays are allocated:                                         *
 *                                                                         *
 * m_rows     (n)  Number of rows                                          *
 * m_cols     (m)  Number of columns                                       *
 * m_data     (Ax) Holds all 'elements' non-zero values, in column order   *
 * m_colstart (Ap) Holds the index of the first element of each column     *
 * m_rowinx   (Ai) Holds the row indices for all elements                  *
 *                                                                         *
 * Column 'col' covers therefore [m_colstart[j], ..., m_colstart[j+1]-1]   *
 * Note that 'm_colstart' has m_cols+1 elements.                           *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GMatrixSparse.cpp
 * @brief Sparse matrix class implementation
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
#include "GSparseSymbolic.hpp"
#include "GSparseNumeric.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS1                       "GMatrixSparse::operator(int&, int&)"
#define G_ACCESS2                       "GMatrixSparse::operator(int&, int&)"
#define G_OP_MUL_VEC                     "GMatrixSparse::operator*(GVector&)"
#define G_OP_ADD                  "GMatrixSparse::operator+=(GMatrixSparse&)"
#define G_OP_SUB                  "GMatrixSparse::operator-=(GMatrixSparse&)"
#define G_OP_MAT_MUL              "GMatrixSparse::operator*=(GMatrixSparse&)"
#define G_INVERT                                "GMatrixSparse::invert(void)"
#define G_ADD_COL                    "GMatrixSparse::add_col(GVector&, int&)"
#define G_ADD_COL2         "GMatrixSparse::add_col(double*, int*, int, int&)"
#define G_CHOL_DECOMP               "GMatrixSparse::cholesky_decompose(bool)"
#define G_CHOL_SOLVE         "GMatrixSparse::cholesky_solver(GVector&, bool)"
#define G_EXTRACT_ROW                      "GMatrixSparse::extract_row(int&)"
#define G_EXTRACT_COL                      "GMatrixSparse::extract_col(int&)"
#define G_INSERT_COL              "GMatrixSparse::insert_col(GVector&, int&)"
#define G_INSERT_COL2   "GMatrixSparse::insert_col(double*, int*, int, int&)"
#define G_STACK_INIT                  "GMatrixSparse::stack_init(int&, int&)"
#define G_STACK_PUSH  "GMatrixSparse::stack_push_column(double*, int*, int&,"\
                                                                     " int&)"
#define G_STACK_FLUSH                      "GMatrixSparse::stack_flush(void)"
#define G_COPY_MEMBERS          "GMatrixSparse::copy_members(GMatrixSparse&)"
#define G_ALLOC_MEMBERS      "GMatrixSparse::alloc_members(int&, int&, int&)"
#define G_ALLOC                   "GMatrixSparse::alloc_elements(int&, int&)"
#define G_FREE                     "GMatrixSparse::free_elements(int&, int&)"
#define G_REMOVE_ZERO              "GMatrixSparse::remove_zero_row_col(void)"
#define G_SYMPERM                    "cs_symperm(GMatrixSparse*, int*, int&)"
#define G_TRANSPOSE                       "cs_transpose(GMatrixSparse*, int)"

/* __ Macros _____________________________________________________________ */
#define G_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define G_MAX(a,b) (((a) > (b)) ? (a) : (b))

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_SPARSE_PENDING                    // Analyse pending values
//#define G_DEBUG_SPARSE_INSERTION                 // Analyse value insertion
//#define G_DEBUG_SPARSE_ADDITION                   // Analyse value addition
//#define G_DEBUG_SPARSE_COMPRESSION   // Analyse zero row/column compression
//#define G_DEBUG_SPARSE_MALLOC                  // Analyse memory management
//#define G_DEBUG_SPARSE_STACK_PUSH           // Analyse stack column pushing
//#define G_DEBUG_SPARSE_STACK_FLUSH         // Analyse stack column flushing


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void matrix constructor
 *
 * This method will allocate a sparse matrix object without any elements.
 * The number of rows and columns of the matrix will be zero.
 *
 * @todo Verify that the class is save against empty matrix objects.
 ***************************************************************************/
GMatrixSparse::GMatrixSparse(void) : GMatrixBase()
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
 * @param[in] elements Number of allocated elements (default: 0).
 *
 * This method allocates a sparse matrix object with the specified number
 * of rows and columns. The elements parameter allows to specify how much
 * physical memory should be allocated initially. By default, no memory
 * will be allocated.
 ***************************************************************************/
GMatrixSparse::GMatrixSparse(const int& rows, const int& cols,
                             const int& elements) : GMatrixBase()
{
    // Initialise private members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(rows, cols, elements);

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrix to GMatrixSparse storage class convertor
 *
 * @param[in] matrix Generic matrix (GMatrix).
 *
 * This constructor converts a generic matrix (of type GMatrix) into a
 * sparse matrix. 
 ***************************************************************************/
GMatrixSparse::GMatrixSparse(const GMatrix& matrix) : GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Construct matrix
    alloc_members(matrix.rows(), matrix.cols());

    // Fill matrix column by column
    for (int col = 0; col < matrix.cols(); ++col) {
        GVector vector = matrix.extract_col(col);
        this->insert_col(vector, col);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrixSymmetric to GMatrixSparse storage class convertor
 *
 * @param[in] matrix Symmetric matrix (GMatrixSymmetric).
 *
 * This constructor converts a symmetric matrix (of type GSymMatrix) into a
 * sparse matrix. 
 ***************************************************************************/
GMatrixSparse::GMatrixSparse(const GMatrixSymmetric& matrix) : GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(matrix.rows(), matrix.cols());

    // Fill matrix column by column
    for (int col = 0; col < matrix.cols(); ++col) {
        GVector vector = matrix.extract_col(col);
        this->insert_col(vector, col);
    }

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
GMatrixSparse::GMatrixSparse(const GMatrixSparse& matrix) : GMatrixBase(matrix)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members. Note that the copy operation does not fill the
    // pending element.
    copy_members(matrix);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GMatrixSparse::~GMatrixSparse(void)
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
GMatrixSparse& GMatrixSparse::operator=(const GMatrixSparse& matrix)
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

    // Return
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
double& GMatrixSparse::operator()(const int& row, const int& col)
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_ACCESS1, row, col, m_rows, m_cols);
    }
    #endif

    // Get element
    fill_pending();
    int inx = get_index(row,col);
    double* value;
    if (inx < 0) {
        value      = &m_fill_val;
        m_fill_row = row;
        m_fill_col = col;
    }
    else {
        value = &(m_data[inx]);
    }

    // Return element
    return *value;
}


/***********************************************************************//**
 * @brief Access operator (const version)
 *
 * @param[in] row Matrix row.
 * @param[in] col Matrix column.
 *
 * @exception GException::out_of_range
 *            Row or column index out of range.
 *
 * We need here the zero element to return also a pointer for 0.0 entry that
 * is not stored. Since we have the const version we don't have to care about
 * modification of this zero value.
 ***************************************************************************/
const double& GMatrixSparse::operator()(const int& row, const int& col) const
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols) {
      throw GException::out_of_range(G_ACCESS2, row, col, m_rows, m_cols);
    }
    #endif

    // Get element
    int inx = get_index(row,col);
    double* value;
    if (inx < 0) {
        value = (double*)&m_zero;
    }
    else if (inx == m_elements) {
        value = (double*)&m_fill_val;
    }
    else {
        value = (double*)&(m_data[inx]);
    }

    // Return element
    return *value;
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
 *
 * The method includes any pending fill.
 ***************************************************************************/
GVector GMatrixSparse::operator*(const GVector& vector) const
{
    // Raise an exception if the matrix and vector dimensions are not compatible
    if (m_cols != vector.size()) {
        throw GException::matrix_vector_mismatch(G_OP_MUL_VEC, vector.size(),
                                                 m_rows, m_cols);
    }

    // Initialise result vector
    GVector result(m_rows);

    // Multiply only if there are elements in matrix
    if (m_elements > 0) {

        // Perform vector multiplication
        for (int col = 0; col < m_cols; ++col) {
            int     i_start    = m_colstart[col];
            int     i_stop     = m_colstart[col+1];
            double* ptr_data   = m_data   + i_start;
            int*    ptr_rowinx = m_rowinx + i_start;
            for (int i = i_start; i < i_stop; ++i) {
                result[*ptr_rowinx++] += *ptr_data++ * vector[col];
            }
        }

        // If fill is pending then add-in also this element into the product
        if (m_fill_val != 0.0) {
            result[m_fill_row] += m_fill_val * vector[m_fill_col];
            #if defined(G_DEBUG_SPARSE_PENDING)
            std::cout << G_OP_MUL_VEC << ": pending value " << m_fill_val
                      << " for location (" <<  m_fill_row << "," << m_fill_col
                      << ") has been used." << std::endl;
            #endif
        }

    } // endif: there were elements in matrix

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Equalty operator
 *
 * @param[in] matrix Matrix.
 *
 * This operator checks if two matrices are identical. Two matrices are
 * considered identical if they have the same dimensions and identicial
 * elements.
 *
 * @todo Implement native sparse code
 ***************************************************************************/
bool GMatrixSparse::operator==(const GMatrixSparse &matrix) const
{
    // Initalise the result to 'equal matrices'
    bool result = true;

    // Perform comparison (only if matrix dimensions are identical)
    if (m_rows == matrix.m_rows && m_cols == matrix.m_cols) {

        // Loop over all matrix elements
        for (int row = 0; row < m_rows; ++row) {
            for (int col = 0; col < m_cols; ++col) {
                if ((*this)(row,col) != matrix(row,col)) {
                    result = false;
                    break;
                }
            }
            if (!result) {
                break;
            }
        }
    }
    else {
        result = false;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] matrix Matrix.
 *
 * This operator checks if two matrices are not identical. Two matrices are
 * considered not identical if they differ in their dimensions or if at
 * least one element differs.
 ***************************************************************************/
bool GMatrixSparse::operator!=(const GMatrixSparse &matrix) const
{
    // Get negated result of equality operation
    bool result = !(this->operator==(matrix));

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
 *
 * @todo Implement native sparse code
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator+=(const GMatrixSparse& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows || m_cols != matrix.m_cols) {
        throw GException::matrix_mismatch(G_OP_ADD,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Perform inplace matrix addition using vectors
    for (int col = 0; col < m_cols; ++col) {
        GVector v_result  = extract_col(col);
        GVector v_operand = matrix.extract_col(col);
        v_result += v_operand;
        insert_col(v_result, col);
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
 *
 * @todo Implement native sparse code
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator-=(const GMatrixSparse& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows || m_cols != matrix.m_cols) {
        throw GException::matrix_mismatch(G_OP_SUB,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Perform inplace matrix subtraction
    for (int col = 0; col < m_cols; ++col) {
        GVector v_result  = extract_col(col);
        GVector v_operand = matrix.extract_col(col);
        v_result -= v_operand;
        insert_col(v_result, col);
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
 * @todo Implement native sparse code
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator*=(const GMatrixSparse& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_cols != matrix.m_rows) {
        throw GException::matrix_mismatch(G_OP_MAT_MUL,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Allocate result matrix
    GMatrixSparse result(m_rows, matrix.m_cols);

    // Multiply only if there are elements in both matrices
    if (m_elements > 0 && matrix.m_elements > 0) {

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

    }

    // Assign result
    *this = result;

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
void GMatrixSparse::clear(void)
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
GMatrixSparse* GMatrixSparse::clone(void) const
{
    // Clone this image
    return new GMatrixSparse(*this);
}


/***********************************************************************//**
 * @brief Transpose matrix
 *
 * The transpose operation exchanges the number of rows against the number
 * of columns.
 ***************************************************************************/
void GMatrixSparse::transpose(void)
{
    // Fill pending element
    fill_pending();
    
    // Compute the transpose
    *this = cs_transpose(*this, 1);
    
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
void GMatrixSparse::invert(void)
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
 * Adds the contect of a vector to a matrix column.
 *
 * This is the main driver routine to add data to a matrix. It handles both
 * normal and stack-based filled. Note that there is another instance of this
 * method that takes a compressed array.
 ***************************************************************************/
void GMatrixSparse::add_col(const GVector& vector, const int& col)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_ADDITION)
    std::cout << "GMatrixSparse::add_col([" << v << "], " << col << "):" << std::endl;
    std::cout << " In Data : ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_data[i] << " ";
    }
    std::cout << std::endl << " In Row .: ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_rowinx[i] << " ";
    }
    std::cout << std::endl << " In Col .: ";
    for (int i = 0; i < m_cols+1; ++i) {
        std::cout << m_colstart[i] << " ";
    }
    std::cout << std::endl;
    #endif

    // Initialise number of non-zero elements to 0
    int non_zero = 0;

    // If we have a stack then try to push vector on stack first. Note that
    // stack_push_column does its own argument verifications, so to avoid
    // double checking we don't do anything before this call ...
    if (m_stack_data != NULL) {
        non_zero = stack_push_column(vector, col);
        if (non_zero == 0) {
            return;
        }
    }

    // ... otherwise check first the arguments and determine the number of
    // non-zero elements in vector
    else {

        // Raise an exception if the column index is invalid
        #if defined(G_RANGE_CHECK)
        if (col < 0 || col >= m_cols) {
            throw GException::out_of_range(G_ADD_COL, col, 0, m_cols-1);
        }
        #endif

        // Raise an exception if the matrix and vector dimensions are incompatible
        if (m_rows != vector.size()) {
            throw GException::matrix_vector_mismatch(G_ADD_COL, vector.size(),
                                                     m_rows, m_cols);
        }

        // Determine number of non-zero elements in vector
        non_zero = vector.non_zeros();
    }

    // Extract vector for column, add elements, and re-insert vector (only if
    // vector to insert has non-zeros)
    if (non_zero > 0) {

        // Copy input vector
        GVector column = vector;

        // Add elements to vector
        for (int i = m_colstart[col]; i < m_colstart[col+1]; ++i) {
            column[m_rowinx[i]] += m_data[i];
        }

        // If there is a pending element then put it in the vector
        if (m_fill_val != 0.0 && m_fill_col == col) {
            column[m_fill_row] += m_fill_val;
        }

        // Insert vector into matrix
        insert_col(column, col);

    }

    // Debugging: show sparse matrix after addition
    #if defined(G_DEBUG_SPARSE_ADDITION)
    std::cout << " Out Data: ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_data[i] << " ";
    }
    std::cout << std::endl << " Out Row : ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_rowinx[i] << " ";
    }
    std::cout << std::endl << " Out Col : ";
    for (int i = 0; i < m_cols+1; ++i) {
        std::cout << m_colstart[i] << " ";
    }
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add compressed array into matrix column
 *
 * @param[in] values Compressed array.
 * @param[in] rows Row indices of array.
 * @param[in] number Number of elements in array.
 * @param[in] col Column index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * Adds the content of a compressed array into a matrix column.
 *
 * This is the main driver routine to add data to a matrix. It handles both
 * normal and stack-based filled. Note that there is another instance of this
 * method that takes a vector.
 ***************************************************************************/
void GMatrixSparse::add_col(const double* values, const int* rows, 
                            int number, const int& col)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_ADDITION)
    std::cout << "GMatrixSparse::add_col(v, i, n, " << col << "):" << std::endl;
    std::cout << " Matrix Data : ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_data[i] << " ";
    }
    std::cout << std::endl << " Matrix Row .: ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_rowinx[i] << " ";
    }
    std::cout << std::endl << " Matrix Col .: ";
    for (int i = 0; i < m_cols+1; ++i) {
        std::cout << m_colstart[i] << " ";
    }
    std::cout << std::endl;
    std::cout << " Array Data .: ";
    for (int i = 0; i < number; ++i) {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl << " Array Row ..: ";
    for (int i = 0; i < number; ++i) {
        std::cout << rows[i] << " ";
    }
    std::cout << std::endl;
    #endif

    // If we have a stack then try to push elements on stack first. Note that
    // stack_push_column does its own argument verifications, so to avoid
    // double checking we don't do anything before this call ...
    if (m_stack_data != NULL) {
        number = stack_push_column(values, rows, number, col);
        if (number == 0) {
            return;
        }
    }

    // ... otherwise check the arguments
    else {
        // If the array is empty there is nothing to do
        if (!values || !rows || (number < 1)) {
            return;
        }

        // Raise an exception if the column index is invalid
        #if defined(G_RANGE_CHECK)
        if (col < 0 || col >= m_cols) {
            throw GException::out_of_range(G_ADD_COL2, col, 0, m_cols-1);
        }
        #endif

        // Raise an exception if the matrix and vector dimensions are incompatible
        if (rows[number-1] >= m_rows) {
            throw GException::matrix_vector_mismatch(G_ADD_COL2, rows[number-1],
                                                     m_rows, m_cols);
        }

    } // endelse: there was no stack

    // Get indices of column in matrix
    int i_start = m_colstart[col];
    int i_stop  = m_colstart[col+1];

    // Case A: the column exists in the matrix, so mix new elements with existing
    // data
    if (i_start < i_stop) {

        // Fill pending element before the merge (it will be lost otherwise)
        fill_pending();

        // Allocate workspace to hold combined column
        int     wrk_size   = number + i_stop - i_start;
        double* wrk_double = new double[wrk_size];
        int*    wrk_int    = new int[wrk_size];

        // Mix matrix column with specified data
        int num_mix;
        mix_column(&(m_data[i_start]), &(m_rowinx[i_start]), i_stop-i_start,
                   values, rows, number,
                   wrk_double, wrk_int, &num_mix);

        // Insert mixed column
        insert_col(wrk_double, wrk_int, num_mix, col);

        // Free workspace
        delete [] wrk_int;
        delete [] wrk_double;

    } // endif: Case A

    // Case B: the column does not yet exist in the matrix, so just insert it
    else {
        insert_col(values, rows, number, col);
    }

    // Debugging: show sparse matrix after insertion
    #if defined(G_DEBUG_SPARSE_ADDITION)
    std::cout << " Out Data: ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_data[i] << " ";
    std::cout << std::endl << " Out Row : ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_rowinx[i] << " ";
    std::cout << std::endl << " Out Col : ";
    for (int i = 0; i < m_cols+1; ++i)
        std::cout << m_colstart[i] << " ";
    std::cout << std::endl;
    #endif

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
 *
 * This is the main driver routine to insert data into a matrix. Note that
 * there is another instance of this function that takes a compressed array.
 ***************************************************************************/
void GMatrixSparse::insert_col(const GVector& vector, const int& col)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_INSERTION)
    std::cout << "GMatrixSparse::insert_col([" << v << "], " << col << "):" << std::endl;
    std::cout << " In Data : ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_data[i] << " ";
    std::cout << std::endl << " In Row .: ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_rowinx[i] << " ";
    std::cout << std::endl << " In Col .: ";
    for (int i = 0; i < m_cols+1; ++i)
        std::cout << m_colstart[i] << " ";
    std::cout << std::endl;
    #endif

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

    // If there is a pending element for this column then delete it since
    // the vector overwrites this element
    if (m_fill_val != 0.0 && m_fill_col == col) {
        #if defined(G_DEBUG_SPARSE_PENDING)
        std::cout << G_INSERT_COL << ": pending value " << m_fill_val << 
                     " for location (" << m_fill_row << "," << m_fill_col << 
                     ") became obsolete" << std::endl;
        #endif
        m_fill_val = 0.0;
        m_fill_row = 0;
        m_fill_col = 0;
    }

    // Determine the number of non-zero elements in the vector
    int n_vector = 0;
    for (int row = 0; row < m_rows; ++row) {
        if (vector[row] != 0.0) {
            n_vector++;
        }
    }

    // Get the start and stop indices of the actual column and compute
    // the number of exisiting elements in the column
    int i_start = m_colstart[col];
    int i_stop  = m_colstart[col+1];
    int n_exist = i_stop - i_start;

    // Compute the size difference for the new matrix. It is positive if
    // the number of non-zero entries in the vector is larger than the
    // number of non-zero entries in the matrix (in this case we have to
    // increase the matrix size).
    int n_diff = n_vector - n_exist;

    // If we need space then allocate it, if we have to much space then free it
    if (n_diff > 0) {
        alloc_elements(i_start, n_diff);
        #if defined(G_DEBUG_SPARSE_INSERTION)
        std::cout << " Insert .: " << n_diff << " elements at index " << i_start << std::endl;
        #endif
    }
    else if (n_diff < 0) {
        free_elements(i_start, -n_diff);
        #if defined(G_DEBUG_SPARSE_INSERTION)
        std::cout << " Remove .: " << -n_diff << " elements at index " << i_start << std::endl;
        #endif
    }

    // Insert the vector elements in the matrix
    if (n_vector > 0) {
        for (int row = 0, i = i_start; row < m_rows; ++row) {
            if (vector[row] != 0.0) {
                m_data[i]   = vector[row];
                m_rowinx[i] = row;
                i++;
            }
        }
    }

    // Update column start indices
    for (int i = col+1; i <= m_cols; ++i) {
        m_colstart[i] += n_diff;
    }

    // Debugging: show sparse matrix after insertion
    #if defined(G_DEBUG_SPARSE_INSERTION)
    std::cout << " Out Data: ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_data[i] << " ";
    std::cout << std::endl << " Out Row : ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_rowinx[i] << " ";
    std::cout << endl << " Out Col : ";
    for (int i = 0; i < m_cols+1; ++i)
        std::cout << m_colstart[i] << " ";
    std::cout << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert compressed array into matrix column
 *
 * @param[in] values Compressed array.
 * @param[in] rows Row indices of array.
 * @param[in] number Number of elements in array.
 * @param[in] col Column index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * Inserts the content of a copressed array into a matrix column. Any 
 * previous content in the matrix column will be overwritted.
 *
 * This is the main driver routine to insert data into a matrix. Note that
 * there is another instance of this function that takes a vector.
 ***************************************************************************/
void GMatrixSparse::insert_col(const double* values, const int* rows, 
                               int number, const int& col)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_INSERTION)
    std::cout << "GMatrixSparse::insert_col(v, i, n, " << col << "):" << std::endl;
    std::cout << " Matrix Data : ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_data[i] << " ";
    }
    std::cout << std::endl << " Matrix Row .: ";
    for (int i = 0; i < m_elements; ++i) {
        std::cout << m_rowinx[i] << " ";
    }
    std::cout << std::endl << " Matrix Col .: ";
    for (int i = 0; i < m_cols+1; ++i) {
        std::cout << m_colstart[i] << " ";
    }
    std::cout << std::endl;
    std::cout << " Array Data .: ";
    for (int i = 0; i < number; ++i) {
        std::cout << values[i] << " ";
    }
    std::cout << std::endl << " Array Row ..: ";
    for (int i = 0; i < number; ++i) {
        std::cout << rows[i] << " ";
    }
    std::cout << std::endl;
    #endif

    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_INSERT_COL2, col, 0, m_cols-1);
    }
    #endif

    // Raise an exception if the index array seems incompatible with matrix 
    // dimensions
    if (rows[number-1] >= m_rows) {
        throw GException::matrix_vector_mismatch(G_INSERT_COL2, rows[number-1],
                                                 m_rows, m_cols);
    }

    // If there is a pending element for this column then delete it since
    // the vector overwrites this element
    if (m_fill_val != 0.0 && m_fill_col == col) {
        #if defined(G_DEBUG_SPARSE_PENDING)
        std::cout << G_INSERT_COL2 << ": pending value " << m_fill_val << 
                " for location (" << m_fill_row << "," << m_fill_col << 
                ") became obsolete" << std::endl;
        #endif
        m_fill_val = 0.0;
        m_fill_row = 0;
        m_fill_col = 0;
    }

    // Get the start and stop indices of the actual column and compute
    // the number of exisiting elements in the column
    int i_start = m_colstart[col];
    int i_stop  = m_colstart[col+1];
    int n_exist = i_stop - i_start;

    // If the array is empty then make sure that the number of elements is 0 
    // (we then just delete the existing column)
    if (!values || !rows) {
        number = 0;
    }

    // Compute the size difference for the new matrix. It is positive if
    // the number of non-zero entries in the array is larger than the
    // number of non-zero entries in the matrix (in this case we have to
    // increase the matrix size).
    int n_diff = number - n_exist;

    // If we need space then allocate it, if we have to much space then free it
    if (n_diff > 0) {
        alloc_elements(i_start, n_diff);
        #if defined(G_DEBUG_SPARSE_INSERTION)
        std::cout << " Insert .: " << n_diff << " elements at index " << i_start << std::endl;
        #endif
    }
    else if (n_diff < 0) {
        free_elements(i_start, -n_diff);
        #if defined(G_DEBUG_SPARSE_INSERTION)
        std::cout << " Remove .: " << -n_diff << " elements at index " << i_start << std::endl;
        #endif
    }

    // Insert the array elements into the matrix
    if (number > 0) {
        for (int row = 0, i = i_start; row < number; ++row, ++i) {
            m_data[i]   = values[row];
            m_rowinx[i] = rows[row];
        }
    }

    // Update column start indices
    for (int i = col+1; i <= m_cols; ++i) {
        m_colstart[i] += n_diff;
    }

    // Debugging: show sparse matrix after insertion
    #if defined(G_DEBUG_SPARSE_INSERTION)
    std::cout << " Out Data: ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_data[i] << " ";
    std::cout << std::endl << " Out Row : ";
    for (int i = 0; i < m_elements; ++i)
        std::cout << m_rowinx[i] << " ";
    std::cout << std::endl << " Out Col : ";
    for (int i = 0; i < m_cols+1; ++i)
        std::cout << m_colstart[i] << " ";
    std::cout << std::endl;
    #endif

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
GVector GMatrixSparse::extract_row(const int& row) const
{
    // Raise an exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_EXTRACT_ROW, row, 0, m_rows-1);
    }
    #endif

    // Create result vector
    GVector result(m_cols);

    // Loop over all columns to extract data
    for (int col = 0; col < m_cols; ++col) {

        // Get the start and stop of the elements
        int i_start = m_colstart[col];
        int i_stop  = m_colstart[col+1];

        // Search requested row in elements
        int i;
        for (i = i_start; i < i_stop; ++i) {
            if (m_rowinx[i] == row) {
                break;
            }
        }

        // Copy element if we found one
        if (i < i_stop) {
            result[col] = m_data[i];
        }

    } // endfor: looped over all columns

    // If there is a pending element then put it in the vector
    if (m_fill_val != 0.0 && m_fill_row == row) {
        result[m_fill_col] = m_fill_val;
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
GVector GMatrixSparse::extract_col(const int& col) const
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_EXTRACT_COL, col, 0, m_cols-1);
    }
    #endif

    // Create result vector
    GVector result(m_rows);

    // Get the start and stop of the elements
    int i_start = m_colstart[col];
    int i_stop  = m_colstart[col+1];

    // Extract elements into vector
    for (int i = i_start; i < i_stop; ++i) {
        result[m_rowinx[i]] = m_data[i];
    }

    // If there is a pending element then put it in the vector
    if (m_fill_val != 0.0 && m_fill_col == col) {
        result[m_fill_row] = m_fill_val;
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Perform a Cholesky decomposition
 *
 * @param[in] compress Use zero-row/column compression (default: true).
 *
 * Cholesky decomposition of a sparse matrix. The decomposition is stored
 * within the GMatrixSparse without destroying the original matrix.
 ***************************************************************************/
void GMatrixSparse::cholesky_decompose(bool compress)
{
    // Save original matrix size
    int matrix_rows = m_rows;
    int matrix_cols = m_cols;

    // Delete any existing symbolic and numeric analysis object and reset
    // pointers
    if (m_symbolic != NULL) delete m_symbolic;
    if (m_numeric  != NULL) delete m_numeric;
    m_symbolic = NULL;
    m_numeric  = NULL;

    // Allocate symbolic analysis object
    GSparseSymbolic* symbolic = new GSparseSymbolic();

    // Declare numeric analysis object. We don't allocate one since we'll
    // throw it away at the end of the function (the L matrix will be copied
    // in this object)
    GSparseNumeric numeric;

    // Fill pending element into matrix
    fill_pending();

    // Remove rows and columns containing only zeros if matrix compression
    // has been selected
    if (compress) {
        remove_zero_row_col();
    }

    // Ordering an symbolic analysis of matrix. This sets up an array 'pinv'
    // which contains the fill-in reducing permutations
    symbolic->cholesky_symbolic_analysis(1, *this);

    // Store symbolic pointer in sparse matrix object
    m_symbolic = symbolic;

    // Perform numeric Cholesky decomposition
    numeric.cholesky_numeric_analysis(*this, *symbolic);

    // Copy L matrix into this object
    free_elements(0, m_elements);
    alloc_elements(0, numeric.m_L->m_elements);
    for (int i = 0; i < m_elements; ++i) {
        m_data[i]   = numeric.m_L->m_data[i];
        m_rowinx[i] = numeric.m_L->m_rowinx[i];
    }
    for (int col = 0; col <= m_cols; ++col) {
        m_colstart[col] = numeric.m_L->m_colstart[col];
    }

    // Insert zero rows and columns if they have been removed previously.
    if (compress) {
        insert_zero_row_col(matrix_rows, matrix_cols);
    }

  // Return
  return;
}


/***********************************************************************//**
 * @brief Cholesky solver
 *
 * @param[in] vector Solution vector.
 * @param[in] compress Request matrix compression.
 *
 * @exception GException::matrix_vector_mismatch
 *            Matrix and vector do not match.
 * @exception GException::matrix_not_factorised
 *            Matrix has not been factorised.
 * 
 * Solves the linear equation A*x=b using a Cholesky decomposition of A.
 * This function is to be applied on a GMatrixSparse matrix for which a
 * Choleksy factorization has been produced using 'cholesky_decompose'.
 ***************************************************************************/
GVector GMatrixSparse::cholesky_solver(const GVector& vector, bool compress)
{
    // Dump header
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << "GMatrixSparse::cholesky_solver" << std::endl;
    std::cout << " Input vector .....: " << vector << std::endl;
    #endif

    // Raise an exception if the matrix and vector dimensions are incompatible
    if (m_rows != vector.size()) {
        throw GException::matrix_vector_mismatch(G_CHOL_SOLVE, vector.size(),
                                                 m_rows, m_cols);
    }
    
    // Raise an exception if there is no symbolic pointer
    if (!m_symbolic) {
        throw GException::matrix_not_factorised(G_CHOL_SOLVE, 
                                                "Cholesky decomposition");
    }

    // Raise an exception if there is no permutation
    if (!m_symbolic->m_pinv) {
        throw GException::matrix_not_factorised(G_CHOL_SOLVE, 
                                                "Cholesky decomposition");
    }

    // Flag row and column compression
    int row_compressed = (m_rowsel != NULL && m_num_rowsel < m_rows);
    int col_compressed = (m_colsel != NULL && m_num_colsel < m_cols);

    // Decide if we need a compression algorithm or not
    int no_zero = !(compress && (row_compressed || col_compressed));

    // Allocate vector for permutation and result vector
    GVector result(m_cols);

    // Setup pointers to L matrix and x vector
    int*    Lp = m_colstart;
    int*    Li = m_rowinx; 
    double* Lx = m_data;

    // Case A: no zero-row/col compression needed
    if (no_zero) {

        // Declare working vector
        GVector x(m_rows);

        // Perform inverse vector permutation
        for (int i = 0; i < vector.size(); ++i) {
            x[m_symbolic->m_pinv[i]] = vector[i];
        }

        // Inplace solve L\x=x
        for (int col = 0; col < m_cols; ++col) {            // loop over columns
            x[col] /= Lx[Lp[col]];                          // divide by diag.
            for (int p = Lp[col]+1; p < Lp[col+1]; p++)     // loop over elements
                x[Li[p]] -= Lx[p] * x[col];
        }

        // Inplace solve L'\x=x
        for (int col = m_cols-1; col >= 0; --col) {         // loop over columns
            for (int p = Lp[col]+1; p < Lp[col+1]; p++)     // loop over elements
                x[col] -= Lx[p] * x[Li[p]];
                x[col] /= Lx[Lp[col]];
        }

        // Perform vector permutation
        for (int i = 0; i < m_cols; ++i) {
            result[i] = x[m_symbolic->m_pinv[i]];
        }

    } // endif: Case A

    // Case B: zero-row/column compression requested
    else {

        // Allocate row and column mapping arrays
        int* row_map = new int[m_rows];
        int* col_map = new int[m_cols];

        // Setup row mapping array that maps original matrix rows into compressed
        // matrix rows. An entry of -1 indicates that the row should be dropped.
        // If no selection exists then setup an identity map.
        if (row_compressed) {
            for (int row = 0; row < m_rows; ++row) {
                row_map[row] = -1;
            }
            for (int c_row = 0; c_row < m_num_rowsel; ++c_row) {
                row_map[m_rowsel[c_row]] = c_row;
            }
        }
        else {
            for (int row = 0; row < m_rows; ++row) {
                row_map[row] = row;
            }
        }
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Row mapping ......:";
        for (int row = 0; row < m_rows; ++row) {
            std::cout << " " << row_map[row];
        }
        std::cout << std::endl;
        #endif

        // Setup column mapping array that maps original matrix column into compressed
        // matrix columns. An entry of -1 indicates that the column should be dropped
        // If no selection exists then setup an identity map.
        if (col_compressed) {
            for (int col = 0; col < m_cols; ++col) {
                col_map[col] = -1;
            }
            for (int c_col = 0; c_col < m_num_colsel; ++c_col) {
                col_map[m_colsel[c_col]] = c_col;
            }
        }
        else {
            for (int col = 0; col < m_cols; ++col) {
                col_map[col] = col;
            }
        }
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Column mapping ...:";
        for (int col = 0; col < m_cols; ++col) {
            std::cout << " " << col_map[col];
        }
        std::cout << std::endl;
        #endif

        // Declare working vector
        GVector x(row_compressed ? m_num_rowsel : m_rows);

        // Compress input vector v -> c_v if required
        if (m_rowsel != NULL && m_num_rowsel < m_rows) {
            for (int c_row = 0; c_row < m_num_rowsel; ++c_row) {
                x[c_row] = vector[m_rowsel[c_row]];
            }
        }
        else {
            x = vector;
        }
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Compressed vector : " << x << std::endl;
        #endif

        // Perform inverse permutation
        x = iperm(x, m_symbolic->m_pinv);
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Permutated vector : " << x << std::endl;
        #endif

        // Inplace solve L\x=x. The column and row maps are just use to see which
        // columns or rows should be skipped in the calculations.
        for (int col = 0; col < m_cols; ++col) {              // loop over columns
            int c_col = col_map[col];
            if (c_col >= 0) {                                 // use only non-zero cols
                x[c_col] /= Lx[Lp[col]];                      // divide by diag.
                for (int p = Lp[col]+1; p < Lp[col+1]; p++) { // loop over elements
                    int c_row = row_map[Li[p]];
                    if (c_row >= 0) {                         // use only non-zero rows
                        x[c_row] -= Lx[p] * x[c_col];
                    }
                }
            }
        }
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Solve Lx=x .......: " << x << std::endl;
        #endif

        // Inplace solve L'\x=x. The column and row maps are just use to see which
        // columns or rows should be skipped in the calculations.
        for (int col = m_cols-1; col >= 0; --col) {           // loop over columns
            int c_col = col_map[col];
            if (c_col >= 0) {                                 // use only non-zero cols
                for (int p = Lp[col]+1; p < Lp[col+1]; p++) { // loop over elements
                    int c_row = row_map[Li[p]];
                    if (c_row >= 0) {                         // use only non-zero rows
                        x[c_col] -= Lx[p] * x[c_row];
                    }
                }
                x[c_col] /= Lx[Lp[col]];
            }
        }
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Solve L'x=x ......: " << x << std::endl;
        #endif

        // Perform vector permutation
        x = perm(x, m_symbolic->m_pinv);
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Permutated vector : " << x << std::endl;
        #endif

        // If column compression has been performed the expand the result vector
        // accordingly
        if (m_colsel != NULL && m_num_colsel < m_cols) {
            for (int c_col = 0; c_col < m_num_colsel; ++c_col) {
                result[m_colsel[c_col]] = x[c_col];
            }
        }
        else {
            result = x;
        }
        #if defined(G_DEBUG_SPARSE_COMPRESSION)
        std::cout << " Restored vector ..: " << result << std::endl;
        #endif

        // Free mapping arrays
        delete [] row_map;
        delete [] col_map;

    } // endelse: Case B

    // Return result vector
    return result;
}


/***********************************************************************//**
 * @brief Cholesky invert
 *
 * @param[in] compress Use zero-row/column compression (default: true).
 *
 * @exception GException::matrix_zero
 *            All matrix elements are zero.
 *
 * Inplace matrix inversion using the Cholesky decomposition.
 ***************************************************************************/
void GMatrixSparse::cholesky_invert(bool compress)
{
    // Generate Cholesky decomposition of matrix
    cholesky_decompose(compress);

    // Allocate result matrix and unit vector
    GMatrixSparse result(m_rows, m_cols);
    GVector       unit(m_rows);

    // Column-wise solving of the problem
    for (int col = 0; col < m_cols; ++col) {

        // Set unit vector
        unit[col] = 1.0;

        // Solve for column
        GVector x = cholesky_solver(unit, compress);

        // Insert column in matrix
        result.insert_col(x, col);

        // Clear unit vector for next round
        unit[col] = 0.0;

    }

    // Assign result
    *this = result;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns fill of matrix
 *
 * The fill of a matrix is defined as the number non-zero elements devided
 * by the number of total elements. By definiton, the fill is comprised
 * in the interval [0,..,1]. The fill of an undefined matrix is defined to
 * be 0.
 ***************************************************************************/
double GMatrixSparse::fill(void) const
{
    // Initialise result
    double result = 0.0;

    // Compute matrix size
    int size = m_rows*m_cols;

    // Continue only if matrix has elements
    if (size > 0) {

        // Determine number of elements in matrix
        int num = (m_fill_val == 0.0) ? m_elements : m_elements + 1;

        // Compute fill
        result = double(num) / double(size);

    } // endif: there were elements in matrix

    // Return fill
    return result;
}


/***********************************************************************//**
 * @brief Return minimum matrix element
 ***************************************************************************/
double GMatrixSparse::min(void) const
{
    // Initialise minimum with fill value
    double result = m_fill_val;

    // Search all elements for the smallest one
    for (int i = 0; i < m_elements; ++i) {
        if (m_data[i] < result) {
            result = m_data[i];
        }
    }

    // If minimum is > 0.0 and there are zero elements then set the minimum
    // to 0.0
    if ((result > 0.0) && (m_elements < (m_rows*m_cols))) {
        result = 0.0;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return maximum matrix element
 ***************************************************************************/
double GMatrixSparse::max(void) const
{
    // Initialise maximum with fill value
    double result = m_fill_val;

    // Search all elements for the largest one
    for (int i = 0; i < m_elements; ++i) {
        if (m_data[i] > result) {
            result = m_data[i];
        }
    }

    // If maximum is < 0.0 and there are zero elements then set the maximum
    // to 0.0
    if ((result < 0.0) && (m_elements < (m_rows*m_cols))) {
        result = 0.0;
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Sum matrix elements
 ***************************************************************************/
double GMatrixSparse::sum(void) const
{
    // Initialise matrix sum with fill value
    double result = m_fill_val;

    // Add all elements
    for (int i = 0; i < m_elements; ++i) {
        result += m_data[i];
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print matrix
 ***************************************************************************/
std::string GMatrixSparse::print(void) const
{
    // Initialise result string
    std::string result;

    // Determine number of elements
    int nonzero;
    int elements = m_elements;
    if (m_fill_val == 0.0) {
        nonzero  = (m_colstart != NULL) ? m_colstart[m_cols] : 0;
    }
    else {
        nonzero  = (m_colstart != NULL) ? m_colstart[m_cols]+1 : 0;
        elements++;
    }

    // Append header
    result.append("=== GMatrixSparse ===");
    result.append("\n"+parformat("Number of rows")+str(m_rows));
    if (m_rowsel != NULL) {
        result.append(" (compressed "+str(m_num_rowsel)+")");
    }
    result.append("\n"+parformat("Number of columns")+str(m_cols));
    if (m_colsel != NULL) {
        result.append(" (compressed "+str(m_num_colsel)+")");
    }
    result.append("\n"+parformat("Number of nonzero elements")+str(nonzero));
    if (m_fill_val != 0.0) {
        result.append("\n"+parformat("Pending element"));
        result.append("("+str(m_fill_row)+","+str(m_fill_col)+")=");
        result.append(str(m_fill_val));
    }
    result.append("\n"+parformat("Number of allocated cells")+str(m_alloc));
    result.append("\n"+parformat("Memory block size")+str(m_mem_block));
    result.append("\n"+parformat("Sparse matrix fill")+str(fill()));
    result.append("\n"+parformat("Pending element")+str(m_fill_val));
    result.append("\n"+parformat("Fill stack size")+str(m_stack_size));
    if (m_stack_data == NULL) {
        result.append(" (none)");
    }
    
    // Append elements and compression schemes
    result.append(print_elements());
    result.append(print_row_compression());
    result.append(print_col_compression());

    // Append symbolic decomposition if available
    //if (m_symbolic != NULL) {
    //    result.append(m_symbolic->print());
    //}

    // Append numeric decomposition if available
    //if (m_numeric != NULL) {
    //    result.append(m_numeric->print());
    //}

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Initialises matrix filling stack
 *
 * @param[in] size Stack size.
 * @param[in] entries Maximum number of entries.
 *
 * The matrix filling stack is used to allow for a fast column-wise filling
 * of a sparse matrix. Columns are successively appended to a stack which is
 * regularily flushed when it is full. This reduces memory copies and
 * movements and increases filling speed.
 ***************************************************************************/
void GMatrixSparse::stack_init(const int& size, const int& entries)
{
    // Free exisiting stack
    free_stack_members();

    // Initialise stack members
    init_stack_members();
    m_stack_max_entries = (entries > 0) ? entries : m_cols;
    m_stack_size        = (size    > 0) ? size    : G_SPARSE_MATRIX_DEFAULT_STACK_SIZE;

    // Allocate stack memory. Raise an exception if allocation fails
    m_stack_colinx = new int[m_stack_max_entries];
    m_stack_start  = new int[m_stack_max_entries+1];
    m_stack_data   = new double[m_stack_size];
    m_stack_rowinx = new int[m_stack_size];
    m_stack_work   = new int[m_cols];
    m_stack_rows   = new int[m_cols];
    m_stack_values = new double[m_cols];

    // Initialise next free stack location to the first stack element
    m_stack_start[0] = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Push a vector column on the matrix stack
 *
 * @param[in] vector Vector.
 * @param[in] col Column index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * Adds the contect of a vector to the matrix stack. This method is
 * identical to the GMatrixSparse::stack_push_column method that uses a
 * compressed array, yet it takes a full vector. Internal working arrays
 * are used to convert the full column vector in a compressed array and to
 * hand it over to the compressed array version.
 ***************************************************************************/
int GMatrixSparse::stack_push_column(const GVector& vector, const int& col)
{
    // Initialise number of non-zero elements
    int non_zero = 0;

    // Compress vector in the buffer and set-up row index array
    for (int i = 0; i < vector.size(); ++i) {
        if (vector[i] != 0.0) {
            m_stack_values[non_zero] = vector[i];
            m_stack_rows[non_zero]   = i;
            non_zero++;
        }
    }

    // Hand over to compressed array method
    int remaining = stack_push_column(m_stack_values, m_stack_rows, non_zero, col);

    // Return remaining number of elements
    return remaining;
}


/***********************************************************************//**
 * @brief Push a compressed array on the matrix stack
 *
 * @param[in] values Compressed array.
 * @param[in] rows Row indices of array.
 * @param[in] number Number of elements in array.
 * @param[in] col Column index (starting from 0).
 *
 * @exception GException::out_of_range
 *            Invalid column index specified.
 * @exception GException::matrix_vector_mismatch
 *            Matrix dimension mismatches the vector size.
 *
 * The method puts new data on the stack while assuring that each column
 * is present only once. If an already existing column is encountered
 * the data from the existing and new column are mixed and put into a new
 * entry one the stack; the old entry is signalled as obsolete.
 *
 * On return the method indicates the number of elements in the input
 * array that have not been put onto the stack (due to memory limits).
 * This value can be either '0' (i.e. all elements have been put on the
 * stack) or 'number' (i.e. none of the elements have been put on the
 * stack). Columns are not partially put onto the stack.
 ***************************************************************************/
int GMatrixSparse::stack_push_column(const double* values, const int* rows,
                                     const int& number, const int& col)
{
    // Initialise return value
    int remaining = number; 

    // Single loop for common exit point
    do {

        // If the array is empty there is nothing to do
        if (!values || !rows || (number < 1)) {
            continue;
        }

        // If there is no stack or the stack can not hold the requested 
        // number of elements then report number of array elements to the
        // caller
        if (m_stack_data == NULL || number > m_stack_size) {
            continue;
        }

        // Raise an exception if the column index is invalid
        #if defined(G_RANGE_CHECK)
        if (col < 0 || col >= m_cols) {
            throw GException::out_of_range(G_STACK_PUSH, col, 0, m_cols-1);
        }
        #endif

        // Raise an exception if the matrix and vector dimensions are
        // incompatible
        if (rows[number-1] >= m_rows) {
            throw GException::matrix_vector_mismatch(G_STACK_PUSH, rows[number-1],
                                                     m_rows, m_cols);
        }

        // Debug header
        #if defined(G_DEBUG_SPARSE_STACK_PUSH)
        std::cout << "GMatrixSparse::stack_push_column(v, i, n, col=" << col << ")";
        std::cout << std::endl;
        std::cout << " Data to push on stack ...:";
        for (int i = 0; i < number; ++i) {
            std::cout << " " << values[i];
        }
        std::cout << std::endl;
        std::cout << " Row indices of data .....:";
        for (int i = 0; i < number; ++i) {
            std::cout << " " << rows[i];
        }
        std::cout << std::endl;
        #endif

        // If the stack is full then flush it before pushing the array onto it.
        // There are 2 conditions that may fill the stack: all entries are
        // occupied or there is not enough space to hold more elements
        if ((m_stack_entries >= m_stack_max_entries) ||
            (number          >= (m_stack_size - m_stack_start[m_stack_entries]))) {
            stack_flush();
        }

        // If the specified column is already on the stack then mix new column
        // with old one and invalidate old column (by setting its column index
        // to -1).
        // Since we loop here over all stack entries we don't want to have too
        // many entries in the stack. So don't set the maximum number of stack
        // entries too large ...
        for (int entry = 0; entry < m_stack_entries; ++entry) {

            // Is the specified column already on the stack?
            if (col == m_stack_colinx[entry]) {

                // Get start and stop indices of existing column in stack
                int i_start = m_stack_start[entry];
                int i_stop  = m_stack_start[entry+1];

                // Allocate variables to hold mixing results
                int num_1;    // Number of elements that are only in stack column
                int num_2;    // Number of elements that are only in new column
                int num_mix;  // Number of elements that are in both columns

                // Determine how many elements are requested to hold the combined
                // column
                mix_column_prepare(&(m_stack_rowinx[i_start]), i_stop-i_start,
                                   rows, number, &num_1, &num_2, &num_mix);
                int num_request = num_1 + num_2 + num_mix;

                // If there is not enough space in the stack to hold the combined
                // column we flush the stack and exit the loop now. In this case
                // the new column is put as the first column in the empty (flushed
                // stack)
                if (num_request >= (m_stack_size - m_stack_start[m_stack_entries])) {
                    stack_flush();
                    break;
                }

                // There is enough space, so combine both columns into a fresh one
                int inx = m_stack_start[m_stack_entries];
                int num_new;
                mix_column(&(m_stack_data[i_start]), &(m_stack_rowinx[i_start]),
                           i_stop-i_start,
                           values, rows, number,
                           &(m_stack_data[inx]), &(m_stack_rowinx[inx]), &num_new);

                // Store entry information and initialise start of next entry
                m_stack_colinx[m_stack_entries] = col;           // Store column index for entry
                m_stack_entries++;                               // Increase the # of entries
                m_stack_start[m_stack_entries] = inx + num_new;  // Set start pointer for next entry

                // Invalidate old column
                m_stack_colinx[entry] = -1;

                // Fall through to end
                remaining = 0;
                break;

            } // endif: we found an existing column
        } // endfor: looped over all columns

        // If column has already been inserted then fall through ...
        if (remaining == 0) {
            continue;
        }

        // Otherwise, push the array on the stack
        int inx = m_stack_start[m_stack_entries];
        for (int i = 0; i < number; ++i) {
            m_stack_data[inx]   = values[i];
            m_stack_rowinx[inx] = rows[i];
            inx++;
        }

        // Store entry information and initialise start of next entry
        m_stack_colinx[m_stack_entries] = col;  // Store column index for entry
        m_stack_entries++;                      // Increase the # of entries
        m_stack_start[m_stack_entries] = inx;   // Set start pointer for next entry

        // Signal success
        remaining = 0;

    } while (0); // End of main do-loop

    // Debug: show stack information  
    #if defined(G_DEBUG_SPARSE_STACK_PUSH)
    std::cout << " Number of stack entries .: " << m_stack_entries << std::endl;
    std::cout << " Stack entry columns .....:";
    for (int i = 0; i < m_stack_entries; ++i) {
        std::cout << " " << m_stack_colinx[i];
    }
    std::cout << std::endl;
    std::cout << " Stack entry starts ......:";
    for (int i = 0; i < m_stack_entries; ++i) {
        std::cout << " " << m_stack_start[i];
    }
    std::cout << " (next at " << m_stack_start[m_stack_entries] << ")" << std::endl;
    std::cout << " Stack data ..............:";
    for (int i = 0; i < m_stack_start[m_stack_entries]; ++i) {
        std::cout << " " << m_stack_data[i];
    }
    std::cout << std::endl;
    std::cout << " Stack rows ..............:";
    for (int i = 0; i < m_stack_start[m_stack_entries]; ++i) {
        std::cout << " " << m_stack_rowinx[i];
    }
    std::cout << std::endl;
    #endif

    // Return remaining number of elements
    return remaining;
}


/***********************************************************************//**
 * @brief Flush matrix stack
 *
 * Adds the stack to the actual matrix. First, the total number of matrix
 * elements is determined. Then new memory is allocated to hold all
 * elements. Finally, all elements are filled into a working array that is
 * then compressed into the matrix.
 *
 * The method uses the internal working array m_stack_work.
 *
 * NOTE: The present algorithm assumes that each column occurs only once
 * in the stack!
 ***************************************************************************/
void GMatrixSparse::stack_flush(void)
{
    // Do nothing if there is no stack
    if (m_stack_data == NULL) {
        return;
    }

    // Fill pending value
    fill_pending();

    // Debug header
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    std::cout << "GMatrixSparse::stack_flush" << std::endl;
    std::cout << " Number of stack entries .: " << m_stack_entries << std::endl;
    std::cout << " Number of stack elements : " << m_stack_start[m_stack_entries] << std::endl;
    std::cout << " Number of matrix elements: " << m_elements << std::endl;
    #endif

    // Use stack working array to flag all columns that exist already in the
    // matrix by 1 and all non-existing columns by 0
    for (int col = 0; col < m_cols; ++col) {
        if (m_colstart[col] < m_colstart[col+1]) {
            m_stack_work[col] = 1;
        }
        else {
            m_stack_work[col] = 0;
        }
    }

    // Initialise some element counters
    int new_elements = 0;
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    int num_matrix   = 0;    // Number of elements only in matrix
    int num_stack    = 0;    // Number of elements only on stack
    int num_both     = 0;    // Number of elements in matrix and on stack
    #endif

    // Loop over all stack entries and gather the number of elements that are
    // new with respect to the initial matrix. For each column set a flag with
    // the following meanings:
    //  -(entry+2): Column exists and is also found in stack=entry
    //  +(entry+2): Column exists in stack=entry only
    //           1: Column exists in matrix only
    //           0: Column does neither exist in matrix or stack  
    for (int entry = 0; entry < m_stack_entries; ++entry) {

        // Get column for actual entry
        int col = m_stack_colinx[entry];

        // Consider only valid entries
        if (col >= 0) {

            // If column exists already in matrix then determine total number
            // of additional elements
            if (m_stack_work[col] == 1) {

                // Flag that this column is a mixed column
                m_stack_work[col] = -(entry+2);

                // Setup index boundaries
                int i_start = m_colstart[col];
                int i_stop  = m_colstart[col+1];
                int k_start = m_stack_start[entry];
                int k_stop  = m_stack_start[entry+1];

                // Allocate output counters
                int num_1;
                int num_2;
                int num_mix;

                // Analyse column mixing
                mix_column_prepare(&(m_rowinx[i_start]), i_stop-i_start,
                                   &(m_stack_rowinx[k_start]), k_stop-k_start,
                                   &num_1, &num_2, &num_mix);

                // Update element counters
                new_elements += num_2;
                #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
                num_matrix   += num_1;
                num_stack    += num_2;
                num_both     += num_mix;
                #endif

            } // endif: column existed in the matrix

            // If column did not exists in the matrix then consider all
            // elements as new
            else {
                m_stack_work[col]  = (entry+2);
                new_elements      += (m_stack_start[entry+1] - m_stack_start[entry]);
            }
        } // endif: entry was valid
    } // endfor: looped over all entries
    int elements = m_elements + new_elements;

    // Dump number of elements in new matrix 
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    std::cout << " New elements ............: " << new_elements << std::endl;
    #endif

    // Allocate memory for new matrix (always keep some elbow room)
    m_alloc = elements + m_mem_block;
    double* new_data   = new double[m_alloc];
    int*    new_rowinx = new int[m_alloc];

    // Fill new matrix. For this purpose we loop over all matrix columns
    // and perform the operation that was identified in the previous scan
    int index = 0;
    for (int col = 0; col < m_cols; ++col) {

        // If column does not exist then skip
        if (m_stack_work[col] == 0) {
            m_colstart[col] = index;
            continue;
        }

        // If column exists in the matrix but not in the stack we copy the
        // existing column into the new matrix
        if (m_stack_work[col] == 1) {
            for (int i = m_colstart[col]; i < m_colstart[col+1]; ++i) {
                new_data[index]   = m_data[i];
                new_rowinx[index] = m_rowinx[i];
                index++;
                #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
                num_matrix++;
                #endif
            }
        }

        // If column is new (i.e. it did not exist in the matrix before)
        // we copy the stack into the new matrix
        else if (m_stack_work[col] > 1) {
            int entry = m_stack_work[col] - 2;
            for (int i = m_stack_start[entry]; i < m_stack_start[entry+1]; ++i) {
                new_data[index]   = m_stack_data[i];
                new_rowinx[index] = m_stack_rowinx[i];
                index++;
                #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
                num_stack++;
                #endif
            }
        }

        // If column exists in the matrix and in the stack then we have to mix both
        else {

            // Get stack entry for column mix
            int entry = -(m_stack_work[col] + 2);

            // Setup index boundaries
            int i_start = m_colstart[col];
            int i_stop  = m_colstart[col+1];
            int k_start = m_stack_start[entry];
            int k_stop  = m_stack_start[entry+1];

            // Allocate output counter
            int num;

            // Perform column mixing
            mix_column(&(m_data[i_start]), &(m_rowinx[i_start]), i_stop-i_start,
                 &(m_stack_data[k_start]), &(m_stack_rowinx[k_start]), k_stop-k_start,
                 &(new_data[index]), &(new_rowinx[index]), &num);

            // Increment element index
            index += num;

        } // endelse: column mixing required

        // Store actual index in column start array
        m_colstart[col] = index;

    } // endfor: looped over all columns

    // Dump number of elements in new matrix after addition
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    std::cout << " Added elements ..........: " << index << " (should be " << elements << ")" << std::endl;
    std::cout << " - Matrix only ...........: " << num_matrix << std::endl;
    std::cout << " - Stack only ............: " << num_stack << std::endl;
    std::cout << " - Matrix & Stack ........: " << num_both << std::endl;
    #endif

    // Correct columns start array
    for (int col = m_cols; col > 0; --col) {
        m_colstart[col] = m_colstart[col-1];
    }
    m_colstart[0] = 0;

    // Delete old matrix memory
    if (m_data   != NULL) delete [] m_data;
    if (m_rowinx != NULL) delete [] m_rowinx;

    // Update pointers to new memory and update element counter
    m_data     = new_data;
    m_rowinx   = new_rowinx;
    m_elements = elements;

    // Stack is empty now, so reset stack counters
    m_stack_entries  = 0;
    m_stack_start[0] = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destroy matrix stack
 *
 * Flush and destroy matrix stack
 ***************************************************************************/
void GMatrixSparse::stack_destroy(void)
{
    // Flush stack first
    stack_flush();

    // Free stack members
    free_stack_members();

    // Initialise stack (no entries)
    init_stack_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class mambers
 ***************************************************************************/
void GMatrixSparse::init_members(void)
{
    // Initialise sparse matrix members
    m_rowinx    = NULL;
    m_mem_block = G_SPARSE_MATRIX_DEFAULT_MEM_BLOCK;
    m_zero      = 0.0;
    m_fill_val  = 0.0;
    m_fill_row  = 0;
    m_fill_col  = 0;
    m_symbolic  = NULL;
    m_numeric   = NULL;

    // Initialise stack members
    init_stack_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] matrix Matrix to be copied.
 ***************************************************************************/
void GMatrixSparse::copy_members(const GMatrixSparse& matrix)
{
    // Copy GMatrixSparse members
    m_mem_block = matrix.m_mem_block;
    m_zero      = matrix.m_zero;
    m_fill_val  = matrix.m_fill_val;
    m_fill_row  = matrix.m_fill_row;
    m_fill_col  = matrix.m_fill_col;

    // Copy row indices if they exist
    if (matrix.m_rowinx != NULL && matrix.m_alloc > 0) {
        m_rowinx = new int[matrix.m_alloc];
        for (int i = 0; i < matrix.m_elements; ++i) {
            m_rowinx[i] = matrix.m_rowinx[i];
        }
    }

    // Clone symbolic decomposition if it exists
    if (matrix.m_symbolic != NULL) {
        m_symbolic  = new GSparseSymbolic();
        *m_symbolic = *matrix.m_symbolic;
    }

    // Copy numeric decomposition if it exists
    if (matrix.m_numeric != NULL) {
        m_numeric  = new GSparseNumeric();
        *m_numeric = *matrix.m_numeric;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMatrixSparse::free_members(void)
{
    // De-allocate only if memory has indeed been allocated
    if (m_numeric  != NULL) delete m_numeric;
    if (m_symbolic != NULL) delete m_symbolic;
    if (m_rowinx   != NULL) delete [] m_rowinx;

    // Properly mark members as free
    m_rowinx   = NULL;
    m_symbolic = NULL;
    m_numeric  = NULL;

    // Free stack members
    free_stack_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate matrix
 *
 * @param[in] rows Number of rows (>1).
 * @param[in] cols Number of columns (>1).
 * @param[in] elements Number of matrix elements to be physically allocated.
 *
 * @exception GException::empty
 *            Attempt to allocate zero size matrix.
 *
 * This is the main constructor code that allocates and initialises memory
 * for matrix elements.
 ***************************************************************************/
void GMatrixSparse::alloc_members(const int& rows, const int& cols,
                                  const int& elements)
{
    // Throw exception if requested matrix size is zero
    if (rows == 0 || cols == 0) {
        throw GException::empty(G_ALLOC_MEMBERS);
    }

    // Allocate column start array. This is the only array that we can
    // allocate at this time. The other arrays can only be allocated during
    // filling of the matrix
    m_colstart = new int[cols+1];

    // Store (logical) matrix size
    m_rows = rows;
    m_cols = cols;

    // Initialise column start indices to 0
    for (int col = 0; col <= m_cols; ++col) {
        m_colstart[col] = 0;
    }

    // Optionally allocate memory for matrix elements
    if (elements > 0) {
        alloc_elements(0, elements);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise fill stack
 *
 * The fill stack is used to fill the sparse matrix without any prior know-
 * ledge about the number of location of the non-zero matrix elements.
 * Values are first filled into the stack and flushed into the sparse
 * matrix once the stack is full.
 ***************************************************************************/
void GMatrixSparse::init_stack_members(void)
{
    // Initialise stack members
    m_stack_max_entries = 0;
    m_stack_size        = 0;
    m_stack_entries     = 0;
    m_stack_colinx      = NULL;
    m_stack_start       = NULL;
    m_stack_data        = NULL;
    m_stack_rowinx      = NULL;
    m_stack_work        = NULL;
    m_stack_rows        = NULL;
    m_stack_values      = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete fill-stack class members
 *
 * Deletes all memory that has been allocated using the init_stack_members()
 * method.
 ***************************************************************************/
void GMatrixSparse::free_stack_members(void)
{
    // Free stack members
    if (m_stack_colinx != NULL) delete [] m_stack_colinx;
    if (m_stack_start  != NULL) delete [] m_stack_start;
    if (m_stack_data   != NULL) delete [] m_stack_data;
    if (m_stack_rowinx != NULL) delete [] m_stack_rowinx;
    if (m_stack_work   != NULL) delete [] m_stack_work;
    if (m_stack_rows   != NULL) delete [] m_stack_rows;
    if (m_stack_values != NULL) delete [] m_stack_values;

    // Properly mark members as free
    m_stack_colinx = NULL;
    m_stack_start  = NULL;
    m_stack_data   = NULL;
    m_stack_rowinx = NULL;
    m_stack_work   = NULL;
    m_stack_rows   = NULL;
    m_stack_values = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Determines element index for (row,col)
 *
 * @param[in] row Row index.
 * @param[in] col Column index.
 *
 * Returns the index in the compressed array for (row,col). The following
 * special results exist:
 * -1: Requested index does not exist in the matrix.
 * m_elements: Requested index is the pending element.
 ***************************************************************************/
int GMatrixSparse::get_index(const int& row, const int& col) const
{
    // Initialise element to 'not found'
    int index = -1;

    // If we have a pending element then check if this element is requested
    if (m_fill_val != 0.0) {
        if (row == m_fill_row && col == m_fill_col) {
            return m_elements;
        }
    }

    // If requested element is not the pending element then check if it exists
    // in the matrix. Only if it is found its index is returned. Otherwise
    // the default index is -1, signalling that the element is absent.
    if (m_elements > 0) {
        int* ptr_colstart = m_colstart + col;
        int  i_start      = *ptr_colstart++;
        int  i_stop       = *ptr_colstart;
        int* ptr_rowinx   = m_rowinx + i_start;
        for (int i = i_start; i < i_stop; ++i) {
            int row_test = *ptr_rowinx++;
            if (row_test == row) {
                index = i;
                break;
            }
        }
    }

    // Return index
    return index;
}


/***********************************************************************//**
 * @brief Fills pending matrix element
 *
 * If 'm_fill_val' is non-zero a pending matrix element exists that should
 * be filled into (row,col)=(m_fill_row,m_fill_col). This routine performs
 * the filling of the matrix with this element and resets 'm_fill_val' to
 * zero. This routine allows for element-by-element filling of a sparse
 * matrix. This is, of course, very time consuming and should in general
 * be avoided. However, it allows to design a sparse matrix class that
 * hides the matrix sparsity completely to the user.
 ***************************************************************************/
void GMatrixSparse::fill_pending(void)
{
    // If we have a pending element then fill it into the matrix
    if (m_fill_val != 0.0) {

        // Debugging
        #if defined(G_DEBUG_SPARSE_PENDING) || defined(G_DEBUG_SPARSE_INSERTION)
        std::cout << "GMatrixSparse::fill_pending(): pending value " <<
                     m_fill_val << " will be filled in location (" <<
                     m_fill_row << "," << m_fill_col <<  ")" << std::endl;
        #endif

        // If there are so far no elements in the matrix then append element ...
        int inx_insert;
        if (m_elements == 0) {
            inx_insert = 0;
        }

        // ... otherwise search for index to insert
        else {
            int* ptr_colstart = m_colstart + m_fill_col;
            int  i_start      = *ptr_colstart++;
            int  i_stop       = *ptr_colstart;
            int* ptr_rowinx   = m_rowinx + i_start;
            for (inx_insert = i_start; inx_insert < i_stop; ++inx_insert) {
                int row_test = *ptr_rowinx++;
                if (row_test > m_fill_row) {
                    break;
                }
            }
        }

        // Allocate memory for new element
        alloc_elements(inx_insert, 1);

        // Insert element
        m_data[inx_insert]   = m_fill_val;
        m_rowinx[inx_insert] = m_fill_row;

        // Update column start indices
        for (int col = m_fill_col+1; col <= m_cols; ++col) {
            m_colstart[col] += 1;
        }

        // Reset fill value
        m_fill_val = 0.0;
        m_fill_row = 0;
        m_fill_col = 0;

        // Debugging: show sparse matrix after filling
        #if defined(G_DEBUG_SPARSE_PENDING) || defined(G_DEBUG_SPARSE_INSERTION)
        std::cout << " Data: ";
        for (int i = 0; i < m_elements; ++i)
            std::cout << m_data[i] << " ";
        std::cout << std::endl << " Row.: ";
        for (int i = 0; i < m_elements; ++i)
            std::cout << m_rowinx[i] << " ";
        std::cout << std::endl << " Col.: ";
        for (int i = 0; i < m_cols+1; ++i)
            std::cout << m_colstart[i] << " ";
        std::cout << std::endl;
        #endif

    } // endif: a pending matrix element was found

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate memory for new matrix elements
 *
 * @param[in] start Index of first allocated element.
 * @param[in] num Number of elements to be allocated.
 *
 * Inserts a memory allocation for 'num' elements at the index 'start'.
 * The new elements are filled with 0.0 and the corresponding row indices
 * are set to 0.
 *
 * NOTE: This method does not take care of updating 'm_colstart'. This has
 * to be done by the client.
 ***************************************************************************/
void GMatrixSparse::alloc_elements(int start, const int& num)
{
    // Dump header
    #if defined(G_DEBUG_SPARSE_MALLOC)
    std::cout << "GMatrixSparse::alloc_elements(start=" << start << ", num=" << 
                 num << ")" << std::endl;
    std::cout << " Before allocation : " << m_elements << " " << m_alloc << std::endl;
    #endif

    // Continue only if we need memory
    if (num > 0) {

        // If start is after the end then append memory
        if (start > m_elements) {
            start = m_elements;
        }

        // Determine the requested new logical size of the matrix
        int new_size = m_elements + num;

        // Case A: the requested memory is already available, so just move the
        // data to make space for new elements and initialise the new cells
        if (new_size <= m_alloc) {

            // Move up all elements after index to insert
            for (int i = m_elements - 1; i >= start; --i) {
                m_data[i+num]   = m_data[i];
                m_rowinx[i+num] = m_rowinx[i];
            }

            // Clear new elements (zero content for row 0)
            for (int i = start; i < start+num; ++i) {
                m_data[i]   = 0.0;
                m_rowinx[i] = 0;
            }

            // Update element counter
            m_elements += num;

        } // endif: Case A: memory already existed

        // Case B: more memory is needed, so allocate it, copy the existing
        // content, and initialise new cells
        else {

            // Make sure that enough memory is allocated
            int new_propose = m_alloc + m_mem_block;
            m_alloc = (new_size > new_propose) ? new_size : new_propose;

            // Allocate memory for new elements
            double* new_data   = new double[m_alloc];
            int*    new_rowinx = new int[m_alloc];

            // Copy all elements before index to insert
            for (int i = 0; i < start; ++i) {
                new_data[i]   = m_data[i];
                new_rowinx[i] = m_rowinx[i];
            }

            // Clear new elements (zero content for row 0)
            for (int i = start; i < start+num; ++i) {
                new_data[i]   = 0.0;
                new_rowinx[i] = 0;
            }

            // Copy all elements after index to insert
            for (int i = start; i < m_elements; ++i) {
                new_data[i+num]   = m_data[i];
                new_rowinx[i+num] = m_rowinx[i];
            }

            // Delete old memory
            if (m_data   != NULL) delete [] m_data;
            if (m_rowinx != NULL) delete [] m_rowinx;

            // Update pointers to new memory and update element counter
            m_data      = new_data;
            m_rowinx    = new_rowinx;
            m_elements += num;

        } // endelse: Case B: more memory was needed

    } // endif: needed new memory

    // Dump new memory size
    #if defined(G_DEBUG_SPARSE_MALLOC)
    std::cout << " After allocation .: " << m_elements << " " << m_alloc << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Free memory for obsolete matrix elements
 *
 * @param[in] start Index of first freed element.
 * @param[in] num Number of elements to be freed.
 *
 * Free memory for 'num' elements starting from index 'start'.
 *
 * NOTE: This method does not take care of updating 'm_colstart'. This has
 * to be done by the client.
 ***************************************************************************/
void GMatrixSparse::free_elements(const int& start, const int& num)
{
    // Dump header
    #if defined(G_DEBUG_SPARSE_MALLOC)
    std::cout << "GMatrixSparse::free_elements(start=" << start << ", num=" << 
                 num << ")" << std::endl;
    #endif

    // Continue only if we need to free memory and if start is within the
    // range
    if (num > 0 && start < m_elements) {

        // Determine the requested new logical size of the matrix
        int new_size = m_elements - num;

        // If there are no elements then simply delete all matrix elements ...
        if (new_size < 1) {
            if (m_data   != NULL) delete [] m_data;
            if (m_rowinx != NULL) delete [] m_rowinx;
            m_data     = NULL;
            m_rowinx   = NULL;
            m_elements = 0;
            m_alloc    = 0;
        }

        // ... otherwise shrink the array
        else {

            // Case A: If at least one entire memory block has been liberated then 
            // physically shrink the matrix array
            if (m_alloc - new_size > m_mem_block) {

                // Shrink, but leave a memory block for possible future filling
                m_alloc = new_size + m_mem_block;

                // Allocate new memory
                double* new_data   = new double[m_alloc];
                int*    new_rowinx = new int[m_alloc];

                // Copy all elements before the starting index
                for (int i = 0; i < start; ++i) {
                    new_data[i]   = m_data[i];
                    new_rowinx[i] = m_rowinx[i];
                }

                // Copy all elements after the starting index
                for (int i = start; i < new_size; ++i) {
                    new_data[i]   = m_data[i+num];
                    new_rowinx[i] = m_rowinx[i+num];
                }

                // Delete old memory
                if (m_data   != NULL) delete [] m_data;
                if (m_rowinx != NULL) delete [] m_rowinx;

                // Update pointers to new memory and update element counter
                m_data     = new_data;
                m_rowinx   = new_rowinx;
                m_elements = new_size;

            } // endif: Case A: memory shrinkage performed

            // Case B: we keep the memory and just move the elements
            else {

                // Move all elements after the starting index
                for (int i = start; i < new_size; ++i) {
                    m_data[i]   = m_data[i+num];
                    m_rowinx[i] = m_rowinx[i+num];
                }

                // Update element counter
                m_elements = new_size;

            } // endelse: Case B

        } // endif: array shrinkage needed
    } // endif: needed new memory

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove rows and columns with zeros
 *
 * @exception GException::matrix_zero
 *            All matrix elements are zero.
 *
 * Remove all rows and columns that contain only zeros from matrix. This
 * function is needed for compressed matrix factorisation. The resulting
 * matrix has reduced size (number of rows and columns) and dimensions
 * (number of elements). Note that the physical memory is not reduced by
 * the method.
 ***************************************************************************/
void GMatrixSparse::remove_zero_row_col(void)
{
    // Dump header
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << "GMatrixSparse::remove_zero_row_col" << std::endl;
    #endif

    // Fill pending value in matrix
    fill_pending();

    // Select non-zero rows and columns of matrix
    select_non_zero();

    // Stop if there is no compression
    if (m_rows == m_num_rowsel && m_cols == m_num_colsel) {
        return;
    }

    // Raise exception if all matrix elements are zero
    if (m_num_rowsel < 1 || m_num_colsel < 1) {
        throw GException::matrix_zero(G_REMOVE_ZERO);
    }

    // Allocate row mapping array
    int* row_map = new int[m_rows];

    // Setup row mapping array that maps original matrix rows into compressed
    // matrix rows. An entry of -1 indicates that the row should be dropped
    for (int row = 0; row < m_rows; ++row) {
        row_map[row] = -1;
    }
    for (int c_row = 0; c_row < m_num_rowsel; ++c_row) {
        row_map[m_rowsel[c_row]] = c_row;
    }

    // Initialise pointers to compressed array
    double* d_data   = m_data;
    int*    d_rowinx = m_rowinx;

    // Initialise column start of first column to zero
    m_colstart[0] = 0;

    // Dump column start array
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << " before compression (" << m_colstart[m_cols] << "):";
    for (int col = 0; col <= m_cols; ++col)
        std::cout << " " << m_colstart[col];
    std::cout << std::endl;
    #endif

    // Loop over all columns of compressed matrix
    for (int c_col = 0; c_col < m_num_colsel; ++c_col) {

        // Get index range for elements in original matrix
        int i_start = m_colstart[m_colsel[c_col]];
        int i_stop  = m_colstart[m_colsel[c_col]+1];

        // Initialise number of element counter for the compressed column
        int num = 0;

        // Loop over all elements in original column
        int c_row;
        for (int i = i_start; i < i_stop; ++i) {
            if ((c_row = row_map[m_rowinx[i]]) >= 0) {
                *d_data++   = m_data[i];
                *d_rowinx++ = c_row;
                num++;
            }
        }

        // Update column stop
        m_colstart[c_col+1] = m_colstart[c_col] + num;

    } // endfor: looped over all columns of compressed matrix

    // Dump column start array
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << " after compression (" << m_colstart[m_num_colsel] << ") :";
    for (int c_col = 0; c_col <= m_num_colsel; ++c_col)
        std::cout << " " << m_colstart[c_col];
    std::cout << std::endl;
    #endif

    // Free row mapping array
    delete [] row_map;

    // Update matrix size
    m_rows = m_num_rowsel;
    m_cols = m_num_colsel;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert zero rows and columns
 *
 * @param[in] rows Number of rows.
 * @param[in] cols Number of columns.
 *
 * Insert zero rows and columns into matrix. Since for a sparse matrix
 * this does not require any allocation of additional memory, the data are
 * not moved by this function, but the pointers are re-arranged.
 ***************************************************************************/
void GMatrixSparse::insert_zero_row_col(const int& rows, const int& cols)
{
    // Dump header
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << "GMatrixSparse::insert_zero_row_col(" << rows << "," << cols << 
                 ")" << std::endl;
    #endif

    // Fill pending value in matrix
    fill_pending();

    // Stop if there is no insertion
    if (m_rows == rows && m_cols == cols) {
        return;
    }

    // If row selection exists then restore row indices
    if (m_rowsel != NULL) {
        for (int i = 0; i < m_elements; ++i) {
            m_rowinx[i] = m_rowsel[m_rowinx[i]];
        }
    }

    // Dump column start array
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << " before restoration (" << m_colstart[m_num_colsel] << "):";
    for (int c_col = 0; c_col <= m_num_colsel; ++c_col)
        std::cout << " " << m_colstart[c_col];
    std::cout << std::endl;
    #endif

    // If column selection exists then restore column counters
    if (m_colsel != NULL) {

        // Start insertion from the last original column
        int col_stop = cols - 1;

        // Loop over all compressed columns
        for (int c_col = m_num_colsel-1; c_col > 0; --c_col) {
            int col_start = m_colsel[c_col-1] + 1;
            int col_value = m_colstart[c_col];
            for (int col = col_start; col <= col_stop; ++col) {
                m_colstart[col] = col_value;
            }
            col_stop = col_start - 1;
        }

        // Set first columns if they are not yet set
        for (int col = 0; col <= col_stop; ++col) {
            m_colstart[col] = 0;
        }

        // Restore the number of elements
        m_colstart[cols] = m_elements;
    }

    // Dump column start array
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << " after restoration (" << m_colstart[cols] << ") :";
    for (int col = 0; col <= cols; ++col) {
        std::cout << " " << m_colstart[col];
    }
    std::cout << std::endl;
    #endif

    // Update matrix size
    m_rows = rows;
    m_cols = cols;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Prepare mix of sparse columns
 *
 * @param[in] src1_row Row index array [0...src1_num-1] of first column.
 * @param[in] src1_num Number of elements in first column.
 * @param[in] src2_row Row index array [0...src2_num-1] of second column.
 * @param[in] src2_num Number of elements in second column.
 * @param[out] num_1 Number of elements only found in column 1.
 * @param[out] num_2 Number of elements only found in column 2.
 * @param[out] num_mix Number of elements found in both columns.
 *
 * This method prepares the mix of two sparse matrix columns into a single
 * column. 
 ***************************************************************************/
void GMatrixSparse::mix_column_prepare(const int* src1_row, int src1_num,
                                       const int* src2_row, int src2_num,
                                       int* num_1, int* num_2, int* num_mix)
{
    // Initialise element counters
    *num_1   = 0;
    *num_2   = 0;
    *num_mix = 0;

    // Initialise indices and row indices of both columns
    int inx_1 = 0;                    // Column 1 element index
    int inx_2 = 0;                    // Column 2 element index
    int row_1 = src1_row[inx_1];      // Column 1 first row index
    int row_2 = src2_row[inx_2];      // Column 2 first row index

    // Mix elements of both columns while both contain still elements
    while (inx_1 < src1_num && inx_2 < src2_num) {

        // Case A: the element exist in both columns
        if (row_1 == row_2) {
            row_1 = src1_row[++inx_1];
            row_2 = src2_row[++inx_2];
            (*num_mix)++;
        }

        // Case B: the element exists only in first column
        else if (row_1 < row_2) {
            row_1 = src1_row[++inx_1];
            (*num_1)++;
        }

        // Case C: the element exists only in second column
        else {
            row_2 = src2_row[++inx_2];
            (*num_2)++;
        }

    } // endwhile: mixing

    // At this point either the first or the second column expired of elements
    // In the case that there are still elements remaining in the first column we
    // count them now ...
    if (inx_1 < src1_num) {
        *num_1 += (src1_num - inx_1);
    }

    // ... or in the case that there are still elements remaining in the second
    // column we count them now
    if (inx_2 < src2_num) {
        *num_2 += (src2_num - inx_2);
    }

    // We're done
    return;
}		


/***********************************************************************//**
 * @brief Mix of sparse columns
 *
 * @param[in] src1_data Data array [0...src1_num-1] of first column.
 * @param[in] src1_row Row index array [0...src1_num-1] of first column.
 * @param[in] src1_num Number of elements in first column.
 * @param[in] src2_data Data array [0...src2_num-1] of second column.
 * @param[in] src2_row Row index array [0...src2_num-1] of second column.
 * @param[in] src2_num Number of elements in second column.
 * @param[in] dst_data Data array [0...dst_num-1] of result column.
 * @param[in] dst_row Row index array [0...dst_num-1] of result column.
 * @param[in] dst_num Number of elements in result column.
 *
 * This method mixes two sparse matrix columns into a single column.
 ***************************************************************************/
void GMatrixSparse::mix_column(const double* src1_data, const int* src1_row,
                               int src1_num,
                               const double* src2_data, const int* src2_row,
                               int src2_num,
                               double* dst_data, int* dst_row, int* dst_num)
{
    // Initialise indices and row indices of both columns
    int inx_1 = 0;                    // Column 1 element index
    int inx_2 = 0;                    // Column 2 element index
    int inx   = 0;                    // Result column element index
    int row_1 = src1_row[inx_1];
    int row_2 = src2_row[inx_2];

    // Mix elements of both columns while both contain still elements
    while (inx_1 < src1_num && inx_2 < src2_num) {

        // Case A: the element exists in both columns, so we add up the values
        if (row_1 == row_2) {
            dst_data[inx] = src1_data[inx_1] + src2_data[inx_2];
            dst_row[inx]  = row_1;
            row_1         = src1_row[++inx_1];
            row_2         = src2_row[++inx_2];
        }

        // Case B: the element exists only in first column, so we copy the element
        // from the first column
        else if (row_1 < row_2) {
            dst_data[inx] = src1_data[inx_1];
            dst_row[inx]  = row_1;
            row_1         = src1_row[++inx_1];
        }

        // Case C: the element exists only in second column, so we copy the element
        // from the second column
        else {
            dst_data[inx] = src2_data[inx_2];
            dst_row[inx]  = row_2;
            row_2         = src2_row[++inx_2];
        }

        // Update the destination index since we added a element
        inx++;

    } // endwhile: mixing

    // At this point either the first or the second column expired of elements
    // In the case that there are still elements remaining in the first column we
    // add them now ...
    for (int i = inx_1; i < src1_num; ++i, ++inx) {
        dst_data[inx] = src1_data[i];
        dst_row[inx]  = src1_row[i];
    }

    // ... or in the case that there are still elements remaining in the second
    // column we add them now
    for (int i = inx_2; i < src2_num; ++i, ++inx) {
        dst_data[inx] = src2_data[i];
        dst_row[inx]  = src2_row[i];
    }

    // Now store the number of columns in the second column
    *dst_num = inx;

    // We're done
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
 * @param[in] matrix Matrix for which absolute values are to be returned.
 ***************************************************************************/
GMatrixSparse abs(const GMatrixSparse& matrix)
{
    // Define result matrix
    GMatrixSparse result = matrix;

    // Fill pending element
    result.fill_pending();

    // Convert all elements to absolute values
    for (int i = 0; i < result.m_elements; ++i) {
        result.m_data[i] = std::abs(result.m_data[i]);
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief cs_symperm
 *
 * @param[in] matrix Matrix.
 * @param[in] pinv TBD.
 *
 * Returns matrix(p,p) where matrix and matrix(p,p) are symmetric the upper
 * part stored.
 ***************************************************************************/
GMatrixSparse cs_symperm(const GMatrixSparse& matrix, const int* pinv)
{
    // Declare loop variables
    int i, j, p, q, i2, j2;

    // Assign matrix attributes
    int     n  = matrix.m_cols;
    int*    Ap = matrix.m_colstart;
    int*    Ai = matrix.m_rowinx; 
    double* Ax = matrix.m_data;

    // Allocate result matrix
    GMatrixSparse C(n, n, Ap[n]);

    // Allocate and initialise workspace
    int  wrk_size = n;
    int* wrk_int  = new int[wrk_size];
    for (i = 0; i < wrk_size; ++i) {
        wrk_int[i] = 0;
    }

    // Assign result matrix attributes
    int*    Cp = C.m_colstart;
    int*    Ci = C.m_rowinx; 
    double* Cx = C.m_data;

    // Count entries in each column of C
    for (j = 0; j < n; j++) {

        // Column j of A is column j2 of C
        j2 = pinv ? pinv[j] : j;

        // Loop over entries in column j
        for (p = Ap[j]; p < Ap[j+1]; p++) {
            i = Ai [p];
            if (i > j) continue;              // skip lower triangular part of A
            i2 = pinv ? pinv[i] : i;          // row i of A is row i2 of C
            wrk_int[G_MAX(i2, j2)]++;         // column count of C
        }
    }

    // Compute column pointers of C
    cs_cumsum(Cp, wrk_int, n);

    // Loop over all columns of A
    for (j = 0 ; j < n ; j++) {

        // Column j of A is column j2 of C
        j2 = pinv ? pinv[j] : j;

        // Loop over entries in column j
        for (p = Ap[j]; p < Ap[j+1]; p++) {
            i = Ai [p] ;
            if (i > j) continue;              // skip lower triangular part of A
            i2    = pinv ? pinv[i] : i;       // row i of A is row i2 of C
            Ci[q  = wrk_int[G_MAX(i2,j2)]++] = G_MIN(i2,j2);
            if (Cx) Cx[q] = Ax[p];
        }
    }

    // Free workspace
    delete [] wrk_int;
  
    // Rectify the number of elements in matrix C
    C.free_elements(Cp[n], (C.m_elements-Cp[n]));

    // Return result
    return C;
}


/***************************************************************************
 * @brief Compute transpose matrix
 *
 * @param[in] matrix Matrix.
 * @param[in] value Flag that signals if values should be copied.
 *
 * The flag 'values' allows to avoid copying the actual data values. This
 * allows to perform a logical matrix transposition, as needed by the
 * symbolic matrix analysis class.
 *
 * Note that this method does not support pending elements (they have to be
 * filled before).
 ***************************************************************************/
GMatrixSparse cs_transpose(const GMatrixSparse& matrix, int values)
{
    // Declare and allocate result matrix 
    GMatrixSparse result(matrix.m_cols, matrix.m_rows, matrix.m_elements);

    // Allocate and initialise workspace
    int  wrk_size = matrix.m_rows;
    int* wrk_int  = new int[wrk_size];
    for (int i = 0; i < wrk_size; ++i) {
        wrk_int[i] = 0;
    }

    // Setup the number of non-zero elements in each row
    // for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;
    // Ap[n] = m.m_colstart[col]
    //     n = m.m_cols
    // Ai[p] = m.m_rowinx[p]
    for (int p = 0; p < matrix.m_colstart[matrix.m_cols]; ++p) {
        wrk_int[matrix.m_rowinx[p]]++;
    }

    // Set row pointers. To use a GSparseSymbolic function we have to
    // allocate and object (but this does not take memory)
    cs_cumsum(result.m_colstart, wrk_int, matrix.m_rows);

    // Case A: Normal transponse, including assignment of values
    if (values) {
        for (int col = 0; col < matrix.m_cols; ++col) {
            for (int p = matrix.m_colstart[col]; p < matrix.m_colstart[col+1] ; ++p) {
                int i              = wrk_int[matrix.m_rowinx[p]]++;
                result.m_rowinx[i] = col;
                result.m_data[i]   = matrix.m_data[p] ;
            }
        }
    }

    // Case B: Logical transponse, no assignment of values is performed
    else {
        for (int col = 0; col < matrix.m_cols; ++col) {
            for (int p = matrix.m_colstart[col]; p < matrix.m_colstart[col+1] ; ++p) {
                result.m_rowinx[wrk_int[matrix.m_rowinx[p]]++] = col;
            }
        }
    }
    
    // Free workspace
    delete [] wrk_int;

    // Return transponse matrix
    return result;
}


/***********************************************************************//**
 * @brief cs_cumsum
 *
 * @param[out] p Integer array (n+1 elements).
 * @param[in] c Integer array (n elements).
 * @param[in] n Number of elements.
 *
 * Evaluate p[0..n] = cumulative sum of c[0..n-1]
 ***************************************************************************/
double cs_cumsum(int* p, int* c, int n)
{
    // Signal error if one of the input pointers is NULL
    if (!p || !c) return (-1);

    // Initialise sums (integer and double)
    int    nz  = 0;
    double nz2 = 0.0;

    // Evaluate p[0..n] = cumulative sum of c[0..n-1]
    for (int i = 0; i < n; ++i) {
        p[i]  = nz ;
        nz   += c[i];
        nz2  += c[i];    // also in double to avoid int overflow
        c[i]  = p[i];    // also copy p[0..n-1] back into c[0..n-1]
    }
    p[n] = nz ;

    // Return cumulative sum of c[0..n-1]
    return nz2;
}
