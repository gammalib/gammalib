/***************************************************************************
 *                   GMatrixSparse.cpp - Sparse matrix class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2021 by Juergen Knoedlseder                         *
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
#define G_CONSTRUCTOR        "GMatrixSparse::GMatrixSparse(int&, int&, int&)"
#define G_OP_MUL_VEC                     "GMatrixSparse::operator*(GVector&)"
#define G_OP_ADD                  "GMatrixSparse::operator+=(GMatrixSparse&)"
#define G_OP_SUB                  "GMatrixSparse::operator-=(GMatrixSparse&)"
#define G_OP_MAT_MUL              "GMatrixSparse::operator*=(GMatrixSparse&)"
#define G_AT                                  "GMatrixSparse::at(int&, int&)"
#define G_EXTRACT_ROW                              "GMatrixSparse::row(int&)"
#define G_SET_ROW                        "GMatrixSparse::row(int&, GVector&)"
#define G_EXTRACT_COLUMN                        "GMatrixSparse::column(int&)"
#define G_SET_COLUMN                  "GMatrixSparse::column(int&, GVector&)"
#define G_SET_COLUMN2       "GMatrixSparse::column(int&, double*, int*, int)"
#define G_ADD_TO_COLUMN        "GMatrixSparse::add_to_column(int&, GVector&)"
#define G_ADD_TO_COLUMN2        "GMatrixSparse::add_to_column(int&, double*,"\
                                                                " int*, int)"
#define G_CHOL_DECOMP               "GMatrixSparse::cholesky_decompose(bool)"
#define G_CHOL_SOLVE         "GMatrixSparse::cholesky_solver(GVector&, bool)"
#define G_STACK_INIT                  "GMatrixSparse::stack_init(int&, int&)"
#define G_STACK_PUSH1      "GMatrixSparse::stack_push_column(GVector&, int&)"
#define G_STACK_PUSH2 "GMatrixSparse::stack_push_column(double*, int*, int&,"\
                                                                     " int&)"
#define G_STACK_FLUSH                      "GMatrixSparse::stack_flush(void)"
#define G_COPY_MEMBERS          "GMatrixSparse::copy_members(GMatrixSparse&)"
#define G_ALLOC_MEMBERS      "GMatrixSparse::alloc_members(int&, int&, int&)"
#define G_GET_INDEX                    "GMatrixSparse::get_index(int&, int&)"
#define G_ALLOC                   "GMatrixSparse::alloc_elements(int&, int&)"
#define G_FREE                     "GMatrixSparse::free_elements(int&, int&)"
#define G_REMOVE_ZERO              "GMatrixSparse::remove_zero_row_col(void)"
#define G_INSERT_ROW       "GMatrixSparse::insert_row(int&, GVector&, bool&)"
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
 * Constructs empty sparse matrix. The number of rows and columns of the
 * matrix will be zero.
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
 * @param[in] rows Number of rows [>=0].
 * @param[in] columns Number of columns [>=0].
 * @param[in] elements Number of allocated elements.
 *
 * @exception GException::invalid_argument
 *            Number of rows or columns is negative.
 *
 * Constructs sparse matrix of dimension @p rows times @p columns. The
 * optional @p elements argument specifies the size of the physical
 * memory that should be allocated. By default, no memory will be allocated
 * and memory allocation will be performed on-the-fly. If the amount of
 * required memory is larger than the size specified by @p elements,
 * additional momeory will be allocated automatically on-the-fly.
 ***************************************************************************/
GMatrixSparse::GMatrixSparse(const int& rows,
                             const int& columns,
                             const int& elements) :
               GMatrixBase()
{
    // Compile option: raise an exception if number of rows or columns is
    // negative
    #if defined(G_RANGE_CHECK)
    if (rows < 0) {
        std::string msg = "Number of rows "+gammalib::str(rows)+" is negative. "
                          "Please specify a non-negative number of rows.";
        throw GException::invalid_argument(G_CONSTRUCTOR, msg);
    }
    if (columns < 0) {
        std::string msg = "Number of columns "+gammalib::str(columns)+" is "
                          "negative. Please specify a non-negative number of "
                          "columns.";
        throw GException::invalid_argument(G_CONSTRUCTOR, msg);
    }
    #endif

    // Initialise private members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(rows, columns, elements);

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrix to GMatrixSparse storage class convertor
 *
 * @param[in] matrix Generic matrix (GMatrix).
 *
 * Constructs a sparse matrix by using the number of rows and columns of
 * a generic matrix and by assigning the elements of the generic matrix
 * to the sparse matrix.
 ***************************************************************************/
GMatrixSparse::GMatrixSparse(const GMatrix& matrix) : GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Construct matrix
    alloc_members(matrix.rows(), matrix.columns());

    // Fill matrix column by column
    for (int col = 0; col < m_cols; ++col) {
        GVector vector = matrix.column(col);
        this->column(col, vector);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrixSymmetric to GMatrixSparse storage class convertor
 *
 * @param[in] matrix Symmetric matrix (GMatrixSymmetric).
 *
 * Constructs a sparse matrix by using the number of rows and columns of
 * a symmetric matrix and by assigning the elements of the symmetric matrix
 * to the sparse matrix.
 ***************************************************************************/
GMatrixSparse::GMatrixSparse(const GMatrixSymmetric& matrix) : GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(matrix.rows(), matrix.columns());

    // Fill matrix column by column
    for (int col = 0; col < m_cols; ++col) {
        GVector vector = matrix.column(col);
        this->column(col, vector);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] matrix Matrix.
 *
 * Constructs matrix by copying an existing matrix.
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
 * @brief Matrix assignment operator
 *
 * @param[in] matrix Matrix.
 * @return Matrix.
 *
 * Assigns the content of another matrix to the actual matrix instance.
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
 * @brief Value assignment operator
 *
 * @param[in] value Value.
 * @return Matrix.
 *
 * Assigns the specified @p value to all elements of the matrix.
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator=(const double& value)
{
    // Continue only if matrix is not empty
    if (!is_empty()) {

        // Fill any pending element to have a non-pending state
        fill_pending();

        // If value is 0 then simply reinitialize column start indices
        if (value == 0) {

            // Initialise column start indices to 0
            for (int col = 0; col <= m_cols; ++col) {
                m_colstart[col] = 0;
            }

        }

        // ... otherwise fill column-wise
        else {

            // Set column vector
            GVector column(m_rows);
            column = value;

            // Column-wise setting
            for (int col = 0; col < m_cols; ++col) {
                this->column(col, column);
            }

        }

    } // endif: matrix was not empty

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Return reference to matrix element
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] column Matrix column [0,...,columns()-1].
 * @return Reference to matrix element.
 *
 * Returns the reference to the matrix element at @p row and @p column. If
 * the matrix element does not yet exist, a reference to the pending element
 * with a value of 0.0 is returned.
 *
 * If the matrix row or column are invalid an exception is thrown by the
 * get_index() method.
 ***************************************************************************/
double& GMatrixSparse::operator()(const int& row, const int& column)
{
    // Fill pending element. This will set the value of the pending element
    // to 0.0.
    fill_pending();

    // Get element
    int inx = get_index(row, column);
    double* value;
    if (inx < 0) {
        value      = &m_fill_val;
        m_fill_row = row;
        m_fill_col = column;
    }
    else {
        value = &(m_data[inx]);
    }

    // Return element
    return *value;
}


/***********************************************************************//**
 * @brief Return reference to matrix element (const version)
 *
 * @param[in] row Matrix row [0,...,rows()-1].
 * @param[in] column Matrix column [0,...,columns()-1].
 * @return Const reference to matrix element.
 *
 * Returns a const reference to the matrix element at @p row and @p column.
 * If the matrix element does not yet exist, a reference to the zero element
 * is returned. If the matrix element corresponds to the pending element,
 * a reference to the pending element is returned. Otherwise, a reference
 * to the matrix elements is returned.
 *
 * If the matrix row or column are invalid an exception is thrown by the
 * get_index() method.
 ***************************************************************************/
const double& GMatrixSparse::operator()(const int& row,
                                        const int& column) const
{
    // Get element. We need here the zero element to return also a pointer
    // for 0.0 entry that is not stored. Since we have the const version we
    // don't have to care about modification of this zero value.
    int inx = get_index(row, column);
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
 * @exception GException::invalid_argument
 *            Vector size differs from number of columns in matrix.
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
    // Throw exception if the matrix and vector dimensions are not compatible
    if (m_cols != vector.size()) {
        std::string msg = "Vector size "+gammalib::str(vector.size())+" does "
                          "not match "+gammalib::str(m_cols)+" matrix columns. "
                          "Please specify a vector of size "+
                          gammalib::str(m_cols)+".";
        throw GException::invalid_argument(G_OP_MUL_VEC, msg);
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
 * @brief Negate matrix elements
 *
 * @return Matrix with negated elements.
 ***************************************************************************/
GMatrixSparse GMatrixSparse::operator-(void) const
{
    // Copy matrix
    GMatrixSparse matrix = *this;

    // Fill pending element
    matrix.fill_pending();

    // Negate all matrix elements
    double* ptr = matrix.m_data;
    for (int i = 0; i < matrix.m_elements; ++i, ++ptr) {
        *ptr = -(*ptr);
    }

    // Return matrix
    return matrix;
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
 * @exception GException::invalid_argument
 *            Incompatible number of matrix rows or columns.
 *
 * This method performs a matrix addition. The operation can only succeed
 * when the dimensions of both matrices are identical.
 *
 * @todo Implement native sparse code
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator+=(const GMatrixSparse& matrix)
{
    // Throw an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows) {
        std::string msg = "Number of matrix rows "+gammalib::str(matrix.m_rows)+
                          " differs from the expected number of "+
                          gammalib::str(m_rows)+" rows. Please specify a "
                          "matrix with "+gammalib::str(m_rows)+" rows.";
        throw GException::invalid_argument(G_OP_ADD, msg);
    }
    if (m_cols != matrix.m_cols) {
        std::string msg = "Number of matrix columns "+gammalib::str(matrix.m_cols)+
                          " differs from the expected number of "+
                          gammalib::str(m_cols)+" columns. Please specify a "
                          "matrix with "+gammalib::str(m_cols)+" columns.";
        throw GException::invalid_argument(G_OP_ADD, msg);
    }

    // Perform inplace matrix addition using vectors
    for (int col = 0; col < m_cols; ++col) {
        GVector v_result  = column(col);
        GVector v_operand = matrix.column(col);
        v_result += v_operand;
        column(col, v_result);
    }

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Unary matrix scalar addition operator
 *
 * @param[in] scalar Scalar.
 *
 * Adds a @p scalar to each matrix element.
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator+=(const double& scalar)
{
    // Perform inplace matrix addition using vectors
    for (int col = 0; col < m_cols; ++col) {
        GVector vector = column(col);
        vector += scalar;
        column(col, vector);
    }

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Unary matrix subtraction operator
 *
 * @param[in] matrix Matrix.
 *
 * @exception GException::invalid_argument
 *            Incompatible matrix size.
 *
 * This method performs a matrix addition. The operation can only succeed
 * when the dimensions of both matrices are identical.
 *
 * @todo Implement native sparse code
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator-=(const GMatrixSparse& matrix)
{
    // Throw an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows) {
        std::string msg = "Number of matrix rows "+gammalib::str(matrix.m_rows)+
                          " differs from the expected number of "+
                          gammalib::str(m_rows)+" rows. Please specify a "
                          "matrix with "+gammalib::str(m_rows)+" rows.";
        throw GException::invalid_argument(G_OP_SUB, msg);
    }
    if (m_cols != matrix.m_cols) {
        std::string msg = "Number of matrix columns "+gammalib::str(matrix.m_cols)+
                          " differs from the expected number of "+
                          gammalib::str(m_cols)+" columns. Please specify a "
                          "matrix with "+gammalib::str(m_cols)+" columns.";
        throw GException::invalid_argument(G_OP_SUB, msg);
    }

    // Perform inplace matrix subtraction
    for (int col = 0; col < m_cols; ++col) {
        GVector v_result  = column(col);
        GVector v_operand = matrix.column(col);
        v_result -= v_operand;
        column(col, v_result);
    }

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Unary matrix scalar subtraction operator
 *
 * @param[in] scalar Scalar.
 *
 * Subtracts a @p scalar from each matrix element.
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator-=(const double& scalar)
{
    // Perform inplace matrix subtraction
    for (int col = 0; col < m_cols; ++col) {
        GVector vector = column(col);
        vector -= scalar;
        column(col, vector);
    }

    // Return result
    return *this;
}


/***********************************************************************//**
 * @brief Unary matrix multiplication operator
 *
 * @param[in] matrix Matrix.
 *
 * @exception GException::invalid_argument
 *            Number of rows in @p matrix is different from number of
 *            matrix columns.
 *
 * This method performs a matrix multiplication. The operation can only
 * succeed when the dimensions of both matrices are compatible.
 *
 * @todo Implement native sparse code
 ***************************************************************************/
GMatrixSparse& GMatrixSparse::operator*=(const GMatrixSparse& matrix)
{
    // Throw exception if the matrix dimensions are not compatible
    if (m_cols != matrix.m_rows) {
        std::string msg = "Number of "+gammalib::str(m_cols)+" columns in "
                          "the first matrix differs from number of "+
                          gammalib::str(matrix.m_rows)+" rows in the second "
                          "matrix. Please specify a second matrix with "+
                          gammalib::str(m_cols)+" rows.";
        throw GException::invalid_argument(G_OP_MAT_MUL, msg);
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
 * @brief Clear matrix
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
 * @brief Clone matrix
 *
 * @return Pointer to deep copy of matrix.
 ***************************************************************************/
GMatrixSparse* GMatrixSparse::clone(void) const
{
    // Clone matrix
    return new GMatrixSparse(*this);
}


/***********************************************************************//**
 * @brief Return reference to matrix element
 *
 * @param[in] row Matrix row [0,...,rows()[.
 * @param[in] column Matrix column [0,...,columns()[.
 * @return Reference to matrix element.
 *
 * @exception GException::out_of_range
 *            Matrix row or column out of range.
 ***************************************************************************/
double& GMatrixSparse::at(const int& row, const int& column)
{
    // Throw exception if row or column index is out of range
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_AT, "Matrix row", row, m_rows);
    }
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_AT, "Matrix column", column, m_cols);
    }

    // Return element
    return ((*this)(row, column));
}


/***********************************************************************//**
 * @brief Return reference to matrix element (const version)
 *
 * @param[in] row Matrix row [0,...,rows()[.
 * @param[in] column Matrix column [0,...,columns()[.
 * @return Reference to matrix element.
 *
 * @exception GException::out_of_range
 *            Matrix row or column out of range.
 ***************************************************************************/
const double& GMatrixSparse::at(const int& row, const int& column) const
{
    // Throw exception if row or column index is out of range
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_AT, "Matrix row", row, m_rows);
    }
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_AT, "Matrix column", column, m_cols);
    }

    // Return element
    return ((*this)(row, column));
}


/***********************************************************************//**
 * @brief Extract row as vector from matrix
 *
 * @param[in] row Matrix row [0,...,rows()[.
 * @return Vector of matrix row elements (columns() elements).
 *
 * @exception GException::out_of_range
 *            Invalid matrix row specified.
 *
 * This method extracts a matrix row into a vector.
 ***************************************************************************/
GVector GMatrixSparse::row(const int& row) const
{
    // Throw exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_EXTRACT_ROW, "Matrix row", row, m_rows);
    }
    #endif

    // Create result vector
    GVector result(m_cols);

    // Loop over all columns to extract data
    for (int col = 0; col < m_cols; ++col) {

        // Get the start and stop of the elements
        int i_start = m_colstart[col];
        int i_stop  = m_colstart[col+1];

        // If column is not empty then search requested row by bisection
        if (i_stop > i_start) {

            // Find row by bisection
            while ((i_stop - i_start) > 1) {
                int mid = (i_start+i_stop) / 2;
                if (m_rowinx[mid] > row) {
                    i_stop = mid;
                }
                else {
                    i_start = mid;
                }
            }

            // Copy element if we found the row
            if (m_rowinx[i_start] == row) {
                result[col] = m_data[i_start];
            }

        } // endif: column was not empty

    } // endfor: looped over all columns

    // If there is a pending element then put it in the vector
    if (m_fill_val != 0.0 && m_fill_row == row) {
        result[m_fill_col] = m_fill_val;
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Set row in matrix
 *
 * @param[in] row Matrix row [0,...,rows()[.
 * @param[in] vector Vector of matrix row elements (columns() elements).
 ***************************************************************************/
void GMatrixSparse::row(const int& row, const GVector& vector)
{
    // Insert row
    insert_row(row, vector, false);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Extract column as vector from matrix
 *
 * @param[in] column Matrix column [0,...,columns()[.
 * @return Vector of matrix column elements (rows() elements).
 *
 * @exception GException::out_of_range
 *            Invalid matrix column specified.
 *
 * This method extracts a matrix column into a vector.
 ***************************************************************************/
GVector GMatrixSparse::column(const int& column) const
{
    // Throw exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_EXTRACT_COLUMN, "Matrix column",
                                       column, m_cols);
    }
    #endif

    // Create result vector
    GVector result(m_rows);

    // Get the start and stop of the elements
    int i_start = m_colstart[column];
    int i_stop  = m_colstart[column+1];

    // Extract elements into vector
    for (int i = i_start; i < i_stop; ++i) {
        result[m_rowinx[i]] = m_data[i];
    }

    // If there is a pending element then put it in the vector
    if (m_fill_val != 0.0 && m_fill_col == column) {
        result[m_fill_row] = m_fill_val;
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Insert vector column into matrix
 *
 * @param[in] column Matrix column [0,...,columns()[.
 * @param[in] vector Vector.
 *
 * @exception GException::out_of_range
 *            Invalid matrix column specified.
 * @exception GException::invalid_argument
 *            Vector size does not match number of matrix rows.
 *
 * Inserts the content of a vector into a matrix column. Any previous
 * content in the matrix column will be overwritted.
 *
 * This is the main driver routine to insert data into a matrix. Note that
 * there is another instance of this function that takes a compressed array.
 ***************************************************************************/
void GMatrixSparse::column(const int& column, const GVector& vector)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_INSERTION)
    std::cout << "GMatrixSparse::column(";
    std::cout << column << ", [" << vector << "]):" << std::endl;
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

    // Throw exception if the column is invalid
    #if defined(G_RANGE_CHECK)
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_SET_COLUMN, "Matrix column",
                                       column, m_cols);
    }
    #endif

    // Throw exception if the vector size does not match the number of rows
    // in the matrix
    if (m_rows != vector.size()) {
        std::string msg = "Vector size "+gammalib::str(vector.size())+
                          " does not match "+gammalib::str(m_rows)+" matrix "
                          "rows. Please specify a vector of size "+
                          gammalib::str(m_rows)+".";
        throw GException::invalid_argument(G_SET_COLUMN, msg);
    }

    // If there is a pending element for this column then delete it since
    // the vector overwrites this element
    if (m_fill_val != 0.0 && m_fill_col == column) {
        #if defined(G_DEBUG_SPARSE_PENDING)
        std::cout << G_SET_COLUMN << ": pending value " << m_fill_val << 
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
    int i_start = m_colstart[column];
    int i_stop  = m_colstart[column+1];
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
    for (int i = column+1; i <= m_cols; ++i) {
        m_colstart[i] += n_diff;
    }

    // Debugging: show sparse matrix after insertion
    #if defined(G_DEBUG_SPARSE_INSERTION)
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
 * @brief Insert compressed array into matrix column
 *
 * @param[in] column Matrix column [0,...,columns()[.
 * @param[in] values Compressed array.
 * @param[in] rows Row numbers of array.
 * @param[in] number Number of elements in array.
 *
 * @exception GException::out_of_range
 *            Invalid matrix column specified.
 * @exception GException::invalid_argument
 *            Row numbers are not smaller than number of matrix rows
 *
 * Inserts the content of a copressed array into a matrix column. Any 
 * previous content in the matrix column will be overwritted.
 *
 * This is the main driver routine to insert data into a matrix. Note that
 * there is another instance of this function that takes a vector.
 ***************************************************************************/
void GMatrixSparse::column(const int& column, const double* values,
                           const int* rows, int number)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_INSERTION)
    std::cout << "GMatrixSparse::column(";
    std::cout << column << ", values, rows, " << number << "):" << std::endl;
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

    // Throw exception if the column is invalid
    #if defined(G_RANGE_CHECK)
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_SET_COLUMN2, "Matrix column",
                                       column, m_cols);
    }
    #endif

    // Throw exception if the index array seems incompatible with matrix
    // dimensions
    if (rows[number-1] >= m_rows) {
        std::string msg = "Row number "+gammalib::str(rows[number-1])+" for "
                          "element "+gammalib::str(number)+" is not smaller "
                          "than number "+gammalib::str(m_rows)+" of matrix "
                          "rows. Please specify only row numbers smaller than "+
                          gammalib::str(m_rows)+".";
        throw GException::invalid_argument(G_SET_COLUMN2, msg);
    }

    // If there is a pending element for this column then delete it since
    // the vector overwrites this element
    if (m_fill_val != 0.0 && m_fill_col == column) {
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
    int i_start = m_colstart[column];
    int i_stop  = m_colstart[column+1];
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
    for (int i = column+1; i <= m_cols; ++i) {
        m_colstart[i] += n_diff;
    }

    // Debugging: show sparse matrix after insertion
    #if defined(G_DEBUG_SPARSE_INSERTION)
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
 * @brief Add row to matrix elements
 *
 * @param[in] row Matrix row [0,...,row()[.
 * @param[in] vector Vector.
 ***************************************************************************/
void GMatrixSparse::add_to_row(const int& row, const GVector& vector)
{
    // Add row
    insert_row(row, vector, true);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Add vector column into matrix
 *
 * @param[in] column Matrix column [0,...,columns()[.
 * @param[in] vector Vector.
 *
 * @exception GException::out_of_range
 *            Invalid matrix column specified.
 * @exception GException::invalid_argument
 *            Matrix dimension mismatches the vector size.
 *
 * Adds the contect of a vector to a matrix column.
 *
 * This is the main driver routine to add data to a matrix. It handles both
 * normal and stack-based filled. Note that there is another instance of this
 * method that takes a compressed array.
 ***************************************************************************/
void GMatrixSparse::add_to_column(const int& column, const GVector& vector)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_ADDITION)
    std::cout << "GMatrixSparse::add_col(";
    std::cout << column << ", [" << v << "]):" << std::endl;
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
        non_zero = stack_push_column(vector, column);
        if (non_zero == 0) {
            return;
        }
    }

    // ... otherwise check first the arguments and determine the number of
    // non-zero elements in vector
    else {

        // Throw exception if the column is invalid
        #if defined(G_RANGE_CHECK)
        if (column < 0 || column >= m_cols) {
            throw GException::out_of_range(G_ADD_TO_COLUMN, "Matrix column",
                                           column, m_cols);
        }
        #endif

        // Throw exception if the vector size does not match the number of rows
        // in the matrix
        if (m_rows != vector.size()) {
            std::string msg = "Vector size "+gammalib::str(vector.size())+
                              " does not match "+gammalib::str(m_rows)+" matrix "
                              "rows. Please specify a vector of size "+
                              gammalib::str(m_rows)+".";
            throw GException::invalid_argument(G_ADD_TO_COLUMN, msg);
        }

        // Determine number of non-zero elements in vector
        non_zero = vector.non_zeros();
    }

    // Extract vector for column, add elements, and re-insert vector (only if
    // vector to insert has non-zeros)
    if (non_zero > 0) {

        // Copy input vector
        GVector v_column = vector;

        // Add elements to vector
        for (int i = m_colstart[column]; i < m_colstart[column+1]; ++i) {
            v_column[m_rowinx[i]] += m_data[i];
        }

        // If there is a pending element then put it in the vector
        if (m_fill_val != 0.0 && m_fill_col == column) {
            v_column[m_fill_row] += m_fill_val;
        }

        // Set vector into matrix
        this->column(column, v_column);

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
 * @param[in] column Matrix column [0,...,columns()[.
 * @param[in] values Compressed array.
 * @param[in] rows Row indices of array.
 * @param[in] number Number of elements in array.
 *
 * @exception GException::out_of_range
 *            Invalid matrix column specified.
 * @exception GException::invalid_argument
 *            Matrix dimension mismatches the vector size.
 *
 * Adds the content of a compressed array into a matrix column.
 *
 * This is the main driver routine to add data to a matrix. It handles both
 * normal and stack-based filled. Note that there is another instance of this
 * method that takes a vector.
 ***************************************************************************/
void GMatrixSparse::add_to_column(const int& column, const double* values,
                                  const int* rows, int number)
{
    // Debug header
    #if defined(G_DEBUG_SPARSE_ADDITION)
    std::cout << "GMatrixSparse::add_col(";
    std::cout << column << ", values, rows, " << number << "):" << std::endl;
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
        number = stack_push_column(values, rows, number, column);
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

        // Throw exception if the column is invalid
        #if defined(G_RANGE_CHECK)
        if (column < 0 || column >= m_cols) {
            throw GException::out_of_range(G_ADD_TO_COLUMN2, "Matrix column",
                                           column, m_cols);
        }
        #endif

        // Throw exception if the index array seems incompatible with matrix
        // dimensions
        if (rows[number-1] >= m_rows) {
            std::string msg = "Row number "+gammalib::str(rows[number-1])+" for "
                              "element "+gammalib::str(number)+" is not smaller "
                              "than number "+gammalib::str(m_rows)+" of matrix "
                              "rows. Please specify only row numbers smaller than "+
                              gammalib::str(m_rows)+".";
            throw GException::invalid_argument(G_ADD_TO_COLUMN2, msg);
        }

    } // endelse: there was no stack

    // Get indices of column in matrix
    int i_start = m_colstart[column];
    int i_stop  = m_colstart[column+1];

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
        this->column(column, wrk_double, wrk_int, num_mix);

        // Free workspace
        delete [] wrk_int;
        delete [] wrk_double;

    } // endif: Case A

    // Case B: the column does not yet exist in the matrix, so just insert it
    else {
        this->column(column, values, rows, number);
    }

    // Debugging: show sparse matrix after insertion
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
 * @brief Return transposed matrix
 *
 * @return Transposed matrix.
 *
 * Returns transposed matrix of the matrix.
 ***************************************************************************/
GMatrixSparse GMatrixSparse::transpose(void) const
{
    // Copy matrix
    GMatrixSparse matrix = *this;

    // Fill pending element
    matrix.fill_pending();

    // Compute the transpose
    matrix = cs_transpose(matrix, 1);

    // Return matrix
    return matrix;
}


/***********************************************************************//**
 * @brief Return inverted matrix
 *
 * @return Inverted matrix.
 *
 * Returns inverse of matrix. Inversion is done for the moment using Cholesky
 * decomposition. This does not work on any kind of matrix.
 *
 * @todo Specify in documentation for which kind of matrix the method works.
 ***************************************************************************/
GMatrixSparse GMatrixSparse::invert(void) const
{
    // Copy matrix
    GMatrixSparse matrix = *this;

    // Fill pending element
    matrix.fill_pending();

    // Invert matrix
    matrix = matrix.cholesky_invert(true);

    // Return matrix
    return matrix;
}


/***********************************************************************//**
 * @brief Solves linear matrix equation
 *
 * @param[in] vector Solution vector.
 * @return Solutions of matrix equation.
 * 
 * Solves the linear equation
 *
 * \f[M \times {\tt solution} = {\tt vector} \f]
 *
 * where \f$M\f$ is the matrix, \f${\tt vector}\f$ is the result, and
 * \f${\tt solution}\f$ is the solution. Solving is done using Cholesky
 * decomposition. This does not work on any kind of matrix.
 *
 * If the matrix is empty and the vector has a zero length the method
 * returns an empty vector.
 *
 * @todo Specify in documentation for which kind of matrix the method works.
 ***************************************************************************/
GVector GMatrixSparse::solve(const GVector& vector) const
{
    // Initialise result with an empty vector
    GVector result;

    // Continue only if matrix is not empty or vector size is not zero
    if (!is_empty() || vector.size() != 0) {

        // Get Cholesky decomposition of matrix
        GMatrixSparse decomposition = cholesky_decompose(true);

        // Solve linear equation
        result = decomposition.cholesky_solver(vector);

    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Return absolute of matrix
 *
 * @return Absolute of matrix
 *
 * Returns matrix where all elements of the matrix have been replaced by
 * their absolute values.
 ***************************************************************************/
GMatrixSparse GMatrixSparse::abs(void) const
{
    // Copy matrix
    GMatrixSparse matrix = *this;

    // Fill pending element
    matrix.fill_pending();

    // Take the absolute value of all matrix elements
    double* ptr = matrix.m_data;
    for (int i = 0; i < matrix.m_elements; ++i, ++ptr) {
        *ptr = std::abs(*ptr);
    }

    // Return matrix
    return matrix;
}


/***********************************************************************//**
 * @brief Returns fill of matrix
 *
 * The fill of a matrix is defined as the number non-zero elements devided
 * by the number of total elements. By definition, the fill is comprised
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
 * @brief Return Cholesky decomposition
 *
 * @param[in] compress Use zero-row/column compression (defaults to true).
 * @return Cholesky decomposition of matrix
 *
 * Returns the Cholesky decomposition of a sparse matrix. The decomposition
 * is stored within a GMatrixSparse object.
 ***************************************************************************/
GMatrixSparse GMatrixSparse::cholesky_decompose(const bool& compress) const
{
    // Create copy of matrix
    GMatrixSparse matrix = *this;

    // Continue only if matrix is not empty
    if (matrix.m_rows > 0 && matrix.m_cols > 0) {

        // Save original matrix size
        int matrix_rows = matrix.m_rows;
        int matrix_cols = matrix.m_cols;

        // Delete any existing symbolic and numeric analysis object and reset
        // pointers
        if (matrix.m_symbolic != NULL) delete matrix.m_symbolic;
        if (matrix.m_numeric  != NULL) delete matrix.m_numeric;
        matrix.m_symbolic = NULL;
        matrix.m_numeric  = NULL;

        // Allocate symbolic analysis object
        GSparseSymbolic* symbolic = new GSparseSymbolic();

        // Declare numeric analysis object. We don't allocate one since we'll
        // throw it away at the end of the function (the L matrix will be copied
        // in this object)
        GSparseNumeric numeric;

        // Fill pending element into matrix
        matrix.fill_pending();

        // Remove rows and columns containing only zeros if matrix compression
        // has been selected
        if (compress) {
            matrix.remove_zero_row_col();
        }

        // Ordering an symbolic analysis of matrix. This sets up an array 'pinv'
        // which contains the fill-in reducing permutations
        symbolic->cholesky_symbolic_analysis(1, matrix);

        // Store symbolic pointer in sparse matrix object
        matrix.m_symbolic = symbolic;

        // Perform numeric Cholesky decomposition
        numeric.cholesky_numeric_analysis(matrix, *symbolic);

        // Copy L matrix into this object
        matrix.free_elements(0, matrix.m_elements);
        matrix.alloc_elements(0, numeric.m_L->m_elements);
        for (int i = 0; i < matrix.m_elements; ++i) {
            matrix.m_data[i]   = numeric.m_L->m_data[i];
            matrix.m_rowinx[i] = numeric.m_L->m_rowinx[i];
        }
        for (int col = 0; col <= matrix.m_cols; ++col) {
            matrix.m_colstart[col] = numeric.m_L->m_colstart[col];
        }

        // Insert zero rows and columns if they have been removed previously.
        if (compress) {
            matrix.insert_zero_row_col(matrix_rows, matrix_cols);
        }

    } // endif: matrix was not empty

    // Return matrix
    return matrix;
}


/***********************************************************************//**
 * @brief Cholesky solver
 *
 * @param[in] vector Solution vector.
 * @param[in] compress Request matrix compression.
 *
 * @exception GException::invalid_argument
 *            Matrix and vector do not match.
 * @exception GException::invalid_value
 *            Matrix has not been factorised.
 * 
 * Solves the linear equation A*x=b using a Cholesky decomposition of A.
 * This function is to be applied on a GMatrixSparse matrix for which a
 * Choleksy factorization has been produced using 'cholesky_decompose'.
 ***************************************************************************/
GVector GMatrixSparse::cholesky_solver(const GVector& vector,
                                       const bool& compress) const
{
    // Dump header
    #if defined(G_DEBUG_SPARSE_COMPRESSION)
    std::cout << "GMatrixSparse::cholesky_solver" << std::endl;
    std::cout << " Input vector .....: " << vector << std::endl;
    #endif

    // Throw exception if the matrix and vector dimensions are incompatible
    if (m_rows != vector.size()) {
        std::string msg = "Vector size "+gammalib::str(vector.size())+
                          " does not match "+gammalib::str(m_rows)+" matrix "
                          "rows. Please specify a vector of size "+
                          gammalib::str(m_rows)+".";
        throw GException::invalid_argument(G_CHOL_SOLVE, msg);
    }

    // Throw exception if there is no symbolic pointer
    if (!m_symbolic) {
        std::string msg = "Matrix not factorised. Please call method "
                          "GMatrixSparse::cholesky_decompose() before calling "
                          "this method.";
        throw GException::invalid_value(G_CHOL_SOLVE, msg);
    }

    // Throw exception if there is no permutation
    if (!m_symbolic->m_pinv) {
        std::string msg = "Matrix not factorised. Please call method "
                          "GMatrixSparse::cholesky_decompose() before calling "
                          "this method.";
        throw GException::invalid_value(G_CHOL_SOLVE, msg);
    }

    // Flag row and column compression
    int row_compressed = (m_rowsel != NULL && m_num_rowsel < m_rows);
    int col_compressed = (m_colsel != NULL && m_num_colsel < m_cols);

    // Decide if we need a compression algorithm or not
    int no_zero = !(compress && (row_compressed || col_compressed));

    // Allocate vector for permutation and result vector
    GVector result(m_cols);

    // Continue only if matrix is not empty
    if (m_rows > 0 && m_cols > 0) {

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
                for (int p = Lp[col]+1; p < Lp[col+1]; p++) {   // loop over elements
                    x[Li[p]] -= Lx[p] * x[col];
                }
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

            // Setup column mapping array that maps original matrix column into
            // compressed matrix columns. An entry of -1 indicates that the
            // column should be dropped. If no selection exists then setup an
            // identity map.
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

            // Inplace solve L\x=x. The column and row maps are just use to see
            // which columns or rows should be skipped in the calculations.
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

            // Inplace solve L'\x=x. The column and row maps are just use to see
            // which columns or rows should be skipped in the calculations.
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

            // If column compression has been performed the expand the result
            // vector accordingly
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

    } // endif: matrix was not empty

    // Return result vector
    return result;
}


/***********************************************************************//**
 * @brief Invert matrix using a Cholesky decomposition
 *
 * @param[in] compress Use zero-row/column compression.
 * @return Inverted matrix.
 *
 * Inverts the matrix using a Cholesky decomposition.
 ***************************************************************************/
GMatrixSparse GMatrixSparse::cholesky_invert(const bool& compress) const
{
    // Generate Cholesky decomposition of matrix
    GMatrixSparse decomposition = cholesky_decompose(compress);

    // Allocate result matrix and unit vector
    GMatrixSparse matrix(m_rows, m_cols);
    GVector       unit(m_rows);

    // Column-wise solving of the problem
    for (int col = 0; col < m_cols; ++col) {

        // Set unit vector
        unit[col] = 1.0;

        // Solve for column
        GVector x = decomposition.cholesky_solver(unit, compress);

        // Insert column in matrix
        matrix.column(col, x);

        // Clear unit vector for next round
        unit[col] = 0.0;

    }

    // Return matrix
    return matrix;
}


/***********************************************************************//**
 * @brief Print matrix
 *
 * @param[in] chatter Chattiness.
 * @return String containing matrix information
 ***************************************************************************/
std::string GMatrixSparse::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Determine number of elements
        int nonzero  = 0;
        int elements = m_elements;
        if (m_cols > 0) {
            if (m_fill_val == 0.0) {
                nonzero  = (m_colstart != NULL) ? m_colstart[m_cols] : 0;
            }
            else {
                nonzero  = (m_colstart != NULL) ? m_colstart[m_cols]+1 : 0;
                elements++;
            }
        }

        // Append header
        result.append("=== GMatrixSparse ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of rows"));
        result.append(gammalib::str(m_rows));
        if (m_rowsel != NULL) {
            result.append(" (compressed "+gammalib::str(m_num_rowsel)+")");
        }
        result.append("\n"+gammalib::parformat("Number of columns"));
        result.append(gammalib::str(m_cols));
        if (m_colsel != NULL) {
            result.append(" (compressed "+gammalib::str(m_num_colsel)+")");
        }
        result.append("\n"+gammalib::parformat("Number of nonzero elements"));
        result.append(gammalib::str(nonzero));
        if (m_fill_val != 0.0) {
            result.append("\n"+gammalib::parformat("Pending element"));
            result.append("("+gammalib::str(m_fill_row)+
                          ","+gammalib::str(m_fill_col)+")=");
            result.append(gammalib::str(m_fill_val));
        }
        result.append("\n"+gammalib::parformat("Number of allocated cells"));
        result.append(gammalib::str(m_alloc));
        result.append("\n"+gammalib::parformat("Memory block size"));
        result.append(gammalib::str(m_mem_block));
        result.append("\n"+gammalib::parformat("Sparse matrix fill"));
        result.append(gammalib::str(fill()));
        result.append("\n"+gammalib::parformat("Pending element"));
        result.append(gammalib::str(m_fill_val));
        result.append("\n"+gammalib::parformat("Fill stack size"));
        result.append(gammalib::str(m_stack_size));
        if (m_stack_data == NULL) {
            result.append(" (none)");
        }

        // Append elements and compression schemes
        result.append(print_elements(chatter));
        result.append(print_row_compression(chatter));
        result.append(print_col_compression(chatter));

    } // endif: chatter was not silent

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
 * @param[in] col Matrix column [0,...,columns()[.
 *
 * @exception GException::out_of_range
 *            Invalid matrix column specified.
 * @exception GException::invalid_argument
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
    // Raise an exception if the column is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_STACK_PUSH1, "Matrix column",
                                       col, m_cols);
    }
    #endif

    // Raise an exception if the matrix and vector dimensions are incompatible
    if (vector.size() != m_rows) {
        std::string msg = "Vector size "+gammalib::str(vector.size())+
                          " does not match "+gammalib::str(m_rows)+" matrix "
                          "rows. Please specify a vector of size "+
                          gammalib::str(m_rows)+".";
        throw GException::invalid_argument(G_STACK_PUSH1, msg);
    }

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
 * @param[in] rows Matrix rows of array elements.
 * @param[in] number Number of elements in array.
 * @param[in] col Matrix column [0,...,columns()[.
 *
 * @exception GException::out_of_range
 *            Invalid matrix column specified.
 * @exception GException::invalid_argument
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
    // Raise an exception if the column is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_STACK_PUSH2, "Matrix column",
                                       col, m_cols);
    }
    #endif

    // Initialise return value
    int remaining = number; 

    // Single loop for common exit point
    do {

        // If the array is empty there is nothing to do
        if (!values || !rows || (number < 1)) {
            continue;
        }

        // Raise an exception if the matrix and vector dimensions are
        // incompatible. This test needs to be done after the test on a
        // positive number!
        if (rows[number-1] >= m_rows) {
            std::string msg = "Row number "+gammalib::str(rows[number-1])+" for "
                              "element "+gammalib::str(number)+" is not smaller "
                              "than number "+gammalib::str(m_rows)+" of matrix "
                              "rows. Please specify only row numbers smaller than "+
                              gammalib::str(m_rows)+".";
            throw GException::invalid_argument(G_STACK_PUSH2, msg);
        }

        // If there is no stack or the stack can not hold the requested 
        // number of elements then report number of array elements to the
        // caller
        if (m_stack_data == NULL || number > m_stack_size) {
            continue;
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
 * Adds the content of the stack to the actual matrix. First, the total
 * number of matrix elements is determined. Then new memory is allocated to
 * hold all elements. Finally, all elements are filled into the new memory
 * that is then replacing the old matrix memory.
 *
 * The method uses the internal working array m_stack_work.
 *
 * NOTE: The present algorithm assumes that each column occurs only once
 * in the stack!
 ***************************************************************************/
void GMatrixSparse::stack_flush(void)
{
    // Do nothing if there is no stack, no stack entries or the matrix is
    // empty
    if (m_stack_data == NULL || m_stack_entries == 0 || is_empty()) {
        return;
    }

    // Fill pending value
    fill_pending();

    // Debug header
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    std::cout << "GMatrixSparse::stack_flush" << std::endl;
    std::cout << " Number of columns on stack : " << m_stack_entries << std::endl;
    std::cout << " Number of elements on stack: " << m_stack_start[m_stack_entries] << std::endl;
    std::cout << " Number of matrix elements .: " << m_elements << std::endl;
    std::cout << " Col.start at end of matrix : " << m_colstart[m_cols] << std::endl;
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
    int new_elements = 0;    // Number of new elements to add to matrix
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    int num_columns  = 0;    // Number of valid columns in stack
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

            // Debug option: update number of valid columns
            #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
            num_columns++;
            #endif

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
                int num_1;    // Number of elements only found in matrix column
                int num_2;    // Number of elements only found in stack column
                int num_mix;  // Number of elements found in both columns

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

            // If column did not exist in the matrix then consider all
            // elements as new
            else {
                m_stack_work[col]  = (entry+2);
                new_elements      += (m_stack_start[entry+1] - m_stack_start[entry]);
            }

        } // endif: entry was valid

    } // endfor: looped over all entries

    // Debug option: Log number of valid columns and elements in new matrix
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    std::cout << " Valid columns on stack ....: " << num_columns << std::endl;
    std::cout << " Valid elements on stack ...: " << new_elements << std::endl;
    #endif

    // Compute total number of matrix elements
    int elements = m_elements + new_elements;

    // Allocate memory for new matrix (always keep some elbow room)
    m_alloc = elements + m_mem_block;
    double* new_data   = new double[m_alloc];
    int*    new_rowinx = new int[m_alloc];

    // Fill new matrix. For this purpose we loop over all matrix columns
    // and perform the operation that was identified in the previous scan
    int index = 0;
    for (int col = 0; col < m_cols; ++col) {

        // If column does not exist then skip. Note that the index is stored
        // in [col] instead of [col+1] to not disturb the existing matrix.
        // The column start indices are corrected later.
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

        // If column exists in the matrix and in the stack then we have to
        // mix both
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

        // Store actual index in column start array. Note that the index is
        // stored in [col] instead of [col+1] to not disturb the existing
        // matrix. The column start indices are corrected later.
        m_colstart[col] = index;

    } // endfor: looped over all columns

    // Dump number of elements in new matrix after stack flushing
    #if defined(G_DEBUG_SPARSE_STACK_FLUSH)
    std::cout << " Added elements ............: " << index;
    std::cout << " (should be " << elements << ")" << std::endl;
    std::cout << " - Matrix only .............: " << num_matrix << std::endl;
    std::cout << " - Stack only ..............: " << num_stack << std::endl;
    std::cout << " - Matrix & Stack ..........: " << num_both << std::endl;
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
 * @param[in] rows Number of rows.
 * @param[in] columns Number of columns.
 * @param[in] elements Number of matrix elements to be physically allocated.
 *
 * This is the main constructor code that allocates and initialises memory
 * for matrix elements. The method only allocates elements if both @p rows
 * and @p columns are positive. Otherwise the method does nothing and will
 * set the m_rows and m_cols attributes to zero.
 ***************************************************************************/
void GMatrixSparse::alloc_members(const int& rows, const int& columns,
                                  const int& elements)
{
    // Continue only if rows and columns are valid
    if (rows > 0 && columns > 0) {

        // Free any existing memory
        if (m_colstart != NULL) delete [] m_colstart;

        // Allocate column start array. This is the only array that we can
        // allocate at this time. The other arrays can only be allocated
        // during filling of the matrix
        m_colstart = new int[columns+1];

        // Store (logical) matrix size
        m_rows = rows;
        m_cols = columns;

        // Initialise column start indices to 0
        for (int col = 0; col <= m_cols; ++col) {
            m_colstart[col] = 0;
        }

        // Optionally allocate memory for matrix elements
        if (elements > 0) {
            alloc_elements(0, elements);
        }

    } // endif: number of elements was positive

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
 * @brief Determines element index for (row,column)
 *
 * @param[in] row Matrix row [0,...,m_rows[.
 * @param[in] column Matrix column [0,...,m_cols[.
 * @return Element index.
 *
 * @exception GException::out_of_range
 *            Matrix row or column out of range
 *
 * Returns the index in the compressed array for (row,col). The following
 * special results exist:
 *
 *       -1: Requested index does not exist in the matrix.
 *       m_elements: Requested index is the pending element.
 ***************************************************************************/
int GMatrixSparse::get_index(const int& row, const int& column) const
{
    // Raise an exception if the row or column index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_GET_INDEX, "Matrix row", row, m_rows);
    }
    if (column < 0 || column >= m_cols) {
        throw GException::out_of_range(G_GET_INDEX, "Matrix column", column, m_cols);
    }
    #endif

    // Initialise element to 'not found'
    int index = -1;

    // If we have a pending element and the pending element is requested
    // then return m_elements as index
    if ((m_fill_val != 0.0) && (row == m_fill_row && column == m_fill_col)) {
        index = m_elements;
    }

    // ... otherwise if there are elements in the matrix then search for
    // the requested row and column.
    else if (m_elements > 0) {

        // Get the start and stop of the elements
        int* ptr_colstart = m_colstart + column;
        int  i_start      = *ptr_colstart++;
        int  i_stop       = *ptr_colstart;

        // If column is not empty then search requested row by bisection
        if (i_stop > i_start) {

            // Find row by bisection
            while ((i_stop - i_start) > 1) {
                int mid = (i_start+i_stop) / 2;
                if (m_rowinx[mid] > row) {
                    i_stop = mid;
                }
                else {
                    i_start = mid;
                }
            }

            // If we found the row then store the index
            if (m_rowinx[i_start] == row) {
                index = i_start;
            }

        } // endif: column was not empty

    } // endif: matrix contained elements

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
        for (int i = 0; i < m_elements; ++i) {
            std::cout << m_data[i] << " ";
        }
        std::cout << std::endl << " Row.: ";
        for (int i = 0; i < m_elements; ++i) {
            std::cout << m_rowinx[i] << " ";
        }
        std::cout << std::endl << " Col.: ";
        for (int i = 0; i < m_cols+1; ++i) {
            std::cout << m_colstart[i] << " ";
        }
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
 * @exception GException::runtime_error
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

    // Throw exception if all matrix elements are zero
    if (m_num_rowsel < 1 || m_num_colsel < 1) {
        std::string msg = "All matrix elements are zero.";
        throw GException::runtime_error(G_REMOVE_ZERO, msg);
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
        for (int i = i_start; i < i_stop; ++i) {
            int c_row = row_map[m_rowinx[i]];
            if (c_row >= 0) {
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
 * column. It counts the number of elements in two columns that have the
 * same row and puts this number into @p num_mix. The number of elements
 * that are only present in the first column will be put in @p num_1, the
 * number of elements that are only present in the second column will be put
 * in @p num_2.
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
    int inx_1 = 0;                  // Column 1 element index
    int inx_2 = 0;                  // Column 2 element index
    int row_1 = src1_row[inx_1];    // Column 1 first row index
    int row_2 = src2_row[inx_2];    // Column 2 first row index

    // Count number of elements with same row or elements that exist only in
    // either the first or the second column while both columns contain still
    // elements
    while (inx_1 < src1_num && inx_2 < src2_num) {

        // Case A: the element exist in both columns
        if (row_1 == row_2) {
            (*num_mix)++;
            inx_1++;
            inx_2++;
            if (inx_1 < src1_num) {
                row_1 = src1_row[inx_1];
            }
            if (inx_2 < src2_num) {
                row_2 = src2_row[inx_2];
            }
        }

        // Case B: the element exists only in first column
        else if (row_1 < row_2) {
            (*num_1)++;
            inx_1++;
            if (inx_1 < src1_num) {
                row_1 = src1_row[inx_1];
            }
        }

        // Case C: the element exists only in second column
        else {
            (*num_2)++;
            inx_2++;
            if (inx_2 < src2_num) {
                row_2 = src2_row[inx_2];
            }
        }

    } // endwhile: counted elements

    // At this point at least one column expired of elements. If there are
    // still elements remaining in the first column then count them now.
    if (inx_1 < src1_num) {
        *num_1 += (src1_num - inx_1);
    }

    // If there are still elements remaining in the second column then count
    // them now.
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
 * @param[in] dst_data Data array [0...dst_num-1] of mixed column.
 * @param[in] dst_row Row index array [0...dst_num-1] of mixed column.
 * @param[out] dst_num Number of elements in mixed column.
 *
 * This method mixes two sparse matrix columns into a single column.
 ***************************************************************************/
void GMatrixSparse::mix_column(const double* src1_data,
                               const int*    src1_row,
                               int           src1_num,
                               const double* src2_data,
                               const int*    src2_row,
                               int           src2_num,
                               double*       dst_data,
                               int*          dst_row,
                               int*          dst_num)
{
    // Initialise indices and row indices of both columns
    int inx   = 0;                    // Mixed column element index
    int inx_1 = 0;                    // Column 1 element index
    int inx_2 = 0;                    // Column 2 element index
    int row_1 = src1_row[inx_1];      // Current row in column 1
    int row_2 = src2_row[inx_2];      // Current row in column 2

    // Mix elements of both columns while both contain still elements
    while (inx_1 < src1_num && inx_2 < src2_num) {

        // Case A: the element exists in both columns, so we add the values
        if (row_1 == row_2) {
            dst_data[inx] = src1_data[inx_1++] + src2_data[inx_2++];
            dst_row[inx]  = row_1;
            if (inx_1 < src1_num) {
                row_1 = src1_row[inx_1];
            }
            if (inx_2 < src2_num) {
                row_2 = src2_row[inx_2];
            }
        }

        // Case B: the element exists only in first column, so we copy the
        // element from the first column
        else if (row_1 < row_2) {
            dst_data[inx] = src1_data[inx_1++];
            dst_row[inx]  = row_1;
            if (inx_1 < src1_num) {
                row_1 = src1_row[inx_1];
            }
        }

        // Case C: the element exists only in second column, so we copy the
        // element from the second column
        else {
            dst_data[inx] = src2_data[inx_2++];
            dst_row[inx]  = row_2;
            if (inx_2 < src2_num) {
                row_2 = src2_row[inx_2];
            }
        }

        // Update the destination index since we added a element
        inx++;

    } // endwhile: mixing

    // At this point either the first or the second column expired of elements
    // In the case that there are still elements remaining in the first column
    // we add them now ...
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

    // Now store the number of elements in the mixed column
    *dst_num = inx;

    // We're done
    return;
}


/***********************************************************************//**
 * @brief Insert row in matrix
 *
 * @param[in] row Matrix row [0,...,row()[.
 * @param[in] vector Vector.
 * @param[in] add Add vector to existing elements?
 *
 * @exception GException::out_of_range
 *            Invalid matrix row specified.
 * @exception GException::invalid_argument
 *            Vector size incompatible with number of matrix columns.
 *
 * @todo Remove elements that are empty after addition
 ***************************************************************************/
void GMatrixSparse::insert_row(const int&     row,
                               const GVector& vector,
                               const bool&    add)
{
    // Raise an exception if the row index is invalid
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows) {
        throw GException::out_of_range(G_SET_ROW, "Matrix row", row, m_rows);
    }
    #endif

    // Raise an exception if the vector length is invalid
    if (vector.size() != m_cols) {
        std::string msg = "Vector size "+gammalib::str(vector.size())+
                          " does not match "+gammalib::str(m_cols)+" matrix "
                          "columns. Please specify a vector of size "+
                          gammalib::str(m_cols)+".";
        throw GException::invalid_argument(G_INSERT_ROW, msg);
    }

    // Fill pending element into matrix
    fill_pending();

    // Determine indices of new columns that need to be allocated and set
    // or add all vector elements that do not need a new memory allocation
    std::vector<int> new_columns;
    std::vector<int> k_insert;
    for (int i = m_cols-1; i >= 0; --i) {
        if (vector[i] != 0.0) {
            int  k          = m_colstart[i];
            int  istop      = m_colstart[i+1];
            bool add_column = true;
            for (; k < istop; ++k) {
                if (m_rowinx[k] == row) {
                    if (add) {
                        m_data[k] += vector[i];
                    }
                    else {
                        m_data[k] = vector[i];
                    }
                    add_column = false;
                    break;
                }
                else if (m_rowinx[k] > row) {
                    break;
                }
            }
            if (add_column) {
                new_columns.push_back(i);
                k_insert.push_back(k);
            }
        }
    }

    // If new columns are needed then allocate new memory and insert
    // new columns
    if (!new_columns.empty()) {

        // Determine the requested new logical size of the matrix
        int new_size = m_elements + new_columns.size();

        // Initialise pointers to destination
        bool    new_memory = false;
        double* new_data   = m_data;
        int*    new_rowinx = m_rowinx;

        // If requested size is smaller than allocated size then allocate
        // new memory
        if (new_size > m_alloc) {

            // Propose a new memory size
            int new_propose = m_alloc + m_mem_block;

            // Make sure that enough memory is allocated
            m_alloc = (new_size > new_propose) ? new_size : new_propose;

            // Allocate memory for new elements
            new_data   = new double[m_alloc];
            new_rowinx = new int[m_alloc];

            // Signal that new memory was allocated
            new_memory = true;

        } // endif: memory allocation required

        // Loop over all new columns
        int n_insert = new_columns.size();
        int k_last   = m_elements;
        for (int i = 0; i < new_columns.size(); ++i) {

            // Move data back by n_insert positions
            for (int k = k_insert[i]; k < k_last; ++k) {
                new_data[k+n_insert]   = m_data[k];
                new_rowinx[k+n_insert] = m_rowinx[k];
            }

            // Insert vector element just before the block that was moved
            // back
            new_data[k_insert[i]+n_insert-1]   = vector[new_columns[i]];
            new_rowinx[k_insert[i]+n_insert-1] = row;

            // Update number of positions to move back and end of block
            n_insert--;
            k_last = k_insert[i];

        } // endfor: looped over all new columns

        // Update column start and number of elements
        int i_insert = 0;
        n_insert = new_columns.size();
        for (int i = m_cols; i > 0; --i) {
            if (new_columns[i_insert] == i) {
                n_insert--;
                i_insert++;
                if (n_insert < 1) {
                    break;
                }
            }
            m_colstart[i] += n_insert;
        }

        // Update number of elements
        m_elements = new_size;

        // If memory was allocated then free old memory and attach new
        // memory
        if (new_memory) {

            // Delete old memory
            if (m_data   != NULL) delete [] m_data;
            if (m_rowinx != NULL) delete [] m_rowinx;

            // Update pointers to new memory and update element counter
            m_data      = new_data;
            m_rowinx    = new_rowinx;

        } // endif: memory was allocated

    } // endif: new columns were needed

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           Friend functions                              =
 =                                                                         =
 ==========================================================================*/

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
 * @param[in] values Flag that signals if values should be copied.
 *
 * The flag 'values' allows to avoid copying the actual data values. This
 * allows to perform a logical matrix transposition, as needed by the
 * symbolic matrix analysis class.
 *
 * Note that this method does not support pending elements (they have to be
 * filled before).
 *
 * The method does nothing on empty matrices.
 ***************************************************************************/
GMatrixSparse cs_transpose(const GMatrixSparse& matrix, int values)
{
    // Declare and allocate result matrix 
    GMatrixSparse result(matrix.m_cols, matrix.m_rows, matrix.m_elements);

    // Transpose only if matrix is not empty
    if (!matrix.is_empty()) {

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

    } // endif: matrix was not empty

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
