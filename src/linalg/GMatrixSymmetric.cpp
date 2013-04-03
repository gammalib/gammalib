/***************************************************************************
 *                GMatrixSymmetric.cpp - Symmetric matrix class            *
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
 * @file GMatrixSymmetric.cpp
 * @brief Symmetric matrix class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cmath>
#include "GTools.hpp"
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrix.hpp"
#include "GMatrixSparse.hpp"
#include "GMatrixSymmetric.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_CAST_MATRIX          "GMatrixSymmetric::GMatrixSymmetric(GMatrix&)"
#define G_CAST_SPARSEMATRIX "GMatrixSymmetric::GMatrixSymmetric(GSparseMatrix&)"
#define G_ACCESS1                    "GMatrixSymmetric::operator(int&, int&)"
#define G_ACCESS2                     "GMatrixSymmetric::operator(int&,int&)"
#define G_OP_ADD            "GMatrixSymmetric::operator+=(GMatrixSymmetric&)"
#define G_OP_SUB            "GMatrixSymmetric::operator-=(GMatrixSymmetric&)"
#define G_OP_MUL_VEC                  "GMatrixSymmetric::operator*(GVector&)"
#define G_OP_MAT_MUL        "GMatrixSymmetric::operator*=(GMatrixSymmetric&)"
#define G_INVERT                             "GMatrixSymmetric::invert(void)"
#define G_ADD_COL                 "GMatrixSymmetric::add_col(GVector&, int&)"
#define G_INSERT_COL           "GMatrixSymmetric::insert_col(GVector&, int&)"
#define G_CHOL_DECOMP            "GMatrixSymmetric::cholesky_decompose(int&)"
#define G_CHOL_SOLVE      "GMatrixSymmetric::cholesky_solver(GVector&, int&)"
#define G_CHOL_INVERT               "GMatrixSymmetric::cholesky_invert(int&)"
#define G_EXTRACT_ROW                   "GMatrixSymmetric::extract_row(int&)"
#define G_EXTRACT_COL                   "GMatrixSymmetric::extract_col(int&)"
#define G_COPY_MEMBERS    "GMatrixSymmetric::copy_members(GMatrixSymmetric&)"
#define G_ALLOC_MEMBERS         "GMatrixSymmetric::alloc_members(int&, int&)"


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Initialises an empty matrix with no rows and columns
 *
 * @todo Check the entire matrix class to see whether handling with an
 *       empty matrix may lead to some conflicts.
 ***************************************************************************/
GMatrixSymmetric::GMatrixSymmetric(void) : GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}

/***********************************************************************//**
 * @brief Matrix constructor
 *
 * @param[in] rows Number of rows in matrix.
 * @param[in] cols Number of columns in matrix.
 ***************************************************************************/
GMatrixSymmetric::GMatrixSymmetric(const int& rows, const int& cols) :
                  GMatrixBase()
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(rows, cols);

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrix to GMatrixSymmetric storage class convertor
 *
 * @param[in] matrix General matrix (GMatrix).
 *
 * @exception GException::matrix_not_symmetric
 *            Matrix is not symmetric.
 *
 * Converts a general matrix into the symmetric storage class. If the input
 * matrix is not symmetric, an exception is thrown.
 ***************************************************************************/
GMatrixSymmetric::GMatrixSymmetric(const GMatrix& matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(matrix.rows(), matrix.cols());

    // Fill matrix
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = col; row < matrix.rows(); ++row) {
            double value_ll = matrix(row,col);
            double value_ur = matrix(col,row);
            if (value_ll != value_ur) {
                throw GException::matrix_not_symmetric(G_CAST_MATRIX,
                                                       matrix.rows(),
                                                       matrix.cols());
            }
            (*this)(row, col) = matrix(row, col);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief GMatrixSparse to GMatrixSymmetric storage class convertor
 *
 * @param[in] matrix Sparse matrix (GMatrixSparse).
 *
 * @exception GException::matrix_not_symmetric
 *            Sparse matrix is not symmetric.
 *
 * Converts a sparse matrix into the symmetric storage class. If the input
 * matrix is not symmetric, an exception is thrown.
 ***************************************************************************/
GMatrixSymmetric::GMatrixSymmetric(const GMatrixSparse& matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Allocate matrix memory
    alloc_members(matrix.rows(), matrix.cols());

    // Fill matrix
    for (int col = 0; col < matrix.cols(); ++col) {
        for (int row = col; row < matrix.rows(); ++row) {
            double value_ll = matrix(row,col);
            double value_ur = matrix(col,row);
            if (value_ll != value_ur) {
                throw GException::matrix_not_symmetric(G_CAST_SPARSEMATRIX,
                                                       matrix.rows(),
                                                       matrix.cols());
            }
            (*this)(row, col) = matrix(row, col);
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
GMatrixSymmetric::GMatrixSymmetric(const GMatrixSymmetric& matrix) :
                  GMatrixBase(matrix)
{
    // Initialise private members for clean destruction
    init_members();

    // Copy members
    copy_members(matrix);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GMatrixSymmetric::~GMatrixSymmetric(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                               Operators                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
GMatrixSymmetric& GMatrixSymmetric::operator=(const GMatrixSymmetric& matrix)
{
    // Execute only if object is not identical
    if (this != &matrix) {

        // Assign base class members. Note that this method will also perform
        // the allocation of the matrix memory and the copy of the matrix
        // attributes.
        this->GMatrixBase::operator=(matrix);

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
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
double& GMatrixSymmetric::operator()(const int& row, const int& col)
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_ACCESS1, row, col, m_rows, m_cols);
    }
    #endif

    // Get element index
    int inx = (row >= col) ? m_colstart[col]+(row-col) : m_colstart[row]+(col-row);

    // Return element
    return m_data[inx];
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
const double& GMatrixSymmetric::operator()(const int& row, const int& col) const
{
    // Compile option: perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
        throw GException::out_of_range(G_ACCESS2, row, col, m_rows, m_cols);
    #endif

    // Get element index
    int inx = (row >= col) ? m_colstart[col]+(row-col) : m_colstart[row]+(col-row);

    // Return element
    return m_data[inx];
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
GVector GMatrixSymmetric::operator*(const GVector& vector) const
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
GMatrixSymmetric& GMatrixSymmetric::operator+=(const GMatrixSymmetric& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows || m_cols != matrix.m_cols) {
        throw GException::matrix_mismatch(G_OP_ADD,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Add matrices
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
GMatrixSymmetric& GMatrixSymmetric::operator-= (const GMatrixSymmetric& matrix)
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_rows != matrix.m_rows || m_cols != matrix.m_cols) {
        throw GException::matrix_mismatch(G_OP_SUB,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

    // Subtract matrices
    subtraction(matrix);

    // Return result
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                              Public methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 ***************************************************************************/
void GMatrixSymmetric::clear(void)
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
GMatrixSymmetric* GMatrixSymmetric::clone(void) const
{
    // Clone this image
    return new GMatrixSymmetric(*this);
}


/***********************************************************************//**
 * @brief Invert matrix
 *
 * @exception GException::feature_not_implemented
 *            Feature not yet implemented.
 *
 * @todo Needs to be implemented.
 ***************************************************************************/
void GMatrixSymmetric::invert(void)
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
void GMatrixSymmetric::add_col(const GVector& vector, const int& col)
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

    // Insert column into vector
    for (int row = 0; row < m_rows; ++row) {
        (*this)(row,col) += vector[row];
    }

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
void GMatrixSymmetric::insert_col(const GVector& vector, const int& col)
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

    // Insert column into vector
    for (int row = 0; row < m_rows; ++row) {
        (*this)(row,col) = vector[row];
    }

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
GVector GMatrixSymmetric::extract_row(const int& row) const
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
    for (int col = 0; col < m_cols; ++col) {
        result[col] = (*this)(row,col);
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
GVector GMatrixSymmetric::extract_col(const int& col) const
{
    // Raise an exception if the column index is invalid
    #if defined(G_RANGE_CHECK)
    if (col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_EXTRACT_COL, col, 0, m_cols-1);
    }
    #endif

    // Create result vector
    GVector result(m_rows);

    // Extract column into vector
    for (int row = 0; row < m_rows; ++row) {
        result[row] = (*this)(row,col);
    }

    // Return vector
    return result;
}


/***********************************************************************//**
 * @brief Extract lower triangle of matrix
 *
 * This method extracts the lower triangle of a matrix into another matrix.
 * (including the diagonal elements).
 * All remaining matrix elements will be zero.
 ***************************************************************************/
GMatrix GMatrixSymmetric::extract_lower_triangle(void) const
{
    // Define result matrix
    GMatrix result(m_rows, m_cols);

    // Extract all elements
    for (int row = 0; row < m_rows; ++row) {
        for (int col = 0; col <= row; ++col) {
            result(row,col) = m_data[m_colstart[col]+(row-col)];
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Extract upper triangle of matrix
 *
 * This method extracts the upper triangle of a matrix into another matrix.
 * (including the diagonal elements).
 * All remaining matrix elements will be zero.
 ***************************************************************************/
GMatrix GMatrixSymmetric::extract_upper_triangle(void) const
{
    // Define result matrix
    GMatrix result(m_rows, m_cols);

    // Extract all elements
    for (int row = 0; row < m_rows; ++row) {
        for (int col = row; col < m_cols; ++col) {
            result(row,col) = m_data[m_colstart[row]+(col-row)];
        }
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Perform a Cholesky decomposition
 *
 * @param[in] compress Use zero-row/column compression (default: true).
 *
 * @exception GException::matrix_not_pos_definite
 *            Matrix is not positive definite.
 * @exception GException::matrix_zero
 *            All matrix elements are zero.
 *
 * Inplace Cholesky decomposition inspired by Numerical Recipes algorithm.
 * The decomposition, which is a matrix occupying only the lower triange,
 * is stored in the elements of the symmetric matrix. To visualise the
 * matrix one has to use 'lower_triangle()' to extract the relevant part.
 * Case A operates on a full matrix, Case B operates on a (logically)
 * compressed matrix where zero rows/columns have been removed.
 ***************************************************************************/
void GMatrixSymmetric::cholesky_decompose(bool compress)
{
    // Set-up incides of non zero rows if matrix compression is requested
    if (compress) {
        set_inx();
    }

    // Check if zero-row/col compression is needed  
    int no_zeros = ((compress && (m_num_inx == m_rows)) || !compress);

    // Case A: no zero-row/col compression needed
    if (no_zeros) {

        // Loop over upper triangle (col >= row)
        double diag = 0.0;
        for (int row = 0; row < m_rows; ++row) {
            double* ptr = m_data + m_colstart[row];
            for (int col = row; col < m_cols; ++col, ++ptr) {
                double sum = *ptr;                         // sum = M(row,col)
                for (int k = 0; k < row; ++k) {
                    int offset = m_colstart[k] - k;        // is always positive
                    sum -= m_data[offset+row] * m_data[offset+col]; // sum -= M(row,k)*M(col,k)
                }
                if (row == col) {
                    if (sum <= 0.0) {
                        throw GException::matrix_not_pos_definite(G_CHOL_DECOMP, row, sum);
                    }
                    *ptr = sqrt(sum);                      // M(row,row) = sqrt(sum)
                    diag = 1.0/(*ptr);
                }
                else
                    *ptr = sum*diag;                       // M(row,col) = sum/M(row,row)
            }
        }
    } // endif: there were no zero rows/cols in matrix

    // Case B: zero-row/col compression needed
    else if (m_num_inx > 0) {

        // Allocate loop variables and pointers
        int  row;
        int  col;
        int  k;
        int* row_ptr;
        int* col_ptr;
        int* k_ptr;

        // Loop over upper triangle (col >= row)
        double diag = 0.0;
        for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
            double* ptr_0 = m_data + m_colstart[*row_ptr] - *row_ptr;
            for (col = row, col_ptr = m_inx + row; col < m_num_inx; ++col, ++col_ptr) {
                double* ptr = ptr_0 + *col_ptr;
                double  sum = *ptr;                                  // sum = M(row,col)
                for (k = 0, k_ptr = m_inx; k < row; ++k, ++k_ptr) {
                    int offset = m_colstart[*k_ptr] - *k_ptr;        // is always positive
                    sum -= m_data[offset+*row_ptr] * m_data[offset+*col_ptr];
                                                                     // sum -= M(row,k)*M(col,k)
                }
                if (*row_ptr == *col_ptr) {
                    if (sum <= 0.0) {
                        throw GException::matrix_not_pos_definite(G_CHOL_DECOMP, *row_ptr, sum);
                    }
                    *ptr = sqrt(sum);                                // M(row,row) = sqrt(sum)
                    diag = 1.0/(*ptr);
                }
                else
                    *ptr = sum*diag;                                 // M(row,col) = sum/M(row,row)
            }
        }
    } // endelse: zero-row/col compression needed

    // Case C: all matrix elements are zero
    else {
        throw GException::matrix_zero(G_CHOL_DECOMP);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Cholesky solver
 *
 * @param[in] vector Vector for which should be solved.
 * @param[in] compress Use zero-row/column compression (default: true).
 *
 * @exception GException::matrix_vector_mismatch
 *            Matrix and vector do not match.
 * @exception GException::matrix_zero
 *            All matrix elements are zero.
 *
 * Solves the linear equation A*x=b using a Cholesky decomposition of A.
 * This function is to be applied on a decomposition GMatrixSymmetric matrix
 * that is produced by 'cholesky_decompose'. Case A operates on a full
 * matrix, Case B on a zero rows/columns (logically) compressed matrix.
 ***************************************************************************/
GVector GMatrixSymmetric::cholesky_solver(const GVector& vector, bool compress)
{
    // Raise an exception if the matrix and vector dimensions are not compatible
    if (m_rows != vector.size()) {
        throw GException::matrix_vector_mismatch(G_CHOL_SOLVE, vector.size(),
                                                 m_rows, m_cols);
    }

    // Allocate result vector
    GVector x(m_rows);

    // Check if zero-row/col compression is needed  
    int no_zeros = ((compress && (m_num_inx == m_rows)) || !compress);

    // Case A: no zero-row/col compression needed
    if (no_zeros) {

        // Solve L*y=b, storing y in x (row>k)
        for (int row = 0; row < m_rows; ++row) {
            double sum = vector[row];
            for (int k = 0; k < row; ++k) {
                sum -= m_data[m_colstart[k]+(row-k)] * x[k]; // sum -= M(row,k) * x(k)
            }
            x[row] = sum/m_data[m_colstart[row]];            // x(row) = sum/M(row,row)
        }

        // Solve trans(L)*x=y (k>row)
        for (int row = m_rows-1; row >= 0; --row) {
            double  sum = x[row];
            double* ptr = m_data + m_colstart[row] + 1;
            for (int k = row+1; k < m_rows; ++k) {
                sum -= *ptr++ * x[k];               // sum -= M(k,row) * x(k)
            }
            x[row] = sum/m_data[m_colstart[row]];   // x(row) = sum/M(row,row)
        }
    } // endif: no zero-row/col compression needed

    // Case B: zero-row/col compression needed
    else if (m_num_inx > 0) {

        // Allocate loop variables and pointers
        int  row;
        int  k;
        int* row_ptr;
        int* k_ptr;

        // Solve L*y=b, storing y in x (row>k)
        for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
            double  sum = vector[*row_ptr];
            double* ptr = m_data + *row_ptr;
            for (k = 0, k_ptr = m_inx; k < row; ++k, ++k_ptr) {
                sum -= *(ptr + m_colstart[*k_ptr] - *k_ptr) * x[*k_ptr]; // sum -= M(row,k) * x(k)
            }
            x[*row_ptr] = sum/m_data[m_colstart[*row_ptr]];              // x(row) = sum/M(row,row)
        }

        // Solve trans(L)*x=y (k>row)
        for (row = m_num_inx-1, row_ptr = m_inx+m_num_inx-1; row >= 0; --row, --row_ptr) {
            double  sum      = x[*row_ptr];
            double* ptr_diag = m_data + m_colstart[*row_ptr];
            double* ptr      = ptr_diag - *row_ptr;
            for (k = row+1, k_ptr = m_inx+row+1; k < m_num_inx; ++k, ++k_ptr) {
                sum -= *(ptr + *k_ptr) * x[*k_ptr];              // sum -= M(k,row) * x(k)
            }
            x[*row_ptr] = sum/(*ptr_diag);                       // x(row) = sum/M(row,row)
        }
    } // endelse: zero-row/col compression needed

    // Case C: all matrix elements are zero
    else {
        throw GException::matrix_zero(G_CHOL_SOLVE);
    }

  // Return result vector
  return x;
}


/***********************************************************************//**
 * @brief Cholesky invert
 *
 * @param[in] compress Use zero-row/column compression (default: true).
 *
 * @exception GException::matrix_zero
 *            All matrix elements are zero.
 *
 * Inplace matrix inversion using the Cholesky decomposition. Case A
 * operates on a full matrix while Case B operates on a (logically)
 * compressed matrix where all zero rows/columns are skipped.
 ***************************************************************************/
void GMatrixSymmetric::cholesky_invert(bool compress)
{
    // Generate Cholesky decomposition of matrix
    this->cholesky_decompose(compress);

    // Check if zero-row/col compression is needed
    int no_zeros = ((compress && (m_num_inx == m_rows)) || !compress);

    // Case A: no zero-row/col compression needed
    if (no_zeros) {

        // Generate inverse of Cholesky decomposition (col>row)
        for (int row = 0; row < m_rows; ++row) {
            double* ptr = m_data + m_colstart[row];
            *ptr        = 1.0/(*ptr);                       // M(row,row) = 1/M(row,row)
            for (int col = row+1; col < m_cols; ++col) {
                double   sum = 0.0;
                double* ptr1 = m_data + col - row;
                double* ptr2 = ptr;
                for (int k = row; k < col; ++k) {
                    sum -= *(ptr1-- + m_colstart[k]) * *ptr2++; // sum -= M(col,k)*M(k,row)
                }
                *(ptr+col-row) = sum/m_data[m_colstart[col]];   // M(col,row) = sum/M(col,col)
            }
        }

        // Matrix multiplication (col>=row)
        for (int row = 0; row < m_rows; ++row) {
            double* ptr = m_data + m_colstart[row];
            for (int col = row; col < m_cols; ++col) {
                double   sum = 0.0;
                double* ptr1 = ptr + col - row;
                double* ptr2 = m_data + m_colstart[col];
                for (int k = col; k < m_cols; ++k) {
                    sum += *ptr1++ * *ptr2++;                 // sum += M(row,k)*M(k,col)
                }
                *(ptr+col-row) = sum;                         // M(row,col) = sum
            }
        }
    } // endif: no zero-row/col compression needed

    // Case B: zero-row/col compression needed
    else if (m_num_inx > 0) {

        // Allocate loop variables and pointers
        int  row;
        int  col;
        int  k;
        int* row_ptr;
        int* col_ptr;
        int* k_ptr;

        // Generate inverse of Cholesky decomposition (col>row)
        for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
            double* ptr_diag = m_data + m_colstart[*row_ptr];
            double* ptr_2    = ptr_diag - *row_ptr;
            *ptr_diag = 1.0/(*ptr_diag);                         // M(row,row) = 1/M(row,row)
            for (col = row+1, col_ptr = m_inx+row+1; col < m_num_inx; ++col, ++col_ptr) {
                double  sum   = 0.0;
                double* ptr_1 = m_data + *col_ptr;
                for (k = row, k_ptr = m_inx+row; k < col; ++k, ++k_ptr) {
                    sum -= *(ptr_1 + m_colstart[*k_ptr] - *k_ptr) *
                           *(ptr_2 + *k_ptr);                    // sum -= M(col,k)*M(k,row)
                }
                *(ptr_2 + *col_ptr) = sum/m_data[m_colstart[*col_ptr]];
                                                                 // M(col,row) = sum/M(col,col)
            }
        }

        // Matrix multiplication (col>=row)
        for (row = 0, row_ptr = m_inx; row < m_num_inx; ++row, ++row_ptr) {
            double* ptr_diag = m_data + m_colstart[*row_ptr];
            double* ptr_1    = ptr_diag - *row_ptr;
            for (col = row, col_ptr = m_inx+row; col < m_num_inx; ++col, ++col_ptr) {
                double  sum   = 0.0;
                double* ptr_2 = m_data + m_colstart[*col_ptr] - *col_ptr;
                for (k = col, k_ptr = m_inx+col; k < m_num_inx; ++k, ++k_ptr) {
                    sum += *(ptr_1 + *k_ptr) * *(ptr_2 + *k_ptr); // sum += M(row,k)*M(k,col)
                }
                *(ptr_1 + *col_ptr) = sum;                        // M(row,col) = sum
            }
        }
    } // endelse: zero-row/col compression needed

    // Case C: all matrix elements are zero
    else {
        throw GException::matrix_zero(G_CHOL_INVERT);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Determine fill of matrix
 *
 * The fill of a matrix is defined as the number of non-zero elements
 * devided by the total number of matrix elements.
 ***************************************************************************/
double GMatrixSymmetric::fill(void) const
{
    // Determine the number of zero elements
    int zero = 0;
    for (int col = 0, i = 0; col < m_cols; ++col) {
        if (m_data[i++] == 0.0) {                    // Count diag. once
            zero++;
        }
        for (int row = col+1; row < m_rows; ++row) {
            if (m_data[i++] == 0.0) {                // Count off-diag. twice
                zero +=2;
            }
        }
    }

    // Return the fill
    return (1.0-double(zero)/double(m_elements));
}


/***********************************************************************//**
 * @brief Sum matrix elements
 ***************************************************************************/
double GMatrixSymmetric::sum(void) const
{
    // Initialise matrix sums (diagonal and off-diagonal)
    double diag     = 0.0;
    double off_diag = 0.0;

    // Calulate sum over diagonal elements
    for (int row = 0; row < m_rows; ++row) {
        diag += m_data[m_colstart[row]];
    }

    // Calulate sum over off-diagonal elements
    for (int row = 0; row < m_rows; ++row) {
        for (int col = row+1; col < m_cols; ++col) {
            off_diag += m_data[m_colstart[row]+(col-row)];
        }
    }

    // Calculate total
    double result = diag + 2.0 * off_diag;

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print matrix
 ***************************************************************************/
std::string GMatrixSymmetric::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GMatrixSymmetric ===");
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
 * @brief Initialise class mambers
 ***************************************************************************/
void GMatrixSymmetric::init_members(void)
{
    // Initialise members
    m_num_inx = 0;
    m_inx     = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
void GMatrixSymmetric::copy_members(const GMatrixSymmetric& matrix)
{
    // Copy attributes
    m_num_inx = matrix.m_num_inx;

    // Copy index selection
    if (m_cols > 0) {
        m_inx     = new int[m_cols];
        for (int i = 0; i < m_cols; ++i) {
            m_inx[i] = matrix.m_inx[i];
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMatrixSymmetric::free_members(void)
{
    // De-allocate only if memory has indeed been allocated by derived class
    if (m_inx != NULL) delete [] m_inx;

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
 *            Attempt to allocate zero size matrix.
 * @exception GException::matrix_not_symmetric
 *            Matrix is not symmetric.
 *
 * This method is the main constructor code that allocates and initialises
 * memory for matrix elements. The method assumes that no memory has been
 * allocated for the matrix elements, the column start index array and the
 * index array.
 * The method allocates the memory for matrix elements, the column start
 * indices and the index array, sets all matrix elements to 0.0, and sets
 * the column start indices. The content of the index array is undefined.
 *
 * @todo Verify if the index array m_inx should be initialized.
 ***************************************************************************/
void GMatrixSymmetric::alloc_members(const int& rows, const int& cols)
{
    // Determine number of physical elements in matrix
    int elements = rows*(rows+1)/2;

    // Throw exception if requested matrix size is zero
    if (elements == 0) {
        throw GException::empty(G_ALLOC_MEMBERS);
    }

    // Throw exception if number of rows and columns is not identical
    if (rows != cols) {
      throw GException::matrix_not_symmetric(G_ALLOC_MEMBERS, rows, cols);
    }

    // Allocate matrix array and column start index array.
    m_data     = new double[elements];
    m_colstart = new int[cols+1];
    m_inx      = new int[cols];

    // Store matrix size (logical and physical)
    m_rows     = rows;
    m_cols     = cols;
    m_elements = elements;
    m_alloc    = elements;

    // Set-up column start indices
    m_colstart[0]   = 0;
    int offset = rows;
    for (int col = 1; col <= m_cols; ++col) {
        m_colstart[col] = m_colstart[col-1] + offset--;
    }

    // Initialise matrix elements to 0.0
    for (int i = 0; i < m_elements; ++i) {
        m_data[i] = 0.0;
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set index selection
 *
 * Determines the non-zero rows/columns in matrix and set up index array
 * that points to these rows/columns. This array is used for compressed
 * matrix calculations (Case B in the above routines).
 ***************************************************************************/
void GMatrixSymmetric::set_inx(void)
{
    // Allocate loop variables and pointers
    int row;
    int col;

    // Find all nonzero rows/cols
    m_num_inx = 0;
    for (row = 0; row < m_rows; ++row) {

        // Check first for zero diagonal. If we have one then check the rest of
        // the row
        if (m_data[m_colstart[row]] == 0.0) {
            for (col = 0; col < row; ++col) {
                if (m_data[m_colstart[col]+(row-col)] != 0.0) {
                    break;
                }
            }
            // Found a non-zero element
            if (col < row) {
                m_inx[m_num_inx++] = row;
            }
            else {
                for (col = row+1; col < m_cols; ++col) {
                    if (m_data[m_colstart[row]+(col-row)] != 0.0) {
                        break;
                    }
                }
                // Found a non-zero element
                if (col < m_cols) {
                    m_inx[m_num_inx++] = row;
                }
            }
        }
        else {
            m_inx[m_num_inx++] = row;
        }
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
 * @brief Return matrix of absolute values
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
GMatrixSymmetric abs(const GMatrixSymmetric& matrix)
{
    // Define result matrix
    GMatrixSymmetric result = matrix;

    // Convert all elements to absolute values
    for (int i = 0; i < result.m_elements; ++i) {
        result.m_data[i] = std::abs(result.m_data[i]);
    }

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Binary matrix multiplication operator
 *
 * @param[in] matrix Matrix to be multiplied.
 *
 * @exception GException::matrix_mismatch
 *            Incompatible matrix size.
 *
 * This method performs a matrix multiplication. Since the product of two
 * symmetric matrices is not necessarily symmetric, this method returns a
 * GMatrix object.
 *
 * The operation can only succeed when the dimensions of both matrices are
 * compatible.
 ***************************************************************************/
GMatrix GMatrixSymmetric::operator*(const GMatrixSymmetric& matrix) const
{
    // Raise an exception if the matrix dimensions are not compatible
    if (m_cols != matrix.m_rows) {
        throw GException::matrix_mismatch(G_OP_MAT_MUL,
                                          m_rows, m_cols,
                                          matrix.m_rows, matrix.m_cols);
    }

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
    
    // Return result
    return result;
}
