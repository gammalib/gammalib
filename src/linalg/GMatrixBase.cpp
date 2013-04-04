/***************************************************************************
 *                GMatrixBase.cpp - Abstract matrix base class             *
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
 * @file GMatrixBase.cpp
 * @brief Abstract matrix base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "GTools.hpp"
#include "GException.hpp"
#include "GVector.hpp"
#include "GMatrixBase.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_COPY_MEMBERS               "GMatrixBase::copy_members(GMatrixBase)"
#define G_SELECT_NON_ZERO                    "GMatrixBase::select_non_zero()"

/* __ Debugging definitions ______________________________________________ */
//#define G_DEBUG_COPY_MEMBERS                      // Dump copy_members info


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GMatrixBase::GMatrixBase(void)
{
    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] matrix Matrix.
 *
 * Constructs matrix by copying information from another matrix. The
 * constructor is sufficiently generic to provide the base constructor for
 * all derived classes, including sparse matrices.
 ***************************************************************************/
GMatrixBase::GMatrixBase(const GMatrixBase& matrix)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(matrix);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GMatrixBase::~GMatrixBase(void)
{
    // Free class members
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
 * @return Matrix.
 ***************************************************************************/
GMatrixBase& GMatrixBase::operator=(const GMatrixBase& matrix)
{
    // Execute only if object is not identical
    if (this != &matrix) {

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
 * @brief Equalty operator
 *
 * @param[in] matrix Matrix.
 * @return True if both matrices are identical.
 *
 * Checks if two matrices are identical. Two matrices are considered
 * identical if they have the same dimensions and identical elements.
 ***************************************************************************/
bool GMatrixBase::operator==(const GMatrixBase& matrix) const
{
    // Initalise result to true (are identical)
    bool result = true;

    // Test for difference. The loop will stop once then first difference
    // is encountered.
    if (m_rows     == matrix.m_rows &&
        m_cols     == matrix.m_cols &&
        m_elements == matrix.m_elements) {
        for (int i = 0; i < m_elements; ++i) {
            if (m_data[i] != matrix.m_data[i]) {
                result = false;
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
 * @return True if both matrices are different.
 *
 * Checks if two matrices are different. Two matrices are considered
 * different if they either have different dimensions or at least one element
 * that differs.
 ***************************************************************************/
bool GMatrixBase::operator!=(const GMatrixBase& matrix) const
{
    // Get negated result of equality operation
    bool result = !(this->operator==(matrix));
	
    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                            Protected methods                            =
 =                                                                         =
 ==========================================================================*/
 
/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GMatrixBase::init_members(void)
{
    // Initialise members
    m_rows       = 0;      // Number of rows
    m_cols       = 0;      // Number of columns
    m_elements   = 0;      // Logically used number of elements
    m_alloc      = 0;      // Allocated # of elements (>= m_elements)
    m_num_rowsel = 0;      // Number of selected rows (for comp. decomp.)
    m_num_colsel = 0;      // Number of selected columns (for comp. decomp.)
    m_data       = NULL;   // Matrix data
    m_colstart   = NULL;   // Column start indices (m_cols+1)
    m_rowsel     = NULL;   // Row selection (for compressed decomposition)
    m_colsel     = NULL;   // Column selection (for compressed decomposition)

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] matrix Matrix.
 ***************************************************************************/
void GMatrixBase::copy_members(const GMatrixBase& matrix)
{
    // Copy matrix attributes
    m_rows       = matrix.m_rows;
    m_cols       = matrix.m_cols;
    m_elements   = matrix.m_elements;
    m_alloc      = matrix.m_alloc;
    m_num_rowsel = matrix.m_num_rowsel;
    m_num_colsel = matrix.m_num_colsel;

    // Allocate only memory if we have rows and columns and data to copy
    if (m_rows > 0 && m_cols > 0) {

        // Allocate memory for column start array and copy content
        if (matrix.m_colstart != NULL) {
            m_colstart = new int[m_cols+1];
            for (int i = 0; i <= m_cols; ++i) {
                m_colstart[i] = matrix.m_colstart[i];
            }
        }

        // Allocate memory for elements and copy them
        if (matrix.m_data != NULL && m_alloc > 0) {
            m_data = new double[m_alloc];
            for (int i = 0; i < m_elements; ++i) {
                m_data[i] = matrix.m_data[i];
            }
        }

        // If there is a row selection then copy it
        if (matrix.m_rowsel != NULL && m_num_rowsel > 0) {
            m_rowsel = new int[m_num_rowsel];
            for (int i = 0; i < m_num_rowsel; ++i) {
                m_rowsel[i] = matrix.m_rowsel[i];
            }
        }

        // If there is a column selection then copy it
        if (matrix.m_colsel != NULL && m_num_colsel > 0) {
            m_colsel = new int[m_num_colsel];
            for (int i = 0; i < m_num_colsel; ++i) {
                m_colsel[i] = matrix.m_colsel[i];
            }
        }

    } // endif: there were data to copy

    // Optionally show debug information
    #if defined(G_DEBUG_COPY_MEMBERS)
    std::cout << "GMatrixBase::copy_members:"
              << " m_colstart=" << matrix.m_colstart << "->" << m_colstart
              << " m_data="     << matrix.m_data     << "->" << m_data
              << " m_rowsel="   << matrix.m_rowsel   << "->" << m_rowsel
              << " m_colsel="   << matrix.m_colsel   << "->" << m_colsel
              << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GMatrixBase::free_members(void)
{
    // De-allocate only if memory has indeed been allocated
    if (m_colsel   != NULL) delete [] m_colsel;
    if (m_rowsel   != NULL) delete [] m_rowsel;
    if (m_data     != NULL) delete [] m_data;
    if (m_colstart != NULL) delete [] m_colstart;

    // Properly mark members as free
    m_colstart = NULL;
    m_data     = NULL;
    m_rowsel   = NULL;
    m_colsel   = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Setup compressed matrix factorisation
 *
 * Determines the non-zero rows and columns in matrix and set up index
 * arrays that point to these rows/columns. These arrays are used for
 * compressed matrix factorisations.
 ***************************************************************************/
void GMatrixBase::select_non_zero(void)
{
    // Free existing selection arrays
    if (m_rowsel != NULL) delete [] m_rowsel;
    if (m_colsel != NULL) delete [] m_colsel;
  
    // Allocate selection arrays
    m_rowsel = new int[m_rows];
    m_colsel = new int[m_cols];

    // Initialise non-zero row and column counters
    m_num_rowsel = 0;
    m_num_colsel = 0;
  
    // Declare loop variables
    int row;
    int col;
  
    // Find all non-zero rows
    for (row = 0; row < m_rows; ++row) {
        for (col = 0; col < m_cols; ++col) {
            if ((*this)(row,col) != 0.0) {
                break;
            }
        }
        // Found a non-zero element in row
        if (col < m_cols) {
            m_rowsel[m_num_rowsel++] = row;
        }
    }

    // Find all non-zero columns
    for (col = 0; col < m_cols; ++col) {
        for (row = 0; row < m_rows; ++row) {
            if ((*this)(row,col) != 0.0)
                break;
        }
        // Found a non-zero element in column
        if (row < m_rows) {
            m_colsel[m_num_colsel++] = col;
        }
    }
  
    // Return
    return;
}



/***********************************************************************//**
 * @brief Scale all matrix elements with a scalar
 *
 * @param[in] scalar Scalar.
 *
 * Multiply all matrix elements with a scalar. There are three cases:
 * - the multiplier is 0, then reset all elements to 0,
 * - the multiplier is +/-1, then do nothing or negate,
 * - in any other case, multiply by multiplier.
 ***************************************************************************/
void GMatrixBase::scale_elements(const double& scalar)
{
    // Case 1: If multiplicator is 0 then set entire matrix to 0
    if (scalar == 0.0) {
        for (int i = 0; i < m_elements; ++i) {
            m_data[i] = 0.0;
        }
    }
  
    // Case 3: If multiplicator is not +/- 1 then perform multiplication
    else if (std::abs(scalar) != 1.0) {
        for (int i = 0; i < m_elements; ++i) {
            m_data[i] *= scalar;
        }
    }
  
    // Case 2: If multiplication is -1 then negate
    else if (scalar == -1.0) {
        for (int i = 0; i < m_elements; ++i) {
            m_data[i] = -m_data[i];
        }
    }
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set all elements to a specific value
 *
 * @param[in] value Value.
 *
 * Sets all matrix elements to a specific value.
 ***************************************************************************/
void GMatrixBase::set_all_elements(const double& value)
{
    // Set all matrix elements
    for (int i = 0; i < m_elements; ++i) {
        m_data[i] = value;
    }
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Return minimum matrix element
 *
 * @return Minimum element in matrix.
 *
 * Returns the minimum matrix element. If the matrix is empty, returns 0.
 ***************************************************************************/
double GMatrixBase::get_min_element(void) const
{
    // Initialise result
    double result = 0.0;

    // Continue only if there are elements
    if (m_elements > 0) {

        // Set actual minimum to first elements
        result = m_data[0];
  
        // Search all elements for the smallest one
        for (int i = 1; i < m_elements; ++i) {
            if (m_data[i] < result) {
                result = m_data[i];
            }
        }
    
    } // endif: there were elements
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns maximum matrix element 
 *
 * @return Maximum element in matrix.
 *
 * Returns the maximum matrix element. If the matrix is empty, returns 0.
 ***************************************************************************/
double GMatrixBase::get_max_element(void) const
{
    // Initialise result
    double result = 0.0;

    // Continue only if there are elements
    if (m_elements > 0) {

        // Set actual maximum to first elements
        result = m_data[0];
  
        // Search all elements for the largest one
        for (int i = 1; i < m_elements; ++i) {
            if (m_data[i] > result) {
                result = m_data[i];
            }
        }
    
    } // endif: there were elements
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns sum over all matrix elements
 ***************************************************************************/
double GMatrixBase::get_element_sum(void) const
{
    // Initialise matrix sum
    double result = 0.0;
  
    // Add all elements  
    for (int i = 0; i < m_elements; ++i) {
        result += m_data[i];
    }
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print all matrix elements
 *
 * @param[in] num Maximum number of rows and columns to print (default: 10)
 *
 * Prints all matrix elements into a string. The parameter max_elements
 * allows to control the maximum number of matrix elements that should be
 * printed. If set to 0, all elements will be printed. Otherwise, the number
 * of rows and columns will be limited by ommitting the central values.
 ***************************************************************************/
std::string GMatrixBase::print_elements(const int& num) const
{
    // Initialise result string
    std::string result;

    // Set row and column limits
    int row_stop  = 0;
    int row_start = 0;
    int col_stop  = 0;
    int col_start = 0;
    if (num > 0) {
        if (m_rows > num) {
            row_stop  = num/2;
            row_start = m_rows - row_stop;
        }
        if (m_cols > num) {
            col_stop  = num/2;
            col_start = m_cols - col_stop;
        }
    }

    // Print matrix elements row by row using the access function
    for (int row = 0; row < row_stop; ++row) {
        result += "\n ";
        for (int col = 0; col < col_stop; ++col) {
            result += str((*this)(row,col));
            if (col != m_cols-1) {
                result += ", ";
            }
        }
        if (col_start > col_stop) {
            result += "... ";
        }
        for (int col = col_start; col < m_cols; ++col) {
            result += str((*this)(row,col));
            if (col != m_cols-1) {
                result += ", ";
            }
        }
    }
    if (row_start > row_stop) {
        result += "\n ";
        for (int col = 0; col < col_stop; ++col) {
            result += "... ";
        }
        if (col_start > col_stop) {
            result += "... ";
        }
        for (int col = col_start; col < m_cols; ++col) {
            result += "... ";
        }
    }
    for (int row = row_start; row < m_rows; ++row) {
        result += "\n ";
        for (int col = 0; col < col_stop; ++col) {
            result += str((*this)(row,col));
            if (col != m_cols-1) {
                result += ", ";
            }
        }
        if (col_start > col_stop) {
            result += "... ";
        }
        for (int col = col_start; col < m_cols; ++col) {
            result += str((*this)(row,col));
            if (col != m_cols-1) {
                result += ", ";
            }
        }
    }
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print row compression scheme if it exists
 ***************************************************************************/
std::string GMatrixBase::print_row_compression(void) const
{
    // Initialise result string
    std::string result;

    // If there is a row compression the print the scheme
    if (m_rowsel != NULL) {
        result.append("\n"+parformat("Row selection"));
        for (int row = 0; row < m_num_rowsel; ++row) {
            result.append(" ");
            result.append(str(m_rowsel[row]));
        }
    }
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Print column compression scheme if it exists
 ***************************************************************************/
std::string GMatrixBase::print_col_compression(void) const
{
    // Initialise result string
    std::string result;

    // If there is a row compression the print the scheme
    if (m_colsel != NULL) {
        result.append("\n"+parformat("Column selection"));
        for (int col = 0; col < m_num_colsel; ++col) {
            result.append(" ");
            result.append(str(m_colsel[col]));
        }
    }
  
    // Return result
    return result;
}
