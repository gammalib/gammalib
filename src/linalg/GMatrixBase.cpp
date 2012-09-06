/***************************************************************************
 *               GMatrixBase.cpp  -  Abstract matrix base class            *
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
 *
 * Allocates an empty matrix that contains no elements.
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
 * The copy constructor is sufficiently general to provide the base
 * constructor for all derived classes, including sparse matrices.
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
 *
 * This operator checks if two matrices are identical. Two matrices are
 * considered identical if they have the same dimensions and identicial
 * elements.
 ***************************************************************************/
bool GMatrixBase::operator==(const GMatrixBase& matrix) const
{
    // Initalise result to true (are identical)
    bool result = true;

    // Perform test for non-identity. At the first non-identify we can
    // stop.
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
 *
 * This operator checks if two matrices are not identical. Two matrices are
 * considered not identical if they differ in their dimensions or if at
 * least one element differs.
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
    // Initialise GMatrixBase members
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
 * @brief Negate matrix elements
 *
 * This method inverts the sign of all matrix elements.
 ***************************************************************************/
void GMatrixBase::negation(void)
{
    // Inverts the sign for all matrix elements
    for (int i = 0; i < m_elements; ++i) {
        m_data[i] = -m_data[i];
    }
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Add matrix elements
 *
 * @param[in] matrix Matrix.
 *
 * Add all matrix elements.
 ***************************************************************************/
void GMatrixBase::addition(const GMatrixBase& matrix)
{
    // Add all elements of matrix
    for (int i = 0; i < m_elements; ++i) {
        m_data[i] += matrix.m_data[i];
    }
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Subtract matrix elements
 *
 * @param[in] matrix Matrix.
 *
 * Subtract all matrix elements.
 ***************************************************************************/
void GMatrixBase::subtraction(const GMatrixBase& matrix)
{
    // Add all elements of matrix
    for (int i = 0; i < m_elements; ++i) {
        m_data[i] -= matrix.m_data[i];
    }
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Multiply matrix elements with scalar
 *
 * @param[in] scalar Scalar.
 *
 * Multiply all matrix elements with a scalar. There are three cases:
 * (1) the multiplier is 0, then reset all elements to 0,
 * (2) the multiplier is +/-1, then do nothing or negate,
 * (3) in any other case, multiply by multiplier.
 ***************************************************************************/
void GMatrixBase::multiplication(const double& scalar)
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
 * @brief Returns minimum matrix element 
 ***************************************************************************/
double GMatrixBase::get_min_element(void) const
{
    // Initialise minimum with first element
    double result = m_data[0];
  
    // Search all elements for the smallest one
    for (int i = 1; i < m_elements; ++i) {
        if (m_data[i] < result) {
            result = m_data[i];
        }
    }
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns maximum matrix element 
 ***************************************************************************/
double GMatrixBase::get_max_element(void) const
{
    // Initialise maximum with first element
    double result = m_data[0];
  
    // Search all elements for the largest one
    for (int i = 1; i < m_elements; ++i) {
        if (m_data[i] > result) {
            result = m_data[i];
        }
    }
  
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
 ***************************************************************************/
std::string GMatrixBase::print_elements(void) const
{
    // Initialise result string
    std::string result;

    // Print matrix elements row by row using the access function
    for (int row = 0; row < m_rows; ++row) {
        result += "\n ";
        for (int col = 0; col < m_cols; ++col) {
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
