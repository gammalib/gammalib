/***************************************************************************
 *               GMatrixBase.cpp  -  Abstract matrix base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2006-2011 by Jurgen Knodlseder                           *
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
 * @brief GMatrixBase class implementation.
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#include <cmath>
#include "GMatrixBase.hpp"
#include "GException.hpp"
#include "GVector.hpp"

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
 * @param[in] m Matrix from which class should be instantiated.
 *
 * The copy constructor is sufficiently general to provide the base
 * constructor for all derived classes, including sparse matrices.
 ***************************************************************************/
GMatrixBase::GMatrixBase(const GMatrixBase& m)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(m);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GMatrixBase::~GMatrixBase()
{
    // Free class members
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
 * @param[in] m Matrix.
 ***************************************************************************/
GMatrixBase& GMatrixBase::operator= (const GMatrixBase& m)
{
    // Execute only if object is not identical
    if (this != &m) {

        // Free members
        free_members();

        // Initialise members
        init_members();

        // Copy members
        copy_members(m);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Equalty operator
 *
 * @param[in] m Matrix.
 *
 * Two matrixes are considered equal if they have the same dimensions and
 * identicial elements.
 ***************************************************************************/
int GMatrixBase::operator== (const GMatrixBase& m) const
{
    // Initalise the result to 'equal matrices'
    int result = 1;
  
    // Perform comparison (only if matrix dimensions and number of physical
    // elements are identical; otherwise set result to 'false')
    if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
        for (int i = 0; i < m_elements; ++i) {
            if (m_data[i] != m.m_data[i]) {
                result = 0;
                break;
            }
        }
    }
    else
        result = 0;

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Non-equality operator
 *
 * @param[in] m Matrix.
 *
 * Two matrixes are considered as non-equal (or different) if the differ
 * in their dimensions or if at least one element differs.
 ***************************************************************************/
int GMatrixBase::operator!= (const GMatrixBase& m) const
{
    // Initalise the result to 'equal matrices'
    int result = 0;

    // Perform comparison (only if matrix dimensions and number of physical
    // elements are identical; otherwise set result to 'true')
    if (m_rows == m.m_rows && m_cols == m.m_cols && m_elements == m.m_elements) {
        for (int i = 0; i < m_elements; ++i) {
            if (m_data[i] != m.m_data[i]) {
                result = 1;
                break;
            }
        }
    }
    else
        result = 1;
	
    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Protected methods                           =
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
 * @param[in] m Matrix from which members should be copied.
 ***************************************************************************/
void GMatrixBase::copy_members(const GMatrixBase& m)
{
    // Copy matrix attributes
    m_rows       = m.m_rows;
    m_cols       = m.m_cols;
    m_elements   = m.m_elements;
    m_alloc      = m.m_alloc;
    m_num_rowsel = m.m_num_rowsel;
    m_num_colsel = m.m_num_colsel;

    // Allocate only memory if we have rows and columns and data to copy
    if (m_rows > 0 && m_cols > 0) {

        // Allocate memory for column start array and copy content
        if (m.m_colstart != NULL) {
            m_colstart = new int[m_cols+1];
            for (int i = 0; i <= m_cols; ++i)
                m_colstart[i] = m.m_colstart[i];
        }

        // Allocate memory for elements and copy them
        if (m.m_data != NULL && m_alloc > 0) {
            m_data = new double[m_alloc];
            for (int i = 0; i < m_elements; ++i)
                m_data[i] = m.m_data[i];
        }

        // If there is a row selection then copy it
        if (m.m_rowsel != NULL && m_num_rowsel > 0) {
            m_rowsel = new int[m_num_rowsel];
            for (int i = 0; i < m_num_rowsel; ++i)
                m_rowsel[i] = m.m_rowsel[i];
        }

        // If there is a column selection then copy it
        if (m.m_colsel != NULL && m_num_colsel > 0) {
            m_colsel = new int[m_num_colsel];
            for (int i = 0; i < m_num_colsel; ++i)
                m_colsel[i] = m.m_colsel[i];
        }

    } // endif: there were data to copy

    // Optionally show debug information
    #if defined(G_DEBUG_COPY_MEMBERS)
    std::cout << "GMatrixBase::copy_members:"
              << " m_colstart=" << m.m_colstart << "->" << m_colstart
              << " m_data="     << m.m_data     << "->" << m_data
              << " m_rowsel="   << m.m_rowsel   << "->" << m_rowsel
              << " m_colsel="   << m.m_colsel   << "->" << m_colsel
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
    for (int i = 0; i < m_elements; ++i)
        m_data[i] = -m_data[i];
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Add matrix elements
 *
 * @param[in] m Matrix.
 *
 * Add all matrix elements.
 ***************************************************************************/
void GMatrixBase::addition(const GMatrixBase& m)
{
    // Add all elements of matrix
    for (int i = 0; i < m_elements; ++i)
        m_data[i] += m.m_data[i];
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Subtract matrix elements
 *
 * @param[in] m Matrix.
 *
 * Subtract all matrix elements.
 ***************************************************************************/
void GMatrixBase::subtraction(const GMatrixBase& m)
{
    // Add all elements of matrix
    for (int i = 0; i < m_elements; ++i)
        m_data[i] -= m.m_data[i];
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Multiply matrix elements with scalar
 *
 * @param[in] s Scalar.
 *
 * Multiply all matrix elements with a scalar. There are three cases:
 * (1) the multiplier is 0, then reset all elements to 0,
 * (2) the multiplier is +/-1, then do nothing or negate,
 * (3) in any other case, multiply by multiplier.
 ***************************************************************************/
void GMatrixBase::multiplication(const double& s)
{
    // Case A: If multiplicator is 0 then set entire matrix to 0
    if (s == 0.0) {
        for (int i = 0; i < m_elements; ++i)
            m_data[i] = 0.0;
    }
  
    // Case C: If multiplicator is not +/- 1 then perform multiplication
    else if (std::abs(s) != 1.0) {
        for (int i = 0; i < m_elements; ++i)
            m_data[i] *= s;
    }
  
    // Case B: If multiplication is -1 then negate
    else if (s == -1.0) {
        for (int i = 0; i < m_elements; ++i)
            m_data[i] = -m_data[i];
    }
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Set all elements to a specific value
 *
 * @param[in] s Value.
 ***************************************************************************/
void GMatrixBase::set_all_elements(const double& s)
{
    // Set all matrix elements
    for (int i = 0; i < m_elements; ++i)
        m_data[i] = s;
	
    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns minimum matrix element 
 ***************************************************************************/
double GMatrixBase::get_min_element() const
{
    // Initialise minimum with first element
    double result = m_data[0];
  
    // Search all elements for the smallest one
    for (int i = 1; i < m_elements; ++i) {
        if (m_data[i] < result)
            result = m_data[i];
    }
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns maximum matrix element 
 ***************************************************************************/
double GMatrixBase::get_max_element() const
{
    // Initialise maximum with first element
    double result = m_data[0];
  
    // Search all elements for the largest one
    for (int i = 1; i < m_elements; ++i) {
        if (m_data[i] > result)
            result = m_data[i];
    }
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Returns sum over all matrix elements
 ***************************************************************************/
double GMatrixBase::get_element_sum() const
{
    // Initialise matrix sum
    double result = 0.0;
  
    // Add all elements  
    for (int i = 0; i < m_elements; ++i)
        result += m_data[i];
  
    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Dump all matrix elements in ostream
 *
 * @param[in] os Output stream.
 ***************************************************************************/
void GMatrixBase::dump_elements(std::ostream& os) const
{
    // Dump matrix elements row by row using the access function
    for (int row = 0; row < m_rows; ++row) {
        os << " ";
        for (int col = 0; col < m_cols; ++col) {
            os << (*this)(row,col);
            if (col != m_cols-1)
                os << ", ";
        }
        if (row != m_rows-1)
            os << std::endl;
    }
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Dump row compression scheme in ostream if it exists
 *
 * @param[in] os Output stream.
 ***************************************************************************/
void GMatrixBase::dump_row_comp(std::ostream& os) const
{
    // If there is a row compression the show scheme
    if (m_rowsel != NULL) {
        os << std::endl << " Row selection ..:";
        for (int row = 0; row < m_num_rowsel; ++row)
            os << " " << m_rowsel[row];
    }
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Dump column compression scheme in ostream if it exists
 *
 * @param[in] os Output stream.
 ***************************************************************************/
void GMatrixBase::dump_col_comp(std::ostream& os) const
{
    // If there is a column compression the show scheme
    if (m_colsel != NULL) {
        os << std::endl << " Column selection:";
        for (int col = 0; col < m_num_colsel; ++col)
            os << " " << m_colsel[col];
    }
  
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
            if ((*this)(row,col) != 0.0)
                break;
        }
        if (col < m_cols)      // found a non-zero element in row
        m_rowsel[m_num_rowsel++] = row;
    }

    // Find all non-zero columns
    for (col = 0; col < m_cols; ++col) {
        for (row = 0; row < m_rows; ++row) {
            if ((*this)(row,col) != 0.0)
                break;
        }
        if (row < m_rows)      // found a non-zero element in column
        m_colsel[m_num_colsel++] = col;
    }
  
    // Return
    return;
}
