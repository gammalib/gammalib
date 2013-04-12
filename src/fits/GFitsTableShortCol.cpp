/***************************************************************************
 *     GFitsTableShortCol.cpp  - FITS table short integer column class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GFitsTableShortCol.cpp
 * @brief FITS table short integer column class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <string>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsTableShortCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INSERT                     "GFitsTableShortCol::insert(int&, int&)"
#define G_REMOVE                     "GFitsTableShortCol::remove(int&, int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsTableShortCol::GFitsTableShortCol(void) : GFitsTableCol()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] name Name of column.
 * @param[in] length Length of column.
 * @param[in] size Vector size of column.
 ***************************************************************************/
GFitsTableShortCol::GFitsTableShortCol(const std::string& name,
                                       const int&         length,
                                       const int&         size)
                                       : GFitsTableCol(name, length, size, 2)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] column Table column.
 ***************************************************************************/
GFitsTableShortCol::GFitsTableShortCol(const GFitsTableShortCol& column)
                                       : GFitsTableCol(column)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(column);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GFitsTableShortCol::~GFitsTableShortCol(void)
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
 * @param[in] column Table column.
 ***************************************************************************/
GFitsTableShortCol& GFitsTableShortCol::operator= (const GFitsTableShortCol& column)
{
    // Execute only if object is not identical
    if (this != &column) {

        // Copy base class members
        this->GFitsTableCol::operator=(column);

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(column);

    } // endif: object was not identical

    // Return this object
    return *this;
}


/***********************************************************************//**
 * @brief Column data access operator
 *
 * @param[in] row Row of column to access.
 * @param[in] inx Vector index in column row to access.
 *
 * Provides access to data in a column.
 ***************************************************************************/
short& GFitsTableShortCol::operator() (const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Return data bin
    return m_data[offset(row, inx)];
}


/***********************************************************************//**
 * @brief Column data access operator (const variant)
 *
 * @param[in] row Row of column to access.
 * @param[in] inx Vector index in column row to access.
 *
 * Provides access to data in a column.
 ***************************************************************************/
const short& GFitsTableShortCol::operator() (const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) {
        const_cast<GFitsTableShortCol*>(this)->fetch_data();
    }

    // Return data bin
    return m_data[offset(row, inx)];
}


/*==========================================================================
 =                                                                         =
 =                               Public methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GFitsTableShortCol::clear(void)
{
    // Free class members (base and derived classes, derived class first)
    free_members();
    this->GFitsTableCol::free_members();

    // Initialise members
    this->GFitsTableCol::init_members();
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone column
 ***************************************************************************/
GFitsTableShortCol* GFitsTableShortCol::clone(void) const
{
    return new GFitsTableShortCol(*this);
}


/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableShortCol::string(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert short into string
    std::ostringstream s_value;
    s_value << m_data[offset(row,inx)];

    // Return value
    return s_value.str();
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as double precision.
 ***************************************************************************/
double GFitsTableShortCol::real(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert short into double
    double value = (double)m_data[offset(row,inx)];

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Get integer value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as integer.
 ***************************************************************************/
int GFitsTableShortCol::integer(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert short into int
    int value = (int)m_data[offset(row,inx)];

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Insert rows in column
 *
 * @param[in] rownum Row after which rows should be inserted (0=first row).
 * @param[in] nrows Number of rows to be inserted.
 *
 * @exception GException::fits_invalid_row
 *            Specified rownum is invalid.
 *
 * This method inserts rows into a FITS table. This implies that all columns
 * will be loaded into memory.
 ***************************************************************************/
void GFitsTableShortCol::insert(const int& rownum, const int& nrows)
{
    // Make sure that rownum is valid
    if (rownum < 0 || rownum > m_length) {
        throw GException::fits_invalid_row(G_INSERT, rownum, m_length);
    }
    
    // Continue only if there are rows to be inserted
    if (nrows > 0) {
    
        // If we have no rows yet then simply set the length to the
        // number of rows to be inserted
        if (m_length == 0) {
            m_length = nrows;
        }
        
        // ... otherwise fetch data, allocate new data and copy over
        // the existing items
        else {

            // If data are not available then load them now
            if (m_data == NULL) fetch_data();

            // Compute new column length
            int length = m_length + nrows;

            // Calculate size of new memory
            m_size = m_number * length;
        
            // Allocate new data to hold the column
            short* new_data = new short[m_size];

            // Compute the number of elements before the insertion point,
            // the number of elements that get inserted, and the total
            // number of elements after the insertion point
            int n_before = m_number * rownum;
            int n_insert = m_number * nrows;
            int n_after  = m_number * (m_length - rownum);
        
            // Copy and initialise data
            short* src = m_data;
            short* dst = new_data;
            for (int i = 0; i < n_before; ++i) {
                *dst++ = *src++;
            }
            for (int i = 0; i < n_insert; ++i) {
                *dst++ = 0;
            }
            for (int i = 0; i < n_after; ++i) {
                *dst++ = *src++;
            }
        
            // Free old data
            if (m_data != NULL) delete [] m_data;
            
            // Set pointer to new data and store length
            m_data   = new_data;
            m_length = length;
        
        } // endelse: there were already data
    
    } // endfor: there were rows to be inserted

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove rows from column
 *
 * @param[in] rownum Row after which rows should be removed (0=first row).
 * @param[in] nrows Number of rows to be removed.
 *
 * @exception GException::fits_invalid_row
 *            Specified rownum is invalid.
 * @exception GException::fits_invalid_nrows
 *            Invalid number of rows specified.
 *
 * This method removes rows from a FITS table. This implies that the column
 * will be loaded into memory.
 ***************************************************************************/
void GFitsTableShortCol::remove(const int& rownum, const int& nrows)
{
    // Make sure that rownum is valid
    if (rownum < 0 || rownum >= m_length) {
        throw GException::fits_invalid_row(G_REMOVE, rownum, m_length-1);
    }
    
    // Make sure that we don't remove beyond the limit
    if (nrows < 0 || nrows > m_length-rownum) {
        throw GException::fits_invalid_nrows(G_REMOVE, nrows, m_length-rownum);
    }
    
    // Continue only if there are rows to be removed
    if (nrows > 0) {
    
        // If data are not available then load them now
        if (m_data == NULL) fetch_data();

        // Compute new column length
        int length = m_length - nrows;
        
        // Calculate size of new memory
        m_size = m_number * length;

        // If we have rows remaining then allocate new data to hold
        // the column
        if (m_size > 0) {
        
            // Allocate new data to hold the column
            short* new_data = new short[m_size];

            // Compute the number of elements before the removal point,
            // the number of elements that get removed, and the total
            // number of elements after the removal point
            int n_before = m_number * rownum;
            int n_remove = m_number * nrows;
            int n_after  = m_number * (length - rownum);

            // Copy data
            short* src = m_data;
            short* dst = new_data;
            for (int i = 0; i < n_before; ++i) {
                *dst++ = *src++;
            }
            src += n_remove;
            for (int i = 0; i < n_after; ++i) {
                *dst++ = *src++;
            }
        
            // Free old data
            if (m_data != NULL) delete [] m_data;
            
            // Set pointer to new data and store length
            m_data   = new_data;
            m_length = length;
        
        } // endif: there are still elements after removal
        
        // ... otherwise just remove all data
        else {

            // Free old data
            if (m_data != NULL) delete [] m_data;

            // Set pointer to new data and store length
            m_data   = NULL;
            m_length = length;
        }
    
    } // endfor: there were rows to be removed

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set nul value
 *
 * @param[in] value Nul value.
 *
 * @todo To correctly reflect the nul value in the data, the column should
 * be reloaded. However, the column may have been changed, so in principle
 * saving is needed. However, we may not want to store the data, hence saving
 * is also not desired. We thus have to develop a method to update the
 * column information for a new nul value in place ...
 ***************************************************************************/
void GFitsTableShortCol::nulval(const short* value)
{
    // Allocate nul value
    alloc_nulval(value);

    // Update column
//    if (m_data != NULL) {
//        save();
//        load();
//    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTableShortCol::init_members(void)
{
    // Initialise members
    m_type   = __TSHORT;
    m_data   = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Table column.
 *
 * Sets the content of the vector column by copying from another column.
 * If the code is compiled with the small memory option, and if the source
 * column has not yet been loaded, then we only load the column temporarily
 * for copying purposes and release it again once copying is finished.
 ***************************************************************************/
void GFitsTableShortCol::copy_members(const GFitsTableShortCol& column)
{
    // Fetch column data if not yet fetched. The casting circumvents the
    // const correctness
    bool not_loaded = (column.m_data == NULL);
    if (not_loaded) {
        const_cast<GFitsTableShortCol*>(&column)->fetch_data();
    }

    // Copy attributes
    m_type = column.m_type;
    m_size = column.m_size;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        alloc_data();
        for (int i = 0; i < m_size; ++i) {
            m_data[i] = column.m_data[i];
        }
    }

    // Copy NULL value
    alloc_nulval(column.m_nulval);

    // Small memory option: release column if it was fetch above
    #if defined(G_SMALL_MEMORY)
    if (not_loaded) {
        const_cast<GFitsTableShortCol*>(&column)->release_data();
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableShortCol::free_members(void)
{
    // Free memory
    if (m_data   != NULL) delete [] m_data;
    if (m_nulval != NULL) delete m_nulval;

    // Mark memory as freed
    m_data   = NULL;
    m_nulval = NULL;

    // Reset load flag
    m_size = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns format string of ASCII table
 ***************************************************************************/
std::string GFitsTableShortCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("I6");

    // Set width
    //format.append(gammalib::str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Returns format string of binary table
 ***************************************************************************/
std::string GFitsTableShortCol::binary_format(void) const
{
    // Initialize format string
    std::string format;

    // Set number of elements
    format.append(gammalib::str(m_number));

    // Set type code
    format.append("I");

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Allocates column data
 ***************************************************************************/
void GFitsTableShortCol::alloc_data(void)
{
    // Free any existing memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free
    m_data = NULL;

    // Allocate new data
    if (m_size > 0) {
        m_data = new short[m_size];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Release column data
 ***************************************************************************/
void GFitsTableShortCol::release_data(void)
{
    // Free any existing memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free and reset loaded vector size
    m_data = NULL;
    m_size = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocates null value
 *
 * @param[in] value Nul value.
 ***************************************************************************/
void GFitsTableShortCol::alloc_nulval(const short* value)
{
    // Free any existing memory
    if (m_nulval != NULL) delete m_nulval;

    // Mark pointer as free
    m_nulval = NULL;

    // If we have valid value, allocate and set nul value
    if (value != NULL) {
        m_nulval  = new short;
        *m_nulval = *value;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise column data
 ***************************************************************************/
void GFitsTableShortCol::init_data(void)
{
    // Initialise data if they exist
    if (m_data != NULL) {
        for (int i = 0; i < m_size; ++i) {
            m_data[i] = 0;
        }
    }

    // Return
    return;
}
