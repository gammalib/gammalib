/***************************************************************************
 * GFitsTableUShortCol.cpp - FITS table unsigned short integer column class*
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GFitsTableUShortCol.cpp
 * @brief FITS table unsigned short integer column class implementation
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
#include "GFitsTableUShortCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INSERT                    "GFitsTableUShortCol::insert(int&, int&)"
#define G_REMOVE                    "GFitsTableUShortCol::remove(int&, int&)"

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
GFitsTableUShortCol::GFitsTableUShortCol(void) : GFitsTableCol()
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
 * @param[in] nrows Number of rows in column.
 * @param[in] size Vector size of column.
 ***************************************************************************/
GFitsTableUShortCol::GFitsTableUShortCol(const std::string& name,
                                         const int&         nrows,
                                         const int&         size) :
                     GFitsTableCol(name, nrows, size, 2)
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
GFitsTableUShortCol::GFitsTableUShortCol(const GFitsTableUShortCol& column) :
                     GFitsTableCol(column)
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
GFitsTableUShortCol::~GFitsTableUShortCol(void)
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
GFitsTableUShortCol& GFitsTableUShortCol::operator=(const GFitsTableUShortCol& column)
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
 * @param[in] inx Vector index in column row to access
 *
 * Provides access to data in a column.
 ***************************************************************************/
unsigned short& GFitsTableUShortCol::operator()(const int& row, const int& inx)
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
 * @param[in] inx Vector index in column row to access
 *
 * Provides access to data in a column.
 ***************************************************************************/
const unsigned short& GFitsTableUShortCol::operator()(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

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
void GFitsTableUShortCol::clear(void)
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
GFitsTableUShortCol* GFitsTableUShortCol::clone(void) const
{
    return new GFitsTableUShortCol(*this);
}


/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableUShortCol::string(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Return value
    return (gammalib::str(m_data[offset(row,inx)]));
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as double precision.
 ***************************************************************************/
double GFitsTableUShortCol::real(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert unsigned short into double
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
int GFitsTableUShortCol::integer(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Convert unsigned short into int
    int value = (int)m_data[offset(row,inx)];

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Insert rows in column
 *
 * @param[in] row Row after which rows should be inserted (0=first row).
 * @param[in] nrows Number of rows to be inserted.
 *
 * @exception GException::out_of_range
 *            Specified row is invalid.
 *
 * Inserts rows into a FITS table. This implies that all columns will be
 * loaded into memory.
 ***************************************************************************/
void GFitsTableUShortCol::insert(const int& row, const int& nrows)
{
    // Make sure that row is valid
    if (row < 0 || row > m_length) {
        throw GException::out_of_range(G_INSERT, "FITS table row number",
                                       row, m_length+1);
    }
    
    // Continue only if there are rows to be inserted
    if (nrows > 0) {
    
        // If we have no rows yet then simply set the length to the
        // number of rows to be inserted
        if (m_length == 0) {
            m_length = nrows;
            m_size   = m_number * m_length;
            alloc_data();
            init_data();
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
            unsigned short* new_data = new unsigned short[m_size];

            // Compute the number of elements before the insertion point,
            // the number of elements that get inserted, and the total
            // number of elements after the insertion point
            int n_before = m_number * row;
            int n_insert = m_number * nrows;
            int n_after  = m_number * (m_length - row);
        
            // Copy and initialise data
            unsigned short* src = m_data;
            unsigned short* dst = new_data;
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
 * @param[in] row Row after which rows should be removed (0=first row).
 * @param[in] nrows Number of rows to be removed.
 *
 * @exception GException::out_of_range
 *            Specified row is invalid.
 *            Invalid number of rows specified.
 *
 * This method removes rows from a FITS table. This implies that the column
 * will be loaded into memory.
 ***************************************************************************/
void GFitsTableUShortCol::remove(const int& row, const int& nrows)
{
    // Make sure that row is valid
    if (row < 0 || row >= m_length) {
        throw GException::out_of_range(G_REMOVE, "FITS table row number",
                                       row, m_length);
    }
    
    // Make sure that we don't remove beyond the limit
    if (nrows < 0 || nrows > m_length-row) {
        throw GException::out_of_range(G_REMOVE, "Number of FITS table rows",
                                       nrows, m_length-row+1);
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
            unsigned short* new_data = new unsigned short[m_size];

            // Compute the number of elements before the removal point,
            // the number of elements that get removed, and the total
            // number of elements after the removal point
            int n_before = m_number * row;
            int n_remove = m_number * nrows;
            int n_after  = m_number * (length - row);

            // Copy data
            unsigned short* src = m_data;
            unsigned short* dst = new_data;
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
void GFitsTableUShortCol::nulval(const unsigned short* value)
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
void GFitsTableUShortCol::init_members(void)
{
    // Initialise members
    m_type   = __TUSHORT;
    m_data   = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Column.
 *
 * Sets the content of the vector column by copying from another column.
 * If the code is compiled with the small memory option, and if the source
 * column has not yet been loaded, then we only load the column temporarily
 * for copying purposes and release it again once copying is finished.
 ***************************************************************************/
void GFitsTableUShortCol::copy_members(const GFitsTableUShortCol& column)
{
    // Fetch data if necessary
    bool not_loaded = (!column.is_loaded());
    if (not_loaded) {
        column.fetch_data();
    }

    // Copy attributes
    m_type     = column.m_type;
    m_size     = column.m_size;
    m_varlen   = column.m_varlen;
    m_rowstart = column.m_rowstart;

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
        const_cast<GFitsTableUShortCol*>(&column)->release_data();
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableUShortCol::free_members(void)
{
    // Free memory
    if (m_data   != NULL) delete [] m_data;
    if (m_nulval != NULL) delete m_nulval;

    // Mark memory as freed
    m_data   = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns format string of ASCII table
 ***************************************************************************/
std::string GFitsTableUShortCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("I5");

    // Set width
    //format.append(gammalib::str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Allocates column data
 ***************************************************************************/
void GFitsTableUShortCol::alloc_data(void)
{
    // Free any existing memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free
    m_data = NULL;

    // Allocate new data
    if (m_size > 0) {
        m_data = new unsigned short[m_size];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fetch column data
 ***************************************************************************/
void GFitsTableUShortCol::fetch_data(void) const
{
    // Load column (circumvent const correctness)
    const_cast<GFitsTableUShortCol*>(this)->load_column();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Resize column data
 *
 * @param[in] index Start index.
 * @param[in] number Number of elements to add/remove.
 *
 * Adds or removes elements from specified index on. Adding is done if
 * @p number is a positive number, removing if @p number is negative.
 * Note that the method does not change the validity of the arguments.
 * This needs to be done by the client.
 ***************************************************************************/
void GFitsTableUShortCol::resize_data(const int& index, const int& number)
{
    // Continue only if number of elements is not zero
    if (number != 0) {

        // If data are not available then load them now
        if (m_data == NULL) fetch_data();

        // If elements should be removed then do not allocate new memory
        // but just move elements forward and change the logical size of
        // memory. Only if all elements should be removed the memory is
        // released.
        if (number < 0) {
            int             left = index - number;
            unsigned short* dst  = m_data + index;
            unsigned short* src  = m_data + left;
            int             num  = m_size - left;
            for (int i = 0; i < num; ++i) {
                *dst++ = *src++;
            }
            m_size += number;
            if (m_size < 1) {
                release_data();
            }
        }

        // If elements should be added then allocate new memory, copy over
        // the old data and initialise the new elements
        else {
            int left                 = m_size - index;
            m_size                  += number;
            unsigned short* new_data = new unsigned short[m_size];
            unsigned short* dst      = new_data;
            unsigned short* src      = m_data;
            for (int i = 0; i < index; ++i) {
                *dst++ = *src++;
            }
            for (int i = 0; i < number; ++i) {
                *dst++ = 0;
            }
            for (int i = 0; i < left; ++i) {
                *dst++ = *src++;
            }
            if (m_data != NULL) delete [] m_data;
            m_data = new_data;
        }

    } // endif: number was non-zero

    // Return
    return;
}


/***********************************************************************//**
 * @brief Release column data
 ***************************************************************************/
void GFitsTableUShortCol::release_data(void)
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
 ***************************************************************************/
void GFitsTableUShortCol::alloc_nulval(const unsigned short* value)
{
    // Free any existing memory
    if (m_nulval != NULL) delete m_nulval;

    // Mark pointer as free
    m_nulval = NULL;

    // If we have valid value, allocate and set nul value
    if (value != NULL) {
        m_nulval  = new unsigned short;
        *m_nulval = *value;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise column data
 ***************************************************************************/
void GFitsTableUShortCol::init_data(void)
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
