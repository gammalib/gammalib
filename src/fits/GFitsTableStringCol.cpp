/***************************************************************************
 *         GFitsTableStringCol.cpp - FITS table string column class        *
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
 * @file GFitsTableStringCol.cpp
 * @brief FITS table string column class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstring>            // std::strncpy, etc...
#include <string>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"
#include "GFitsTableStringCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INSERT                    "GFitsTableStringCol::insert(int&, int&)"
#define G_REMOVE                    "GFitsTableStringCol::remove(int&, int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                        Constructors/destructors                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsTableStringCol::GFitsTableStringCol(void) : GFitsTableCol()
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
 * @param[in] width Length of individual string.
 * @param[in] size Number of strings in each column.
 ***************************************************************************/
GFitsTableStringCol::GFitsTableStringCol(const std::string& name,
                                         const int&         nrows,
                                         const int&         width,
                                         const int&         size) :
                     GFitsTableCol(name, nrows, size, width)
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
GFitsTableStringCol::GFitsTableStringCol(const GFitsTableStringCol& column) :
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
GFitsTableStringCol::~GFitsTableStringCol(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Operators                               =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] column Table column.
 * @return Table column.
 ***************************************************************************/
GFitsTableStringCol& GFitsTableStringCol::operator=(const GFitsTableStringCol& column)
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
std::string& GFitsTableStringCol::operator()(const int& row, const int& inx)
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
const std::string& GFitsTableStringCol::operator()(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Return data bin
    return m_data[offset(row, inx)];
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GFitsTableStringCol::clear(void)
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
GFitsTableStringCol* GFitsTableStringCol::clone(void) const
{
    return new GFitsTableStringCol(*this);
}


/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableStringCol::string(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Return value
    return m_data[offset(row,inx)];
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as double precision.
 ***************************************************************************/
double GFitsTableStringCol::real(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Assign string to double
    double value = gammalib::todouble(m_data[offset(row,inx)]);

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
int GFitsTableStringCol::integer(const int& row, const int& inx) const
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Assign string to int
    int value = gammalib::toint(m_data[offset(row,inx)]);

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
void GFitsTableStringCol::insert(const int& row, const int& nrows)
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
            std::string* new_data = new std::string[m_size];

            // Compute the number of elements before the insertion point,
            // the number of elements that get inserted, and the total
            // number of elements after the insertion point
            int n_before = m_number * row;
            int n_insert = m_number * nrows;
            int n_after  = m_number * (m_length - row);
        
            // Copy and initialise data
            std::string* src = m_data;
            std::string* dst = new_data;
            for (int i = 0; i < n_before; ++i) {
                *dst++ = *src++;
            }
            for (int i = 0; i < n_insert; ++i) {
                dst->clear();
                dst++;
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
void GFitsTableStringCol::remove(const int& row, const int& nrows)
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
            std::string* new_data = new std::string[m_size];

            // Compute the number of elements before the removal point,
            // the number of elements that get removed, and the total
            // number of elements after the removal point
            int n_before = m_number * row;
            int n_remove = m_number * nrows;
            int n_after  = m_number * (length - row);

            // Copy data
            std::string* src = m_data;
            std::string* dst = new_data;
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
 * @param[in] value Nul value string.
 *
 * Allows the specification of the FITS table nul string. The default nul
 * strings is an empty string which results in NULL entries in case that
 * no data has been written.
 *
 * @todo To correctly reflect the nul value in the data, the column should
 * be reloaded. However, the column may have been changed, so in principle
 * saving is needed. However, we may not want to store the data, hence saving
 * is also not desired. We thus have to develop a method to update the
 * column information for a new nul value in place ...
 ***************************************************************************/
void GFitsTableStringCol::nulval(const std::string& value)
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
 =                            Private methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTableStringCol::init_members(void)
{
    // Initialise members
    m_type   = __TSTRING;
    m_data   = NULL;
    m_buffer = NULL;
    m_nulval = NULL;

    // Initialise nul value
    alloc_nulval("");

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
void GFitsTableStringCol::copy_members(const GFitsTableStringCol& column)
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
        if (m_data != NULL) delete [] m_data;
        m_data = new std::string[m_size];
        for (int i = 0; i < m_size; ++i) {
            m_data[i] = column.m_data[i];
        }
    }

    // Copy NULL value
    if (column.m_nulval != NULL) {
        if (m_nulval != NULL) delete [] m_nulval;
        m_nulval = new char[m_width+1];
        std::strncpy(m_nulval, column.m_nulval, m_width);
    }

    // Small memory option: release column if it was fetch above
    #if defined(G_SMALL_MEMORY)
    if (not_loaded) {
        const_cast<GFitsTableStringCol*>(&column)->release_data();
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableStringCol::free_members(void)
{
    // Free memory
    if (m_data   != NULL) delete [] m_data;
    if (m_nulval != NULL) delete [] m_nulval;

    // Mark memory as freed
    m_data   = NULL;
    m_nulval = NULL;

    // Reset load flag
    m_size = 0;

    // Free buffer
    free_buffer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * The data is transferred from the internal string vector into a character
 * transfer buffer that can be handled by the cfitsio routines.
 *
 * The table column is only saved if it is linked to a FITS file and if the
 * data are indeed present in the class instance. This avoids saving of data
 * that have not been modified.
 *
 * Refer to GFitsTableCol::save_column() for more information.
 ***************************************************************************/
void GFitsTableStringCol::save(void)
{
    // Free buffer
    free_buffer();

    // Allocate buffer. This also fills the entire buffer with '\0'
    // characters.
    alloc_buffer();

    // Transfer string into buffer
    for (int i = 0; i < m_size; ++i) {
        if (m_data[i].length() > 0) {
            std::strncpy(m_buffer[i], m_data[i].c_str(), m_width);
        }
    }

    // Save column
    save_column();

    // Free buffer
    free_buffer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns format string of ASCII table
 ***************************************************************************/
std::string GFitsTableStringCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("A");

    // Set width
    format.append(gammalib::str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Allocates column data
 ***************************************************************************/
void GFitsTableStringCol::alloc_data(void)
{
    // Free memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free
    m_data = NULL;

    // Allocate new data
    if (m_size > 0) {
        m_data = new std::string[m_size];
    }

    // Allocate also buffer
    free_buffer();
    alloc_buffer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Fetch column data
 *
 * If a FITS file is attached to the column the data are loaded into memory
 * from the FITS file. If no FITS file is attached, memory is allocated
 * to hold the column data and all cells are set to 0.
 *
 * Refer to GFitsTableCol::load_column for more information.
 ***************************************************************************/
void GFitsTableStringCol::fetch_data(void) const
{
    // Calculate size of memory
    //m_size = m_number * m_length;

    // Free old buffer memory
    //free_buffer();

    // Allocate buffer memory
    //alloc_buffer();

    // Load column
    const_cast<GFitsTableStringCol*>(this)->load_column();

    // Extract string from buffer
    for (int i = 0; i < m_size; ++i) {
        if (m_buffer[i] != NULL) {
            m_data[i].assign(m_buffer[i]);
        }
    }

    // Free buffer memory
    free_buffer();

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
void GFitsTableStringCol::resize_data(const int& index, const int& number)
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
            int          left = index - number;
            std::string* dst  = m_data + index;
            std::string* src  = m_data + left;
            int          num  = m_size - left;
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
            int left              = m_size - index;
            m_size               += number;
            std::string* new_data = new std::string[m_size];
            std::string* dst      = new_data;
            std::string* src      = m_data;
            for (int i = 0; i < index; ++i) {
                *dst++ = *src++;
            }
            for (int i = 0; i < number; ++i) {
                *dst++ = "";
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
void GFitsTableStringCol::release_data(void)
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
 * @brief Allocate nul value
 *
 * @param[in] value Nul value string.
 ***************************************************************************/
void GFitsTableStringCol::alloc_nulval(const std::string& value)
{
    // Free any existing memory
    if (m_nulval != NULL) delete [] m_nulval;

    // Mark pointer as free
    m_nulval = NULL;

    // Allocate and initialise nul value
    m_nulval = new char[m_width+1];
    for (int j = 0; j <= m_width; ++j) {
        m_nulval[j] = '\0';
    }

    // Set nul value
    std::strncpy(m_nulval, value.c_str(), m_width);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise column data
 ***************************************************************************/
void GFitsTableStringCol::init_data(void)
{
    // Initialise data if they exist
    if (m_data != NULL) {
        for (int i = 0; i < m_size; ++i) {
            m_data[i].clear();
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate CFITSIO transfer buffer
 *
 * The CFITSIO transfer buffer allows transparent conversion from the CFITSIO
 * storage format to a vector of strings.
 ***************************************************************************/
void GFitsTableStringCol::alloc_buffer(void) const
{
    // Allocate buffer memory
    if (m_size > 0) {
        m_buffer = new char*[m_size];
        for (int i = 0; i < m_size; ++i) {
            m_buffer[i] = new char[m_width+1];
        }
    }

    // Initialise buffer
    if (m_buffer != NULL) {
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j <= m_width; ++j) {
                (m_buffer[i])[j] = '\0';
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Free CFITSIO transfer buffer
 *
 * Release memory that has been allocated for the CFITSIO transfer buffer.
 ***************************************************************************/
void GFitsTableStringCol::free_buffer(void) const
{
    // If there was a buffer allocated then free it
    if (m_buffer != NULL) {
        for (int i = 0; i < m_size; ++i) {
            if (m_buffer[i] != NULL) delete [] m_buffer[i];
        }
        delete [] m_buffer;
        m_buffer = NULL;
    }

    // Return
    return;
}
