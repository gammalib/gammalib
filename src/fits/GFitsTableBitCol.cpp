/***************************************************************************
 *           GFitsTableBitCol.cpp  - FITS table Bit column class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @file GFitsTableBitCol.cpp
 * @brief FITS table bit column class implementation
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
#include "GFitsTableBitCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_INSERT                       "GFitsTableBitCol::insert(int&, int&)"
#define G_REMOVE                       "GFitsTableBitCol::remove(int&, int&)"
#define G_LOAD_COLUMN                       "GFitsTableBitCol::load_column()"
#define G_SAVE_COLUMN                       "GFitsTableBitCol::save_column()"
#define G_GET_BIT                      "GFitsTableBitCol::get_bit(int&,int&)"

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
GFitsTableBitCol::GFitsTableBitCol(void) : GFitsTableCol()
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
GFitsTableBitCol::GFitsTableBitCol(const std::string& name,
                                   const int&         length,
                                   const int&         size)
                                   : GFitsTableCol(name, length, size, 1)
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
GFitsTableBitCol::GFitsTableBitCol(const GFitsTableBitCol& column) 
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
GFitsTableBitCol::~GFitsTableBitCol(void)
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
 * @param[in] column Table column.
 ***************************************************************************/
GFitsTableBitCol& GFitsTableBitCol::operator= (const GFitsTableBitCol& column)
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
bool& GFitsTableBitCol::operator() (const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Set any pending Bit
    set_pending();

    // Get Bit
    get_bit(row, inx);

    // Signal that a Bit is pending. We need this here since the non-const
    // operator allows changing the Bit after exiting the method, hence
    // we have to signal that the actual value of 'm_bit_value' could have
    // been modified and needs to be written back into the data array.
    m_bit_pending = true;

    // Return Bit
    return m_bit_value;
}


/***********************************************************************//**
 * @brief Column data access operator (const variant)
 *
 * @param[in] row Row of column to access.
 * @param[in] inx Vector index in column row to access
 *
 * Provides access to data in a column.
 ***************************************************************************/
const bool& GFitsTableBitCol::operator() (const int& row, const int& inx) const
{
    // Circumvent const correctness
    GFitsTableBitCol* ptr = const_cast<GFitsTableBitCol*>(this);
    
    // If data are not available then load them now
    if (m_data == NULL) ptr->fetch_data();

    // Set any pending Bit
    ptr->set_pending();

    // Get Bit
    ptr->get_bit(row, inx);

    // Return data bin
    return m_bit_value;
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
void GFitsTableBitCol::clear(void)
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
GFitsTableBitCol* GFitsTableBitCol::clone(void) const
{
    return new GFitsTableBitCol(*this);
}


/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableBitCol::string(const int& row, const int& inx) const
{
    // Get Bit value
    bool bit = (*this)(row, inx);

    // Convert bit into string
    std::string result = (bit) ? "T" : "F";

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as double precision.
 ***************************************************************************/
double GFitsTableBitCol::real(const int& row, const int& inx) const
{
    // Get Bit value
    bool bit = (*this)(row, inx);

    // Convert bit into double
    double result = (bit) ? 1.0 : 0.0;

    // Return result
    return result;
}


/***********************************************************************//**
 * @brief Get integer value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as integer.
 ***************************************************************************/
int GFitsTableBitCol::integer(const int& row, const int& inx) const
{
    // Get Bit value
    bool bit = (*this)(row, inx);

    // Convert bit into double
    int result = (bit) ? 1 : 0;

    // Return result
    return result;
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
 * This method inserts rows into a FITS table. This implies that the column
 * will be loaded into memory.
 ***************************************************************************/
void GFitsTableBitCol::insert(const int& rownum, const int& nrows)
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

            // Set any pending Bit
            set_pending();

            // Compute new column length
            int length = m_length + nrows;

            // Compute total number of Bits in column
            m_bits = m_number * length;

            // Compute number of Bytes and Bits per row
            m_bytes_per_row = (m_number > 0) ? ((m_number-1) / 8) + 1 : 0;
            m_bits_per_row  = m_bytes_per_row * 8;

            // Compute length of memory array
            m_size = m_bytes_per_row * length;
        
            // Allocate new data to hold the column
            unsigned char* new_data = new unsigned char[m_size];

            // Compute the number of elements before the insertion point,
            // the number of elements that get inserted, and the total
            // number of elements after the insertion point
            int n_before = m_bytes_per_row * rownum;
            int n_insert = m_bytes_per_row * nrows;
            int n_after  = m_bytes_per_row * (m_length - rownum);

            // Copy and initialise data
            unsigned char* src = m_data;
            unsigned char* dst = new_data;
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
void GFitsTableBitCol::remove(const int& rownum, const int& nrows)
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

        // Set any pending Bit
        set_pending();

        // Compute new column length
        int length = m_length - nrows;
        
        // Compute total number of Bits in column
        m_bits = m_number * length;

        // Compute number of Bytes and Bits per row
        m_bytes_per_row = (m_number > 0) ? ((m_number-1) / 8) + 1 : 0;
        m_bits_per_row  = m_bytes_per_row * 8;

        // Compute length of memory array
        m_size = m_bytes_per_row * length;
        
        // If we have rows remaining then allocate new data to hold
        // the column
        if (m_size > 0) {
        
            // Allocate new data to hold the column
            unsigned char* new_data = new unsigned char[m_size];

            // Compute the number of elements before the removal point,
            // the number of elements that get removed, and the total
            // number of elements after the removal point
            int n_before = m_bytes_per_row * rownum;
            int n_remove = m_bytes_per_row * nrows;
            int n_after  = m_bytes_per_row * (length - rownum);

            // Copy data
            unsigned char* src = m_data;
            unsigned char* dst = new_data;
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
void GFitsTableBitCol::nulval(const unsigned char* value)
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
void GFitsTableBitCol::init_members(void)
{
    // Initialise members
    m_type          = __TBIT;
    m_bits          = 0;
    m_bytes_per_row = 0;
    m_bits_per_row  = 0;
    m_data          = NULL;
    m_nulval        = NULL;
    m_bit_pending   = false;
    m_bit_value     = false;
    m_bit_byte      = 0;
    m_bit_mask      = 0;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Column for which members should be copied.
 ***************************************************************************/
void GFitsTableBitCol::copy_members(const GFitsTableBitCol& column)
{
    // Fetch column data if not yet fetched. The casting circumvents the
    // const correctness
    bool not_loaded = (column.m_data == NULL);
    if (not_loaded) {
        const_cast<GFitsTableBitCol*>(&column)->fetch_data();
    }

    // Copy attributes
    m_type          = column.m_type;
    m_size          = column.m_size;
    m_bits          = column.m_bits;
    m_bytes_per_row = column.m_bytes_per_row;
    m_bits_per_row  = column.m_bits_per_row;
    m_bit_pending   = column.m_bit_pending;
    m_bit_value     = column.m_bit_value;
    m_bit_byte      = column.m_bit_byte;
    m_bit_mask      = column.m_bit_mask;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        alloc_data();
        for (int i = 0; i < m_size; ++i)
            m_data[i] = column.m_data[i];
    }

    // Copy NULL value
    alloc_nulval(column.m_nulval);

    // Small memory option: release column if it was fetch above
    #if defined(G_SMALL_MEMORY)
    if (not_loaded) {
        const_cast<GFitsTableBitCol*>(&column)->release_data();
    }
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableBitCol::free_members(void)
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
std::string GFitsTableBitCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("I");

    // Set width
    format.append(str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Returns format string of binary table
 ***************************************************************************/
std::string GFitsTableBitCol::binary_format(void) const
{
    // Initialize format string
    std::string format;

    // Set number of elements
    format.append(str(m_number));

    // Set type code
    format.append("X");

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Allocates column data
 ***************************************************************************/
void GFitsTableBitCol::alloc_data(void)
{
    // Free any existing memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free
    m_data = NULL;

    // Allocate new data
    if (m_size > 0) {
        m_data = new unsigned char[m_size];
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Release column data
 ***************************************************************************/
void GFitsTableBitCol::release_data(void)
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
void GFitsTableBitCol::alloc_nulval(const unsigned char* value)
{
    // Free any existing memory
    if (m_nulval != NULL) delete m_nulval;

    // Mark pointer as free
    m_nulval = NULL;

    // If we have valid value, allocate and set nul value
    if (value != NULL) {
        m_nulval  = new unsigned char;
        *m_nulval = *value;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise column data
 ***************************************************************************/
void GFitsTableBitCol::init_data(void)
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


/***********************************************************************//**
 * @brief Load table column from FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 * @exception GException::fits_error
 *            An error occured while loading column data from FITS file.
 *
 * Load Bit (vector) column into memory by reading 8 Bits at once.
 ***************************************************************************/
void GFitsTableBitCol::load_column(void)
{
    // Compute total number of Bits in column
    m_bits = m_number * m_length;

    // Compute number of Bytes and Bits per row
    m_bytes_per_row = (m_number > 0) ? ((m_number-1) / 8) + 1 : 0;
    m_bits_per_row  = m_bytes_per_row * 8;

    // Compute length of memory array
    m_size = m_bytes_per_row * m_length;

    // Load only if the column has a positive size
    if (m_size > 0) {

        // Allocate and initialise fresh memory
        alloc_data();
        init_data();

        // If a FITS file is attached then load column data from the FITS
        // file
        if (FPTR(m_fitsfile)->Fptr != NULL) {

            // Move to the HDU
            int status = 0;
            status     = __ffmahd(FPTR(m_fitsfile),
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  NULL, &status);
            if (status != 0) {
                throw GException::fits_hdu_not_found(G_LOAD_COLUMN,
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  status);
            }

            // Load data 8 Bits at once
            status = __ffgcv(FPTR(m_fitsfile), __TBYTE, m_colnum, 1, 1, m_size,
                             m_nulval, m_data, &m_anynul, &status);
            if (status != 0) {
                throw GException::fits_error(G_LOAD_COLUMN, status,
                                  "for column \""+m_name+"\".");
            }
        }

    } // endif: column has a positive size

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 * @exception GException::fits_error
 *            Error occured during writing of the column data.
 *
 * Save Bit (vector) column into FITS file by writing 8 Bits at once.
 ***************************************************************************/
void GFitsTableBitCol::save_column(void)
{
    // Continue only if a FITS file is connected and data have been loaded
    if (FPTR(m_fitsfile)->Fptr != NULL && m_colnum > 0 && m_data != NULL) {

        // Set any pending Bit
        set_pending();

        // Move to the HDU
        int status = 0;
        status     = __ffmahd(FPTR(m_fitsfile),
                              (FPTR(m_fitsfile)->HDUposition)+1, NULL,
                              &status);
        if (status != 0) {
            throw GException::fits_hdu_not_found(G_SAVE_COLUMN,
                              (FPTR(m_fitsfile)->HDUposition)+1,
                              status);
        }

        // Save data 8 Bits at once
        status = __ffpcn(FPTR(m_fitsfile), __TBYTE, m_colnum, 1, 1,
                         m_size, m_data, m_nulval, &status);
        if (status != 0) {
            throw GException::fits_error(G_SAVE_COLUMN, status);
        }

    } // endif: FITS file was connected

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get Bit for boolean access
 *
 * @param[in] row Row of column.
 * @param[in] inx Vector index in column row.
 *
 * @exception GException::fits_invalid_row
 *            Table row out of valid range.
 * @exception GException::out_of_range
 *            Table vector index out of valid range.
 *
 * Set the Bit for boolean data access. Note that this method assumes that
 * the data have already been loaded.
 ***************************************************************************/
void GFitsTableBitCol::get_bit(const int& row, const int& inx)
{
    // Check row value
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_length) {
        throw GException::fits_invalid_row(G_GET_BIT, row, m_length-1);
    }
    #endif

    // Check inx value
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_number) {
        throw GException::out_of_range(G_GET_BIT, inx, 0, m_number-1);
    }
    #endif

    // Compute Byte and Bit mask
    m_bit_byte = row * m_bytes_per_row + inx / 8;
    m_bit_mask = 1 << (7 - (inx % 8));

    // Set Bit value
    m_bit_value = (m_data[m_bit_byte] & m_bit_mask);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Set pending Bit
 *
 * Write the pending Bit into the data. Note that this method assumes that
 * the data have already been loaded.
 ***************************************************************************/
void GFitsTableBitCol::set_pending(void)
{
    // Continue only if we have a pending Bit
    if (m_bit_pending) {

        // Set or unset Bit
        if (m_bit_value) {
            m_data[m_bit_byte] = m_data[m_bit_byte] | m_bit_mask;
        }
        else {
            m_data[m_bit_byte] = m_data[m_bit_byte] & ~m_bit_mask;
        }

        // Signal that no more Bit is pending
        m_bit_pending = false;

    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/
