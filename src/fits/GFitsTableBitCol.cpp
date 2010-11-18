/***************************************************************************
 *           GFitsTableBitCol.cpp  - FITS table Bit column class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

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
 * @param[in] column Column from which class instance should be built.
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
 * @param[in] column Column which should be assigned
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
 * @param[in] inx Vector index in column row to access
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
    // If data are not available then load them now (circumvent const
    // correctness)
    if (m_data == NULL) ((GFitsTableBitCol*)this)->fetch_data();

    // Set any pending Bit (circumvent const correctness)
    ((GFitsTableBitCol*)this)->set_pending();

    // Get Bit (circumvent const correctness)
    ((GFitsTableBitCol*)this)->get_bit(row, inx);

    // Return data bin
    return m_bit_value;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableBitCol::string(const int& row, const int& inx)
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
double GFitsTableBitCol::real(const int& row, const int& inx)
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
int GFitsTableBitCol::integer(const int& row, const int& inx)
{
    // Get Bit value
    bool bit = (*this)(row, inx);

    // Convert bit into double
    int result = (bit) ? 1 : 0;

    // Return result
    return result;
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
    if (not_loaded) ((GFitsTableBitCol*)(&column))->fetch_data();

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
    if (not_loaded) ((GFitsTableBitCol*)(&column))->release_data();
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
 * @brief Clone column
 ***************************************************************************/
GFitsTableBitCol* GFitsTableBitCol::clone(void) const
{
    return new GFitsTableBitCol(*this);
}


/***********************************************************************//**
 * @brief Returns format string of ASCII table
 ***************************************************************************/
std::string GFitsTableBitCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("X");

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
    if (m_size > 0)
        m_data = new unsigned char[m_size];

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
        for (int i = 0; i < m_size; ++i)
            m_data[i] = 0;
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
            if (status != 0)
                throw GException::fits_hdu_not_found(G_LOAD_COLUMN,
                                  (FPTR(m_fitsfile)->HDUposition)+1,
                                  status);

            // Load data 8 Bits at once
            status = __ffgcv(FPTR(m_fitsfile), __TBYTE, m_colnum, 1, 1, m_size,
                             m_nulval, m_data, &m_anynul, &status);
            if (status != 0)
                throw GException::fits_error(G_LOAD_COLUMN, status,
                                  "for column \""+m_name+"\".");
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
        if (status != 0)
            throw GException::fits_hdu_not_found(G_SAVE_COLUMN,
                              (FPTR(m_fitsfile)->HDUposition)+1,
                              status);

        // Save data 8 Bits at once
        status = __ffpcn(FPTR(m_fitsfile), __TBYTE, m_colnum, 1, 1,
                         m_size, m_data, m_nulval, &status);
        if (status != 0)
            throw GException::fits_error(G_SAVE_COLUMN, status);

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
 * @exception GException::out_of_range
 *            Table row or vector index are out of valid range.
 *
 * Set the Bit for boolean data access. Note that this method assumes that
 * the data have already been loaded.
 ***************************************************************************/
void GFitsTableBitCol::get_bit(const int& row, const int& inx)
{
    // Check row value
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_GET_BIT, row, 0, m_length-1);
    #endif

    // Check inx value
    #if defined(G_RANGE_CHECK)
    if (inx < 0 || inx >= m_number)
        throw GException::out_of_range(G_GET_BIT, inx, 0, m_number-1);
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
        if (m_bit_value)
            m_data[m_bit_byte] = m_data[m_bit_byte] | m_bit_mask;
        else
            m_data[m_bit_byte] = m_data[m_bit_byte] & ~m_bit_mask;

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
