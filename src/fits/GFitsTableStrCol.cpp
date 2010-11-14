/***************************************************************************
 *          GFitsTableStrCol.cpp  - FITS table string column class         *
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
#include <iostream>
#include <string.h>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsTableStrCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_STRING           "GFitsTableStrCol::string(const int&, const int&)"
#define G_REAL               "GFitsTableStrCol::real(const int&, const int&)"
#define G_INTEGER         "GFitsTableStrCol::integer(const int&, const int&)"

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
GFitsTableStrCol::GFitsTableStrCol() : GFitsTableCol()
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
 * @param[in] width Length of individual string.
 * @param[in] size Number of strings in each column.
 ***************************************************************************/
GFitsTableStrCol::GFitsTableStrCol(const std::string& name,
                                   const int&         length,
                                   const int&         width,
                                   const int&         size)
                                   : GFitsTableCol(name, length, size, width)
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
GFitsTableStrCol::GFitsTableStrCol(const GFitsTableStrCol& column) 
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
GFitsTableStrCol::~GFitsTableStrCol()
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
 * @param[in] column Column which should be assigned
 ***************************************************************************/
GFitsTableStrCol& GFitsTableStrCol::operator= (const GFitsTableStrCol& column)
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
 * Provides access to data in a column. No range checking is performed.
 * Use one of
 *   GFitsTableStrCol::string(ix,iy),
 *   GFitsTableStrCol::real(ix;iy) or
 *   GFitsTableStrCol::integer(ix;iy)
 * if range checking is required.
 ***************************************************************************/
std::string& GFitsTableStrCol::operator() (const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Calculate pixel offset
    int offset = row * m_number + inx;

    // Return data bin
    return m_data[offset];
}


/***********************************************************************//**
 * @brief Column data access operator (const variant)
 *
 * @param[in] row Row of column to access.
 * @param[in] inx Vector index in column row to access
 *
 * Provides access to data in a column. No range checking is performed.
 * Use one of
 *   GFitsTableStrCol::string(ix,iy),
 *   GFitsTableStrCol::real(ix;iy) or
 *   GFitsTableStrCol::integer(ix;iy)
 * if range checking is required.
 ***************************************************************************/
const std::string& GFitsTableStrCol::operator() (const int& row, const int& inx)
                                                 const
{
    // If data are not available then load them now
    if (m_data == NULL) ((GFitsTableStrCol*)this)->fetch_data();

    // Calculate pixel offset
    int offset = row * m_number + inx;

    // Return data bin
    return m_data[offset];
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
 * @exception GException::out_of_range
 *            Table row or vector index are out of valid range.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableStrCol::string(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_STRING, row, 0, m_length-1);

    // Check inx value
    if (inx < 0 || inx >= m_number)
        throw GException::out_of_range(G_STRING, inx, 0, m_number-1);

    // Get index
    int offset = row * m_number + inx;

    // Return value
    return m_data[offset];
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * @exception GException::out_of_range
 *            Table row or vector index are out of valid range.
 *
 * Returns value of specified row and vector index as double precision.
 ***************************************************************************/
double GFitsTableStrCol::real(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_REAL, row, 0, m_length-1);

    // Check inx value
    if (inx < 0 || inx >= m_number)
        throw GException::out_of_range(G_REAL, inx, 0, m_number-1);

    // Get index
    int offset = row * m_number + inx;

    // Assign C string to double
    double value = atof(m_data[offset].c_str());

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Get integer value
 *
 * @param[in] row Table row.
 * @param[in] inx Table column vector index.
 *
 * @exception GException::out_of_range
 *            Table row or vector index are out of valid range.
 *
 * Returns value of specified row and vector index as integer.
 ***************************************************************************/
int GFitsTableStrCol::integer(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_INTEGER, row, 0, m_length-1);

    // Check inx value
    if (inx < 0 || inx >= m_number)
        throw GException::out_of_range(G_INTEGER, inx, 0, m_number-1);

    // Get index
    int offset = row * m_number + inx;

    // Assign C string to int
    int value = int(atof(m_data[offset].c_str()));

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Return pointer to std::string column
 ***************************************************************************/
std::string* GFitsTableStrCol::data(void)
{
    return m_data;
}



/***********************************************************************//**
 * @brief Set NULL string
 *
 * @param[in] value Pointer on NULL string
 *
 * Allows the specification of the FITS table NULL string. If the string
 * is empty the data will not be screened for NULL values.
 ***************************************************************************/
void GFitsTableStrCol::set_nullval(const std::string& string)
{
    // If NULL string is empty then reset the NULL string
    if (string.empty()) {
        if (m_nulval != NULL) delete [] m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy string into NULL string
    else {
        if (m_nulval != NULL) delete [] m_nulval;
        m_nulval = new char[m_width+1];
        strncpy(m_nulval, string.c_str(), m_width);
    }

    // Re-load column
/*
    if (m_data != NULL) {
        save();
        load();
    }
*/

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
void GFitsTableStrCol::init_members(void)
{
    // Initialise members
    m_type   = __TSTRING;
    m_data   = NULL;
    m_buffer = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Column for which members should be copied.
 ***************************************************************************/
void GFitsTableStrCol::copy_members(const GFitsTableStrCol& column)
{
    // Fetch column data if not yet fetched. The casting circumvents the
    // const correctness
    if (column.m_data == NULL)
        ((GFitsTableStrCol*)(&column))->fetch_data();

    // Copy attributes
    m_type = column.m_type;
    m_size = column.m_size;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        alloc_data();
        for (int i = 0; i < m_size; ++i)
            m_data[i] = column.m_data[i];
    }

    // Copy nulval data
    if (column.m_nulval != NULL) {
        m_nulval = new char[m_width+1];
        strncpy(m_nulval, column.m_nulval, m_width);
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableStrCol::free_members(void)
{
    // Free memory
    if (m_data   != NULL) delete [] m_data;
    if (m_nulval != NULL) delete m_nulval;

    // Mark memory as freed
    m_data   = NULL;
    m_nulval = NULL;

    // Free buffer
    free_buffer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * The table column is only saved if it is linked to a FITS file and if the
 * data are indeed present in the class instance. This avoids saving of data
 * that have not been modified.
 *
 * Refer to GFitsTableCol::save_column() for more information.
 ***************************************************************************/
void GFitsTableStrCol::save(void)
{
    // Free buffer
    free_buffer();

    // Allocate buffer
    alloc_buffer();

    // Transfer string into buffer
    for (int i = 0; i < m_size; ++i)
       strncpy(m_buffer[i], m_data[i].c_str(), m_width);

    // Save column
    save_column();

    // Free buffer
    free_buffer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone column
 ***************************************************************************/
GFitsTableStrCol* GFitsTableStrCol::clone(void) const
{
    return new GFitsTableStrCol(*this);
}


/***********************************************************************//**
 * @brief Returns format string of ASCII table
 ***************************************************************************/
std::string GFitsTableStrCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("A");

    // Set width
    format.append(str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Returns format string of binary table
 ***************************************************************************/
std::string GFitsTableStrCol::binary_format(void) const
{
    // Initialize format string
    std::string format;

    // Set number of elements
    format.append(str(m_repeat));

    // Set type code
    format.append("A");

    // If there are substrings then add width of substring
    if (m_repeat > m_width)
        format.append(str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Allocates column data
 ***************************************************************************/
void GFitsTableStrCol::alloc_data(void)
{
    // Free memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free
    m_data = NULL;

    // Allocate new data
    if (m_size > 0)
        m_data = new std::string[m_size];

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise column data
 ***************************************************************************/
void GFitsTableStrCol::init_data(void)
{
    // Initialise data if they exist
    if (m_data != NULL) {
        for (int i = 0; i < m_size; ++i)
            m_data[i].clear();
    }

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
void GFitsTableStrCol::fetch_data(void)
{
    // Calculate size of memory
    m_size = m_number * m_length;

    // Free old buffer memory
    free_buffer();

    // Allocate buffer memory
    alloc_buffer();

    // Save column
    load_column();

    // Extract string from buffer
    for (int i = 0; i < m_size; ++i) {
       if (m_buffer[i] != NULL)
           m_data[i].assign(m_buffer[i]);
    }

    // Free buffer memory
    free_buffer();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Allocate buffer memory
 ***************************************************************************/
void GFitsTableStrCol::alloc_buffer(void)
{
    // Allocate buffer memory
    if (m_size > 0) {
        m_buffer = new char*[m_size];
        for (int i = 0; i < m_size; ++i)
            m_buffer[i] = new char[m_width+1];
    }

    // Initialise buffer
    if (m_buffer != NULL) {
        for (int i = 0; i < m_size; ++i) {
            for (int j = 0; j <= m_width; ++j)
                (m_buffer[i])[j] = '\0';
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Free buffer memory
 ***************************************************************************/
void GFitsTableStrCol::free_buffer(void)
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


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] column Column to put in output stream.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GFitsTableStrCol& column)
{
    // Dump column in output stream
    column.dump_column(os, column.m_data);

    // Return output stream
    return os;
}
