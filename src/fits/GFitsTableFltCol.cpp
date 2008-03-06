/***************************************************************************
 *           GFitsTableFltCol.cpp  - FITS table float column class         *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2008 by Jurgen Knodlseder                   *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsTableFltCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_STRING  "GFitsTableFltCol::string(const int&, const int&)"
#define G_REAL    "GFitsTableFltCol::real(const int&, const int&)"
#define G_INTEGER "GFitsTableFltCol::integer(const int&, const int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GFitsTableFltCol constructors/destructors              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsTableFltCol::GFitsTableFltCol() : GFitsTableCol()
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
GFitsTableFltCol::GFitsTableFltCol(const std::string& name,
                                   const int&         length,
                                   const int&         size)
                                   : GFitsTableCol(name, length, size, 4)
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
GFitsTableFltCol::GFitsTableFltCol(const GFitsTableFltCol& column)
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
GFitsTableFltCol::~GFitsTableFltCol()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GFitsTableFltCol operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] column Column which should be assigned
 ***************************************************************************/
GFitsTableFltCol& GFitsTableFltCol::operator= (const GFitsTableFltCol& column)
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
 *   GFitsTableFltCol::string(ix,iy),
 *   GFitsTableFltCol::real(ix;iy) or
 *   GFitsTableFltCol::integer(ix;iy)
 * if range checking is required.
 ***************************************************************************/
float& GFitsTableFltCol::operator() (const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Calculate pixel offset
    int offset = row * m_number + inx;

    // Return image pixel
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
 *   GFitsTableFltCol::string(ix,iy),
 *   GFitsTableFltCol::real(ix;iy) or
 *   GFitsTableFltCol::integer(ix;iy)
 * if range checking is required.
 ***************************************************************************/
const float& GFitsTableFltCol::operator() (const int& row, const int& inx)
                                           const
{
    // If data are not available then load them now
    if (m_data == NULL) ((GFitsTableFltCol*)this)->fetch_data();

    // Calculate pixel offset
    int offset = row * m_number + inx;

    // Return image pixel
    return m_data[offset];
}


/*==========================================================================
 =                                                                         =
 =                      GFitsTableFltCol public methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * The table column is only saved if it is linked to a FITS file and if the
 * data are indeed present in the class instance. This avoids saving of data
 * that have not been modified.
 *
 * Refer to GFitsTableCol::save_column() for more information.
 ***************************************************************************/
void GFitsTableFltCol::save(void)
{
    // Save column
    save_column();

    // Return
    return;
}


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
std::string GFitsTableFltCol::string(const int& row, const int& inx)
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

    // Convert double into string
    ostringstream s_value;
    s_value << scientific << m_data[offset];

    // Return value
    return s_value.str();
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
double GFitsTableFltCol::real(const int& row, const int& inx)
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

    // Convert float into double
    double value = (double)m_data[offset];

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
int GFitsTableFltCol::integer(const int& row, const int& inx)
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

    // Convert float into int
    int value = (int)m_data[offset];

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Clone column
 ***************************************************************************/
GFitsTableFltCol* GFitsTableFltCol::clone(void) const
{
    return new GFitsTableFltCol(*this);
}


/***********************************************************************//**
 * @brief Return pointer to floating point column
 ***************************************************************************/
float* GFitsTableFltCol::data(void)
{
    return m_data;
}


/***********************************************************************//**
 * @brief Set NULL value
 *
 * @param[in] value Pointer on NULL value
 *
 * Allows the specification of the FITS table NULL value. If value=NULL the
 * data will not be screened for NULL values.
 ***************************************************************************/
void GFitsTableFltCol::set_nullval(const float* value)
{
    // If NULL value is empty then reset the NULL value
    if (value == NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy value into NULL value
    else {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new float;
        *m_nulval = *value;
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
 =                      GFitsTableFltCol private methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTableFltCol::init_members(void)
{
    // Initialise members
    m_type   = __TFLOAT;
    m_data   = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] column Column for which members should be copied.
 ***************************************************************************/
void GFitsTableFltCol::copy_members(const GFitsTableFltCol& column)
{
    // Copy attributes

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        if (m_data != NULL) delete [] m_data;
        m_data = new float[m_size];
        memcpy(m_data, column.m_data, m_size*sizeof(float));
    }

    // Copy NULL value
    if (column.m_nulval != NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new float;
        *m_nulval = *column.m_nulval;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableFltCol::free_members(void)
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
std::string GFitsTableFltCol::ascii_format(void) const
{
    // Initialize format string
    std::string format;

    // Set type code
    format.append("F");

    // Set width
    format.append(str(m_width));

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Returns format string of binary table
 ***************************************************************************/
std::string GFitsTableFltCol::binary_format(void) const
{
    // Initialize format string
    std::string format;

    // Set number of elements
    format.append(str(m_number));

    // Set type code
    format.append("E");

    // Return format
    return format;
}


/***********************************************************************//**
 * @brief Allocates column data
 ***************************************************************************/
void GFitsTableFltCol::alloc_data(void)
{
    // Free any existing memory
    if (m_data != NULL) delete [] m_data;

    // Mark pointer as free
    m_data = NULL;

    // Allocate new data
    if (m_size > 0)
        m_data = new float[m_size];

    // Return
    return;
}


/***********************************************************************//**
 * @brief Initialise column data
 ***************************************************************************/
void GFitsTableFltCol::init_data(void)
{
    // Initialise data if they exist
    if (m_data != NULL) {
        for (int i = 0; i < m_size; ++i)
            m_data[i] = 0.0;
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
void GFitsTableFltCol::fetch_data(void)
{
    // Save column
    load_column();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsTableFltCol friends                       =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] column Column to put in output stream
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFitsTableFltCol& column)
{
    // Dump column in output stream
    column.dump_column(os, column.m_data);

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsTableFltCol               =
 =                                                                         =
 ==========================================================================*/
