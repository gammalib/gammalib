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
#define G_SAVE       "GFitsTableFltCol::save()"
#define G_STRING     "GFitsTableFltCol::string(const int&, const int&)"
#define G_REAL       "GFitsTableFltCol::real(const int&, const int&)"
#define G_INTEGER    "GFitsTableFltCol::integer(const int&, const int&)"
#define G_PTR_FLOAT  "GFitsTableFltCol::ptr_float()"
#define G_PTR_DOUBLE "GFitsTableFltCol::ptr_double()"
#define G_PTR_SHORT  "GFitsTableFltCol::ptr_short()"
#define G_PTR_LONG   "GFitsTableFltCol::ptr_long()"
#define G_PTR_INT    "GFitsTableFltCol::ptr_int()"
#define G_LOAD       "GFitsTableFltCol::load()"

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
GFitsTableFltCol::GFitsTableFltCol(const GFitsTableFltCol& column) :
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


/*==========================================================================
 =                                                                         =
 =                      GFitsTableFltCol public methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 ***************************************************************************/
void GFitsTableFltCol::save(void)
{
    // Continue only if a FITS file is connected
    if (m_fitsfile.Fptr != NULL) {

        // Move to the HDU
        int status = 0;
        status     = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL,
                              &status);
        if (status != 0)
            throw GException::fits_hdu_not_found(G_SAVE, 
                                                 (m_fitsfile.HDUposition)+1,
                                                 status);

        // Save the column data
        // TBD

    } // endif: FITS file was connected

    // Return
    return;
}


/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] col Table column vector index.
 *
 * Returns value of specified row and vector index as string.
 ***************************************************************************/
std::string GFitsTableFltCol::string(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_STRING, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_number)
        throw GException::out_of_range(G_STRING, col, 0, m_number-1);

    // Get index
    int inx = row * m_repeat + col;

    // Convert double into string
    ostringstream s_value;
    s_value << scientific << m_data[inx];

    // Return value
    return s_value.str();
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] col Table column vector index.
 *
 * Returns value of specified row and vector index as double precision.
 ***************************************************************************/
double GFitsTableFltCol::real(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_REAL, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_number)
        throw GException::out_of_range(G_REAL, col, 0, m_number-1);

    // Get index
    int inx = row * m_repeat + col;

    // Convert float into double
    double value = (double)m_data[inx];

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Get integer value
 *
 * @param[in] row Table row.
 * @param[in] col Table column vector index.
 *
 * Returns value of specified row and vector index as integer.
 ***************************************************************************/
int GFitsTableFltCol::integer(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_INTEGER, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_number)
        throw GException::out_of_range(G_INTEGER, col, 0, m_number-1);

    // Get index
    int inx = row * m_repeat + col;

    // Convert float into int
    int value = (int)m_data[inx];

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
    m_size   = 0;
    m_anynul = 0;
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
    m_size   = column.m_size;
    m_anynul = column.m_anynul;

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
 * @brief Load column data
 ***************************************************************************/
void GFitsTableFltCol::load(void)
{
    // Calculate size of memory
    m_size = m_number * m_length;

    // Allocate memory
    if (m_data != NULL) delete [] m_data;
    m_data = new float[m_size];

    // Load column data
    int status = 0;
    status = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);
    status = __ffgcv(&m_fitsfile, __TFLOAT, m_colnum, 1, 1, m_size,
                     m_nulval, m_data, &m_anynul, &status);
    if (status != 0)
        throw GException::fits_error(G_LOAD, status);

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


/*==========================================================================
 =                                                                         =
 =                          GFitsTableFltCol friends                       =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsTableFltCol               =
 =                                                                         =
 ==========================================================================*/
