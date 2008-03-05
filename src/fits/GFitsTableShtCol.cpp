/***************************************************************************
 *           GFitsTableShtCol.cpp  - FITS table short column class         *
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
#include "GFitsTableShtCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_SAVE       "GFitsTableShtCol::save()"
#define G_STRING     "GFitsTableShtCol::string(const int&, const int&)"
#define G_REAL       "GFitsTableShtCol::real(const int&, const int&)"
#define G_INTEGER    "GFitsTableShtCol::integer(const int&, const int&)"
#define G_PTR_FLOAT  "GFitsTableShtCol::ptr_float()"
#define G_PTR_DOUBLE "GFitsTableShtCol::ptr_double()"
#define G_PTR_SHORT  "GFitsTableShtCol::ptr_short()"
#define G_PTR_LONG   "GFitsTableShtCol::ptr_long()"
#define G_PTR_INT    "GFitsTableShtCol::ptr_int()"
#define G_LOAD       "GFitsTableShtCol::load()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GFitsTableShtCol constructors/destructors              =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableShtCol::GFitsTableShtCol() : GFitsTableCol()
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
GFitsTableShtCol::GFitsTableShtCol(const std::string& name,
                                   const int&         length,
                                   const int&         size)
                                   : GFitsTableCol(name, length, size, 2)
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***************************************************************************
 *                              Copy constructor                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableShtCol::GFitsTableShtCol(const GFitsTableShtCol& column) : GFitsTableCol(column)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(column);

    // Return
    return;
}


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableShtCol::~GFitsTableShtCol()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GFitsTableShtCol operators                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableShtCol& GFitsTableShtCol::operator= (const GFitsTableShtCol& column)
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
 =                      GFitsTableShtCol public methods                    =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Save table column into FITS file
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 ***************************************************************************/
void GFitsTableShtCol::save(void)
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

/***************************************************************************
 *                             Get string value                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsTableShtCol::string(const int& row, const int& col)
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

    // Convert short into string
    ostringstream s_value;
    s_value << m_data[inx];

    // Return value
    return s_value.str();
}


/***************************************************************************
 *                              Get real value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GFitsTableShtCol::real(const int& row, const int& col)
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

    // Convert short into double
    double value = (double)m_data[inx];

    // Return value
    return value;
}


/***************************************************************************
 *                               Get int value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
int GFitsTableShtCol::integer(const int& row, const int& col)
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

    // Convert short into int
    int value = (int)m_data[inx];

    // Return value
    return value;
}


/***************************************************************************
 *                      Access to invalid data pointers                    *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
float* GFitsTableShtCol::ptr_float(void)
{
    throw GException::fits_invalid_type(G_PTR_FLOAT,
                              "No <float> pointer allowed to <short> array");
}
double* GFitsTableShtCol::ptr_double(void)
{
    throw GException::fits_invalid_type(G_PTR_DOUBLE,
                             "No <double> pointer allowed to <short> array");
}
long* GFitsTableShtCol::ptr_long(void)
{
    throw GException::fits_invalid_type(G_PTR_LONG, 
                               "No <long> pointer allowed to <short> array");
}
int* GFitsTableShtCol::ptr_int(void)
{
    throw GException::fits_invalid_type(G_PTR_INT,
                                "No <int> pointer allowed to <short> array");
}


/***************************************************************************
 *                              Set NULL string                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableShtCol::set_nullval(const short* value)
{
    // If NULL value is empty then reset the NULL value
    if (value == NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy value into NULL value
    else {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new short;
        *m_nulval = *value;
    }

    // Re-load column 
    // NOTE: THIS WILL LEAD TO A LOSS OF MODIFICATIONS; ISSUE SAVE BEFORE !!!
    if (m_data != NULL) {
        //save();
        load();
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                      GFitsTableShtCol private methods                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableShtCol::init_members(void)
{
    // Initialise members
    m_type   = __TSHORT;
    m_size   = 0;
    m_anynul = 0;
    m_data   = NULL;
    m_nulval = NULL;

    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableShtCol::copy_members(const GFitsTableShtCol& column)
{
    // Copy attributes
    m_size   = column.m_size;
    m_anynul = column.m_anynul;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        if (m_data != NULL) delete [] m_data;
        m_data = new short[m_size];
        memcpy(m_data, column.m_data, m_size*sizeof(short));
    }

    // Copy NULL value
    if (column.m_nulval != NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new short;
        *m_nulval = *column.m_nulval;
    }

    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableShtCol::free_members(void)
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


/***************************************************************************
 *                             Load column data                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableShtCol::load(void)
{
    // Calculate size of memory
    m_size = m_number * m_length;

    // Allocate memory
    if (m_data != NULL) delete [] m_data;
    m_data = new short[m_size];

    // Load column data
    int status = 0;
    status = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);
    status = __ffgcv(&m_fitsfile, __TSHORT, m_colnum, 1, 1, m_size,
                     m_nulval, m_data, &m_anynul, &status);
    if (status != 0)
        throw GException::fits_error(G_LOAD, status);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns format string of ASCII table
 ***************************************************************************/
std::string GFitsTableShtCol::ascii_format(void) const
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
std::string GFitsTableShtCol::binary_format(void) const
{
    // Initialize format string
    std::string format;

    // Set number of elements
    format.append(str(m_number));

    // Set type code
    format.append("I");

    // Return format
    return format;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsTableShtCol friends                       =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsTableShtCol               =
 =                                                                         =
 ==========================================================================*/
