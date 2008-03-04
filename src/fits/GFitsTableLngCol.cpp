/***************************************************************************
 *           GFitsTableLngCol.cpp  - FITS table long column class          *
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
#include "GException.hpp"
#include "GFitsTableLngCol.hpp"
#include <iostream>                           // cout, cerr

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_STRING     "GFitsTableLngCol::string(const int&, const int&)"
#define G_REAL       "GFitsTableLngCol::real(const int&, const int&)"
#define G_INTEGER    "GFitsTableLngCol::integer(const int&, const int&)"
#define G_FETCH_DATA "GFitsTableLngCol::fetch_data()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GFitsTableLngCol constructors/destructors              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsTableLngCol::GFitsTableLngCol() : GFitsTableCol()
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
GFitsTableLngCol::GFitsTableLngCol(const GFitsTableLngCol& column) :
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
GFitsTableLngCol::~GFitsTableLngCol()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GFitsTableLngCol operators                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] column Column which should be assigned
 ***************************************************************************/
GFitsTableLngCol& GFitsTableLngCol::operator= (const GFitsTableLngCol& column)
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
 *   GFitsTableLngCol::string(ix,iy),
 *   GFitsTableLngCol::real(ix;iy) or
 *   GFitsTableLngCol::integer(ix;iy)
 * if range checking is required.
 ***************************************************************************/
long& GFitsTableLngCol::operator() (const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Calculate pixel offset
    int offset = row * m_repeat + inx;

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
 *   GFitsTableLngCol::string(ix,iy),
 *   GFitsTableLngCol::real(ix;iy) or
 *   GFitsTableLngCol::integer(ix;iy)
 * if range checking is required.
 ***************************************************************************/
const long& GFitsTableLngCol::operator() (const int& row, const int& inx)
                                                                        const
{
    // If data are not available then load them now
    if (m_data == NULL) ((GFitsTableLngCol*)this)->fetch_data();

    // Calculate pixel offset
    int offset = row * m_repeat + inx;

    // Return image pixel
    return m_data[offset];
}


/*==========================================================================
 =                                                                         =
 =                      GFitsTableLngCol public methods                    =
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
std::string GFitsTableLngCol::string(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_STRING, row, 0, m_length-1);

    // Check inx value
    if (inx < 0 || inx >= m_repeat)
        throw GException::out_of_range(G_STRING, inx, 0, m_repeat-1);

    // Get index
    int offset = row * m_repeat + inx;

    // Convert long into string
    ostringstream s_value;
    s_value << m_data[offset];

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
double GFitsTableLngCol::real(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_REAL, row, 0, m_length-1);

    // Check inx value
    if (inx < 0 || inx >= m_repeat)
        throw GException::out_of_range(G_REAL, inx, 0, m_repeat-1);

    // Get index
    int offset = row * m_repeat + inx;

    // Convert long into double
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
int GFitsTableLngCol::integer(const int& row, const int& inx)
{
    // If data are not available then load them now
    if (m_data == NULL) fetch_data();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_INTEGER, row, 0, m_length-1);

    // Check inx value
    if (inx < 0 || inx >= m_repeat)
        throw GException::out_of_range(G_INTEGER, inx, 0, m_repeat-1);

    // Get index
    int offset = row * m_repeat + inx;

    // Convert long into int
    int value = (int)m_data[offset];

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Clone column
 ***************************************************************************/
GFitsTableLngCol* GFitsTableLngCol::clone(void) const
{
    return new GFitsTableLngCol(*this);
}


/***********************************************************************//**
 * @brief Return pointer to long integer column
 ***************************************************************************/
long* GFitsTableLngCol::data(void)
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
void GFitsTableLngCol::set_nullval(const long* value)
{
    // If NULL value is empty then reset the NULL value
    if (value == NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy value into NULL value
    else {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new long;
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
 =                      GFitsTableLngCol private methods                   =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTableLngCol::init_members(void)
{
    // Initialise members
    m_type   = __TLONG;
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
void GFitsTableLngCol::copy_members(const GFitsTableLngCol& column)
{
    // Copy attributes
    m_size   = column.m_size;
    m_anynul = column.m_anynul;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        if (m_data != NULL) delete [] m_data;
        m_data = new long[m_size];
        memcpy(m_data, column.m_data, m_size*sizeof(long));
    }

    // Copy NULL value
    if (column.m_nulval != NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new long;
        *m_nulval = *column.m_nulval;
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GFitsTableLngCol::free_members(void)
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
 * @brief Fetch column data
 *
 * @exception GException::fits_error
 *            An error occured while loading column data from FITS file.
 *
 * If a FITS file is attached to the column the data are loaded into memory
 * from the FITS file. If no FITS file is attached, memory is allocated
 * to hold the column data and all cells are set to 0.
 ***************************************************************************/
void GFitsTableLngCol::fetch_data(void)
{
    // Calculate size of memory
    m_size = m_repeat * m_length;

    // Load only if the column has a positive size
    if (m_size > 0) {

        // Allocate fresh memory
        if (m_data != NULL) delete [] m_data;
        m_data = new long[m_size];

        // If a FITS file is attached then load column data from the FITS
        // file
        if (m_fitsfile.Fptr != NULL) {
            int status = 0;
            status = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL,
                              &status);
            status = __ffgcv(&m_fitsfile, __TLONG, m_colnum, 1, 1, m_size,
                             m_nulval, m_data, &m_anynul, &status);
            if (status != 0)
                throw GException::fits_error(G_FETCH_DATA, status);
        }

        // ... otherwise initialise all column values to 0
        else {
            for (int i = 0; i < m_size; ++i)
                m_data[i] = 0;
        }

    } // endif: column has a positive size

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsTableLngCol friends                       =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsTableLngCol               =
 =                                                                         =
 ==========================================================================*/
