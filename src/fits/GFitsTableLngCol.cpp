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
#define G_PTR_FLOAT  "GFitsTableLngCol::ptr_float()"
#define G_PTR_DOUBLE "GFitsTableLngCol::ptr_double()"
#define G_PTR_SHORT  "GFitsTableLngCol::ptr_short()"
#define G_PTR_LONG   "GFitsTableLngCol::ptr_long()"
#define G_PTR_INT    "GFitsTableLngCol::ptr_int()"
#define G_LOAD       "GFitsTableLngCol::load()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GFitsTableLngCol constructors/destructors              =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableLngCol::GFitsTableLngCol() : GFitsTableCol()
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
GFitsTableLngCol::GFitsTableLngCol(const GFitsTableLngCol& column) : GFitsTableCol(column)
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

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
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


/*==========================================================================
 =                                                                         =
 =                      GFitsTableLngCol public methods                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Get string value                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsTableLngCol::string(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_STRING, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_repeat)
        throw GException::out_of_range(G_STRING, col, 0, m_repeat-1);

    // Get index
    int inx = row * m_repeat + col;

    // Convert long into string
    ostringstream s_value;
    s_value << m_data[inx];

    // Return value
    return s_value.str();
}


/***************************************************************************
 *                              Get real value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GFitsTableLngCol::real(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_REAL, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_repeat)
        throw GException::out_of_range(G_REAL, col, 0, m_repeat-1);

    // Get index
    int inx = row * m_repeat + col;

    // Convert long into double
    double value = (double)m_data[inx];

    // Return value
    return value;
}


/***************************************************************************
 *                               Get int value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
int GFitsTableLngCol::integer(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_INTEGER, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_repeat)
        throw GException::out_of_range(G_INTEGER, col, 0, m_repeat-1);

    // Get index
    int inx = row * m_repeat + col;

    // Convert long into int
    int value = (int)m_data[inx];

    // Return value
    return value;
}


/***************************************************************************
 *                      Access to invalid data pointers                    *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
float* GFitsTableLngCol::ptr_float(void)
{
    throw GException::fits_invalid_type(G_PTR_FLOAT,
                               "No <float> pointer allowed to <long> array");
}
double* GFitsTableLngCol::ptr_double(void)
{
    throw GException::fits_invalid_type(G_PTR_DOUBLE,
                              "No <double> pointer allowed to <long> array");
}
short* GFitsTableLngCol::ptr_short(void)
{
    throw GException::fits_invalid_type(G_PTR_SHORT,
                               "No <short> pointer allowed to <long> array");
}
int* GFitsTableLngCol::ptr_int(void)
{
    throw GException::fits_invalid_type(G_PTR_INT,
                                 "No <int> pointer allowed to <long> array");
}


/***************************************************************************
 *                              Set NULL string                            *
 * ----------------------------------------------------------------------- *
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
 =                      GFitsTableLngCol private methods                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableLngCol::init_members(void)
{
    // Initialise members
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


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
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


/***************************************************************************
 *                             Load column data                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableLngCol::load(void)
{
    // Calculate size of memory
    m_size = m_repeat * m_length;

    // Allocate memory
    if (m_data != NULL) delete [] m_data;
    m_data = new long[m_size];

    // Load column data
    int status = 0;
    status     = __ffgcv(m_fitsfile, __TLONG, m_colnum, 1, 1, m_size,
                         m_nulval, m_data, &m_anynul, &status);
    if (status != 0)
        throw GException::fits_error(G_LOAD, status);

    // Return
    return;
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
