/***************************************************************************
 *          GFitsTableDblCol.cpp  - FITS table double column class         *
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
#include "GFitsTableDblCol.hpp"
#include <iostream>                           // cout, cerr

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_STRING     "GFitsTableDblCol::string(const int&, const int&)"
#define G_REAL       "GFitsTableDblCol::real(const int&, const int&)"
#define G_INTEGER    "GFitsTableDblCol::integer(const int&, const int&)"
#define G_PTR_FLOAT  "GFitsTableDblCol::ptr_float()"
#define G_PTR_DOUBLE "GFitsTableDblCol::ptr_double()"
#define G_PTR_SHORT  "GFitsTableDblCol::ptr_short()"
#define G_PTR_LONG   "GFitsTableDblCol::ptr_long()"
#define G_PTR_INT    "GFitsTableDblCol::ptr_int()"
#define G_LOAD       "GFitsTableDblCol::load()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GFitsTableDblCol constructors/destructors              =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableDblCol::GFitsTableDblCol() : GFitsTableCol()
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
GFitsTableDblCol::GFitsTableDblCol(const GFitsTableDblCol& column) : GFitsTableCol(column)
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
GFitsTableDblCol::~GFitsTableDblCol()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GFitsTableDblCol operators                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableDblCol& GFitsTableDblCol::operator= (const GFitsTableDblCol& column)
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
 =                      GFitsTableDblCol public methods                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Get string value                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsTableDblCol::string(const int& row, const int& col)
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

    // Convert double into string
    ostringstream s_value;
    s_value << scientific << m_data[inx];

    // Return value
    return s_value.str();
}


/***************************************************************************
 *                              Get real value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GFitsTableDblCol::real(const int& row, const int& col)
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

    // Return value
    return m_data[inx];
}


/***************************************************************************
 *                               Get int value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
int GFitsTableDblCol::integer(const int& row, const int& col)
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

    // Convert double into int
    int value = (int)m_data[inx];

    // Return value
    return value;
}


/***************************************************************************
 *                      Access to invalid data pointers                    *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
float* GFitsTableDblCol::ptr_float(void)
{
    throw GException::fits_invalid_type(G_PTR_FLOAT,
                             "No <float> pointer allowed to <double> array");
}
short* GFitsTableDblCol::ptr_short(void) 
{
    throw GException::fits_invalid_type(G_PTR_SHORT,
                             "No <short> pointer allowed to <double> array");
}
long* GFitsTableDblCol::ptr_long(void) 
{
    throw GException::fits_invalid_type(G_PTR_LONG,
                              "No <long> pointer allowed to <double> array");
}
int* GFitsTableDblCol::ptr_int(void) 
{
    throw GException::fits_invalid_type(G_PTR_INT,
                               "No <int> pointer allowed to <double> array");
}


/***************************************************************************
 *                              Set NULL string                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableDblCol::set_nullval(const double* value)
{
    // If NULL value is empty then reset the NULL value
    if (value == NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval = NULL;
    }

    // ... otherwise copy value into NULL value
    else {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new double;
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
 =                      GFitsTableDblCol private methods                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableDblCol::init_members(void)
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
void GFitsTableDblCol::copy_members(const GFitsTableDblCol& column)
{
    // Copy attributes
    m_size   = column.m_size;
    m_anynul = column.m_anynul;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        if (m_data != NULL) delete [] m_data;
        m_data = new double[m_size];
        memcpy(m_data, column.m_data, m_size*sizeof(double));
    }

    // Copy NULL value
    if (column.m_nulval != NULL) {
        if (m_nulval != NULL) delete m_nulval;
        m_nulval  = new double;
        *m_nulval = *column.m_nulval;
    }

    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableDblCol::free_members(void)
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
void GFitsTableDblCol::load(void)
{
    // Calculate size of memory
    m_size = m_repeat * m_length;

    // Allocate memory
    if (m_data != NULL) delete [] m_data;
    m_data = new double[m_size];

    // Load column data
    int status = 0;
    status = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);
    status = __ffgcv(&m_fitsfile, __TDOUBLE, m_colnum, 1, 1, m_size,
                     m_nulval, m_data, &m_anynul, &status);
    if (status != 0)
        throw GException::fits_error(G_LOAD, status);

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsTableDblCol friends                       =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsTableDblCol               =
 =                                                                         =
 ==========================================================================*/
