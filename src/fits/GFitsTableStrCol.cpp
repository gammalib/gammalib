/***************************************************************************
 *          GFitsTableStrCol.cpp  - FITS table string column class         *
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
#include "GFitsTableStrCol.hpp"
#include <iostream>                           // cout, cerr

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_STRING     "GFitsTableStrCol::string(const int&, const int&)"
#define G_REAL       "GFitsTableStrCol::real(const int&, const int&)"
#define G_INTEGER    "GFitsTableStrCol::integer(const int&, const int&)"
#define G_PTR_STRING "GFitsTableStrCol::ptr_string()"
#define G_PTR_FLOAT  "GFitsTableStrCol::ptr_float()"
#define G_PTR_DOUBLE "GFitsTableStrCol::ptr_double()"
#define G_PTR_SHORT  "GFitsTableStrCol::ptr_short()"
#define G_PTR_LONG   "GFitsTableStrCol::ptr_long()"
#define G_PTR_INT    "GFitsTableStrCol::ptr_int()"
#define G_LOAD       "GFitsTableStrCol::load()"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                  GFitsTableStrCol constructors/destructors              =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableStrCol::GFitsTableStrCol() : GFitsTableCol()
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
GFitsTableStrCol::GFitsTableStrCol(const GFitsTableStrCol& column) : GFitsTableCol(column)
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
GFitsTableStrCol::~GFitsTableStrCol()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                       GFitsTableStrCol operators                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
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


/*==========================================================================
 =                                                                         =
 =                      GFitsTableStrCol public methods                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Get string value                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
std::string GFitsTableStrCol::string(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_STRING, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_num_subs)
        throw GException::out_of_range(G_STRING, col, 0, m_num_subs-1);

    // Get string index
    int inx = row * m_num_subs + col;

    // Assign C string to C++ std::string
    std::string value;
    value.assign(m_data[inx]);

    // Return value
    return value;
}


/***************************************************************************
 *                              Get real value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
double GFitsTableStrCol::real(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_STRING, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_num_subs)
        throw GException::out_of_range(G_STRING, col, 0, m_num_subs-1);

    // Get string index
    int inx = row * m_num_subs + col;

    // Assign C string to double
    double value = atof(m_data[inx]);

    // Return value
    return value;
}


/***************************************************************************
 *                               Get int value                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
int GFitsTableStrCol::integer(const int& row, const int& col)
{
    // Make sure that data are loaded
    if (m_data == NULL)
        load();

    // Check row value
    if (row < 0 || row >= m_length)
        throw GException::out_of_range(G_STRING, row, 0, m_length-1);

    // Check col value
    if (col < 0 || col >= m_num_subs)
        throw GException::out_of_range(G_STRING, col, 0, m_num_subs-1);

    // Get string index
    int inx = row * m_num_subs + col;

    // Assign C string to int
    int value = atoi(m_data[inx]);

    // Return value
    return value;
}


/***************************************************************************
 *                      Access to invalid data pointers                    *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
float* GFitsTableStrCol::ptr_float(void) 
{
    throw GException::fits_invalid_type(G_PTR_FLOAT, "No <float> pointer allowed to <string> array");
}
double* GFitsTableStrCol::ptr_double(void) 
{
    throw GException::fits_invalid_type(G_PTR_FLOAT, "No <double> pointer allowed to <string> array");
}
short* GFitsTableStrCol::ptr_short(void) 
{
    throw GException::fits_invalid_type(G_PTR_SHORT, "No <short> pointer allowed to <string> array");
}
long* GFitsTableStrCol::ptr_long(void) 
{
    throw GException::fits_invalid_type(G_PTR_LONG, "No <long> pointer allowed to <string> array");
}
int* GFitsTableStrCol::ptr_int(void) 
{
    throw GException::fits_invalid_type(G_PTR_INT, "No <int> pointer allowed to <string> array");
}


/***************************************************************************
 *                              Set NULL string                            *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableStrCol::set_nullstr(const std::string string)
{
    // If NULL string is empty then reset the NULL string
    if (string.empty()) {
        if (m_nulstr != NULL) delete [] m_nulstr;
        m_nulstr = NULL;
    }

    // ... otherwise copy string into NULL string
    else {
        if (m_nulstr != NULL) delete [] m_nulstr;
        m_nulstr = new char[m_width+1];
        strncpy(m_nulstr, string.c_str(), m_width);
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
 =                      GFitsTableStrCol private methods                   =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableStrCol::init_members(void)
{
    // Initialise members
    m_type   = __TSTRING;
    m_size     = 0;
    m_num_subs = 0;
    m_data     = NULL;
    m_nulstr   = NULL;
    m_anynul   = 0;

    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 * Note: it is assumed that no memory has been allocated so far for the    *
 *       actual instance.                                                  *
 ***************************************************************************/
void GFitsTableStrCol::copy_members(const GFitsTableStrCol& column)
{
    // Copy string size
    m_size     = column.m_size;
    m_num_subs = column.m_num_subs;
    m_anynul   = column.m_anynul;

    // Copy column data
    if (column.m_data != NULL && m_size > 0) {
        m_data = new char*[m_size];
        for (int i = 0; i < m_size; ++i) {
            if (column.m_data[i] != NULL) {
                m_data[i] = new char[m_width+1];
                memcpy(m_data[i], column.m_data[i], m_width);
            }
            else
                m_data[i] = NULL;
        }
    }

    // Copy nulstr data
    if (column.m_nulstr != NULL) {
        m_nulstr = new char[m_width+1];
        strncpy(m_nulstr, column.m_nulstr, m_width);
    }

    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsTableStrCol::free_members(void)
{
    // Free memory
    if (m_data != NULL) {
        for (int i = 0; i < m_size; ++i) {
            if (m_data[i] != NULL) delete [] m_data[i];
        }
        delete [] m_data;
    }

    // Mark memory as freed
    m_data = NULL;

    // Return
    return;
}


/***************************************************************************
 *                             Load column data                            *
 * ----------------------------------------------------------------------- *
 * The string column is read in a list of character buffers. The length of *
 * each character buffer is 'm_width+1' to hold the final '\0'. The number *
 * of strings is given by the number of rows of the table times the number *
 * of substrings.                                                          *
 ***************************************************************************/
void GFitsTableStrCol::load(void)
{
    // Calculate number of substrings
    if (m_repeat == 1)             // ASCII tables
        m_num_subs = 1;
    else                           // Binary tables
        m_num_subs = m_repeat / m_width;

    // Calculate total number of strings
    m_size = m_num_subs * m_length;

    // Free memory
    free_members();

    // Reset number of NULLs
    m_anynul = 0;

    // Load data only if some are available
    if (m_size > 0) {

        // Allocate memory
        m_data = new char*[m_size];
        for (int i = 0; i < m_size; ++i)
            m_data[i] = new char[m_width+1];

        // Load column data
        int status = 0;
        status = __ffmahd(&m_fitsfile, (m_fitsfile.HDUposition)+1, NULL, &status);
        status = __ffgcvs(&m_fitsfile, m_colnum, 1, 1, m_size, m_nulstr, m_data, 
                          &m_anynul, &status);
        if (status != 0)
            throw GException::fits_error(G_LOAD, status);

    } // endif: data were available

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsTableStrCol friends                       =
 =                                                                         =
 ==========================================================================*/


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsTableStrCol               =
 =                                                                         =
 ==========================================================================*/
