/***************************************************************************
 *             GFitsAsciiTable.cpp  - FITS ASCII table class               *
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
#include "GFitsAsciiTable.hpp"
#include "GFitsTableStrCol.hpp"
#include "GFitsTableShtCol.hpp"
#include "GFitsTableLngCol.hpp"
#include "GFitsTableFltCol.hpp"
#include "GFitsTableDblCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */



/*==========================================================================
 =                                                                         =
 =                  GFitsAsciiTable constructors/destructors               =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                                Constructor                              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable() : GFitsData()
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
GFitsAsciiTable::GFitsAsciiTable(const GFitsAsciiTable& table) : GFitsData(table)
{
    // Initialise class members for clean destruction
    init_members();

    // Copy members
    copy_members(table);

    // Return
    return;
}


/***************************************************************************
 *                               Destructor                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsAsciiTable::~GFitsAsciiTable()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                        GFitsAsciiTable operators                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                            Assignment operator                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsAsciiTable& GFitsAsciiTable::operator= (const GFitsAsciiTable& table)
{
    // Execute only if object is not identical
    if (this != &table) {
  
        // Copy base class members
        this->GFitsData::operator=(table);

        // Free members
        free_members();
  
        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(table);
	
    } // endif: object was not identical

    // Return this object
    return *this;
}


/*==========================================================================
 =                                                                         =
 =                      GFitsAsciiTable public methods                     =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                               Open Table                                *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
/// NOT YET IMPLEMENTED
void GFitsAsciiTable::open(__fitsfile* fptr)
{
    cout << "open ASCII table" << endl;
    // Return
    return;
}


/***************************************************************************
 *                               Close Table                               *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsAsciiTable::close(void)
{
    // Free members
    free_members();
  
    // Initialise members
    init_members();
    
    // Return
    return;
}


/***************************************************************************
 *               Return pointer to column with specified name              *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableCol* GFitsAsciiTable::column(const std::string colname)
{
    // Initialise pointer
    GFitsTableCol* ptr = NULL;
    
    // If there are columns then search for the specified name
    if (m_columns != NULL) {
        for (int i = 0; i < m_cols; ++i) {
            if (m_columns[i]->name() == colname) {
                ptr = m_columns[i];
                break;
            }
        }
    }
    
    // Return column pointer
    return ptr;
}


/***************************************************************************
 *                   Return pointer to column with number                  *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
GFitsTableCol* GFitsAsciiTable::column(const int colnum)
{
    // Initialise pointer
    GFitsTableCol* ptr = NULL;
    
    // If there are columns then search for the specified name
    if (m_columns != NULL) {
        ptr = m_columns[colnum];
    }
    
    // Return column pointer
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                      GFitsAsciiTable private methods                    =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                         Initialise class members                        *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsAsciiTable::init_members(void)
{
    // Initialise members
    m_rows    = 0;
    m_cols    = 0;
    m_columns = NULL;

    // Return
    return;
}


/***************************************************************************
 *                            Copy class members                           *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsAsciiTable::copy_members(const GFitsAsciiTable& table)
{
    // Copy attributes
    m_rows = table.m_rows;
    m_cols = table.m_cols;
    
    // Copy column definition
    if (table.m_columns != NULL && m_cols > 0) {
        m_columns = new GFitsTableCol*[m_cols];
        for (int i = 0; i < m_cols; ++i)
            m_columns[i] = table.m_columns[i]->clone();
    }
    
    // Return
    return;
}


/***************************************************************************
 *                           Delete class members                          *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
void GFitsAsciiTable::free_members(void)
{
    // Free memory
    for (int i = 0; i < m_cols; ++i) {
        if (m_columns[i] != NULL) delete m_columns[i];
    }
    if (m_columns != NULL) delete [] m_columns;
    
    // Mark memory as freed
    m_columns = NULL;
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsAsciiTable friends                        =
 =                                                                         =
 ==========================================================================*/

/***************************************************************************
 *                             Output operator                             *
 * ----------------------------------------------------------------------- *
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFitsAsciiTable& table)
{
    // Put header in stream
    os << "=== GFitsAsciiTable ===" << endl;
    os << " Number of rows ............: " << table.m_rows << endl;
    os << " Number of columns .........: " << table.m_cols << endl;
    for (int i = 0; i < table.m_cols; ++i)
        os << " " << *(table.m_columns[i]);

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                  Other functions used by GFitsAsciiTable                =
 =                                                                         =
 ==========================================================================*/
