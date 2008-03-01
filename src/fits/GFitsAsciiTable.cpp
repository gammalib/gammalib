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

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsAsciiTable::GFitsAsciiTable() : GFitsData()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param table Table to be used to initialise instance
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


/***********************************************************************//**
 * @brief Destructor
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

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param table Table to be assigned
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

/***********************************************************************//**
 * @brief Open Table
 *
 * NOT YET IMPLEMENTED
 ***************************************************************************/
void GFitsAsciiTable::open(__fitsfile* fptr)
{
    // Return
    return;
}


/***********************************************************************//**
 * @brief Save ASCII table
 *
 * NOT YET IMPLEMENTED
 ***************************************************************************/
void GFitsAsciiTable::save(void)
{
    cout << "GFitsAsciiTable::save entry" << endl;
    // Return
    return;
}


/***********************************************************************//**
 * @brief Close Table
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


/***********************************************************************//**
 * @brief Return pointer to column with specified name
 *
 * @param colname Name of ASCII table column
 ***************************************************************************/
GFitsTableCol* GFitsAsciiTable::column(const std::string& colname)
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


/***********************************************************************//**
 * @brief Return pointer to column with number
 *
 * @param colnum Number of ASCII column (starting from 0)
 ***************************************************************************/
GFitsTableCol* GFitsAsciiTable::column(const int& colnum)
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

/***********************************************************************//**
 * @brief Initialise class members
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


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param table Table to be copied
 ***************************************************************************/
void GFitsAsciiTable::copy_members(const GFitsAsciiTable& table)
{
    // Copy attributes
    m_rows = table.m_rows;
    m_cols = table.m_cols;
    
    // Copy column definition
    if (table.m_columns != NULL && m_cols > 0) {
        m_columns = new GFitsTableCol*[m_cols];
        for (int i = 0; i < m_cols; ++i) {
            if (table.m_columns[i] != NULL)
                m_columns[i] = table.m_columns[i]->clone();
        }
    }
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
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


/***********************************************************************//**
 * @brief Connect ASCII table to FITS file
 *
 * @param fptr FITS file pointer to which the binary table should be connected
 *
 * The connection of the ASCII table is done by connecting all columns.
 ***************************************************************************/
void GFitsAsciiTable::connect(__fitsfile* fptr)
{
    // First connect ASCII table
    
    // Then connect all columns
    for (int i = 0; i < m_cols; ++i) {
        if (m_columns[i] != NULL) m_columns[i]->connect(fptr);
    }
    
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                          GFitsAsciiTable friends                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param os Stream into which the output is done
 * @param table Table to be dumped
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
