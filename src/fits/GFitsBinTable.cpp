/***************************************************************************
 *              GFitsBinTable.cpp  - FITS binary table class               *
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
#include "GFitsBinTable.hpp"
#include "GFitsTableStrCol.hpp"
#include "GFitsTableShtCol.hpp"
#include "GFitsTableLngCol.hpp"
#include "GFitsTableFltCol.hpp"
#include "GFitsTableDblCol.hpp"

/* __ Namespaces _________________________________________________________ */

/* __ Method name definitions ____________________________________________ */
#define G_OPEN    "GFitsBinTable::open(fitsfile*)"
#define G_COLUMN1 "GFitsBinTable::column(const std::string&)"
#define G_COLUMN2 "GFitsBinTable::column(const int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                   GFitsBinTable constructors/destructors                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Constructor
 ***************************************************************************/
GFitsBinTable::GFitsBinTable() : GFitsData()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Constructor
 *
 * @param[in] nrows Number of rows in table
 * @param[in] ncols Number of columns in table
 ***************************************************************************/
GFitsBinTable::GFitsBinTable(int nrows, int ncols) : GFitsData()
{
    // Initialise class members for clean destruction
    init_members();

    // Store table dimension
    m_rows = nrows;
    m_cols = ncols;

    // If there are columns then allocate memory for column pointers
    if (m_cols > 0) {

        // Allocate pointer memory
        m_columns = new GFitsTableCol*[m_cols];

        // Initialise pointers
        for (int i = 0; i < m_cols; ++i)
            m_columns[i] = NULL;

    } // endfor: looped over all columns

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] table Table which will be used to construct GFitsBinTable
 *                  instance
 ***************************************************************************/
GFitsBinTable::GFitsBinTable(const GFitsBinTable& table) : GFitsData(table)
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
GFitsBinTable::~GFitsBinTable()
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                         GFitsBinTable operators                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] table Table which will be assigned
 ***************************************************************************/
GFitsBinTable& GFitsBinTable::operator= (const GFitsBinTable& table)
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
 =                       GFitsBinTable public methods                      =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Open Table
 *
 * @param[in] fptr FITS file pointer
 *
 * Builds a description of the binary table in memory. Columns are not loaded
 * but column descriptors are allocated. The column data will only be loaded
 * once it needs to be accessed.
 ***************************************************************************/
void GFitsBinTable::open(__fitsfile* fptr)
{
    // Move to HDU
    int status = 0;
    status     = __ffmahd(fptr, (fptr->HDUposition)+1, NULL, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Determine number of rows in table
    long nrows  = 0;
    status      = __ffgnrw(fptr, &nrows, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);
    else
        m_rows = (int)nrows;

    // Determine number of columns in table
    status = __ffgncl(fptr, &m_cols, &status);
    if (status != 0)
        throw GException::fits_error(G_OPEN, status);

    // Allocate memory for column pointers
    if (m_columns != NULL) delete [] m_columns;
    m_columns = new GFitsTableCol*[m_cols];

    // Get table column information
    int  typecode = 0;
    long repeat   = 0;
    long width    = 0;
    for (int i = 0; i < m_cols; ++i) {

        // Get column name
        char keyname[10];
        char value[80];
        sprintf(keyname, "TTYPE%d", i+1);
        status = __ffgkey(fptr, keyname, value, NULL, &status);
        if (status != 0)
            throw GException::fits_error(G_OPEN, status);
        value[strlen(value)-1] = '\0';

        // Get column definition
        status = __ffgtcl(fptr, i+1, &typecode, &repeat, &width, &status);
        if (status != 0)
            throw GException::fits_error(G_OPEN, status);

        // Allocate column
        switch (typecode) {
        case __TSTRING:
            m_columns[i] = new GFitsTableStrCol();
            break;
        case __TSHORT:
            m_columns[i] = new GFitsTableShtCol();
            break;
        case __TLONG:
            m_columns[i] = new GFitsTableLngCol();
            break;
        case __TFLOAT:
            m_columns[i] = new GFitsTableFltCol();
            break;
        case __TDOUBLE:
            m_columns[i] = new GFitsTableDblCol();
            break;
        default:
            throw GException::fits_unknown_coltype(G_OPEN, typecode);
            break;
        }

        // Store column definition
        m_columns[i]->set_name(strip_whitespace(&(value[1])));
        m_columns[i]->set_colnum(i+1);
        m_columns[i]->set_type(typecode);
        m_columns[i]->set_repeat(repeat);
        m_columns[i]->set_width(width);
        m_columns[i]->set_length(m_rows);
        m_columns[i]->set_fitsfile(fptr);

    } // endfor: looped over all columns

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save binary table
 *
 * NOT YET IMPLEMENTED
 ***************************************************************************/
void GFitsBinTable::save(void)
{
    cout << "GFitsBinTable::save entry" << endl;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close binary table
 ***************************************************************************/
void GFitsBinTable::close(void)
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
 * @param[in] colname Name of column which should be accessed
 *
 * @exception GException::fits_no_data
 *            There are no column data in the table or there is not data
 *            for this column
 * @exception GException::fits_column_not_found
 *            Requested column has not been found in table.
 ***************************************************************************/
GFitsTableCol* GFitsBinTable::column(const std::string& colname)
{
    // If there is no data then throw an exception
    if (m_columns == NULL)
        throw GException::fits_no_data(G_COLUMN1, "No column data in table");

    // Initialise pointer
    GFitsTableCol* ptr = NULL;

    // If there are columns then search for the specified name
    if (m_columns != NULL) {
        for (int i = 0; i < m_cols; ++i) {
            if (m_columns[i]->name() == colname) {
                ptr = m_columns[i];
                if (ptr == NULL)
                    throw GException::fits_no_data(G_COLUMN1, 
                                                   "No data for this column");
                break;
            }
        }
    }

    // If column has not been found throw an exception
    if (ptr == NULL)
        throw GException::fits_column_not_found(G_COLUMN1, colname);

    // Return column pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Return pointer to column with number
 *
 * @param[in] colnum Number of requested column (starting from 0)
 *
 * @exception GException::fits_no_data
 *            There are no column data in the table or there is not data
 *            for this column
 * @exception GException::out_of_range
 *            Requested column has not been found in table.
 ***************************************************************************/
GFitsTableCol* GFitsBinTable::column(const int& colnum)
{
    // If there is no data then throw an exception
    if (m_columns == NULL)
        throw GException::fits_no_data(G_COLUMN2, "No column data in table");

    // If column number is out of range then throw an exception
    if (colnum < 0 || colnum >= m_cols)
        throw GException::out_of_range(G_COLUMN2, colnum, 0, m_cols-1);

    // Get column pointer
    GFitsTableCol* ptr = ptr = m_columns[colnum];
    if (ptr == NULL)
        throw GException::fits_no_data(G_COLUMN2, "No data for this column");

    // Return column pointer
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                       GFitsBinTable private methods                     =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsBinTable::init_members(void)
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
 * @param[in] table Table to copy
 ***************************************************************************/
void GFitsBinTable::copy_members(const GFitsBinTable& table)
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
void GFitsBinTable::free_members(void)
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
 * @brief Connect binary table to FITS file
 *
 * @param[in] fptr FITS file pointer to which the binary table should be 
 *                 connected
 *
 * The connection of the binary table is done by connecting all columns.
 ***************************************************************************/
void GFitsBinTable::connect(__fitsfile* fptr)
{
    // First connect binary table

    // Then connect all columns
    for (int i = 0; i < m_cols; ++i) {
        if (m_columns[i] != NULL) m_columns[i]->connect(fptr);
    }

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                           GFitsBinTable friends                         =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream
 * @param[in] table Table to put in output stream
 ***************************************************************************/
ostream& operator<< (ostream& os, const GFitsBinTable& table)
{
    // Put header in stream
    os << "=== GFitsBinTable ===" << endl;
    os << " Number of rows ............: " << table.m_rows << endl;
    os << " Number of columns .........: " << table.m_cols << endl;
    for (int i = 0; i < table.m_cols; ++i)
        os << " " << *(table.m_columns[i]);

    // Return output stream
    return os;
}


/*==========================================================================
 =                                                                         =
 =                   Other functions used by GFitsBinTable                 =
 =                                                                         =
 ==========================================================================*/
