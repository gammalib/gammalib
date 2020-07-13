/***************************************************************************
 *                  GFitsTable.cpp - FITS table base class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GFitsTable.cpp
 * @brief FITS table abstract base class implementation
 * @author Juergen Knoedlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstring>
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include "GException.hpp"
#include "GTools.hpp"
#include "GFitsCfitsio.hpp"
#include "GFits.hpp"
#include "GFitsTable.hpp"
#include "GFitsTableBitCol.hpp"
#include "GFitsTableByteCol.hpp"
#include "GFitsTableBoolCol.hpp"
#include "GFitsTableStringCol.hpp"
#include "GFitsTableUShortCol.hpp"
#include "GFitsTableShortCol.hpp"
#include "GFitsTableULongCol.hpp"
#include "GFitsTableLongCol.hpp"
#include "GFitsTableLongLongCol.hpp"
#include "GFitsTableFloatCol.hpp"
#include "GFitsTableDoubleCol.hpp"
#include "GFitsTableCFloatCol.hpp"
#include "GFitsTableCDoubleCol.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS1                              "GFitsTable::operator[](int&)"
#define G_ACCESS2                      "GFitsTable::operator[](std::string&)"
#define G_SET1                        "GFitsTable::set(int&, GFitsTableCol&)"
#define G_SET2                "GFitsTable::set(std::string&, GFitsTableCol&)"
#define G_INSERT1                   "GFitsTable::insert(int, GFitsTableCol&)"
#define G_INSERT2          "GFitsTable::insert(std::string&, GFitsTableCol&)"
#define G_REMOVE1                                  "GFitsTable::remove(int&)"
#define G_REMOVE2                          "GFitsTable::remove(std::string&)"
#define G_INSERT_ROWS                   "GFitsTable::insert_rows(int&, int&)"
#define G_REMOVE_ROWS                   "GFitsTable::remove_rows(int&, int&)"
#define G_DATA_OPEN                            "GFitsTable::data_open(void*)"
#define G_DATA_SAVE                                 "GFitsTable::data_save()"
#define G_GET_TFORM                             "GFitsTable::get_tform(int&)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */
//#define G_DEBUG_SAVE                          //!< Debug data_save() method


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 *
 * Construct empty FITS table.
 ***************************************************************************/
GFitsTable::GFitsTable(void) : GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Table constructor
 *
 * @param[in] nrows Number of rows in table.
 *
 * Construct FITS table with @p nrows table rows.
 ***************************************************************************/
GFitsTable::GFitsTable(const int& nrows) : GFitsHDU()
{
    // Initialise class members for clean destruction
    init_members();

    // Store numnber of rows
    m_rows = nrows;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] table FITS table.
 ***************************************************************************/
GFitsTable::GFitsTable(const GFitsTable& table) : GFitsHDU(table)
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
GFitsTable::~GFitsTable(void)
{
    // Free members
    free_members();

    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                Operators                                =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Assignment operator
 *
 * @param[in] table FITS table.
 * @return FITS table.
 ***************************************************************************/
GFitsTable& GFitsTable::operator=(const GFitsTable& table)
{
    // Execute only if object is not identical
    if (this != &table) {

        // Copy base class members
        this->GFitsHDU::operator=(table);

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


/***********************************************************************//**
 * @brief Returns pointer to table column
 *
 * @param[in] colnum Column number [0,...,ncols()-1].
 * @return Table column pointer.
 *
 * @exception GException::fits_no_data
 *            No data found in table.
 * @exception GException::out_of_range
 *            Column number is out of range.
 ***************************************************************************/
GFitsTableCol* GFitsTable::operator[](const int& colnum)
{
    // If there is no data then throw an exception
    if (m_columns == NULL) {
        throw GException::fits_no_data(G_ACCESS1, "No columns in table.");
    }

    // Compile option: raise exception if column number is out of range
    #if defined(G_RANGE_CHECK)
    if (colnum < 0 || colnum >= m_cols) {
        throw GException::out_of_range(G_ACCESS1, colnum, 0, m_cols-1);
    }
    #endif

    // Get column pointer
    GFitsTableCol* ptr = m_columns[colnum];
    if (ptr == NULL) {
        throw GException::fits_no_data(G_ACCESS1, "No data for this column.");
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Returns pointer to table column (const version)
 *
 * @param[in] colnum Column number [0,...,ncols()-1].
 * @return Table column pointer.
 *
 * @exception GException::fits_no_data
 *            No data found in table.
 * @exception GException::out_of_range
 *            Column number is out of range.
 ***************************************************************************/
const GFitsTableCol* GFitsTable::operator[](const int& colnum) const
{
    // If there is no data then throw an exception
    if (m_columns == NULL) {
        throw GException::fits_no_data(G_ACCESS1, "No columns in table.");
    }

    // Compile option: raise exception if column number is out of range
    #if defined(G_RANGE_CHECK)
    if (colnum < 0 || colnum >= m_cols) {
        throw GException::out_of_range(G_ACCESS1, colnum, 0, m_cols-1);
    }
    #endif

    // Get column pointer
    const GFitsTableCol* ptr = m_columns[colnum];
    if (ptr == NULL) {
        throw GException::fits_no_data(G_ACCESS1, "No data for this column.");
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Returns pointer to table column
 *
 * @param[in] colname Column name.
 * @return Table column pointer.
 *
 * @exception GException::fits_no_data
 *            No data found in table.
 * @exception GException::fits_column_not_found
 *            Column name not found.
 ***************************************************************************/
GFitsTableCol* GFitsTable::operator[](const std::string& colname)
{
    // If there is no data then throw an exception
    if (m_columns == NULL) {
        throw GException::fits_no_data(G_ACCESS2, "No columns in table.");
    }

    // Get column pointer
    GFitsTableCol* ptr = ptr_column(colname);

    // If column has not been found throw an exception
    if (ptr == NULL) {
        throw GException::fits_column_not_found(G_ACCESS2, colname);
    }

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Returns pointer to table column (const version)
 *
 * @param[in] colname Column name.
 * @return Table column pointer.
 *
 * @exception GException::fits_no_data
 *            No data found in table.
 * @exception GException::fits_column_not_found
 *            Column name not found.
 ***************************************************************************/
const GFitsTableCol* GFitsTable::operator[](const std::string& colname) const
{
    // If there is no data then throw an exception
    if (m_columns == NULL) {
        throw GException::fits_no_data(G_ACCESS2, "No columns in table.");
    }

    // Get column pointer
    GFitsTableCol* ptr = ptr_column(colname);

    // If column has not been found throw an exception
    if (ptr == NULL) {
        throw GException::fits_column_not_found(G_ACCESS2, colname);
    }

    // Return pointer
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Set column
 *
 * @param[in] colnum Column number [0,...,ncols()-1].
 * @param[in] column Table column.
 * @return Pointer to table column that has been set.
 *
 * Sets the column of a table by making a deep copy of the @p column
 * provided.
 ***************************************************************************/
GFitsTableCol* GFitsTable::set(const int& colnum, const GFitsTableCol& column)
{
    // Check if column number is valid
    if (colnum < 0 || colnum >= ncols()) {
        throw GException::out_of_range(G_SET1, "Column number", colnum, ncols());
    }

    // Free existing column only if it differs from current column. This
    // prevents unintential deallocation of the argument
    if ((m_columns[colnum] != NULL) && (m_columns[colnum] != &column)) {
        delete m_columns[colnum];
    }

    // Clone column
    m_columns[colnum] = column.clone();

    // Update header
    update_header();

    // Return pointer to column
    return m_columns[colnum];
}


/***********************************************************************//**
 * @brief Set column
 *
 * @param[in] colname Column name.
 * @param[in] column Table column.
 * @return Pointer to table column that has been set.
 *
 * Sets the column of a table by making a deep copy of the @p column
 * provided.
 ***************************************************************************/
GFitsTableCol* GFitsTable::set(const std::string&   colname,
                               const GFitsTableCol& column)
{
    // Get column number
    int colnum = this->colnum(colname);

    // If column has not been found throw an exception
    if (colnum < 0) {
        throw GException::fits_column_not_found(G_SET2, colname);
    }

    // Set column and return pointer to column
    return (set(colnum, column));
}



/***********************************************************************//**
 * @brief Insert column into the table
 *
 * @param[in] colnum Column number [0,...,ncols()].
 * @param[in] column Table column.
 * @return Pointer to table column that has been appended.
 *
 * @exception GException::fits_bad_col_length
 *            The length of the column is incompatible with the number of
 *            rows in the table.
 *
 * A column will be inserted at position 'colnum' of the table. If the
 * position is beyond the end of the table the column will be appended.
 *
 * If the table is empty and has 0 rows, the number of rows will be set to
 * the length of the column.
 *
 * The length of the column to be inserted has to be identical to the number
 * of rows in the table.
 ***************************************************************************/
GFitsTableCol* GFitsTable::insert(int colnum, const GFitsTableCol& column)
{
    // Check if column number is valid
    if (colnum < 0 || colnum > ncols()) {
        throw GException::out_of_range(G_INSERT1, "Column number", colnum, ncols()+1);
    }

    // If the table is empty and has 0 rows then set the number of rows in
    // the table to the length of the column
    if (m_columns == NULL && m_rows == 0) {
        m_rows = column.nrows();
    }

    // Throw exception if the column length is incompatible with number of
    // rows in the table
    if (m_rows != column.nrows()) {
        throw GException::fits_bad_col_length(G_INSERT1,
                                              column.nrows(), m_rows);
    }

    // If no column data exist then allocate them now
    if (m_columns == NULL) {
        m_columns    = new GFitsTableCol*[1];
        m_cols       = 1;
        m_columns[0] = NULL;
        colnum       = 0;     // We insert at the first place
    }

    // ... otherwise make space to insert column
    else {

        // Allocate fresh memory
        GFitsTableCol** tmp = new GFitsTableCol*[m_cols+1];

        // Copy over old column pointers. Leave some space at the position
        // where we want to insert the column
        for (int src = 0, dst = 0; dst < m_cols+1; ++dst) {
            if (dst == colnum) {
                tmp[dst] = NULL;
            }
            else {
                tmp[dst] = m_columns[src];
                src++;
            }
        }

        // Free old column pointer array. Do not free the columns since they
        // are still alive in the new column pointer array.
        delete [] m_columns;

        // Connect new memory
        m_columns = tmp;

        // Increment column counter
        m_cols++;

    } // endelse: we made space to insert a column

    // Copy column by cloning it
    m_columns[colnum] = column.clone();

    // Reset column number since column does not already exist in FITS
    // file
    m_columns[colnum]->colnum(0);

    // Update header
    update_header();

    // Return column pointer
    return (m_columns[colnum]);
}


/***********************************************************************//**
 * @brief Insert column into the table
 *
 * @param[in] colname Column name.
 * @param[in] column Table column.
 * @return Pointer to table column that has been appended
 *
 * Inserts the column at the position given by the specified column name.
 ***************************************************************************/
GFitsTableCol* GFitsTable::insert(const std::string&   colname,
                                  const GFitsTableCol& column)
{
    // Get column number
    int colnum = this->colnum(colname);

    // If column has not been found throw an exception
    if (colnum < 0) {
        throw GException::fits_column_not_found(G_INSERT2, colname);
    }

    // Insert column and return pointer to column
    return (insert(colnum, column));
}


/***********************************************************************//**
 * @brief Remove column from the table
 *
 * @param[in] colnum Column number [0,...,ncols()-1].
 *
 * @exception GException::fits_bad_col_length
 *            The length of the column is incompatible with the number of
 *            rows in the table.
 *
 * Remove the column at position @p colnum from the table.
 ***************************************************************************/
void GFitsTable::remove(const int& colnum)
{
    // Check if column number is valid
    if (colnum < 0 || colnum >= ncols()) {
        throw GException::out_of_range(G_REMOVE1, "Column number", colnum, ncols());
    }

    // At this point we should have at least one column in the table and we
    // can now reduce the table by removing the column. If there is just a
    // single column in the table then delete all columns and set the number
    // of columns to zero ...
    if (m_cols == 1) {
        free_columns();
        m_cols = 0;
    }

    // ... otherwise remove the specified column
    else {

        // Allocate fresh memory
        GFitsTableCol** tmp = new GFitsTableCol*[m_cols-1];

        // Copy over old column pointers, skipping the column that should be
        // removed
        for (int src = 0, dst = 0; src < m_cols; ++src) {
            if (src != colnum) {
                tmp[dst] = m_columns[src];
                tmp[dst]->colnum(dst);
                dst++;
            }
        }

        // Free deleted column and old column pointer array. Do not free the
        // other columns since they are still alive in the new column pointer
        // array.
        if (m_columns[colnum] != NULL) {
            delete m_columns[colnum];
        }
        delete [] m_columns;

        // Connect new memory
        m_columns = tmp;

        // Decerement column counter
        m_cols--;

    } // endelse: we made space to insert a column

    // Update header
    update_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove column from the table
 *
 * @param[in] colname Column name.
 *
 * Remove the column with name @p colname from the table.
 ***************************************************************************/
void GFitsTable::remove(const std::string& colname)
{
    // Get column number
    int colnum = this->colnum(colname);

    // If column has not been found throw an exception
    if (colnum < 0) {
        throw GException::fits_column_not_found(G_REMOVE2, colname);
    }

    // Remove column
    remove(colnum);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Append rows to the table
 *
 * @param[in] nrows Number of rows to be appended.
 *
 * Appends @p nrows rows at the end of the FITS table. The method calls the
 * insert() method to do the job. Note that a call of this method will lead
 * to loading all table data into memory.
 ***************************************************************************/
void GFitsTable::append_rows(const int& nrows)
{
    // Set row number for insertion to end of the file
    int rownum = this->nrows();

    // Insert rows
    insert_rows(rownum, nrows);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Insert rows into the table
 *
 * @param[in] row Row at which rows are inserted [0,...,nrows()].
 * @param[in] nrows Number of rows to be inserted.
 *
 * @exception GException::fits_invalid_row
 *            Specified rownum is invalid.
 *
 * Inserts @p nrows table rows at the specified @p row into the table. If the
 * @p row index is set to the number of existing rows, nrows(), @p nrows
 * will be appened at the end of the table.
 *
 * Note that this method will load all table data into memory.
 ***************************************************************************/
void GFitsTable::insert_rows(const int& row, const int& nrows)
{
    // Make sure that row is valid
    if (row < 0 || row > m_rows) {
        throw GException::fits_invalid_row(G_INSERT_ROWS, row, m_rows);
    }

    // Continue only if there are rows to be inserted
    if (nrows > 0) {

        // Insert rows for all columns
        for (int icol = 0; icol < m_cols; ++icol) {
            m_columns[icol]->insert(row, nrows);
        }

        // Increment number of rows in table
        m_rows += nrows;

    } // endfor: there were rows to be inserted

    // Update header
    update_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Remove rows from the table
 *
 * @param[in] row Row from which on rows are removed [0,...,nrows()-1].
 * @param[in] nrows Number of rows to be removed.
 *
 * @exception GException::fits_invalid_row
 *            Specified rownum is invalid.
 * @exception GException::fits_invalid_nrows
 *            Invalid number of rows specified.
 *
 * Removes @p nrows table rows from the specified @p row on from the table.
 *
 * Note that this method will load all column data into memory.
 ***************************************************************************/
void GFitsTable::remove_rows(const int& row, const int& nrows)
{
    // Make sure that row is valid
    if (row < 0 || row >= m_rows) {
        throw GException::fits_invalid_row(G_REMOVE_ROWS, row, m_rows-1);
    }

    // Make sure that we don't remove beyond the limit
    if (nrows < 0 || nrows > m_rows-row) {
        throw GException::fits_invalid_nrows(G_REMOVE_ROWS, nrows, m_rows-row);
    }

    // Continue only if there are rows to be removed
    if (nrows > 0) {

        // Remove rows for all columns
        for (int icol = 0; icol < m_cols; ++icol) {
            m_columns[icol]->remove(row, nrows);
        }

        // Decrement number of rows in table
        m_rows -= nrows;

    } // endfor: there were rows to be removed

    // Update header
    update_header();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Checks the presence of a column in table
 *
 * @param[in] colname Column name.
 * @return True if the column exists, false otherwise.
 ***************************************************************************/
bool GFitsTable::contains(const std::string& colname) const
{
    // Get pointer in column
    GFitsTableCol* ptr = ptr_column(colname);

    // Return state
    return (ptr != NULL);
}


/***********************************************************************//**
 * @brief Print table information
 *
 * @param[in] chatter Chattiness.
 * @return String containing table information.
 ***************************************************************************/
std::string GFitsTable::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GFitsTable ===");

        // Append HDU information
        result.append("\n"+print_hdu(chatter));

        // Append table type
        result.append("\n"+gammalib::parformat("Table type"));
        switch (m_type) {
        case GFitsHDU::HT_ASCII_TABLE:
            result.append("ASCII table");
            break;
        case GFitsHDU::HT_BIN_TABLE:
            result.append("Binary table");
            break;
        default:
            result.append("Unknown");
            break;
        }

        // Append table dimensions
        result.append("\n"+gammalib::parformat("Number of rows"));
        result.append(gammalib::str(m_rows));
        result.append("\n"+gammalib::parformat("Number of columns"));
        result.append(gammalib::str(m_cols));

        // NORMAL: Append header information
        if (chatter >= NORMAL) {
            result.append("\n"+m_header.print(gammalib::reduce(chatter)));
        }

        // NORMAL: Append table columns
        if (chatter >= VERBOSE) {
            if (m_columns != NULL) {
                for (int i = 0; i < m_cols; ++i) {
                    result.append("\n");
                    if (m_columns[i] != NULL) {
                        result.append(m_columns[i]->print(gammalib::reduce(chatter)));
                    }
                    else {
                        result.append(" Column "+gammalib::str(i)+" undefined");
                    }
                }
            }
            else {
                result.append("\n Table columns undefined");
            }
        }

    } // endif: chatter was not silent

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Protected methods                           =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GFitsTable::init_members(void)
{
    // Initialise members
    m_type    = -1;
    m_rows    = 0;
    m_cols    = 0;
    m_columns = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] table FITS Table.
 *
 * Copies FITS table members from @p table.
 ***************************************************************************/
void GFitsTable::copy_members(const GFitsTable& table)
{
    // Copy attributes
    m_type = table.m_type;
    m_rows = table.m_rows;
    m_cols = table.m_cols;

    // Copy column definition
    if (table.m_columns != NULL && m_cols > 0) {
        m_columns = new GFitsTableCol*[m_cols];
        for (int i = 0; i < m_cols; ++i) {
            if (table.m_columns[i] != NULL) {
                m_columns[i] = table.m_columns[i]->clone();
            }
            else {
                m_columns[i] = NULL;
            }
        }
    }

    // Return
    return;
}


/***********************************************************************//**
 * @brief Free class members
 *
 * Free all class members.
 ***************************************************************************/
void GFitsTable::free_members(void)
{
    // Free columns
    free_columns();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Free column pointers
 *
 * De-allocates all column pointers.
 ***************************************************************************/
void GFitsTable::free_columns(void)
{
    // Free memory
    if (m_columns != NULL) {
        for (int i = 0; i < m_cols; ++i) {
            if (m_columns[i] != NULL) delete m_columns[i];
        }
        delete [] m_columns;
    }

    // Mark memory as freed
    m_columns = NULL;

    // Return
    return;
}


/***********************************************************************//**
 * @brief Update header after row or column manipulations
 *
 * Updates the header after row or column manipulations so that the header
 * is kept up-to-date with the current table structure.
 *
 * The following keywords are updated:
 *     - TFIELDS
 *     - NAXIS1
 *     - NAXIS2
 *     - TTYPEi
 *     - TFORMi
 *     - TUNITi
 *     - TZEROi (only for unsigned integer columns)
 *     - TSCALi (only for unsigned integer columns)
 *     - TDIMi (only for columns with dimensions)
 *
 * This method should be called after adding or removing table rows or
 * columns.
 ***************************************************************************/
void GFitsTable::update_header(void)
{
    // Get non-const reference to header for manipulations
    GFitsHeader& hdr = const_cast<GFitsHeader&>(header());

    // Initialise TFIELDS, NAXIS1 and NAXIS2 keywords
    int tfields = 0;
    int naxis1  = 0;
    int naxis2  = 0;

    // Compute NAXIS1 keyword value
    for (int i = 0; i < m_cols; ++i) {
        if (m_columns[i] != NULL) {
            naxis1 += m_columns[i]->number() * m_columns[i]->width();
        }
    }

    // Set TFIELDS, NAXIS1 and NAXIS2 keywords. In case that there are no
    // rows, columns or the table has a zero width, all keywords will be
    // set to zero for compliance with cfitsio.
    if ((m_cols > 0) && (m_rows > 0) && (naxis1 > 0)) {
        tfields = m_cols;
        naxis2  = m_rows;
    }
    else {
        tfields = 0;
        naxis1  = 0;
        naxis2  = 0;
    }

    // Update TFIELDS, NAXIS1 and NAXIS2 keywords
    card("TFIELDS", tfields, "number of table fields");
    card("NAXIS1",  naxis1,  "width of table in bytes");
    card("NAXIS2",  naxis2,  "number of rows in table");

    // Loop over all columns
    for (int i = 0; i < m_cols; ++i) {

        // Setup keywords
        std::string key_ttype = "TTYPE" + gammalib::str(i+1);
        std::string key_tform = "TFORM" + gammalib::str(i+1);
        std::string key_tunit = "TUNIT" + gammalib::str(i+1);
        std::string key_tzero = "TZERO" + gammalib::str(i+1);
        std::string key_tscal = "TSCAL" + gammalib::str(i+1);
        std::string key_tdim  = "TDIM"  + gammalib::str(i+1);

        // Continue only if the column is valid
        if (m_columns[i] != NULL) {

            // Get values
            char* ttype = get_ttype(i);
            char* tform = get_tform(i);
            char* tunit = get_tunit(i);

            // Convert into strings
            std::string val_ttype(ttype);
            std::string val_tform(tform);
            std::string val_tunit(tunit);

            // Replace "U" by "I" and "V" by "J" in TFORM keyword. This is
            // needed due to a cfitsio internal translation of these
            // characters.
            val_tform = gammalib::replace_segment(val_tform, "U", "I");
            val_tform = gammalib::replace_segment(val_tform, "V", "J");

            // Delete values
            if (ttype != NULL) delete [] ttype;
            if (tform != NULL) delete [] tform;
            if (tunit != NULL) delete [] tunit;

            // Build TDIM value
            std::vector<int> tdim     = m_columns[i]->dim();
            std::string      val_tdim = "";
            if (tdim.size() > 0) {
                val_tdim.append("("+gammalib::str(tdim[0]));
                for (int k = 1; k < tdim.size(); ++k) {
                    val_tdim.append(","+gammalib::str(tdim[k]));
                }
                val_tdim.append(")");
            }

            // Setup comments
            std::string com_ttype = "label for field   " + gammalib::str(i+1);
            std::string com_tform = "data format of field " + gammalib::str(i+1);
            std::string com_tunit = "physical unit of field" + gammalib::str(i+1);
            std::string com_tdim  = "dimensions of field " + gammalib::str(i+1);

            // Set keywords
            card(key_ttype, val_ttype, com_ttype);
            card(key_tform, val_tform, com_tform);
            card(key_tunit, val_tunit, com_tunit);
            if (!val_tdim.empty()) {
                card(key_tdim, val_tdim, com_tdim);
            }
            else {
                if (hdr.contains(key_tdim)) {
                    hdr.remove(key_tdim);
                }
            }

            // Set TZERO and TSCAL keywords for unsigned columns
            if (std::abs(m_columns[i]->type()) == __TUSHORT) {
                card(key_tzero, 0, "offset for unsigned integers");
                card(key_tscal, 0, "data are not scaled");
                hdr[key_tzero].value((unsigned long)(32768)); // Use appropriate value method
                hdr[key_tscal].value(1);
            }
            else if ((std::abs(m_columns[i]->type()) == __TULONG) ||
                     (std::abs(m_columns[i]->type()) == __TUINT)) {
                card(key_tzero, 0, "offset for unsigned integers");
                card(key_tscal, 0, "data are not scaled");
                hdr[key_tzero].value((unsigned long)(2147483648)); // Use appropriate value method
                hdr[key_tscal].value(1);
            }

        } // endif: column was valid

        // ... otherwise remove keywords
        else {
            if (hdr.contains(key_ttype)) hdr.remove(key_ttype);
            if (hdr.contains(key_tform)) hdr.remove(key_tform);
            if (hdr.contains(key_tunit)) hdr.remove(key_tunit);
            if (hdr.contains(key_tzero)) hdr.remove(key_tzero);
            if (hdr.contains(key_tscal)) hdr.remove(key_tscal);
            if (hdr.contains(key_tdim))  hdr.remove(key_tdim);
        }

    } // endif: loop over all columns

    // Loop over possible columns remaining in header
    for (int i = m_cols; i < 100; ++i) {

        // Setup keywords
        std::string key_ttype = "TTYPE" + gammalib::str(i+1);
        std::string key_tform = "TFORM" + gammalib::str(i+1);
        std::string key_tunit = "TUNIT" + gammalib::str(i+1);
        std::string key_tzero = "TZERO" + gammalib::str(i+1);
        std::string key_tscal = "TSCAL" + gammalib::str(i+1);
        std::string key_tdim  = "TDIM"  + gammalib::str(i+1);

        // Remove keywords
        if (hdr.contains(key_ttype)) hdr.remove(key_ttype);
        if (hdr.contains(key_tform)) hdr.remove(key_tform);
        if (hdr.contains(key_tunit)) hdr.remove(key_tunit);
        if (hdr.contains(key_tzero)) hdr.remove(key_tzero);
        if (hdr.contains(key_tscal)) hdr.remove(key_tscal);
        if (hdr.contains(key_tdim))  hdr.remove(key_tdim);

    } // endif: loop over possible columns remaining in header

    // Return
    return;
}


/***********************************************************************//**
 * @brief Open Table
 *
 * @param[in] vptr FITS file pointer.
 *
 * @exception GException::fits_hdu_not_found
 *            Specified HDU not found in FITS file.
 * @exception GException::fits_error
 *            A CFITSIO error occured during loading the table.
 * @exception GException::fits_unknown_coltype
 *            FITS column of unsupported type has been found in the FITS file.
 * @exception GException::fits_inconsistent_tdim
 *            The TDIM information provided in the header is inconsistent
 *            with the size of the column.
 *
 * This method loads a description of the table into memory. Column data are
 * not loaded at this point to save memory. They will be loaded on request
 * later. 
 ***************************************************************************/
void GFitsTable::data_open(void* vptr)
{
    // Move to HDU
    gammalib::fits_move_to_hdu(G_DATA_OPEN, vptr);

    // Save FITS file pointer
    FPTR_COPY(m_fitsfile, vptr);

    // Determine number of rows in table
    int status = 0;
    long nrows = 0;
    status     = __ffgnrw(FPTR(m_fitsfile), &nrows, &status);
    if (status != 0) {
        throw GException::fits_error(G_DATA_OPEN, status);
    }
    else {
        m_rows = (int)nrows;
    }

    // Determine number of columns in table
    status = __ffgncl(FPTR(m_fitsfile), &m_cols, &status);
    if (status != 0) {
        throw GException::fits_error(G_DATA_OPEN, status);
    }

    // Allocate and initialise memory for column pointers. Note that this
    // initialisation is needed to allow for a clean free_members() call
    // in case of any exception.
    free_columns();
    m_columns = new GFitsTableCol*[m_cols];
    for (int i = 0; i < m_cols; ++i) {
        m_columns[i] = NULL;
    }

    // Get table column information
    int  typecode = 0;
    long repeat   = 0;
    long width    = 0;
    for (int i = 0; i < m_cols; ++i) {

        // Get column name
        char keyname[10];
        char value[80];
        sprintf(keyname, "TTYPE%d", i+1);
        status = __ffgkey(FPTR(m_fitsfile), keyname, value, NULL, &status);
        if (status != 0) {
            throw GException::fits_error(G_DATA_OPEN, status);
        }
        value[strlen(value)-1] = '\0';

        // Get column definition
        status = __ffgtcl(FPTR(m_fitsfile), i+1, &typecode, &repeat, &width,
                          &status);
        if (status != 0) {
            throw GException::fits_error(G_DATA_OPEN, status);
        }

        // Check for unsigned columns
        unsigned long offset = 0;
        sprintf(keyname, "TZERO%d", i+1);
        status = __ffgky(FPTR(m_fitsfile), __TULONG, keyname, &offset, NULL, &status);
        if (status == 0) {
            if (typecode == __TSHORT && offset == 32768u) {
                typecode = __TUSHORT;
            }
            else if (typecode == -__TSHORT && offset == 32768u) {
                typecode = -__TUSHORT;
            }
            else if (typecode == __TLONG && offset == 2147483648u) {
                typecode = __TULONG;
            }
            else if (typecode == -__TLONG && offset == 2147483648u) {
                typecode = -__TULONG;
            }
            else if (typecode == __TINT && offset == 2147483648u) {
                typecode = __TUINT;
            }
            else if (typecode == -__TINT && offset == 2147483648u) {
                typecode = -__TUINT;
            }
            else {
                std::string msg = "But column '"+std::string(value)+"' has"
                                  " typecode "+gammalib::str(typecode)+" and"
                                  " unexpected associated TZERO="+
                                  gammalib::str(offset)+".";
                throw GException::fits_error(G_DATA_OPEN, 0, msg);
            }
        }
        else {
            status = 0;
        }

        // Get column unit (optional, leave blank if not found)
        char unit[80];
        sprintf(keyname, "TUNIT%d", i+1);
        status = __ffgkey(FPTR(m_fitsfile), keyname, unit, NULL, &status);
        if (status != 0) {
            status = 0;
            unit[0] = '\0';
            unit[1] = '\0';
        }
        else {
            unit[strlen(unit)-1] = '\0';
        }

        // Get column dimension (optional, leave blank if not found)
        char dim[80];
        sprintf(keyname, "TDIM%d", i+1);
        status = __ffgkey(FPTR(m_fitsfile), keyname, dim, NULL, &status);
        if (status != 0) {
            status = 0;
            dim[0] = '\0';
            dim[1] = '\0';
        }
        else {
            dim[strlen(dim)-1] = '\0';
        }

        // If found, extract column dimension into vector array of integers
        std::vector<int> vdim;
        std::string      sdim = gammalib::strip_chars(gammalib::strip_whitespace(&(dim[1])),"()");
        if (sdim.length() > 0) {
            std::vector<std::string> elements = gammalib::split(sdim, ",");
            for (int k = 0; k < elements.size(); ++k) {
                vdim.push_back(gammalib::toint(elements[k]));
            }
        }

        // Allocate column
        m_columns[i] = alloc_column(typecode);
        if (m_columns[i] == NULL) {
            std::ostringstream colname;
            colname << value;
            throw GException::fits_unknown_coltype(G_DATA_OPEN, colname.str(), typecode);
        }

        // Store column definition
        m_columns[i]->name(gammalib::strip_whitespace(&(value[1])));
        m_columns[i]->unit(gammalib::strip_whitespace(&(unit[1])));
        m_columns[i]->dim(vdim);
        m_columns[i]->colnum(i+1);
        m_columns[i]->type(typecode);
        m_columns[i]->repeat(repeat);
        m_columns[i]->width(width);
        m_columns[i]->nrows(m_rows);
        m_columns[i]->is_variable(typecode < 0);
        m_columns[i]->connect(FPTR(m_fitsfile));

        // Extract column vector size
        if (m_columns[i]->repeat() == 1) { // ASCII tables
            m_columns[i]->number(1);
        }
        else {                             // Binary tables
            if (typecode == __TSTRING) {
                m_columns[i]->number(m_columns[i]->repeat() /
                                     m_columns[i]->width());
            }
            else {
                m_columns[i]->number(m_columns[i]->repeat());
            }
        }

        // If TDIM information was set then check its consistency
        if (!vdim.empty()) {

            // Compute expectation
            int num = vdim[0];
            for (int k = 1; k < vdim.size(); ++k) {
                num *= vdim[k];
            }

            // Compare with real size
            if (num != m_columns[i]->number()) {
                throw GException::fits_inconsistent_tdim(G_DATA_OPEN,
                                                         vdim,
                                                         m_columns[i]->number());
            }

        } // endif: Valid TDIM information was found

    } // endfor: looped over all columns

    // Return
    return;
}


/***********************************************************************//**
 * @brief Save table into FITS file
 *
 * @exception GException::fits_error
 *            A CFITSIO error occured in this method.
 * @exception GException::fits_bad_col_length
 *            Table columns have inconsistent lengths.
 *
 * Saves the FITS table into the FITS file.
 *
 * In case that no table HDU exists so far in the FITS file, a new table
 * will be appended to the FITS file. In case that a different HDU type
 * exists at the requested extension number, the existing HDU will be
 * deleted and be replaced by the requested table type.
 *
 * The method also verifies the consistency of all table columns. Table
 * columns need to have identical lengths to be saved into a FITS table.
 * All columns with a length of zero will be excluded from saving, and if
 * they exist in the FITS file, they will be removed from the file.
 *
 * @todo This method should also update the header. Even easier, this method
 *       should save the header into the file using the m_header.save()
 *       method.
 *       Only this assures coherence between the files !!!! Once this has
 *       been implemented (also in the GFitsImage method) we should delete
 *       the m_header.save() call in GFitsHDU::save.
 ***************************************************************************/
void GFitsTable::data_save(void)
{
    // Debug definition: Dump method entry
    #if defined(G_DEBUG_SAVE)
    std::cout << "GFitsTable::save: entry" << std::endl;
    for (int i = 0; i < m_cols; ++i) {
        std::cout << m_columns[i]->print() << std::endl;
    }
    #endif

    // Make sure that column lengths are consistent with table length.
    // Columns with zero length will not be considered (why?)
    for (int i = 0; i < m_cols; ++i) {
        if (m_columns[i] != NULL && m_columns[i]->nrows() > 0) {
            if (m_columns[i]->nrows() != m_rows) {
                throw GException::fits_bad_col_length(G_DATA_SAVE,
                                                      m_columns[i]->nrows(),
                                                      m_rows);
            }
        }
    }

    // Move to HDU
    int status = 0;
    int type   = 0;
    status     = __ffmahd(FPTR(m_fitsfile), m_hdunum+1, &type, &status);

    // If move was successful but HDU type in file differs from HDU type
    // of object then replace the HDU in the file
    bool replace = (status == 0 && type != m_type);

    // If HDU does not exist in file or should be replaced then create or
    // replace it now
    if (status == 107 || replace) {

        // Reset status
        status = 0;

        // Initialise number of fields
        int tfields = 0;

        // Setup cfitsio column definition arrays
        char** ttype = NULL;
        char** tform = NULL;
        char** tunit = NULL;
        if (m_cols > 0) {
            ttype = new char*[m_cols];
            tform = new char*[m_cols];
            tunit = new char*[m_cols];
            for (int i = 0; i < m_cols; ++i) {
                ttype[i] = NULL;
                tform[i] = NULL;
                tunit[i] = NULL;
            }
            for (int i = 0; i < m_cols; ++i) {
                ttype[tfields] = get_ttype(i);
                tform[tfields] = get_tform(i);
                tunit[tfields] = get_tunit(i);
                if (ttype[tfields] != NULL && tform[tfields] != NULL && 
                    tunit[tfields] != NULL)
                    tfields++;
            }
        }

        // Replace FITS HDU by table
        if (replace) {

            // Delete current FITS HDU
            status = __ffdhdu(FPTR(m_fitsfile), NULL, &status);
            //if (status != 0) {
            //    throw GException::fits_error(G_DATA_SAVE, status);
            //}
            if (status == 0) {

                // Insert either ASCII or Binary table at current HDU position
                if (exttype() == GFitsHDU::HT_ASCII_TABLE) {
                    long tbcol  = 0;
                    long rowlen = 0;
                    status = __ffgabc(tfields, tform, 1, &rowlen, &tbcol, &status);
                    status = __ffitab(FPTR(m_fitsfile), rowlen, m_rows, tfields, ttype,
                                      &tbcol, tform, tunit, NULL, &status);
                }
                else {
                    status = __ffibin(FPTR(m_fitsfile), m_rows, tfields, ttype, tform,
                                      tunit, NULL, 0, &status);
                }
                //if (status != 0) {
                //    throw GException::fits_error(G_DATA_SAVE, status);
                //}

            } // endif: deletion of current FITS HDU successful

        } // endif: replacement of FITS table requested

        // ... otherwise create FITS table
        else {
            status = __ffcrtb(FPTR(m_fitsfile), m_type, m_rows, tfields,
                              ttype, tform, tunit, NULL, &status);
            //if (status != 0) {
            //    throw GException::fits_error(G_DATA_SAVE, status);
            //}
        }

        // De-allocate column definition arrays
        if (m_cols > 0) {
            for (int i = 0; i < m_cols; ++i) {
                if (ttype[i] != NULL) delete [] ttype[i];
                if (tform[i] != NULL) delete [] tform[i];
                if (tunit[i] != NULL) delete [] tunit[i];
            }
            if (ttype != NULL) delete [] ttype;
            if (ttype != NULL) delete [] tform;
            if (ttype != NULL) delete [] tunit;
        }

        // If status is not okay then throw an exception
        if (status != 0) {
            throw GException::fits_error(G_DATA_SAVE, status);
        }

        // Connect all existing columns to FITS table
        if (m_columns != NULL) {
            for (int i = 0; i < m_cols; ++i) {
                if (m_columns[i] != NULL) {
                    FPTR_COPY(m_columns[i]->m_fitsfile, m_fitsfile);
                    m_columns[i]->colnum(i+1);
                }
            }
        }

        // Debug option: Signal table creation
        #if defined(G_DEBUG_SAVE)
        std::cout << "GFitsTable::save: created new table" << std::endl;
        #endif

    }

    // ... otherwise we signal a FITS error
    else if (status != 0) {
        throw GException::fits_error(G_DATA_SAVE, status);
    }

    // Determine number of columns in table
    int num_cols = 0;
    status = __ffgncl(FPTR(m_fitsfile), &num_cols, &status);
    if (status != 0) {
        throw GException::fits_error(G_DATA_SAVE, status);
    }

    // Debug option: Log number of columns in FITS table
    #if defined(G_DEBUG_SAVE)
    std::cout << "GFitsTable::save: FITS table contains ";
    std::cout << num_cols << " columns." << std::endl;
    #endif

    // If we have no columns in the table then delete all columns from
    // FITS table
    if (m_columns == NULL && num_cols > 0) {
        for (int i = 0; i < num_cols; ++i) {
            status = __ffdcol(FPTR(m_fitsfile), i, &status);
            if (status != 0) {
                throw GException::fits_error(G_DATA_SAVE, status);
            }
        }
    }

    // ... otherwise update the FITS table
    else {

        // Determine number of rows in table
        long num_rows = 0;
        status        = __ffgnrw(FPTR(m_fitsfile), &num_rows, &status);
        if (status != 0) {
            throw GException::fits_error(G_DATA_SAVE, status);
        }

        // Debug option: Log number of rows in FITS table
        #if defined(G_DEBUG_SAVE)
        std::cout << "GFitsTable::save: FITS table contains ";
        std::cout << num_rows << " rows." << std::endl;
        #endif

        // If the table length differs from number of rows in the FITS file
        // then readjust FITS table length. We do this by adding or
        // deleting rows at the end of the table as we anyways rewrite the
        // entire table later
        if (m_rows > num_rows) {

            // Insert rows at the end of the table
            long long firstrow = num_rows;
            long long nrows    = m_rows - num_rows;
            status = __ffirow(FPTR(m_fitsfile), firstrow, nrows, &status);
            if (status != 0) {
                throw GException::fits_error(G_DATA_SAVE, status);
            }

            // Debug option: Log row adding
            #if defined(G_DEBUG_SAVE)
            std::cout << "GFitsTable::save: Added " << nrows;
            std::cout << " rows to FITS table." << std::endl;
            #endif

        }
        else if (m_rows < num_rows) {

            // Delete rows at the end of the table
            long long firstrow = num_rows;
            long long nrows    = num_rows - m_rows;
            status = __ffdrow(FPTR(m_fitsfile), firstrow, nrows, &status);
            if (status != 0) {
                throw GException::fits_error(G_DATA_SAVE, status);
            }

            // Debug option: Log row adding
            #if defined(G_DEBUG_SAVE)
            std::cout << "GFitsTable::save: Deleted " << nrows;
            std::cout << " rows from FITS table." << std::endl;
            #endif

        }

        // Debug option: Show where we are
        #if defined(G_DEBUG_SAVE)
        std::cout << "GFitsTable::save: Now update all columns." << std::endl;
        #endif

        // Update all columns. The 'm_colnum' field specifies where in the
        // FITS file the column resides. If 'm_colnum=0' then we have a new
        // column that does not yet exist. In this case we append a new column
        // to the FITS file.
        for (int i = 0; i < m_cols; ++i) {

            // Only consider valid columns
            if (m_columns[i] != NULL) {

                // If column has no correspondance than add new column in
                // FITS table and link column to table.
                if (m_columns[i]->colnum() == 0) {

                    // Increment number of columns in FITS file
                    num_cols++;

                    // Append column to FITS file
                    status = __fficol(FPTR(m_fitsfile), num_cols, get_ttype(i),
                                      get_tform(i), &status);
                    if (status != 0) {
                        throw GException::fits_error(G_DATA_SAVE, status);
                    }

                    // Connect all column to FITS table by copying over the
                    // FITS file pointer.
                    FPTR_COPY(m_columns[i]->m_fitsfile, m_fitsfile);
                    m_columns[i]->colnum(num_cols);

                } // endif: column appended to FITS file

                // Now write column into FITS file (only if length is positive)
                if (m_columns[i]->nrows() > 0) {
                    // Debug option: Show which column we're going to write
                    #if defined(G_DEBUG_SAVE)
                    std::cout << "GFitsTable::save: Write column " << i;
                    std::cout << "." << std::endl;
                    #endif
                    m_columns[i]->save();
                }

            } // endif: column was valid
        } // endfor: looped over all table columns

        // Debug option: Show where we are
        #if defined(G_DEBUG_SAVE)
        std::cout << "GFitsTable::save: Now delete all obsolete columns.";
        std::cout << std::endl;
        #endif

        // Delete all obsolete columns from FITS file. We do this from last to
        // first so that the column numbers remain valid. Also note that
        // FITS column counting starts from 1.
        for (int colnum = num_cols; colnum > 0; --colnum) {

            // Get column name from FITS file
            char keyname[10];
            char value[80];
            sprintf(keyname, "TTYPE%d", colnum);
            status = __ffgkey(FPTR(m_fitsfile), keyname, value, NULL, &status);
            if (status != 0) {
                throw GException::fits_error(G_DATA_SAVE, status);
            }
            value[strlen(value)-1] = '\0';
            std::string colname = gammalib::strip_whitespace(&(value[1]));

            // Check if this column is actually in our list of columns
            bool used = false;
            for (int i = 0; i < m_cols; ++i) {
                if (m_columns[i]          != NULL &&
                    m_columns[i]->nrows() > 0 &&
                    m_columns[i]->name()  == colname) {
                    used = true;
                    break;
                }
            }

            // If column is not used then delete it now from FITS table
            if (!used) {

                // Delete column
                status = __ffdcol(FPTR(m_fitsfile), colnum, &status);
                if (status != 0) {
                    throw GException::fits_error(G_DATA_SAVE, status);
                }

                // Debug option: Log column deletion
                #if defined(G_DEBUG_SAVE)
                std::cout << "GFitsTable::save: Deleted obsolete column ";
                std::cout << value << " from FITS table." << std::endl;
                #endif

            } // endif: deleted column

        } // endfor: Looped over all FITS columns

    } // endelse: FITS table has been updated

    // Debug option: Show where we are
    #if defined(G_DEBUG_SAVE)
    std::cout << "GFitsTable::save: Now update the FITS header for all columns.";
    std::cout << std::endl;
    #endif

    // Now update the header (IS THIS REALLY NEEDED HERE???)
    update_header();

    // Debug definition: Dump method exit
    #if defined(G_DEBUG_SAVE)
    std::cout << "GFitsTable::save: exit" << std::endl;
    #endif

    // Return
    return;
}


/***********************************************************************//**
 * @brief Close table
 ***************************************************************************/
void GFitsTable::data_close(void)
{
    // Free members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Connect table data to FITS file
 *
 * @param[in] vptr FITS file pointer.
 *
 * Connects the table columns to the file specified by the FITS file pointer.
 * This method does nothing if the file pointer in not valid.
 ***************************************************************************/
void GFitsTable::data_connect(void* vptr)
{
    // Continue only if file pointer is valid
    if (vptr != NULL) {

        // Connect all columns
        if (m_columns != NULL) {
            for (int i = 0; i < m_cols; ++i) {
                if (m_columns[i] != NULL) m_columns[i]->connect(vptr);
            }
        }

    } // endif: file pointer was valid

    // Return
    return;
}


/***********************************************************************//**
 * @brief Returns pointer to column type
 *
 * @param[in] colnum Column number for which type is to be returned
 *
 * This methods allocates memory for the character string that holds the
 * column type. The client has to de-allocate this memory after usage.
 * In case that the column does not exist a NULL pointer is returned.
 ***************************************************************************/
char* GFitsTable::get_ttype(const int& colnum) const
{
    // Initialise result with NULL pointer
    char* ptr = NULL;

    // Get type only if column exists
    if (m_columns != NULL && colnum >=0 && colnum < m_cols && 
        m_columns[colnum] != NULL) {
        int size = m_columns[colnum]->name().length();
        ptr      = new char[size+1];
        std::strncpy(ptr, m_columns[colnum]->name().c_str(), size);
        ptr[size] = '\0';
   }

    // Return result
    return ptr;
}


/***********************************************************************//**
 * @brief Returns pointer to column format
 *
 * @param[in] colnum Column number (starting from 0).
 *
 * @exception GException::fits_unknown_tabtype
 *            Table is neither ASCII nor Binary.
 *
 * This methods allocates memory for the character string that holds the
 * column format. The client has to de-allocate this memory after usage.
 * In case that the column does not exist a NULL pointer is returned.
 ***************************************************************************/
char* GFitsTable::get_tform(const int& colnum) const
{
    // Initialise result with NULL pointer
    char* ptr = NULL;

    // Get type only if column exists
    if (m_columns != NULL && colnum >=0 && colnum < m_cols &&
        m_columns[colnum] != NULL) {

        // Get table type specific format
        int size;
        switch (m_type) {
        case GFitsHDU::HT_ASCII_TABLE:
            size = m_columns[colnum]->ascii_format().length();
            ptr  = new char[size+1];
            std::strncpy(ptr, m_columns[colnum]->ascii_format().c_str(), size);
            break;
        case GFitsHDU::HT_BIN_TABLE:
            size = m_columns[colnum]->tform_binary().length();
            ptr  = new char[size+1];
            std::strncpy(ptr, m_columns[colnum]->tform_binary().c_str(), size);
            break;
        default:
            throw GException::fits_unknown_tabtype(G_GET_TFORM, m_type);
        }
        ptr[size] = '\0';
    }

    // Return result
    return ptr;
}


/***********************************************************************//**
 * @brief Returns pointer to column unit
 *
 * @param[in] colnum Column number (starting from 0).
 *
 * This methods allocates memory for the character string that holds the
 * column unit. The client has to de-allocate this memory after usage.
 * In case that the column does not exist a NULL pointer is returned.
 ***************************************************************************/
char* GFitsTable::get_tunit(const int& colnum) const
{
    // Initialise result with NULL pointer
    char* ptr = NULL;

    // Get type only if column exists
    if (m_columns != NULL && colnum >=0 && colnum < m_cols && 
        m_columns[colnum] != NULL) {
        int size = m_columns[colnum]->unit().length();
        ptr      = new char[size+1];
        std::strncpy(ptr, m_columns[colnum]->unit().c_str(), size);
        ptr[size] = '\0';
    }

    // Return result
    return ptr;
}


/*==========================================================================
 =                                                                         =
 =                              Private methods                            =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Allocates column
 *
 * @param[in] typecode cfitsio type code
 *
 * Allocates a table column depending on the cfitsio type code. If the type
 * code is not found then return a NULL pointer.
 ***************************************************************************/
GFitsTableCol* GFitsTable::alloc_column(int typecode) const
{
    // Initialise table column pointer
    GFitsTableCol* ptr = NULL;

    // Allocate column
    switch (std::abs(typecode)) {
    case __TBIT:
        ptr = new GFitsTableBitCol;
        break;
    case __TBYTE:
        ptr = new GFitsTableByteCol;
        break;
    case __TLOGICAL:
        ptr = new GFitsTableBoolCol;
        break;
    case __TSTRING:
        ptr = new GFitsTableStringCol;
        break;
    case __TUSHORT:
        ptr = new GFitsTableUShortCol;
        break;
    case __TSHORT:
        ptr = new GFitsTableShortCol;
        break;
    case __TULONG:
        ptr = new GFitsTableULongCol;
        break;
    case __TLONG:
        ptr = new GFitsTableLongCol;
        break;
    case __TFLOAT:
        ptr = new GFitsTableFloatCol;
        break;
    case __TLONGLONG:
        ptr = new GFitsTableLongLongCol;
        break;
    case __TDOUBLE:
        ptr = new GFitsTableDoubleCol;
        break;
    case __TCOMPLEX:
        ptr = new GFitsTableCFloatCol;
        break;
    case __TDBLCOMPLEX:
        ptr = new GFitsTableCDoubleCol;
        break;
    default:
        break;
    }

    // Return
    return ptr;
}


/***********************************************************************//**
 * @brief Returns pointer of column with given name
 *
 * @param[in] colname Name of column.
 * @return Pointer to table column (NULL if column has not been found).
 *
 * Returns a pointer to the column with the specified name. If more columns
 * with the same name exist, a pointer to the first of these columns is
 * returned. If no column with the specified name is found, a NULL pointer
 * will be returned.
 ***************************************************************************/
GFitsTableCol* GFitsTable::ptr_column(const std::string& colname) const
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

    // Return pointer
    return ptr;
}


/***********************************************************************//**
 * @brief Returns column number of a given column name
 *
 * @param[in] colname Name of column.
 * @return Column number (-1 if column has not been found).
 *
 * Returns the column number of the column with the specified name. If more
 * columns with the same name exist, the first of these columns is returned.
 ***************************************************************************/
int GFitsTable::colnum(const std::string& colname) const
{
    // Initialise column number
    int colnum = -1;

    // If there are columns then search for the specified name
    if (m_columns != NULL) {
        for (int i = 0; i < m_cols; ++i) {
            if (m_columns[i]->name() == colname) {
                colnum = i;
                break;
            }
        }
    }

    // Return column number
    return colnum;
}
