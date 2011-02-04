/***************************************************************************
 *             GCsv.cpp - Column separated values table class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GCsv.cpp
 * @brief Column separated values table class implementation
 * @author J. Knodlseder
 */

/* __ Includes ___________________________________________________________ */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <cstdio>            // std::fopen, std::fgets, std::fclose, etc...
#include "GCsv.hpp"
#include "GTools.hpp"
#include "GException.hpp"

/* __ Method name definitions ____________________________________________ */
#define G_ACCESS                              "GCsv::operator() (int&, int&)"
#define G_LOAD                        "GCsv::load(std::string&, std::string)"

/* __ Macros _____________________________________________________________ */

/* __ Coding definitions _________________________________________________ */

/* __ Debug definitions __________________________________________________ */


/*==========================================================================
 =                                                                         =
 =                         Constructors/destructors                        =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Void constructor
 ***************************************************************************/
GCsv::GCsv(void)
{
    // Initialise private members
    init_members();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Filename.
 * @param[in] sep Column separator (default is whitespace).
 ***************************************************************************/
GCsv::GCsv(const std::string& filename, std::string sep)
{ 
    // Initialise private
    init_members();

    // Load CSV table
    load(filename, sep);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy constructor
 *
 * @param[in] csv Column separated values table.
 ***************************************************************************/
GCsv::GCsv(const GCsv& csv)
{ 
    // Initialise private
    init_members();

    // Copy members
    copy_members(csv);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Destructor
 ***************************************************************************/
GCsv::~GCsv(void)
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
 * @param[in] csv Column separated values table.
 ***************************************************************************/
GCsv& GCsv::operator= (const GCsv& csv)
{ 
    // Execute only if object is not identical
    if (this != &csv) {

        // Free members
        free_members();

        // Initialise private members for clean destruction
        init_members();

        // Copy members
        copy_members(csv);

    } // endif: object was not identical
  
    // Return
    return *this;
}


/***********************************************************************//**
 * @brief Table element access operator
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 ***************************************************************************/
std::string& GCsv::operator() (const int& row, const int& col)
{
    // Perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
        throw GException::out_of_range(G_ACCESS, row, col, m_rows-1, m_cols-1);
    #endif

    // Return element
    return m_data[row][col];
}


/***********************************************************************//**
 * @brief Table element access operator (const version)
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 ***************************************************************************/
const std::string& GCsv::operator() (const int& row, const int& col) const
{
    // Perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols)
        throw GException::out_of_range(G_ACCESS, row, col, m_rows, m_cols);
    #endif

    // Return element
    return m_data[row][col];
}


/*==========================================================================
 =                                                                         =
 =                             Public methods                              =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Clear instance
 *
 * This method properly resets the object to an initial state.
 ***************************************************************************/
void GCsv::clear(void)
{
    // Free class members
    free_members();

    // Initialise members
    init_members();

    // Return
    return;
}


/***********************************************************************//**
 * @brief Clone instance
 **************************************************************************/
GCsv* GCsv::clone(void) const
{
    return new GCsv(*this);
}


/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 *
 * Returns value of specified row and column as string.
 ***************************************************************************/
std::string GCsv::string(const int& row, const int& col) const
{
    // Return value
    return (*this)(row, col);
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 *
 * Returns value of specified row and column as double precision.
 ***************************************************************************/
double GCsv::real(const int& row, const int& col) const
{
    // Convert element into double
    double value = todouble((*this)(row, col));

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Get integer value
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 *
 * Returns value of specified row and column as integer.
 ***************************************************************************/
int GCsv::integer(const int& row, const int& col) const
{
    // Convert element into int
    int value = toint((*this)(row, col));

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Load CSV table
 *
 * @param[in] filename Filename.
 * @param[in] sep Column separator (default is whitespace).
 *
 * @exception GException::file_not_found
 *            CSV table file not found.
 * @exception GException::csv_bad_columns
 *            Inconsistent columns encountered in CSV table file.
 **************************************************************************/
void GCsv::load(const std::string& filename, std::string sep)
{
    // Clear instance
    clear();

    // If no seperator is given then assume a whitespace
    if (sep.length() == 0)
        sep = " ";

    // Allocate line buffer
    const int n = 10000; 
    char  line[n];

    // Open CSV table (read-only
    FILE* fptr = std::fopen(filename.c_str(), "r");
    if (fptr == NULL)
        throw GException::file_not_found(G_LOAD, filename);

    // Read lines
    int iline = 0;
    while (std::fgets(line, n, fptr) != NULL) {

        // Increment line counter
        iline++;

        // Get line with leading and trailing whitespace removed
        std::string sline = strip_chars(strip_whitespace(std::string(line)),"\n");

        // Skip line if empty
        if (sline.length() == 0)
            continue;

        // Split line in elements
        std::vector<std::string> elements = split(sline, sep);
        for (int i = 0; i < elements.size(); ++i)
            elements[i] = strip_whitespace(elements[i]);

        // If this is the first valid line then simply store the elements
        if (m_data.size() == 0) {
            m_data.push_back(elements);
            m_cols = elements.size();
        }

        // ... otherwise check table consistency and add elements
        else {
            // Check table consistency
            if (m_cols != elements.size())
                throw GException::csv_bad_columns(G_LOAD, filename,
                                  iline, m_cols, elements.size());

            // Append elements
            m_data.push_back(elements);
        }

        // Increment number of rows
        m_rows++;

    } // endwhile: looped over lines

    // Close file
    std::fclose(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print class information
 ***************************************************************************/
std::string GCsv::print(void) const
{
    // Initialise result string
    std::string result;

    // Append header
    result.append("=== GCsv ===");
    result.append("\n"+parformat("Number of columns")+str(m_cols));
    result.append("\n"+parformat("Number of rows")+str(m_rows));

    // Return result
    return result;
}


/*==========================================================================
 =                                                                         =
 =                             Private methods                             =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Initialise class members
 ***************************************************************************/
void GCsv::init_members(void)
{
    // Initialise members
    m_cols = 0;
    m_rows = 0;
    m_data.clear();
  
    // Return
    return;
}


/***********************************************************************//**
 * @brief Copy class members
 *
 * @param[in] csv Column separated values table.
 ***************************************************************************/
void GCsv::copy_members(const GCsv& csv)
{
    // Copy members
    m_cols = csv.m_cols;
    m_rows = csv.m_rows;
    m_data = csv.m_data;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Delete class members
 ***************************************************************************/
void GCsv::free_members(void)
{
    // Return
    return;
}


/*==========================================================================
 =                                                                         =
 =                                 Friends                                 =
 =                                                                         =
 ==========================================================================*/

/***********************************************************************//**
 * @brief Output operator
 *
 * @param[in] os Output stream.
 * @param[in] csv Column separated values table.
 ***************************************************************************/
std::ostream& operator<< (std::ostream& os, const GCsv& csv)
{
     // Write CSV in output stream
    os << csv.print();

    // Return output stream
    return os;
}


/***********************************************************************//**
 * @brief Log operator
 *
 * @param[in] log Logger.
 * @param[in] csv Column separated values table.
 ***************************************************************************/
GLog& operator<< (GLog& log, const GCsv& csv)
{
    // Write CSV into logger
    log << csv.print();

    // Return logger
    return log;
}
