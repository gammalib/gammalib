/***************************************************************************
 *             GCsv.cpp - Column separated values table class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GCsv.cpp
 * @brief Column separated values table class implementation
 * @author Juergen Knoedlseder
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
#define G_ACCESS                               "GCsv::operator()(int&, int&)"
#define G_LOAD                       "GCsv::load(std::string&, std::string&)"
#define G_SAVE                "GCsv::save(std::string&, std::string&, bool&)"

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
 * @brief Table constructor
 *
 * @param[in] nrows Number of rows.
 * @param[in] ncols Number of columns.
 ***************************************************************************/
GCsv::GCsv(const int& nrows, const int& ncols)
{
    // Initialise private members
    init_members();

    // Allocate empty row
    std::vector<std::string> row(ncols, "");

    // Allocate table
    m_data.assign(nrows, row);

    // Set table dimensions
    m_cols = ncols;
    m_rows = nrows;

    // Return
    return;
}


/***********************************************************************//**
 * @brief File constructor
 *
 * @param[in] filename Filename.
 * @param[in] sep Column separator (default is whitespace).
 ***************************************************************************/
GCsv::GCsv(const std::string& filename, const std::string& sep)
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
 * @return Column separated values table.
 ***************************************************************************/
GCsv& GCsv::operator=(const GCsv& csv)
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
std::string& GCsv::operator()(const int& row, const int& col)
{
    // Perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_ACCESS, row, col, m_rows-1, m_cols-1);
    }
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
const std::string& GCsv::operator()(const int& row, const int& col) const
{
    // Perform range check
    #if defined(G_RANGE_CHECK)
    if (row < 0 || row >= m_rows || col < 0 || col >= m_cols) {
        throw GException::out_of_range(G_ACCESS, row, col, m_rows, m_cols);
    }
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
 * @brief Clear CSV table
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
 * @brief Clone CSV table
 *
 * @return Pointer to deep copy of CSV table.
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
    double value = gammalib::todouble((*this)(row, col));

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
    int value = gammalib::toint((*this)(row, col));

    // Return value
    return value;
}


/***********************************************************************//**
 * @brief Get string value
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 * @param[in] value String value.
 *
 * Set value of specified row and column as string.
 ***************************************************************************/
void GCsv::string(const int& row, const int& col, const std::string& value)
{
    // Set value
    (*this)(row, col) = value;
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get double precision value
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 * @param[in] value Double precision floating point value.
 *
 * Set value of specified row and column as double precision floating point.
 ***************************************************************************/
void GCsv::real(const int& row, const int& col, const double& value)
{
    // Set value
    (*this)(row, col) = gammalib::str(value, m_precision);
    
    // Return
    return;
}


/***********************************************************************//**
 * @brief Get integer value
 *
 * @param[in] row Table row.
 * @param[in] col Table column.
 * @param[in] value Integer value.
 *
 * Set value of specified row and column as integer.
 ***************************************************************************/
void GCsv::integer(const int& row, const int& col, const int& value)
{
    // Set value
    (*this)(row, col) = gammalib::str(value);
    
    // Return
    return;
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
 *
 * Load CSV table from ASCII file. Any environment variable present in the
 * filename will be expanded.
 **************************************************************************/
void GCsv::load(const std::string& filename, const std::string& sep)
{
    // Clear instance
    clear();

    // Allocate line buffer
    const int n = 10000; 
    char  line[n];

    // Expand environment variables
    std::string fname = gammalib::expand_env(filename);

    // Open CSV table (read-only)
    FILE* fptr = std::fopen(fname.c_str(), "r");
    if (fptr == NULL) {
        throw GException::file_not_found(G_LOAD, fname);
    }

    // Read lines
    int iline = 0;
    while (std::fgets(line, n, fptr) != NULL) {

        // Increment line counter
        iline++;

        // Get line with leading and trailing whitespace removed
        std::string sline =
            gammalib::strip_chars(gammalib::strip_whitespace(std::string(line)),"\n");

        // Skip line if empty
        if (sline.length() == 0) {
            continue;
        }

        // Split line in elements
        std::vector<std::string> elements = gammalib::split(sline, sep);
        for (int i = 0; i < elements.size(); ++i) {
            elements[i] = gammalib::strip_whitespace(elements[i]);
        }

        // If this is the first valid line then simply store the elements
        if (m_data.empty()) {
            m_data.push_back(elements);
            m_cols = elements.size();
        }

        // ... otherwise check table consistency and add elements
        else {
            // Check table consistency
            if (m_cols != elements.size()) {
                throw GException::csv_bad_columns(G_LOAD, fname,
                                  iline, m_cols, elements.size());
            }

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
 * @brief Save CSV table
 *
 * @param[in] filename Filename.
 * @param[in] sep Column separator (default: whitespace).
 * @param[in] bool Overwrite existing file? (default: false).
 *
 * @exception GException::invalid_value
 *            Attempt to overwrite existing file.
 * @exception GException::file_error
 *            Unable to create file.
 *
 * Save CSV table into ASCII file. Any environment variable present in the
 * filename will be expanded.
 **************************************************************************/
void GCsv::save(const std::string& filename, const std::string& sep,
                const bool& clobber) const
{
    // Throw exception if file exists but clobber flag is false
    if (!clobber && gammalib::file_exists(filename)) {
        std::string msg = "File \""+filename+"\" exists already but the "
                          "clobber flag is set to \"false\". Set the "
                          "clobber flag to true to overwrite the existing "
                          "file or specify another file name.";
        throw GException::invalid_value(G_SAVE, msg);
    }

    // Expand environment variables
    std::string fname = gammalib::expand_env(filename);

    // Open CSV table (write-only)
    FILE* fptr = std::fopen(fname.c_str(), "w");
    if (fptr == NULL) {
        std::string msg = "Unable to create file \""+filename+"\".";
        throw GException::file_error(G_SAVE, msg);
    }

    // Loop over the rows
    for (int row = 0; row < m_rows; ++row) {

        // Write columns
        for (int col = 0; col < m_cols; ++col) {
            std::fputs(m_data[row][col].c_str(), fptr);
            std::fputs(sep.c_str(), fptr);
        }

        // Write linebreak
        std::fputs("\n", fptr);

    } // endfor: looped over rows

    // Close file
    std::fclose(fptr);

    // Return
    return;
}


/***********************************************************************//**
 * @brief Print column separated values information
 *
 * @param[in] chatter Chattiness (defaults to NORMAL).
 * @return String containing column separated values information.
 ***************************************************************************/
std::string GCsv::print(const GChatter& chatter) const
{
    // Initialise result string
    std::string result;

    // Continue only if chatter is not silent
    if (chatter != SILENT) {

        // Append header
        result.append("=== GCsv ===");

        // Append information
        result.append("\n"+gammalib::parformat("Number of columns"));
        result.append(gammalib::str(m_cols));
        result.append("\n"+gammalib::parformat("Number of rows"));
        result.append(gammalib::str(m_rows));
        result.append("\n"+gammalib::parformat("Floating point precision"));
        if (m_precision == 0) {
            result.append("default");
        }
        else {
        result.append(gammalib::str(m_precision));
        }

    } // endif: chatter was not silent

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
    m_precision = 0;
  
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
    m_cols      = csv.m_cols;
    m_rows      = csv.m_rows;
    m_data      = csv.m_data;
    m_precision = csv.m_precision;
    
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
