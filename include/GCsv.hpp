/***************************************************************************
 *             GCsv.hpp - Column separated values table class              *
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
 * @file GCsv.hpp
 * @brief Column separated values table class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCSV_HPP
#define GCSV_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GCsv
 *
 * @brief Column separated values table class
 *
 * This class implements a table of std::string elements that is loaded
 * from a column separated value ASCII file. The column separation string
 * can be specified upon loading of the file (by default the class assumes
 * that elements are separated by a white space).
 *
 * The class provides operators for string element access, and methods for
 * conversion of the string values:
 *
 *    double      real    = csv.real(row,col);
 *    int         integer = csv.integer(row,col);
 *    std::string string  = csv.string(row,col);
 *  
 ***************************************************************************/
class GCsv : public GBase {

public:
    // Constructors and destructors
    GCsv(void);
    GCsv(const int& nrows, const int& ncols);
    GCsv(const std::string& filename, const std::string& sep = " ");
    GCsv(const GCsv& csv);
    virtual ~GCsv(void);
 
    // Operators
    GCsv&              operator=(const GCsv& csv);
    std::string&       operator()(const int& row, const int& col);
    const std::string& operator()(const int& row, const int& col) const;

    // Methods
    void        clear(void);
    GCsv*       clone(void) const;
    std::string classname(void) const;
    int         size(void) const;
    const int&  ncols(void) const;
    const int&  nrows(void) const;
    const int&  precision(void) const;
    void        precision(const int& precision);
    std::string string(const int& row, const int& col) const;
    double      real(const int& row, const int& col) const;
    int         integer(const int& row, const int& col) const;
    void        string(const int& row, const int& col, const std::string& value);
    void        real(const int& row, const int& col, const double& value);
    void        integer(const int& row, const int& col, const int& value);
    void        load(const std::string& filename, const std::string& sep = " ");
    void        save(const std::string& filename, const std::string& sep = " ",
                     const bool& clobber = false) const;
    std::string print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCsv& csv);
    void free_members(void);

    // Protected data members
    int                                    m_cols;      //!< Number of columns
    int                                    m_rows;      //!< Number of rows
    std::vector<std::vector<std::string> > m_data;      //!< CSV table data
    int                                    m_precision; //!< Precision for floats
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCsv").
 ***************************************************************************/
inline
std::string GCsv::classname(void) const
{
    return ("GCsv");
}


/***********************************************************************//**
 * @brief Return table size (columns times rows)
 *
 * @return Table size.
 ***************************************************************************/
inline
int GCsv::size(void) const
{
    return m_rows*m_cols;
}


/***********************************************************************//**
 * @brief Return number of columns
 *
 * @return Number of columns.
 ***************************************************************************/
inline
const int& GCsv::ncols(void) const
{
    return m_cols;
}


/***********************************************************************//**
 * @brief Return number of rows
 *
 * @return Number of rows.
 ***************************************************************************/
inline
const int& GCsv::nrows(void) const
{
    return m_rows;
}


/***********************************************************************//**
 * @brief Return fixed field floating point precision
 *
 * @return Fixed field floating point precision.
 *
 * Returns the precision for floating point values when setting values using
 * the real() method. Any value >0 indicates the number of decimal places
 * that the floating point value will have.
 ***************************************************************************/
inline
const int& GCsv::precision(void) const
{
    return m_precision;
}


/***********************************************************************//**
 * @brief Set fixed field floating point precision
 *
 * @param[in] precision Fixed field floating point precision.
 *
 * Set the precision for floating point values when setting values using
 * the real() method. Any value >0 indicates the number of decimal places
 * that the floating point value will have.
 ***************************************************************************/
inline
void GCsv::precision(const int& precision)
{
    m_precision = precision;
    return;
}

#endif /* GCSV_HPP */
