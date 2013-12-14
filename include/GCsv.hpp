/***************************************************************************
 *             GCsv.hpp - Column separated values table class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
    int         size(void) const;
    const int&  ncols(void) const;
    const int&  nrows(void) const;
    std::string string(const int& row, const int& col) const;
    double      real(const int& row, const int& col) const;
    int         integer(const int& row, const int& col) const;
    void        load(const std::string& filename, const std::string& sep = " ");
    std::string print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCsv& csv);
    void free_members(void);

    // Protected data members
    int                                    m_cols;  //!< Number of columns
    int                                    m_rows;  //!< Number of rows
    std::vector<std::vector<std::string> > m_data;  //!< CSV table data
};


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

#endif /* GCSV_HPP */
