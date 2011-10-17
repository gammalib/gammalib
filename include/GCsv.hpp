/***************************************************************************
 *             GCsv.hpp - Column separated values table class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
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
 * @author J. Knodlseder
 */

#ifndef GCsv_HPP
#define GCsv_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GLog.hpp"


/***********************************************************************//**
 * @class GCsv
 *
 * @brief Column separated values table class definition
 ***************************************************************************/
class GCsv {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCsv& csv);
    friend GLog&         operator<< (GLog& log, const GCsv& csv);

public:
    // Constructors and destructors
    GCsv(void);
    GCsv(const std::string& filename, std::string sep = " ");
    GCsv(const GCsv& csv);
    virtual ~GCsv(void);
 
    // Operators
    GCsv&              operator= (const GCsv& csv);
    std::string&       operator() (const int& row, const int& col);
    const std::string& operator() (const int& row, const int& col) const;

    // Methods
    void        clear(void);
    GCsv*       clone(void) const;
    std::string string(const int& row, const int& col) const;
    double      real(const int& row, const int& col) const;
    int         integer(const int& row, const int& col) const;
    void        load(const std::string& filename, std::string sep = " ");
    int         ncols(void) const { return m_cols; }
    int         nrows(void) const { return m_rows; }
    int         size(void) const { return m_rows*m_cols; }
    std::string print(void) const;
  
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

#endif /* GCsv_HPP */
