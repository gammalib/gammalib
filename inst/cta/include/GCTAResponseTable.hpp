/***************************************************************************
 *            GCTAResponseTable.hpp  -  CTA response table class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCTAResponseTable.hpp
 * @brief CTA response table class definition
 * @author J. Knoedlseder
 */

#ifndef GCTARESPONSETABLE_HPP
#define GCTARESPONSETABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include <iostream>
#include "GLog.hpp"
#include "GNodeArray.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GCTAResponseTable
 *
 * @brief Interface for the CTA response table class
 *
 * A response table contains response parameters in multi-dimensional vector
 * column format. Each dimension is described by axes columns. 
 ***************************************************************************/
class GCTAResponseTable {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAResponseTable& table);
    friend GLog&         operator<< (GLog& log,        const GCTAResponseTable& table);

public:
    // Constructors and destructors
    GCTAResponseTable(void);
    GCTAResponseTable(const GCTAResponseTable& table);
    GCTAResponseTable(const GFitsTable* hdu);
    virtual ~GCTAResponseTable(void);

    // Operators
    GCTAResponseTable&  operator= (const GCTAResponseTable & table);
    std::vector<double> operator()(const double& arg) const;
    std::vector<double> operator()(const double& arg1, const double& arg2) const;

    // Methods
    void               clear(void);
    GCTAResponseTable* clone(void) const;
    int                size(void) const { return m_colname_lo.size(); }
    int                axis(const int& index) const;
    double             axis_lo(const int& index, const int& bin) const;
    double             axis_hi(const int& index, const int& bin) const;
    void               axis_linear(const int& index);
    void               axis_log10(const int& index);
    void               read(const GFitsTable* hdu);
    void               write(GFitsTable* hdu) const;
    std::string        print(void) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAResponseTable& table);
    void free_members(void);
    void read_colnames(const GFitsTable* hdu);
    void read_axes(const GFitsTable* hdu);
    void read_pars(const GFitsTable* hdu);

    // Table information
    std::vector<std::string>          m_colname_lo;  //!< Column names for lower boundaries
    std::vector<std::string>          m_colname_hi;  //!< Column names for upper boundaries
    std::vector<std::string>          m_colname_par; //!< Column names for parameters
    std::vector<std::vector<double> > m_axis_lo;     //!< Axes lower boundaries
    std::vector<std::vector<double> > m_axis_hi;     //!< Axes upper boundaries
    std::vector<GNodeArray>           m_axis_nodes;  //!< Axes node arrays
    std::vector<std::vector<double> > m_pars;        //!< Parameters
};

#endif /* GCTARESPONSETABLE_HPP */
