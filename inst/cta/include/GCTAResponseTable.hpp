/***************************************************************************
 *             GCTAResponseTable.hpp - CTA response table class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2016 by Juergen Knoedlseder                         *
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
 * @author Juergen Knoedlseder
 */

#ifndef GCTARESPONSETABLE_HPP
#define GCTARESPONSETABLE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GBase.hpp"
#include "GNodeArray.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GCTAResponseTable
 *
 * @brief CTA response table class
 *
 * A CTA response table holds n-dimensional data cubes that describe a
 * component of the instrumental response function.
 *
 * Each dimension of the n-dimensional cube is describe by two axes vectors,
 * providing the lower and upper bin boundaries of each axis.
 *
 * Each n-dimensional data cube is called a table. Table elements can be
 * accessed by element index, or through linear, bilinear or trilinear
 * interpolation operators.
 ***************************************************************************/
class GCTAResponseTable : public GBase {

public:
    // Constructors and destructors
    GCTAResponseTable(void);
    GCTAResponseTable(const GCTAResponseTable& table);
    explicit GCTAResponseTable(const GFitsTable& hdu);
    virtual ~GCTAResponseTable(void);

    // Operators
    GCTAResponseTable&  operator=(const GCTAResponseTable& table);
    std::vector<double> operator()(const double& arg) const;
    std::vector<double> operator()(const double& arg1, const double& arg2) const;
    std::vector<double> operator()(const double& arg1, const double& arg2,
                                   const double& arg3) const;
    const double&       operator()(const int& element) const;
    double&             operator()(const int& element);
    const double&       operator()(const int& table, const int& element) const;
    double&             operator()(const int& table, const int& element);
    double              operator()(const int& table, const double& arg) const;
    double              operator()(const int& table, const double& arg1,
                                   const double& arg2) const;
    double              operator()(const int& index, const double& arg1,
                                   const double& arg2, const double& arg3) const;

    // Methods
    void               clear(void);
    GCTAResponseTable* clone(void) const;
    std::string        classname(void) const;
    int                tables(void) const;
    const int&         elements(void) const;
    const std::string& unit(const int& table) const;
    void               scale(const int& table, const double& scale);
    void               append_table(const std::string& name,
                                    const std::string& unit);
    const int&         axes(void) const;
    int                axis(const std::string& name) const;
    int                axis_bins(const int& axis) const;
    const double&      axis_lo(const int& axis, const int& bin) const;
    const double&      axis_hi(const int& axis, const int& bin) const;
    const GNodeArray&  axis_nodes(const int& axis) const;
    const std::string& axis_lo_name(const int& axis) const;
    const std::string& axis_hi_name(const int& axis) const;
    const std::string& axis_lo_unit(const int& axis) const;
    const std::string& axis_hi_unit(const int& axis) const;
    void               axis_linear(const int& axis);
    void               axis_log10(const int& axis);
    void               axis_radians(const int& axis);
    void               append_axis(const std::vector<double>& axis_lo, 
                                   const std::vector<double>& axis_hi,
                                   const std::string&         name,
                                   const std::string&         unit);    
    void               read(const GFitsTable& table);
    void               write(GFitsTable& table) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAResponseTable& table);
    void free_members(void);
    void read_colnames(const GFitsTable& hdu);
    void read_axes(const GFitsTable& hdu);
    void read_tables(const GFitsTable& hdu);
    void update(const double& arg) const;
    void update(const double& arg1, const double& arg2) const;
    void update(const double& arg1, const double& arg2,
                const double& arg3) const;

    // Table information
    int                               m_naxes;         //!< Number of axes
    int                               m_ntables;       //!< Number of tables
    int                               m_nelements;     //!< Number of elements per table
    std::vector<std::string>          m_colname_lo;    //!< Column names for lower boundaries
    std::vector<std::string>          m_colname_hi;    //!< Column names for upper boundaries
    std::vector<std::string>          m_colname_table; //!< Column names for table
    std::vector<std::vector<double> > m_axis_lo;       //!< Axes lower boundaries
    std::vector<std::vector<double> > m_axis_hi;       //!< Axes upper boundaries
    std::vector<std::string>          m_units_lo;      //!< Lower boundaries units
    std::vector<std::string>          m_units_hi;      //!< Upper boundaries units
    std::vector<std::string>          m_units_table;   //!< Parameter units
    std::vector<GNodeArray>           m_axis_nodes;    //!< Axes node arrays
    std::vector<std::vector<double> > m_tables;        //!< Tables

    // Response table computation cache for 1D access
    mutable int    m_inx_left;        //!< Index of left node
    mutable int    m_inx_right;       //!< Index of right node
    mutable double m_wgt_left;        //!< Weight of left node
    mutable double m_wgt_right;       //!< Weight of right node

    // Response table computation cache for 2D access
    mutable int    m_inx1;            //!< Index of upper left node
    mutable int    m_inx2;            //!< Index of lower left node
    mutable int    m_inx3;            //!< Index of upper right node
    mutable int    m_inx4;            //!< Index of lower right node
    mutable double m_wgt1;            //!< Weight of upper left node
    mutable double m_wgt2;            //!< Weight of lower left node
    mutable double m_wgt3;            //!< Weight of upper right node
    mutable double m_wgt4;            //!< Weight of lower right node

    // Response table computation cache for 3D access
    mutable int    m_inx5;            //!< Index of upper left node
    mutable int    m_inx6;            //!< Index of lower left node
    mutable int    m_inx7;            //!< Index of upper right node
    mutable int    m_inx8;            //!< Index of lower right node
    mutable double m_wgt5;            //!< Weight of upper left node
    mutable double m_wgt6;            //!< Weight of lower left node
    mutable double m_wgt7;            //!< Weight of upper right node
    mutable double m_wgt8;            //!< Weight of lower right node
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAResponseTable").
 ***************************************************************************/
inline
std::string GCTAResponseTable::classname(void) const
{
    return ("GCTAResponseTable");
}


/***********************************************************************//**
 * @brief Return number of tables
 *
 * @return Number of tables.
 *
 * Returns the number of tables.
 ***************************************************************************/
inline
int GCTAResponseTable::tables(void) const
{
    return (m_ntables);
}


/***********************************************************************//**
 * @brief Return number of elements per table
 *
 * @return Number of elements per table.
 *
 * Returns the number of elements per table.
 ***************************************************************************/
inline
const int& GCTAResponseTable::elements(void) const
{
    return (m_nelements);
}


/***********************************************************************//**
 * @brief Return number of axes of the tables
 *
 * @return Number of axes of tables.
 *
 * Returns the number of axes of tables.
 ***************************************************************************/
inline
const int& GCTAResponseTable::axes(void) const
{
    // Return number of axes
    return (m_naxes);
}

#endif /* GCTARESPONSETABLE_HPP */
