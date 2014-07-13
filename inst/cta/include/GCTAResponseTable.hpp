/***************************************************************************
 *             GCTAResponseTable.hpp - CTA response table class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2014 by Juergen Knoedlseder                         *
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
 * @brief Interface for the CTA response table class
 *
 * A response table contains response parameters in multi-dimensional vector
 * column format. Each dimension is described by axes columns. 
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
    const double&       operator()(const int& index, const int& element) const;
    double&             operator()(const int& index, const int& element);
    double              operator()(const int& index, const double& arg) const;
    double              operator()(const int& index, const double& arg1,
                                   const double& arg2) const;
    double              operator()(const int& index, const double& arg1,
                                   const double& arg2, const double& arg3) const;

    // Methods
    void               clear(void);
    GCTAResponseTable* clone(void) const;
    int                size(void) const;
    const int&         elements(void) const;
    const int&         axes(void) const;
    int                axis(const int& index) const;
    double             axis_lo(const int& index, const int& bin) const;
    double             axis_hi(const int& index, const int& bin) const;
    void               axis_linear(const int& index);
    void               axis_log10(const int& index);
    void               axis_radians(const int& index);
    std::string        axis_lo_name(const int& index) const;
    std::string        axis_hi_name(const int& index) const;
    std::string        axis_lo_unit(const int& index) const;
    std::string        axis_hi_unit(const int& index) const;
    std::string        unit(const int& index) const;
    const GNodeArray&  nodes(const int& index) const;
    void               scale(const int& index, const double& scale);
    void               read(const GFitsTable& hdu);
    void               write(GFitsTable& hdu) const;
    std::string        print(const GChatter& chatter = NORMAL) const;
    void               add_axis(std::vector<double> axis_lo, 
				std::vector<double> axis_hi,
				std::string name_lo, std::string name_hi,
				std::string unit);    
    void               add_par(std::string name, std::string unit);

private:
    // Methods
    void init_members(void);
    void copy_members(const GCTAResponseTable& table);
    void free_members(void);
    void read_colnames(const GFitsTable& hdu);
    void read_axes(const GFitsTable& hdu);
    void read_pars(const GFitsTable& hdu);
    void update(const double& arg) const;
    void update(const double& arg1, const double& arg2) const;
    void update(const double& arg1, const double& arg2,
                const double& arg3) const;

    // Table information
    int                               m_naxes;       //!< Number of axes
    int                               m_npars;       //!< Number of parameters
    int                               m_nelements;   //!< Number of elements per parameter
    std::vector<std::string>          m_colname_lo;  //!< Column names for lower boundaries
    std::vector<std::string>          m_colname_hi;  //!< Column names for upper boundaries
    std::vector<std::string>          m_colname_par; //!< Column names for parameters
    std::vector<std::vector<double> > m_axis_lo;     //!< Axes lower boundaries
    std::vector<std::vector<double> > m_axis_hi;     //!< Axes upper boundaries
    std::vector<std::string>          m_units_lo;    //!< Lower boundaries units
    std::vector<std::string>          m_units_hi;    //!< Upper boundaries units
    std::vector<std::string>          m_units_par;   //!< Parameter units
    std::vector<GNodeArray>           m_axis_nodes;  //!< Axes node arrays
    std::vector<std::vector<double> > m_pars;        //!< Parameters

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
 * @brief Return number of parameters in response table
 *
 * @return Number of parameters in response table.
 *
 * Returns the number of parameters in response table.
 ***************************************************************************/
inline
int GCTAResponseTable::size(void) const
{
    return m_npars;
}


/***********************************************************************//**
 * @brief Return number of elements per parameter
 *
 * @return Number of elements per parameter.
 *
 * Returns the number of elements per parameter.
 ***************************************************************************/
inline
const int& GCTAResponseTable::elements(void) const
{
    return m_nelements;
}


/***********************************************************************//**
 * @brief Return number of axes in response table
 *
 * @return Number of axes in response table.
 *
 * Returns the number of axes in response table.
 ***************************************************************************/
inline
const int& GCTAResponseTable::axes(void) const
{
    // Return number of axes
    return m_naxes;
}

#endif /* GCTARESPONSETABLE_HPP */
