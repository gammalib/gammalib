/***************************************************************************
 *              GCTAResponseTable.i - CTA response table class             *
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
 * @file GCTAResponseTable.i
 * @brief CTA response table class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GTools.hpp"
#include "GCTAResponseTable.hpp"
%}

/* Define std::vector<double> as valid return type (otherwise a memory leak
   occurs. */
%include "std_vector.i"
%template(VecDouble) std::vector<double>;


/***********************************************************************//**
 * @class GCTAResponseTable
 *
 * @brief Interface for the CTA response table class
 ***************************************************************************/
class GCTAResponseTable : public GBase {

public:
    // Constructors and destructors
    GCTAResponseTable(void);
    GCTAResponseTable(const GCTAResponseTable& table);
    explicit GCTAResponseTable(const GFitsTable& hdu);
    virtual ~GCTAResponseTable(void);

    // Interpolation operators
    std::vector<double> operator()(const double& arg) const;
    std::vector<double> operator()(const double& arg1, const double& arg2) const;
    std::vector<double> operator()(const double& arg1, const double& arg2,
                                   const double& arg3) const;
    double&             operator()(const int& element);
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
    std::string        axis_lo_unit(const int& index) const;
    std::string        axis_hi_unit(const int& index) const;
    std::string        unit(const int& index) const;
    const GNodeArray&  nodes(const int& index) const;
    void               scale(const int& index, const double& scale);
    void               read(const GFitsTable& hdu);
    void               write(GFitsTable& hdu) const;
};


/***********************************************************************//**
 * @brief GCTAResponse class extension
 ***************************************************************************/
%extend GCTAResponseTable {
    GCTAResponseTable copy() {
        return (*self);
    }
};
