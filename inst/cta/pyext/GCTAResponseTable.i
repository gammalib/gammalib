/***************************************************************************
 *              GCTAResponseTable.i - CTA response table class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012-2020 by Juergen Knoedlseder                         *
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

    // Operators
    std::vector<double> operator()(const double& arg) const;
    std::vector<double> operator()(const double& arg1, const double& arg2) const;
    std::vector<double> operator()(const double& arg1, const double& arg2,
                                   const double& arg3) const;
    double              operator()(const int& index, const double& arg) const;
    double              operator()(const int& index, const double& arg1,
                                   const double& arg2) const;
    double              operator()(const int& index, const double& arg1,
                                   const double& arg2, const double& arg3) const;

    // Methods
    void               clear(void);
    GCTAResponseTable* clone(void) const;
    std::string        classname(void) const;
    bool               has_table(const std::string& name) const;
    bool               has_axis(const std::string& name) const;
    const int&         axes(void) const;
    const int&         tables(void) const;
    const int&         elements(void) const;
    int                axis(const std::string& name) const;
    int                table(const std::string& name) const;
    const std::string& unit(const int& table) const;
    void               scale(const int& table, const double& scale);
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
    const std::string& telescope(void) const;
    void               telescope(const std::string& telescope);
    void               append_axis(const std::vector<double>& axis_lo, 
                                   const std::vector<double>& axis_hi,
                                   const std::string&         name,
                                   const std::string&         unit);    
    void               append_table(const std::string& name,
                                    const std::string& unit);
    void               read(const GFitsTable& table);
    void               write(GFitsTable& table) const;
};


/***********************************************************************//**
 * @brief GCTAResponse class extension
 ***************************************************************************/
%extend GCTAResponseTable {
    double __getitem__(int GTuple1D2D[]) {
        if (GTuple1D2D[0] == 1) {
            return (*self)(GTuple1D2D[1]);
        }
        else {
            return (*self)(GTuple1D2D[1], GTuple1D2D[2]);
        }
    }
    void __setitem__(int GTuple1D2D[], double value) {
        if (GTuple1D2D[0] == 1) {
            (*self)(GTuple1D2D[1]) = value;
        }
        else {
            (*self)(GTuple1D2D[1], GTuple1D2D[2]) = value;
        }
    } 
    GCTAResponseTable copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        hdu = gammalib.GFitsBinTable()
        self.write(hdu)
        state = (hdu,)
        return state
    def __setstate__(self, state):
        self.__init__()
        self.read(state[0])
}
};
