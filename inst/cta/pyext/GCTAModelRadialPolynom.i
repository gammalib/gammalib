/***************************************************************************
 *          GCTAModelRadialPolynom.i - Radial Polynom model class          *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011-2018 by Juergen Knoedlseder                         *
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
 * @file GCTAModelRadialPolynom.i
 * @brief Radial Polynom model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialPolynom.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialPolynom
 *
 * @brief Radial Polynom CTA model class
 ***************************************************************************/
class GCTAModelRadialPolynom : public GCTAModelRadial {
public:
    // Constructors and destructors
    GCTAModelRadialPolynom(void);
    explicit GCTAModelRadialPolynom(const std::vector<double>& coeffs);
    explicit GCTAModelRadialPolynom(const GXmlElement& xml);
    GCTAModelRadialPolynom(const GCTAModelRadialPolynom& model);
    virtual ~GCTAModelRadialPolynom(void);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelRadialPolynom* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual double                  eval(const double& offset,
                                         const bool& gradients = false) const;
    virtual GCTAInstDir             mc(GRan& ran) const;
    virtual double                  mc_max_value(const GCTAObservation& obs) const;
    virtual double                  omega(void) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;

    // Other methods
    //double coeff(void) const;
    //void   coeff(const double& value);
};


/***********************************************************************//**
 * @brief GCTAModelRadialPolynom class extension
 ***************************************************************************/
%extend GCTAModelRadialPolynom {
    GCTAModelRadialPolynom copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        xml = gammalib.GXmlElement()
        self.write(xml)
        state = (xml,)
        return state
    def __setstate__(self, state):
        if state[0].elements('parameter') == 0:
            self.__init__()
        else:
            self.__init__(state[0])
}
};
