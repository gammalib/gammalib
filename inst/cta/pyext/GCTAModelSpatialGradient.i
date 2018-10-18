/***************************************************************************
 *       GCTAModelSpatialGradient.i - Spatial gradient CTA model class     *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018 by Juergen Knoedlseder                              *
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
 * @file GCTAModelSpatialGradient.i
 * @brief Spatial gradient CTA interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelSpatialGradient.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelSpatialGradient
 *
 * @brief Spatial gradient CTA model class
 ***************************************************************************/
class GCTAModelSpatialGradient  : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialGradient(void);
    explicit GCTAModelSpatialGradient(const GXmlElement& xml);
    GCTAModelSpatialGradient(const GCTAModelSpatialGradient& model);
    virtual ~GCTAModelSpatialGradient(void);

    // Operators
    virtual GCTAModelSpatialGradient& operator=(const GCTAModelSpatialGradient& model);

    // Implemented pure virtual methods
    virtual void                      clear(void);
    virtual GCTAModelSpatialGradient* clone(void) const;
    virtual std::string               classname(void) const;
    virtual std::string               type(void) const;
    virtual double                    eval(const GCTAInstDir& dir,
                                           const GEnergy&     energy,
                                           const GTime&       time,
                                           const bool&        gradients = false) const;
    virtual GCTAInstDir               mc(const GEnergy& energy,
                                         const GTime&   time,
                                         GRan& ran) const;
    virtual void                      read(const GXmlElement& xml);
    virtual void                      write(GXmlElement& xml) const;

    // Other methods
    double detx_gradient(void) const;
    double dety_gradient(void) const;
    void   detx_gradient(const double& detx_gradient);
    void   dety_gradient(const double& dety_gradient);
};


/***********************************************************************//**
 * @brief GCTAModelSpatialGradient class extension
 ***************************************************************************/
%extend GCTAModelSpatialGradient {
    GCTAModelSpatialGradient copy() {
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
