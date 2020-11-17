/***************************************************************************
 *   GCTAModelSpatialGaussSpectrum.i - Spatial energy dependent Gaussian   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2019-2020 by Juergen Knoedlseder                         *
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
 * @file GCTAModelSpatialGaussSpectrum.i
 * @brief Spatial energy dependent Gaussian interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelSpatialGaussSpectrum.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelSpatialGaussSpectrum
 *
 * @brief Spatial energy dependent Gaussian model class
 ***************************************************************************/
class GCTAModelSpatialGaussSpectrum  : public GCTAModelSpatial {

public:
    // Constructors and destructors
    GCTAModelSpatialGaussSpectrum(void);
    explicit GCTAModelSpatialGaussSpectrum(const double& sigma);
    explicit GCTAModelSpatialGaussSpectrum(const GModelSpectral& sigma);
    explicit GCTAModelSpatialGaussSpectrum(const GXmlElement& xml);
    GCTAModelSpatialGaussSpectrum(const GCTAModelSpatialGaussSpectrum& model);
    virtual ~GCTAModelSpatialGaussSpectrum(void);

    // Implemented pure virtual methods
    virtual void                           clear(void);
    virtual GCTAModelSpatialGaussSpectrum* clone(void) const;
    virtual std::string                    classname(void) const;
    virtual std::string                    type(void) const;
    virtual double                         eval(const GCTAInstDir& dir,
                                                const GEnergy&     energy,
                                                const GTime&       time,
                                                const bool&        gradients = false) const;
    virtual double                         mc_max_value(const GCTAObservation& obs) const;
    virtual void                           read(const GXmlElement& xml);
    virtual void                           write(GXmlElement& xml) const;

    // Other methods
    const GModelSpectral* sigma(void) const;
    void                  sigma(const double& sigma);
    void                  sigma(const GModelSpectral& sigma);
};


/***********************************************************************//**
 * @brief GCTAModelSpatialGaussSpectrum class extension
 ***************************************************************************/
%extend GCTAModelSpatialGaussSpectrum {
    GCTAModelSpatialGaussSpectrum copy() {
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
