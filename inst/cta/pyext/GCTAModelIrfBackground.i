/***************************************************************************
 *        GCTAModelIrfBackground.i - CTA IRF background model class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2014-2018 by Juergen Knoedlseder                         *
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
 * @file GCTAModelIrfBackground.i
 * @brief CTA IRF background model class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelIrfBackground.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelIrfBackground
 *
 * @brief CTA IRF background model class
 ***************************************************************************/
class GCTAModelIrfBackground : public GModelData {
public:
    // Constructors and destructors
    GCTAModelIrfBackground(void);
    explicit GCTAModelIrfBackground(const GXmlElement& xml);
    explicit GCTAModelIrfBackground(const GModelSpectral& spectral);
    GCTAModelIrfBackground(const GCTAModelIrfBackground& bgd);
    virtual ~GCTAModelIrfBackground(void);

    // Implemented pure virtual methods
    virtual void                    clear(void);
    virtual GCTAModelIrfBackground* clone(void) const;
    virtual std::string             classname(void) const;
    virtual std::string             type(void) const;
    virtual bool                    is_constant(void) const;
    virtual double                  eval(const GEvent& event,
                                         const GObservation& obs,
                                         const bool& gradients = false) const;
    virtual double                  npred(const GEnergy& obsEng,
                                          const GTime& obsTime,
                                          const GObservation& obs) const;
    virtual GCTAEventList*          mc(const GObservation& obs, GRan& ran) const;
    virtual void                    read(const GXmlElement& xml);
    virtual void                    write(GXmlElement& xml) const;

    // Other methods
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;
    void            spectral(const GModelSpectral* spectral);
    void            temporal(const GModelTemporal* temporal);
};

/***********************************************************************//**
* @brief GCTAModelIrfBackground class extension
***************************************************************************/
%extend GCTAModelIrfBackground {
    GCTAModelIrfBackground copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        args = gammalib.GModelData.__getstate__(self), self.spectral(), \
               self.temporal()
        return args
    def __setstate__(self, state):
        gammalib.GModelData.__setstate__(self, state[0])
        self.spectral(state[1])
        self.temporal(state[2])
}
};
