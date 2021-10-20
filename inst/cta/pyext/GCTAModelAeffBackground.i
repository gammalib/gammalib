/***************************************************************************
 *        GCTAModelAeffBackground.i - CTA Aeff background model class      *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2015-2021 by Michael Mayer                               *
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
 * @file GCTAModelAeffBackground.i
 * @brief CTA IRF background model class definition
 * @author Michael Mayer
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelAeffBackground.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelAeffBackground
 *
 * @brief CTA IRF background model class
 ***************************************************************************/
class GCTAModelAeffBackground : public GModelData {
public:
    // Constructors and destructors
    GCTAModelAeffBackground(void);
    explicit GCTAModelAeffBackground(const GXmlElement& xml);
    explicit GCTAModelAeffBackground(const GModelSpectral& spectral);
    GCTAModelAeffBackground(const GModelSpectral& spectral,
                            const GModelTemporal& temporal);
    GCTAModelAeffBackground(const GCTAModelAeffBackground& bgd);
    virtual ~GCTAModelAeffBackground(void);

    // Implemented pure virtual methods
    virtual void                     clear(void);
    virtual GCTAModelAeffBackground* clone(void) const;
    virtual std::string              classname(void) const;
    virtual std::string              type(void) const;
    virtual bool                     is_constant(void) const;
    virtual double                   eval(const GEvent&       event,
                                          const GObservation& obs,
                                          const bool&         gradients = false) const;
    virtual double                   npred(const GEnergy&       obsEng,
                                           const GTime&         obsTime,
                                           const GPolarization& obsPol,
                                           const GObservation&  obs) const;
    virtual GCTAEventList*           mc(const GObservation& obs, GRan& ran) const;
    virtual void                     read(const GXmlElement& xml);
    virtual void                     write(GXmlElement& xml) const;

    // Other methods
    GModelSpectral* spectral(void) const;
    GModelTemporal* temporal(void) const;
    void            spectral(const GModelSpectral* spectral);
    void            temporal(const GModelTemporal* temporal);
};

/***********************************************************************//**
* @brief GCTAModelAeffBackground class extension
***************************************************************************/
%extend GCTAModelAeffBackground {
    GCTAModelAeffBackground copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.spectral(), self.temporal(),
                 gammalib.GModelData.__getstate__(self))
        return state
    def __setstate__(self, state):
        if state[0] == None and state[1] == None:
            self.__init__()
        elif state[1] == None:
            self.__init__(state[0])
        else:
            self.__init__(state[0], state[1])
        gammalib.GModelData.__setstate__(self, state[2])
}
};
