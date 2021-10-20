/***************************************************************************
 *            GCTAModelBackground.i - Background model class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2018-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAModelBackground.i
 * @brief Background model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelBackground.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelBackground
 *
 * @brief Background model class
 ***************************************************************************/
class GCTAModelBackground : public GModelData {

public:
    // Constructors and destructors
    GCTAModelBackground(void);
    explicit GCTAModelBackground(const GXmlElement& xml);
    GCTAModelBackground(const GCTAModelSpatial& spatial,
                        const GModelSpectral&   spectral);
    GCTAModelBackground(const GCTAModelSpatial& spatial,
                        const GModelSpectral&   spectral,
                        const GModelTemporal&   temporal);
    GCTAModelBackground(const GCTAModelBackground& model);
    virtual ~GCTAModelBackground(void);

    // Implemented pure virtual methods
    virtual void                 clear(void);
    virtual GCTAModelBackground* clone(void) const;
    virtual std::string          classname(void) const;
    virtual std::string          type(void) const;
    virtual bool                 is_constant(void) const;
    virtual double               eval(const GEvent&       event,
                                      const GObservation& obs,
                                      const bool&         gradients = false) const;
    virtual double               npred(const GEnergy&       energy,
                                       const GTime&         time,
                                       const GPolarization& polarization,
                                       const GObservation&  obs) const;
    virtual GCTAEventList*       mc(const GObservation& obs, GRan& ran) const;
    virtual void                 read(const GXmlElement& xml);
    virtual void                 write(GXmlElement& xml) const;

    // Other methods
    GCTAModelSpatial* spatial(void) const;
    GModelSpectral*   spectral(void) const;
    GModelTemporal*   temporal(void) const;
    void              spatial(const GCTAModelSpatial* spatial);
    void              spectral(const GModelSpectral* spectral);
    void              temporal(const GModelTemporal* temporal);
};

/***********************************************************************//**
* @brief GCTAModelBackground class extension
***************************************************************************/
%extend GCTAModelBackground {
    GCTAModelBackground copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.spatial(), self.spectral(), self.temporal(),
                 gammalib.GModelData.__getstate__(self))
        return state
    def __setstate__(self, state):
        if state[0] != None and state[1] != None and state[2] != None:
            self.__init__(state[0], state[1], state[2])
        elif state[0] != None and state[1] != None:
            self.__init__(state[0], state[1])
        else:
            self.__init__()
        gammalib.GModelData.__setstate__(self, state[3])
}
};
