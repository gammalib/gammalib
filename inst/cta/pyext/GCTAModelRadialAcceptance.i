/***************************************************************************
 *       GCTAModelRadialAcceptance.i - Radial acceptance model class       *
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
 * @file GCTAModelRadialAcceptance.i
 * @brief Radial acceptance model class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCTAModelRadialAcceptance.hpp"
%}


/***********************************************************************//**
 * @class GCTAModelRadialAcceptance
 *
 * @brief Radial acceptance model class
 ***************************************************************************/
class GCTAModelRadialAcceptance : public GModelData {

public:
    // Constructors and destructors
    GCTAModelRadialAcceptance(void);
    explicit GCTAModelRadialAcceptance(const GXmlElement& xml);
    GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                              const GModelSpectral&  spectral);
    GCTAModelRadialAcceptance(const GCTAModelRadial& radial,
                              const GModelSpectral&  spectral,
                              const GModelTemporal&  temporal);
    GCTAModelRadialAcceptance(const GCTAModelRadialAcceptance& model);
    virtual ~GCTAModelRadialAcceptance(void);

    // Implemented pure virtual methods
    virtual void                       clear(void);
    virtual GCTAModelRadialAcceptance* clone(void) const;
    virtual std::string                classname(void) const;
    virtual std::string                type(void) const;
    virtual bool                       is_constant(void) const;
    virtual double                     eval(const GEvent& event,
                                            const GObservation& obs,
                                            const bool& gradients = false) const;
    virtual double                     npred(const GEnergy& obsEng, const GTime& obsTime,
                                             const GObservation& obs) const;
    virtual GCTAEventList*             mc(const GObservation& obs, GRan& ran) const;
    virtual void                       read(const GXmlElement& xml);
    virtual void                       write(GXmlElement& xml) const;

    // Other methods
    GCTAModelRadial* radial(void)   const;
    GModelSpectral*  spectral(void) const;
    GModelTemporal*  temporal(void) const;
    void             radial(const GCTAModelRadial* radial);
    void             spectral(const GModelSpectral* spectral);
    void             temporal(const GModelTemporal* temporal);
};


/***********************************************************************//**
 * @brief GCTAModelRadialAcceptance class extension
 ***************************************************************************/
%extend GCTAModelRadialAcceptance {
    GCTAModelRadialAcceptance copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.radial(), self.spectral(), self.temporal(),
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
