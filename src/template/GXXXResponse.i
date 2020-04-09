/***************************************************************************
 *                GXXXResponse.i - [INSTRUMENT] response class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXResponse.i
 * @brief [INSTRUMENT] instrument response function class definition
 * @author [AUTHOR]
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXResponse.hpp"
%}


/***********************************************************************//**
 * @class GXXXResponse
 *
 * @brief [INSTRUMENT] instrument response function class
 ***************************************************************************/
class GXXXResponse : public GResponse {

public:
    // Constructors and destructors
    GXXXResponse(void);
    GXXXResponse(const GXXXResponse& rsp);
    virtual ~GXXXResponse(void);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GXXXResponse* clone(void) const;
    virtual std::string   classname(void) const;
    virtual bool          use_edisp(void) const;
    virtual bool          use_tdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        nroi(const GModelSky&    model,
                               const GEnergy&      obsEng,
                               const GTime&        obsTime,
                               const GObservation& obs) const;
    virtual GEbounds      ebounds(const GEnergy& obsEnergy) const;

    // Other Methods
    // TODO: Copy methods from GXXXObservation.hpp file
};


/***********************************************************************//**
 * @brief GXXXResponse class extension
 ***************************************************************************/
%extend GXXXResponse {
    GXXXResponse copy() {
        return (*self);
    }
%pythoncode {
    def __getstate__(self):
        state = (self.classname()) # TODO: Replace by appropriate class members
        return state
    def __setstate__(self, state):
        self.__init__()
}
};
