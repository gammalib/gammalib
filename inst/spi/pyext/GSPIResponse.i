/***************************************************************************
 *               GSPIResponse.i - INTEGRAL/SPI response class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2020 by Juergen Knoedlseder                              *
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
 * @file GSPIResponse.i
 * @brief INTEGRAL/SPI instrument response function class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GSPIResponse.hpp"
%}


/***********************************************************************//**
 * @class GSPIResponse
 *
 * @brief INTEGRAL/SPI instrument response function class
 ***************************************************************************/
class GSPIResponse : public GResponse {

public:
    // Constructors and destructors
    GSPIResponse(void);
    GSPIResponse(const GSPIResponse& rsp);
    explicit GSPIResponse(const GFilename& rspname);
    virtual ~GSPIResponse(void);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GSPIResponse* clone(void) const;
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
    void             rspname(const GFilename& rspname);
    const GFilename& rspname(void) const;
    bool             is_precomputed(void) const;
    const double&    energy_keV(void) const;
    const double&    dlogE(void) const;
    const double&    gamma(void) const;
    void             set(const GSPIObservation& obs,
                         const GEnergy&         energy = GEnergy());
    double           irf_value(const GSkyDir&      srcDir,
                               const GSPIEventBin& bin,
                               const int&          ireg) const;
    double           zenith(const int& ipt, const GSkyDir& dir) const;
    double           azimuth(const int& ipt, const GSkyDir& dir) const;
    void             read(const GFits& fits);
    void             write(GFits& fits) const;
    void             load(const GFilename& filename);
    void             save(const GFilename& filename,
                          const bool&      clobber = false) const;
};


/***********************************************************************//**
 * @brief GSPIResponse class extension
 ***************************************************************************/
%extend GSPIResponse {
    GSPIResponse copy() {
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
