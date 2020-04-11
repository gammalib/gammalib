/***************************************************************************
 *                 GLATResponse.i - Fermi/LAT Response class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2020 by Juergen Knoedlseder                         *
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
 * @file GLATResponse.i
 * @brief Fermi/LAT Response class definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATResponse.hpp"
%}


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Fermi/LAT Response class
 ***************************************************************************/
class GLATResponse : public GResponse {
public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    virtual ~GLATResponse(void);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GLATResponse* clone(void) const;
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

    // Overloaded virtual base class methods
    virtual double        irf_spatial(const GEvent&       event,
                                      const GSource&      source,
                                      const GObservation& obs) const;

    // Other Methods
    int                size(void) const;
    const std::string& rspname(void) const;
    void               load(const std::string& rspname);
    void               save(const std::string& rspname) const;
    const bool&        force_mean(void) const;
    void               force_mean(const bool& value);
    GLATAeff*          aeff(const int& index) const;
    GLATPsf*           psf(const int& index) const;
    GLATEdisp*         edisp(const int& index) const;

    // Reponse methods
    double irf_spatial_atom(const GLATEventAtom& event,
                            const GSource&       source,
                            const GObservation&  obs) const;
    double irf_spatial_bin(const GLATEventBin& event,
                           const GSource&      source,
                           const GObservation& obs) const;
};


/***********************************************************************//**
 * @brief GLATResponse class extension
 ***************************************************************************/
%extend GLATResponse {
    GLATResponse copy() {
        return (*self);
    }
};
