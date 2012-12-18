/***************************************************************************
 *                GLATResponse.i - Fermi-LAT response class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2008-2012 by Juergen Knoedlseder                         *
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
 * @brief Fermi-LAT response class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GLATResponse.hpp"
%}


/***********************************************************************//**
 * @class GLATResponse
 *
 * @brief Fermi-LAT response class
 ***************************************************************************/
class GLATResponse : public GResponse {

public:
    // Constructors and destructors
    GLATResponse(void);
    GLATResponse(const GLATResponse& rsp);
    virtual ~GLATResponse(void);

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GLATResponse* clone(void) const;
    virtual bool          hasedisp(void) const;
    virtual bool          hastdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;

    // Implemented virtual methods
    virtual double irf(const GEvent&       event,
                       const GSource&      source,
                       const GObservation& obs) const;

    // Other Methods
    void        caldb(const std::string& caldb);
    std::string caldb(void) const;
    std::string rspname(void) const;
    void        load(const std::string& rspname);
    int         size(void) const;
    GLATAeff*   aeff(const int& index) const;
    GLATPsf*    psf(const int& index) const;
    GLATEdisp*  edisp(const int& index) const;
    void        save(const std::string& rspname) const;
    bool        force_mean(void);
    void        force_mean(const bool& value);

    // Reponse methods
    double irf(const GLATEventAtom& event,
               const GSource&       source,
               const GObservation&  obs) const;
    double irf(const GLATEventBin& event,
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
