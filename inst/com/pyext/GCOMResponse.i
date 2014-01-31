/***************************************************************************
 *                 GCOMResponse.i  -  COMPTEL Response class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2012 by Juergen Knoedlseder                              *
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
 * @file GCOMResponse.hpp
 * @brief COMPTEL instrument response function class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GCOMResponse.hpp"
%}


/***********************************************************************//**
 * @class GCOMResponse
 *
 * @brief Interface for the COMPTEL instrument response function
 ***************************************************************************/
class GCOMResponse : public GResponse {
public:
    // Constructors and destructors
    GCOMResponse(void);
    GCOMResponse(const GCOMResponse& rsp);
    explicit GCOMResponse(const std::string& iaqname, const std::string& caldb = "");
    virtual ~GCOMResponse(void);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GCOMResponse* clone(void) const;
    virtual bool          use_edisp(void) const;
    virtual bool          use_tdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;

    // Other Methods
    void        caldb(const std::string& caldb);
    std::string caldb(void) const;
    std::string iaqname(void) const;
    void        load(const std::string& iaqname);
    void        read(const GFitsImage& hdu);
};


/***********************************************************************//**
 * @brief GCOMResponse class extension
 ***************************************************************************/
%extend GCOMResponse {
    GCOMResponse copy() {
        return (*self);
    }
};
