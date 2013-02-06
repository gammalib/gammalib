/***************************************************************************
 *                   GXXXResponse.i  -  XXX Response class                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2013 by Juergen Knoedlseder                              *
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
 * @file GXXXResponse.hpp
 * @brief XXX instrument response function class interface definition
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GXXXResponse.hpp"
%}


/***********************************************************************//**
 * @class GXXXResponse
 *
 * @brief Interface for the XXX instrument response function
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
    virtual bool          hasedisp(void) const;
    virtual bool          hastdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;
};


/***********************************************************************//**
 * @brief GXXXResponse class extension
 ***************************************************************************/
%extend GXXXResponse {
    GXXXResponse copy() {
        return (*self);
    }
};
