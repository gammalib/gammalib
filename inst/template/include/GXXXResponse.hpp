/***************************************************************************
 *                  GXXXResponse.hpp  -  XXX Response class                *
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

#ifndef GXXXRESPONSE_HPP
#define GXXXRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declaration ________________________________________________ */


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

    // Operators
    virtual GXXXResponse& operator= (const GXXXResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GXXXResponse* clone(void) const;
    virtual bool          hasedisp(void) const { return false; }
    virtual bool          hastdisp(void) const { return false; }
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;
    virtual std::string   print(void) const;

private:
    // Private methods
    void init_members(void);
    void copy_members(const GXXXResponse& rsp);
    void free_members(void);

    // Private data members
};

#endif /* GXXXRESPONSE_HPP */
