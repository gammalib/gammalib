/***************************************************************************
 *                  GSPIResponse.hpp  -  SPI Response class                *
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
 * @file GSPIResponse.hpp
 * @brief SPI instrument response function class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIRESPONSE_HPP
#define GSPIRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEvent.hpp"
#include "GPhoton.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declaration ________________________________________________ */


/***********************************************************************//**
 * @class GSPIResponse
 *
 * @brief Interface for the SPI instrument response function
 ***************************************************************************/
class GSPIResponse : public GResponse {

public:
    // Constructors and destructors
    GSPIResponse(void);
    GSPIResponse(const GSPIResponse& rsp);
    virtual ~GSPIResponse(void);

    // Operators
    virtual GSPIResponse& operator= (const GSPIResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GSPIResponse* clone(void) const;
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
    void copy_members(const GSPIResponse& rsp);
    void free_members(void);

    // Private data members
};

#endif /* GSPIRESPONSE_HPP */
