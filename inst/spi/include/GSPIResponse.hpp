/***************************************************************************
 *              GSPIResponse.hpp - INTEGRAL/SPI response class             *
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
 * @file GSPIResponse.hpp
 * @brief INTEGRAL/SPI instrument response function class definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIRESPONSE_HPP
#define GSPIRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GResponse.hpp"

/* __ Forward declaration ________________________________________________ */
class GEvent;
class GPhoton;
class GEnergy;
class GTime;
class GSource;
class GObservation;
class GModelSky;
class GEbounds;

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GSPIResponse
 *
 * @brief INTEGRAL/SPI instrument response function class
 *
 * The INTEGRAL/SPI instrument response function class defines the function
 * that translates from physical quantities to measured quantities.
 *
 * @todo Complete the class description.
 ***************************************************************************/
class GSPIResponse : public GResponse {

public:
    // Constructors and destructors
    GSPIResponse(void);
    GSPIResponse(const GSPIResponse& rsp);
    virtual ~GSPIResponse(void);

    // Operators
    virtual GSPIResponse& operator=(const GSPIResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GSPIResponse* clone(void) const;
    virtual std::string   classname(void) const;
    virtual bool          use_edisp(void) const;
    virtual bool          use_tdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        irf(const GEvent&       event,
                              const GSource&      source,
                              const GObservation& obs) const;
    virtual double        nroi(const GModelSky&    model,
                               const GEnergy&      obsEng,
                               const GTime&        obsTime,
                               const GObservation& obs) const;
    virtual GEbounds      ebounds(const GEnergy& obsEnergy) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Other Methods
    // TODO: Add any further methods that are needed

private:
    // Private methods
    void init_members(void);
    void copy_members(const GSPIResponse& rsp);
    void free_members(void);

    // Private data members
    // TODO: Add any data members that are necessary. Note that the events
    // are stored in the GObservation base class
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GSPIResponse").
 ***************************************************************************/
inline
std::string GSPIResponse::classname(void) const
{
    return ("GSPIResponse");
}


/***********************************************************************//**
 * @brief Signal if energy dispersion will be used
 *
 * @return False.
 *
 * @todo Implement method as needed.
 ***************************************************************************/
inline
bool GSPIResponse::use_edisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if time dispersion will be used
 *
 * @return False.
 *
 * @todo Implement method as needed.
 ***************************************************************************/
inline
bool GSPIResponse::use_tdisp(void) const
{
    return false;
}

#endif /* GSPIRESPONSE_HPP */
