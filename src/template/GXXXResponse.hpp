/***************************************************************************
 *               GXXXResponse.hpp - [INSTRUMENT] response class            *
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
 * @file GXXXResponse.hpp
 * @brief [INSTRUMENT] instrument response function class definition
 * @author [AUTHOR]
 */

#ifndef GXXXRESPONSE_HPP
#define GXXXRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GResponse.hpp"
//#include "GCaldb.hpp"

/* __ Type definitions ___________________________________________________ */

/* __ Forward declaration ________________________________________________ */
class GEvent;
class GPhoton;
class GEnergy;
class GTime;
class GSource;
class GObservation;
class GModelSky;
class GEbounds;


/***********************************************************************//**
 * @class GXXXResponse
 *
 * @brief [INSTRUMENT] instrument response function class
 *
 * The [INSTRUMENT] instrument response function class defines the function
 * that translates from physical quantities to measured quantities.
 *
 * @todo Complete the class description.
 ***************************************************************************/
class GXXXResponse : public GResponse {

public:
    // Constructors and destructors
    GXXXResponse(void);
    GXXXResponse(const GXXXResponse& rsp);
    virtual ~GXXXResponse(void);

    // Operators
    virtual GXXXResponse& operator=(const GXXXResponse & rsp);

    // Implement pure virtual base class methods
    virtual void          clear(void);
    virtual GXXXResponse* clone(void) const;
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
    void copy_members(const GXXXResponse& rsp);
    void free_members(void);

    // Private data members
    // TODO: Add any data members that are necessary. Note that the events
    // are stored in the GObservation base class
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXXXResponse").
 ***************************************************************************/
inline
std::string GXXXResponse::classname(void) const
{
    return ("GXXXResponse");
}


/***********************************************************************//**
 * @brief Signal if energy dispersion will be used
 *
 * @return False.
 *
 * @todo Implement method as needed.
 ***************************************************************************/
inline
bool GXXXResponse::use_edisp(void) const
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
bool GXXXResponse::use_tdisp(void) const
{
    return false;
}

#endif /* GXXXRESPONSE_HPP */
