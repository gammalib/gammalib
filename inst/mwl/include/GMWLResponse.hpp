/***************************************************************************
 *          GMWLResponse.hpp  -  Multi-wavelength response class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2013 by Juergen Knoedlseder                         *
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
 * @file GMWLResponse.hpp
 * @brief Multi-wavelength response class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GMWLRESPONSE_HPP
#define GMWLRESPONSE_HPP

/* __ Includes ___________________________________________________________ */
#include "GResponse.hpp"


/***********************************************************************//**
 * @class GMWLResponse
 *
 * @brief Multi-wavelength response class
 *
 * This class implements a dummy response class for multi-wavelength
 * observations. Since the multi-wavelength instrument classes handles data
 * that a provided in photon space, no instrument response is in fact needed.
 * The dummy response implemented by this class provides a simple diagonal
 * response matrix that allows integration of multi-wavelength observations
 * using the standard instrument specific interface.
 ***************************************************************************/
class GMWLResponse : public GResponse {

public:
    // Constructors and destructors
    GMWLResponse(void);
    GMWLResponse(const GMWLResponse& rsp);
    virtual ~GMWLResponse(void);

    // Operators
    virtual GMWLResponse& operator= (const GMWLResponse & rsp);

    // Implemented pure virtual methods
    virtual void          clear(void);
    virtual GMWLResponse* clone(void) const;
    virtual bool          hasedisp(void) const;
    virtual bool          hastdisp(void) const;
    virtual double        irf(const GEvent&       event,
                              const GPhoton&      photon,
                              const GObservation& obs) const;
    virtual double        npred(const GPhoton&      photon,
                                const GObservation& obs) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLResponse& pnt);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Signal if response supports energy dispersion
 *
 * @return True if response supports energy dispersion.
 ***************************************************************************/
inline
bool GMWLSpectrum::hasedisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Signal if response supports time dispersion
 *
 * @return True if response supports time dispersion.
 ***************************************************************************/
inline
bool GMWLSpectrum::hastdisp(void) const
{
    return false;
}


/***********************************************************************//**
 * @brief Return instrument response function
 *
 * @param[in] event Event.
 * @param[in] photon Photon.
 * @param[in] obs Observation.
 *
 * @return Instrument response function (always 1).
 ***************************************************************************/
inline
double GMWLSpectrum::irf(const GEvent& event, const GPhoton& photon,
                         const GObservation& obs) const
{
    return 1.0;
}


/***********************************************************************//**
 * @brief Return predicted number of events
 *
 * @param[in] photon Photon.
 * @param[in] obs Observation.
 *
 * @return Instrument response function (always 1).
 ***************************************************************************/
inline
double GMWLSpectrum::npred(const GPhoton& photon, const GObservation& obs) const
{
    return 1.0;
}

#endif /* GMWLRESPONSE_HPP */
