/***************************************************************************
 *             GSPIObservation.hpp  -  SPI observation class               *
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
 * @file GSPIObservation.hpp
 * @brief SPI observation class interface definition
 * @author Juergen Knoedlseder
 */

#ifndef GSPIOBSERVATION_HPP
#define GSPIOBSERVATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GTime.hpp"
#include "GXmlElement.hpp"
#include "GObservation.hpp"
#include "GResponse.hpp"
#include "GSPIResponse.hpp"
#include "GSPIPointing.hpp"


/***********************************************************************//**
 * @class GSPIObservation
 *
 * @brief Interface class for SPI observations
 ***************************************************************************/
class GSPIObservation : public GObservation {

public:
    // Constructors and destructors
    GSPIObservation(void);
    GSPIObservation(const GSPIObservation& obs);
    virtual ~GSPIObservation(void);

    // Operators
    virtual GSPIObservation& operator=(const GSPIObservation& obs);

    // Implement pure virtual methods
    virtual void             clear(void);
    virtual GSPIObservation* clone(void) const;
    virtual void             response(const GResponse& rsp);
    virtual GSPIResponse*    response(void) const;
    virtual GSPIPointing*    pointing(void) const;
    virtual std::string      instrument(void) const { return m_instrument; }
    virtual double           ontime(void) const { return m_ontime; }
    virtual double           livetime(void) const { return m_livetime; }
    virtual double           deadc(const GTime& time) const { return m_deadc; }
    virtual void             read(const GXmlElement& xml);
    virtual void             write(GXmlElement& xml) const;
    virtual std::string      print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GSPIObservation& obs);
    void free_members(void);

    // Protected members
    std::string   m_instrument;  //!< Instrument name
    double        m_ontime;      //!< Ontime (sec)
    double        m_livetime;    //!< Livetime (sec)
    double        m_deadc;       //!< Deadtime correction
    GSPIPointing* m_pointing;    //!< Pointer to pointing direction
    GSPIResponse* m_response;    //!< Pointer to response functions
};

#endif /* GSPIOBSERVATION_HPP */
