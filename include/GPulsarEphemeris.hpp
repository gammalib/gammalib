/***************************************************************************
 *              GPulsarEphemeris.hpp - Pulsar ephemeris class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GPulsarEphemeris.hpp
 * @brief Pulsar ephemeris class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPULSAREPHEMERIS_HPP
#define GPULSAREPHEMERIS_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GTime.hpp"
#include "GSkyDir.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GPulsarEphemeris
 *
 * @brief Pulsar ephemeris class
 *
 * This class implements an ephemeris for a pulsar.
 ***************************************************************************/
class GPulsarEphemeris : public GBase {

public:
    // Constructors and destructors
    GPulsarEphemeris(void);
    GPulsarEphemeris(const GPulsarEphemeris& ephemeris);
    virtual ~GPulsarEphemeris(void);

    // Operators
    GPulsarEphemeris& operator=(const GPulsarEphemeris& ephemeris);

    // Implemented pure virtual base class methods
    virtual void              clear(void);
    virtual GPulsarEphemeris* clone(void) const;
    virtual std::string       classname(void) const;
    virtual std::string       print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const std::string& name(void) const;
    void               name(const std::string& name);
    const GSkyDir&     dir(void) const;
    void               dir(const GSkyDir& dir);
    const GTime&       tstart(void) const;
    void               tstart(const GTime& tstart);
    const GTime&       tstop(void) const;
    void               tstop(const GTime& tstop);
    const std::string& timesys(void) const;
    void               timesys(const std::string& timesys);
    GTime              t0(void) const;
    void               t0(const GTime& t0);
    double             phase(void) const;
    void               phase(const double& phase);
    double             f0(void) const;
    void               f0(const double& f0);
    double             f1(void) const;
    void               f1(const double& f1);
    double             f2(void) const;
    void               f2(const double& f2);
    double             phase(const GTime&       time,
                             const std::string& timesys) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPulsarEphemeris& ephemeris);
    void free_members(void);
    
    // Protected members
    std::string m_name;    //!< Pulsar name
    GTime       m_tstart;  //!< Validity start time
    GTime       m_tstop;   //!< Validity stop time
    GSkyDir     m_dir;     //!< Pulsar sky direction
    std::string m_timesys; //!< Time system of pulsar ephemeris
    GTime       m_t0;      //!< Reference epoch of pulsar ephemeris
    double      m_phase;   //!< Pulse phase
    double      m_f0;      //!< Pulsar frequency (Hz)
    double      m_f1;      //!< Pulsar frequency derivative (s^-2)
    double      m_f2;      //!< Pulsar second frequency derivative (s^-3)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPulsarEphemeris").
 ***************************************************************************/
inline
std::string GPulsarEphemeris::classname(void) const
{
    return ("GPulsarEphemeris");
}


/***********************************************************************//**
 * @brief Returns pulsar name
 *
 * @return Pulsar name.
 ***************************************************************************/
inline
const std::string& GPulsarEphemeris::name(void) const
{
    // Return
    return m_name;
}


/***********************************************************************//**
 * @brief Set pulsar name
 *
 * @param[in] name Pulsar name.
 ***************************************************************************/
inline
void GPulsarEphemeris::name(const std::string& name)
{
    m_name = name;
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar sky direction
 *
 * @return Pulsar sky direction.
 ***************************************************************************/
inline
const GSkyDir& GPulsarEphemeris::dir(void) const
{
    // Return
    return m_dir;
}


/***********************************************************************//**
 * @brief Set pulsar sky direction
 *
 * @param[in] dir Pulsar sky direction.
 ***************************************************************************/
inline
void GPulsarEphemeris::dir(const GSkyDir& dir)
{
    m_dir = dir;
    return;
}


/***********************************************************************//**
 * @brief Returns validity start time
 *
 * @return Validity start time.
 ***************************************************************************/
inline
const GTime& GPulsarEphemeris::tstart(void) const
{
    // Return
    return m_tstart;
}


/***********************************************************************//**
 * @brief Set validity start time
 *
 * @param[in] tstart Validity start time.
 ***************************************************************************/
inline
void GPulsarEphemeris::tstart(const GTime& tstart)
{
    m_tstart = tstart;
    return;
}


/***********************************************************************//**
 * @brief Returns validity stop time
 *
 * @return Validity stop time.
 ***************************************************************************/
inline
const GTime& GPulsarEphemeris::tstop(void) const
{
    // Return
    return m_tstop;
}


/***********************************************************************//**
 * @brief Set validity stop time
 *
 * @param[in] tstop Validity stop time.
 ***************************************************************************/
inline
void GPulsarEphemeris::tstop(const GTime& tstop)
{
    m_tstop = tstop;
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar ephemeris time system
 *
 * @return Pulsar ephemeris time system.
 *
 * The pulsar ephemeris time system is one of "TT", "TAI" or "UTC".
 ***************************************************************************/
inline
const std::string& GPulsarEphemeris::timesys(void) const
{
    // Return
    return (m_timesys);
}


/***********************************************************************//**
 * @brief Set pulsar ephemeris time system
 *
 * @param[in] timesys Pulsar ephemeris time system.
 *
 * The pulsar ephemeris time system is one of "TT", "TAI" or "UTC".
 ***************************************************************************/
inline
void GPulsarEphemeris::timesys(const std::string& timesys)
{
    m_timesys = timesys;
    return;
}


/***********************************************************************//**
 * @brief Returns reference epoch of pulsar ephemeris
 *
 * @return Reference epoch of pulsar ephemeris.
 ***************************************************************************/
inline
GTime GPulsarEphemeris::t0(void) const
{
    return (m_t0);
}


/***********************************************************************//**
 * @brief Set reference epoch of pulsar ephemeris
 *
 * @param[in] t0 Reference epoch of pulsar ephemeris.
 ***************************************************************************/
inline
void GPulsarEphemeris::t0(const GTime& t0)
{
    m_t0 = t0;
    return;
}


/***********************************************************************//**
 * @brief Returns pulse phase
 *
 * @return Pulse phase.
 ***************************************************************************/
inline
double GPulsarEphemeris::phase(void) const
{
    return (m_phase);
}


/***********************************************************************//**
 * @brief Set pulse phase
 *
 * @param[in] phase Pulse phase.
 ***************************************************************************/
inline
void GPulsarEphemeris::phase(const double& phase)
{
    m_phase = phase;
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar frequency (Hz)
 *
 * @return Pulsar frequency (Hz).
 ***************************************************************************/
inline
double GPulsarEphemeris::f0(void) const
{
    return (m_f0);
}


/***********************************************************************//**
 * @brief Set pulsar frequency (Hz)
 *
 * @param[in] f0 Pulsar frequency (Hz).
 ***************************************************************************/
inline
void GPulsarEphemeris::f0(const double& f0)
{
    m_f0 = f0;
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar frequency derivative (s^-2)
 *
 * @return Pulsar frequency derivative (s^-2).
 ***************************************************************************/
inline
double GPulsarEphemeris::f1(void) const
{
    return (m_f1);
}


/***********************************************************************//**
 * @brief Set pulsar frequency derivative (s^-2)
 *
 * @param[in] f1 Pulsar frequency derivative (s^-2).
 ***************************************************************************/
inline
void GPulsarEphemeris::f1(const double& f1)
{
    m_f1 = f1;
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar second frequency derivative (s^-3)
 *
 * @return Pulsar second frequency derivative (s^-3).
 ***************************************************************************/
inline
double GPulsarEphemeris::f2(void) const
{
    return (m_f2);
}


/***********************************************************************//**
 * @brief Set pulsar second frequency derivative (s^-3)
 *
 * @param[in] f2 Pulsar second frequency derivative (s^-3).
 ***************************************************************************/
inline
void GPulsarEphemeris::f2(const double& f2)
{
    m_f2 = f2;
    return;
}

#endif /* GPULSAREPHEMERIS_HPP */
