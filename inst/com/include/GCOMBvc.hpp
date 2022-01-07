/***************************************************************************
 *         GCOMBvc.hpp - COMPTEL Solar System Barycentre Data class        *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knodlseder                               *
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
 * @file GCOMBvc.hpp
 * @brief COMPTEL Solar System Barycentre Data class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMBVC_HPP
#define GCOMBVC_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GTime.hpp"
#include "GVector.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMBvc
 *
 * @brief COMPTEL Solar System Barycentre Data class
 *
 * The class holds one record of a COMPTEL Solar System Barycentre Data file.
 ***************************************************************************/
class GCOMBvc : public GBase {

public:
    // Constructors and destructors
    GCOMBvc(void);
    GCOMBvc(const GCOMBvc& bvc);
    virtual ~GCOMBvc(void);

    // Operators
    GCOMBvc& operator=(const GCOMBvc& bvc);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMBvc*    clone(void) const;
    virtual std::string classname(void) const;
    virtual std::string print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GTime&   tstart(void) const;
    void           tstart(const GTime& tstart);
    const GTime&   tstop(void) const;
    void           tstop(const GTime& tstop);
    const int&     tjd(void) const;
    void           tjd(const int& tjd);
    const int&     tics(void) const;
    void           tics(const int& tics);
    const GVector& ssb(void) const;
    void           ssb(const GVector& ssb);
    const double&  tdelta(void) const;
    void           tdelta(const double& tdelta);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMBvc& bvc);
    void free_members(void);

    // Protected members
    GTime   m_tstart;  //!< Start time of superpacket
    GTime   m_tstop;   //!< Stop time of superpacket
    int     m_tjd;     //!< TJD of BVC record
    int     m_tics;    //!< Tics of BVC record
    GVector m_ssb;     //!< Solar System Barycentre vector in celestial system (km)
    double  m_tdelta;  //!< Time difference TDB-UTC (sec)
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMBvc").
 ***************************************************************************/
inline
std::string GCOMBvc::classname(void) const
{
    return ("GCOMBvc");
}


/***********************************************************************//**
 * @brief Return start time of superpacket
 *
 * @return Start time of superpacket.
 *
 * Returns the start time of the superpacket.
 ***************************************************************************/
inline
const GTime& GCOMBvc::tstart(void) const
{
    return (m_tstart);
}


/***********************************************************************//**
 * @brief Set start time of superpacket
 *
 * @param[in] tstart Start time of superpacket.
 *
 * Set the start time of the superpacket.
 ***************************************************************************/
inline
void GCOMBvc::tstart(const GTime& tstart)
{
    m_tstart = tstart;
    return;
}


/***********************************************************************//**
 * @brief Return stop time of superpacket
 *
 * @return Stop time of superpacket.
 *
 * Returns the stop time of the superpacket. The stop time is defined as the
 * start time plus 131071 tics, since the length of one superpacket is
 * 16.384 secs, i.e. 16.384 * 8000 = 131072 ticks.
 ***************************************************************************/
inline
const GTime& GCOMBvc::tstop(void) const
{
    return (m_tstop);
}


/***********************************************************************//**
 * @brief Set stop time of superpacket
 *
 * @param[in] tstop Stop time of superpacket.
 *
 * Set the stop time of the superpacket.
 ***************************************************************************/
inline
void GCOMBvc::tstop(const GTime& tstop)
{
    m_tstop = tstop;
    return;
}


/***********************************************************************//**
 * @brief Return Truncated Julian Days of Orbit Aspect Record
 *
 * @return Truncated Julian Days of Orbit Aspect Record.
 *
 * Returns the Truncated Julian Days of the Orbit Aspect Record.
 ***************************************************************************/
inline
const int& GCOMBvc::tjd(void) const
{
    return (m_tjd);
}


/***********************************************************************//**
 * @brief Set Truncated Julian Days of Orbit Aspect Record
 *
 * @param[in] tjd Truncated Julian Days of Orbit Aspect Record.
 *
 * Set the Truncated Julian Days of the Orbit Aspect Record.
 ***************************************************************************/
inline
void GCOMBvc::tjd(const int& tjd)
{
    m_tjd = tjd;
    return;
}


/***********************************************************************//**
 * @brief Return tics of Orbit Aspect Record
 *
 * @return Tics of Orbit Aspect Record.
 *
 * Returns the tics of the Orbit Aspect Record.
 ***************************************************************************/
inline
const int& GCOMBvc::tics(void) const
{
    return (m_tics);
}


/***********************************************************************//**
 * @brief Set tics of Orbit Aspect Record
 *
 * @param[in] tics Tics of Orbit Aspect Record.
 *
 * Set the tics of the Orbit Aspect Record.
 ***************************************************************************/
inline
void GCOMBvc::tics(const int& tics)
{
    m_tics = tics;
    return;
}


/***********************************************************************//**
 * @brief Return Solar System Barycentre vector
 *
 * @return Solar System Barycentre vector (km).
 *
 * Returns the Solar System Barycentre vector in km. The vector is given in
 * the celestial system.
 ***************************************************************************/
inline
const GVector& GCOMBvc::ssb(void) const
{
    return (m_ssb);
}


/***********************************************************************//**
 * @brief Set Solar System Barycentre vector
 *
 * @param[in] ssb Solar System Barycentre vector (km).
 *
 * Set the Solar System Barycentre vector in km. The vector is defined in
 * the celestial system.
 ***************************************************************************/
inline
void GCOMBvc::ssb(const GVector& ssb)
{
    m_ssb = ssb;
    return;
}


/***********************************************************************//**
 * @brief Return TDB-UTC time difference
 *
 * @return TDB-UTC time difference (s).
 *
 * Returns the TDB-UTC time difference in seconds.
 ***************************************************************************/
inline
const double& GCOMBvc::tdelta(void) const
{
    return (m_tdelta);
}


/***********************************************************************//**
 * @brief Set TDB-UTC time difference
 *
 * @param[in] tdelta TDB-UTC time difference (s).
 *
 * Set the TDB-UTC time difference.
 ***************************************************************************/
inline
void GCOMBvc::tdelta(const double& tdelta)
{
    m_tdelta = tdelta;
    return;
}

#endif /* GCOMBVC_HPP */
