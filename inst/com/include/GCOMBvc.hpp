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
    const GTime&   time(void) const;
    void           time(const GTime& time);
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
    GTime   m_time;    //!< Time for Solar System Barycentre Data
    int     m_tjd;     //!< TJD of Solar System Barycentre Data
    int     m_tics;    //!< Tics of Solar System Barycentre Data
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
 * @brief Return time of Solar System Barycentre Data
 *
 * @return Time of Solar System Barycentre Data.
 *
 * Returns the time of the Solar System Barycentre Data. The time is defined
 * as the mid-point of the corresponding superpacket.
 ***************************************************************************/
inline
const GTime& GCOMBvc::time(void) const
{
    return (m_time);
}


/***********************************************************************//**
 * @brief Set time of Solar System Barycentre Data
 *
 * @param[in] time Time of Solar System Barycentre Data.
 *
 * Set the time of the Solar System Barycentre Data. The time is defined as
 * the mid-point of the corresponding superpacket.
 ***************************************************************************/
inline
void GCOMBvc::time(const GTime& time)
{
    m_time = time;
    return;
}


/***********************************************************************//**
 * @brief Return Truncated Julian Days of Solar System Barycentre Data
 *
 * @return Truncated Julian Days of Solar System Barycentre Data.
 *
 * Returns the Truncated Julian Days of Solar System Barycentre Data.
 ***************************************************************************/
inline
const int& GCOMBvc::tjd(void) const
{
    return (m_tjd);
}


/***********************************************************************//**
 * @brief Set Truncated Julian Days of Solar System Barycentre Data
 *
 * @param[in] tjd Truncated Julian Days of Solar System Barycentre Data.
 *
 * Set the Truncated Julian Days of the Solar System Barycentre Data.
 ***************************************************************************/
inline
void GCOMBvc::tjd(const int& tjd)
{
    m_tjd = tjd;
    return;
}


/***********************************************************************//**
 * @brief Return tics of Solar System Barycentre Data
 *
 * @return Tics of Solar System Barycentre Data.
 *
 * Returns the tics of the Solar System Barycentre Data.
 ***************************************************************************/
inline
const int& GCOMBvc::tics(void) const
{
    return (m_tics);
}


/***********************************************************************//**
 * @brief Set tics of Solar System Barycentre Data
 *
 * @param[in] tics Tics of Solar System Barycentre Data.
 *
 * Set the tics of the Solar System Barycentre Data.
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
