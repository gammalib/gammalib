/***************************************************************************
 *              GCOMOad.hpp - COMPTEL Orbit Aspect Data class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knodlseder                               *
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
 * @file GCOMOad.hpp
 * @brief COMPTEL Orbit Aspect Data class definition
 * @author Juergen Knodlseder
 */

#ifndef GCOMOAD_HPP
#define GCOMOAD_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GTime.hpp"
#include "GSkyDir.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GCOMOad
 *
 * @brief COMPTEL Orbit Aspect Data class
 *
 * The class holds one record of a COMPTEL Orbit Aspect Data file.
 ***************************************************************************/
class GCOMOad : public GBase {

public:
    // Constructors and destructors
    GCOMOad(void);
    GCOMOad(const GCOMOad& oad);
    virtual ~GCOMOad(void);

    // Operators
    GCOMOad& operator=(const GCOMOad& oad);

    // Implemented pure virtual base class methods
    virtual void        clear(void);
    virtual GCOMOad*    clone(void) const;
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
    const float&   gcaz(void) const;
    void           gcaz(const float& gcaz);
    const float&   gcel(void) const;
    void           gcel(const float& gcel);
    const float&   georad(void) const;
    void           georad(const float& georad);
    const GSkyDir& zaxis(void) const;
    void           zaxis(const GSkyDir& zaxis);
    const GSkyDir& xaxis(void) const;
    void           xaxis(const GSkyDir& xaxis);
    double         theta(const GSkyDir& sky) const;
    double         phi(const GSkyDir& sky) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMOad& oad);
    void free_members(void);

    // Protected members
    GTime   m_tstart;  //!< Start time of superpacket
    GTime   m_tstop;   //!< Stop time of superpacket
    GSkyDir m_zaxis;   //!< Telescope z-axis
    GSkyDir m_xaxis;   //!< Telescope x-axis
    int     m_tjd;     //!< TJD of OAD record
    int     m_tics;    //!< Tics of OAD record
    float   m_gcaz;    //!< Geocentre azimuth angle (deg)
    float   m_gcel;    //!< Geocentre zenith angle (deg)
    float   m_georad;  //!< Apparent radius of Earth (deg)

    // Precomputation cache
    mutable double m_posang; //!< X-axis position angle in COMPTEL system
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMOad").
 ***************************************************************************/
inline
std::string GCOMOad::classname(void) const
{
    return ("GCOMOad");
}


/***********************************************************************//**
 * @brief Return start time of superpacket
 *
 * @return Start time of superpacket.
 *
 * Returns the start time of the superpacket.
 ***************************************************************************/
inline
const GTime& GCOMOad::tstart(void) const
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
void GCOMOad::tstart(const GTime& tstart)
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
const GTime& GCOMOad::tstop(void) const
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
void GCOMOad::tstop(const GTime& tstop)
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
const int& GCOMOad::tjd(void) const
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
void GCOMOad::tjd(const int& tjd)
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
const int& GCOMOad::tics(void) const
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
void GCOMOad::tics(const int& tics)
{
    m_tics = tics;
    return;
}


/***********************************************************************//**
 * @brief Return Geocentre azimuth angle
 *
 * @return Geocentre azimuth angle (deg).
 *
 * Returns the Geocentre azimuth angle in degrees.
 ***************************************************************************/
inline
const float& GCOMOad::gcaz(void) const
{
    return (m_gcaz);
}


/***********************************************************************//**
 * @brief Set Geocentre azimuth angle
 *
 * @param[in] gcaz Geocentre azimuth angle (deg).
 *
 * Set the Geocentre azimuth angle.
 ***************************************************************************/
inline
void GCOMOad::gcaz(const float& gcaz)
{
    m_gcaz = gcaz;
    return;
}


/***********************************************************************//**
 * @brief Return Geocentre zenith angle
 *
 * @return Geocentre zenith angle (deg).
 *
 * Returns the Geocentre zenith angle in degrees.
 ***************************************************************************/
inline
const float& GCOMOad::gcel(void) const
{
    return (m_gcel);
}


/***********************************************************************//**
 * @brief Set Geocentre zenith angle
 *
 * @param[in] gcel Geocentre zenith angle (deg).
 *
 * Set the Geocentre zenith angle.
 ***************************************************************************/
inline
void GCOMOad::gcel(const float& gcel)
{
    m_gcel = gcel;
    return;
}


/***********************************************************************//**
 * @brief Return apparent radius of Earth
 *
 * @return Apparent radius of Earth (deg).
 *
 * Returns the apparent radius of Earth in degrees.
 ***************************************************************************/
inline
const float& GCOMOad::georad(void) const
{
    return (m_georad);
}


/***********************************************************************//**
 * @brief Set apparent radius of Earth
 *
 * @param[in] georad Apparent radius of Earth (deg).
 *
 * Set the apparent radius of Earth.
 ***************************************************************************/
inline
void GCOMOad::georad(const float& georad)
{
    m_georad = georad;
    return;
}


/***********************************************************************//**
 * @brief Return telescope Z-axis
 *
 * @return Telescope Z-axis.
 *
 * Returns the telescope Z-axis.
 ***************************************************************************/
inline
const GSkyDir& GCOMOad::zaxis(void) const
{
    return (m_zaxis);
}


/***********************************************************************//**
 * @brief Set telescope Z-axis
 *
 * @param[in] zaxis Telescope Z-axis.
 *
 * Set the telescope Z-axis.
 ***************************************************************************/
inline
void GCOMOad::zaxis(const GSkyDir& zaxis)
{
    m_posang = 1.0e30; // To assure initialisation of position angle
    m_zaxis  = zaxis;
    return;
}


/***********************************************************************//**
 * @brief Return telescope X-axis
 *
 * @return Telescope X-axis.
 *
 * Returns the telescope X-axis.
 ***************************************************************************/
inline
const GSkyDir& GCOMOad::xaxis(void) const
{
    return (m_xaxis);
}


/***********************************************************************//**
 * @brief Set telescope X-axis
 *
 * @param[in] xaxis Telescope X-axis.
 *
 * Set the telescope X-axis.
 ***************************************************************************/
inline
void GCOMOad::xaxis(const GSkyDir& xaxis)
{
    m_posang = 1.0e30; // To assure initialisation of position angle
    m_xaxis  = xaxis;
    return;
}


/***********************************************************************//**
 * @brief Return zenith angle of sky direction in COMPTEL coordinates
 *
 * @param[in] sky Sky direction.
 * @return Zenith angle of sky direction in COMPTEL coordinates (deg).
 *
 * Returns the zenith angle of a sky direction in COMPTEL coordinates.
 ***************************************************************************/
inline
double GCOMOad::theta(const GSkyDir& sky) const
{
    return (m_zaxis.dist_deg(sky));
}


/***********************************************************************//**
 * @brief Return azimuth angle of sky direction in COMPTEL coordinates
 *
 * @param[in] sky Sky direction.
 * @return Azimuth angle of sky direction in COMPTEL coordinates (deg).
 *
 * Returns the azimuth angle of a sky direction in COMPTEL coordinates.
 ***************************************************************************/
inline
double GCOMOad::phi(const GSkyDir& sky) const
{
    // If position angle has not be initialised the do it now
    if (m_posang > 1.0e20) {
        m_posang = m_zaxis.posang_deg(m_xaxis);
    }
    //return (m_zaxis.posang_deg(m_xaxis) - m_zaxis.posang_deg(sky));
    return (m_posang - m_zaxis.posang_deg(sky));
}

#endif /* GCOMOAD_HPP */
