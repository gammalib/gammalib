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
    const GTime& tstart(void) const;
    void         tstart(const GTime& tstart);
    const GTime& tstop(void) const;
    void         tstop(const GTime& tstop);
    const int&   tjd(void) const;
    void         tjd(const int& tjd);
    const int&   tics(void) const;
    void         tics(const int& tics);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMOad& oad);
    void free_members(void);
    
    // Protected members
    GTime m_tstart; //!< Start time of superpacket
    GTime m_tstop;  //!< Stop time of superpacket
    int   m_tjd;    //!< TJD of OAD record
    int   m_tics;   //!< Tics of OAD record
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

#endif /* GCOMOAD_HPP */
