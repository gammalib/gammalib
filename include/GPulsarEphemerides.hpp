/***************************************************************************
 *            GPulsarEphemerides.hpp - Pulsar ephemerides class            *
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
 * @file GPulsarEphemerides.hpp
 * @brief Pulsar ephemerides class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPULSAREPHEMERIDES_HPP
#define GPULSAREPHEMERIDES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GTime.hpp"
#include "GSkyDir.hpp"
#include "GModelTemporalPhaseCurve.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GPulsarEphemerides
 *
 * @brief Pulsar ephemerides class
 *
 * @todo Add class description.
 ***************************************************************************/
class GPulsarEphemerides : public GBase {

public:
    // Constructors and destructors
    GPulsarEphemerides(void);
    GPulsarEphemerides(const GPulsarEphemerides& ephemerides);
    virtual ~GPulsarEphemerides(void);

    // Operators
    GPulsarEphemerides& operator=(const GPulsarEphemerides& ephemerides);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GPulsarEphemerides* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const GTime& tstart(void) const;
    const GTime& tstop(void) const;
    GTime        t0(void) const;
    void         t0(const GTime& t0);
    double       f0(void) const;
    void         f0(const double& f0);
    double       f1(void) const;
    void         f1(const double& f1);
    double       f2(void) const;
    void         f2(const double& f2);
    double       phase(const GTime& time) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPulsarEphemerides& ephemerides);
    void free_members(void);
    
    // Protected members
    std::string              m_name;         //!< Pulsar name
    GTime                    m_tstart;       //!< Validity start time
    GTime                    m_tstop;        //!< Validity stop time
    GSkyDir                  m_dir;          //!< Pulsar sky direction
    GModelTemporalPhaseCurve m_phase_curve;  //!< Pulsar phase curve
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPulsarEphemerides").
 ***************************************************************************/
inline
std::string GPulsarEphemerides::classname(void) const
{
    return ("GPulsarEphemerides");
}


/***********************************************************************//**
 * @brief Returns validity start time
 *
 * @return Validity start time.
 ***************************************************************************/
inline
const GTime& GPulsarEphemerides::tstart(void) const
{
    // Return
    return m_tstart;
}


/***********************************************************************//**
 * @brief Returns validity stop time
 *
 * @return Validity stop time.
 ***************************************************************************/
inline
const GTime& GPulsarEphemerides::tstop(void) const
{
    // Return
    return m_tstop;
}


/***********************************************************************//**
 * @brief Returns time of phase 0
 *
 * @return Reference time of phase 0.
 ***************************************************************************/
inline
GTime GPulsarEphemerides::t0(void) const
{
    return (m_phase_curve.mjd());
}


/***********************************************************************//**
 * @brief Set time of phase 0
 *
 * @param[in] t0 Reference time of phase 0.
 ***************************************************************************/
inline
void GPulsarEphemerides::t0(const GTime& t0)
{
    m_phase_curve.mjd(t0);
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar frequency (Hz)
 *
 * @return Pulsar frequency (Hz).
 ***************************************************************************/
inline
double GPulsarEphemerides::f0(void) const
{
    return (m_phase_curve.f0());
}


/***********************************************************************//**
 * @brief Set pulsar frequency (Hz)
 *
 * @param[in] f0 Pulsar frequency (Hz).
 ***************************************************************************/
inline
void GPulsarEphemerides::f0(const double& f0)
{
    m_phase_curve.f0(f0);
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar frequency derivative (s^-2)
 *
 * @return Pulsar frequency derivative (s^-2).
 ***************************************************************************/
inline
double GPulsarEphemerides::f1(void) const
{
    return (m_phase_curve.f1());
}


/***********************************************************************//**
 * @brief Set pulsar frequency derivative (s^-2)
 *
 * @param[in] f1 Pulsar frequency derivative (s^-2).
 ***************************************************************************/
inline
void GPulsarEphemerides::f1(const double& f1)
{
    m_phase_curve.f1(f1);
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar second frequency derivative (s^-3)
 *
 * @return Pulsar second frequency derivative (s^-3).
 ***************************************************************************/
inline
double GPulsarEphemerides::f2(void) const
{
    return (m_phase_curve.f2());
}


/***********************************************************************//**
 * @brief Set pulsar second frequency derivative (s^-3)
 *
 * @param[in] f2 Pulsar second frequency derivative (s^-3).
 ***************************************************************************/
inline
void GPulsarEphemerides::f2(const double& f2)
{
    m_phase_curve.f2(f2);
    return;
}


/***********************************************************************//**
 * @brief Returns pulsar phase
 *
 * @return Pulsar phase.
 ***************************************************************************/
inline
double GPulsarEphemerides::phase(const GTime& time) const
{
    return (m_phase_curve.phase(time));
}

#endif /* GPULSAREPHEMERIDES_HPP */
