/***************************************************************************
 *                 GCTAEventAtom.hpp - CTA event atom class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2021 by Juergen Knoedlseder                         *
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
 * @file GCTAEventAtom.hpp
 * @brief CTA event atom class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCTAEVENTATOM_HPP
#define GCTAEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventAtom.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPolarization.hpp"
#include "GCTAInstDir.hpp"


/***********************************************************************//**
 * @class GCTAEventAtom
 *
 * @brief CTA event atom class
 *
 * This class implements a CTA event atom. It collects all the relevant event
 * information needed for CTA unbinned analysis.
 ***************************************************************************/
class GCTAEventAtom : public GEventAtom {

    // Friend classes
    friend class GCTAEventList;

public:
    // Constructors and destructors
    GCTAEventAtom(void);
    GCTAEventAtom(const GCTAInstDir& dir, const GEnergy& energy, const GTime& time);
    GCTAEventAtom(const GCTAEventAtom& atom);
    virtual ~GCTAEventAtom(void);

    // Operators
    GCTAEventAtom& operator=(const GCTAEventAtom& atom);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GCTAEventAtom*       clone(void) const;
    virtual std::string          classname(void) const;
    virtual const GCTAInstDir&   dir(void) const;
    virtual const GEnergy&       energy(void) const;
    virtual const GTime&         time(void) const;
    virtual const GPolarization& polarization(void) const;
    virtual std::string          print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const int&           index(void) const;
    const unsigned long& event_id(void) const;
    const int&           mc_id(void) const;
    const float&         phase(void) const;
    void                 index(const int& index);
    void                 event_id(const unsigned long& id);
    void                 mc_id(const int& id);
    void                 phase(const float& phase);
    void                 dir(const GCTAInstDir& dir);
    void                 energy(const GEnergy& energy);
    void                 time(const GTime& time);
    void                 polarization(const GPolarization& polarization);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAEventAtom& atom);
    void free_members(void);

    // Protected members
    int           m_index;          //!< Index in list
    GCTAInstDir   m_dir;            //!< Event direction
    GEnergy       m_energy;         //!< Event energy
    GTime         m_time;           //!< Event time
    GPolarization m_polarization;   //!< Event polarization
    unsigned long m_event_id;       //!< Event identifier
    int           m_mc_id;          //!< Monte Carlo identifier
    float         m_phase;          //!< Optional phase
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCTAEventAtom").
 ***************************************************************************/
inline
std::string GCTAEventAtom::classname(void) const
{
    return ("GCTAEventAtom");
}


/***********************************************************************//**
 * @brief Return instrument direction
 *
 * @return Instrument direction.
 ***************************************************************************/
inline
const GCTAInstDir& GCTAEventAtom::dir(void) const
{
    return (m_dir);
}


/***********************************************************************//**
 * @brief Return energy
 *
 * @return Energy.
 ***************************************************************************/
inline
const GEnergy& GCTAEventAtom::energy(void) const
{
    return (m_energy);
}


/***********************************************************************//**
 * @brief Return time
 *
 * @return Time.
 ***************************************************************************/
inline
const GTime& GCTAEventAtom::time(void) const
{
    return (m_time);
}


/***********************************************************************//**
 * @brief Return polarization
 *
 * @return Polarization.
 ***************************************************************************/
inline
const GPolarization& GCTAEventAtom::polarization(void) const
{
    return (m_polarization);
}


/***********************************************************************//**
 * @brief Set instrument direction
 *
 * @param[in] dir Instrument direction.
 ***************************************************************************/
inline
void GCTAEventAtom::dir(const GCTAInstDir& dir)
{
    m_dir = dir;
    return;
}


/***********************************************************************//**
 * @brief Set energy
 *
 * @param[in] energy Energy.
 ***************************************************************************/
inline
void GCTAEventAtom::energy(const GEnergy& energy)
{
    m_energy = energy;
    return;
}


/***********************************************************************//**
 * @brief Set time
 *
 * @param[in] time Time.
 ***************************************************************************/
inline
void GCTAEventAtom::time(const GTime& time)
{
    m_time = time;
    return;
}


/***********************************************************************//**
 * @brief Set polarization
 *
 * @param[in] polarization Polarization.
 ***************************************************************************/
inline
void GCTAEventAtom::polarization(const GPolarization& polarization)
{
    m_polarization = polarization;
    return;
}


/***********************************************************************//**
 * @brief Return event index in list
 *
 * @return Index.
 *
 * Returns the index of the event in case it is part of a list.
 ***************************************************************************/
inline
const int& GCTAEventAtom::index(void) const
{
    return (m_index);
}


/***********************************************************************//**
 * @brief Return event identifier
 *
 * @return Event identifier.
 ***************************************************************************/
inline
const unsigned long& GCTAEventAtom::event_id(void) const
{
    return (m_event_id);
}


/***********************************************************************//**
 * @brief Return Monte Carlo identifier
 *
 * @return Monte Carlo identifier.
 ***************************************************************************/
inline
const int& GCTAEventAtom::mc_id(void) const
{
    return (m_mc_id);
}


/***********************************************************************//**
 * @brief Return event phase
 *
 * @return Event phase.
 ***************************************************************************/
inline
const float& GCTAEventAtom::phase(void) const
{
    return (m_phase);
}


/***********************************************************************//**
 * @brief Set event index
 *
 * @param[in] index Event index.
 ***************************************************************************/
inline
void GCTAEventAtom::index(const int& index)
{
    m_index = index;
    return;
}


/***********************************************************************//**
 * @brief Set event identifier
 *
 * @param[in] id Event identifier.
 ***************************************************************************/
inline
void GCTAEventAtom::event_id(const unsigned long& id)
{
    m_event_id = id;
    return;
}


/***********************************************************************//**
 * @brief Set Monte Carlo identifier
 *
 * @param[in] id Monte Carlo identifier.
 ***************************************************************************/
inline
void GCTAEventAtom::mc_id(const int& id)
{
    m_mc_id = id;
    return;
}


/***********************************************************************//**
 * @brief Set event phase
 *
 * @param[in] phase Event phase.
 ***************************************************************************/
inline
void GCTAEventAtom::phase(const float& phase)
{
    m_phase = phase;
    return;
}

#endif /* GCTAEVENTATOM_HPP */
