/***************************************************************************
 *                 GCTAEventAtom.hpp - CTA event atom class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
    GCTAEventAtom(const GCTAEventAtom& atom);
    virtual ~GCTAEventAtom(void);

    // Operators
    GCTAEventAtom& operator=(const GCTAEventAtom& atom);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCTAEventAtom*     clone(void) const;
    const GCTAInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
    void               dir(const GCTAInstDir& dir);
    void               energy(const GEnergy& energy);
    void               time(const GTime& time);
    std::string        print(const GChatter& chatter = NORMAL) const;

    // Other methods
    const int&           index(void) const;
    const unsigned long& event_id(void) const;
    const unsigned long& obs_id(void)   const;
    void                 event_id(const unsigned long& id);
    void                 obs_id(const unsigned long& id);

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
    unsigned long m_event_id;       //!< Event identifier
    unsigned long m_obs_id;         //!< Observation identifier
    int           m_multip;         //!< Multiplicity
    char          m_telmask;        //!< Telescope mask
    float         m_dir_err;        //!< Error on event direction
    float         m_alt;            //!< Event altitude
    float         m_az;             //!< Event azimuth
    float         m_corex;          //!< Position on ground (m)
    float         m_corey;          //!< Position on ground (m)
    float         m_core_err;       //!< Error on core reconstruction
    float         m_xmax;           //!< Position of shower max (g/cm2)
    float         m_xmax_err;       //!< Error on shower max (g/cm2)
    float         m_shwidth;        //!< Shower width (m)
    float         m_shlength;       //!< Shower length (m)
    float         m_energy_err;     //!< Error on event energy (MeV)
    float         m_hil_msw;        //!< Hillas MSW
    float         m_hil_msw_err;    //!< Hillas MSW error
    float         m_hil_msl;        //!< Hillas MSL
    float         m_hil_msl_err;    //!< Hillas MSL error
};


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
 * @brief Return observation identifier
 *
 * @return Observation identifier.
 ***************************************************************************/
inline
const unsigned long& GCTAEventAtom::obs_id(void) const
{
    return (m_obs_id);
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
 * @brief Set observation identifier
 *
 * @param[in] id Observation identifier.
 ***************************************************************************/
inline
void GCTAEventAtom::obs_id(const unsigned long& id)
{
    m_obs_id = id;
    return;
}

#endif /* GCTAEVENTATOM_HPP */
