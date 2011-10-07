/***************************************************************************
 *                GCTAEventAtom.hpp  -  CTA event atom class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
 * @brief GCTAEventAtom class interface definition.
 * @author J. Knodlseder
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
 * @brief GCTAEventAtom class interface defintion
 *
 * This class implement a CTA event atom. It collects all the relevant event
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
    GCTAEventAtom& operator= (const GCTAEventAtom& atom);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCTAEventAtom*     clone(void) const;
    const GCTAInstDir& dir(void) const { return m_dir; }
    const GEnergy&     energy(void) const { return m_energy; }
    const GTime&       time(void) const { return m_time; }
    void               dir(const GCTAInstDir& dir) { m_dir=dir; }
    void               energy(const GEnergy& energy) { m_energy=energy; }
    void               time(const GTime& time) { m_time=time; }
    std::string        print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAEventAtom& atom);
    void free_members(void);

    // Protected members
    GCTAInstDir   m_dir;            //!< Event direction
    GEnergy       m_energy;         //!< Event energy
    GTime         m_time;           //!< Event time
    unsigned long m_event_id;       //!< Event identifier
    unsigned long m_obs_id;         //!< Observation identifier
    int           m_multip;         //!< Multiplicity
    char          m_telmask;        //!< Telescope mask
    float         m_dir_err;        //!< Error on event direction
    float         m_detx;           //!< Tangential coordinate in nominal sys
    float         m_dety;           //!< Tangential coordinate in nominal sys
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
};

#endif /* GCTAEVENTATOM_HPP */
