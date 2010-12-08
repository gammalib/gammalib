/***************************************************************************
 *                GCTAEventAtom.hpp  -  CTA event atom class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
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
    std::string        print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAEventAtom& atom);
    void free_members(void);

    // Protected members
    GCTAInstDir m_dir;            //!< Event direction
    GEnergy     m_energy;         //!< Event energy
    GTime       m_time;           //!< Event time
    long        m_event_id;       //!< Event identifier
    char        m_flags;          //!< Flags
    int         m_multip;         //!< Multiplicity
    char        m_telmask;        //!< Telescope mask
    float       m_dir_err;        //!< Error on event direction
    float       m_detx;           //!<
    float       m_dety;           //!<
    float       m_alt_pnt;        //!< Pointing altitude
    float       m_az_pnt;         //!< Pointing azimuth
    float       m_alt;            //!< Event altitude
    float       m_az;             //!< Event azimuth
    float       m_corex;          //!<
    float       m_corey;          //!<
    float       m_core_err;       //!<
    float       m_xmax;           //!<
    float       m_xmax_err;       //!<
    float       m_energy_err;     //!< Error on event energy (MeV)
    float       m_hil_msw;        //!< Hillas width
    float       m_hil_msw_err;    //!< Error on Hillas width
    float       m_hil_msl;        //!< Hillas length
    float       m_hil_msl_err;    //!< Error on Hillas length
};

#endif /* GCTAEVENTATOM_HPP */
