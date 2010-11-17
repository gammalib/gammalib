/***************************************************************************
 *                GCTAEventCube.hpp  -  CTA event cube class               *
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
 * @file GCTAEventCube.hpp
 * @brief GCTAEventCube class interface definition.
 * @author J. Knodlseder
 */

#ifndef GCTAEVENTCUBE_HPP
#define GCTAEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventCube.hpp"
#include "GCTAEventBin.hpp"
#include "GCTAObservation.hpp"
#include "GSkymap.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAPointing.hpp"
#include "GFitsTable.hpp"
#include "GFitsImage.hpp"


/***********************************************************************//**
 * @class GCTAEventCube
 *
 * @brief GCTAEventCube class interface defintion.
 *
 * The GCTAEventCube class holds event information for the CTA binned
 * analysis. This information is handled on a bin-by-bin basis by the
 * GCTAEventBin class.
 ***************************************************************************/
class GCTAEventCube : public GEventCube {

    // Friend classes
    friend class GCTAObservation;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAEventCube& cube);

public:
    // Constructors and destructors
    GCTAEventCube(void);
    GCTAEventCube(const GCTAEventCube& cube);
    virtual ~GCTAEventCube(void);

    // Operators
    GCTAEventCube& operator= (const GCTAEventCube& cube);

    // Methods
    void             clear(void);
    void             load(const std::string& filename);
    GCTAEventBin*    pointer(int index);
    int              number(void) const;
    int              nx(void) const { return m_map.nx(); }
    int              ny(void) const { return m_map.ny(); }
    int              npix(void) const { return m_map.npix(); }
    int              ebins(void) const { return m_map.nmaps(); }
    void             obs(GCTAObservation* ptr) { m_obs=ptr; return; }
    GCTAObservation* obs(void) const { return m_obs; }

protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GCTAEventCube& cube);
    void           free_members(void);
    GCTAEventCube* clone(void) const;
    void           read_cntmap(GFitsImage* hdu);
    void           read_ebds(GFitsTable* hdu);
    void           read_gti(GFitsTable* hdu);
    void           set_directions(void);
    void           set_energies(void);
    void           set_time(void);

    // Protected fundamental data
    GSkymap          m_map;            //!< Counts map stored as sky map
    GEbounds         m_ebds;           //!< Energy boundaries
    GGti             m_gti;            //!< Good Time Intervals

    // Protected derived data
    GCTAEventBin     m_bin;            //!< Actual event bin
    double*          m_counts;         //!< Pointer to skymap pixels
    GCTAInstDir*     m_dirs;           //!< Array of event directions
    double*          m_omega;          //!< Array of solid angles (sr)
    GEnergy*         m_energies;       //!< Array of log mean energies
    GEnergy*         m_ewidth;         //!< Array of energy bin widths
    GTime            m_time;           //!< Event cube mean time
    double           m_ontime;         //!< Event cube ontime (sec)
    GCTAPointing     m_pnt;            //!< CTA pointing
    GCTAObservation* m_obs;            //!< Points back to CTA observation
};

#endif /* GCTAEVENTCUBE_HPP */
