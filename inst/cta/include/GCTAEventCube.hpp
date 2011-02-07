/***************************************************************************
 *                GCTAEventCube.hpp  -  CTA event cube class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
#include <string>
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

public:
    // Constructors and destructors
    GCTAEventCube(void);
    GCTAEventCube(const GSkymap& map);
    GCTAEventCube(const GCTAEventCube& cube);
    virtual ~GCTAEventCube(void);

    // Operators
    GCTAEventCube& operator= (const GCTAEventCube& cube);

    // Implemented pure virtual base class methods
    void            clear(void);
    GCTAEventCube*  clone(void) const;
    int             size(void) const;
    int             dim(void) const;
    int             naxis(int axis) const;
    void            load(const std::string& filename);
    GCTAEventBin*   pointer(int index);
    int             number(void) const;
    std::string     print(void) const;


    // Other methods
    void            write(GFits* file) const;
    void            map(const GSkymap& map) { m_map=map; }
    void            ebounds(const GEbounds& ebds) { m_ebds=ebds; }
    void            gti(const GGti& gti) { m_gti=gti; }
    const GSkymap&  map(void) const { return m_map; }
    const GEbounds& ebounds(void) const { return m_ebds; }
    const GGti&     gti(void) const { return m_gti; }
    int             nx(void) const { return m_map.nx(); }
    int             ny(void) const { return m_map.ny(); }
    int             npix(void) const { return m_map.npix(); }
    int             ebins(void) const { return m_map.nmaps(); }

protected:
    // Protected methods
    void         init_members(void);
    void         copy_members(const GCTAEventCube& cube);
    void         free_members(void);
    void         read_cntmap(GFitsImage* hdu);
    void         read_ebds(GFitsTable* hdu);
    void         read_gti(GFitsTable* hdu);
    void         set_directions(void);
    void         set_energies(void);
    void         set_time(void);

    // Protected fundamental data
    GSkymap      m_map;            //!< Counts map stored as sky map
    GEbounds     m_ebds;           //!< Energy boundaries
    GGti         m_gti;            //!< Good Time Intervals

    // Protected derived data
    GCTAEventBin m_bin;            //!< Actual event bin
    GCTAInstDir* m_dirs;           //!< Array of event directions
    double*      m_omega;          //!< Array of solid angles (sr)
    GEnergy*     m_energies;       //!< Array of log mean energies
    GEnergy*     m_ewidth;         //!< Array of energy bin widths
    GTime        m_time;           //!< Event cube mean time
    double       m_ontime;         //!< Event cube ontime (sec)
};

#endif /* GCTAEVENTCUBE_HPP */
