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
#include "GSkymap.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
//#include "GFits.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCTAInstDir.hpp"
#include "GCTAPointing.hpp"
#include "GCTAResponse.hpp"


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
	void          load(const std::string& filename);
    GCTAEventBin* pointer(int index);
    int           number(void) const;
    int           ebins(void) const { return m_map.nmaps(); }
    int           pixels(void) const { return m_map.npix(); }
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GCTAEventCube& cube);
    void           free_members(void);
    GCTAEventCube* clone(void) const;
    void           load_cntmap(GFitsHDU* hdu);
    void           load_ebds(GFitsHDU* hdu);
    void           load_gti(GFitsHDU* hdu);

    // Protected data area
    GCTAEventBin m_bin;         //!< Actual event bin
    GSkymap      m_map;         //!< Counts map stored as sky map
    double*      m_counts;      //!< Pointer to skymap pixels
    GCTAInstDir* m_dirs;        //!< Array of events directions
    GEnergy*     m_energies;    //!< Array of log mean energies
    GTime        m_time;        //!< Event cube mean time
    GCTAPointing m_pnt;         //!< CTA pointing
    GCTAResponse m_rsp;         //!< CTA instrument response function
    GEbounds     m_ebds;        //!< Energy boundaries
    GGti         m_gti;         //!< Good Time Intervals
};

#endif /* GCTAEVENTCUBE_HPP */
