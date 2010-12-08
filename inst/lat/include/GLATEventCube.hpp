/***************************************************************************
 *                GLATEventCube.hpp  -  LAT event cube class               *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GLATEventCube.hpp
 * @brief GLATEventCube class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENTCUBE_HPP
#define GLATEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <vector>
#include "GEventCube.hpp"
#include "GLATObservation.hpp"
#include "GLATInstDir.hpp"
#include "GLATEventBin.hpp"
#include "GEbounds.hpp"
#include "GGti.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GFits.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"
#include "GSkymap.hpp"
#include "GNodeArray.hpp"


/***********************************************************************//**
 * @class GLATEventCube
 *
 * @brief GLATEventCube class interface defintion.
 ***************************************************************************/
class GLATEventCube : public GEventCube {

public:
    // Constructors and destructors
    GLATEventCube(void);
    GLATEventCube(const GLATEventCube& cube);
    virtual ~GLATEventCube(void);

    // Operators
    GLATEventCube& operator= (const GLATEventCube& cube);

    // Implemented pure virtual base class methods
    void           clear(void);
    GLATEventCube* clone(void) const;
    void           load(const std::string& filename);
    GLATEventBin*  pointer(int index);
    int            number(void) const;
    std::string    print(void) const;

    // Other methods
    GEbounds&   ebds(void) { return m_ebds; }
    GGti&       gti(void) { return m_gti; }
    GTime&      time(void) { return m_time; }
    double      ontime(void) { return m_ontime; }
    int         nx(void) const { return m_map.nx(); }
    int         ny(void) const { return m_map.ny(); }
    int         npix(void) const { return m_map.npix(); }
    int         ebins(void) const { return m_map.nmaps(); }
    int         ndiffrsp(void) const { return m_srcmap.size(); }
    std::string diffname(const int& index) const;
    GSkymap*    diffrsp(const int& index) const;
    GNodeArray* enodes(void) { return &m_enodes; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATEventCube& cube);
    void free_members(void);
    void read_cntmap(GFitsImage* hdu);
    void read_srcmap(GFitsImage* hdu);
    void read_ebds(GFitsTable* hdu);
    void read_gti(GFitsTable* hdu);
    void set_directions(void);
    void set_energies(void);
    void set_time(void);

    // Protected data area
    GLATEventBin             m_bin;          //!< Actual energy bin
    GSkymap                  m_map;          //!< Counts map stored as sky map
    GEbounds                 m_ebds;         //!< Energy boundaries
    GGti                     m_gti;          //!< Good Time Intervals
    GTime                    m_time;         //!< Event cube mean time
    double                   m_ontime;       //!< Event cube ontime (sec)
    double*                  m_counts;       //!< Pointer to skymap pixels
    GLATInstDir*             m_dirs;         //!< Array of event directions
    double*                  m_omega;        //!< Array of solid angles (sr)
    GEnergy*                 m_energies;     //!< Array of log mean energies
    GEnergy*                 m_ewidth;       //!< Array of energy bin widths
    std::vector<GSkymap*>    m_srcmap;       //!< Pointers to source maps
    std::vector<std::string> m_srcmap_names; //!< Source map names
    GNodeArray               m_enodes;       //!< Energy nodes
};

#endif /* GLATEVENTCUBE_HPP */
