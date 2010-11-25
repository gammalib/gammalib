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
#include "GEventCube.hpp"
#include "GLATEventBin.hpp"
#include "GLATObservation.hpp"
#include "GEbounds.hpp"
#include "GFits.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GLATInstDir.hpp"
#include "GLATPointing.hpp"
#include "GLATResponse.hpp"
#include "GFitsImage.hpp"
#include "GFitsTable.hpp"


/***********************************************************************//**
 * @class GLATEventCube
 *
 * @brief GLATEventCube class interface defintion.
 ***************************************************************************/
class GLATEventCube : public GEventCube {

    // Friend classes
    friend class GLATObservation;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATEventCube& cube);

public:
    // Constructors and destructors
    GLATEventCube(void);
    GLATEventCube(const GLATEventCube& cube);
    virtual ~GLATEventCube(void);

    // Operators
    GLATEventCube& operator= (const GLATEventCube& cube);

    // Implemented pure virtul methods
    void           clear(void);
    GLATEventCube* clone(void) const;
    void           load(const std::string& filename);
    GLATEventBin*  pointer(int index);
    int            number(void) const;

    // Other methods
    int ebins(void) const { return m_ebins; }
    int pixels(void) const { return m_pixels; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATEventCube& cube);
    void free_members(void);
    void load_cntmap(GFitsImage* hdu);
    void load_ebds(GFitsTable* hdu);

    // Protected data area
    GLATEventBin m_bin;            //!< Actual energy bin
    int          m_nx;             //!< Number of pixels in x
    int          m_ny;             //!< Number of pixels in y
    int          m_pixels;         //!< Number of pixels in x and y
    int          m_ebins;          //!< Number of energy bins
    double*      m_counts;         //!< Pointer to counts array
    GEnergy*     m_energies;       //!< Pointer to energies
    GTime        m_time;           //!< Event cube mean time
    GLATInstDir* m_dirs;           //!< Pointer to event directions
    GEbounds     m_ebds;           //!< Energy boundaries

private:
};

#endif /* GLATEVENTCUBE_HPP */
