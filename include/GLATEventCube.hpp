/***************************************************************************
 *                GLATEventCube.hpp  -  LAT event cube class               *
 * ----------------------------------------------------------------------- *
 *  copyright            : (C) 2009 by Jurgen Knodlseder                   *
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
#include "GFits.hpp"


/***********************************************************************//**
 * @class GLATEventCube
 *
 * @brief GLATEventCube class interface defintion.
 ***************************************************************************/
class GLATEventCube : public GEventCube {

public:
    // Constructors and destructors
    GLATEventCube();
    GLATEventCube(const GLATEventCube& cube);
    virtual ~GLATEventCube();

    // Operators
    GLATEventCube& operator= (const GLATEventCube& cube);

    // Methods
	void          load(const std::string& filename);
    void          load(GFitsHDU* hdu);
    GLATEventBin* pointer(int index) const;
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GLATEventCube& cube);
    void           free_members(void);
    GLATEventCube* clone(void) const;

    // Protected data area
    GLATEventBin* m_bins;             //!< Pointer to bins

private:
};

#endif /* GLATEVENTCUBE_HPP */
