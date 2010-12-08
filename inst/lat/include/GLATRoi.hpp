/***************************************************************************
 *               GLATRoi.hpp  -  LAT region of interest class              *
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
 * @file GLATRoi.hpp
 * @brief GLATRoi class definition.
 * @author J. Knodlseder
 */

#ifndef GLATROI_HPP
#define GLATROI_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GRoi.hpp"
#include "GLATInstDir.hpp"


/***********************************************************************//**
 * @class GLATRoi
 *
 * @brief Interface for the LAT region of interest class.
 *
 * The LAT region of interest class defines the region of photon arrival
 * directions that is used for unbinned data analysis. A circular LAT has
 * been implemented.
 ***************************************************************************/
class GLATRoi : public GRoi {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATRoi& roi);
    friend GLog&         operator<< (GLog& log, const GLATRoi& roi);

public:
    // Constructors and destructors
    GLATRoi(void);
    GLATRoi(const GLATRoi& roi);
    virtual ~GLATRoi(void);

    // Operators
    GLATRoi& operator= (const GLATRoi& roi);

    // Implemented pure virtual base class methods
    void        clear(void);
    GLATRoi*    clone(void) const;
    std::string print(void) const;

    // Other methods
    GLATInstDir centre(void) const { return m_centre; }
    double      radius(void) const { return m_radius; }
    void        centre(const GLATInstDir& centre) { m_centre=centre; return; }
    void        radius(const double& radius) { m_radius=radius; return; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATRoi& roi);
    void free_members(void);
    
    // Protected members
    GLATInstDir m_centre;   //!< Centre of ROI in instrument coordinates
    double      m_radius;   //!< Radius of ROI in degrees
};

#endif /* GLATROI_HPP */
