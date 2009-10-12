/***************************************************************************
 *                GLATEventAtom.hpp  -  LAT event atom class               *
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
 * @file GLATEventAtom.hpp
 * @brief GLATEventAtom class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENTATOM_HPP
#define GLATEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventAtom.hpp"


/***********************************************************************//**
 * @class GLATEventAtom
 *
 * @brief GLATEventAtom class interface defintion.
 ***************************************************************************/
class GLATEventAtom : public GEventAtom {

    // Friend classes
    friend class GLATEventList;
//    friend class GLATEventCube;

public:
    // Constructors and destructors
    GLATEventAtom();
    GLATEventAtom(const GLATEventAtom& atom);
    virtual ~GLATEventAtom();

    // Operators
    GLATEventAtom& operator= (const GLATEventAtom& atom);

    // Methods
    
protected:
    // Protected methods
    void           init_members(void);
    void           copy_members(const GLATEventAtom& atom);
    void           free_members(void);
    GLATEventAtom* clone(void) const;

    // Protected data area (defines all LAT specific event attributes)
    float  m_theta;                //!< Zenith angle in instrument system
    float  m_phi;                  //!< Azimuth angle in instrument system
    float  m_zenith_angle;         //!< Zenith angle in Earth system
    float  m_earth_azimuth_angle;  //!< Azimuth angle in Earth system
    short  m_event_class;          //!< Event class

private:
};

#endif /* GLATEVENTATOM_HPP */
