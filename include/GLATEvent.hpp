/***************************************************************************
 *                   GLATEvent.hpp  -  LAT Event class                     *
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
 * @file GLATEvent.hpp
 * @brief GLATEvent class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENT_HPP
#define GLATEVENT_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvent.hpp"


/***********************************************************************//**
 * @class GLATEvent
 *
 * @brief GLATEvent class interface defintion.
 ***************************************************************************/
class GLATEvent : public GEvent {

    // Friend classes
    friend class GLATEventList;
    friend class GLATEventCube;

public:
    // Constructors and destructors
    GLATEvent();
    GLATEvent(const GLATEvent& event);
    virtual ~GLATEvent();

    // Operators
    GLATEvent& operator= (const GLATEvent& event);

    // Methods
    
protected:
    // Protected methods
    void       init_members(void);
    void       copy_members(const GLATEvent& event);
    void       free_members(void);
    GLATEvent* clone(void) const;

    // Protected data area (defines all LAT specific event attributes)
    float  m_theta;                //!< Zenith angle in instrument system
    float  m_phi;                  //!< Azimuth angle in instrument system
    float  m_zenith_angle;         //!< Zenith angle in Earth system
    float  m_earth_azimuth_angle;  //!< Azimuth angle in Earth system
    short  m_event_class;          //!< Event class

private:
};

#endif /* GLATEVENTS_HPP */
