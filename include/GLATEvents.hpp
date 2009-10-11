/***************************************************************************
 *                  GLATEvents.hpp  -  LAT Events class                    *
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
 * @file GLATEvents.hpp
 * @brief GLATEvents class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENTS_HPP
#define GLATEVENTS_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvents.hpp"
#include "GLATEvent.hpp"
#include "GFits.hpp"


/***********************************************************************//**
 * @class GLATEvents
 *
 * @brief GLATEvents class interface defintion.
 ***************************************************************************/
class GLATEvents : public GEvents {

public:
    // Constructors and destructors
    GLATEvents();
    GLATEvents(const GLATEvents& events);
    virtual ~GLATEvents();

    // Operators
    GLATEvents& operator= (const GLATEvents& events);

    // Methods
	void       load(const std::string& filename);
    void       load(GFitsHDU* hdu);
    GLATEvent* pointer(int index) const;
    
protected:
    // Protected methods
    void        init_members(void);
    void        copy_members(const GLATEvents& events);
    void        free_members(void);
    GLATEvents* clone(void) const;

    // Protected data area

private:
};

#endif /* GLATEVENTS_HPP */
