/***************************************************************************
 *                 GEvents.hpp  -  Events abstract base class              *
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
 * @file GEvents.hpp
 * @brief GEvents abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTS_HPP
#define GEVENTS_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GEvents
 *
 * @brief Abstract interface for the events classes.
 ***************************************************************************/
class GEvents {

	// Friend classes
    friend class GObservation;

public:
    // Constructors and destructors
    GEvents();
    GEvents(const GEvents& events);
    virtual ~GEvents();

    // Operators
    virtual GEvents& operator= (const GEvents& events);

    // Methods
  
protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GEvents& events);
    void    free_members(void);
    virtual GEvents* clone(void) const = 0;

    // Protected data area

private:
};

#endif /* GEVENTS_HPP */
