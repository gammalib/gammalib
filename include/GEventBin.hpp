/***************************************************************************
 *               GEventBin.hpp  -  Event bin abstract base class           *
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
 * @file GEventBin.hpp
 * @brief GEventBin abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTBIN_HPP
#define GEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvent.hpp"


/***********************************************************************//**
 * @class GEventBin
 *
 * @brief Abstract interface for the event bin class
 ***************************************************************************/
class GEventBin : public GEvent {

    // Friend classes
    friend class GEvents;

public:
    // Constructors and destructors
    GEventBin();
    GEventBin(const GEventBin& bin);
    virtual ~GEventBin();

    // Operators
    virtual GEventBin& operator= (const GEventBin& bin);

    // Virtual methods
    virtual std::string string(void) const = 0;
    
    // Implemented methods
    
protected:
    // Protected methods
    void               init_members(void);
    void               copy_members(const GEventBin& bin);
    void               free_members(void);
    virtual GEventBin* clone(void) const = 0;

    // Protected data area

private:
};

#endif /* GEVENTBIN_HPP */
