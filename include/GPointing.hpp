/***************************************************************************
 *              GPointing.hpp  -  Pointing abstract base class             *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GPointing.hpp
 * @brief Pointing abstract base class definition
 * @author J. Knodlseder
 */

#ifndef GPOINTING_HPP
#define GPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GSkyDir.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GPointing
 *
 * @brief Abstract interface for the pointing classes
 *
 * The pointing class holds information about the time dependent telescope
 * pointing.
 ***************************************************************************/
class GPointing {

    // I/O friends
    friend std::ostream& operator<<(std::ostream& os, const GPointing& pnt);
    friend GLog&         operator<<(GLog& log,        const GPointing& pnt);

public:
    // Constructors and destructors
    GPointing(void);
    GPointing(const GPointing& pnt);
    virtual ~GPointing(void);

    // Operators
    virtual GPointing& operator= (const GPointing& pnt);

    // Pure virtual methods
    virtual void           clear(void) = 0;
    virtual GPointing*     clone(void) const = 0;
    virtual const GSkyDir& dir(void) const = 0;
    virtual std::string    print(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPointing& pnt);
    void free_members(void);
};

#endif /* GPOINTING_HPP */
