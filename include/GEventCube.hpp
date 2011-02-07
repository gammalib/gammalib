/***************************************************************************
 *          GEventCube.hpp  -  Abstract event cube container class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEventCube.hpp
 * @brief GEventCube container class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTCUBE_HPP
#define GEVENTCUBE_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEvents.hpp"
#include "GEventBin.hpp"


/***********************************************************************//**
 * @class GEventCube
 *
 * @brief GEventCube container class interface defintion.
 *
 * @todo Remove islist() and iscube() and use dynamic_cast instead. These
 *       methods have also to be removed in the GEvents class and all classes
 *       that derive from GEvents.
 * @todo Replace pointer() method by access operators. Also here we have to
 *       go back to the GEvents class and change this consistently
 *       thoughout all classes.
 ***************************************************************************/
class GEventCube : public GEvents {

public:
    // Constructors and destructors
    GEventCube(void);
    GEventCube(const GEventCube& cube);
    virtual ~GEventCube(void);

    // Operators
    virtual GEventCube& operator= (const GEventCube& cube);

    // Pure virtual methods
    virtual void        clear(void) = 0;
    virtual GEventCube* clone(void) const = 0;
    virtual int         size(void) const = 0;
    virtual int         dim(void) const = 0;
    virtual int         naxis(int axis) const = 0;
    virtual void        load(const std::string& filename) = 0;
    virtual GEventBin*  pointer(int index) = 0;
    virtual int         number(void) const = 0;
    virtual std::string print(void) const = 0;

    // Implemented pure virtual base class methods
    bool islist(void) const { return false; }
    bool iscube(void) const { return true; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEventCube& cube);
    void free_members(void);
};

#endif /* GEVENTCUBE_HPP */
