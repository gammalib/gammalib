/***************************************************************************
 *                 GLATEventBin.hpp  -  LAT event bin class                *
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
 * @file GLATEventBin.hpp
 * @brief GLATEventBin class interface definition.
 * @author J. Knodlseder
 */

#ifndef GLATEVENTBIN_HPP
#define GLATEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GEventBin.hpp"


/***********************************************************************//**
 * @class GLATEventBin
 *
 * @brief GLATEventBin class interface defintion.
 ***************************************************************************/
class GLATEventBin : public GEventBin {

    // Friend classes
    friend class GLATEventCube;

public:
    // Constructors and destructors
    GLATEventBin();
    GLATEventBin(const GLATEventBin& bin);
    virtual ~GLATEventBin();

    // Operators
    GLATEventBin& operator= (const GLATEventBin& bin);

    // Methods
    
protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GLATEventBin& bin);
    void          free_members(void);
    GLATEventBin* clone(void) const;

    // Protected data area (defines all LAT specific event attributes)

private:
};

#endif /* GLATEVENTBIN_HPP */
