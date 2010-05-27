/***************************************************************************
 *                 GLATPointing.hpp  -  LAT pointing class                 *
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
 * @file GLATPointing.hpp
 * @brief GLATPointing class definition.
 * @author J. Knodlseder
 */

#ifndef GLATPOINTING_HPP
#define GLATPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include "GPointing.hpp"


/***********************************************************************//**
 * @class GLATPointing
 *
 * @brief Interface for the LAT pointing class.
 *
 * The LAT pointing class contains information for a specific LAT pointing.
 ***************************************************************************/
class GLATPointing : public GPointing {

public:
    // Constructors and destructors
    GLATPointing(void);
    GLATPointing(const GLATPointing& pnt);
    ~GLATPointing(void);

    // Operators
    GLATPointing& operator= (const GLATPointing& pnt);

    // Methods

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GLATPointing& pnt);
    void          free_members(void);
    GLATPointing* clone(void) const;
};

#endif /* GLATPOINTING_HPP */
