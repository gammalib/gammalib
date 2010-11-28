/***************************************************************************
 *              GPointing.hpp  -  Pointing abstract base class             *
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
 * @file GPointing.hpp
 * @brief GPointing abstract base class definition.
 * @author J. Knodlseder
 */

#ifndef GPOINTING_HPP
#define GPOINTING_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GPointing
 *
 * @brief Abstract interface for the pointing classes.
 *
 * The pointing classes hold instrument specific information about the
 * pointing of the telescope.
 ***************************************************************************/
class GPointing {

public:
    // Constructors and destructors
    GPointing(void);
    GPointing(const GPointing& pnt);
    virtual ~GPointing(void);

    // Operators
    virtual GPointing& operator= (const GPointing& pnt);

    // Pure virtual methods
    virtual GPointing* clone(void) const = 0;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPointing& pnt);
    void free_members(void);
};

#endif /* GPOINTING_HPP */
