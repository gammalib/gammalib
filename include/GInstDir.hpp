/***************************************************************************
 *         GInstDir.hpp  -  Instrument direction abstract base class       *
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
 * @file GInstDir.hpp
 * @brief GInstDir abstract base class definition.
 * @author J. Knodlseder
 */

#ifndef GINSTDIR_HPP
#define GINSTDIR_HPP

/* __ Includes ___________________________________________________________ */


/***********************************************************************//**
 * @class GInstDir
 *
 * @brief Abstract interface for the instrument direction of an event.
 *
 * The instrument direction of an event is the equivalent of the sky
 * direction (implemented by GSkyDir) but in the instrument data space.
 * The instrument direction may be any kind of position or direction
 * information encoded in the data space, such as incident event
 * reconstructions for imaging devices or detector numbers etc. for
 * non-imaging devices.
 ***************************************************************************/
class GInstDir {

public:
    // Constructors and destructors
    GInstDir(void);
    GInstDir(const GInstDir& dir);
    virtual ~GInstDir(void);

    // Operators
    virtual GInstDir& operator= (const GInstDir& dir);

    // Pure virtual methods

protected:
    // Protected methods
    void              init_members(void);
    void              copy_members(const GInstDir& dir);
    void              free_members(void);
    virtual GInstDir* clone(void) const = 0;
};

#endif /* GINSTDIR_HPP */
