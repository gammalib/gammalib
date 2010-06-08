/***************************************************************************
 *                 GCTAPointing.hpp  -  CTA pointing class                 *
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
 * @file GCTAPointing.hpp
 * @brief GCTAPointing class definition.
 * @author J. Knodlseder
 */

#ifndef GCTAPOINTING_HPP
#define GCTAPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GPointing.hpp"


/***********************************************************************//**
 * @class GCTAPointing
 *
 * @brief Interface for the CTA pointing class.
 *
 * The CTA pointing class contains information for a specific CTA pointing.
 ***************************************************************************/
class GCTAPointing : public GPointing {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GCTAPointing& pnt);

public:
    // Constructors and destructors
    GCTAPointing(void);
    GCTAPointing(const GCTAPointing& pnt);
    ~GCTAPointing(void);

    // Operators
    GCTAPointing& operator= (const GCTAPointing& pnt);

    // Methods
    void clear(void);

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GCTAPointing& pnt);
    void          free_members(void);
    GCTAPointing* clone(void) const;
};

#endif /* GCTAPOINTING_HPP */
