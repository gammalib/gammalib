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
#include <string>
#include "GPointing.hpp"
#include "GSkyDir.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GCTAPointing
 *
 * @brief Interface for the CTA pointing class.
 *
 * This class implements a CTA pointing. For the time being it is assumed
 * that the pointing direction is time-independent.
 ***************************************************************************/
class GCTAPointing : public GPointing {

public:
    // Constructors and destructors
    GCTAPointing(void);
    explicit GCTAPointing(const GSkyDir& dir);
    GCTAPointing(const GCTAPointing& pnt);
    virtual ~GCTAPointing(void);

    // Operators
    GCTAPointing& operator= (const GCTAPointing& pnt);

    // Implemented pure virtual methods
    void           clear(void);
    GCTAPointing*  clone(void) const;
    const GSkyDir& dir(void) const { return m_dir; }
    std::string    print(void) const;

    // Other methods
    void dir(const GSkyDir& dir);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAPointing& pnt);
    void free_members(void);

    // Protected members
    GSkyDir m_dir;  //!< Pointing direction
};

#endif /* GCTAPOINTING_HPP */
