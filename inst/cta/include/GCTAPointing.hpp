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
#include "GMatrix.hpp"


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
    virtual GCTAPointing& operator=(const GCTAPointing& pnt);

    // Implemented pure virtual methods
    virtual void           clear(void);
    virtual GCTAPointing*  clone(void) const;
    virtual const GSkyDir& dir(void) const { return m_dir; }
    virtual std::string    print(void) const;

    // Other methods
    void  dir(const GSkyDir& dir);
    const GMatrix& rot(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCTAPointing& pnt);
    void free_members(void);
    void update(void) const;

    // Protected members
    GSkyDir         m_dir;        //!< Pointing direction

    // Cached members
    mutable bool    m_has_cache;  //!< Has transformation cache
    mutable GMatrix m_Rback;      //!< Rotation matrix
};

#endif /* GCTAPOINTING_HPP */
