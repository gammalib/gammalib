/***************************************************************************
 *                 GLATEventBin.hpp  -  LAT event bin class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
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
#include <ostream>
#include "GEventBin.hpp"
#include "GModels.hpp"
#include "GVector.hpp"
#include "GLATInstDir.hpp"
#include "GLATPointing.hpp"
#include "GLATResponse.hpp"


/***********************************************************************//**
 * @class GLATEventBin
 *
 * @brief GLATEventBin class interface defintion.
 ***************************************************************************/
class GLATEventBin : public GEventBin {

    // Friend classes
    friend class GLATEventCube;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GLATEventBin& bin);

public:
    // Constructors and destructors
    GLATEventBin(void);
    GLATEventBin(const GLATEventBin& bin);
    virtual ~GLATEventBin(void);

    // Operators
    GLATEventBin& operator= (const GLATEventBin& bin);

    // Methods
    double              model(GModels& models, GVector* gradient = NULL) const;
    double              size(void) const;
    const GLATInstDir*  dir(void) const { return m_dir; }
    const GLATPointing* pnt(void) const { return m_pnt; }
    const GLATResponse* rsp(void) const { return m_rsp; }

protected:
    // Protected methods
    void          init_members(void);
    void          copy_members(const GLATEventBin& bin);
    void          free_members(void);
    GLATEventBin* clone(void) const;

    // LAT specific event attributes
    GLATInstDir*  m_dir;     //!< Pointer to event direction
    GLATPointing* m_pnt;     //!< Pointer to instrument pointing
    GLATResponse* m_rsp;     //!< Pointer to instrument response function

};

#endif /* GLATEVENTBIN_HPP */
