/***************************************************************************
 *          GMWLPointing.hpp  -  Multi-wavelength pointing class           *
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
 * @file GMWLPointing.hpp
 * @brief GMWLPointing class definition.
 * @author J. Knodlseder
 */

#ifndef GMWLPOINTING_HPP
#define GMWLPOINTING_HPP

/* __ Includes ___________________________________________________________ */
#include "GPointing.hpp"


/***********************************************************************//**
 * @class GMWLPointing
 *
 * @brief Interface for the Multi-wavelength pointing class
 *
 * The Multi-wavelength pointing class is a dummy class that is needed but
 * not used for the implementation of the Multi-wavelength interface as an
 * instrument class.
 ***************************************************************************/
class GMWLPointing : public GPointing {

public:
    // Constructors and destructors
    GMWLPointing(void);
    GMWLPointing(const GMWLPointing& pnt);
    ~GMWLPointing(void);

    // Operators
    GMWLPointing& operator= (const GMWLPointing& pnt);

    // Methods
    void           clear(void);
    GMWLPointing*  clone(void) const;
    const GSkyDir& dir(void) const { return m_dir; }
    std::string    print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLPointing& pnt);
    void free_members(void);

    // Protected members
    GSkyDir m_dir;  //!< Pointing direction
};

#endif /* GMWLPOINTING_HPP */
