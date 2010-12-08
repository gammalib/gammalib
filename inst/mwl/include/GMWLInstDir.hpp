/***************************************************************************
 *     GMWLInstDir.hpp  -  Multi-wavelength instrument direction class     *
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
 * @file GMWLInstDir.hpp
 * @brief GMWLInstDir class definition.
 * @author J. Knodlseder
 */

#ifndef GMWLINSTDIR_HPP
#define GMWLINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"


/***********************************************************************//**
 * @class GMWLInstDir
 *
 * @brief Interface for the Multi-wavelength instrument direction class.
 *
 * The Multi-wavelength instrument direction class is a dummy class that is
 * needed but not used for the implementation of the Multi-wavelength
 * interface as an instrument class.
 ***************************************************************************/
class GMWLInstDir : public GInstDir {

public:
    // Constructors and destructors
    GMWLInstDir(void);
    GMWLInstDir(const GMWLInstDir& dir);
    virtual ~GMWLInstDir(void);

    // Operators
    GMWLInstDir& operator= (const GMWLInstDir& dir);

    // Methods
    void         clear(void);
    GMWLInstDir* clone(void) const;
    std::string  print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLInstDir& dir);
    void free_members(void);
};

#endif /* GMWLINSTDIR_HPP */
