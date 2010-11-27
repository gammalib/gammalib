/***************************************************************************
 *          GMWLSpectrum.hpp  -  Multi-wavelength spectrum class           *
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
 * @file GMWLSpectrum.hpp
 * @brief GMWLSpectrum class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMWLSPECTRUM_HPP
#define GMWLSPECTRUM_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GEventCube.hpp"
#include "GMWLDatum.hpp"


/***********************************************************************//**
 * @class GMWLSpectrum
 *
 * @brief GMWLSpectrum class interface defintion
 *
 * This class defines a multi-wavelength spectrum and is a container for
 * spectral points of type GMWLDatum. It derives from the abstract
 * GEventCube base class.
 ***************************************************************************/
class GMWLSpectrum : public GEventCube {

    // Friend classes
    friend class GMWLObservation;

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GMWLSpectrum& spec);
    friend GLog&         operator<< (GLog& log, const GMWLSpectrum& spec);

public:
    // Constructors and destructors
    GMWLSpectrum(void);
    GMWLSpectrum(const std::string& filename);
    GMWLSpectrum(const GMWLSpectrum& spec);
    virtual ~GMWLSpectrum(void);

    // Operators
    GMWLSpectrum& operator= (const GMWLSpectrum& spec);

    // Implemented pure virtul methods
    void          clear(void);
    GMWLSpectrum* clone(void) const;
    void          load(const std::string& filename);
    GMWLDatum*    pointer(int index);
    int           number(void) const;

    // Other methods
    int         size(void) const;
    std::string print(void) const;
    void        load_fits(const std::string& filename, int extno = 0);
    void        load_fits(const std::string& filename, const std::string& extname);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLSpectrum& spec);
    void free_members(void);
    void read_fits(const GFitsTable* table);

    // Protected members
    std::vector<GMWLDatum> m_data;  //!< Spectral data
};

#endif /* GMWLSPECTRUM_HPP */
