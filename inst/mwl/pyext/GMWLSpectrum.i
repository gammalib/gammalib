/***************************************************************************
 *       GMWLSpectrum.i  -  Multi-wavelength spectrum class python I/F     *
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
 * @file GMWLSpectrum.i
 * @brief GMWLSpectrum class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMWLSpectrum.hpp"
%}


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
public:
    // Constructors and destructors
    GMWLSpectrum(void);
    GMWLSpectrum(const std::string& filename);
    GMWLSpectrum(const GMWLSpectrum& spec);
    virtual ~GMWLSpectrum(void);

    // Implemented pure virtul methods
    void          clear(void);
    GMWLSpectrum* clone(void) const;
    void          load(const std::string& filename);
    GMWLDatum*    pointer(int index);
    int           number(void) const;

    // Other methods
    int         size(void) const;
    std::string telescope(void) const { return m_telescope; }
    std::string instrument(void) const { return m_instrument; }
    GEbounds    ebounds(void) const;
    void        load(const std::string& filename, const std::string& extname);
    void        load_fits(const std::string& filename, int extno = 0);
};


/***********************************************************************//**
 * @brief GMWLSpectrum class extension
 ***************************************************************************/
%extend GMWLSpectrum {
    GMWLSpectrum copy() {
        return (*self);
    }
};
