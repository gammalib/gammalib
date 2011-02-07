/***************************************************************************
 *          GMWLSpectrum.hpp  -  Multi-wavelength spectrum class           *
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
 * @file GMWLSpectrum.hpp
 * @brief GMWLSpectrum class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMWLSPECTRUM_HPP
#define GMWLSPECTRUM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEventCube.hpp"
#include "GMWLDatum.hpp"
#include "GEbounds.hpp"
#include "GFitsTable.hpp"
#include "GEnergy.hpp"


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
    explicit GMWLSpectrum(const std::string& filename);
    GMWLSpectrum(const GMWLSpectrum& spec);
    virtual ~GMWLSpectrum(void);

    // Operators
    GMWLSpectrum& operator= (const GMWLSpectrum& spec);

    // Implemented pure virtual methods
    void          clear(void);
    GMWLSpectrum* clone(void) const;
    int           size(void) const { return m_data.size(); }
    int           dim(void) const { return 1; }
    int           naxis(int axis) const { return m_data.size(); }
    void          load(const std::string& filename);
    GMWLDatum*    pointer(int index);
    int           number(void) const;
    std::string   print(void) const;

    // Other methods
    std::string telescope(void) const { return m_telescope; }
    std::string instrument(void) const { return m_instrument; }
    GEbounds    ebounds(void) const;
    void        load(const std::string& filename, const std::string& extname);
    void        load_fits(const std::string& filename, int extno = 0);

protected:
    // Protected methods
    void    init_members(void);
    void    copy_members(const GMWLSpectrum& spec);
    void    free_members(void);
    void    read_fits(const GFitsTable* table);
    GEnergy conv_energy(const double& energy, const std::string& unit);
    double  conv_flux(const GEnergy& energy, const double& flux, const std::string& unit);

    // Protected members
    std::string            m_telescope;   //!< Telescope name
    std::string            m_instrument;  //!< Instrument name
    std::vector<GMWLDatum> m_data;        //!< Spectral data
};

#endif /* GMWLSPECTRUM_HPP */
