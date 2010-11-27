/***************************************************************************
 *         GMWLDatum.hpp  -  Multi-wavelength spectral point class         *
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
 * @file GMWLDatum.hpp
 * @brief GMWLDatum class interface definition.
 * @author J. Knodlseder
 */

#ifndef GMWLDATUM_HPP
#define GMWLDATUM_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GLog.hpp"
#include "GEventBin.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GMWLInstDir.hpp"


/***********************************************************************//**
 * @class GMWLDatum
 *
 * @brief GMWLDatum class interface defintion
 *
 * This class defines a spectral data point for the multi-wavelength
 * interface. It derives from the abstract GEventBin base class.
 ***************************************************************************/
class GMWLDatum : public GEventBin {

    // Friend classes
    friend class GMWLSpectrum;   //!< Needs access to load data

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GMWLDatum& datum);
    friend GLog&         operator<< (GLog& log, const GMWLDatum& datum);

public:
    // Constructors and destructors
    GMWLDatum(void);
    GMWLDatum(const GMWLDatum& datum);
    virtual ~GMWLDatum(void);

    // Operators
    GMWLDatum& operator= (const GMWLDatum& datum);

    // Event access methods
    const GEnergy&  energy(void) const { return m_eng; }
    const GInstDir& dir(void) const { return m_dir; }
    const GTime&    time(void) const { return m_time; }
    double          counts(void) const { return m_flux; }
    double          error(void) const { return m_flux_err; }

    // Other methods
    void        clear(void);
    GMWLDatum*  clone(void) const;
    double      size(void) const { return 1.0; }
    GEnergy     energy_err(void) const { return m_eng_err; }
    double      flux(void) const { return m_flux; }
    double      flux_err(void) const { return m_flux_err; }
    std::string print(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLDatum& datum);
    void free_members(void);

    // Protected members
    GMWLInstDir m_dir;       //!< Instrument direction of spectral point (not used)
    GTime       m_time;      //!< Time of spectral point (not used)
    GEnergy     m_eng;       //!< Energy of spectral point
    GEnergy     m_eng_err;   //!< Uncertainty in energy
    double      m_flux;      //!< Flux of spectral point
    double      m_flux_err;  //!< Uncertainty in flux

};

#endif /* GMWLDATUM_HPP */
