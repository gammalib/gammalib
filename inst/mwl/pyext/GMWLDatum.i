/***************************************************************************
 *     GMWLDatum.i  -  Multi-wavelength spectral point class python I/F    *
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
 * @file GMWLDatum.i
 * @brief GMWLDatum class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GMWLDatum.hpp"
%}


/***********************************************************************//**
 * @class GMWLDatum
 *
 * @brief GMWLDatum class interface defintion
 *
 * This class defines a spectral data point for the multi-wavelength
 * interface. It derives from the abstract GEventBin base class.
 ***************************************************************************/
class GMWLDatum : public GEventBin {
public:
    // Constructors and destructors
    GMWLDatum(void);
    GMWLDatum(const GMWLDatum& datum);
    virtual ~GMWLDatum(void);

    // Event access methods
    const GEnergy&  energy(void) const { return m_eng; }
    const GInstDir& dir(void) const { return m_dir; }
    const GTime&    time(void) const { return m_time; }
    double          counts(void) const { return m_flux; }
    double          error(void) const;

    // Other methods
    void        clear(void);
    GMWLDatum*  clone(void) const;
    double      size(void) const { return 1.0; }
    GEnergy     energy_err(void) const { return m_eng_err; }
    double      flux(void) const { return m_flux; }
    double      flux_err(void) const { return m_flux_err; }
};


/***********************************************************************//**
 * @brief GMWLDatum class extension
 ***************************************************************************/
%extend GMWLDatum {
    char *__str__() {
        static std::string result = self->print();
        return ((char*)result.c_str());
    }
    GMWLDatum copy() {
        return (*self);
    }
};
