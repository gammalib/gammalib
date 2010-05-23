/***************************************************************************
 *                        GEnergy.hpp - Energy class                       *
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
 * @file GEnergy.hpp
 * @brief Energy value class definition.
 * @author J. Knodlseder
 */

#ifndef GENERGY_HPP
#define GENERGY_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>


/***********************************************************************//**
 * @class GEnergy
 *
 * @brief Class that handles energies in a unit independent way.
 *
 * The GEnergy class stores an energy value in units of MeV and implements
 * methods that provide automatic conversion of the energy values in other
 * units. This makes instrument specific implementations more robust and
 * reduces the risk of unit errors.
 ***************************************************************************/
class GEnergy {

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GEnergy& eng);

public:
    // Constructors and destructors
    GEnergy(void);
    GEnergy(const GEnergy& eng);
    ~GEnergy(void);
 
    // Operators
    GEnergy& operator= (const GEnergy& eng);

    // Methods
    void   clear(void) { m_energy = 0.0; }
    double keV(void) const;
    double MeV(void) const;
    double GeV(void) const;
    double TeV(void) const;
    void   keV(const double& eng);
    void   MeV(const double& eng);
    void   GeV(const double& eng);
    void   TeV(const double& eng);
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEnergy& eng);
    void free_members(void);

    // Protected data members
    double m_energy;          //!< Energy in MeV
};

#endif /* GENERGY_HPP */
