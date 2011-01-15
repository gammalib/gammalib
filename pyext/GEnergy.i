/***************************************************************************
 *                   GEnergy.i  -  Energy class python I/F                 *
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
 * @file GEnergy.i
 * @brief GEnergy class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEnergy.hpp"
#include "GTools.hpp"
%}


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
    // Operator friends
    /*
    friend GEnergy operator+ (const GEnergy &a, const GEnergy &b);
    friend GEnergy operator- (const GEnergy &a, const GEnergy &b);
    friend GEnergy operator* (const double &a, const GEnergy &b);
    friend GEnergy operator* (const GEnergy &a, const double &b);
    friend GEnergy operator/ (const GEnergy &a, const double &b);
    friend bool    operator== (const GEnergy &a, const GEnergy &b);
    friend bool    operator!= (const GEnergy &a, const GEnergy &b);
    friend bool    operator< (const GEnergy &a, const GEnergy &b);
    friend bool    operator<= (const GEnergy &a, const GEnergy &b);
    friend bool    operator> (const GEnergy &a, const GEnergy &b);
    friend bool    operator>= (const GEnergy &a, const GEnergy &b);
    */

public:
    // Constructors and destructors
    GEnergy(void);
    GEnergy(const GEnergy& eng);
    ~GEnergy(void);
 
    // Operators
    GEnergy& operator+= (const GEnergy& eng);
    GEnergy& operator-= (const GEnergy& eng);

    // Methods
    void        clear(void) { m_energy = 0.0; }
    double      erg(void) const;
    double      keV(void) const;
    double      MeV(void) const;
    double      GeV(void) const;
    double      TeV(void) const;
    void        erg(const double& eng);
    void        keV(const double& eng);
    void        MeV(const double& eng);
    void        GeV(const double& eng);
    void        TeV(const double& eng);
};


/***********************************************************************//**
 * @brief GEnergy class extension
 ***************************************************************************/
%extend GEnergy {
    char *__str__() {
        return tochar(self->print());
    }
    GEnergy copy() {
        return (*self);
    }
};
