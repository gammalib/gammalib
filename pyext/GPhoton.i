/***************************************************************************
 *                   GPhoton.i  -  Photon class python I/F                 *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2011 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GPhoton.i
 * @brief GPhoton class python interface
 * @author J. Knodlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GPhoton.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GPhoton
 *
 * @brief Class that handles photons.
 *
 * The GPhoton class stores the physical attributes of a photon such as the
 * photon arrival direction, its energy and its arrival time. This class is
 * mainly used for Monte Carlo simulations.
 ***************************************************************************/
class GPhoton {
public:
    // Constructors and destructors
    GPhoton(void);
    GPhoton(const GPhoton& ph);
    ~GPhoton(void);
 
    // Methods
    void        clear(void);
    GSkyDir&    dir(void) { return m_dir; }
    GEnergy&    energy(void) { return m_energy; }
    GTime&      time(void) { return m_time; }
    void        dir(const GSkyDir& dir) { m_dir=dir; }
    void        energy(const GEnergy& energy) { m_energy=energy; }
    void        time(const GTime& time) { m_time=time; }
};


/***********************************************************************//**
 * @brief GPhoton class extension
 ***************************************************************************/
%extend GPhoton {
    char *__str__() {
        return tochar(self->print());
    }
    bool __eq__(const GPhoton& ph) const {
        return ((*self) == ph);
    }
    bool __ne__(const GPhoton& ph) const {
        return ((*self) != ph);
    }
    GPhoton copy() {
        return (*self);
    }
};


/***********************************************************************//**
 * @brief Define GPhotons vector
 ***************************************************************************/
//#pragma SWIG nowarn=512
//%warnfilter(512) GTimes;
%include "std_vector.i"
namespace std {
    %template(GPhotons) vector<GPhoton>;
}
