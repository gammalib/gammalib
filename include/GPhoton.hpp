/***************************************************************************
 *                        GPhoton.hpp - Photon class                       *
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
 * @file GPhoton.hpp
 * @brief Photon class definition.
 * @author J. Knodlseder
 */

#ifndef GPHOTON_HPP
#define GPHOTON_HPP

/* __ Includes ___________________________________________________________ */
#include <vector>
#include <string>
#include <iostream>
#include "GLog.hpp"
#include "GSkyDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


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

    // I/O friends
    friend std::ostream& operator<< (std::ostream& os, const GPhoton& ph);
    friend GLog&         operator<< (GLog& log, const GPhoton& ph);

    // Operator friends
    friend bool operator== (const GPhoton &a, const GPhoton &b);
    friend bool operator!= (const GPhoton &a, const GPhoton &b);

public:
    // Constructors and destructors
    GPhoton(void);
    GPhoton(const GPhoton& ph);
    virtual ~GPhoton(void);
 
    // Operators
    GPhoton& operator= (const GPhoton& ph);

    // Methods
    void        clear(void);
    GSkyDir&    dir(void) { return m_dir; }
    GEnergy&    energy(void) { return m_energy; }
    GTime&      time(void) { return m_time; }
    int         mcid(void) { return m_mc_id; }
    void        dir(const GSkyDir& dir) { m_dir=dir; }
    void        energy(const GEnergy& energy) { m_energy=energy; }
    void        time(const GTime& time) { m_time=time; }
    void        mcid(const int& mcid) { m_mc_id=mcid; }
    std::string print(void) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPhoton& ph);
    void free_members(void);

    // Protected data members
    GSkyDir m_dir;      //!< Photon arrival direction
    GEnergy m_energy;   //!< Photon energy
    GTime   m_time;     //!< Photon arrival time
    int     m_mc_id;    //!< Monte Carlo simulation origin
};


/***************************************************************************
 *                               Inline friends                            *
 ***************************************************************************/
inline
bool operator== (const GPhoton &a, const GPhoton &b)
{
    return (a.m_energy == b.m_energy && a.m_time == b.m_time &&
            a.m_dir.dist(b.m_dir) == 0.0);
}
inline
bool operator!= (const GPhoton &a, const GPhoton &b)
{
    return (a.m_energy != b.m_energy ||  a.m_time != b.m_time ||
            a.m_dir.dist(b.m_dir) > 0.0);
}


/***************************************************************************
 *                                 Typedefs                                *
 ***************************************************************************/
typedef std::vector<GPhoton> GPhotons;

#endif /* GPHOTON_HPP */
