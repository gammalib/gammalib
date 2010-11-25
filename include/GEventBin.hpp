/***************************************************************************
 *               GEventBin.hpp  -  Event bin abstract base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2010 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEventBin.hpp
 * @brief GEventBin abstract base class interface definition.
 * @author J. Knodlseder
 */

#ifndef GEVENTBIN_HPP
#define GEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GEvent.hpp"
#include "GModels.hpp"
#include "GVector.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GInstDir.hpp"
#include "GPointing.hpp"
#include "GResponse.hpp"


/***********************************************************************//**
 * @class GEventBin
 *
 * @brief Abstract interface for the event bin class.
 *
 * An event bin is a collection of event atoms with similar properties.
 * It has three generic attributes (m_counts, m_time and m_energy) that can
 * be accessed through the counts(), energy() and time() methods. The
 * counts() method returns a double precision value, while the energy() and
 * time() methods return const pointers.
 *
 * The event bin class does not actually allocate memory for event bins
 * but handles pointers that point the information that is relevant for
 * analysis. Filling of the information is handled by the pointer() method of
 * the corresponding event cube.
 ***************************************************************************/
class GEventBin : public GEvent {

    // Friend classes
    friend class GEvents;

public:
    // Constructors and destructors
    GEventBin(void);
    GEventBin(const GEventBin& bin);
    virtual ~GEventBin(void);

    // Operators
    virtual GEventBin& operator= (const GEventBin& bin);

    // Pure virtual methods
    virtual double          size(void) const = 0;
    virtual const GInstDir* dir(void) const = 0;
    virtual GEventBin*      clone(void) const = 0;

    // Implemented methods
    bool           isatom(void) const { return false; }
    bool           isbin(void) const { return true; }
    double         counts(void) const { return *m_counts; }
    const GEnergy* energy(void) const { return m_energy; }
    const GTime*   time(void) const { return m_time; }

protected:
    // Protected data area
    double*  m_counts;      //!< Pointer to number of counts
    GTime*   m_time;        //!< Pointer to bin time
    GEnergy* m_energy;      //!< Pointer to bin energy

private:
    // Provate methods
    void init_members(void);
    void copy_members(const GEventBin& bin);
    void free_members(void);

};

#endif /* GEVENTBIN_HPP */
