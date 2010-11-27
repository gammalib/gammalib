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
#include "GInstDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GEventBin
 *
 * @brief Abstract interface for the event bin class.
 *
 * An event bin is a collection of event atoms with similar properties.
 * Event bins are used for binned analysis.
 *
 * Each event has 3 attributes: energy, instrument direction and time.
 * These attributes can be accessed and changed through the energy(),
 * dir(), and time() methods.
 *
 * The counts() and error() methods return the number of events within an
 * event bin and the uncertainty in this number, which is typically the
 * square root of the number of events.
 *
 * The size() method returns the size of an event bin, which is the
 * quantity that has to be multiplied by the probability for an event to
 * occur to predict the number of events in a bin). The size is the solid
 * angle of the event bin times the energy width times the ontime interval
 * covered by the events.
 *
 * The GEventBin class does not hold any data members. Data members are
 * stored in the derived classes.
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

    // Event access methods
    virtual const GInstDir& dir(void) const = 0;
    virtual const GEnergy&  energy(void) const = 0;
    virtual const GTime&    time(void) const = 0;
    virtual double          counts(void) const = 0;
    virtual double          error(void) const = 0;

    // Other methods
    virtual void       clear(void) = 0;
    virtual double     size(void) const = 0;
    virtual GEventBin* clone(void) const = 0;
    bool               isatom(void) const { return false; }
    bool               isbin(void) const { return true; }

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEventBin& bin);
    void free_members(void);
};

#endif /* GEVENTBIN_HPP */
