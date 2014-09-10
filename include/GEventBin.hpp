/***************************************************************************
 *                GEventBin.hpp - Abstract event bin base class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2014 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *  This program is free software: you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation, either version 3 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
 *                                                                         *
 ***************************************************************************/
/**
 * @file GEventBin.hpp
 * @brief Abstract event bin base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GEVENTBIN_HPP
#define GEVENTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
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
 * event bin and the uncertainty in this number, which by default is
 * given by the square root of the number of events (this is the default
 * implementation provided by this class).
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
    virtual GEventBin& operator=(const GEventBin& bin);

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GEvent*         clone(void) const = 0;
    virtual std::string     classname(void) const = 0;
    virtual double          size(void) const = 0;
    virtual const GInstDir& dir(void) const = 0;
    virtual const GEnergy&  energy(void) const = 0;
    virtual const GTime&    time(void) const = 0;
    virtual double          counts(void) const = 0;
    virtual double          error(void) const = 0;
    virtual void            counts(const double& counts) = 0;
    virtual std::string     print(const GChatter& chatter = NORMAL) const = 0;

    // Other methods
    bool is_atom(void) const;
    bool is_bin(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEventBin& bin);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Signal if event is an atom
 *
 * @return False.
 ***************************************************************************/
inline
bool GEventBin::is_atom(void) const
{
    return (false);
}


/***********************************************************************//**
 * @brief Signal if event is a bin
 *
 * @return True.
 ***************************************************************************/
inline
bool GEventBin::is_bin(void) const
{
    return (true);
}

#endif /* GEVENTBIN_HPP */
