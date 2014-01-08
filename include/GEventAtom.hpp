/***************************************************************************
 *               GEventAtom.hpp - Abstract event atom base class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2009-2013 by Jurgen Knodlseder                           *
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
 * @file GEventAtom.hpp
 * @brief Abstract event atom base class definition
 * @author Juergen Knoedlseder
 */

#ifndef GEVENTATOM_HPP
#define GEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GEvent.hpp"
#include "GInstDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"


/***********************************************************************//**
 * @class GEventAtom
 *
 * @brief Abstract interface for the event atom class.
 *
 * An event atom is a single event occuring in an instrument. Event atoms
 * are used for unbinned analysis.
 *
 * Each event has 3 attributes: energy, instrument direction and time.
 * These attributes can be accessed and changed through the energy(),
 * dir(), and time() methods.
 *
 * The counts() and error() methods return the number of counts and the
 * error in this number for each event. For event atoms the number of
 * counts is 1 and the error is 0.
 *
 * The size() method returns the size of an event bin ,which is the
 * quantity that has to be multiplied by the probability for an event to
 * occur to predict the number of events in a bin. For event atoms this
 * quantity is by definition 1.
 *
 * The GEventAtom class does not hold any data members. Data members are
 * stored in the derived classes.
 ***************************************************************************/
class GEventAtom : public GEvent {

    // Friend classes
    friend class GEvents;

public:
    // Constructors and destructors
    GEventAtom(void);
    GEventAtom(const GEventAtom& atom);
    virtual ~GEventAtom(void);

    // Operators
    virtual GEventAtom& operator=(const GEventAtom& atom);

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GEvent*         clone(void) const = 0;
    virtual double          size(void) const;
    virtual const GInstDir& dir(void) const = 0;
    virtual const GEnergy&  energy(void) const = 0;
    virtual const GTime&    time(void) const = 0;
    virtual double          counts(void) const;
    virtual double          error(void) const;
    virtual std::string     print(const GChatter& chatter = NORMAL) const = 0;

    // Other methods
    bool is_atom(void) const;
    bool is_bin(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEventAtom& atom);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return size of an event atom
 *
 * @return Event size, fixed to 1.0 for an event atom.
 *
 * Returns the size of an event atom. The size is only useful for event bins,
 * for which it is defined as the quantity that needs to be multiplied by the
 * event probability to give the predicted number of events in a bin. For
 * event atoms, the size is fixed to 1.0.
 ***************************************************************************/
inline
double GEventAtom::size(void) const
{
    return (1.0);
}


/***********************************************************************//**
 * @brief Return number of counts in event atom
 *
 * @return Number of counts, fixed to 1.0 for an event atom.
 ***************************************************************************/
inline
double GEventAtom::counts(void) const
{
    return (1.0);
}


/***********************************************************************//**
 * @brief Return error in number of counts in event atom
 *
 * @return Error in number of counts, fixed to 0.0 for an event atom.
 ***************************************************************************/
inline
double GEventAtom::error(void) const
{
    return (0.0);
}


/***********************************************************************//**
 * @brief Signal if event is an atom
 *
 * @return True.
 ***************************************************************************/
inline
bool GEventAtom::is_atom(void) const
{
    return (true);
}


/***********************************************************************//**
 * @brief Signal if event is a bin
 *
 * @return False.
 ***************************************************************************/
inline
bool GEventAtom::is_bin(void) const
{
    return (false);
}

#endif /* GEVENTATOM_HPP */
