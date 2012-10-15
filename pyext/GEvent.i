/***************************************************************************
 *               GEvent.i  -  Abstract event class python I/F              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2012 by Juergen Knoedlseder                         *
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
 * @file GEvent.i
 * @brief GEvent class python interface
 * @author Juergen Knoedlseder
 */
%{
/* Put headers and other declarations here that are needed for compilation */
#include "GEvent.hpp"
#include "GTools.hpp"
%}


/***********************************************************************//**
 * @class GEvent
 *
 * @brief Abstract interface for the event classes.
 *
 * This class provides an abstract interface to an event. An event can be
 * either a physical event occuring in the instrument (called an event
 * atom) or a collection of events with similar properties (called an
 * event bin). While event atoms are used for unbinned analysis, event bins
 * are used for binned analysis. The methods isatom() and isbin() inform
 * whether an event is an atom or a bin.
 *
 * Each event has 3 attributes: energy, instrument direction and time.
 * These attributes can be accessed and changed through the energy(),
 * dir(), and time() methods.
 * 
 * The counts() and error() methods return the number of counts and the
 * error in this number for each event. For event atoms the number of
 * counts is 1 and the error is 0. The event bins, the number of counts
 * is the number of events within a bin, and error is the uncertainty in
 * this number (typically the square root of the number, yet also other
 * schemes may be implemented).
 *
 * The size() method returns the size of an event bin, which is the
 * quantity that has to be multiplied by the probability for an event to
 * occur to predict the number of events in a bin. For event atoms this
 * quantity is by definition 1. For event bins, the size is the solid
 * angle of the event bin times the energy width times the ontime interval
 * covered by the events.
 *
 * The GEvent class does not hold any data members. Data members are stored
 * in the derived classes.
 ***************************************************************************/
class GEvent : public GBase {

public:
    // Constructors and destructors
    GEvent(void);
    GEvent(const GEvent& event);
    virtual ~GEvent(void);

    // Pure virtual methods
    virtual void            clear(void) = 0;
    virtual GEvent*         clone(void) const = 0;
    virtual double          size(void) const = 0;
    virtual const GInstDir& dir(void) const = 0;
    virtual const GEnergy&  energy(void) const = 0;
    virtual const GTime&    time(void) const = 0;
    virtual double          counts(void) const = 0;
    virtual double          error(void) const = 0;
    virtual bool            isatom(void) const = 0;
    virtual bool            isbin(void) const = 0;
};


/***********************************************************************//**
 * @brief GEvent class extension
 ***************************************************************************/
%extend GEvent {
    char *__str__() {
        return tochar(self->print());
    }
};
