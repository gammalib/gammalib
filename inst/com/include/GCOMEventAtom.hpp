/***************************************************************************
 *               GCOMEventAtom.hpp - COMPTEL event atom class              *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2017 by Juergen Knoedlseder                              *
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
 * @file GCOMEventAtom.hpp
 * @brief COMPTEL event atom class definition
 * @author Juergen Knoedlseder
 */

#ifndef GCOMEVENTATOM_HPP
#define GCOMEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GEventAtom.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GCOMInstDir.hpp"

/* __ Forward declarations _______________________________________________ */


/***********************************************************************//**
 * @class GCOMEventAtom
 *
 * @brief COMPTEL event atom class
 ***************************************************************************/
class GCOMEventAtom : public GEventAtom {

public:
    // Constructors and destructors
    GCOMEventAtom(void);
    GCOMEventAtom(const GCOMEventAtom& atom);
    virtual ~GCOMEventAtom(void);

    // Operators
    GCOMEventAtom& operator=(const GCOMEventAtom& atom);

    // Implemented pure virtual base class methods
    void               clear(void);
    GCOMEventAtom*     clone(void) const;
    std::string        classname(void) const;
    const GCOMInstDir& dir(void) const;
    const GEnergy&     energy(void) const;
    const GTime&       time(void) const;
    std::string        print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GCOMEventAtom& atom);
    void free_members(void);

    // Protected members
    GCOMInstDir m_dir;    //!< Event direction
    GEnergy     m_energy; //!< Event energy
    GTime       m_time;   //!< Event time
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GCOMEventAtom").
 ***************************************************************************/
inline
std::string GCOMEventAtom::classname(void) const
{
    return ("GCOMEventAtom");
}


/***********************************************************************//**
 * @brief Return event instrument direction
 *
 * @return Event instrument direction.
 *
 * Returns the direction of the event.
 *
 * @todo Specify what actually an event direction is. You instrument may
 * not be an imaging instrument.
 ***************************************************************************/
inline
const GCOMInstDir& GCOMEventAtom::dir(void) const
{
    return m_dir;
}


/***********************************************************************//**
 * @brief Return event energy
 *
 * @return Event energy.
 *
 * Returns the reconstructed energy of the photon on the sky.
 ***************************************************************************/
inline
const GEnergy& GCOMEventAtom::energy(void) const
{
    return m_energy;
}


/***********************************************************************//**
 * @brief Return event time
 *
 * @return Event time.
 *
 * Returns the event triggering time.
 ***************************************************************************/
inline
const GTime& GCOMEventAtom::time(void) const
{
    return m_time;
}

#endif /* GCOMEVENTATOM_HPP */
