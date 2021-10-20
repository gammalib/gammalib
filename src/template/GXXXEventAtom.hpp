/***************************************************************************
 *             GXXXEventAtom.hpp - [INSTRUMENT] event atom class           *
 * ----------------------------------------------------------------------- *
 *  copyright (C) [YEAR] by [AUTHOR]                                       *
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
 * @file GXXXEventAtom.hpp
 * @brief [INSTRUMENT] event atom class definition
 * @author [AUTHOR]
 */

#ifndef GXXXEVENTATOM_HPP
#define GXXXEVENTATOM_HPP

/* __ Includes ___________________________________________________________ */
#include <iostream>
#include "GEventAtom.hpp"
#include "GXXXInstDir.hpp"
#include "GEnergy.hpp"
#include "GTime.hpp"
#include "GPolarization.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GXXXEventAtom
 *
 * @brief [INSTRUMENT] event atom class
 ***************************************************************************/
class GXXXEventAtom : public GEventAtom {

public:
    // Constructors and destructors
    GXXXEventAtom(void);
    GXXXEventAtom(const GXXXEventAtom& atom);
    virtual ~GXXXEventAtom(void);

    // Operators
    GXXXEventAtom& operator=(const GXXXEventAtom& atom);

    // Implemented pure virtual base class methods
    virtual void                 clear(void);
    virtual GXXXEventAtom*       clone(void) const;
    virtual std::string          classname(void) const;
    virtual const GXXXInstDir&   dir(void) const;
    virtual const GEnergy&       energy(void) const;
    virtual const GTime&         time(void) const;
    virtual const GPolarization& polarization(void) const;
    virtual std::string          print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GXXXEventAtom& atom);
    void free_members(void);

    // Protected members
    GXXXInstDir   m_dir;          //!< Event direction
    GEnergy       m_energy;       //!< Event energy
    GTime         m_time;         //!< Event time
    GPolarization m_polarization; //!< Event polarization
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GXXXEventAtom").
 ***************************************************************************/
inline
std::string GXXXEventAtom::classname(void) const
{
    return ("GXXXEventAtom");
}


/***********************************************************************//**
 * @brief Return event instrument direction
 *
 * @return Event instrument direction.
 *
 * Returns the direction of the event.
 *
 * @todo Specify what actually an event direction is. Your instrument may
 * not be an imaging instrument.
 ***************************************************************************/
inline
const GXXXInstDir& GXXXEventAtom::dir(void) const
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
const GEnergy& GXXXEventAtom::energy(void) const
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
const GTime& GXXXEventAtom::time(void) const
{
    return m_time;
}


/***********************************************************************//**
 * @brief Return event polarization
 *
 * @return Event polarization.
 *
 * Returns the event polarization.
 ***************************************************************************/
inline
const GPolarization& GXXXEventAtom::polarization(void) const
{
    return m_polarization;
}

#endif /* GXXXEVENTATOM_HPP */
