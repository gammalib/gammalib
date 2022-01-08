/***************************************************************************
 *            GPulsarEphemerides.hpp - Pulsar ephemerides class            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2022 by Juergen Knoedlseder                              *
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
 * @file GPulsarEphemerides.hpp
 * @brief Pulsar ephemerides class definition
 * @author Juergen Knoedlseder
 */

#ifndef GPULSAREPHEMERIDES_HPP
#define GPULSAREPHEMERIDES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GPulsarEphemerides
 *
 * @brief Pulsar ephemerides class
 *
 * @todo Add class description.
 ***************************************************************************/
class GPulsarEphemerides : public GBase {

public:
    // Constructors and destructors
    GPulsarEphemerides(void);
    GPulsarEphemerides(const GPulsarEphemerides& ephemerides);
    virtual ~GPulsarEphemerides(void);

    // Operators
    GPulsarEphemerides& operator=(const GPulsarEphemerides& ephemerides);

    // Implemented pure virtual base class methods
    virtual void                clear(void);
    virtual GPulsarEphemerides* clone(void) const;
    virtual std::string         classname(void) const;
    virtual std::string         print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPulsarEphemerides& ephemerides);
    void free_members(void);
    
    // Protected members
    // TODO: Add any data members that are necessary
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPulsarEphemerides").
 ***************************************************************************/
inline
std::string GPulsarEphemerides::classname(void) const
{
    return ("GPulsarEphemerides");
}

#endif /* GPULSAREPHEMERIDES_HPP */
