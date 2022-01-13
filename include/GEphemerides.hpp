/***************************************************************************
 *                   GEphemerides.hpp - Ephemerides class                  *
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
 * @file GEphemerides.hpp
 * @brief Ephemerides class definition
 * @author Juergen Knoedlseder
 */

#ifndef GEPHEMERIDES_HPP
#define GEPHEMERIDES_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"
#include "GFilename.hpp"

/* __ Forward declarations _______________________________________________ */

/* __ Constants __________________________________________________________ */


/***********************************************************************//**
 * @class GEphemerides
 *
 * @brief Ephemerides class
 *
 * This class implements the JPL ephemerides.
 ***************************************************************************/
class GEphemerides : public GBase {

public:
    // Constructors and destructors
    GEphemerides(void);
    GEphemerides(const GEphemerides& ephemerides);
    virtual ~GEphemerides(void);

    // Operators
    GEphemerides& operator=(const GEphemerides& ephemerides);

    // Implemented pure virtual base class methods
    virtual void          clear(void);
    virtual GEphemerides* clone(void) const;
    virtual std::string   classname(void) const;
    virtual std::string   print(const GChatter& chatter = NORMAL) const;

    // Other methods
    // TODO: Add any further methods that are needed

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GEphemerides& ephemerides);
    void free_members(void);
    
    // Protected members
    std::string m_name;     //!< Ephemerides (e.g. DE200)
    GFilename   m_filename; //!< Ephemerides filename
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GEphemerides").
 ***************************************************************************/
inline
std::string GEphemerides::classname(void) const
{
    return ("GEphemerides");
}

#endif /* GEPHEMERIDES_HPP */
