/***************************************************************************
 *                   GPolarization.hpp - Polarization class                *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2021 by Juergen Knoedlseder                              *
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
 * @file GPolarization.hpp
 * @brief Polarization value class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GPOLARIZATION_HPP
#define GPOLARIZATION_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GBase.hpp"


/***********************************************************************//**
 * @class GPolarization
 *
 * @brief Class that handles polarizations.
 ***************************************************************************/
class GPolarization : public GBase {

public:
    // Constructors and destructors
    GPolarization(void);
    GPolarization(const GPolarization& polarization);
    virtual ~GPolarization(void);
 
    // Operators
    GPolarization& operator=(const GPolarization& polarization);

    // Methods
    void           clear(void);
    GPolarization* clone(void) const;
    std::string    classname(void) const;
    std::string    print(const GChatter& chatter = NORMAL) const;
  
protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GPolarization& polarization);
    void free_members(void);

    // Protected data members

};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GPolarization").
 ***************************************************************************/
inline
std::string GPolarization::classname(void) const
{
    return ("GPolarization");
}

#endif /* GPOLARIZATION_HPP */
