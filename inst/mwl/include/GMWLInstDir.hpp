/***************************************************************************
 *     GMWLInstDir.hpp - Multi-wavelength instrument direction class       *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2020 by Juergen Knoedlseder                         *
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
 * @file GMWLInstDir.hpp
 * @brief GMWLInstDir class definition.
 * @author Juergen Knoedlseder
 */

#ifndef GMWLINSTDIR_HPP
#define GMWLINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include <cstdint>
#include "GInstDir.hpp"


/***********************************************************************//**
 * @class GMWLInstDir
 *
 * @brief Interface for the Multi-wavelength instrument direction class.
 *
 * The Multi-wavelength instrument direction class is a dummy class that is
 * needed but not used for the implementation of the Multi-wavelength
 * interface as an instrument class.
 ***************************************************************************/
class GMWLInstDir : public GInstDir {

public:
    // Constructors and destructors
    GMWLInstDir(void);
    GMWLInstDir(const GMWLInstDir& dir);
    virtual ~GMWLInstDir(void);

    // Operators
    GMWLInstDir& operator= (const GMWLInstDir& dir);

    // Methods
    virtual void         clear(void);
    virtual GMWLInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual uint64_t     hash(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GMWLInstDir& dir);
    void free_members(void);
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GMWLInstDir").
 ***************************************************************************/
inline
std::string GMWLInstDir::classname(void) const
{
    return ("GMWLInstDir");
}


/***********************************************************************//**
 * @brief Return instrument direction hash value
 *
 * @return Hash value.
 *
 * Returns a hash value that can be used in the response cache.
 ***************************************************************************/
inline
uint64_t GMWLInstDir::hash(void) const
{
    return 0;
}

#endif /* GMWLINSTDIR_HPP */
