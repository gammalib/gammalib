/***************************************************************************
 *          GLATInstDir.hpp - Fermi/LAT instrument direction class         *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2014 by Juergen Knoedlseder                         *
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
 * @file GLATInstDir.hpp
 * @brief Fermi/LAT instrument direction class definition
 * @author Juergen Knoedlseder
 */

#ifndef GLATINSTDIR_HPP
#define GLATINSTDIR_HPP

/* __ Includes ___________________________________________________________ */
#include <string>
#include "GInstDir.hpp"
#include "GSkyDir.hpp"


/***********************************************************************//**
 * @class GLATInstDir
 *
 * @brief Fermi/LAT instrument direction class
 *
 * The LAT instrument direction is an encapsulation of a sky direction
 * as LAT is an imaging device.
 ***************************************************************************/
class GLATInstDir : public GInstDir {

public:
    // Constructors and destructors
    GLATInstDir(void);
    explicit GLATInstDir(const GSkyDir& dir);
    GLATInstDir(const GLATInstDir& dir);
    virtual ~GLATInstDir(void);

    // Operators
    GLATInstDir& operator=(const GLATInstDir& dir);

    // Implemented pure virtual base class methods
    virtual void         clear(void);
    virtual GLATInstDir* clone(void) const;
    virtual std::string  classname(void) const;
    virtual std::string  print(const GChatter& chatter = NORMAL) const;

    // Other methods
    void           dir(const GSkyDir& dir);
    GSkyDir&       dir(void);
    const GSkyDir& dir(void) const;

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const GLATInstDir& dir);
    void free_members(void);
    
    // Data members
    GSkyDir m_dir;  //!< Observed incident direction of event
};


/***********************************************************************//**
 * @brief Return class name
 *
 * @return String containing the class name ("GLATInstDir").
 ***************************************************************************/
inline
std::string GLATInstDir::classname(void) const
{
    return ("GLATInstDir");
}


/***********************************************************************//**
 * @brief Returns reference to sky direction
 *
 * @return Reference to sky direction.
 *
 * Returns reference to the sky direction.
 ***************************************************************************/
inline
GSkyDir& GLATInstDir::dir(void)
{
    return (m_dir);
}


/***********************************************************************//**
 * @brief Returns reference to sky direction (const version)
 *
 * @return Reference to sky direction.
 *
 * Returns reference to the sky direction.
 ***************************************************************************/
inline
const GSkyDir& GLATInstDir::dir(void) const
{
    return (m_dir);
}


/***********************************************************************//**
 * @brief Set sky direction
 *
 * @param[in] dir Sky direction.
 *
 * Set the sky direction.
 ***************************************************************************/
inline
void GLATInstDir::dir(const GSkyDir& dir)
{
    m_dir = dir;
    return;
}

#endif /* GLATINSTDIR_HPP */
